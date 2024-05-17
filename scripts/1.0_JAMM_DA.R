library(tidyverse)
library(plyranges)
library(glue)
library(parallel)
library(DESeq2)
library(BiocParallel)
library(dplyr)
library(glue)
library(patchwork)

#setwd("/data/share/htp/pleiotropy/paper_data/")

# get coverage from bed files -------------------------------
# the input files for this first step are too large to upload anywhere
# everything that follows will be reproducible
jamm_path<-"/data/ngs/epigenetics/dhs/bed_files_for_jamm"
tissues <- list.files(jamm_path)

tmp<- readRDS("roadmap_DHS_summaries/region_summary/jamm_region_info_full.rds") %>%  
  dplyr::rename(seqnames = chromosome) %>% 
  as_granges

tmp  %>%  write_bed("DHS.bed")

all_cnt<- mclapply(tissues, function(tt){
  print(tt)
  sapply( list.files(path=glue("{jamm_path}/{tt}/"), pattern="SRX"), 
          function(srx){
            print(srx)
            system( glue("bedtools coverage -counts -a DHS.bed -b {jamm_path}/{tt}/{srx} | cut -f7"), intern =T )
          })}, mc.cores = 12 )

names(all_cnt)<-tissues

saveRDS(all_cnt, "DA_JAMM/dhs_counts.RDS")
system("rm DHS.bed")






# Make count matrix -------------------------------------------------------

jamm<-readRDS("roadmap_DHS_summaries/region_summary/jamm_region_info_full.rds") 
all_cnt<-readRDS("DA_JAMM/dhs_counts.RDS")

DHScnts<-list()

DHScnts$inf<-lapply(1:length(all_cnt),function(i){
  tibble( tissue=names(all_cnt)[i],
          sample_id=colnames(all_cnt[[i]])) }) %>% bind_rows()

DHScnts$mat<- apply( cbind(all_cnt$adrenal_gland ,
                           all_cnt$brain ,
                           all_cnt$heart ,
                           all_cnt$kidney ,
                           all_cnt$large_intestine ,
                           all_cnt$lung ,
                           all_cnt$muscle ,
                           all_cnt$stomach ,
                           all_cnt$thymus ),2, as.numeric)
rownames(DHScnts$mat)<-paste0("region_id_",jamm$region_id)

DHScnts$inf<-data.frame(DHScnts$inf)
rownames(DHScnts$inf)<-DHScnts$inf$sample_id
DHScnts$inf<-DHScnts$inf[colnames(DHScnts$mat),]



# Run DESeq2 using 0+tissue in the design ---------------------------------

bigdds<-DESeqDataSetFromMatrix(countData = DHScnts$mat,
                               colData = DHScnts$inf,
                               design = ~0+tissue)

register(MulticoreParam(10))
bigdds<-DESeq(bigdds, parallel = T)
saveRDS(bigdds, "DA_JAMM/dds.rds")


bigdds<-readRDS("DA_JAMM/dds.rds")
resultsNames(bigdds)


# res<-results(bigdds, contrast=c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,-1,0.125), independentFiltering = FALSE, parallel = TRUE)
# plotMA(res)


# Do DA testing --------------------------------------------------------------

# since we have 0 and 1 possible at each position, to get all combinations, we need to generate a grid of 2^9
comb <- expand.grid(rep(list(0:1), 9))
names(comb)<-c("adrenal_gland","brain", "heart", "kidney", "large_intestine", "lung", "muscle", "stomach", "thymus")
# exclude complete 0s and 1s.
comb<-comb[2:(dim(comb)[1]-1),]



# Reoccuring summarization where CREs coming from certain tissue combination are selected and their DA status summarized

summarize_DA<-function(res, jamm, comb_i, collect=F){
  
  reswrangled<-
    res %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("region_id") %>% 
    dplyr::mutate(region_id = as.double(gsub("region_id_","",region_id))) %>% 
    dplyr::inner_join(jamm %>% dplyr::select(region_id, total, adrenal_gland, brain, 
                                             heart,  kidney, large_intestine, lung, 
                                             muscle, stomach, thymus)) %>% 
    # select only peaks of interest for the specific context
    dplyr::filter(total==sum(comb_i),
                  adrenal_gland == comb_i$adrenal_gland,
                  brain == comb_i$brain, 
                  heart == comb_i$heart, 
                  kidney == comb_i$kidney, 
                  large_intestine == comb_i$large_intestine, 
                  lung == comb_i$lung, 
                  muscle == comb_i$muscle, 
                  stomach == comb_i$stomach, 
                  thymus == comb_i$thymus) %>% 
    # do p-value adjustment
    dplyr::mutate(padj_adj = p.adjust(pvalue, "BH")) 
  
  if (collect){
    wrong<-reswrangled %>% filter(padj_adj>0.1) 
    saveRDS(wrong, paste0("DA_JAMM/to_filter_out/combi_",paste(comb_i, collapse = ""),".rds"))
  }
  
  reswrangled %>% 
    plyranges::group_by(total, adrenal_gland, brain, heart, kidney, 
                        large_intestine, lung, muscle, stomach, thymus) %>% 
    dplyr::summarise(sign_lfc_pos = length(region_id[padj_adj<0.1 & log2FoldChange>0]),
                     nonsign_lfc_pos = length(region_id[padj_adj>=0.1 & log2FoldChange>=0]),
                     sign_lfc_neg = length(region_id[padj_adj<0.1 & log2FoldChange<0]),
                     nonsign_lfc_neg = length(region_id[padj_adj>=0.1 & log2FoldChange<=0]),
                     n_regions = length(region_id))
  
}




DA_out<-lapply(1:nrow(comb), function(i){
  print(i)
  comb_i<-comb[i,]
  fract1<-1/sum(comb_i)
  fract0<-1/(length(resultsNames(bigdds))-sum(comb_i))
  
  vec<-sapply(comb_i, function(x){if (x==1){return(fract1)} else {return(-fract0)}}) 
  
  register(MulticoreParam(10))
  res<-results(bigdds, contrast=vec, independentFiltering = FALSE, parallel = TRUE, altHypothesis = "greater")
  
  # the reviewer wants us to only test the CREs that were assigned a specific label
  summarize_DA(res, jamm, comb_i, collect=T)
  
})

saveRDS(DA_out, "DA_JAMM/DA_out_1sided.rds")






# Do the same 10 times with shuffled labels -------------------------------

list_of_shuffles<-list()
for (i in 1:10){
  list_of_shuffles[[i]]<-lapply(1:nrow(comb), function(i){
    print(i)
    comb_i<-comb[i,]
    fract1<-1/sum(comb_i)
    fract0<-1/(length(resultsNames(bigdds))-sum(comb_i))
    
    vec<-sapply(comb_i, function(x){if (x==1){return(fract1)} else {return(-fract0)}})
    shuffled<-sample(vec)
    
    register(MulticoreParam(10))
    res<-results(bigdds, contrast=shuffled, independentFiltering = FALSE, parallel = TRUE, altHypothesis = "greater")
    
    summarize_DA(res, jamm, comb_i)
  }) %>%
    bind_rows()
}

saveRDS(list_of_shuffles, "DA_JAMM/DA_shuffled_1sided.rds")






# Collect ambiguous CREs and plot ----------------------------------------------

# collect ambiguous CREs (the output rds is available)
to_remove<-list.files("DA_JAMM/to_filter_out/", full.names = T)
to_remove_df<-lapply(to_remove, function(x){
  readRDS(x)
}) %>% bind_rows()
saveRDS(to_remove_df %>% dplyr::select(region_id, total), "DA_JAMM/to_filter_out.rds")


# plot error rates
DA_out<-readRDS("DA_JAMM/DA_out_1sided.rds") %>% 
  bind_rows() %>% 
  rownames_to_column("nn") %>% 
  dplyr::mutate(prop_sign_lfc_pos = sign_lfc_pos/n_regions,
                prop_lfc_pos = (sign_lfc_pos+nonsign_lfc_pos)/n_regions,
                prop_sign_lfc_neg = sign_lfc_neg/n_regions,
                total = as.factor(total)) 

DA_out<-DA_out %>% 
  left_join(DA_out %>% 
              pivot_longer(3:11, names_to = "tissue", values_to = "OPEN") %>% 
              filter(OPEN == 1) %>% 
              group_by(nn) %>% 
              reframe( open_in= glue_collapse(tissue,sep=","))
  ) %>% 
  left_join(DA_out %>% 
              pivot_longer(3:11, names_to = "tissue", values_to = "OPEN") %>% 
              filter(OPEN == 0) %>% 
              group_by(nn) %>% 
              reframe( closed_in= glue_collapse(tissue,sep=","))
  ) %>% 
  dplyr::mutate(label = ifelse(as.numeric(as.character(total))<5,open_in, closed_in))


DA_out<-DA_out %>% 
  group_by(total) %>% 
  dplyr::mutate(median_sign_pos_prop = median(prop_sign_lfc_pos)) %>% 
  ungroup() %>% 
  left_join(DA_out %>% 
              group_by(total) %>% 
              top_n(2, wt = -prop_sign_lfc_pos) %>% 
              dplyr::transmute(nn, label2=label))




p1<-DA_out %>% 
  group_by(total) %>% 
  dplyr::summarise(`LFC>0, p-adj<0.1` = sum(sign_lfc_pos)/sum(n_regions),
                   other = 1 - `LFC>0, p-adj<0.1`) %>% 
  pivot_longer(2:3) %>% 
  dplyr::mutate(name = factor(name, levels=c("other","LFC>0, p-adj<0.1"))) %>% 
  ggplot(aes(x=total, y=value, fill=name))+
  geom_col()+
  ylab("Proportion per PD")+
  scale_fill_manual(values=c(`LFC>0, p-adj<0.1` = "#003049",
                             other = "#669bbc"))+
  theme_bw()+
  theme(legend.title = element_blank(), 
        legend.position = "none",
        axis.text = element_text(color="black"))+
  xlab("PD")



p2<-DA_out %>% 
  dplyr::summarise(`LFC>0, p-adj<0.1` = sum(sign_lfc_pos)/sum(n_regions),
                   other = 1 - `LFC>0, p-adj<0.1`) %>% 
  pivot_longer(1:2) %>% 
  dplyr::mutate(name = factor(name, levels=c("other","LFC>0, p-adj<0.1"))) %>% 
  ggplot(aes(x=1, y=value, fill=name))+
  geom_col()+
  ylab("Overall proportion")+
  scale_fill_manual(values=c(`LFC>0, p-adj<0.1` = "#003049",
                             other = "#669bbc"))+
  theme_bw()+
  theme(legend.title = element_blank(),
        axis.text = element_text(color="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())


#cowplot::plot_grid(p1,p2, rel_widths = c(1,0.6), align = "hv", axis ="lrtb")
p1+p2 + plot_layout(widths = c(2, 0.3)) + plot_annotation(tag_levels = "A")
ggsave("figures/JAMM_DA_final.pdf", height = 4, width = 7)
# exactly 5% error rate