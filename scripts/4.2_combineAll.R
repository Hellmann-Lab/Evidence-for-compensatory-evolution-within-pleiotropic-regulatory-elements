library(broom)
# library(MAST)
library(viridis)
library(tidyverse)
library(cowplot)
library(data.table)
library(DESeq2)

tissueColors <- c( "#9E0142" ,"#D53E4F", "#F46D43", "#FDAE61", "#FFD92F" ,"#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
names(tissueColors) <-  c("adrenal_gland", "brain", "heart", "kidney", "large_intestine", "lung", "muscle", "stomach", "thymus")

specificityColors <- c( "#A3753B", "#CC9B57", "#E7CF97", "#F8EDD0", "#F7F7F7", "#D2EEEA", "#99D7CE", "#5DACA5", "#33847E")
names(specificityColors) <-  c(1:9)


#setwd("/data/share/htp/pleiotropy/paper_data/")
source("ATACseq/scripts/TFBS_analysis_functions_ATAC.R")


# PUT ALL INFO TOGETHER ####
# chromosomal positions
jamm_hg19<-readRDS("roadmap_DHS_summaries/region_summary/jamm_region_info_full.rds") %>%
  mutate(width=end-start)

# enhancer / promoter assignment
DHS_to_gene<-readRDS("CRE_to_Gene/DHS_to_gene.rds") %>%
  distinct(region_id, assignment, .keep_all = T) %>%
  mutate(assignment=gsub("enh","enhancer",assignment),
         assignment=gsub("prom","promoter",assignment))



# CGI EVO -----------------------------------------------------------------

jamm_inATAC<-readRDS("ATACseq/RDS/jamm_inATAC.rds") %>%
  filter(similar_width=="yes", Ns_hg38==0, Ns_macFas6==0) %>%
  dplyr::mutate(GC_distance = abs(GC_content.mac - GC_content.hg),
                CpG_obs_exp_distance = abs(CpG_obs_exp.mac - CpG_obs_exp.hg)) 


# SEQUENCE EVO ------------------------------------------------------------

phastCons_DHS<-readRDS("RDS/phastCons_DHS.rds") %>% mutate(total=as.factor(as.character(total)))
phyloP_DHS<-readRDS("RDS/phyloP_DHS.rds") %>% mutate(total=as.factor(as.character(total)))


# TFBS CONSERVATION -------------------------------------------------------
dd_10perc_exprTissues<-data.table::fread("ATACseq/cbust/top10perc_motifs_exprTissues_dt.csv") %>%
  calculate_per_region_TF_stats() %>%
  mutate(region_id = as.factor(region_id)) 


dd_10perc_exprNPC<-data.table::fread("ATACseq/cbust/top10perc_motifs_exprNPC_dt.csv") %>% 
  calculate_per_region_TF_stats() %>% 
  mutate(region_id = as.factor(region_id)) 


# still should remove ones with bad aln quality
TF_dists<-readRDS("ATACseq/cbust/RDS/TFBS_position/summarized_TFBS_distances_exprTissues_10perc.rds") %>%
  mutate(indel_perc=indel_bp/aln_length,
         mm_perc=mismatches/aln_length)

quantile(TF_dist$aln_length, probs=c(0.05,0.5,0.95))
quantile(TF_dist$indel_perc, probs=c(0.05,0.5,0.95))
quantile(TF_dist$mm_perc, probs=c(0.05,0.5,0.95))
quantile(TF_dist$nindel, probs=c(0.05,0.5,0.95))


# DEG, DA, OPENNESS -----------------------------------------------------------------

exprdds <- readRDS("RNAseq/RDS/dds_clean.rds")
DE<-lapply( c("NPC"), function(i){
  tmp <- colData(exprdds)
  tmp <- tmp[tmp$Differentiation == i,]
  tmpdds <- DESeqDataSetFromMatrix( colData = tmp,
                                    countData = DESeq2::counts(exprdds)[,rownames(tmp)],
                                    design = ~ Species)
  tmpdds <-DESeq(tmpdds)
  results(tmpdds) %>% data.frame %>% rownames_to_column("gene_id") %>% 
    mutate(celltype = i)
}) %>% bind_rows %>%
  dplyr::mutate(DEG=ifelse(padj<0.1,T,F))




atacDDS <- readRDS("ATACseq/RDS/atacDDS.rds")
DA<-lapply( c("NPC"), function(i){
  tmp <- colData(atacDDS)
  tmp <- tmp[tmp$cell_type == i,]
  tmpdds <- DESeqDataSetFromMatrix( colData = tmp,
                                    countData = DESeq2::counts(atacDDS)[,rownames(tmp)],
                                    design = ~ species)
  tmpdds <-DESeq(tmpdds)
  results(tmpdds) %>% data.frame %>% rownames_to_column("region_id") %>% 
    mutate(celltype = i) %>% 
    separate_wider_delim(region_id, ".", names = c(NA, "region_id")) %>% 
    mutate(region_id = as.numeric(region_id)) %>% 
    inner_join(DHS_to_gene)
}) %>% bind_rows %>%
  dplyr::mutate(DA=ifelse(padj<0.1,T,F))




jamm_forATAC_NPC_all<-readRDS("ATACseq/RDS/jammPeaks_vs_ATACPeaks_NPC_all.rds")
stringent_set_NPC<-readRDS("ATACseq/RDS/jammPeaks_vs_ATACPeaks_NPC_stringent.rds") 





# PUT IT ALL TOGETHER -----------------------------------------------------
# inner join for metrics that could be done on all sequences
# left join for NPC specific parameters: DA stuff


seq_vs_TFBS<-jamm_hg19 %>% 
  transmute(seqnames = chromosome, start, end, 
            region_id = as.factor(as.character(region_id)), CGI, total) %>%
  left_join(DHS_to_gene %>% 
              dplyr::transmute(region_id = as.factor(as.character(region_id)), 
                               assignment) %>% distinct()) %>%
  inner_join(jamm_inATAC %>% 
               filter(!is.nan(CpG_obs_exp_distance)) %>%
               dplyr::transmute(region_id = as.factor(as.character(region_id)), 
                                CpG_obs_exp_cons = round(1-CpG_obs_exp_distance, digits = 4),
                                GC_cons = round(1-GC_distance, digits = 4))) %>%
  inner_join(phastCons_DHS %>% ungroup() %>% 
               dplyr::transmute(region_id = as.factor(as.character(region_id)),
                                meanPhastCons = round(meanPhastCons, digits = 4)))  %>%
  inner_join(phyloP_DHS %>% ungroup() %>% 
               dplyr::transmute(region_id = as.factor(as.character(region_id)),
                                meanPhyloP = round(meanPhyloP, digits = 4))) %>%
  inner_join(dd_10perc_exprTissues %>% 
               dplyr::transmute(region_id = as.factor(as.character(region_id)), 
                                canb_cons.exprTissues = round(1-canb_corr, digits = 4))) %>%
  inner_join(dd_10perc_exprNPC %>% 
               dplyr::transmute(region_id = as.factor(as.character(region_id)), 
                                canb_cons.exprNPC = round(1-canb_corr, digits = 4))) %>%
  inner_join(TF_dists %>% 
               dplyr::transmute(region_id = as.factor(as.character(region_id)), 
                                meanBindingAgreement = round(mean_agreement, digits = 4))) %>%
  left_join(DA %>% dplyr::transmute(region_id = as.factor(as.character(region_id)), 
                                     DA, LFC.DA = log2FoldChange)) %>%
  left_join(jamm_forATAC_NPC_all %>% 
               dplyr::transmute(region_id = as.factor(as.character(region_id)),
                                openness_all = openness)) %>%
  left_join(stringent_set_NPC %>% 
               dplyr::transmute(region_id = as.factor(as.character(region_id)),
                                openness_stringent = openness)) 

saveRDS(seq_vs_TFBS, "ATACseq/cbust/RDS/seq_vs_TFBS_10perc.rds")






# Spider / radar plot -------------------------------------------------------------

library(fmsb)
seq_vs_TFBS<-readRDS("ATACseq/cbust/RDS/seq_vs_TFBS_10perc.rds")


chart<-seq_vs_TFBS %>%
  pivot_longer(cols=c(CpG_obs_exp_cons, meanPhyloP, canb_cons.exprTissues, meanBindingAgreement, LFC.DA)) %>%
  group_by(name, total) %>%
  dplyr::summarise(mean = mean(abs(value), na.rm=T)) %>%
  pivot_wider(id_cols = total, names_from = "name", values_from = "mean") %>%
  # Add downstream gene expression
  left_join(readRDS("RNAseq/RDS/DE_DA_table.rds") %>%
              filter(celltype=="NPC", !is.na(region_id)) %>%
              group_by(total) %>%
              dplyr::summarise(LFC.DE = mean(abs(log2FoldChange)))) %>%
  # add proportion of conserved peaks
  left_join(seq_vs_TFBS %>% 
              filter(openness_stringent %in% 
                       c("Human-only", "Macaque-only", "Always open")) %>%
              group_by(total) %>%
              dplyr::summarise(conserved_peaks = length(total[openness_stringent=="Always open"])/length(total))) %>%
  mutate(index = total+2,
         total= paste0("PD", total)) 


chart2<-chart %>%
  dplyr::select(-conserved_peaks) %>%
  add_row(total = as.character(1), 
          CpG_obs_exp_cons = max(chart$CpG_obs_exp_cons), 
          LFC.DA = min(chart$LFC.DA),
          canb_cons.exprTissues = max(chart$canb_cons.exprTissues), 
          meanBindingAgreement = max(chart$meanBindingAgreement), 
          meanPhyloP = max(chart$meanPhyloP),
          LFC.DE = min(chart$LFC.DE),
          index = 1) %>%
  add_row(total = as.character(2), 
          CpG_obs_exp_cons = min(chart$CpG_obs_exp_cons), 
          LFC.DA = max(chart$LFC.DA),
          canb_cons.exprTissues = min(chart$canb_cons.exprTissues), 
          meanBindingAgreement = min(chart$meanBindingAgreement), 
          meanPhyloP = min(chart$meanPhyloP),
          LFC.DE = max(chart$LFC.DE),
          index = 2) %>%
  arrange(index) %>%
  dplyr::select(-index) %>%
  relocate(total, CpG_obs_exp_cons, LFC.DE, LFC.DA, 
           canb_cons.exprTissues, meanBindingAgreement,
           meanPhyloP) %>%
  column_to_rownames("total")


chart2_relocated<-chart2 %>%
  relocate(LFC.DE, LFC.DA, CpG_obs_exp_cons, canb_cons.exprTissues, 
           meanBindingAgreement, meanPhyloP)

pdf("figures/combi_all.pdf", height=6.5, width=6.5)
par(xpd = TRUE, mfrow = c(1,1), mar = c(1, 1, 1, 1))
radarchart(chart2_relocated, pcol=specificityColors, plty=1,  cglwd=1.5,
           seg=1, plwd = 4.5, cglcol = "grey40", vlcex=1.5,
           vlabels=c("Downstream\nexpression", "Accessibility","CpG\nobs/exp",
                     "TFBS\nrepertoire", "TFBS\nposition", "Sequence"))
dev.off()
