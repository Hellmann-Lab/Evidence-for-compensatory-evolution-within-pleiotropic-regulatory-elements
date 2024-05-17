# best of chipseq vs our peaks
library(plyranges)
library(tidyverse)
library(JASPAR2020)
library(TFBSTools)
library(DESeq2)

#setwd("/data/share/htp/pleiotropy/paper_data/")

specificityColors <- c( "#A3753B", "#CC9B57", "#E7CF97", "#F8EDD0", "#F7F7F7", "#D2EEEA", "#99D7CE", "#5DACA5", "#33847E")
names(specificityColors) <-  c(1:9)

speciesCols<-c("human"="#8B8BA7","macaque"="#D1D1DC")

source("ATACseq/scripts/TFBS_analysis_functions_ATAC.R")
source("ATACseq/scripts/TFBS_plotting_functions_ATAC.R")




# LOAD CRE INFO -----------------------------------------------------------

# need the coords in hg38 of the region_ids that have ORTHOLOGUES
jamm_inATAC<-readRDS("ATACseq/cbust/RDS/jamm_inATAC.rds") %>%
  filter(similar_width=="yes", Ns_macFas6==0, Ns_hg38==0)

# coordinate info hg19 (phyloP etc relevant)
jamm_hg19<-readRDS("roadmap_DHS_summaries/region_summary/jamm_region_info_full.rds") %>% 
  transmute(seqnames = chromosome, start, end, region_id, total)

# get the chromosomes and files where these region_ids can be found (cbust-relevant)
locfile<-setNames(read.table("ATACseq/cbust/RDS/TFBS_position/selected_ids.txt"), c("region_id","human_file", "macaque_file"))


# CRE characterization
seq_vs_TFBS_NPCs<-readRDS("ATACseq/cbust/RDS/seq_vs_TFBS_10perc.rds") %>%
  filter(openness_stringent %in% c("Human-only", "Macaque-only", "Always open"))




# EXPRESSION & CRE-GENE ASSOCIATION ---------------------------------------

# expressed TFs in npcs
expressedTFs<-readRDS("ATACseq/cbust/expressedTFs.rds") %>%
  dplyr::select(-IC_rank, -IC, -class) %>% dplyr::rename(TF_gene_id = gene_id)

# get annotation gtf: i'm only using symbols and gene_ids 
gene_to_symbol <- rtracklayer::import("gtf/gencode.v19.annotation.gtf.gz") %>% 
  S4Vectors::mcols() %>% data.frame() %>% 
  dplyr::select(gene_id, gene_name) %>% distinct() %>% 
  mutate(gene_id = sub("(.*)\\..*", "\\1", gene_id)) 


# expression
exprdds <- readRDS("expression_conservation/RDS/dds_clean.rds")
DE<-lapply( c("NPC"), function(i){
  tmp <- colData(exprdds)
  tmp <- tmp[tmp$Differentiation == i,]
  tmpdds <- DESeqDataSetFromMatrix( colData = tmp,
                                    countData = DESeq2::counts(exprdds)[,rownames(tmp)],
                                    design = ~ Species)
  tmpdds <-DESeq(tmpdds)
  results(tmpdds) %>% data.frame %>% rownames_to_column("gene_id") %>% mutate(celltype = i)
}) %>% bind_rows %>%
  dplyr::mutate(DEG=ifelse(padj<0.1,T,F)) %>%  
  left_join(gene_to_symbol) 
  


# identify the subset of CREs that have been associated with genes
DHS_genes <- readRDS("CRE_to_Gene/DHS_to_gene.rds") %>%
  dplyr::select(-tissue) %>%
  distinct(region_id, gene_id, assignment, total, distance) %>%
  group_by(gene_id) %>%  
  mutate(promoter_number_per_gene = length(gene_id)) %>% 
  ungroup() %>%
  # join with DE
  inner_join(DE %>% dplyr::select(gene_id, DEG, gene_name)) %>%
  dplyr::filter(region_id %in% jamm_inATAC$region_id)


# load TF filter for being among the top 10% for that CRE
dd_10perc<-data.table::fread("ATACseq/cbust/top10perc_motifs_exprNPC_dt.csv")



# ADDITIONAL EXTERNAL INFORMATION -----------------------------------------

# metadata of GTRD database with ids and celltypes
metdat<-read.table("ATACseq/chipseq_gtrd/metadata_chipseq.txt", 
                   header = T, fill = T, sep="\t") %>%
  dplyr::mutate(id=as.character(id))

# topGO neurogenesis category
neuro<-read.table("ATACseq/chipseq_gtrd/neurogenesis_GO", fill = TRUE)

length(unique(neuro$V3))




# RANKING AND SELECTING CREs ----------------------------------------------

seq_vs_TFBS_rank<-seq_vs_TFBS_NPCs %>%
  arrange(-CpG_obs_exp_cons) %>%
  mutate(rankCpG=1:nrow(.)) %>%
  arrange(-canb_cons.exprTissues) %>%
  mutate(rankCanberra=1:nrow(.)) %>%
  arrange(-meanPhastCons) %>%
  mutate(rankPhastCons=1:nrow(.)) %>%
  arrange(-meanBindingAgreement) %>%
  mutate(rankBinding=1:nrow(.)) %>%
  # lowest LFC = most cons
  arrange(abs(LFC.DA)) %>%
  mutate(rankDA=1:nrow(.),
         region_id = as.character(region_id))


divCREs<-DHS_genes %>% 
  #filter(gene_name %in% neuro$V3) %>%
  filter(DEG==FALSE, promoter_number_per_gene<5, assignment=="prom") %>%
  inner_join(seq_vs_TFBS_rank %>% 
               dplyr::transmute(region_id = as.double(region_id), rankCanberra, rankPhastCons, rankBinding)) %>%
  rowwise() %>%
  mutate(dist_ranks = abs(rankCanberra-rankPhastCons),
         distTF_ranks = abs(rankCanberra-rankBinding)) %>%
  filter(rankCanberra < rankPhastCons, rankCanberra < rankBinding, total==9) %>%
  filter(rankCanberra < 14000, distTF_ranks>10000, dist_ranks>10000) %>%
  arrange(-distTF_ranks)





# POTENTIAL EXAMPLES ------------------------------------------------------

# ATXN3

lab = "ATXN3"
r<-divCREs$region_id[divCREs$gene_name==lab]
dt_tfs<-dd_10perc[region_id == r,]

# also when excluding neuro_topGO_filter, the result is ok
genetab<-get_chip(name=lab, dnase = DHS_genes, metadata_filter = "neur", old = F,
                  neuro_topGO_filter = neuro$V3) 

if (sum(grepl(genetab$motifs[genetab$region_id==r], dt_tfs$cluster_id_or_motif_name))>0){
  
  ext_grey_alpha<-plot_targeted_forFig(bw = "bw/primates.phastCons46way.bw", 
                                       regid = r, 
                                       coord_file_4phastCons = jamm_hg19,
                                       locfile_4cbust = locfile,
                                       padding_bp = 50, 
                                       conc_TFs = genetab$motifs[genetab$region_id==r],
                                       motif_filter = dt_tfs$cluster_id_or_motif_name,
                                       TF_score_cutoff = 3,
                                       cluster_score_cutoff = 0,
                                       highlight1 = "#D4D8D8",
                                       highlight1_alpha = 0.3, 
                                       plot_labels = c("F","G","H"))
  
  ranks<-plot_ranks(r=r, lab=lab, plot_labels = c("A","B","C","D","E"))
  
  mycn<- ggdraw() + 
    cowplot::draw_image(magick::image_read_pdf("figures/logo_MYCN.pdf"))+
    theme(plot.margin = unit(c(0, 0,0,0), "cm"))
  pou3f2<- ggdraw() + 
    cowplot::draw_image(magick::image_read_pdf("figures/logo_POU3F2.pdf"))+
    theme(plot.margin = unit(c(0, 0,0,0), "cm"))

  
  motif_p<-plot_grid(NULL,mycn, pou3f2, labels = c("","I", "J"), scale=0.9, rel_widths=c(0.1,1,1), ncol=3, label_size=12)
  right_p<-plot_grid(ext_grey_alpha, motif_p, NULL, ncol=1, rel_heights = c(1,0.25,0.03))
  
  ggsave(paste0("figures/new_fig5/",lab,"_allFilters_new.pdf"), plot_grid(ranks, right_p, rel_widths = c(0.78,0.92)), height = 180*0.97, width = 162.5*0.92, units = "mm")
}


pfm_mycn <- getMatrixByID(JASPAR2020, ID = "MA0104.4")
pdf("figures/logos/logo_MYCN.pdf", height = 3.2, width = 5)
seqLogo(toICM(pfm_mycn))
dev.off()

pfm_pou3f2 <- getMatrixByID(JASPAR2020, ID = "MA1114.1")
pdf("figures/logos/logo_POU3F2.pdf", height = 3.2, width = 5)
seqLogo(toICM(pfm_pou3f2))
dev.off()


