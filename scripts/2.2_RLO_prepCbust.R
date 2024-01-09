
library(Biostrings)
library(stringr)
library(tidyverse)
library(plyranges)
library(RMySQL)
library(cowplot)
library(GenomicRanges)
library(DESeq2)
library(pheatmap)
library(UpSetR)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mfascicularis.NCBI.6.0)
library(BSgenome.Hsapiens.NCBI.GRCh38)

# liftover hg19 dhs to hg38
# liftover hg38 to cyno6
# make a bed file for each species
# count reads falling in these


#setwd("/data/share/htp/pleiotropy/paper_data/")
source("ATACseq/scripts/functions.R")
source("ATACseq/scripts/helper_functions.R")
workdir<-"/data/share/htp/pleiotropy/paper_data/"


jamm_region_identifiers<-readRDS("roadmap_DHS_summaries/region_summary/jamm_region_info_full.rds") %>%
  transmute(seqnames, start,end, region_id) %>%
  dplyr::filter(start < end & end - start<=5000)


# RECIPROCAL LIFTOVER -----------------------------------------------------

jamm_hg38_noExt<-translate_jamm(chain_file = "macFas6_liftOver_chains/hg19ToHg38.over.chain",
                                coordinate_file = as_granges(jamm_region_identifiers), 
                                extend = 20,
                                reverse_chain_file = "macFas6_liftOver_chains/hg38ToHg19.over.chain")

jamm_macfas6_noExt<-translate_jamm(chain_file = "macFas6_liftOver_chains/hg38ToMacFas6.over.chain", 
                                   coordinate_file = as_granges(jamm_hg38_noExt), 
                                   extend = 20,
                                   reverse_chain_file = "macFas6_liftOver_chains/macFas6ToHg38.over.chain")


#save both as bed files
write.table(jamm_hg38_noExt %>% 
              transmute(seqnames = gsub("chr","",seqnames), start, end, region_id),
            "ATACseq/bed/hg38_DHS.bed", col.names = F, row.names = F, quote = F, sep="\t")

write.table(jamm_macfas6_noExt %>% 
              transmute(seqnames = gsub("chr","",seqnames), start, end, region_id),
            "ATACseq/bed/macFas6_DHS.bed", col.names = F, row.names = F, quote = F, sep="\t")


# run the counting in these regions
system("ATACseq/scripts/count_table.sh")





# GET SEQUENCES -----------------------------------------------------------


# first, prepare data for cbust: need to save the sequences and the PWMs
# SAVE SEQUENCES
#make a MACAQUE granges object
coords_mac_gr<-jamm_macfas6_noExt %>%
  filter(region_id %in% jamm_hg38_noExt$region_id) %>%
  filter(!grepl("JAA|MU|MT", seqnames)) %>% 
  mutate(seqnames=gsub("chr","",seqnames)) %>%
  filter(start>300) %>%
  plyranges::as_granges()

writeFasta4cbust(gr = coords_mac_gr,
                 genome = BSgenome.Mfascicularis.NCBI.6.0,
                 fasta.file = paste0(workdir,"ATACseq/cbust/fastas/macFas6/macFas6"),
                 padding = 500, 
                 max.seq = 5000,
                 id.col = "region_id")


writeFasta4cbust(gr = coords_mac_gr,
                 genome = BSgenome.Mfascicularis.NCBI.6.0,
                 fasta.file = paste0(workdir,"ATACseq/cbust/fastas/macFas6NoExt/macFas6NoExt"),
                 padding = 0, 
                 max.seq = 5000,
                 id.col = "region_id")


#make a HUMAN granges object, only keep common region_ids
coords_hum_gr<-jamm_hg38_noExt %>% 
  filter(region_id %in% jamm_macfas6_noExt$region_id) %>%
  filter(!grepl("KI|GL|MT", seqnames)) %>% 
  mutate(seqnames=gsub("chr","",seqnames)) %>%
  plyranges::as_granges()

writeFasta4cbust(gr = coords_hum_gr,
                 genome = BSgenome.Hsapiens.NCBI.GRCh38,
                 fasta.file = paste0(workdir,"ATACseq/cbust/fastas/hg38/hg38"),
                 padding = 500, 
                 max.seq = 5000,
                 id.col = "region_id")


writeFasta4cbust(gr = coords_hum_gr,
                 genome = BSgenome.Hsapiens.NCBI.GRCh38,
                 fasta.file = paste0(workdir,"ATACseq/cbust/fastas/hg38NoExt/hg38NoExt"),
                 padding = 0, 
                 max.seq = 5000,
                 id.col = "region_id")



# SANITY CHECKS -----------------------------------------------------------


# quantify width differences and Ns  
upper=1.2
lower=0.8

width_comp<-inner_join(jamm_region_identifiers %>% transmute(region_id, width.hg19 = end-start), 
                   jamm_hg38_noExt %>% transmute(region_id, width.hg38 = end-start)) %>%
  inner_join(jamm_macfas6_noExt %>% transmute(region_id, width.macfas6 = end-start)) %>%
  mutate(hg38 = width.hg38/width.hg19,
         macfas6 = width.macfas6/width.hg19) %>%
  mutate(similar_width=ifelse(hg38>lower & hg38<upper & macfas6>lower & macfas6<upper, "yes","no"))

macFas_Ns<-get_Ns("macFas6")
hg38_Ns<-get_Ns("hg38")





# GC AND CGI --------------------------------------------------------------


# MACAQUE
seq_mac<-data.frame(region_id = coords_mac_gr$region_id,
                    seq = getSeq(BSgenome.Mfascicularis.NCBI.6.0, coords_mac_gr)) %>%
  rowwise() %>%
  mutate(GC = calc_gc(DNAString(seq)),
         CpG_obs_exp = calc_cpg_oe(DNAString(seq))) %>%
  dplyr::transmute(region_id, GC_content.mac = GC, CpG_obs_exp.mac = CpG_obs_exp)


# HUMAN
seq_hum<-data.frame(region_id = coords_hum_gr$region_id,
                    seq = getSeq(BSgenome.Hsapiens.NCBI.GRCh38, coords_hum_gr)) %>%
  rowwise() %>%
  mutate(GC = calc_gc(DNAString(seq)),
         CpG_obs_exp = calc_cpg_oe(DNAString(seq))) %>%
  dplyr::transmute(region_id, GC_content.hg = GC, CpG_obs_exp.hg = CpG_obs_exp)





# COMBINE AND SAVE --------------------------------------------------------


jamm_inATAC<-inner_join(as_tibble(coords_hum_gr),
                        as_tibble(coords_mac_gr),
                        by="region_id",
                        suffix = c("", ".macFas6")) %>%
  mutate(region_id=as.factor(region_id)) %>%
  left_join(width_comp %>% dplyr::transmute(region_id=as.factor(region_id), similar_width)) %>%
  left_join(macFas_Ns %>% dplyr::transmute(region_id, Ns_macFas6 = N)) %>%
  left_join(hg38_Ns %>% dplyr::transmute(region_id, Ns_hg38 = N)) %>%
  left_join(seq_mac %>% ungroup() %>% mutate(region_id = as.factor(region_id))) %>%
  left_join(seq_hum %>% ungroup() %>% mutate(region_id = as.factor(region_id)))

saveRDS(jamm_inATAC, "ATACseq/RDS/jamm_inATAC.rds")





# GET PWMS ----------------------------------------------------------------

# get annotation gtf 
gtf <- rtracklayer::import("gtf/gencode.v19.annotation.gtf.gz")

#get the jaspar2020 motif PFM list
opts <- list()
opts[["collection"]] <- "CORE"
opts[["matrixtype"]] <- "PFM"
opts[["tax_group"]] <- "vertebrates"

PFMatrixList <- getMatrixSet(JASPAR2020, opts)

jaspar_icm <- toICM(PFMatrixList, pseudocounts=0.8, schneider=TRUE)
jaspar_icm_mat <- Matrix(jaspar_icm)
names(jaspar_icm_mat) <- paste0(ID(jaspar_icm),"_",name(jaspar_icm))
jaspar_ic <- jaspar_icm_mat %>% sapply(sum)

jaspar_ic_df <- jaspar_ic %>%
  as.data.frame() %>%
  rownames_to_column(var = "motif") %>%
  dplyr::rename(IC = ".") %>%
  mutate(motif_ID = word(motif,1,1,"_"),
         motif_name = word(motif,2,2,"_")) %>%
  relocate(motif_ID,motif_name) %>%
  dplyr::select(!motif)



# expressed TFs in iPSC-NPC dataset ####

cnts<-readRDS("RNAseq/RDS/cnt_clean.rds")

gene_to_symbol <- gtf %>% S4Vectors::mcols() %>% data.frame() %>% 
  dplyr::select(gene_id, gene_name) %>% distinct() %>% 
  mutate(gene_id = sub("(.*)\\..*", "\\1", gene_id)) %>%
  filter(gene_id %in% rownames(cnts))

# motif ID to gene name info
jaspar_ic_df_diff <-jaspar_ic_df %>%
  mutate(TF = gsub("[(].*","",motif_name),
         TF = toupper(TF)) %>%
  separate_rows(TF,sep='::') %>%
  group_by(motif_ID) %>%
  mutate(nTFs=length(motif_ID)) %>%
  #filter expressed TFs
  inner_join(gene_to_symbol, by=c("TF"="gene_name")) %>%
  group_by(motif_ID) %>%
  mutate(nTFs_expressed=length(motif_ID))

length(unique(jaspar_ic_df_diff$motif_ID)) # 521 motifs
length(unique(jaspar_ic_df_diff$TF)) # 446 TFs

table(jaspar_ic_df_diff$nTFs==jaspar_ic_df_diff$nTFs_expressed)

# ok basically nearly all TFs in all dimers are detected (besides 18); good enough, keep.

# filter PWMs
PFMatrixList_expr<-PFMatrixList[names(PFMatrixList) %in% unique(jaspar_ic_df_diff$motif_ID)]
savePFM(PFMatrixList_obj = PFMatrixList_expr, out.file = "ATACseq/cbust/PWMs/exprIpscNpc")

# also save the TF family and class
class<-bind_rows(lapply(PFMatrixList, function(x){data.frame(class=x@matrixClass)}), .id="motif_ID")


jaspar_ic_df_diff %>% 
  filter(motif_ID %in% names(PFMatrixList)) %>%
  left_join(class) %>%
  dplyr::select(-nTFs, -nTFs_expressed) %>%
  saveRDS("ATACseq/cbust/expressedTFs.rds")

family<-bind_rows(lapply(PFMatrixList, function(x){data.frame(family=x@tags$family)}), .id="motif_ID")
saveRDS(family, "ATACseq/cbust/TF_family.rds")







# expressed in NPCs ####

exprdds <- readRDS("RNAseq/RDS/dds_clean.rds")
tmp <- colData(exprdds)
tmp <- tmp[tmp$Differentiation == "NPC",]
countData = DESeq2::counts(exprdds)[,rownames(tmp)]

freq_umicounts <- apply(countData > 0, 1, mean, na.rm = T)
expressed_genes <- freq_umicounts > 0.25
table(expressed_genes)
umi_counts_filt <- countData[expressed_genes,]

exprTFs_ipsc_npc<-readRDS("ATACseq/cbust/expressedTFs.rds")

exprTFs_ipsc_npc[!exprTFs_ipsc_npc$gene_id %in% rownames(umi_counts_filt),]
# nothing to drop - same set as when taking both ipscs and npcs





# expressed TFs in human tissues ####

summarized_expression_all<-readRDS("roadmap_expression_summaries/summarized_expression_allTissues.rds")

gene_to_symbol_tissues <- gtf %>% S4Vectors::mcols() %>% data.frame() %>% 
  dplyr::select(gene_id, gene_name) %>% distinct() %>% 
  mutate(gene_id = sub("(.*)\\..*", "\\1", gene_id)) %>%
  filter(gene_id %in% summarized_expression_all$gene_id)

# motif ID to gene name info
jaspar_ic_df_tissues <-jaspar_ic_df %>%
  mutate(TF = gsub("[(].*","",motif_name),
         TF = toupper(TF)) %>%
  separate_rows(TF,sep='::') %>%
  group_by(motif_ID) %>%
  mutate(nTFs=length(motif_ID)) %>%
  #filter expressed TFs
  inner_join(gene_to_symbol_tissues, by=c("TF"="gene_name")) %>%
  group_by(motif_ID) %>%
  mutate(nTFs_expressed=length(motif_ID))

length(unique(jaspar_ic_df_tissues$motif_ID)) # 643 motifs
length(unique(jaspar_ic_df_tissues$TF)) # 560 TFs

table(jaspar_ic_df_tissues$nTFs==jaspar_ic_df_tissues$nTFs_expressed)

# ok basically nearly all TFs in all dimers are detected (9); good enough, keep.

# filter PWMs
PFMatrixList_Tissueexpr<-PFMatrixList[names(PFMatrixList) %in% unique(jaspar_ic_df_tissues$motif_ID)]
savePFM(PFMatrixList_obj = PFMatrixList_Tissueexpr, out.file = "ATACseq/cbust/PWMs/exprTissues")

# also save the TF family and class
class<-bind_rows(lapply(PFMatrixList, function(x){data.frame(class=x@matrixClass)}), .id="motif_ID")


jaspar_ic_df_tissues %>% 
  filter(motif_ID %in% names(PFMatrixList)) %>%
  left_join(class) %>%
  dplyr::select(-nTFs, -nTFs_expressed) %>% 
  saveRDS("ATACseq/cbust/expressedTFs_inTissues.rds")

