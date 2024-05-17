
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
library(patchwork)



setwd("/data/share/htp/pleiotropy/paper_data/")
source("ATACseq/scripts/functions.R")
source("ATACseq/scripts/helper_functions.R")


regionColors <- c("#9CA578","#3288BD","#E37063")
names(regionColors) <-  c("Enhancer", "Promoter","Expression")

specificityColors <- c( "#A3753B", "#CC9B57", "#E7CF97", "#F8EDD0", "#F7F7F7", "#D2EEEA", "#99D7CE", "#5DACA5", "#33847E")
names(specificityColors) <-  c(1:9)




# get the original jamm peaks
jamm_gr<-readRDS("roadmap_DHS_summaries/region_summary/jamm_region_info_full.rds") %>%
  dplyr::rename(seqnames=chromosome) %>%
  as_granges()

# keep only similar length and N-free sequences 
jamm_inATAC<-readRDS("ATACseq/RDS/jamm_inATAC.rds") %>% 
  filter(region_id %in% jamm_gr$region_id, similar_width=="yes", Ns_hg38==0, Ns_macFas6==0) 




# Human peaks vs DHS (download these genrich peaks from E-MTAB-13373) ----------------------------------


path_hum<-"/data/share/htp/ATACseq/Novogene_ATAC/02_Preprocessing/bwa_mapping/bwa_cleaning/05_peaks/genrich/"
hg38_peaks_NPC<-read_narrowpeaks(paste0(path_hum,"ATAC_b1_sample05_hg38.bam.narrowPeak")) %>%
  mutate(peak_id.hum = paste0(name,".hum"))

#first find overlaps with DHS
hg38_jamm<-as_granges(jamm_inATAC %>% transmute(seqnames, start, end, region_id)) 


OL_jamm_ATAC_hg38<-calc_OL(hg38_jamm, 
                           hg38_peaks_NPC, 
                           filt = T, 
                           rel_width_OL_cutoff = 0.1, 
                           id1_col = "region_id", 
                           id2_col = "peak_id.hum") %>%
  dplyr::rename(width_jamm.hum = width1,
                width_ATAC.hum = width2,
                fracOverlap1.hum = fracOverlap1,
                fracOverlap2.hum = fracOverlap2,
                OneToOne.hum = OneToOne) %>%
  dplyr::select(-n_Ind1ToInd2,-n_Ind2ToInd1)

# a good proportion (72%) of human NPC peaks are present in DHS
length(hg38_peaks_NPC)
length(unique(OL_jamm_ATAC_hg38$peak_id.hum))/length(unique(hg38_peaks_NPC$peak_id.hum))
length(unique(OL_jamm_ATAC_hg38$region_id))/length(unique(hg38_jamm$region_id))






# Macaque peaks vs DHS ----------------------------------------------------

path_mac<-"/data/share/htp/ATACseq/Novogene_ATAC/02_Preprocessing/bwa_mapping/bwa_cleaning/05_peaks/genrich/"
mac6_peaks_NPC<-read_narrowpeaks(paste0(path_mac,"ATAC_b1_sample07_macFas6.bam.narrowPeak")) %>%
  mutate(peak_id.mac = paste0(name,".mac"))

#first find overlaps with DHS
mac6_jamm<-as_granges(jamm_inATAC %>% transmute(seqnames = seqnames.macFas6, start = start.macFas6, end = end.macFas6, region_id))

# this function filters overlaps with at least 10% and annotates whether 
OL_jamm_ATAC_mac6<-calc_OL(mac6_jamm, 
                           mac6_peaks_NPC, 
                           filt = T, 
                           rel_width_OL_cutoff = 0.1, 
                           id1_col = "region_id", 
                           id2_col = "peak_id.mac") %>%
  dplyr::rename(width_jamm.mac = width1,
                width_ATAC.mac = width2,
                fracOverlap1.mac = fracOverlap1,
                fracOverlap2.mac = fracOverlap2,
                OneToOne.mac = OneToOne) %>%
  dplyr::select(-n_Ind1ToInd2,-n_Ind2ToInd1)

# a good proportion (68%) of macaque NPC peaks are present in human DHS
length(mac6_peaks_NPC)
length(unique(OL_jamm_ATAC_mac6$peak_id.mac))/length(unique(mac6_peaks_NPC$peak_id.mac))
length(unique(OL_jamm_ATAC_mac6$region_id))/length(unique(mac6_jamm$region_id))





# Multimatches ------------------------------------------------------------

# We want to exclude cases with multimatches and require that in a case of 1-to-1 overlap, also the human and macaque NPC peaks need to overlap to at least 10%

multimatches<-bind_rows(OL_jamm_ATAC_hg38 %>% filter(!OneToOne.hum) %>% mutate(species="human"),
                        OL_jamm_ATAC_mac6 %>% filter(!OneToOne.mac) %>% mutate(species="macaque"))





# DHS peaks vs ATAC peaks -------------------------------------------------

OL_jamm_ATAC_hg38_loose<-calc_OL(hg38_jamm, 
                           hg38_peaks_NPC, 
                           filt = T, 
                           rel_width_OL_cutoff = 0, 
                           id1_col = "region_id", 
                           id2_col = "peak_id.hum") 

OL_jamm_ATAC_mac6_loose<-calc_OL(mac6_jamm, 
                           mac6_peaks_NPC, 
                           filt = T, 
                           rel_width_OL_cutoff = 0, 
                           id1_col = "region_id", 
                           id2_col = "peak_id.mac") 

# anly annotate whether the DHS coordinates overlap any NPC peak
jamm_forATAC_NPC_all<-jamm_inATAC %>%
  #loose 550 genes due to this step
  inner_join(readRDS("CRE_to_Gene/DHS_to_gene.rds") %>% 
               distinct(region_id = as.factor(region_id), assignment, total, CGI)) %>%
  dplyr::select(region_id, assignment, total, CGI) %>% 
  dplyr::mutate(OL.hum = ifelse(region_id %in% OL_jamm_ATAC_hg38$region_id, T, F),
                OL.mac = ifelse(region_id %in% OL_jamm_ATAC_mac6$region_id, T, F),
                openness = case_when(OL.hum & OL.mac ~ "Always open",
                                     !OL.hum & !OL.mac ~ "Not open",
                                     !OL.hum & OL.mac ~ "Macaque-only",
                                     OL.hum & !OL.mac ~ "Human-only"))

saveRDS(jamm_forATAC_NPC_all, "ATACseq/RDS/jammPeaks_vs_ATACPeaks_NPC_all.rds")



# EXCLUDE MULTIMATCHES
jamm_forATAC_NPC<-jamm_inATAC %>%
  #loose 550 genes due to this step
  inner_join(readRDS("CRE_to_Gene/DHS_to_gene.rds") %>% 
               distinct(region_id = as.factor(region_id), assignment, total, CGI)) %>%
  #inner_join(readRDS("../expression_conservation/Data/DGE_PD/DHS_to_gene_expr.rds") %>%
  #             distinct(region_id = as.factor(region_id), assignment, total, CGI)) %>% 
  dplyr::select(region_id, assignment, total, CGI) %>% 
  filter(!region_id %in% multimatches$region_id) %>%
  left_join(OL_jamm_ATAC_hg38 %>% 
              filter(!region_id %in% multimatches$region_id, !peak_id.hum %in% multimatches$peak_id.hum), by="region_id") %>%
  left_join(OL_jamm_ATAC_mac6 %>% 
              filter(!region_id %in% multimatches$region_id, !peak_id.mac %in% multimatches$peak_id.mac), by="region_id") %>%
  mutate(openness = case_when(OneToOne.hum & OneToOne.mac ~ "Always open",
                              is.na(OneToOne.hum) & is.na(OneToOne.mac) ~ "Not open",
                              is.na(OneToOne.hum) & OneToOne.mac ~ "Macaque-only",
                              OneToOne.hum & is.na(OneToOne.mac) ~ "Human-only"))

# visualize overlaps
jamm_forATAC_NPC %>% pivot_longer(contains("fracOverlap")) %>% ggplot(aes(x=value))+geom_histogram()+facet_grid(.~name) 
# 25% have 100% overlap






# Peak - peak OLs ---------------------------------------------------------

# how many of these overlap between human and mac peak? Liftover human peaks, that had any overlaps with DHS, to macaque genome 

coord_hummatches<-hg38_peaks_NPC %>% as_tibble() %>% 
  mutate(region_id = peak_id.hum, seqnames = paste0("chr",seqnames)) %>% 
  filter(peak_id.hum %in% jamm_forATAC_NPC$peak_id.hum[jamm_forATAC_NPC$openness=="Always open"]) %>% 
  as_granges()

humNPC_onMacfas6<-translate_jamm(chain_file = "liftOvers/hg38ToMacFas6.over.chain", 
                                 coordinate_file = coord_hummatches, 
                                 extend = 50,
                                 reverse_chain_file = "liftOvers/macFas6ToHg38.over.chain") %>%
  as_tibble() %>% mutate(seqnames = gsub("chr", "", seqnames), peak_id.hum=region_id)


# these could not be liftovered (696)
#noLO<-coord_hummatches %>% filter(!peak_id.hum %in% humNPC_onMacfas6$peak_id.hum)

# now let's see which have at least 10% OL rel to either species peak
macNPC_onMacfas6<-mac6_peaks_NPC %>% 
  filter(peak_id.mac %in% jamm_forATAC_NPC$peak_id.mac[jamm_forATAC_NPC$openness=="Always open"]) %>% 
  as_tibble()

OL_hum_mac_inMac6<-calc_OL(humNPC_onMacfas6, 
                           macNPC_onMacfas6, 
                           filt = T, 
                           rel_width_OL_cutoff = 0.1, 
                           id1_col = "peak_id.hum", 
                           id2_col = "peak_id.mac") 


# these ones might just actually not overlap either due to noLO or noOL
noOL<-jamm_forATAC_NPC %>% 
  filter(openness=="Always open") %>% 
  anti_join(OL_hum_mac_inMac6, by = join_by(peak_id.hum, peak_id.mac)) 

stringent_set<-jamm_forATAC_NPC %>% anti_join(noOL, by = join_by(peak_id.hum, peak_id.mac))

saveRDS(stringent_set, "ATACseq/RDS/jammPeaks_vs_ATACPeaks_NPC_stringent.rds")


#just a checkup: indeed, also no LO is excluded as expected
#stringent_set %>% inner_join(as_tibble(noLO), by="peak_id.hum")






# DA analysis -------------------------------------------------------------

# we only care about the npcs currently
jamm_forATAC_NPC_all<-readRDS("ATACseq/RDS/jammPeaks_vs_ATACPeaks_NPC_all.rds") %>%
  filter(openness!="Not open")


# get the count matrix
# first 4 columns are chr, start, end, region_id
hum_counts<-setNames(read.table("ATACseq/count_tables/hg38_counts.bed") %>% 
                       as_tibble(),
                     c("seqnames", "start", "end", "region_id", 
                       "iPSC_01.hg", "iPSC_02.hg", 
                       "NPC_05.hg", "NPC_06.hg")) 

mac_counts<-setNames(read.table("ATACseq/count_tables/macFas6_counts.bed") %>%
                       as_tibble(), 
                     c("seqnames", "start", "end", "region_id", 
                       "iPSC_03.mac", "iPSC_04.mac", 
                       "NPC_07.mac", "NPC_08.mac"))

# inner join is important, otherwise we have included region ids that did not have recipr LO to macaque
count_matrix<-inner_join(mac_counts %>% dplyr::select(-seqnames,-start,-end),
                         hum_counts %>% dplyr::select(-seqnames,-start,-end)) %>%
  filter(region_id %in% jamm_forATAC_NPC_all$region_id) %>%
  mutate(region_id=paste0("region_id.",region_id)) %>%
  column_to_rownames("region_id")

saveRDS(count_matrix,"ATACseq/count_tables/count_matrix_4dds.rds")


metaData<-read.table("ATACseq/sampleinfo", header = T) 
metaData<-metaData %>%
  mutate(short_number = formatC(short_number, width = 2, format = "d", flag = "0"),
         sample=paste0(cell_type, "_", short_number, ifelse(species=="human", ".hg", ".mac")))
rownames(metaData)<-metaData$sample



# match the order of annotation and count mat
metaData<-metaData[match(colnames(count_matrix), rownames(metaData)),]
count_matrix<-count_matrix[,match(rownames(metaData),colnames(count_matrix))]

# do DESeq2 ####
atacDDS <- DESeqDataSetFromMatrix(count_matrix, metaData, 
                                  design = ~ species + cell_type + species*cell_type)
atacDDS <- DESeq(atacDDS)
saveRDS(atacDDS,"ATACseq/RDS/atacDDS.rds")

