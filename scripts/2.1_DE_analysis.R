library(tidyverse)
library(DESeq2)
library(cowplot)
library(vsn)
library(plyranges)
library(tidyverse)

#setwd("/data/share/htp/pleiotropy/paper_data/")
source("scripts/helper_functions.R")
# to see the zUMIs run, check RNAseq folder

# COLLECT THE DATA --------------------------------------------------------


# Build info data frame (bar code, sample, species) ####
bc_sample_species<-read.csv("RNAseq/sample_annotation.txt")

# Get zUMI expression data ####
## hg38
zumihg38 <- readRDS("RNAseq/zumis/hg38/zUMIs_output/expression/npc_diff_hg38.dgecounts.rds")
exon_mat_hg38<- as.matrix(zumihg38$umicount$exon$all) 

## macFas6
zumimacFas6 <- readRDS("RNAseq/zumis/macFas6/zUMIs_output/expression/npc_diff_macFas6.dgecounts.rds")
exon_mat_macFas6 <- as.matrix(zumimacFas6$umicount$exon$all)


# Curate rows & cols of info df and exon_mat ####
# Now we need to make sure to only keep bc_sample_species rows that are also present in our exon_mat_hg38/macFas6,
# the row names of bc_sample_species must be identical and in the same order than the column names of the count table.

## Check if col names in both exon_mat is the same
table(colnames(exon_mat_hg38) == colnames(exon_mat_macFas6))

## Subset the bc_sample_species by species ####
hSample <- bc_sample_species %>% filter(Species == "human")
macSample <- bc_sample_species %>% filter(Species == "macFas") 


## Subset count matrices by species ####
#select only human counts & make it a df to be able to inner_join
exon_mat_hg38 <- exon_mat_hg38[,hSample$barcode] %>%
  as.data.frame() %>% 
  rownames_to_column("exon")

#select only rheMac10 counts & make it a df to be able to inner_join
exon_mat_macFas6 <- exon_mat_macFas6[,macSample$barcode] %>%  
  as.data.frame() %>% 
  rownames_to_column("exon")

## Create info df for Samples Of Interest colData ####
colDataSOI <-bind_rows(hSample, macSample) %>%
  dplyr::mutate(Time=as.factor(Time),
                Species=as.factor(Species))

rownames(colDataSOI) <- colDataSOI$barcode #coldata Samples of interest


## Create count matrix for Samples Of Interest cntSOI ####
cntSOI <- exon_mat_hg38 %>% 
  inner_join(exon_mat_macFas6, by = "exon") %>% 
  column_to_rownames("exon") %>% 
  as.matrix() 

write_csv2(as.data.frame(cntSOI), "RNAseq/cnt_matrix_full.csv")

## Check rows of colDataSOI match cols of cntSOI ####
table(rownames(colDataSOI) == colnames(cntSOI))

## Change col & row names to sample Name easier to identify ####
colnames(cntSOI) <- colDataSOI$Name
rownames(colDataSOI)<- colDataSOI$Name
table(rownames(colDataSOI) == colnames(cntSOI))



# SELECT TIME POINTS ------------------------------------------------------

colDataDif <- colDataSOI %>% filter(Time != 3) %>% 
  mutate(Differentiation = factor(case_when(Time %in% c("5","7","9")  ~ "NPC",
                                            Time %in% c("0","1") ~ "iPSC")))
cntDif <- cntSOI[,colDataDif$Name] 
table(colnames(cntDif) == rownames(colDataDif))



# FILTER GENES ------------------------------------------------------------

NPCs <- colDataDif %>%  filter(Differentiation == "NPC") %>% pull(Name)

# we are only going to filter based on NPCs as we will not use iPSCs later: require 6 (28.57%) expression presence for the gene to be kept (was the same though if also adding: | sum(x[iPSCs]>=4))
gene_select <- apply(cntDif >= 1, 1, function(x){sum(x[NPCs]) >= 6})
summary(gene_select)

cntDif_filtered <- cntDif[gene_select,] # count matrix only filtered for gene_select



## Drop-out summary ####
par(mfrow = c(1,2))
#Dropouts
quick_zero <- rowMeans(cntDif_filtered == 0)
s<-summary(quick_zero)
boxplot(quick_zero)
dropout_umi = sum(cntDif_filtered == 0)/(nrow(cntDif_filtered) * ncol(cntDif_filtered))
#median: 0.05, mean: 0.33

# average expression
rwm <- rowMeans(cntDif_filtered)
hist(rwm[rwm<100], 
     freq=FALSE, 
     xlab="Average UMI count per gene", 
     col="grey80", 
     cex.lab=1, 
     cex.axis=0.8, 
     breaks = 50, 
     xaxt="n")

axis(side=1, 
     at=seq(0,100, 10), 
     labels=seq(0,100, 10))

length(rwm[rwm<=2])/length(rwm) # 22% genes show <= 2 UMI counts in average
length(rwm[rwm<=10])/length(rwm)  # 54% genes show <= 10 UMI counts in average

dev.copy2pdf(file = "RNAseq/figures/DropOut_AvgUMIcnt.pdf", width = 8, height = 6)


# Create DDS object ####

table(colnames(cntDif_filtered) == rownames(colDataDif))
ddsDif_filtered <- DESeqDataSetFromMatrix(countData = cntDif_filtered,
                                          colData = colDataDif,
                                          design = ~ Species)

ddsDif_filtered <- estimateSizeFactors(ddsDif_filtered)
summary(sizeFactors(ddsDif_filtered))
ddsDif_filtered <- estimateDispersions(ddsDif_filtered)

par(mfrow = c(1,1))
plotDispEsts(ddsDif_filtered)

dev.copy2pdf(file = "RNAseq/figures/dispEsts.pdf", width = 8, height = 6)

vsd_filtered<-varianceStabilizingTransformation(ddsDif_filtered)
vsdMat_filtered<-assay(vsd_filtered)

vsn::meanSdPlot(vsdMat_filtered)
ggsave("RNAseq/figures/variance_stabilized_counts_clean.png", height = 6, width = 8)

k<-reshape2::melt(as.matrix(vsdMat_filtered))
colnames(k)<-c("Ensembl","full_name","var_stabilized_counts")
ggplot(data = k, aes(x=full_name, y=var_stabilized_counts)) + geom_boxplot() + coord_flip()+ theme(axis.text.y = element_text(size=5))

ggsave("RNAseq/figures/normBoxplot_clean.png", height = 6, width = 8)

# SAVE THE DATA ####
saveRDS(colDataDif, file ="RNAseq/RDS/colData_clean.rds")
saveRDS(cntDif_filtered, file ="RNAseq/RDS/cnt_clean.rds")
saveRDS(ddsDif_filtered, file ="RNAseq/RDS/dds_clean.rds")

