
library(methods)
library(broom)
library(tidyverse)

#setwd("/data/share/htp/pleiotropy/paper_data/")

source("scripts/helper_functions.R")
summarized_expression_all<- readRDS("roadmap_expression_summaries/summarized_expression_allTissues.rds")
DHS_to_gene<-readRDS("CRE_to_Gene/DHS_to_gene.rds")


data<- DHS_to_gene %>% 
  mutate(total = factor(total,level=1:9)) %>% 
  left_join(summarized_expression_all %>% dplyr::select(tissue,gene_id,log2_mean_expression)) 

out<-analyse_permutations(df= data ,n = 30)  
#pp<-make_permut_plots(out)
#fout<-list( simout=out, plots=pp)

saveRDS(out, file = "roadmap_expression_summaries/complex_mixed_model_boot.RDS")

