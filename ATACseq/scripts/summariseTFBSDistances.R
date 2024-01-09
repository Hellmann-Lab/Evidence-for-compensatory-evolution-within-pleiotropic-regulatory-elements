
#setwd("/data/share/htp/pleiotropy/paper_data/")
library(tidyverse)


combis<-list.files("ATACseq/cbust/RDS/TFBS_position/combis_exprTissues_10perc/", full.names = TRUE)

TFBS_dist<-lapply(combis, function(x){
  combi<-readRDS(x)
  #print(x)
  if (dim(combi$TF_out)[1]>0){
    
    summarized_TF_dist<-data.frame(mean_agreement=mean(combi$TF_out$n_bases_agreement/combi$TF_out$n_bases),
                                   mean_dist=mean(combi$TF_out$mean_dist, na.rm=TRUE),
                                   total_TFs=length(combi$TF_out$TF),
                                   n_sites_hum_mean=mean(combi$TF_out$n_sites_hum),
                                   n_sites_hum_sum=sum(combi$TF_out$n_sites_hum),
                                   n_sites_mac_mean=mean(combi$TF_out$n_sites_mac),
                                   n_sites_mac_sum=sum(combi$TF_out$n_sites_mac),
                                   region_id=unique(combi$TF_out$region_id)) %>%
      left_join(combi$alnQ)
    return(summarized_TF_dist)
  } else {
    return(data.frame(region_id=gsub("ATACseq/cbust/RDS/TFBS_position/combis_exprTissues_10perc/region_id_","", x)))
  }
})

TFBS_dist<-lapply(TFBS_dist, function(x){x %>% dplyr::mutate(region_id = as.factor(as.character(region_id)))}) %>% bind_rows()

saveRDS(TF_dist, "ATACseq/cbust/RDS/TFBS_position/summarized_TFBS_distances_exprTissues_10perc.rds")

