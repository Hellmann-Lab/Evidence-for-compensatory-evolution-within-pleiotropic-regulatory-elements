
library(tidyverse)
library(data.table)

setwd("/data/share/htp/pleiotropy/paper_data/")

cbust_rnk_clean <- function(file) {
 
  tmp<- fread(file)[ cluster_or_motif == 'motif', ][
                                    ,rnk := frank(-cluster_or_motif_score), by =.(cluster_id,seq_name) ][ ,
                                    .(cscore = sum(cluster_or_motif_score),
                                      minRNK = min(rnk)),
                    by = .(seq_name, cluster_id_or_motif_name) ][
                      , region_id := tstrsplit(seq_name, "@@", fixed = TRUE)[[2]] ][
                      , seq_name := NULL ]  

  return(tmp)
}




# SUMMARIZE TISSUE-EXPRESSED TF DATA --------------------------------------


hum_cbust_files<-list.files("DA_ipsc_npc/cbust/cbust_out/hg38/exprTissues/", pattern="txt", full.names = T)
hum_cbust <-lapply(hum_cbust_files, cbust_rnk_clean) %>% rbindlist()

mac_cbust_files<-list.files("DA_ipsc_npc/cbust/cbust_out/macFas6/exprTissues/", pattern="txt", full.names = T)
mac_cbust <-lapply(mac_cbust_files, cbust_rnk_clean) %>% rbindlist()

combList<- rbindlist( list( M = mac_cbust, H = hum_cbust), idcol = "species")[
  ,minminRNK:= min(minRNK), by = .(cluster_id_or_motif_name, region_id) ][
    ,nMotif := .N, by = .(region_id) ]


data.table::fwrite( combList, file="DA_ipsc_npc/cbust/topNmotifs_exprTissues_unfiltered_dt.csv" )



# SUMMARIZE NPC-EXPRESSED DATA --------------------------------------------



hum_cbust_files<-list.files("DA_ipsc_npc/cbust/cbust_out/hg38/exprIpscNpc/", pattern="txt", full.names = T)
hum_cbust <-lapply(hum_cbust_files, cbust_rnk_clean) %>% rbindlist()

mac_cbust_files<-list.files("DA_ipsc_npc/cbust/cbust_out/macFas6/exprIpscNpc/", pattern="txt", full.names = T)
mac_cbust <-lapply(mac_cbust_files, cbust_rnk_clean) %>% rbindlist()

combList<- rbindlist( list( M = mac_cbust, H = hum_cbust), idcol = "species")[
            ,minminRNK:= min(minRNK), by = .(cluster_id_or_motif_name, region_id) ][
            ,nMotif := .N, by = .(region_id) ]

data.table::fwrite( combList, file="DA_ipsc_npc/cbust/topNmotifs_ipsc_npc_unfiltered_dt.csv" )
