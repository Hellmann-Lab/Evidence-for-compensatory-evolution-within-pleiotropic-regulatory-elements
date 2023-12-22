library(tidyverse)
library(foreach)
library(doParallel)

#try saving per region id (need to delete after all calculations) --> so far very quick! sbatch --cpus-per-task=30 --wrap="Rscript splitCbust.R"
setwd("/data/share/htp/pleiotropy/paper_data/")

#generate output directories
output<-"DA_ipsc_npc/cbust/RDS/TFBS_position/cbust_sep_exprTissues/"
dir.create(paste0(output,"/hg38/"), showWarnings = F, recursive = T)
dir.create(paste0(output,"/macFas6/"), showWarnings = F, recursive = T)


names_cbust<-c("region_id", "genomic_start__bed", "genomic_end__bed", "cluster_id_or_motif_name", "cluster_or_motif_score", "strand", "seq_name", "relative_start__bed", "relative_end__bed", "seq_number", "cluster_or_motif", "cluster_id", "motif_id", "motif_sequence", "motif_type_contribution_score", "extra_info")

#human cbust files
hum<-list.files("DA_ipsc_npc/cbust/cbust_out/hg38/exprTissues", pattern = ".fa.txt", full.names = T)

cl <- makeCluster(30)
registerDoParallel(cl)  # use multicore, set to the number of cores

foreach (i=1:length(hum)) %dopar% {
  library(tidyverse)
  
  Data1<-setNames(read.table(hum[i],  sep="\t"), names_cbust) %>%
    dplyr::select(seq_name, cluster_or_motif, cluster_id_or_motif_name, cluster_or_motif_score, relative_start__bed, relative_end__bed, strand, motif_sequence) %>%
    #dplyr::filter(cluster_or_motif=="motif") %>%
    group_by(seq_name) %>%
    group_split(seq_name)
  
  lapply(Data1, function(x){saveRDS(dplyr::ungroup(x), paste0(output,"/hg38/",stringr::word(x$seq_name[1],2,2,sep="@@"),".rds"))})
}




#macaque cbust files
mac<-list.files("DA_ipsc_npc/cbust/cbust_out/macFas6/exprTissues", pattern = ".fa.txt", full.names = T)

foreach (i=1:length(mac)) %dopar% {
  library(tidyverse)
    Data1<-setNames(read.table(mac[i],  sep="\t"), names_cbust) %>%
    dplyr::select(seq_name, cluster_or_motif, cluster_id_or_motif_name, cluster_or_motif_score, relative_start__bed, relative_end__bed, strand, motif_sequence) %>%
    #dplyr::filter(cluster_or_motif=="motif") %>%
    group_by(seq_name) %>%
    group_split(seq_name)
  
  lapply(Data1, function(x){saveRDS(dplyr::ungroup(x), paste0(output,"/macFas6/",stringr::word(x$seq_name[1],2,2,sep="@@"),".rds"))})
}

stopCluster(cl)





# make a file with region_id and the chromosome on which it lies for both species

setwd("/data/share/htp/pleiotropy/paper_data/")
library(Biostrings)
library(tidyverse)

# region_id vs file  - recorder ####
human_files<-list.files("DA_ipsc_npc/cbust/fastas/hg38", full.names = F)

region_ids_h<-lapply(human_files, function(x){
  strings<-readDNAStringSet(paste0("DA_ipsc_npc/cbust/fastas/hg38/",x))
  return(data.frame(region_id=stringr::word(names(strings),2,2,sep="@@"), human_file=x))
})

# add macaque files
macaque_files<-list.files("DA_ipsc_npc/cbust/fastas/macFas6", full.names = F)

region_ids_m<-lapply(macaque_files, function(x){
  strings<-readDNAStringSet(paste0("DA_ipsc_npc/cbust/fastas/macFas6/",x))
  return(data.frame(region_id=stringr::word(names(strings),2,2,sep="@@"), macaque_file=x))
})

region_ids_decoder<-inner_join(bind_rows(region_ids_h), bind_rows(region_ids_m))

write.table(region_ids_decoder, "DA_ipsc_npc/cbust/RDS/TFBS_position/selected_ids.txt",
            row.names = F, col.names = F, quote = F)


