

# Functions to compute diversity statistics ####################################
calc_entropy <- function(p_i) {
  H = -sum(p_i*log(p_i))
  return(H)
}
calc_H_beta <- function(H_gamma, H_alpha, weight, overlap,
                        n = 2, option = "weighted"){
  if (option == "equal"){
    if (overlap[1] == 0) {
      H_beta <- log(n)
    } else {
      H_beta <- H_gamma - H_alpha
    }
  }
  
  if (option == "weighted") {
    if (weight[1] == 1) {
      H_beta <- 0 #log(n)
    } else {
      if (overlap[1] == 0) {
        H_beta <- -sum(weight*log(weight))
      } else {
        H_beta <- H_gamma - H_alpha
      }
    }
  }
  
  return(H_beta)
}

calc_Shannon_diff <- function(H_beta, weight, overlap,
                              n = 2, option = "weighted"){
  if (option == "equal"){
    Shannon_diff <- H_beta/log(n)
  }
  
  if (option == "weighted") {
    if (weight[1] == 1) {
      Shannon_diff <- 1 #H_beta/log(n)
    } else {
      Shannon_diff <- H_beta/-sum(weight*log(weight))
    }
  }
  return(Shannon_diff)
}




# run wilcox test for enrichment
run_wilcox_per_PD<-function(df, selCol, selFilt, alt = "greater"){
  
  #, allTissue){
  PD_selected<-df %>% mutate(classes = ifelse({{selCol}} == selFilt, 1, 0)) 
  
  lapply(unique(PD_selected$motif_ID), function(x){
    wilcox.test(motif_density ~ classes, data = PD_selected %>% filter(motif_ID == x),
                alternative = alt) %>% tidy() %>% mutate(motif_ID = x)
  })
}








#alignment stats
get_aln_stats<-function(seqs1, seqs2, aln_mafft){
  
  #for these stats could trim the endings; aka I save the full alignment length as well as the end- trimmed as well as the non-gap part!
  all_aln<-PairwiseAlignmentsSingleSubject(aln_mafft)
  
  aln_stats<-data.frame(nameS1 = names(aln_mafft)[1],
                        nameS2 = names(aln_mafft)[2],
                        score = sapply(all_aln, score),
                        aln_length = sapply(all_aln, nchar),
                        s1_length  = width(seqs1),
                        s2_length  = width(seqs2),
                        nindel = sapply(all_aln, function(x){
                          tmp<-nindel(x)
                          tmp@insertion[1,1] + tmp@deletion[1,1] }),
                        indel_bp = sapply(all_aln, function(x){
                          tmp<-nindel(x)
                          tmp@insertion[1,2] + tmp@deletion[1,2] }))
  
  return(list(all_aln, aln_stats))
}



#align using mafft
#it's weird that i give it "i", it should still work if there is just one region_id, human_file, macaque_file, but should change this for the future again!
align_with_mafft_ATAC<-function(region_id, human_file, macaque_file){
  
  print(paste0("region_id ",region_id))
  hum<-readDNAStringSet(paste0("DA_ipsc_npc/cbust/fastas/hg38/",human_file))
  mac<-readDNAStringSet(paste0("DA_ipsc_npc/cbust/fastas/macFas6/",macaque_file))
  
  hum_red<-hum[grepl(paste0("@@",region_id,"$"),names(hum))]
  mac_red<-mac[grepl(paste0("@@",region_id,"$"),names(mac))]
  both<-c(hum_red,mac_red)
  names(both)<-c(paste0(region_id,"_Homo_sapiens"),
                 paste0(region_id,"_Macaca_fascicularis"))
  dir.create("DA_ipsc_npc/cbust/fastas/alns/tmp", showWarnings = F, recursive = T)
  writeXStringSet(both,paste0("DA_ipsc_npc/cbust/fastas/alns/tmp/region_id_",region_id,".fa"))
  writeXStringSet(both,paste0("DA_ipsc_npc/cbust/fastas/alns/tmp_mafft/region_id_",region_id,".fa"))
  
  com<-paste0("/home/zane/mafft7/mafft-linux64/mafft.bat --adjustdirection --maxiterate 1000 --auto /data/share/htp/pleiotropy/paper_data/DA_ipsc_npc/cbust/fastas/alns/tmp/region_id_",region_id,".fa > /data/share/htp/pleiotropy/paper_data/DA_ipsc_npc/cbust/fastas/alns/tmp_mafft/mafft_",region_id,".fa; wait")
  
  system(com)
  
  #in the cases where sequence was reverse, it adds a _R_ to the name. save this info
  checkAln<-readDNAStringSet(paste0("DA_ipsc_npc/cbust/fastas/alns/tmp_mafft/mafft_",region_id,".fa"))
  revseq1<-setNames(grepl("_R_",names(checkAln)), c("Homo_sapiens","Macaca_fascicularis"))
  
  humvsmac_table<-get_aln_stats(seqs1=hum_red,seqs2=mac_red, aln_mafft = checkAln)
  
  #add a warning if the alignment is below 70% of the shorter seq
  if (humvsmac_table[[2]]$aln_length[1] < 0.7* min(humvsmac_table[[2]]$s1_length[1],humvsmac_table[[2]]$s2_length[1])){
    warning(paste0("region_id_", region_id, " alignment shows a <70% length rel to sequence!"))
  }
  
  strings<-setNames(checkAln, c("Homo_sapiens","Macaca_fascicularis"))
  
  return(list(aln_stats=humvsmac_table[[2]], alignment = strings, revseq1=revseq1,
              hum_red = hum_red, mac_red = mac_red))
  
}



# gap-removing count function
fun1 <- function(x) stringi::stri_length(x) - stringi::stri_count_fixed(x, "-")


names_cbust<-c("region_id", "genomic_start__bed", "genomic_end__bed", "cluster_id_or_motif_name", "cluster_or_motif_score", "strand", "seq_name", "relative_start__bed", "relative_end__bed", "seq_number", "cluster_or_motif", "cluster_id", "motif_id", "motif_sequence", "motif_type_contribution_score", "extra_info")



#read in TFBS motif file
get_TFBS_sites_sep_ATAC<-function(species="human", species_latin = "Homo_sapiens", 
                             #file = human_file,
                             region_id = region_id, 
                             revseq1=revseq1, seq_red=hum_red, 
                             cm = "motif", score_cutoff = 5, cluster_cutoff = 3, 
                             folder = "tmp_cbust_sep"){
  
  TFBS_gr <-readRDS(paste0("DA_ipsc_npc/cbust/RDS/TFBS_position/", folder, "/",species,"/",region_id,".rds")) %>%
    dplyr::filter(cluster_or_motif %in% cm) %>%
    dplyr::filter((cluster_or_motif == "cluster" & cluster_or_motif_score>cluster_cutoff) | (cluster_or_motif == "motif" & cluster_or_motif_score>score_cutoff)) %>%
    dplyr::rowwise() %>%
    #if the seq has been turned around, the cbust positions are actually start = seqwidth - position end; end = seqwidth - position start
    dplyr::mutate(start = ifelse(revseq1[species_latin], (width(seq_red)-relative_end__bed+1), relative_start__bed+1),
                  end = ifelse(revseq1[species_latin], (width(seq_red)-relative_start__bed), relative_end__bed),
                  strand = case_when(revseq1[species_latin] & strand == "+" ~ "-",
                                     revseq1[species_latin] & strand == "-" ~ "+",
                                     T ~ strand),
                  middle_point=round(start+(end-start)/2,0),
                  seqnames = stringr::word(seq_name,2,2,sep="@@")) %>%
    #mutate(seqnames=stringr::word(seqnames,2,2,sep="@@")) %>%
    as_granges()
  return(TFBS_gr)
}





quantify_binding_site_overlaps<-function(tf, aln_df, aln_hum_gr, aln_mac_gr, humTF, macTF){
  
  if (tf %in% humTF$cluster_id_or_motif_name){
    
    alignment_df_tf<-aln_df %>% 
      dplyr::mutate(binding_human= ifelse(pos_seq_Homo_sapiens %in% (filter_by_overlaps(aln_hum_gr,humTF %>% filter(cluster_id_or_motif_name == tf) ) %>% as_tibble() %>% pull(start)), 1, 0),
                    middle_point_human=ifelse(pos_seq_Homo_sapiens %in% (as_tibble(humTF) %>% filter(cluster_id_or_motif_name == tf) %>% pull(middle_point)),1,0))
  } else {
    alignment_df_tf<-aln_df %>% dplyr::mutate(binding_human =0, middle_point_human = 0)
  }
  
  #do the same for macaque
  if (tf %in% macTF$cluster_id_or_motif_name){
    
    alignment_df_tf<-alignment_df_tf %>% 
      dplyr::mutate(binding_macaque= ifelse(pos_seq_Macaca_fascicularis %in% (filter_by_overlaps(aln_mac_gr, macTF %>% filter(cluster_id_or_motif_name == tf) ) %>% as_tibble() %>% pull(start)), 1, 0),
                    middle_point_macaque=ifelse(pos_seq_Macaca_fascicularis %in% (as_tibble(macTF) %>% filter(cluster_id_or_motif_name == tf) %>% pull(middle_point)),1,0))
  } else {
    alignment_df_tf<-alignment_df_tf %>% dplyr::mutate(binding_macaque=0, middle_point_macaque = 0)
  }
  return(alignment_df_tf)
}


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}






get_mafft_aln_ATAC<-function(region_id, human_file, macaque_file){
  
  print(paste0("region_id ",region_id))
  hum<-readDNAStringSet(paste0("DA_ipsc_npc/cbust/fastas/hg38/",human_file))
  mac<-readDNAStringSet(paste0("DA_ipsc_npc/cbust/fastas/macFas6",macaque_file))
  
  hum_red<-hum[names(hum)==region_id]
  mac_red<-mac[names(mac)==region_id]
  both<-c(hum_red,mac_red)
  names(both)<-c(paste0(region_id,"_Homo_sapiens"),
                 paste0(region_id,"_Macaca_fascicularis"))
  
  #in the cases where sequence was reverse, it adds a _R_ to the name. save this info
  checkAln<-readDNAStringSet(paste0("DA_ipsc_npc/cbust/fastas/alns/tmp_mafft/mafft_",region_id,".fa"))
  revseq1<-setNames(grepl("_R_",names(checkAln)), c("Homo_sapiens","Macaca_fascicularis"))
  
  humvsmac_table<-get_aln_stats(seqs1=hum_red,seqs2=mac_red, aln_mafft = checkAln)
  
  #add a warning if the alignment is below 70% of the shorter seq
  if (humvsmac_table[[2]]$aln_length[1] < 0.7* min(humvsmac_table[[2]]$s1_length[1],humvsmac_table[[2]]$s2_length[1])){
    warning(paste0("region_id_", region_id, " alignment shows a <70% length rel to sequence!"))
  }
  
  strings<-setNames(checkAln, c("Homo_sapiens","Macaca_fascicularis"))
  
  return(list(aln_stats=humvsmac_table[[2]], alignment = strings, revseq1=revseq1,
              hum_red = hum_red, mac_red = mac_red))
  
}





calculate_per_region_TF_stats<-function(dff) {
  dff %>% dplyr::select(region_id,species,cscore,cluster_id_or_motif_name) %>% 
    pivot_wider( names_from = species, values_from = c("cscore") , values_fill = 0) %>% 
    group_by(region_id) %>% 
    mutate( pi_avg = (M+H) / sum(M+H),
            pi_H = ifelse(H >0, H/ sum(H),0),
            pi_M = ifelse(M >0, M/ sum(M),0) ) %>% 
    summarise( nmotifs = length(cluster_id_or_motif_name),
               d_eucl = sqrt(sum((M-H)^2)),
               d_canb =  sum( abs(M-H)/(M+H) ),
               H_gamma = -sum( pi_avg * log(pi_avg)),
               H_alpha = -sum( pi_M * log(pi_M), na.rm = T) - sum( pi_H * log(pi_H), na.rm = T),
               shannon_diff = H_gamma - H_alpha/2,
               diversityH = exp(   - sum( pi_H * log(pi_H), na.rm = T) ) ,
               canb_corr =d_canb/nmotifs)
}

