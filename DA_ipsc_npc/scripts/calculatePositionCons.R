
# let's do alignments
library(Biostrings)
library(tidyverse)
library(stringi)
library(plyranges)
library(stringr)
library(magrittr)
library(readr)

setwd("/data/share/htp/pleiotropy/paper_data/")
source("DA_ipsc_npc/scripts/TFBS_analysis_functions_ATAC.R")
outdir<-"DA_ipsc_npc/cbust/RDS/TFBS_position/"
#extension<-"_exprTissues"
extension<-"_exprTissues_10perc"
dir.create(paste0(outdir,"regions_noTFs",extension), showWarnings = F, recursive = T)
dir.create(paste0(outdir,"combis",extension), showWarnings = F, recursive = T)



args <- commandArgs(trailingOnly = TRUE)
print(args)

region_id<-args[1]
human_file<-args[2]
macaque_file<-args[3]
if (length(args)>3){
  regid = region_id
  # dd_10perc_regid<-data.table::fread("DA_ipsc_npc/cbust/top10perc_motifs_exprTissues_dt.csv")[region_id == regid,]
  dd_10perc_regid<-data.table::fread(args[4])[region_id == regid,]
}


# example
#region_id<-124
#human_file<-'hg38_1.fa'
#macaque_file<-'macFas6_1.fa'


print(paste0("### Region_id ",region_id))
aln_out<-align_with_mafft_ATAC(region_id, human_file, macaque_file)
alnQ<-aln_out[["aln_stats"]] %>% mutate(region_id = region_id)
strings<-aln_out[["alignment"]]
revseq1<-aln_out[["revseq1"]]


#collect all the positions
position_list<-list()

for (p in 1:length(strings)){
  pos_split<-strsplit(x=as.character(strings[p]), split=character(0)) %>%
    as.data.frame() %>%
    rownames_to_column("pos_alignment")
  
  char_seq<-DNAString(as.character(strings[p]))
  
  for (k in 1:nrow(pos_split)){
    pos_split[k,paste0("pos_seq_",names(strings[p]))]<-fun1(char_seq[1:k])
    position_list[[names(strings[p])]]<-pos_split %>% column_to_rownames("pos_alignment")
  }
}

alignment_df<-bind_cols(position_list) %>% 
  rownames_to_column("pos_alignment") %>%
  mutate(seqnames = region_id)
#we probably only want to look at TFBS sites that theoretically could've been
#included in both species: remove gaps? yes.
aln_dfs<-alignment_df %>% 
  mutate(region_id = region_id) %>%
  filter(Macaca_fascicularis!="-", Homo_sapiens !="-") 

mismatches<-aln_dfs %>% filter(Homo_sapiens !=Macaca_fascicularis)
#pull out mismatches myself since the function does not always work..
#validated on a few examples that this gives the same as mismatch summary
alnQ$mismatches<-dim(mismatches)[1]

#read in TFBS motif file
hum_TFBS_gr<-get_TFBS_sites_sep_ATAC(species="hg38", 
                                     species_latin = "Homo_sapiens",
                                     region_id = region_id,
                                     revseq1 = revseq1, 
                                     seq_red = aln_out[["hum_red"]],
                                     folder = "cbust_sep_exprTissues", 
                                     cm = "motif",
                                     score_cutoff = 3)

if (length(args)>3 & length(hum_TFBS_gr)>0){
  hum_TFBS_gr<-hum_TFBS_gr %>% dplyr::filter(cluster_id_or_motif_name %in% dd_10perc_regid$cluster_id_or_motif_name)
}


mac_TFBS_gr<-get_TFBS_sites_sep_ATAC(species="macFas6", 
                                     species_latin = "Macaca_fascicularis", 
                                     region_id = region_id,
                                     revseq1 = revseq1, 
                                     seq_red = aln_out[["mac_red"]],
                                     folder = "cbust_sep_exprTissues",
                                     cm = "motif",
                                     score_cutoff = 3)

if (length(args)>3 & length(mac_TFBS_gr)>0){
  mac_TFBS_gr<-mac_TFBS_gr %>% dplyr::filter(cluster_id_or_motif_name %in% dd_10perc_regid$cluster_id_or_motif_name)
}

if (length(hum_TFBS_gr)==0 | length(mac_TFBS_gr)==0){
  
  saveRDS(bindingTFs, paste0(outdir,"regions_noTFs",extension,"/region_id_",region_id,".rds"))
  
} else {
  
  bindingTFs<-data.frame(TF=unique(c(hum_TFBS_gr$cluster_id_or_motif_name, mac_TFBS_gr$cluster_id_or_motif_name)))
  
  alignment_df_tf_gr_hum<-alignment_df %>% 
    transmute(seqnames, start = pos_seq_Homo_sapiens, width = 1, strand = "*") %>% 
    as_granges()
  
  alignment_df_tf_gr_mac<-alignment_df %>% 
    transmute(seqnames, start = pos_seq_Macaca_fascicularis, width = 1, strand = "*") %>% 
    as_granges()
  
  for (tf in bindingTFs$TF){
    #print(tf)
    binding_site_agreement<-quantify_binding_site_overlaps(tf = tf, aln_df = alignment_df,
                                                           aln_hum_gr = alignment_df_tf_gr_hum, 
                                                           aln_mac_gr = alignment_df_tf_gr_mac, 
                                                           humTF = hum_TFBS_gr, 
                                                           macTF = mac_TFBS_gr) %>%
      filter(binding_human == 1 | binding_macaque ==1 ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(agree = ifelse(binding_human == 1 & binding_macaque == 1, 1, 0))
    
    
    #add stats
    bindingTFs$mean_agreement[bindingTFs$TF==tf]<-mean(binding_site_agreement$agree)
    bindingTFs$n_bases[bindingTFs$TF==tf]<-length(binding_site_agreement$agree)
    bindingTFs$n_bases_agreement[bindingTFs$TF==tf]<-length(binding_site_agreement$agree[binding_site_agreement$agree==1])
    bindingTFs$n_sites_hum[bindingTFs$TF==tf]<-sum(binding_site_agreement$middle_point_human)
    bindingTFs$n_sites_mac[bindingTFs$TF==tf]<-sum(binding_site_agreement$middle_point_macaque)
    
    bindingTFs$n_unique_sites[bindingTFs$TF==tf]<-sum(binding_site_agreement$n_sites)
    bindingTFs$n_total_sites[bindingTFs$TF==tf]<-sum(binding_site_agreement$n_sites_agreement)
    bindingTFs$n_sites_agreement[bindingTFs$TF==tf]<-length(binding_site_agreement$n_sites_agreement[binding_site_agreement$n_sites_agreement==2])
    
    #project all onto human? or does this introduce some biases?
    if (max(binding_site_agreement$middle_point_human==1) & max(binding_site_agreement$middle_point_macaque==1)){
      
      hum_midp<-binding_site_agreement %>% 
        filter(middle_point_human == 1) %>% 
        as_granges(seqnames = 1, start = pos_seq_Homo_sapiens, width = 1, strand = "*")
      
      mac_midp<-binding_site_agreement %>% 
        filter(middle_point_macaque == 1) %>% 
        as_granges(seqnames = 1, start = pos_seq_Homo_sapiens, width = 1, strand = "*")
      
      dist<-add_nearest_distance(hum_midp, mac_midp, name = "distance")
      bindingTFs$mean_dist[bindingTFs$TF==tf]<-mean(dist$distance)
    }
  }
  aln_list<-bindingTFs %>% dplyr::mutate(region_id = region_id)
  saveRDS(list(TF_out=aln_list, alnQ=alnQ, aln_dfs=aln_dfs), paste0(outdir,"combis",extension,"/region_id_",region_id,".rds"))
}


