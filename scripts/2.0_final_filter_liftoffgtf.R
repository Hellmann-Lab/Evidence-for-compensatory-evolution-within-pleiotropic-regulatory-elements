library(tidyverse)
library(plyranges)


get_transcript_length<-function(gtf){
  gtf %>% as_tibble %>% group_by(transcript_id,gene_name) %>% summarise(total_bp=sum(width))
}

filter_liftoff_gtf <- function( input_gtf , gtf, species, path="/data/share/htp/EBgrant/genome_data_new/liftoff/"){
  
  ogtf <- plyranges::read_gff(paste0(path, input_gtf)) 
  
  good_genes <- ogtf %>% 
    filter(type == "gene" & is.na(partial_mapping) & is.na(low_identity)) %>% 
    as_tibble() %>% 
    pull(gene_id)
    
  ogtf_filt <- ogtf %>% 
    filter(gene_id %in% good_genes)
  
  too_long<- get_transcript_length( gtf %>% filter(type == "exon")  ) %>% 
    left_join( get_transcript_length( ogtf_filt %>% filter(type == "exon") ),
               by="transcript_id",suffix = c("",".lo")) %>% 
    group_by( gene_name) %>% 
    mutate(   delta = abs(total_bp - total_bp.lo),
              ratio = total_bp.lo/total_bp ) %>% 
    filter(delta>100 & ratio>2) %>% pull(transcript_id)
  
  ogtf2 <- ogtf_filt  %>%  filter(!(transcript_id %in% too_long) )
  newgtf <- file.path(path, paste0("genes.gtf") )
  
  #rtracklayer::export.gff3(ogtf2, con=newgtf )
  rtracklayer::export.gff2(ogtf2, con=newgtf )   #changed to gff2 to get transcript id for cellranger compatibility
}


filter_liftoff_gtf( input_gtf="mf6_liftoff_polished.gtf",
                    gtf = plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/remapping/genomes/hg38/genes.gtf"),
                    species ="mf6")