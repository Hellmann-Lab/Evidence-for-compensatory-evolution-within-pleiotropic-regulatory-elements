
# let's do alignments
library(Biostrings)
library(tidyverse)
#library(stringi)
library(plyranges)
library(stringr)
library(magrittr)
library(readr)
library(wesanderson)
library(cowplot)
library(RColorBrewer)

basic_theme_ins<-  theme(axis.title=element_text(size=9.5),
                         axis.text = element_text(color="black", size=8.5),
                         legend.title = element_text(size = 7.5), 
                         legend.text = element_text(size = 8),
                         legend.key.size = unit(0.5, "lines"),
                         legend.margin=margin(0.1,0.1,0.1,0.1),
                         legend.box.margin=margin(0.1,0.1,0.1,0.1),
                         strip.text.x = element_text(size = 10))

#setwd("/data/share/htp/pleiotropy/paper_data/")
#source("/data/share/htp/pleiotropy/paper_data/ATACseq/scripts/TFBS_analysis_functions_ATAC.R")


# the input (separated rds) are already filtered for motif, probably to reduce size
# translate macaque positions to the human sequence
# need to bind using macaque coords (adjusted for alignment) and then pull out the respective human coords
# we want the start and the end coords
# still need to remove 500 bp from each side of the cbust output!

get_TFBS_positions_ATAC<-function(get_region_id, human_file, macaque_file,
                             TF_score_cutoff = 5, cluster_score_cutoff = 0){
  
  aln_out<-align_with_mafft_ATAC(get_region_id, human_file, macaque_file)
  alnQ<-aln_out[["aln_stats"]] %>% mutate(region_id = get_region_id)
  strings<-aln_out[["alignment"]]
  revseq1<-aln_out[["revseq1"]]
  
  genomic_pos<-readRDS("ATACseq/cbust/RDS/jamm_inATAC.rds") %>%
    dplyr::mutate(region_id = as.numeric(as.character(region_id))) %>%
    dplyr::filter(region_id == get_region_id) %>%
    # the extended seqs are aligned here..
    dplyr::transmute(seqnames, 
                     start_DHS = start, 
                     end_DHS = end, 
                     start = ifelse(start-500>=0,start-500,0),
                     # too lazy to control for the chr ends right now..
                     # unlikely we encounter these in a few examples.
                     end = end+500)
  
  #collect all the positions
  position_list<-list()
  
  for (p in 1:length(strings)){
    pos_split<-strsplit(x=as.character(strings[p]), 
                        split=character(0)) %>%
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
    mutate(seqnames = as.factor(get_region_id)) %>%
    rowwise() %>%
    mutate(pos_hg38 = ifelse(revseq1["Homo_sapiens"],
                             genomic_pos$end[1]-pos_seq_Homo_sapiens+1, 
                             genomic_pos$start[1]+pos_seq_Homo_sapiens-1),
           chr=genomic_pos$seqnames[1])
  
  
  mac_TFBS_inHum<-get_TFBS_sites_sep_ATAC(species="macFas6", 
                                     species_latin = "Macaca_fascicularis", 
                                     region_id = get_region_id,
                                     revseq1 = revseq1, 
                                     seq_red = aln_out[["mac_red"]],
                                     cm =c("cluster", "motif"), 
                                     score_cutoff = TF_score_cutoff,
                                     cluster_cutoff = cluster_score_cutoff,
                                     folder = "cbust_sep_exprTissues") %>%
    as_tibble() %>%
    #mutate(seqnames=stringr::word(seqnames,2,2,sep="@@")) %>%
    left_join(alignment_df %>% filter(Macaca_fascicularis !="-") %>%
                dplyr::transmute(seqnames, 
                                 start = pos_seq_Macaca_fascicularis,
                                 pos_alignment, 
                                 pos_seq_Homo_sapiens,
                                 pos_hg38, 
                                 chr)) %>%
    left_join(alignment_df %>% filter(Macaca_fascicularis !="-") %>%
                dplyr::transmute(seqnames, 
                                 end = pos_seq_Macaca_fascicularis, 
                                 pos_alignment, 
                                 pos_seq_Homo_sapiens, 
                                 pos_hg38), 
              by = c("seqnames","end"), suffix=c(".start",".end"))
  
  
  hum_TFBS_inHum<-get_TFBS_sites_sep_ATAC(species="hg38", 
                                     species_latin = "Homo_sapiens",
                                     region_id = get_region_id,
                                     revseq1 = revseq1, 
                                     seq_red = aln_out[["hum_red"]],
                                     cm =c("cluster", "motif"), 
                                     score_cutoff = TF_score_cutoff,
                                     cluster_cutoff = cluster_score_cutoff,
                                     folder = "cbust_sep_exprTissues") %>%
    as_tibble() %>%
    #mutate(seqnames=stringr::word(seqnames,2,2,sep="@@")) %>%
    left_join(alignment_df %>% filter(Homo_sapiens !="-") %>%
                dplyr::transmute(seqnames, 
                                 start = pos_seq_Homo_sapiens, 
                                 pos_hg38, 
                                 chr)) %>%
    left_join(alignment_df %>% filter(Homo_sapiens !="-") %>%
                dplyr::transmute(seqnames, 
                                 end = pos_seq_Homo_sapiens, 
                                 pos_hg38), 
              by = c("seqnames","end"), suffix=c(".start",".end"))
  
  
  # combine
  binding_both<-bind_rows(mac_TFBS_inHum %>%
                            dplyr::transmute(seqnames, 
                                             cluster_or_motif,
                                             cluster_id_or_motif_name, 
                                             cluster_or_motif_score,
                                             strand,
                                             motif_sequence,
                                             chr, 
                                             pos_hg38.start, 
                                             pos_hg38.end, 
                                             species = "macaque"),
                          hum_TFBS_inHum %>%
                            dplyr::transmute(seqnames, 
                                             cluster_or_motif, 
                                             cluster_id_or_motif_name, 
                                             cluster_or_motif_score, 
                                             strand,
                                             motif_sequence, 
                                             chr, 
                                             pos_hg38.start, 
                                             pos_hg38.end, 
                                             species = "human")) %>%
    group_by(cluster_id_or_motif_name) %>%
    mutate(sum_score = sum(cluster_or_motif_score), 
           sum_score) %>%
    ungroup()
  
  return(list(binding_both = binding_both, 
              aln = strings,
              genomic_pos = genomic_pos))
}



plot_bindingSites<-function(outfile = outfile, get_region_id, filt = "top_TFs", 
                            n_filt = 5, TFs = NULL,
                            padding = 0, fill_highlight = "lightyellow",
                            expressedTF_string = NULL){
  
  if (!is.null(expressedTF_string)){
    outfile$binding_both<-outfile$binding_both %>%
      filter(grepl(expressedTF_string, cluster_id_or_motif_name) | cluster_or_motif == "cluster")
  }
  
  if (!is.null(TFs)){
    top<-outfile$binding_both %>% 
      filter(grepl(TFs, cluster_id_or_motif_name))
  } else {
    if (filt == "top_TFs"){
      top<-outfile$binding_both %>% 
        filter(cluster_or_motif == "motif") %>%
        distinct(cluster_id_or_motif_name, sum_score, .keep_all=T) %>%
        slice_max(order_by=sum_score, n=n_filt)
    }
    if (filt == "top_sites"){
      top<-outfile$binding_both %>% 
        filter(cluster_or_motif == "motif") %>%
        group_by(species) %>%
        slice_max(order_by=cluster_or_motif_score, n=n_filt) %>%
        ungroup()
    }}
  
  motifs<-outfile$binding_both %>% 
    filter(cluster_id_or_motif_name %in% top$cluster_id_or_motif_name) %>%
    mutate(cluster_id_or_motif_name = stringr::word(cluster_id_or_motif_name, 2,2,"_")) 
  
  motifs %>%
    ggplot() + 
    # cluster positions
    geom_rect(data=outfile$binding_both %>% filter(cluster_or_motif == "cluster"), 
              mapping=aes(xmin=pos_hg38.start, xmax=pos_hg38.end, 
                          ymin=-max(motifs$cluster_or_motif_score)-0.5,
                          ymax=max(motifs$cluster_or_motif_score)-0.5), 
              fill=fill_highlight, alpha=0.9) +
    geom_hline(yintercept=0) +
    # binding site positions
    geom_rect(aes(xmin=pos_hg38.start, xmax=pos_hg38.end, 
                  ymin=0, ymax=as.numeric(paste0(strand,cluster_or_motif_score)),
                  fill = cluster_id_or_motif_name)) +
    theme_light(base_size = 9)+
    scale_x_continuous(breaks = scales::pretty_breaks(n=3),
                       limits = c(outfile$genomic_pos$start_DHS[1]-padding,
                                  outfile$genomic_pos$end_DHS[1]+padding))+
    coord_cartesian(expand=F)+
    xlab(outfile$genomic_pos$seqnames[1]) +
    ylab("Stranded motif score") +
    facet_grid(species~.)+
    scale_size(guide = "none")+
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(top$cluster_id_or_motif_name))))+
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0),
          legend.key.size = unit(0.5,"line"))
}




plot_bindingSites_alt_noStrand<-function(outfile = outfile, get_region_id,
                                         n_filt = 5, TFs = NULL,
                                         padding = 0, fill_highlight = "lightyellow",
                                         bindingCols=NULL, consCols=NULL){
  
  
  if (!is.null(TFs)){
    top<-outfile$binding_both %>% 
      filter(grepl(TFs, cluster_id_or_motif_name))
  } 
  
  motifs<-outfile$binding_both %>% 
    filter(cluster_id_or_motif_name %in% top$cluster_id_or_motif_name) %>%
    mutate(cluster_id_or_motif_name = stringr::word(cluster_id_or_motif_name, 2,2,"_")) %>%
    group_by(cluster_id_or_motif_name, strand, pos_hg38.start) %>%
    mutate(n_species = length(cluster_id_or_motif_name),
           conserved = ifelse(n_species == 2, "conserved", "species-\nspecific")) %>%
    ungroup() %>%
    group_by(cluster_id_or_motif_name, strand, pos_hg38.end) %>%
    mutate(n_species = length(cluster_id_or_motif_name),
           conserved = ifelse(n_species == 1 & conserved!="conserved", 
                              "species-\nspecific",
                              "conserved")) %>%
    ungroup() %>%
    group_by(cluster_id_or_motif_name) %>%
    mutate(keep = ifelse("species-\nspecific" %in% conserved, T, F)) %>%
    #mutate(keep = ifelse(max(cluster_or_motif_score>=5), T, F)) %>%
    filter(keep) %>%
    # collapse strand info to max 
    group_by(species, cluster_id_or_motif_name, pos_hg38.start, pos_hg38.end) %>%
    mutate(ismax_double = ifelse(length(unique(strand))==2 & cluster_or_motif_score == max(cluster_or_motif_score), T,F),
           rm = ifelse(sum(ismax_double)>1 & strand =="-", T,F)) %>% 
    ungroup() %>%
    filter(!rm)
  
  
  pp_main<-motifs %>%
    mutate(cluster_or_motif_score_species=ifelse(species=="human", cluster_or_motif_score, -cluster_or_motif_score)) %>%
    ggplot() + 
    # cluster positions
    geom_rect(data=outfile$binding_both %>% filter(cluster_or_motif == "cluster"), 
              mapping=aes(xmin=pos_hg38.start, 
                          xmax=pos_hg38.end, 
                          #ymin=0,
                          ymin=-abs(max(motifs$cluster_or_motif_score))-1,
                          ymax=abs(max(motifs$cluster_or_motif_score))+1), 
              fill=fill_highlight, alpha=0.9) +
    geom_hline(yintercept=0) +
    # binding site positions
    geom_rect(aes(xmin=pos_hg38.start, xmax=pos_hg38.end, 
                  ymin=0, 
                  ymax=as.numeric(cluster_or_motif_score_species),
                  fill=cluster_id_or_motif_name),
              linewidth = rel(0.08)) +
    theme_light(base_size = 9)+
    scale_x_continuous(breaks = scales::pretty_breaks(n=3),
                       limits = c(outfile$genomic_pos$start_DHS[1]-padding,
                                  outfile$genomic_pos$end_DHS[1]+padding))+
    scale_y_continuous(limits=c(-9,9), labels = abs)+
    coord_cartesian(expand=F)+
    xlab(outfile$genomic_pos$seqnames[1]) +
    ylab("Motif score") +
    theme(#axis.title.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          strip.text= element_text(size = rel(0.9), color = "black"),
          strip.background = element_rect(fill = "grey85", colour = "grey20"),
          panel.border = element_rect(colour = "grey20", linewidth = rel(0.7)),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(0,0,0,0),
          legend.spacing.y = unit(0.1, "cm"),
          legend.spacing.x = unit(0.1, "cm"),
          legend.key.size = unit(0.35, "cm")) +
    xlab("Alignment position")

  
  
  if (!is.null(bindingCols)){
    pp_main<-pp_main+scale_fill_manual(values=bindingCols)
  }
  if (!is.null(consCols)){
    pp_main<-pp_main+scale_color_manual(values=consCols)
  }
  return(pp_main)
}



plot_bindingSites_alt<-function(outfile = outfile, get_region_id,
                            n_filt = 5, TFs = NULL,
                            padding = 0, fill_highlight = "lightyellow",
                            bindingCols=NULL, consCols=NULL){
  
  
  if (!is.null(TFs)){
    top<-outfile$binding_both %>% 
      filter(grepl(TFs, cluster_id_or_motif_name))
  } 
  
  motifs<-outfile$binding_both %>% 
    filter(cluster_id_or_motif_name %in% top$cluster_id_or_motif_name) %>%
    mutate(cluster_id_or_motif_name = stringr::word(cluster_id_or_motif_name, 2,2,"_")) %>%
    group_by(cluster_id_or_motif_name, strand, pos_hg38.start) %>%
    mutate(n_species = length(cluster_id_or_motif_name),
           conserved = ifelse(n_species == 2, "conserved", "species-\nspecific")) %>%
    ungroup() %>%
    group_by(cluster_id_or_motif_name, strand, pos_hg38.end) %>%
    mutate(n_species = length(cluster_id_or_motif_name),
           conserved = ifelse(n_species == 1 & conserved!="conserved", 
                              "species-\nspecific",
                              "conserved")) %>%
    ungroup() %>%
    group_by(cluster_id_or_motif_name) %>%
    mutate(keep = ifelse("species-\nspecific" %in% conserved, T, F)) %>%
    #mutate(keep = ifelse(max(cluster_or_motif_score>=5), T, F)) %>%
    filter(keep)
  
  
  # just to get the legend.. 
  pp<-motifs %>%
    ggplot() + 
    geom_line(aes(x=pos_hg38.start, 
                  y=0, 
                  color = cluster_id_or_motif_name, 
                  linetype=conserved))+
    theme_light(base_size=9) +
    scale_linetype_manual(values=c(3,1)) +
    theme(legend.title = element_blank(),
          legend.position = "top",
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(0,0,0,0),
          legend.spacing.y = unit(0.1, "cm"),
          legend.spacing.x = unit(0.1, "cm"),
          legend.key.size = unit(0.35, "cm"),
          legend.box = "vertical")
  
  if (!is.null(bindingCols)){
    pp<-pp+scale_color_manual(values=bindingCols)
  }
  
  leg<-get_legend(pp)
  pp_main<-motifs %>%
    #filter(cluster_id_or_motif_name=="ZEB1") %>%
    ggplot() + 
    # cluster positions
    geom_rect(data=outfile$binding_both %>% filter(cluster_or_motif == "cluster"), 
              mapping=aes(xmin=pos_hg38.start, xmax=pos_hg38.end, 
                          ymin=-max(motifs$cluster_or_motif_score)-1,
                          ymax=max(motifs$cluster_or_motif_score)+1), 
              fill=fill_highlight, alpha=0.9) +
    geom_hline(yintercept=0) +
    # binding site positions
    geom_rect(aes(xmin=pos_hg38.start, xmax=pos_hg38.end, 
                  ymin=0, ymax=as.numeric(paste0(strand,cluster_or_motif_score)),
                  #color = cluster_id_or_motif_name, 
                  #fill=conserved,
                  alpha = conserved, 
                  linetype=conserved,
                  fill=cluster_id_or_motif_name),
              color = "grey20",
              linewidth = rel(0.08)) +
    theme_light(base_size = 9)+
    scale_x_continuous(breaks = scales::pretty_breaks(n=3),
                       limits = c(outfile$genomic_pos$start_DHS[1]-padding,
                                  outfile$genomic_pos$end_DHS[1]+padding))+
    coord_cartesian(expand=F)+
    xlab(outfile$genomic_pos$seqnames[1]) +
    ylab("Stranded motif score") +
    facet_grid(species~.)+
    scale_alpha_manual(values=c(0.4,1))+
    scale_linetype_manual(values=c(3,1))+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(), 
          strip.text= element_text(size = rel(0.9), color = "black"),
          strip.background = element_rect(fill = "grey85", colour = "grey20"),
          panel.border = element_rect(colour = "grey20", linewidth = rel(0.7)),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          #legend.position = "left",
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(0,0,0,0),
          legend.spacing.y = unit(0.1, "cm"),
          legend.spacing.x = unit(0.1, "cm"),
          legend.key.size = unit(0.35, "cm"))

  
  if (!is.null(bindingCols)){
    pp_main<-pp_main+scale_fill_manual(values=bindingCols)
  }
  if (!is.null(consCols)){
    pp_main<-pp_main+scale_color_manual(values=consCols)
  }
  return(list(legend=leg, binding_main = pp_main))
}




plot_TFBScoverage_ATAC<-function(outfile = outfile, get_region_id, speciesColors, padding = 0, 
                                 highlight = F, fill_highlight = "lightyellow", fill_alpha=0.5, 
                                 tmpdir="ATACseq/cbust/RDS/TFBS_coverage/tmp/"){
  
  top<-outfile$binding_both
  
  coord_DHS<-data.frame(seqnames=top$chr[1], 
                        start=outfile$genomic_pos$start_DHS[1]-padding, 
                        end=outfile$genomic_pos$end_DHS[1]+padding, 
                        region_id=get_region_id)
  

  dir.create(file.path(tmpdir), showWarnings = FALSE)
  write.table(coord_DHS, paste0(tmpdir,"/DHS_coord.bed"), col.names = F, row.names = F, quote = F, sep = "\t")
  
  # need to do this per species
  for (sp in c("human", "macaque")){
    write.table(top %>% filter(species==sp) %>% dplyr::select(chr, pos_hg38.start, pos_hg38.end), 
                paste0(tmpdir,"/TFBS_coord_",sp,".bed"), 
                col.names = F, row.names = F, quote = F, sep = "\t")
    
    # calculate TFBS density at each position
    system(paste0("bedtools coverage -d -a ",tmpdir,"/DHS_coord.bed -b ",tmpdir,"/TFBS_coord_",sp,".bed > ",tmpdir,"/cov_TFBS_",sp,".bed; wait"))
  }
  
  # read in the data
  cov<-setNames(bind_rows(
    read.table(paste0(tmpdir,"/cov_TFBS_human.bed")) %>%
      mutate(species = "human"),
    read.table(paste0(tmpdir,"/cov_TFBS_macaque.bed")) %>%
      mutate(species = "macaque")),
    c("chr","start","end","region_id","position","coverage","species")) %>% 
    mutate(position=position+coord_DHS$start[1])
  
  # plot
  if (highlight){
    hl<-geom_rect(coord_DHS, 
                  mapping=aes(xmin = start, 
                              xmax = end, 
                              ymin = 0, 
                              ymax = max(cov$coverage)*1.05, 
                              fill=fill_highlight, alpha = fill_alpha))
    
    hl<-geom_rect(data=outfile$binding_both %>% filter(cluster_or_motif == "cluster"), 
                      mapping=aes(xmin=pos_hg38.start, xmax=pos_hg38.end, 
                                  ymin = 0, 
                                  ymax = max(cov$coverage)*1.05), 
                      fill=fill_highlight, alpha=0.9) 
    
    
  } else {
    hl<-NULL
  }
  
  motifs<-outfile$binding_both %>% 
    filter(cluster_id_or_motif_name %in% top$cluster_id_or_motif_name) %>%
    mutate(cluster_id_or_motif_name = stringr::word(cluster_id_or_motif_name, 2,2,"_")) %>%
    group_by(cluster_id_or_motif_name, strand, pos_hg38.start, pos_hg38.end) %>%
    mutate(n_species = length(cluster_id_or_motif_name)) %>%
    ungroup()
  
  ggplot(cov)+
    hl+
    geom_line(aes(x = position, y = coverage, color = species, group = species), linewidth=rel(0.3)) +
    theme_light(base_size = 9) +
    ylab("TFBS density") +
    xlab(coord_DHS$seqnames[1]) +
    theme(panel.grid = element_blank(), 
          legend.position = "none",
          legend.title = element_blank())+
    scale_x_continuous(breaks = scales::pretty_breaks(n=3),
                       limits = c(min(cov$position),
                                  max(cov$position)), 
                       expand = c(0,0))+
    scale_y_continuous(limits = c(0, max(cov$coverage)*1.05), 
                       expand = c(0,0))+
    guides(color = guide_legend(nrow = 1))+
    scale_color_manual(values = speciesColors) +
    facet_grid(species~.)
}





plot_phastCons<-function(bigWigFile, coord, get_region_id, 
                         padding = 0, fill_col = "white",
                         highlight = F, 
                         fill_highlight = "lightyellow", fill_alpha=0.5){
  
  # coordinate files
  coord_DHS<-coord %>% filter(region_id == get_region_id) %>%
    mutate(start = start - padding, end = end + padding)
  
  # get bigwig
  options(ucscChromosomeNames=FALSE)
  phastC<-rtracklayer::import(bigWigFile, 
                              which= as_granges(coord_DHS),
                              as="NumericList") %>%
    data.frame() %>%
    mutate(position=coord_DHS$start[1]:(length(value)-1+coord_DHS$start[1]))
  
  # plot
  if (highlight){
    coord_DHS_orig<-coord %>% filter(region_id == get_region_id) 
    hl<-geom_rect(coord_DHS_orig, 
                  mapping=aes(xmin = start, 
                              xmax = end, 
                              ymin = 0, 
                              ymax = 1), 
                  fill=fill_highlight, alpha = fill_alpha)
  } else {
    hl<-NULL
  }
  ggplot(phastC)+
    hl+
    geom_line(aes(x = position, y = value, group = 1), linewidth=rel(0.3))+
    #geom_ribbon(aes(x = position, y = value, ymin = 0, ymax = value), 
    #            fill = fill_col, position = "identity") +
    theme_light(base_size = 9)+
    ylab("phastCons")+
    xlab(coord_DHS$seqnames[1]) +
    theme(panel.grid = element_blank())+
    scale_x_continuous(breaks = scales::pretty_breaks(n=3),
                       limits = c(min(phastC$position),
                                  max(phastC$position)), 
                       expand = c(0,0))+
    scale_y_continuous(limits = c(0,1), 
                       expand = c(0,0))
}





plot_ATAC<-function(genome, bw_sample, get_region_id, fill_col, padding = 0, 
                    highlight = F, fill_highlight = "lightyellow", fill_alpha=0.5, upper_lim = NULL){
  
  # coordinate files
  if (genome == "hg38"){
    coord = setNames(read.table("ATACseq/bed/hg38_DHS.bed"), c("seqnames", "start", "end", "region_id"))
  }
  if (genome == "macFas6"){
    coord = setNames(read.table("ATACseq/bed/macFas6_DHS.bed"), c("seqnames", "start", "end", "region_id"))
  }
  
  coord_DHS<-coord %>% filter(region_id == get_region_id) %>%
    mutate(start = start - padding, end = end + padding)
  
  # get bigwig (can be easily made from the peak files in E-MTAB-13373)
  bigWigFile= paste0("/data/share/htp/ATACseq/Novogene_ATAC/02_Preprocessing/bwa_mapping/ATAC_b1_sample",bw_sample, "_", genome, ".bam_slidingWindow1.bw")
  options(ucscChromosomeNames=TRUE)
  atac_1<-rtracklayer::import(bigWigFile, 
                              which= as_granges(coord_DHS),
                              as="NumericList") %>%
    data.frame() %>%
    mutate(position=0+coord_DHS$start[1]:(length(value)-1+coord_DHS$start[1]))
  
  if (!is.null(upper_lim)){
    upper = upper_lim
  } else {
    upper = max(atac_1$value)+10
  }
  
  # plot
  if (highlight){
    coord_DHS_orig<-coord %>% filter(region_id == get_region_id) 
    hl<-geom_rect(coord_DHS_orig, 
                  mapping=aes(xmin = start, 
                              xmax = end, 
                              ymin = min(atac_1$value), 
                              ymax = upper), 
                  fill=fill_highlight, alpha = fill_alpha)
  } else {
    hl<-NULL
  }
  
  ggplot(atac_1) +
    hl +
    geom_ribbon(aes(x = position, y = value, ymin = min(value), ymax = value), 
                fill = fill_col, alpha = 0.8, position = "identity", na.rm = T) +
    geom_line(aes(x = position, y = value, group = 1), na.rm = T, lwd = 0.01) +
    theme_light(base_size = 9) +
    ylab("Normalized\ncoverage") +
    xlab(paste0("chr",coord_DHS$seqnames[1])) +
    scale_x_continuous(breaks = scales::pretty_breaks(n=3),
                       limits = c(min(atac_1$position),
                                  max(atac_1$position)), 
                       expand = c(0,0))+
    scale_y_continuous(limits = c(min(atac_1$value), 
                                  upper), 
                       expand = c(0,0))+
    theme(panel.grid = element_blank())
}





# label tissues for the cbust summary ###
label_tissue <- function(df,dnase_jamm_long){
  cbust_subset_list <- list()
  for (i in c("adrenal_gland","brain","thymus","stomach","heart","muscle","kidney","large_intestine","lung")) {
    id_subset <- dnase_jamm_long %>% filter(tissue==i) %>% dplyr::select(region_id, total) %>% distinct # subset regions by tissue
    cbust_subset_list[[i]] <- df %>% inner_join(id_subset) %>% dplyr::mutate(tissue=as.factor(i)) %>% distinct
  }
  # bind the tissue subset df by rows
  cbust_tissue_labelled <- bind_rows(cbust_subset_list) %>%
    dplyr::rename(Tissue=tissue) %>%
    dplyr::mutate(Tissue=fct_relevel(Tissue,"adrenal_gland","brain","heart","kidney","large_intestine","lung","muscle","thymus","stomach"))
  
  return(cbust_tissue_labelled)
}








get_chip<-function(name, jamm=jamm_inATAC, dnase=DHS_genes, exprtf=expressedTFs,
                   metadata_filter = "neur", old=TRUE, neuro_topGO_filter=NULL){
  if (old == TRUE){
    chip<-read.table(paste0("ATACseq/chipseq_gtrd/macs2 peaks_Homo sapiens.interval_target",name,"_allTFs")) 
  } else {
    chip<-read.table(paste0("ATACseq/chipseq_gtrd/macs2 peaks_Homo sapiens.interval_",name)) 
  }
  
  chip<-chip %>%
    dplyr::rename(seqnames=V1, start=V2, end=V3, TF=V4) %>%
    mutate(seqnames=gsub("chr","", seqnames),
           id=stringr::word(V11,2,2,"[.]")) %>%
    left_join(metdat %>% dplyr::select(id, title, specie)) %>%
    filter(specie=="Homo sapiens", grepl(metadata_filter,title)) %>%
    as_granges()
  
  jamm_chip<-find_overlaps(as_granges(jamm), chip) %>% 
    as_tibble %>%
    mutate(region_id = as.factor(region_id)) %>%
    inner_join(dnase %>% mutate(region_id = as.factor(region_id))) %>% 
    distinct(region_id = as.factor(region_id),TF, .keep_all=T) 
  
  if (!is_null(neuro_topGO_filter)){
    jamm_chip<-jamm_chip %>% filter(TF %in% neuro_topGO_filter) 
  }
  
  jamm_chip<-jamm_chip %>% 
    #select the ones we have a PWM for
    inner_join(exprtf %>% dplyr::select(motif_ID, motif_name), by=c("TF"="motif_name"))
    
  
  to_check<-jamm_chip %>%
    group_by(region_id, total) %>%
    dplyr::summarise(distance=min(distance),
                     motifs=paste(motif_ID, collapse = "|"),
                     motif_names=paste(TF, collapse = "|")) %>%
    arrange(distance)
  
  return(to_check)
}





# put it together with the peak stuff - make an adjusted peak function for this

plot_targeted_forFig<-function(bw = "/data/share/htp/HERVH/ESRG/primates.phastCons46way.bw", 
                               regid, PD = NULL, 
                               coord_file_4phastCons, 
                               locfile_4cbust, 
                               padding_bp = 30,
                               conc_TFs = NULL, 
                               motif_filter = NULL,
                               n_TFs = 5,
                               TF_score_cutoff = 5, 
                               cluster_score_cutoff = 0,
                               highlight1 = "white",
                               highlight1_alpha=0.5,
                               plot_labels=NULL){
  
  #"#D1D1DC","#8B8BA7"
  print("Plotting ATAC-seq data")
  pt<-plot_ATAC(genome = "hg38", bw_sample = "05", get_region_id = regid, 
                fill_col = "#8B8BA7", highlight = T, padding = padding_bp, 
                fill_highlight = highlight1, fill_alpha = highlight1_alpha)+
    annotate("text",x = -Inf, y=Inf, label="Human NPCs", hjust = -0.15, vjust = 1.5, size=2.5) +
    theme(panel.border = element_rect(colour = "grey20", linewidth = rel(0.7)),
          axis.title.x = element_text(size = 8))
  
  pb<-plot_ATAC(genome = "macFas6", bw_sample = "07", get_region_id = regid, 
                fill_col = "#D1D1DC", highlight = T, padding = padding_bp, 
                fill_highlight = highlight1, fill_alpha = highlight1_alpha)+#, upper_lim = 500) +
    annotate("text",x = -Inf, y=Inf, label="Macaque NPCs", hjust = -0.15, vjust = 1.5, size=2.5)+
    #ggtitle("Macaque iPSC") +
    theme(panel.border = element_rect(colour = "grey20", linewidth = rel(0.7)),
          axis.title.x = element_text(size = 8))
  
  print("Plotting phastCons")
  pph<-plot_phastCons(bigWigFile = bw, coord = coord_file_4phastCons,
                      get_region_id = regid, padding = padding_bp, 
                      highlight = T, fill_highlight = highlight1, fill_alpha = highlight1_alpha)+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(colour = "grey20", linewidth = rel(0.7)))
  
  print("Plotting TFBS coverage")
  outfile<-get_TFBS_positions_ATAC(get_region_id = regid, 
                                   human_file = locfile_4cbust$human_file[locfile_4cbust$region_id == regid], 
                                   macaque_file = locfile_4cbust$macaque_file[locfile_4cbust$region_id == regid],
                                   TF_score_cutoff = TF_score_cutoff, 
                                   cluster_score_cutoff = cluster_score_cutoff)
  
  if (!is.null(motif_filter)){
    outfile$binding_both<-outfile$binding_both %>% filter(cluster_id_or_motif_name %in% motif_filter)
  }
  
  print("Plotting TFBS sites")
  pbinding3<-plot_bindingSites_alt_noStrand(outfile = outfile, 
                                            get_region_id = regid, 
                                            padding = padding_bp,
                                            TFs = conc_TFs,
                                            bindingCols = c("#A5DDEB","#DDCC68"))+
    annotate("text",x = -Inf, y=8, label="Human", hjust = -0.15, #vjust = 1.5, 
             size=2.5)+
    annotate("text",x = -Inf, y=-8, label="Macaque", hjust = -0.15, #vjust = 28, 
             size=2.5)
  
  plot_grid(pt, pb, 
            pbinding3,
            rel_heights=c(0.7, 0.7, 1),
            ncol=1, align = "v", axis="lr", labels=plot_labels, label_y = 1.03, scale = 0.95, label_x=c(0,0,0.01), label_size = 12)
  
}




plot_ranks<-function(r, lab, tab = seq_vs_TFBS_rank, DE_tab = DE, plot_labels=NULL){
  
  rank_phc<-ggplot(tab, aes(x=rankPhastCons, y=meanPhastCons, group=1))+
    geom_line()+
    geom_area(fill="#DEB8A8")+
    geom_vline(xintercept = tab$rankPhastCons[tab$region_id==r])+
    geom_label(data=subset(tab, tab$region_id==r),
               aes(rankPhastCons,0.45, label=lab), size=2.7)+
    xlab("CRE rank")+
    ylab("Sequence\nconservation")+
    theme_bw(base_size = 9)+
    theme(panel.grid = element_blank(),
          axis.title.x=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  
  
  rank_cgi<-ggplot(tab, aes(x=rankCpG, y=CpG_obs_exp_cons, group=1))+
    geom_line()+
    geom_area(fill="#DEA9A9")+
    geom_vline(xintercept = tab$rankCpG[tab$region_id==r])+
    geom_label(data=subset(tab, tab$region_id==r),
               aes(rankCpG,0.45, label=lab), size=2.7)+
    xlab("CRE rank")+
    ylab("CpG exp/obs\nconservation")+
    theme_bw(base_size = 9)+
    theme(panel.grid = element_blank(),
          axis.title.x=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  
  rank_binding<-ggplot(tab, aes(x=rankBinding, y=meanBindingAgreement, group=1))+
    geom_line()+
    geom_area(fill="#DA8892")+
    geom_vline(xintercept = tab$rankBinding[tab$region_id==r])+
    geom_label(data=subset(tab, tab$region_id==r),
               aes(rankBinding,0.75, label=lab), size=2.7)+
    xlab("CRE rank")+
    ylab("TFBS binding site\nconservation")+
    theme_bw(base_size = 9)+
    theme(panel.grid = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none")

  rank_sh<-ggplot(tab, aes(x=rankCanberra, y=canb_cons.exprTissues, group=1))+
    geom_line()+
    geom_area(fill="#90475A")+ 
    geom_vline(xintercept = tab$rankCanberra[tab$region_id==r])+
    geom_label(data=subset(tab, tab$region_id==r),
               aes(rankCanberra,0.75, label=lab), size=2.7)+
    xlab("CRE rank")+
    ylab("TFBS repertoire\nconservation")+
    theme_bw(base_size = 9)+
    theme(panel.grid = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none")

  # openness
  rank_DA<-ggplot(tab, aes(x=rankDA, y=abs(LFC.DA), group=1))+
    geom_line()+
    geom_area(fill="#457b9d")+
    geom_vline(xintercept = tab$rankDA[tab$region_id==r])+
    geom_label(data=subset(tab, tab$region_id==r),
               aes(rankDA, 7.5, label=lab), size=2.7)+
    xlab("CRE rank")+
    ylab("Openness\n|LFC|")+
    theme_bw(base_size = 9)+
    theme(panel.grid = element_blank(),
          plot.margin = unit(c(5.5, 7, 5.5, 5.8), "points"))
  
  # gene expression
  expr_rank<-DE_tab %>% 
    distinct(gene_id, gene_name, LFC.DE = log2FoldChange) %>%
    arrange(abs(LFC.DE)) %>%
    mutate(rankDE = 1:nrow(.))
  
  
  rank_DE<-ggplot(expr_rank, aes(x=rankDE, y=abs(LFC.DE), group=1))+
    geom_line()+
    geom_area(fill="#1d3557")+
    geom_vline(xintercept = expr_rank$rankDE[expr_rank$gene_name==lab])+
    geom_label(data=subset(expr_rank, expr_rank$gene_name==lab),
               aes(rankDE, 6, label=lab), size=2.7)+
    xlab("Gene rank")+
    ylab("Expression\n|LFC|")+
    theme_bw(base_size = 9)+
    #basic_theme_ins+
    theme(panel.grid = element_blank(),
          plot.margin = unit(c(5.5, 7, 5.5, 5.8), "points"))
  
  
  all_ranks<-plot_grid(rank_phc, rank_binding, rank_sh, rank_DA, rank_DE, 
                       ncol=1, align = "v", axis="lr",
                       rel_heights = c(0.88, 0.88, 0.88, 1, 1),
                       labels = plot_labels, scale = 0.95, label_y = 1.04, label_size = 12)
  
}






