# functions


find_overlapping_region_ids<-function(ranges_obj){
  
  inner_overlaps<-findOverlaps(ranges_obj, minoverlap = 1) #adjacent peaks
  inner_overlaps_df<-inner_overlaps %>% 
    as.data.frame %>% 
    filter(queryHits!=subjectHits) %>% 
    mutate(queryHits=as.character(queryHits)) %>%
    left_join(as_tibble(ranges_obj) %>% 
                rownames_to_column("queryHits"))
  
}


reduce_liftover<-function(LO_ranges, extend = 50){
  
  #first extend all to downstream by 10 bases and see if they overlap now. 
  #then after reduction, remove 10b from the end (this is because large regions sometimes get split up)
  reduced_LO_ranges<-LO_ranges %>%
    as_tibble() %>%
    mutate(end=end+extend) %>%
    dplyr::select(-width) %>%
    as_granges() %>%
    group_by(region_id) %>%
    reduce_ranges() %>%
    as_tibble() %>%
    mutate(end=end-extend) %>%
    #now select only matches in 1 piece
    group_by(region_id) %>%
    mutate(n=length(region_id)) %>%
    ungroup() %>%
    filter(n==1) %>%
    dplyr::select(-width, -n)
  
  return(reduced_LO_ranges)
}



#get the orthologous coordinates
translate_jamm<-function(chain_file, coordinate_file, extend=50,
                         reverse_chain_file, secondary = NULL){
  
  require(rtracklayer)
  
  chain <- import.chain(chain_file)
  
  jamm_forCbust_cross<-coordinate_file %>%
    liftOver(chain) %>% 
    unlist() %>%
    reduce_liftover(extend = extend)
  
  #check if it lands back where it should in the human (reciprocal part)
  chainBack <- import.chain(reverse_chain_file)
  
  jamm_forCbust_return<-as_granges(jamm_forCbust_cross) %>%
    liftOver(chainBack) %>% 
    unlist() %>%
    reduce_liftover(extend = extend)
  
  if (!is.null(secondary)){
    
    #this is typically hg38 to hg19
    chainBackMore <- import.chain(secondary)
    
    jamm_forCbust_return2<-as_granges(jamm_forCbust_return) %>%
      liftOver(chainBackMore) %>% 
      unlist() %>%
      reduce_liftover(extend = extend)
    
    jamm_forCbust_return<-jamm_forCbust_return2
  }
  
  #now check that it overlaps with the same region as it came from
  coord_OL<-find_overlaps(as_granges(coordinate_file),
                          as_granges(jamm_forCbust_return)) %>%
    filter(region_id.x==region_id.y)
  
  out<-jamm_forCbust_cross %>% dplyr::filter(region_id %in% coord_OL$region_id.x) %>%
    #add width
    as_granges() %>%
    as_tibble()
  
  #return the coordinates of the reciprocal matches in the OTHER species
  return(out)
}



# save the selected PFMs
savePFM<-function(PFMatrixList_obj, out.file, append_PFM1=NULL, name_PFM1=NULL, append_PFM2=NULL, name_PFM2=NULL){
  
  system(paste0("rm ",out.file)) #since i'm appending, first need to remove the old version
  
  lapply(PFMatrixList_obj, function(x) {
    pfm_t<-t(x@profileMatrix)
    paste0(">",x@ID,"_",x@name)
    write(paste0(">",x@ID,"_",x@name),file=out.file,append=TRUE)
    write.table( data.frame(pfm_t), out.file, quote=F, row.names = F, col.names = F, append= T, sep=" ")
  })
  
  if (!is.null(append_PFM1)){
    write(paste0(">",name_PFM1),file=out.file,append=TRUE)
    write.table( data.frame(append_PFM1), out.file, quote=F, row.names = F, col.names = F, append= T, sep=" ")
  }
  
  if (!is.null(append_PFM2)){
    write(paste0(">",name_PFM2),file=out.file,append=TRUE)
    write.table( data.frame(append_PFM2), out.file, quote=F, row.names = F, col.names = F, append= T, sep=" ")
  }
}




subset_chipseq_byDNaseReg<-function(tissue_vec=tissue4_vec, remap = remap_tissue_df, 
                                    cbust = cbust_regions, jaspar = jaspar_ic_df){
  
  tissue_tables<-list()

  for (tis in tissue_vec){
  
    remap_tissue_TF<-remap %>%
      filter(tissue==tis) %>%
      filter(region_id %in% cbust$region_id[cbust$tissue==tis]) %>%
      filter(TF %in% jaspar$TF) %>%
      mutate(given=1) %>%
      distinct() %>%
      tidyr::pivot_wider(id_cols = c("tissue","region_id"), names_from = "TF", 
                  values_from="given", values_fill = 0) %>%
      mutate(at_least_1Chip=1) %>%
      right_join(cbust %>% filter(tissue==tis)) %>%
      filter(at_least_1Chip==1)
  
    remap_tissue_TF[is.na(remap_tissue_TF)]<-0
  
    remap_tissue_TF_long<-remap_tissue_TF %>%
      pivot_longer(cols=c(-tissue, -region_id, -at_least_1Chip, -CpG_obs_exp, -GC_content, -total),
                   names_to="TF", values_to = "Class") %>%
      inner_join(jaspar) 
  
    tissue_tables[[tis]]<-remap_tissue_TF_long
  }
  tissue_tables_combined<-bind_rows(tissue_tables, .id="tissue")
  return(tissue_tables_combined)
}




## functions for motif density etc
#library(MASS, lib.loc = "/opt/R/4.1.0/lib/R/library")

fit_boxcox_anova<-function(df, to_predict, var, predictors, group_var=total,  tukey_predictor="total", title_plot="TFBS density", y_lab="Box-Cox[number of\nTFBSs per bp]", x_lab="Number of tissues"){
  
  model<-MASS::boxcox(as.formula(paste(to_predict, "~", paste(predictors, collapse="+"))), data=df) 
  lambda <- model$x[which.max(model$y)]
  # 2.2 ANOVA
  ANOVA <- aov(as.formula(paste("(((", to_predict, ")^lambda-1)/lambda) ~", paste(predictors, collapse="+"))), data = df)

  plot(ANOVA,1) # diagnostics
  plot(ANOVA,2)
  plot(ANOVA,3)
  summary(ANOVA)
  #write.csv(tidy(ANOVA) ,"stat/roadmap_DHS_summaries/tfbs_dens_ANOVA.csv")
  # 2.3 Post-hoc Tukey test
  tukey <- TukeyHSD(ANOVA, which = tukey_predictor)
  
  value_max <- df  %>% dplyr::group_by({{ group_var }}) %>% 
    dplyr::summarise(max_value = max(((({{ var }})^lambda-1)/lambda)))
  

  hsd <- agricolae::HSD.test(ANOVA, trt = tukey_predictor, group = T)
  sig.letters <- hsd$groups[order(row.names(hsd$groups)), ]
  df<-df %>%
    mutate(lambda_norm=((({{ var }})^lambda-1)/lambda))
  
  fig_tfbs <- df %>% 
    #filter(motif_dens>0) %>%
    ggplot(aes(x=total, y= ((({{ var }})^lambda-1)/lambda))) +
    geom_boxplot(notch=T,lwd=1,position = "dodge",aes(fill=species)) +
    geom_text(data = value_max, aes(x=total, y = 0.1 + max_value, label = sig.letters$groups), vjust=-1, size=6)+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+ #0.1
    viridis::scale_fill_viridis(discrete = T, end = 0.4, alpha=0.7,direction=1) +
    labs(title = title_plot)+
    theme_light(base_size =20) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(colour = "grey", size = 0.5),
          panel.border = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.direction = "horizontal",
          plot.title = element_text(face="bold",hjust = 0.5)) +
    ylab(y_lab) + #\u03bb check the lamda
    xlab(x_lab)
  return(list(fig_tfbs=fig_tfbs, df=df))
}




fastaToCpG<-function(path){
  
  seq<-readDNAStringSet(path)
  
  seq_df<-data.frame(region_id = names(seq),
                      seq = seq) %>%
    rowwise() %>%
    mutate(GC_content = calc_gc(DNAString(seq)),
           CpG_obs_exp = calc_cpg_oe(DNAString(seq))) %>%
    dplyr::select(-seq)
  
  return(seq_df)
}

