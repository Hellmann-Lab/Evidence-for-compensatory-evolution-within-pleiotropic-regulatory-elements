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





extract_region_fa<-function(genome_path, coord_file, species){
  
  require(Biostrings)
  gen<-readDNAStringSet(genome_path)
  chrom_sizes<-data.frame(seqnames=names(gen), chrom_width=width(gen))
  
  coord_500bExt<- coord_file %>%
    dplyr::filter(similar=="yes") %>%
    dplyr::left_join(chrom_sizes) %>%
    dplyr::mutate(start_up=ifelse(start-500>=0,start-500,0),
           end_down=ifelse(end+500<=chrom_width,end+500,chrom_width)) %>%
    dplyr::select(seqnames,start_up,end_down,region_id) %>%
    dplyr::rename(start=start_up, end=end_down)
  
  #and cut out the sequences 
  registerDoParallel(30)  # use multicore, set to the number of cores
  foreach (seq=unique(coord_500bExt$seqnames)) %dopar% {
    to_cut_coord<-coord_500bExt %>% filter(seqnames==seq)
    to_cut_fa<-gen[names(gen)==seq]
    
    #very annoying, subseq did not work on all seqs at once so i'm looping
    dhs_string<-DNAStringSet()
    for (i in 1:nrow(to_cut_coord)){
      subs<-subseq(to_cut_fa,start=to_cut_coord$start[i], end=to_cut_coord$end[i])
      names(subs)<-to_cut_coord$region_id[i]
      dhs_string<-c(dhs_string,subs)
    }
    dir.create(file.path(paste0("TFBS/cbust_on_paddedSeq/fastas/",species)), showWarnings = FALSE)
    writeXStringSet(dhs_string,paste0("TFBS/cbust_on_paddedSeq/fastas/",species,"/dhs_",seq,".fa"))
  }
}


extract_region_fa_noExtension<-function(genome_path, coord_file, species){
  
  gen<-readDNAStringSet(genome_path)
  
  coord_500bExt<- coord_file %>%
    filter(similar=="yes") 
  
  #and cut out the sequences 
  registerDoParallel(25)  # use multicore, set to the number of cores
  foreach (seq=unique(coord_500bExt$seqnames)) %dopar% {
    to_cut_coord<-coord_500bExt %>% filter(seqnames==seq)
    to_cut_fa<-gen[names(gen)==seq]
    
    #very annoying, subseq did not work on all seqs at once so i'm looping
    dhs_string<-DNAStringSet()
    for (i in 1:nrow(to_cut_coord)){
      subs<-subseq(to_cut_fa,start=to_cut_coord$start[i], end=to_cut_coord$end[i])
      names(subs)<-to_cut_coord$region_id[i]
      dhs_string<-c(dhs_string,subs)
    }
    writeXStringSet(dhs_string,paste0("TFBS/cbust_on_paddedSeq/fastas_short/",species,"/dhs_",seq,".fa"))
  }
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



#for prediction evaluation

predict_binding <- function(model, data, observation_col){
  
  test_set <- data.frame(obs = data %>% pull({{observation_col}}),
                         Class_0 = predict(model,data, type = "prob")[,"Class_0"],
                         Class_1 = predict(model,data, type = "prob")[,"Class_1"])
  test_set$pred <- factor(ifelse(test_set$Class_1 >= .5, "Class_1", "Class_0"), 
                          levels=c("Class_1", "Class_0"))
  test_set$obs<-factor(test_set$obs, levels=c("Class_1","Class_0"))
  return(test_set)
  
}


evaluate_prediction<-function(test_set){
  
  predobs_df<-data.frame(pred_0=length(test_set$pred[test_set$pred==0 | test_set$pred=="Class_0"]),
                            pred_1=length(test_set$pred[test_set$pred==1 | test_set$pred=="Class_1"]),
                            obs_0=length(test_set$obs[test_set$obs==0 | test_set$obs=="Class_0"]),
                            obs_1=length(test_set$obs[test_set$obs==1 | test_set$obs=="Class_1"]))
  postResample_df <- as.data.frame(t(postResample(pred = test_set$pred, obs = test_set$obs)))
  twoClassSummary_df <- as.data.frame(t(twoClassSummary(test_set, lev=levels(test_set$obs))))
  #prSummary_df <- as.data.frame(t(prSummary(test_set, lev = levels(test_set$obs)))) 
  perf <- cbind(predobs_df,postResample_df,twoClassSummary_df)
  return(perf)
  
}




wrap_evaluate_prediction_automl<-function(outlist, bootstrap=F, seed=NULL, frac){
  
  predictions_automl<-lapply(outlist, function(x){
    print(x)
    out<-read_csv(x) %>% setNames(c("obs", "pred", "pred_prob")) 
    
    if (bootstrap){
      set.seed(as.numeric(seed))
      out<-out %>% dplyr::sample_frac(size = frac)
    }
    
    out %>%
      dplyr::transmute(obs=factor(obs, levels=c(1,0)),
                       pred=factor(pred, levels=c(1,0))) %>%
      as.data.frame() %>%
      evaluate_prediction() %>%
      mutate(ROC=as.numeric(pROC::auc(out$obs, out$pred_prob)))
  })
  predictions_automl_df<-bind_rows(predictions_automl)
  return(predictions_automl_df)
  
}




calc_thresh<-function(outlist, get_max = T, get_similar = F, bootstrap=F, 
                      seed=NULL, upper_lim=NULL, lower_lim=NULL, plotdir=NULL){
  
  predictions_automl<-lapply(outlist, function(x){
    print(x)
    out<-read_csv(x) %>% setNames(c("obs", "pred", "pred_prob")) 
    
    if (bootstrap){
      set.seed(as.numeric(seed))
      out<-out %>% dplyr::sample_frac(size = 0.8)
    }
    
    roc <- pROC::roc(out$obs, out$pred_prob, percent = F)
    
    perf<-data.frame(sensspec=(roc$sensitivities+roc$specificities)/2,
                     sens=roc$sensitivities, 
                     spec=roc$specificities,
                     n=1:length(roc$sensitivities),
                     thres=roc$thresholds) %>% 
      mutate(sp = stringr::str_extract(x[1], "[[:digit:]]+")) 
    
    perf_long<-perf %>%
      pivot_longer(cols = c("sensspec","sens","spec"), names_to = "metric", values_to = "sensspec")
    
    if (get_max){
      sel<-perf[perf$sensspec==max(perf$sensspec),]
      #sel<-perf[perf$sens==max(perf$sens),]
      
    } 
    
    if (get_similar){
      sel<-perf[(perf$sensspec>=lower_lim & perf$sensspec<=upper_lim),]
      sel<-sel[sel$sens==max(sel$sens),]
      sel<-sel[sel$spec==max(sel$spec),]
      #cand<-perf[(perf$sens>=lower_lim & perf$sens<=upper_lim),]
      #sel<-cand[cand$spec==max(cand$spec),]
      
      pl<-ggplot(perf_long, 
                   aes(x=thres, y=sensspec, group=metric, color=metric))+
              geom_line()+
              geom_vline(xintercept = sel$thres)+
              geom_hline(yintercept = upper_lim)+
              ggtitle(paste("Number of tissues:", sel$sp))+
        theme_bw()+
        xlab("Threshold")+
        ylab("1 - Error")
      print(pl)
      ggsave(filename = paste0(plotdir,"/error_",sel$sp,".pdf"),plot=pl,device="pdf", height=3.5, width=5.3)
    }
    sel
  })
  predictions_automl_df<-bind_rows(predictions_automl)
  return(predictions_automl_df)
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
  #write.csv(tidy(ANOVA) ,"stat/general/tfbs_dens_ANOVA.csv")
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




# label tissues for the cbust summary ###
label_tissue <- function(df,dnase_jamm_long){
  cbust_subset_list <- list()
  for (i in c("adrenal_gland","brain","thymus","stomach","heart","muscle","kidney","large_intestine","lung")) {
    id_subset <- dnase_jamm_long %>% filter(tissue==i) %>% dplyr::select(region_id) %>% distinct # subset regions by tissue
    cbust_subset_list[[i]] <- df %>% inner_join(id_subset) %>% filter(!is.na(total)) %>% dplyr::mutate(tissue=as.factor(i)) %>% distinct
  }
  # bind the tissue subset df by rows
  cbust_tissue_labelled <- bind_rows(cbust_subset_list) %>%
    dplyr::rename(Tissue=tissue) %>%
    dplyr::mutate(Tissue=fct_relevel(Tissue,"adrenal_gland","brain","heart","kidney","large_intestine","lung","muscle","thymus","stomach"))
  
  return(cbust_tissue_labelled)
}




get_aln_stats<-function(seqs1, seqs2, aln_mafft){
  
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
                          tmp@insertion[1,2] + tmp@deletion[1,2] }),
                        mismatches = sapply(all_aln, function(x){
                          tmp<-sum(mismatchSummary(x)$pattern$position$Count) })) #,
  #dnds = sapply(all_aln, function(x){
  # dnds(as.DNAbin(c(aligned(subject(x)), aligned(pattern(x))))) })
  
  return(list(all_aln, aln_stats))
}




calc_cpg_oe <- function(seq){
  ff<- Biostrings::alphabetFrequency(seq)
  l <- sum(ff[c("A","C","G","T")])
  eff.length <- l - (sum(ff)-l +1 )
  Biostrings::dinucleotideFrequency(seq)["CG"] / ( ff["G"]* ff["C"]) * l^2 / eff.length
}


calc_gc<-function(seq){
  ff<- Biostrings::alphabetFrequency(seq)
  gc <- sum(ff[c("C","G")])
  l <- sum(ff[c("A","C","G","T")])
  gc/l
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

