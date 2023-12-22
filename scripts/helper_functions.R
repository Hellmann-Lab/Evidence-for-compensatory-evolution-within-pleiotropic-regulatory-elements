# makeTSSgr -----------------------------------
#' @name makeTSSgr
#' @description Uses a txdb object to create a genomic ranges object with transcription start sites
#' @usage makeTSSgr( gene_id, upstream =2000, downstream = 2000, txdbFile)
#' @param gene_id vector of gene ids for which to find TSS gene_ids need to be identical to the ones used in txdb
#' @param upstream how much sequence upstream of the transcript start to use.
#' @param downstream how much sequence downstream of the transcript start to use.
#' @param txdbFile location of the file with the txdb
#' @return A GRanges object with transcript ids as rownames and gene_id in mcols.
#' @importFrom AnnotationDbi loadDb select
#' @importFrom GenomicRanges GRanges
#' @export
makeTSSgr <- function(gene_id , upstream = 2000, downstream = 2000 , txdbFile){
  
  txdb<- AnnotationDbi::loadDb(txdbFile)
  g2T <- AnnotationDbi::select(txdb, keys=bt$ensembl_gene_id,
                               columns = c("GENEID","TXNAME","TXSTRAND","TXCHROM","TXSTART","TXSTRAND","TXEND"),
                               keytype = "GENEID" )

  GRanges( seqnames = g2T$TXCHROM,
           ranges =   IRanges( start = ifelse( g2T$TXSTRAND == "+", 
                                               ifelse(g2T$TXSTART - upstream >0, g2T$TXSTART-upstream, 0), 
                                               ifelse(g2T$TXEND - downstream >0 , g2T$TXEND-downstream ,  0) ),
                               end =  ifelse( g2T$TXSTRAND == "+", 
                                              g2T$TXSTART + downstream,
                                              g2T$TXEND + upstream ),
                               names = g2T$TXNAME),
           strand = g2T$TXSTRAND,
           gene_id = g2T$GENEID )
  
}



#associate dhs to genes

require(tidyverse)
require(plyranges)


find_elements_left_right <- function(x, y , overlaps, dist) {
  left_of_peak <- x %>% filter( !(peak_id %in% overlaps$peak_id) ) %>% 
    join_precede( y ) %>% 
    plyranges::mutate( distance  = tss.start - end,
                       direction = ifelse(tss.strand == "+", "upstream","downstream")) %>% 
    filter(distance< dist )
  
  right_of_peak<- x %>% filter( !(peak_id %in% overlaps$peak_id) ) %>% 
    join_follow( y ) %>% 
    plyranges::mutate( distance  = start - tss.start,
                       direction = ifelse(tss.strand == "-", "upstream","downstream")) %>% 
    filter(distance< dist )
  c(left_of_peak,right_of_peak)
}

annot_peaks<- function(x,y, prom_dist=2000, enh_dist=1e6, gtf_file, rel_dist=10,
                       expr_filter_var_y = F){
  
  # rel_dist only keeps left and right peaks if one is not rel_dist times closer set to 1e6 for no filtering
  
  y<- y %>%  plyranges::mutate( tss.start = start,
                                tss.strand = strand ) 
  
  if( expr_filter_var_y == F){
    print("Evaluate all TSS in y for enhancer annotation.")
    expr_filter_var_y<-"nofilter"
    y$nofilter<- T
  }else if (expr_filter_var_y == T){
    print("Only evaluate TSS for enhancer annotation that overlap with a peak.")
    expr_filter_var_y<-"nofilter"
    y$nofilter <- F
  }else{
    print(paste("Evaluate TSS for enhancer annotation that overlap with a peak or are marked in", expr_filter_var_y))
  }
  
  #find Promoters TSS either falls within the peak or
  overlaps<- x %>% 
    join_overlap_inner( y ) %>% 
    plyranges::mutate( distance = 0,
                       direction = "overlap")
  
  # and the next upstream or downstream TSS if not further away than prom_dist is used
  lr_prom<- find_elements_left_right( x, y, overlaps, prom_dist)
  
  promoter_peaks<- c(overlaps,lr_prom) %>% plyranges::mutate(element_type = "prom")
  
  #### Enhancer consider expressed genes & peaks as putatively expressed
  # remove all promoter peaks each peak will only have at most 2 genes associated
  
  enhancer_peaks <- find_elements_left_right( x = x %>% filter( !( peak_id %in% unique(promoter_peaks$peak_id)) ),
                                              y = y %>% filter(  (  tss_id %in% unique(promoter_peaks$tss_id)) |
                                                                   !!as.name( expr_filter_var_y ) ),
                                              overlaps = overlaps,
                                              dist = enh_dist ) %>% 
    plyranges::mutate(element_type = "enh") %>% 
    group_by(peak_id) %>% 
    filter( distance < 10*min(distance))
  
  #summarise promoter gene associations
  # for ensembl gtfs could add gene_type to grouping variable to 
  allPeaks<-bind_ranges(promoter_peaks, enhancer_peaks) %>% 
    as_tibble %>% 
    dplyr::group_by(seqnames,start,end,peak_id, gene_id, 
                    gene_name,tss.strand,element_type) %>% 
    dplyr::summarise(distance = min(distance),
                     tag      = paste(unique(tag),collapse=","),
                     n_tss    = length(unique(tss_id)),
                     direction = paste(unique(direction),collapse=",")) %>% 
    as_granges
  
  if(!missing(gtf_file)){
    exons<- read_gff(gtf_file) %>% filter(type == "exon") 
    allPeaks <- allPeaks %>% plyranges::mutate(exonOL = count_overlaps(., exons))
  }
  return(allPeaks)
}






annot_peaks2<- function(x,y, prom_dist=2000, enh_dist=1e6, gtf_file, rel_dist=10 ,
                       expr_filter_var_y = F){
  
  # rel_dist only keeps left and right peaks if one is not rel_dist times closer set to 1e6 for no filtering
  
  y<- y %>%  plyranges::mutate( tss.start = start,
                                tss.strand = strand ) 
  
  if( expr_filter_var_y == F){
    print("Evaluate all TSS in y for enhancer annotation.")
    expr_filter_var_y<-"nofilter"
    y$nofilter<- T
  }else if (expr_filter_var_y == T){
    print("Only evaluate TSS for enhancer annotation that overlap with a peak.")
    expr_filter_var_y<-"nofilter"
    y$nofilter <- F
  }else{
    print(paste("Evaluate TSS for enhancer annotation that are marked in", expr_filter_var_y))
  }
  
  #find Promoters TSS either falls within the peak or
  overlaps<- x %>% 
    join_overlap_inner( y ) %>% 
    plyranges::mutate( distance = 0,
                       direction = "overlap")
  
  # and the next upstream or downstream TSS if not further away than prom_dist is used
  lr_prom<- find_elements_left_right( x, y, overlaps, prom_dist)
  
  promoter_peaks<- c(overlaps,lr_prom) %>% plyranges::mutate(element_type = "prom")
  
  #### Enhancer consider expressed genes & peaks as putatively expressed
  # remove all promoter peaks each peak will only have at most 2 genes associated
  
  enhancer_peaks <- find_elements_left_right( x = x %>% filter( !( peak_id %in% unique(promoter_peaks$peak_id)) ),
                                              # this !! selects expressed = T
                                              y = y %>% filter(  #(  tss_id %in% unique(promoter_peaks$tss_id)) |
                                                                   !!as.name( expr_filter_var_y ) ),
                                              overlaps = overlaps,
                                              dist = enh_dist ) %>% 
    plyranges::mutate(element_type = "enh") %>% 
    group_by(peak_id) %>% 
    filter( distance < 10*min(distance))
  
  #summarise promoter gene associations
  # for ensembl gtfs could add gene_type to grouping variable to 
  allPeaks<-bind_ranges(promoter_peaks, enhancer_peaks) %>% 
    as_tibble %>% 
    dplyr::group_by(seqnames,start,end,peak_id, gene_id, 
                    gene_name,tss.strand,element_type) %>% 
    dplyr::summarise(distance = min(distance),
                     tag      = paste(unique(tag),collapse=","),
                     n_tss    = length(unique(tss_id)),
                     direction = paste(unique(direction),collapse=",")) %>% 
    as_granges
  
  if(!missing(gtf_file)){
    exons<- read_gff(gtf_file) %>% filter(type == "exon") 
    allPeaks <- allPeaks %>% plyranges::mutate(exonOL = count_overlaps(., exons))
  }
  return(allPeaks)
}




#plot assignment numbers
summarise_assignments<-function(assignment_table, regids, tissuewise=F){
  
  if (tissuewise){
    
    hg_s<-assignment_table %>%
      as_tibble() %>% 
      group_by(region_id, tissue) %>%
      summarise(n_genes=length(unique(gene_id)),
                minDist=min(distance),
                element_type=element_type[1]) %>%
      group_by(region_id) %>%
      summarise(n_genes=max(n_genes),
                minDist=min(minDist),
                element_type=element_type) %>%
      distinct()
    
  } else {
    
  hg_s<-assignment_table %>%
    as_tibble() %>% 
    group_by(region_id) %>%
    summarise(n_genes=length(unique(gene_id)),
              minDist=min(distance),
              element_type=element_type[1])
  }
  
  unassigned<-unique(regids[!regids %in% hg_s$region_id])
  
  p1<-ggplot(hg_s, aes(x=n_genes, fill=element_type))+geom_bar(width=0.8)+
    #scale_y_log10()+
    theme_bw()+theme(legend.position = "top", legend.title = element_blank())+
    scale_fill_manual(values=c("#9CA578","#3288BD"))+xlab("# of associated genes")+ylab("# of CREs")
  
  p2_1<-ggplot(hg_s %>% filter(element_type=="Enhancer"), aes(x=minDist))+# / 1000))+
    geom_histogram(bins=50, fill="grey70")+
    scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3))+
    scale_x_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3))+
    facet_grid(.~element_type, scale="free")+
    theme_bw(base_size = 7)+
    xlab("Distance to gene")+ylab("# of CREs")
  
  p2_2<-ggplot(hg_s %>% filter(element_type=="Promoter"), aes(x=minDist))+# / 1000))+
    geom_histogram(bins=50, fill="grey70")+
    scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3))+
    scale_x_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3))+
    #scale_y_log10()+
    facet_grid(.~element_type, scale="free")+
    theme_bw(base_size = 7)+
    xlab("Distance to gene")+ylab("# of CREs")
    
    
  sp2<-ggplot(hg_s, aes(x=minDist))+geom_histogram(bins=50, fill="grey70")+
    scale_y_log10()+
    facet_grid(.~element_type, scale="free")+
    theme_bw()+
    ggtitle(paste0("# enhancer ",round(table(hg_s$element_type)[1]/1000,1),"k (",round(table(hg_s$element_type)[1]/length(regids)*100,1),"%), ",
                     "# promoter ", round(table(hg_s$element_type)[2]/1000,1),"k (",round(table(hg_s$element_type)[2]/length(regids)*100,1),"%), ",
                     "# unassigned: ", length(unassigned), " (",round(length(unassigned)/length(regids)*100,1), "%)"))+ 
    theme(plot.title = element_text(size = 12))+
    xlab("Distance to closest gene")+
    ylab("# of CREs")
  
  return(list(p1=p1,p2_1=p2_1,p2_2=p2_2,p2=p2))

}







filter_function <- function(data, var, value, op = `==`) {
  data %>% filter(op(.data[[var]], .env$value))
}

associate_per_tissue<-function(peak_GR=jammGR, tss_GR=tss, tss_GR_tag="expressed", count_table =summarized_expression_all, tissue, prom_or_expr = F){
  
  print(tissue)
  x = peak_GR %>% as_tibble() %>% filter_function(tissue, 1) %>% as_granges()
  #print("here")
  #tis<-deparse(substitute(var))
  y = tss_GR %>% plyranges::mutate(expressed = ifelse(gene_id %in% count_table$gene_id[count_table$tissue==tissue], T,F))
  #print("here")
  if (prom_or_expr){
    peaks_to_gene<-annot_peaks(x = x,
                                y = y,
                                prom_dist = 2000,
                                enh_dist = 1e6,
                                rel_dist = 10,
                                expr_filter_var_y = tss_GR_tag)
  } else {
  peaks_to_gene<-annot_peaks2(x = x,
                              y = y,
                              prom_dist = 2000,
                              enh_dist = 1e6,
                              rel_dist = 10,
                              expr_filter_var_y = tss_GR_tag)
  }
  return(peaks_to_gene)
}




########################
# SEQUENCE ANALYSIS ####
########################

get_insight<-function(dirpath, original=F){
  
  INSIGHT_res_colnames<-c("dataID", "thres", "rho", "rho_stderr", "E.A.", "E.A._stderr", "E.W.",
                          "E.W._stderr", "alpha", "alpha_stderr", "tau", "tau_stderr", "eta", 
                          "eta_stderr", "gamma", "gamma_stderr", "lnLd",
                          "LRT.rho>0.",  "LRT.eta>0.",  "LRT.gamma>0.", "em_status")
  
  if (original){
    
    result_list<-list.files(dirpath,full.names = T,recursive = T, pattern = "insight.results.txt")
    results <- setNames(lapply(result_list, function(x) read.table(x, header = T)) %>% bind_rows,
                        INSIGHT_res_colnames)
    
  } else {
    
    result_list<-list.files(dirpath,full.names = T,recursive = T, pattern = "ins.results.txt")
    results <- setNames(lapply(result_list, function(x) read.table(x)) %>% bind_rows,
                        INSIGHT_res_colnames) %>% distinct()
  }
  return(results)
}




pivot_longer_insight<-function(insight_wide){
  
  insight_long<-data.table::melt(setDT(insight_wide), 
                                        measure.vars = list(c("rho","E.A.", "E.W."), 
                                                            c("rho_stderr", "E.A._stderr", 
                                                              "E.W._stderr")),
                                        variable.name="selection_type", 
                                        value.name=c("selection_coeff", "stderr")) %>%
    dplyr::mutate(selection_type=case_when(selection_type=="1" ~ "rho",
                                    selection_type=="2" ~ "E.A.",
                                    selection_type=="3" ~ "E.W."),
           selection_coeff=as.numeric(selection_coeff),
           stderr=as.numeric(stderr),
           Tissues=as.factor(gsub("sp","",sp)))
  return(insight_long)
}





analysePhastCons <- function( bigWigFile, gr, probcut=0.9){
  phastCons  <- rtracklayer::import(bigWigFile, 
                                    which= gr,
                                    as="NumericList")
  sumP<- sapply(phastCons, function(x){ 
    n<-length(x)
    c(n, sum(x), sum(x>probcut) )
  }) %>% t() %>% data.frame()
  
  names(sumP)<- c("bp","sumPhastCons","neg_n")
  sumP <- as_tibble(phastCons@metadata$ranges) %>% dplyr::select(-strand) %>% bind_cols(sumP) %>% distinct()
  
  gr %>% as_tibble() %>% inner_join(sumP, by = c("seqnames", "start", "end", "width")) %>%
    dplyr::rowwise() %>% mutate(meanPhastCons = sumPhastCons/bp,
                                fracPhastCons = neg_n/bp)

}

# analysePhyloP <- function( bigWigFile, gr){
#   phyloP  <- rtracklayer::import(bigWigFile, 
#                                  which= gr,
#                                  as="NumericList")
#   sumP<- sapply(phyloP, function(x){ 
#     c(length(x), max(x), sum(x>=1.3) )
#   }) %>% t() %>% data.frame()
#   
#   names(sumP)<- c("bp","maxPhyloP","consPhylo")
#   gr %>% as_tibble() %>% bind_cols(sumP)
# }

analysePhyloP <- function( bigWigFile, gr, consCutoff = 1){
  phyloP  <- rtracklayer::import(bigWigFile, 
                                 which= gr,
                                 as="NumericList")
  sumP<- sapply(phyloP, function(x){ 
    c(length(x), max(x), sum(x), sum(x>=consCutoff), sum(x[x>=0]), length(x[x>=0]), sum(x[x<0]), length(x[x<0]) )
  }) %>% t() %>% data.frame()
  
  names(sumP)<- c("bp","maxPhyloP", "sumPhyloP", "consPhyloP_overCutoff", "sum_consPhyloP", "bp_consPhyloP", "sum_posPhyloP", "bp_posPhyloP")
  sumP <- as_tibble(phyloP@metadata$ranges) %>% dplyr::select(-strand) %>% bind_cols(sumP) %>% distinct()
  
  gr %>% as_tibble() %>%  inner_join(sumP, by = c("seqnames", "start", "end", "width")) %>%
    dplyr::rowwise() %>% mutate(meanPhyloP = sumPhyloP/bp,
                                fracPhyloP_overCutoff = consPhyloP_overCutoff/bp,
                                meanConsPhyloP = sum_consPhyloP/bp_consPhyloP,
                                meanConsPhyloP_relToAll = sum_consPhyloP/bp,
                                meanPosSelPhyloP = sum_posPhyloP/bp_posPhyloP,
                                meanPosSelPhyloP_relToAll = sum_posPhyloP/bp)
}




library(tidyverse)
library(lme4)
library(lmerTest)

permutation_fit<- function(df,shuffle=T) {
  
  if(shuffle){ 
    print("shuffeling")
    df2<- df %>% ungroup %>% group_by(assignment) %>%  mutate(total= sample(total), size=length(total)) 
  }else{
    df2<-df
    print("not shuffeling")
  }
  sumdf<-df2 %>% 
    group_by(tissue, gene_id, total, CGI, assignment, log2_mean_expression ) %>% 
    dplyr::summarise(nn = sum(1/log2(distance+2))) %>%  
    pivot_wider( names_from = assignment, values_from = nn) %>% 
    pivot_wider( names_from = CGI, values_from = c(enh,prom), names_sep = "_") %>% 
    pivot_wider( names_from = total, values_from = c(enh_CGI,enh_noCGI, prom_CGI,prom_noCGI ),  names_sep="_") %>% 
    replace(is.na(.), 0) 
  
  fixedVar<-glue::glue_collapse( colnames(sumdf[,-(1:3)]), sep=" + ") 
  f<- glue::glue(" log2_mean_expression ~ {fixedVar} + (1|tissue)")
  ff<- as.formula(f)
  mod <- lme4::lmer( ff, sumdf %>% ungroup)
  
  mr2 <- partR2::partR2(mod, data = sumdf, R2_type = "marginal", nboot = 10)
  
  return( mr2 )
}

analyse_permutations <- function(df, n=100){
  r2<-tibble()
  coef<-tibble()
  
  for(i in 1:n){
    print(i)
    f<-permutation_fit(df)
    r2<-bind_rows(r2,f$R2)
    coef<-bind_rows(coef,f$BW)
  }
  
  real<-permutation_fit(df,shuffle=F)
  
  r2_sum<- r2 %>% 
    dplyr::summarise( mean_r2 = mean(estimate),
                      max_r2 = max(estimate),
                      percentile = c(5,50,95),
                      CI = quantile(estimate,probs=c(0.05,0.5,.95)),
                      min_CI_lower = min(CI_lower),
                      max_CI_upper = max(CI_upper)) %>% 
    pivot_wider(names_from = percentile,names_prefix = "perc_",values_from = CI) %>% 
    mutate(R2=real$R2$estimate,
           CI_lower = real$R2$CI_lower,
           CI_upper = real$R2$CI_upper) 
  
  coef_sum<- coef %>%  group_by(term) %>% 
    dplyr::summarise(mean_beta = mean(estimate),
                     max_beta = max(estimate),
                     percentile = c(5,50,95),
                     CI = quantile(estimate,c(0.05,0.5,.95) )) %>% 
    pivot_wider(names_from = percentile,names_prefix = "perc_",values_from = CI) %>%
    inner_join(real$BW)
  
  return(list(r2=r2_sum, coef=coef_sum, allr2=r2,allcoef=coef))
}

make_permut_plots<- function( permout, title=""){
  
  specificityColors <- c( "#A3753B", "#CC9B57", "#E7CF97", "#F8EDD0", "#F7F7F7", "#D2EEEA", "#99D7CE", "#5DACA5", "#33847E")
  names(specificityColors) <-  c(1:9)
  
  coef_plot<- permout$coef %>% mutate(term=as.factor(gsub("e","",term)))%>% 
    ggplot(aes(x=term, y=estimate, 
               ymin=perc_5, ymax=perc_95,
               col=PD)) +
    geom_pointrange(position =  position_dodge(width=0.5))+
    scale_color_manual(values = specificityColors)+
    xlab("CRE pleiotropic degree")+
    ylab("coefficient")+
    theme_bw()+ theme(legend.position = "None")+ggtitle(title)
  
  r2_plot<-permout$allr2 %>% 
    ggplot(aes(x=PD,y=estimate, fill=PD)) +
    geom_violin()+
    geom_pointrange(data=permout$r2, 
                    aes(x=PD,y=estimate,
                        ymin=CI_lower,
                        ymax=CI_upper), color="darkgrey")+
    xlab("expression pleiotropic degree")+
    ylab(expression("marginal"~R^{2}))+
    scale_fill_manual(values=specificityColors)+
    theme_bw()
  
  cowplot::plot_grid(coef_plot,r2_plot, nrow = 1, axis="tblr",align = "h")
  
}

