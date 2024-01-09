
#helper functions for the cbust and phastcons analyses

writeFasta4cbust<-function(gr, genome, fasta.file, id.col= "peak_id", other.id,
                           padding=500 ,max.seq=200){
  require(BSgenome)
  require(Biostrings)
  
  gr <- gr %>% anchor_center() %>% 
    stretch(2*padding)
  x  <- 1:length(gr)
  dd <- split( x, ceiling(x/max.seq))
  
  for(j in 1:length(dd) ){
    i<- dd[[j]]
    s<-getSeq( genome, gr[i])
    names(s) <- paste0(seqnames(gr[i]),":",start(gr[i]),"-",end(gr[i]),"@@",mcols(gr)[[id.col]][i])
    if(!missing(other.id)){
      names(s)<- paste( names(s), mcols(gr)[[other.id]][i] )
    } 
    writeXStringSet(s, paste0(fasta.file,"_",j,".fa"))
  }
}



extract_region_fa<-function(genome_path, coord_file, genome, extension=300, out_folder, ncores=10){
  
  gen<-readDNAStringSet(genome_path)
  chrom_sizes<-data.frame(seqnames=stringr::word(names(gen),1,1," "), chrom_width=width(gen),
                          seqnames_full=names(gen))
  
  coord_500bExt<- coord_file %>%
    as_tibble() %>% 
    #filter(similar=="yes") %>%
    left_join(chrom_sizes) %>% 
    rowwise() %>%
    dplyr::mutate(start_up=ifelse(start-extension>=0,start-extension,0),
                  end_down=ifelse(end+extension<=chrom_width,end+extension,chrom_width)) %>%
    dplyr::mutate(extended_by=(end_down-start_up) - (end-start)) %>%
    dplyr::select(seqnames_full,start_up,end_down,peak_id, extended_by) %>%
    dplyr::rename(seqnames=seqnames_full,start=start_up, end=end_down) 
  #saveRDS(coord_500bExt,paste0(out_folder,"/",genome,"_coord_",extension,"ext.rds"))
  
  #and cut out the sequences 
  registerDoParallel(ncores)  # use multicore, set to the number of cores --> separate by chromosome
  foreach (seq=unique(coord_500bExt$seqnames)) %dopar% {
    to_cut_coord<-coord_500bExt %>% filter(seqnames==seq)
    to_cut_fa<-gen[names(gen)==seq]
    
    #very annoying, subseq did not work on all seqs at once so i'm looping
    peak_string<-DNAStringSet()
    for (i in 1:nrow(to_cut_coord)){
      print(to_cut_coord$peak_id[i])
      subs<-subseq(to_cut_fa,start=to_cut_coord$start[i], end=to_cut_coord$end[i])
      #names(subs)<-to_cut_coord$peak_id[i]
      names(subs)<-paste0("chr",to_cut_coord$seqnames[i],":",to_cut_coord$start[i],"-",to_cut_coord$end[i],"@@", to_cut_coord$peak_id[i])
      peak_string<-c(peak_string,subs)
    }
    writeXStringSet(peak_string,paste0(out_folder,"/",genome,"/peaks_",seq,".fa"))
  }
}







analysePhastCons <- function( bigWigFile, gr, probcut=0.9){
  phastcons  <- rtracklayer::import(bigWigFile, 
                                    which= gr,
                                    as="NumericList")
  sumP<- sapply(phastcons, function(x){ 
    n<-length(x)
    c(n, sum(x), sum(x>probcut) )
  }) %>% t() %>% data.frame()
  
  names(sumP)<- c("bp","sumP","neg_n")
  gr %>% as_tibble() %>% bind_cols(sumP) %>%
    dplyr::rowwise() %>% mutate(meanCons = sum(sumP)/sum(bp),
                                fracCons = sum(neg_n)/sum(bp),
                                bp = sum(bp))
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

check_baseFreq<-function(ggranges_obj){
  
  out<-alphabetFrequency(ggranges_obj) %>% 
    as.data.frame() %>% #there are currently no Ns, nice
    rowwise() %>%
    mutate(GC=sum(G,C), 
           length=sum(G,C,A,T), 
           prop_GC=GC/length,
           A=A/length,
           T=T/length,
           C=C/length,
           G=G/length) 
  return(out)
}




mergePeaks<- function( genome, ff , spec, cellt ,max.peak.size=5000,canonical=F ){
  gr<-lapply( 1:length(ff), function(i){
    read_narrowpeaks(ff[i]) %>% 
      plyranges::mutate(species = spec[i], 
                        celltype = cellt[i] ) } ) %>% 
    bind_ranges() %>% 
    filter(width <= max.peak.size) %>% 
    reduce_ranges(celltype = paste(unique(celltype), collapse=","),
                  species  = paste(unique(species), collapse = ",")) %>% 
    plyranges::mutate( peak_id = paste(genome,1:length(.),sep = "_") )
  
  if(canonical){
    return(keepStandardChromosomes(gr,pruning.mode = "coarse"))
  }else{
    return(gr)
  }
  
}



get_Ns<-function(fa_folder, base_folder="ATACseq/cbust/fastas/", cut_id = TRUE){
  
  cbust<-c(list.files(paste0(base_folder,fa_folder),pattern = ".fa",recursive = T))
  seqs<-lapply(cbust,function(x){readDNAStringSet(paste0(base_folder,fa_folder,"/",x))})
  names(seqs)<-cbust
  
  Ns<-lapply(seqs, function(x){data.frame(region_id=as.factor(names(x)), width=width(x)) %>% 
      bind_cols(alphabetFrequency(x) %>% as.data.frame() %>% dplyr::select(N))}) %>%
    bind_rows() 
  
  if (cut_id){
    Ns<-Ns %>% dplyr::mutate(region_id = as.factor(word(region_id,2,2,"@@")))
  }
  Ns
}





calc_OL<-function(tab1, tab2, filt = T, rel_width_OL_cutoff = 0.1, id1_col = "region_id", id2_col="peak_id"){
  
  tab1_gr<-tab1 %>%
    as_tibble() %>%
    mutate(index1 = 1:nrow(.), width1 = width) %>%
    as_granges()  
  
  tab2_gr<-tab2 %>%
    as_tibble() %>%
    mutate(index2 = 1:nrow(.), width2 = width) %>%
    as_granges() 
  
  # quantify the overlap
  hits <- findOverlaps(tab1_gr, tab2_gr)
  overlaps <- pintersect(tab1_gr[queryHits(hits)], tab2_gr[subjectHits(hits)])
  hits_overlaps<-data.frame(hits, width(overlaps))
  
  overlaps_df<-join_overlap_left(tab1_gr, tab2_gr) %>%
    as_tibble() %>%
    inner_join(hits_overlaps, by=c("index1"="queryHits", "index2"="subjectHits")) %>%
    dplyr::rename(overlap=`width.overlaps.`) %>%
    mutate(fracOverlap1=overlap/width1,
           fracOverlap2=overlap/width2)
  
  if (filt){
    overlaps_df<-overlaps_df %>% filter(fracOverlap1>rel_width_OL_cutoff, fracOverlap2>rel_width_OL_cutoff)
  }
  
  overlaps_df<-overlaps_df %>%
    group_by(index1) %>%
    mutate(n_Ind1ToInd2 = length(unique(index2))) %>%
    group_by(index2) %>%
    mutate(n_Ind2ToInd1= length(unique(index1))) %>%
    ungroup() %>%
    mutate(OneToOne = ifelse(n_Ind1ToInd2==1 & n_Ind2ToInd1==1, T, F)) %>%
    dplyr::select( {{id1_col}}, {{id2_col}}, fracOverlap1, fracOverlap2, width1, width2, n_Ind1ToInd2, n_Ind2ToInd1, OneToOne)
  
  return(overlaps_df)
}




