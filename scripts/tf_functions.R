
library(org.Hs.eg.db)
library(topGO)
library(ReactomePA)


#do topGO
do_topGO<-function(all_genes_table, ont = "BP", return_table=T){
  rownames(all_genes_table)<-all_genes_table[,1]
  sign<-all_genes_table$padj<0.05
  print(table(sign))
  List<-rep(1,dim(all_genes_table)[1])
  names(List)<-rownames(all_genes_table)
  List[sign]<-0.01
  topDiffGenes<-function(x){return(x<0.05)}
  newClass<-new("topGOdata",description="sign",
                ontology=ont, allGenes=List,
                geneSel=topDiffGenes, annot=annFUN.org,
                nodeSize=10,ID="ENSEMBL",
                mapping="org.Hs.eg.db")
  if(!return_table){
    return(newClass)
  }
  numSigGenes(newClass)
  GO.res<-runTest(newClass,algorithm = "elim",statistic = "fisher")
  sum(score(GO.res)<=0.01)
  topGO_table<-GenTable(newClass,Fisher=GO.res,orderBy="Fisher",ranksOf="Fisher",topNodes=40,numChar=100)
  return(topGO_table)
}


plot_topGO<-function(topGO, p_cutoff = NULL, top_n = NULL, n_size_breaks=2){
  
  topGO$Term <- factor(topGO$Term, levels = unique(topGO$Term[order(topGO$Significant/topGO$Annotated)]))
  topGO$Significant<-as.integer(topGO$Significant)
  topGO$Annotated<-as.integer(topGO$Annotated)
  topGO$Fisher<-as.numeric(topGO$Fisher)
  
  if(!is.null(p_cutoff)){
    topGO<-topGO %>% dplyr::filter(Fisher<=p_cutoff)
  }
  
  if(!is.null(top_n)){
    topGO<-topGO[1:top_n,]
  }
  
  breaks = seq(min(topGO$Significant),max(topGO$Significant), length.out = n_size_breaks)
  breaks_x = round(seq(min(topGO$Significant/topGO$Annotated),max(topGO$Significant/topGO$Annotated), length.out = 3), 2)
  breaks_col = round(seq(min(topGO$Fisher),max(topGO$Fisher), length.out = 2), 2)
  
  ggplot(topGO %>% dplyr::rename(`Fisher's p`=Fisher), 
         aes(reorder(Term,Significant/Annotated),
             Significant/Annotated, 
             fill=`Fisher's p`))+ #, size=Significant
    geom_point(size=2.3,colour="black",pch=21)+
    theme_classic()+
    coord_flip()+
    scale_fill_continuous(low = "grey80", high = "grey20",
                           breaks=breaks_col)+#,breaks=c(0.015, 0.025, 0.035))+
    theme(axis.title = element_text(size=8),
          axis.text = element_text(size=7, color="black"),
          axis.title.y=element_blank(),
          axis.text.y=element_text())+
    guides(fill=guide_colorbar(barwidth = 0.5, barheight=1.5))+
    theme(legend.title = element_text(size = 6.5, margin = margin(b=0)), 
          legend.text = element_text(size = 6, margin=margin(t=0,b=0,l=0,r=0)),
          legend.background = element_rect(fill="transparent",colour=NA),
          legend.spacing.y = unit(1,"mm"),
          legend.spacing.x = unit(1,"mm"),
          legend.position = c(0.83,0.25)) +
    scale_y_continuous(breaks=breaks_x)
  
}





get_topGO_genes<-function(full_DE_table=NULL,topGO_summary=NULL, topGO_table, term_vector, annot, genesOfInterest=NULL){
  
  topGO_summary<-do_topGO(full_DE_table, return_table = F)
  
  all_terms<-list()
  for (i in 1:length(term_vector)){
    
    top_terms<-genesInTerm(topGO_summary, term_vector[[i]])
    GOsig<-intersect(top_terms[[1]], genesOfInterest)
    
    all_terms[[topGO_table$Term[topGO_table$GO.ID==term_vector[i]]]]<-data.frame(gene_id=as.factor(GOsig)) %>%
      inner_join(full_DE_table) %>%
      mutate(GO.ID=as.factor(term_vector[i])) %>%
      left_join(topGO_table) %>%
      left_join(annot)
  }
  return(all_terms)
}


make_cnet_from_topGO<-function(full_DE_table = NULL, topGO_object=NULL, 
                               topGO_table, 
                               term_vector,
                               annot,
                               genesOfInterest=NULL,
                               pointval){
  
  ggs<-get_topGO_genes(full_DE_table=full_DE_table,
                       topGO_table=topGO_table,
                       term_vector=term_vector, 
                       annot = annot,
                       genesOfInterest=genesOfInterest)
  
  
  ggs_red<-lapply(ggs,function(x){x %>% pull(TF)})
  ggs_red2<-bind_rows(ggs) %>% distinct(TF, .keep_all = T)
  lfcs<-ggs_red2 %>% pull({{pointval}})
  names(lfcs)<-ggs_red2$TF

  cnetplot(ggs_red, foldChange = lfcs, showCategory = length(term_vector), 
           cex_label_gene=0.65*0.5, cex_label_category=0.85*0.5, cex_category=0.7, cex_gene=0.7)
  
}




ranking_by_abundance<-function(tis = NULL, reg_id_selection, dd){
  
  print(tis)
  # select only expressed TFs
  expressedTFs_tissues<-readRDS("../ATACseq/cbust/expressedTFs_inTissues.rds")
  
  if (!is.null(tis)){
    expressedInTissue<-readRDS("../roadmap_expression_summaries/summarized_expression_allTissues.rds") %>%
      filter(tissue == tis) %>% 
      inner_join(expressedTFs_tissues %>% dplyr::select(gene_id, motif_ID))
  } else {
    expressedInTissue<-readRDS("../roadmap_expression_summaries/summarized_expression_allTissues.rds") %>%
      distinct(gene_id, .keep_all = T) %>% 
      inner_join(expressedTFs_tissues %>% dplyr::select(gene_id, motif_ID))
  }
  
  
  # generate the ranks based on how many CREs are bound by a TF in each PD (normalized by CREs in that PD)
  ranked<-dd %>%
    as_tibble() %>%
    dplyr::filter(region_id %in% reg_id_selection) %>%
    group_by(cluster_id_or_motif_name, total) %>%
    dplyr::summarise(n_CREs = length(region_id)) %>%
    dplyr::mutate(motif_ID = stringr::word(cluster_id_or_motif_name,1,1,"_")) %>%
    dplyr::filter(motif_ID %in% expressedInTissue$motif_ID) %>%
    inner_join(dd %>%
                 dplyr::filter(region_id %in% reg_id_selection) %>%
                 group_by(total) %>%
                 dplyr::summarise(n_total_CREs = length(total))) %>%
    mutate(prop_CREs=n_CREs/n_total_CREs) %>%
    arrange(-prop_CREs) %>%
    group_by(cluster_id_or_motif_name) %>%
    dplyr::mutate(rank=row_number())
  return(ranked)
}






ranks_to_topGO<-function(tis = NULL, ranked=NULL, PD, FC = 1.5, n_cnet=4, PD_filt){
  
  print(tis)
  # select only expressed TFs
  expressedTFs_tissues<-readRDS("../ATACseq/cbust/expressedTFs_inTissues.rds")
  
  if (!is.null(tis)){
    expressedInTissue<-readRDS("../roadmap_expression_summaries/summarized_expression_allTissues.rds") %>%
      filter(tissue == tis) %>% 
      inner_join(expressedTFs_tissues %>% dplyr::select(gene_id, motif_ID))
  } else {
    expressedInTissue<-readRDS("../roadmap_expression_summaries/summarized_expression_allTissues.rds") %>%
      distinct(gene_id, .keep_all = T) %>% 
      inner_join(expressedTFs_tissues %>% dplyr::select(gene_id, motif_ID))
  }
  
  # calculate FCs and add TF names for cnet annotations
  annot<-expressedTFs_tissues %>% 
    ungroup() %>%
    dplyr::distinct(gene_id, TF, motif_ID, .keep_all = T) %>% 
    inner_join(ranked %>% 
                 ungroup() %>%
                 filter(total == PD) %>% 
                 dplyr::select(motif_ID, prop_CREs, cluster_id_or_motif_name)) %>%
    # add the average abundance across PDs
    inner_join(ranked %>% 
                 ungroup() %>%
                 filter(total != PD) %>%
                 group_by(motif_ID) %>%
                 dplyr::summarise(avr_prop_CREs=mean(prop_CREs)), by="motif_ID") %>%
    dplyr::transmute(gene_id,prop_CREs, avr_prop_CREs, TF, cluster_id_or_motif_name, FC=prop_CREs/avr_prop_CREs)
  
  # prepare for topGO (assign p-values as the threshold of foreground, background)
  motifs_ranks<-ranked %>%
    dplyr::distinct(cluster_id_or_motif_name) %>%
    dplyr::mutate(isPD = ifelse(cluster_id_or_motif_name %in% PD_filt, T, F),
                  padj = ifelse(isPD, 0.01, 1)) %>%
    as_tibble() %>%
    dplyr::mutate(motif_ID = stringr::word(cluster_id_or_motif_name,1,1,"_")) %>%
    left_join(expressedInTissue[,c("motif_ID","gene_id")] %>% distinct(.keep_all=T)) %>%
    relocate(gene_id) %>%
    group_by(gene_id) %>%
    slice_min(padj, n=1, with_ties = F) %>%
    ungroup() %>%
    as.data.frame() 

  print(table(motifs_ranks$isPD))
  
  # Do topGO (padj is the decisive column for foreground / background)
  topGO_PD_TFs_ranked<-do_topGO(motifs_ranks, ont="BP")
  
  # Make cnets
  cnet_PD<-make_cnet_from_topGO(full_DE_table = motifs_ranks,
                                annot = annot,
                                topGO_table = topGO_PD_TFs_ranked,
                                term_vector = topGO_PD_TFs_ranked$GO.ID[1:n_cnet],
                                genesOfInterest = motifs_ranks$gene_id[motifs_ranks$padj<0.05],
                                pointval = "FC")+ 
    scale_color_gradient(name="Odds\nratio",low='#F5E9E2', high='#a53860',
                         limits=c(0.95,1.85), breaks=c(1,1.4,1.8))+
    guides(color=guide_colorbar(barwidth = 0.5, barheight=1.5))+
    theme(plot.title = element_text(size=8),
          legend.title = element_text(size = 7.5), 
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.5, "lines"),
          legend.margin=margin(0.1,0.1,0.1,0.1),
          legend.box.margin=margin(0.1,0.1,0.1,0.1))+
    guides(size="none")
  
  
  return(list(topGO_PD = topGO_PD_TFs_ranked,
              cnet_PD = cnet_PD))
  
}