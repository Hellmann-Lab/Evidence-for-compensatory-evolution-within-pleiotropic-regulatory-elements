
# REANALYSE ROLLER 2021 DATA ----------------------------------------------

# Steps:
# 1) Translate PD region_id positions from hum and cyno to mmul using liftover (same version as used in the paper: mmul10)
# 2) Combine that with the annotations from the other species
# 3) Recheck PD agreement
# 4) Anything different for testis? 
# 5) Check the evo age for different PD 

#setwd("/data/share/htp/pleiotropy/paper_data/")
library(tidyverse)
library(data.table)
library(plyranges)
library(ggtree)
library(ape)
library(Biostrings)


# these are all the helper functions for the project
source("scripts/helper_functions.R")
source("ATACseq/scripts/functions.R")
source("ATACseq/scripts/helper_functions.R")



basic_theme_ins2<-  theme(axis.title=element_text(size=8),
                          axis.text = element_text(color="black", size=7),
                          legend.title = element_text(size = 7.3), 
                          legend.text = element_text(size = 7),
                          legend.key.size = unit(0.5, "lines"),
                          legend.margin=margin(0,0,0,0),
                          legend.box.margin=margin(0,0,0,0),
                          legend.background = element_blank(), 
                          legend.box.background = element_blank(),
                          panel.grid.major = element_line(size = 0.35),
                          panel.grid.minor = element_blank(),
                          strip.text = element_text(size = 8))





# JAMM Liftover to MMUL ---------------------------------------------------


# 1) --> need our PDs; liftover script, mmul coords and genome version

jamm_region_identifiers<-readRDS("roadmap_DHS_summaries/region_summary/jamm_region_info_full.rds") %>%
  transmute(seqnames=chromosome, start,end, region_id, total, brain, muscle)

additional<-readRDS("CRE_to_Gene/DHS_to_gene.rds") %>% 
  distinct(region_id, assignment, CGI)

jamm_hg38_noExt<-translate_jamm(chain_file = "liftOvers/hg19ToHg38.over.chain",
                                coordinate_file = as_granges(jamm_region_identifiers), 
                                extend = 50,
                                reverse_chain_file = "liftOvers/hg38ToHg19.over.chain")

# now translate both to mmul10 --> nope, there are still no macfas6 LO files, just do human

jamm_hg38_to_mmul10<-translate_jamm(chain_file = "liftOvers/hg38ToRheMac10.over.chain",
                                    coordinate_file = as_granges(jamm_hg38_noExt), 
                                    extend = 50,
                                    reverse_chain_file = "liftOvers/rheMac10ToHg38.over.chain") %>% 
  left_join(jamm_region_identifiers %>% dplyr::select(region_id, brain, muscle,total)) %>% 
  left_join(additional) %>% 
  as_granges()

length(jamm_hg38_to_mmul10)
dim(jamm_hg38_noExt)
length(jamm_hg38_to_mmul10)/dim(jamm_hg38_noExt)[1]





# Get all macaque regions, compare to our DHS -----------------------------
# I want to compare: widths, numbers, number of overlaps, number of multi-overlaps, tissue specificity, annotation type

macfiles<-list.files("Roller2021/Chipseq_tables", pattern = "Macaque_regRegions", full.names = T)

# these coordinates can overlap within the file, if the enh/prom role of the CRE is different between tissues
macregs<-lapply(macfiles, function(x){
  data.table::fread(x) %>% 
    dplyr::select(V4, V5, V6)
}) %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  # somehow there are around 2k coords that are doubled.. correct this and write them an email
  mutate(n_separators = str_count(V4, pattern=":|-"),
         coords = ifelse(n_separators>2, str_sub(V4, start=1, end=str_length(V4)/2), V4),
         mmul_coords = coords) %>% 
  distinct(mmul_coords, .keep_all = T) %>% 
  separate(col = "coords", into = c("seqnames","start","end"), sep = ":|-") %>% 
  transmute(seqnames = paste0("chr",seqnames),
            start = as.numeric(start), 
            end = as.numeric(end), 
            mmul_region_id = 1:length(mmul_coords),
            annot_mmul = V5, 
            tissues_mmul = V6) %>% 
  as_granges()

seqlevelsStyle(macregs) <- "UCSC"
macregs<-keepStandardChromosomes(macregs, pruning.mode="coarse")
saveRDS(macregs,"Roller2021/chipseq_macregs.rds")






###############################
# COLLAPSE OL MACAQUE CREs ####
###############################


macregs_collapsed<-join_overlap_left(reduce_ranges(macregs) %>% 
                                       mutate(collapsed_region_id = 1:length(.)), 
                                     macregs) %>% 
  as_tibble() %>% 
  group_by(seqnames, start, end, collapsed_region_id) %>% 
  dplyr::summarise(annot_mmul = paste(sort(annot_mmul), collapse = ","),
                   tissues_mmul = paste(sort(trimws(strsplit(tissues_mmul[1], ',')[[1]])), collapse=',')) %>% 
  dplyr::mutate(total_mmul = str_count(tissues_mmul, pattern = ",")+1)

saveRDS(macregs_collapsed, "Roller2021/chipseq_macregs_collapsed.rds")

OL_jamm_mmul_col<-calc_OL(as_granges(jamm_hg38_to_mmul10), 
                          as_granges(macregs_collapsed), 
                          filt = T, 
                          rel_width_OL_cutoff = 0.1, 
                          id1_col = "region_id", 
                          id2_col = "collapsed_region_id") %>% 
  left_join(jamm_hg38_to_mmul10 %>% 
              as_tibble() %>% 
              dplyr::select(-seqnames,-start,-end)) %>% 
  left_join(macregs_collapsed) %>% 
  dplyr::mutate(brain_mmul = ifelse(grepl("Brain", tissues_mmul), 1, 0),
                muscle_mmul = ifelse(grepl("Muscle", tissues_mmul), 1, 0))

saveRDS(OL_jamm_mmul_col, "Roller2021/chipseq_OL_humJamm_mac_collapsed.rds")


# width comparison
ggplot(OL_jamm_mmul_col, aes(x = width1, y = width2))+
  theme_bw()+
  geom_abline(intercept=0, slope=1)+
  geom_point()+
  xlab("DHS width (bp)")+
  ylab("ChIP-CRE width (bp)")


length(jamm_hg38_to_mmul10)
length(macregs_collapsed)
length(unique(OL_jamm_mmul_col$region_id))/length(jamm_hg38_to_mmul10) # only 80k are present in macaque
length(unique(OL_jamm_mmul_col$collapsed_region_id)) # on
dim(OL_jamm_mmul_col)

# 1-to-1, many-to-1, 1-to-many summaries
nOL2<-OL_jamm_mmul_col %>% 
  group_by(n_Ind1ToInd2,n_Ind2ToInd1) %>% 
  dplyr::summarise(n = length(n_Ind2ToInd1)) 

# 65.6% of OLs are 1-to-1; 33.4% and 0.8% are 1-to-many or many-to-1, with much more multi-DHS vs 1 cipseq
nOL2$n[nOL2$n_Ind1ToInd2==1 & nOL2$n_Ind2ToInd1==1]/sum(nOL2$n)
sum(nOL2$n[nOL2$n_Ind1ToInd2==1 & nOL2$n_Ind2ToInd1>1])/sum(nOL2$n)
sum(nOL2$n[nOL2$n_Ind1ToInd2>1 & nOL2$n_Ind2ToInd1==1])/sum(nOL2$n)
sum(nOL2$n[nOL2$n_Ind1ToInd2>1 & nOL2$n_Ind2ToInd1>1])/sum(nOL2$n) # basically no many to many



# Plot PD agreement between hum and mac (1-to-1 CREs) -----------------------------------

OL_jamm_mmul_1to1<-OL_jamm_mmul_col %>% 
  dplyr::filter(OneToOne)

PDmac<-OL_jamm_mmul_1to1 %>% 
  distinct(region_id,total,total_mmul) %>%
  group_by(total) %>% 
  mutate(tot =length(total)) %>% 
  group_by(total,total_mmul,tot) %>% 
  summarize(cnt=length(total_mmul)) %>%
  ggplot(aes(x=as.factor(total), y=cnt/tot, fill=as.factor(total_mmul)))+
  theme_bw()+
  geom_bar(stat="identity")+ 
  scale_fill_manual(values = specificityColors)+
  ylab("Fraction of CREs")+
  guides(fill=guide_legend(title="PD in\nmacaque"))+
  xlab("PD")+
  basic_theme_ins2+
  theme(legend.title = element_text(size=7),
        axis.text = element_text(color="black"))



# How many brain and muscle CREs overlap? -----------------------------------------

table(OL_jamm_mmul_1to1$brain, OL_jamm_mmul_1to1$brain_mmul)
fisher.test(OL_jamm_mmul_1to1$brain, OL_jamm_mmul_1to1$brain_mmul)

table(OL_jamm_mmul_1to1$muscle, OL_jamm_mmul_1to1$muscle_mmul)
fisher.test(OL_jamm_mmul_1to1$muscle, OL_jamm_mmul_1to1$muscle_mmul)






#############################
# INCLUDE OTHER SPECIES #####
#############################


macToOther_merged_collapsed<-lapply(macfiles, function(x){
  
  macToOther_red<-data.table::fread(x) %>% 
    as_tibble() %>% 
    dplyr::mutate(n_separators = str_count(V4, pattern=":|-"),
                  coords = ifelse(n_separators>2, 
                                  str_sub(V4, start=1, end=str_length(V4)/2), 
                                  V4),
                  mmul_coords = coords) %>% 
  distinct(mmul_coords, .keep_all = T) %>% 
  separate(col = "coords", into = c("seqnames","start","end"), sep = ":|-") %>% 
  dplyr::transmute(seqnames = paste0("chr",seqnames),
            start = as.numeric(start), 
            end = as.numeric(end), 
            annot_other = V11, 
            tissues_other = V12, 
            n_other = str_count(tissues_other, pattern = ",")+1,
            species_other = str_split_i(V13, pattern = "_", 1)) %>% 
  as_granges()
  
  join_overlap_inner(as_granges(macregs_collapsed), macToOther_red) %>% 
    as_tibble() %>% 
    group_by(seqnames, start, end, width, collapsed_region_id, 
             annot_mmul, tissues_mmul, total_mmul) %>% 
    dplyr::summarise(annot_other = paste(sort(annot_other), 
                                         collapse = ","),
                     tissues_other = paste(sort(trimws(strsplit(tissues_other[1], ',')[[1]])),
                                           collapse=',')) %>% 
    dplyr::mutate(total_other = str_count(tissues_other, pattern = ",")+1,
                  species_other = macToOther_red$species_other[1]) %>% 
    dplyr::inner_join(OL_jamm_mmul_col %>% 
                        distinct(region_id, collapsed_region_id, n_Ind1ToInd2, n_Ind2ToInd1, OneToOne))
  
})


saveRDS(macToOther_merged_collapsed, "Roller2021/chipseq_PD.rds")







########################
# HUMAN VS MARMOSET ####
########################

# the used marmoset genome assembly does not seem to have a liftover chain file available. use the PD annotations from Roller via hum <-> mac <-> marmoset
marm<-lapply(macToOther_merged_collapsed, function(x){
  if (x$species_other[1]=="Marmoset"){
    return(x)
  }
}) %>% bind_rows()

marm<-marm %>% left_join(jamm_region_info_full %>% dplyr::select(region_id, total))


PDmarm<-marm %>% 
  filter(OneToOne) %>% 
    distinct(region_id,total,total_other) %>%
    group_by(total) %>% 
    mutate(tot =length(total)) %>% 
    group_by(total,total_other,tot) %>% 
    summarize(cnt=length(total_other)) %>%
    ggplot(aes(x=as.factor(total), y=cnt/tot, fill=as.factor(total_other)))+
    theme_bw()+
    geom_bar(stat="identity")+ 
    scale_fill_manual(values = specificityColors)+
    ylab("Fraction of CREs")+
    guides(fill=guide_legend(title="PD in\nmarmoset"))+
    xlab("PD")+
    basic_theme_ins2+
    theme(legend.title = element_text(size=7),
          axis.text = element_text(color="black"))

plot_grid(PDmac,PDmarm)  






################
# SANKEY PD ####
################


jamm_region_identifiers<-readRDS("roadmap_DHS_summaries/region_summary/jamm_region_info_full.rds") %>%
  transmute(seqnames=chromosome, start,end, region_id, total, brain, muscle)

macToOther_merged_collapsed<-readRDS("Roller2021/chipseq_PD.rds")
macregs_collapsed<-readRDS("Roller2021/chipseq_macregs_collapsed.rds")

allSpecies<-bind_rows(macToOther_merged_collapsed) %>% 
  bind_rows(macregs_collapsed %>% dplyr::inner_join(OL_jamm_mmul_col) %>% 
              dplyr::transmute(region_id, total_other=total_mmul, species_other = 'Macaque', OneToOne)) %>% 
  ungroup() %>% 
  left_join(jamm_region_identifiers %>% dplyr::select(region_id, total)) %>% 
  dplyr::filter(OneToOne) %>% 
  dplyr::select(region_id, total, total_other, species_other)


# lower this one to only 5 species
allSpecies_wide<-allSpecies %>% 
  dplyr::filter(species_other %in% c("Macaque", "Marmoset")) %>%  #, "Mouse", "Rat")) %>% 
  pivot_wider(id_cols = c(region_id,total), 
              names_from = "species_other", 
              values_from = "total_other", 
              values_fn = function(x){min(4,sum(x))}) %>% 
  ungroup() %>% 
  dplyr::select(-region_id) 

allSpecies_wide<-allSpecies_wide[order(allSpecies_wide$total, allSpecies_wide$Macaque, allSpecies_wide$Marmoset),]


links <- allSpecies_wide %>% 
  dplyr::mutate(row = row_number()) %>%  # add a row id
  pivot_longer(-row, names_to = "column", values_to = "source") %>%  # gather all columns
  dplyr::mutate(column = match(column, names(allSpecies_wide))) %>%  # convert col names to col ids
  group_by(row) %>%
  dplyr::mutate(target = lead(source, order_by = column)) %>%  # get target from following node in row
  ungroup() %>% 
  filter(!is.na(target) & !is.na(source)) %>% 
  dplyr::mutate(source = paste0(source, '_', column)) %>%
  dplyr::mutate(target = paste0(target, '_', column + 1)) %>%
  dplyr::select(source, target)

nodes <- data.frame(name = unique(c(links$source, links$target)))
nodes$label <- sub('_[0-9]*$', '', nodes$name) # remove column id from node label

links$source_id <- match(links$source, nodes$name) - 1
links$target_id <- match(links$target, nodes$name) - 1
links$value <- 1
links$group <- as.factor('a')

library(networkD3)
my_color <- 'd3.scaleOrdinal() .domain(["a", "1", "2", "3", "4", "5", "6", "7", "8", "9"]) .range(["#dad7cd", "#a3753b", "#cc9B57", "#e7cf97", "#f8edd0", "#f7f7f7", "#d2eeea", "#99d7ce", "#5daca5","#33847e"])'

sankey_out = sankeyNetwork(Links = links, Nodes = nodes,
                           Source = 'source_id', LinkGroup = 'group',
                           colourScale = my_color, 
                           Target = 'target_id', Value = 'value', NodeGroup='label',
                           NodeID = 'label', iterations = 0, fontSize =20)
saveNetwork(sankey_out, file = "Roller2021/figures/tmp.html")

library(webshot)
webshot( "Roller2021/figures/tmp.html", 
         "Roller2021/figures/humMacMarmoset.png", 
         vwidth = 800*1.1, vheight = 500*1.1)






#########################
# avr PD and CRE age ####
#########################


# Construct a tree --------------------------------------------------------

mammaltree<-read.tree("/data/share/htp/TRNP1/paper_data/Co-evolution-TRNP1-and-GI/protein/trees/mammaltree.txt")
speciesTree<-drop.tip(mammaltree, mammaltree$tip.label[!mammaltree$tip.label %in% c("Macaca_mulatta","Callithrix_jacchus","Felis_catus","Canis_lupus", "Equus_ferus", "Mus_musculus", "Rattus_norvegicus", "Oryctolagus_cuniculus", "Sus_scrofa", "Monodelphis_domestica")])


# Get annotations of region_id vs species ---------------------------------

jamm_region_identifiers<-readRDS("roadmap_DHS_summaries/region_summary/jamm_region_info_full.rds") %>%
  transmute(seqnames=chromosome, start,end, region_id, total, brain, muscle)

macToOther_merged_collapsed<-readRDS("Roller2021/chipseq_PD.rds")
macregs_collapsed<-readRDS("Roller2021/chipseq_macregs_collapsed.rds")
OL_jamm_mmul_col<-readRDS("Roller2021/chipseq_OL_humJamm_mac_collapsed.rds")

allSpecies<-bind_rows(macToOther_merged_collapsed) %>% 
  bind_rows(macregs_collapsed %>% dplyr::inner_join(OL_jamm_mmul_col) %>% 
              dplyr::transmute(region_id, total_other=total_mmul, species_other = 'Macaque', OneToOne)) %>% 
  ungroup() %>% 
  left_join(jamm_region_identifiers %>% dplyr::select(region_id, total)) %>% 
  dplyr::filter(OneToOne) %>% 
  dplyr::select(region_id, total, total_other, species_other) %>% 
  dplyr::mutate(latin=case_when(species_other == "Macaque" ~ "Macaca_mulatta",
                                species_other == "Marmoset" ~ "Callithrix_jacchus",
                                species_other == "Cat" ~ "Felis_catus",
                                species_other == "Dog" ~ "Canis_lupus", 
                                species_other == "Horse" ~ "Equus_ferus", 
                                species_other == "Mouse" ~ "Mus_musculus", 
                                species_other == "Rat" ~ "Rattus_norvegicus",
                                species_other == "Rabbit" ~ "Oryctolagus_cuniculus", 
                                species_other == "Pig" ~ "Sus_scrofa", 
                                species_other == "Opossum" ~ "Monodelphis_domestica"))





# Average PD across the tree ----------------------------------------------
# also quantify the average number of tissues in which it is used across species

meanPD<-allSpecies %>% 
  group_by(region_id, total) %>% 
  summarize(mean_total_other = mean(total_other)) %>% 
  ggplot(aes(x=as.factor(total), y=mean_total_other, fill=as.factor(total)))+
  geom_boxplot(notch = T, outlier.size = 0.1)+
  theme_bw()+
  xlab("PD")+
  ylab(expression(mu [n(tissues)]))+
  scale_fill_manual(values=specificityColors)+
  basic_theme_ins2+
  theme(legend.position = "none",
        axis.text = element_text(color="black"))

ggsave("Roller2021/figures/avrTissues.png", height=3, width=4)






# Make a supplementary figure ---------------------------------------------

sup_roller<-plot_grid(PDmac +theme(legend.position = "top"),
                      PDmarm +theme(legend.position = "top"),
                      meanPD, ncol=3)  
figSankey<- ggdraw() + 
  cowplot::draw_image(magick::image_read("figures/humMacMarmoset_phylopic.png")) +
  theme(plot.margin = unit(c(0, -0.3, 0, 0), "cm"))



sup_roller2<-plot_grid(PDmac,
                       PDmarm,
                       plot_grid(NULL,figSankey, rel_widths=c(0.035,1)),
                       plot_grid(NULL,meanPD,NULL, ncol=3, rel_widths = c(0.05,1,0.24)),
                       scale=c(0.95,0.95,1,0.95),
                       labels=c("A","B","C","D"), label_size=12)#, align = "hv", axis="tblr")  

sup_roller2
ggsave("figures/supplRoller.pdf", width =162.5, height=105, units = "mm")







# CRE phylogenetic age ----------------------------------------------------

#calculate_tree_length(df=allSpecies, regid = "228816", tree = speciesTree)

treeLengths<-allSpecies %>%
  distinct(region_id) %>% 
  rowwise() %>%
  dplyr::mutate(tree_length=calculate_tree_length(df=allSpecies, regid=region_id, tree=speciesTree)) %>%
  ungroup() %>%
  dplyr::mutate(max_tree_length=max(tree_length)) %>%
  rowwise() %>%
  mutate(rel_tree_length=tree_length/max_tree_length) %>%
  ungroup() %>% 
  left_join(jamm_region_identifiers %>% dplyr::select(region_id, total))

saveRDS(treeLengths, "Roller2021/PhyloTreeLengths.rds")

evoage<-ggplot(treeLengths, aes(x=as.factor(total), y=rel_tree_length, fill=as.factor(total)))+
  geom_boxplot(notch = T)+
  theme_bw()+
  xlab("CRE PD")+
  ylab(expression(paste("Evolutionary age (relative ", lambda, ")")))+
  scale_fill_manual(values=specificityColors)+
  basic_theme_ins2+
  theme(legend.position = "none",
        axis.text = element_text(color="black"))

evoage
ggsave("Roller2021/figures/evoAge.png", height=3, width=4)
ggsave("Roller2021/figures/evoAge.pdf", height=3, width=4)



# also plot species tree
#speciesTree$tip.label<-gsub("_"," ",speciesTree$tip.label)

speciesTree_plot<-speciesTree
speciesTree_plot$tip.label<- case_when(speciesTree_plot$tip.label == "Macaca_mulatta" ~ "Macaque",
                                       speciesTree_plot$tip.label == "Callithrix_jacchus" ~ "Marmoset",
                                       speciesTree_plot$tip.label == "Felis_catus" ~ "Cat",
                                       speciesTree_plot$tip.label == "Canis_lupus" ~ "Dog", 
                                       speciesTree_plot$tip.label == "Equus_ferus" ~ "Horse", 
                                       speciesTree_plot$tip.label == "Mus_musculus" ~ "Mouse", 
                                       speciesTree_plot$tip.label == "Rattus_norvegicus" ~ "Rat",
                                       speciesTree_plot$tip.label == "Oryctolagus_cuniculus" ~ "Rabbit", 
                                       speciesTree_plot$tip.label == "Sus_scrofa" ~ "Pig", 
                                       speciesTree_plot$tip.label == "Monodelphis_domestica" ~ "Opossum")

mamtree<-ggtree(speciesTree_plot, lwd=0.3)+
  geom_tiplab(fontface="italic", size=2, geom = "text")+ 
  xlim(0,200) +
  #basic_theme_ins2+
  theme(plot.margin = unit(c(0, 0, 0,0), "cm")) +
  geom_treescale(width=20, x=0.05, y=10, offset=0.05)

mamtree
ggsave("Roller2021/figures/species_tree.pdf", height = 5*0.85, width=7*0.85)

plot_grid(mamtree,evoage, rel_widths = c(1,1.4), scale = c(0.9, 1), labels=c("C","D"), label_size=12, label_x = c(0,-0.05))
ggsave("figures/PD_evoAge.pdf", width =162.5*0.75, height=50, units = "mm")


