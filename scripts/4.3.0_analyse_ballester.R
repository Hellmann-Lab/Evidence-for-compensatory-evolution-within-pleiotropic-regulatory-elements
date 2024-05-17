
# Ballester data analysis:
# take liver CREs from macaque, translate to hg38 and then to hg18
# overlap with TFBS (translate to the same coordinate system)
# check if high PD is more conserved (more 1s)


setwd("/data/share/htp/pleiotropy/paper_data/")
library(tidyverse)
library(data.table)
library(plyranges)
library(readxl)

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


# Strategy:
# 1) Select collapsed_ids for CREs in macaque that are open in the liver in any of the 5 species 
# 2) Identify PD1: the ones that are PD1 in human (i.e. not in our data) AND PD1 in macaque in Roller et al
# 3) Identify PD9: PD9 in human AND PD4 in macaque in Roller et al
# 4) Split by TF etc. Calculate phylogenetic age as for CRE activity


# Get liver CRE coordinates (with PD) -------------------------------------

# need to select collapsed_region_ids from all livers in 4 species that 

# these are macaque collapsed CREs
hummac<-readRDS("Roller2021/chipseq_OL_humJamm_mac_collapsed.rds") %>% 
  dplyr::select(region_id, collapsed_region_id) %>% 
  left_join(readRDS("Roller2021/chipseq_macregs_collapsed.rds")) %>% 
  right_join(readRDS("roadmap_DHS_summaries/region_summary/jamm_region_info_full.rds") %>%
              dplyr::select(region_id, total))


# Select PD9s (us) and PD4s (Roller)
hummac_pleiotropic<-hummac %>% filter(total_mmul == 4 & total == 9)

# Select liver PD1s (Roller) not present in our data
macregs_PD1_liver<-readRDS("Roller2021/chipseq_macregs_collapsed.rds") %>% 
  dplyr::filter(tissues_mmul == "Liver", total_mmul ==1) %>% 
  anti_join(hummac)



# Combine liver-specific and pleiotropic CREs
liverCREs_mmul10<-bind_rows(hummac_pleiotropic %>% 
                              dplyr::select(seqnames, start, end, collapsed_region_id) %>% 
                              dplyr::mutate(type="pleiotropic"),
                            macregs_PD1_liver %>% 
                              dplyr::select(seqnames, start, end, collapsed_region_id) %>% 
                              dplyr::mutate(type="liver-specific"))

dim(liverCREs_mmul10)




# Translate TFBS from hg18 to mmul10 --------------------------------------

tfbs_liver_hg18 <- read_excel("Ballester2014/elife-02626-fig2-data1-v1.xlsx") %>%
  dplyr::transmute(seqnames=chr, start=start_hg18, end=stop_hg18, CRM_or_singleton_name,
                   CRM_status, shared_status_hsap_mmul_mmus_rnor_cfam, region_id = 1:nrow(.))


tfbs_hg19<-translate_jamm(chain_file = "liftOvers/hg18ToHg19.over.chain",
                                coordinate_file = as_granges(tfbs_liver_hg18), 
                                extend = 50,
                                reverse_chain_file = "liftOvers/hg19ToHg18.over.chain")


tfbs_mmul10<-translate_jamm(chain_file = "liftOvers/hg19ToRheMac10.over.chain",
                                    coordinate_file = as_granges(tfbs_hg19), 
                                    extend = 50,
                                    reverse_chain_file = "liftOvers/rheMac10ToHg19.over.chain") %>%
  left_join(tfbs_liver_hg18 %>% 
              dplyr::select(region_id, CRM_or_singleton_name,
                            CRM_status, shared_status_hsap_mmul_mmus_rnor_cfam)) %>% 
  dplyr::select(-region_id)






# Annotate liver CREs where TFBS are within -------------------------------

liverPD_TFBS<-join_overlap_left(as_granges(liverCREs_mmul10), as_granges(tfbs_mmul10)) %>%
  filter(!is.na(shared_status_hsap_mmul_mmus_rnor_cfam)) %>% 
  as_tibble() %>% 
  rowwise() %>% 
  dplyr::mutate(
    n_TFspecies = sum(as.numeric(strsplit(as.character(shared_status_hsap_mmul_mmus_rnor_cfam), "")[[1]])),
    macaque=ifelse(substring(as.character(shared_status_hsap_mmul_mmus_rnor_cfam), 2, 2) == "1",T,F),
    TF = word(CRM_or_singleton_name, 1,1,":"))

saveRDS(liverPD_TFBS, "Ballester2014/summarized_PD1_PD9_TFBS_mmul10.rds")



pProps<-ggplot(liverPD_TFBS %>% 
                 mutate(PD = ifelse(type=="pleiotropic","9","1")),
               aes(x = PD, fill = as.factor(n_TFspecies)))+
  geom_bar(position="fill")+
  theme_bw()+
  scale_fill_manual(values=c("1" = "#ffcdb2",
                             "2" = "#ffb4a2",
                             "3" = "#e5989b",
                             "4" = "#b5838d",
                             "5" = "#6d6875"),
                    name="# of species\nwith conserved\nTF repertoire")+
  xlab("PD") +
  ylab("Fraction")+
  basic_theme_ins2+
  theme(axis.text = element_text(color="black"))
pProps

ggsave("Roller2021/figures/nSpecies_vs_TFrepertoire.png", height=3, width=4.1)
ggsave("Roller2021/figures/nSpecies_vs_TFrepertoire.pdf", height=3, width=4.1)




# Calculate TFBS repertoire phylogenetic age  -----------------------------

liverPD_TFBS_long<-readRDS("Ballester2014/summarized_PD1_PD9_TFBS_mmul10.rds") %>% 
  dplyr::transmute(collapsed_region_id, CRM_or_singleton_name,
                   status = shared_status_hsap_mmul_mmus_rnor_cfam, TF, type) %>% 
  dplyr::mutate(
    Homo_sapiens=ifelse(substring(as.character(status), 1, 1) == "1",1,0),
    Macaca_mulatta=ifelse(substring(as.character(status), 2, 2) == "1",1,0),
    Mus_musculus=ifelse(substring(as.character(status), 3, 3) == "1",1,0),
    Rattus_norvegicus=ifelse(substring(as.character(status), 4, 4) == "1",1,0),
    Canis_lupus=ifelse(substring(as.character(status), 5, 5) == "1",1,0)
  ) %>% 
  pivot_longer(6:10, names_to = "latin") %>% 
  filter(value == 1) %>% 
  dplyr::select(-value)


# generate a tree
BallesterSpeciesTree<-drop.tip(mammaltree, mammaltree$tip.label[!mammaltree$tip.label %in% c("Homo_sapiens","Macaca_mulatta","Canis_lupus", "Mus_musculus", "Rattus_norvegicus")])

BallesterSpeciesTree_plot<-BallesterSpeciesTree
BallesterSpeciesTree_plot$tip.label<- case_when(BallesterSpeciesTree_plot$tip.label == "Homo_sapiens" ~ "Human",
                                                BallesterSpeciesTree_plot$tip.label == "Macaca_mulatta" ~ "Macaque",
                                                BallesterSpeciesTree_plot$tip.label == "Felis_catus" ~ "Cat",
                                                BallesterSpeciesTree_plot$tip.label == "Canis_lupus" ~ "Dog", 
                                                BallesterSpeciesTree_plot$tip.label == "Equus_ferus" ~ "Horse", 
                                                BallesterSpeciesTree_plot$tip.label == "Mus_musculus" ~ "Mouse", 
                                                BallesterSpeciesTree_plot$tip.label == "Rattus_norvegicus" ~ "Rat",
                                                BallesterSpeciesTree_plot$tip.label == "Oryctolagus_cuniculus" ~ "Rabbit", 
                                                BallesterSpeciesTree_plot$tip.label == "Sus_scrofa" ~ "Pig", 
                                                BallesterSpeciesTree_plot$tip.label == "Monodelphis_domestica" ~ "Opossum")

ggtree(BallesterSpeciesTree_plot, color="grey60", lwd=0.75)+geom_tiplab()+xlim(0,150) 

mamtree2<-ggtree(BallesterSpeciesTree_plot, lwd=0.3)+
  geom_tiplab(fontface="italic", size=2.5, geom = "text")+ 
  xlim(0,130) +
  theme(plot.margin = unit(c(0, 0, 0,0), "cm")) 
mamtree2

ggsave("Roller2021/figures/ballester_species_tree.png", height = 5*0.70, width=7*0.70)



# calculate relative tree length
treeLengthsBallester<-liverPD_TFBS_long %>%
  distinct(collapsed_region_id, TF, type, status) %>% 
  rowwise() %>%
  dplyr::mutate(
    tree_length=calculate_tree_length(df = liverPD_TFBS_long, 
                                      regid = collapsed_region_id,
                                      regid_column = collapsed_region_id,
                                      tree = BallesterSpeciesTree)) %>%
  ungroup() %>%
  dplyr::mutate(max_tree_length=max(tree_length)) %>%
  rowwise() %>%
  mutate(rel_tree_length=tree_length/max_tree_length) %>%
  ungroup() 

# recheck where I loose like 500 TF-CRE pairs rel to liverPD_TFBS
saveRDS(treeLengthsBallester, "Ballester2014/PhyloTreeLengths.rds")


# plot
ggplot(treeLengthsBallester %>% 
         mutate(PD = ifelse(type=="pleiotropic","9","1")), 
       aes(x=PD, y=rel_tree_length, fill=PD))+
  geom_boxplot(notch = T, width=0.6)+
  theme_bw()+
  xlab("PD")+
  ylab(expression(paste("TF repertoire age (relative ", lambda, ")")))+
  scale_fill_manual(values=c(specificityColors[[1]],specificityColors[[9]]))+
  theme(legend.position = "none",
        axis.text = element_text(color="black"))

ggsave("Roller2021/figures/ballester_evoAge.png", height=3, width=4)


wilcox.test(
  treeLengthsBallester$rel_tree_length[treeLengthsBallester$type=="pleiotropic"], 
  treeLengthsBallester$rel_tree_length[treeLengthsBallester$type=="liver-specific"])

tflist<-list()
for (tf in unique(treeLengthsBallester$TF)){
  dff<-treeLengthsBallester %>% filter(TF == tf)
  tflist[[tf]]<-wilcox.test(
    dff$rel_tree_length[dff$type=="pleiotropic"], 
    dff$rel_tree_length[dff$type=="liver-specific"]) %>% broom::tidy()
}
TF_df<-tflist %>% bind_rows(.id="TF") %>% mutate(padj = p.adjust(p.value, method = "BH"))


p_age_all<-treeLengthsBallester %>% 
  group_by(type) %>% 
  summarize(mean = mean(rel_tree_length),
            sem = sd(rel_tree_length)/sqrt(length(type))) %>% 
  dplyr::mutate(#TF = factor(TF, levels=c("CEBPA","FOXA1","HNF4A","HNF6","CRM")),
                PD = ifelse(type=="pleiotropic","9","1")) %>% 
  ggplot(aes(x=PD, y=mean, color=PD))+
  geom_point()+
  geom_errorbar(aes(ymin = mean-sem, ymax=mean+sem), width=0.2)+
  theme_bw()+
  xlab("PD")+
  ylab(expression(paste("TF binding conservation (", lambda[i] / lambda, ")")))+
  scale_color_manual(values=c(specificityColors[[1]],specificityColors[[9]]))+
  basic_theme_ins2+
  theme(legend.position = "none",
        axis.text = element_text(color="black"))#+
  #facet_grid(.~TF)



pTFage<-treeLengthsBallester %>% 
  group_by(type, TF) %>% 
  summarize(mean = mean(rel_tree_length),
            sem = sd(rel_tree_length)/sqrt(length(type))) %>% 
  dplyr::mutate(TF = factor(TF, levels=c("CEBPA","FOXA1","HNF4A","HNF6","CRM")),
                PD = ifelse(type=="pleiotropic","9","1")) %>% 
  ggplot(aes(x=PD, y=mean, color=PD))+
  geom_point()+
  geom_errorbar(aes(ymin = mean-sem, ymax=mean+sem), width=0.2)+
  theme_bw()+
  xlab("PD")+
  ylab(expression(paste("TF binding conservation (", lambda[i] / lambda, ")")))+
  scale_color_manual(values=c(specificityColors[[1]],specificityColors[[9]]))+
  basic_theme_ins2+
  theme(legend.position = "none",
        axis.text = element_text(color="black"))+
  facet_grid(.~TF)
pTFage

ggsave("Roller2021/figures/ballester_evoAge_perTF.png", height=3, width=10)




ptop<-plot_grid(mamtree2,pProps, p_age_all, ncol=3,labels=c("A","B","C"), label_size = 12, scale=c(0.9,1,1), rel_widths = c(0.9,1,0.85))
psuppTF<-plot_grid(ptop,pTFage, ncol=1, labels=c("","D"), label_size=12)#, rel_heights=c(1.2,1))
psuppTF
ggsave("figures/TFBS_evoAge.pdf", width =162.5, height=113, units = "mm")

