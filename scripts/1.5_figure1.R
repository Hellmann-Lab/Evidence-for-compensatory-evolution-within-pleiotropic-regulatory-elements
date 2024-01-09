
# Settings ----------------------------------------------------------------

library(DBI)
library(RMySQL)
library(tidyverse)
library(plyranges)

library(methods)
library(broom)
library(tidyverse)
library(patchwork)
library(figpatch)

#setwd("/data/share/htp/pleiotropy/paper_data/")


tissueColors <- c( "#9E0142" ,"#D53E4F", "#F46D43", "#FDAE61", "#FFD92F" ,"#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
tissueColors2 <- tissueColors
names(tissueColors2) <-  c("adrenal gland", "brain", "heart", "kidney", "large intestine", "lung", "muscle", "stomach", "thymus")


regionColors <- c("#9CA578","#3288BD")
names(regionColors) <-  c("Enhancer", "Promoter")

specificityColors <- c( "#A3753B", "#CC9B57", "#E7CF97", "#F8EDD0", "#F7F7F7", "#D2EEEA", "#99D7CE", "#5DACA5", "#33847E")
names(specificityColors) <-  c(1:9)

theme_set( theme_bw( base_size = 8 ) +
             theme(legend.title = element_blank(),
                   legend.position = "bottom",
                   legend.box = "horizontal",
                   legend.direction = "horizontal",
                   legend.key.size = unit(3,"mm")) )




# Sample Plot -------------------------------------------------------------

sampleSize<-readRDS("roadmap_DHS_summaries/srx_info.rds") %>% 
  dplyr::select(srx_id, donor_id, tissue, experiment) %>%
  mutate(tissue=gsub("_", " ", tissue)) %>% 
  filter(experiment == "DNase hypersensitivity")

sample_plot <-  sampleSize %>% filter(tissue %in% names(tissueColors2)) %>% 
  ggplot( aes(x = tissue, fill = tissue)) + geom_bar() + 
  scale_fill_manual(values=tissueColors2) +
  ylab("# of replicates") +
  coord_flip() + 
  theme(legend.position = "None",
        axis.title.y=element_blank(),
        plot.margin=unit(c(0.5,0,0.2,1),"cm"),
        panel.background = element_rect(fill = "transparent"))+
  scale_x_discrete(limits=rev)

sample_plot


# Peaks per tissue --------------------

dhs <- readRDS("roadmap_DHS_summaries/region_summary/jamm_region_info_full.rds")
dhs_long <- dhs %>% pivot_longer( cols = adrenal_gland:thymus, names_to = "tissue", values_to = "open") %>%  filter(open == 1)

DHS_to_gene<-readRDS("CRE_to_Gene/DHS_to_gene.rds")

pleio_tissue <- dhs_long %>%  inner_join( DHS_to_gene ) %>% 
  group_by(tissue,total, assignment) %>%  
  dplyr::summarise( n = dplyr::n(),
                    `% GC`= median(GC_content,na.rm = T),
                    `CpG obs/exp` = median(CpG_obs_exp, na.rm = T),
                     bp=median(total_length,na.rm=T) ) %>% 
  group_by(tissue,assignment) %>% 
  mutate(total_peaks= sum(n), frac_peaks = n/sum(n), `pleiotropic degree` = total) %>% 
  mutate(type = ifelse(assignment == "enh", "Enhancer","Promoter")) 



peaks_per_tissue <- dhs_long %>% dplyr::count(tissue) %>% 
  ggplot(  aes(y=tissue, x=n/1000, col = tissue)) +
  geom_point() +
  scale_color_manual(values=tissueColors) +
  xlab("# of CREs / 1000")+
  theme(legend.position = "None",
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin=unit(c(0.5,1,0.2,0),"cm"),
        panel.background = element_blank())+
  scale_y_discrete(limits=rev)



# Promoter Fraction plot --------------------------------------------------

promFrac<-DHS_to_gene %>% distinct(region_id,total,assignment,CGI ) %>% 
  dplyr::count(total,assignment,CGI) %>% 
  group_by(total) %>% mutate(all = sum(n)) %>% 
  mutate(total=as.factor(total)) %>% 
  mutate(assignment = gsub("enh","Enhancer",assignment),
         assignment = gsub("prom","Promoter",assignment))  %>%
  ggplot(aes(x=total, y=n/all*100, fill=assignment, alpha=CGI)) +
    geom_col() +
   scale_fill_manual(values = regionColors)+
   scale_alpha_manual(values=c(1,0.6)) +
  xlab("PD")+
  ylab('% of CREs')+
  theme(legend.position = "None")


# DHS - number plot -------------------------------------------------------

dhs_no <-  pleio_tissue %>% ggplot(aes(x=`pleiotropic degree`,y=n/1000,col=tissue,linetype=type)) +
  geom_line()+
  scale_color_manual(values = tissueColors)+
  scale_x_continuous(breaks=1:9)+ ylab("# of DHS/1,000")+
  xlab("PD")+
  scale_linetype_manual(values=c(3,1))+
  theme(legend.position = "None")

cre_size <-  pleio_tissue %>% 
 ggplot(aes(x=`pleiotropic degree`,y=bp,col=tissue,linetype=type)) +
  geom_line()+scale_color_manual(values = tissueColors)+
  guides(color="none")+
  scale_x_continuous(breaks=1:9)+ ylab("bp")+
  xlab("PD")+ 
  scale_linetype_manual(values=c(3,1))+
  theme(legend.position = c(0.3,0.75),
        legend.background = element_rect(fill='transparent'),
        legend.title = element_blank())

# Expression pleiotropy ---------------------------------------------------

summarized_expression_all<- readRDS("roadmap_expression_summaries/summarized_expression_allTissues.rds") %>% 
                      group_by(gene_id) %>% 
                      mutate( expr_pleio = length(tissue),
                              expr_pleio = factor(expr_pleio, levels=1:9 ))

# The promoter pleiotropy of a gene is determined by the promoters that are closest to a TSS for those, we then take the one with the max. pleiotropic degree.

expr_DHS_promoter<-DHS_to_gene %>% 
  group_by(gene_id,tissue) %>%
  filter( assignment == "prom" & distance == min(distance) ) %>%
  summarise( PD = max(total)) %>% 
  ungroup %>% 
  inner_join(summarized_expression_all) %>% 
  mutate( PD = factor(PD,level=1:9))


exprPD_vs_promPD<-expr_DHS_promoter %>% distinct(gene_id,PD,expr_pleio) %>%
  group_by(PD) %>% 
  mutate(tot =length(PD)) %>% 
  group_by(PD,expr_pleio,tot) %>% 
  summarize(cnt=length(expr_pleio)) %>%
  ggplot(aes(x=PD, y=cnt/tot, fill=expr_pleio))+
  geom_bar(stat="identity")+ scale_fill_manual(values = specificityColors)+
  ylab("fraction of promoters")+
  guides(fill=guide_legend(title="Expression PD"))+
  xlab("Promoter PD")+theme(legend.title = element_text(size=7) )
 

exprPD_vs_promPD

# impact of element composition on expression -----------------------------

data<- DHS_to_gene %>% 
  mutate(total = factor(total,level=1:9)) %>% 
  left_join(summarized_expression_all %>% 
              dplyr::select(tissue,gene_id,log2_mean_expression)) 

# need to first run this one
#system("Rscript 1.4_run_model_permutations_expression_pleiotropy.R")

pred<-readRDS("roadmap_expression_summaries/complex_mixed_model_boot.RDS")


enh_expr_plot<-pred$coef %>% 
  separate(term, into=c("type","CGI","PD"), sep="_",convert = T) %>% 
  mutate( type = ifelse(type == "enh", "Enhancer","Promoter") ) %>% 
  ggplot(aes(x=PD,y=estimate,col=type, alpha=CGI)) + 
  geom_line( )+
  #geom_point(aes( shape=CGI, 
                  # ymin = CI_lower, 
                  # ymax = CI_upper, 
                  # size=0.25 ))+
  scale_color_manual(values = regionColors)+
  scale_alpha_manual(values=c(1,0.5))+
  scale_x_continuous(breaks=1:9)+
  scale_linetype_manual(values=c(3,1))+
  ylab(expression(beta["scaled"]))

enh_expr_plot





# Assemble Figure Panels --------------------------------------------------

fig_1A<-fig("figures/fig1A_short.pdf", 
            b_margin = ggplot2::margin(-10,0,-10,-10), link_dim = T)+
  theme(plot.margin = margin(0,0,0,0))
  
sample_plot1<-sample_plot+theme(plot.margin = margin(0,0,0,6))
dhs_no1<-dhs_no+theme(plot.margin = margin(0,3,0,6))

p1.1 <- ( sample_plot1 + peaks_per_tissue + plot_layout(widths = c(2, 1.5)) +
            theme(plot.margin = margin(0,0,0,0))) 
p1.2 <- ( dhs_no1 + cre_size+theme(plot.margin = margin(0,0,0,0))) 
# p2 <-  promFrac+theme(legend.position = "bottom") + 
#   enh_expr_plot+theme(legend.position = "none") +  
#   exprPD_vs_promPD+theme(legend.position = "bottom")
                                               
p2 <-  promFrac + enh_expr_plot + exprPD_vs_promPD & theme(legend.position = "bottom", 
                                                           legend.box = "vertical", 
                                                           legend.spacing.y = unit(2, "mm"))

top<-p1.1/p1.2
top[[2]][[1]]<-top[[2]][[1]]+theme(axis.title.y = element_text(vjust = -15))
p1<-fig_1A+top+plot_layout(widths = c(0.7,1))+theme(plot.margin = margin(0,0,0,0))

fig1_new<-p1/p2/guide_area()+
  plot_layout(heights=c(1.15,0.5,0.2))+
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size=12,face = "bold"))

ggsave(file="figures/fig1_v3.pdf", fig1_new, width =162.5, height=138, units = "mm")


# extra
p1.1 <- ( sample_plot1 + peaks_per_tissue + plot_layout(widths = c(2, 1.5)) +
            theme(plot.margin = margin(7,7,7,7))) 
ggsave("figures/Fig1BC.pdf", p1.1, width =162.5*0.55, height=138*0.29, units = "mm")




# Supplementary figure ----------------------------------------------------

#maybe add R2
enh_expr_suppl<-pred$coef %>% separate(term, into=c("type","CGI","PD"), sep="_",convert = T) %>%
  mutate( type = ifelse(type == "enh", "Enhancer","Promoter") ) %>% 
  ggplot(aes(x=PD,y=estimate)) + 
  geom_line(col="darkgrey")+
  geom_ribbon(aes(ymin=perc_5, ymax=perc_95),col="darkgrey",fill="darkgrey",alpha=0.4)+
  geom_point(aes(col=as.factor(PD)),
             size=2)+
  scale_color_manual(values = specificityColors) +
  scale_x_continuous(breaks=1:9)+
  facet_grid(type~CGI,scales = "free")+
  ylab(expression(beta~"weights"))+ guides(colour = "none")+
  theme(legend.title = element_blank(),
        legend.position = "None",
        legend.background = element_rect(fill='transparent'))

enh_expr_suppl

