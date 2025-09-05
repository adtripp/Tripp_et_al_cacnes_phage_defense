
#
#Libraries
library(tidyverse)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(patchwork)
library(readr)

options(scipen=999)
`%notin%` <- Negate(`%in%`)


#Import data

#all genes pastml
cacnes_pastmlout <- read_csv("cacnes_pastmlout.csv", 
                             col_types = cols(...1 = col_skip()))
View(cacnes_pastmlout)

#names of pdes
cacnes_pde_pastmlout_filt <- read_csv("cacnes_pde_pastmlout-filt.csv")
View(cacnes_pde_pastmlout_filt)

#Categorize genes as core, or noncore based on lineage prevalence
#Above 95% is core
#
cacnes_gene_presence_absence_lineage <- read_csv("cacnes_gene_presence_absence_lineage.csv")
View(cacnes_gene_presence_absence_lineage)
gene_categories = cacnes_gene_presence_absence_lineage %>%
  group_by(Lineage) %>%
  reshape2::melt() %>%
  ungroup() %>%
  mutate(gene = variable) %>%
  group_by(gene) %>%
  filter(value == 1) %>%
  summarise(gene_cnt = sum(value, na.rm=F),
            gene_carriage = sum(value, na.rm=F)/(nrow(cacnes_gene_presence_absence_lineage)) ) %>%
  ungroup() %>%
  mutate(pangenome = ifelse(gene_carriage < 0.15, "cloud_gene",
                            ifelse(gene_carriage < 0.95 & gene_carriage >= 0.15, "shell_gene", NA))) 
  
  
#Add the gain/loss info
gainloss = rbind(
  left_join( gene_categories,
             cacnes_pastmlout) %>%
    mutate(event_cnt = gains + losses) %>%
    select(gene, gains, losses, count, event_cnt, pangenome),
  cacnes_pde_pastmlout_filt %>%
    mutate(event_cnt = gains + losses,
           pangenome = "pdes") )
View(gainloss)


#Plot boxplot
gene_order = c("shell_gene","cloud_gene","pdes")

ggplot(gainloss, aes(x=factor(pangenome, levels=gene_order) ,y=event_cnt)) + 
  geom_jitter(alpha=0.3,width=0.05) +
  stat_compare_means(comparisons = list(c("pdes","shell_gene"),
                                        c("pdes","cloud_gene"),
                                        c("shell_gene","cloud_gene"))) +
  theme_pubr() +
  stat_summary(fun.y = mean, geom = "crossbar", width=0.4, color="blue")+
  labs(title="C.acnes gain/loss") +
  theme(axis.title.x = element_blank())

#

