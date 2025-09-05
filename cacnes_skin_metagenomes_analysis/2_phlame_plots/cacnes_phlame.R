



#Libraries
library(tidyverse)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(patchwork)


options(scipen=999)
`%notin%` <- Negate(`%in%`)


########




#########

#Import in PHLAME data set, with prevalence and relative frequency info

#Import data
#phage mapping data
samtools_coverage_stats_cutiphage_update <- read_csv("samtools_coverage_stats-cutiphage-update.csv")
names(samtools_coverage_stats_cutiphage_update)[1] = "ref"
phage = samtools_coverage_stats_cutiphage_update %>%
  mutate(SampleName = gsub("_ref_Pacnes_phage_NC028967_cov.tsv","", filename),
         vir_depth=meandepth,
         vir_cov=coverage,
         vir_reads=numreads) %>%
  select(vir_depth,vir_cov,vir_reads,SampleName) 

#bacteria mapping data
samtools_coverage_stats_cacnes <- read_csv("samtools_coverage_stats-cacnes.csv")
names(samtools_coverage_stats_cacnes)[1] = "ref"
cacnes = samtools_coverage_stats_cacnes %>%
  mutate(SampleName = gsub("_ref_Pacnes_C1_aligned.sorted.bam","", filename),
         bac_depth=meandepth,
         bac_cov=coverage,
         bac_reads=numreads) %>%
  select(bac_depth,bac_cov,bac_reads,SampleName)

#####
#Import metadata
global_skin_metadata <- read_csv("global_skin_metadata.csv")
#Read lengths by study
global_skin_metadata_readlength <- read_csv("global_skin_metadata-readlength.csv")
View(global_skin_metadata_readlength)

#####
#phlame results
# Cacnes_global_phylogroup_frequencies <- read_csv("Cacnes_global_phylogroup_frequencies.csv")
Cacnes_global_subphylogroup_frequencies <- read_csv("Cacnes_global_subphylogroup_frequencies.csv")
######
#import brakcen results
bracken_dcast_cacnes <- read_csv("bracken_dcast-cacnes.csv")

#Merge into final df
vmr = left_join(left_join(left_join(left_join(left_join(phage, 
                                                        cacnes), 
                                              Cacnes_global_subphylogroup_frequencies), 
                                    global_skin_metadata),
                          bracken_dcast_cacnes) %>% 
                  group_by(SampleName) %>%
                  mutate(vmr = vir_depth/bac_depth,
                         vmr = ifelse(vir_cov >= 10, vmr, 0)) %>%
                  select(-c(Instrument,LibraryLayout,LibrarySelection,Site,`Public data?`,`Run Number`,AvgSpotLen,City,Country,`United States, Europe or Canada`,Ethnicity,Age, BioProject))
                ,
                global_skin_metadata_readlength)
View(vmr)


#####
#Plot prevalence & relative abundance of each subphylo
# subphylo_order = c("A.1","A.2","A.3","B.1","C.1","F.1","F.2","F.3","E.1","D.1","H.1","H.2","H.3","K.1.1","K.1.2","K.2","K.3")
subphylo_order = c("A.3","A.1","D.1","F.3","H.2","K.2","E.1","K.1.1","K.1.2","F.2","A.2","K.3","C.1","F.1","H.1","B.1","H.3")

#
#Second prevalence
g<-ggplot(vmr %>% 
         filter(tot_reads >=1000000 & bac_depth >= 5) %>% 
         select(-c(vir_depth, vir_cov, vir_reads, bac_depth, bac_cov, bac_reads, tot_reads, bracken_cacnes, ApproxAge, Region, Sex, StudyName, `Acne?`, BioProject, NumberIndividuals, ReadLength, `L.1`)) %>% 
         group_by(SampleName, vmr) %>% 
         reshape2::melt(id.vars = c("SampleName", "vmr")) %>% 
         select(SampleName, variable, value) %>% 
         filter(value >= 0.01) %>% 
         group_by(variable) %>% 
         summarise(prevalence = n()/471),
       aes(y=factor(variable, levels=rev(subphylo_order)), 
           x= prevalence, 
           fill=ifelse(variable %in% c("D.1","H.1","H.2","H.3","K.1.1","K.1.2","K.2","K.3"), "PDE", "noPDE"))) +
  geom_col(alpha=0.8)  +
  labs(y="Subphylogroup", x="prevalence") +
  theme_pubr() + xlim(0,0.75) +
  scale_fill_manual(values=c("grey50","olivedrab3")) +
  theme(legend.position = "none") 

q<-ggplot(vmr %>%
         filter(tot_reads >=1000000 & bac_depth >= 5) %>%
         select(-c(vir_depth, vir_cov, vir_reads, bac_depth, bac_cov, bac_reads, tot_reads, bracken_cacnes, ApproxAge, Region, Sex, StudyName, `Acne?`, BioProject, NumberIndividuals, ReadLength, `L.1`)) %>%
         group_by(SampleName, vmr) %>%
         reshape2::melt(id.vars = c("SampleName", "vmr")),
       aes(y=factor(variable, levels=rev(subphylo_order)), x= value,
           color=ifelse(variable %in% c("D.1","H.1","H.2","H.3","K.1.1","K.1.2","K.2","K.3"), "PDE", "noPDE"))) +
  geom_point(alpha=0.3, size=3) +
  stat_summary(fun.y = mean,  width=0.75, color="black", size=0.3, geom="crossbar", alpha=0.5) +
  labs(y="Subphylogroup", x="phlame relative abundance") +
  theme_pubr()  +
  scale_color_manual(values=c("grey50","olivedrab3")) +
  theme(legend.position = "none")


grid.arrange(g,q, ncol=2, widths=c(1/3,2/3))


# Save plot to a PDF file
ggsave(
  filename = "phalme_subphylo_rel_abund.pdf",
  plot = grid.arrange(g,q, ncol=2, widths=c(1/3,2/3)),
  device = "pdf",
  width = 11,    # adjust width (in inches) as needed
  height = 6,    # adjust height (in inches) as needed
  units = "in"
)
########





###
#
#Second prevalence
h<-ggplot(vmr %>% 
         filter(tot_reads >=1000000 & bac_depth >= 5) %>% 
         select(-c(vir_depth, vir_cov, vir_reads, bac_depth, bac_cov, bac_reads, tot_reads, 
                   bracken_cacnes, ApproxAge, Region, Sex, StudyName, `Acne?`, BioProject,
                   NumberIndividuals, ReadLength, `L.1`)) %>% 
         group_by(SampleName, vmr) %>% 
         reshape2::melt(id.vars = c("SampleName", "vmr")) %>% 
         select(SampleName, variable, value) %>% 
         filter(value >= 0.01) %>% 
         group_by(variable) %>% 
         summarise(prevalence = n()/471) %>%
         mutate(myfill=ifelse(variable %in% c("D.1","H.1","H.2","H.3","K.1.1","K.1.2","K.2","K.3"), "PDE", "noPDE")),
       aes(x=factor(myfill), 
           y= prevalence, 
           fill=myfill, 
           color=myfill)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions="all", alpha=0.5 ) +
  stat_compare_means(label.y=0.8) +
  scale_fill_manual(values = c("grey50","olivedrab3"))+   
  scale_color_manual(values = c("grey50","olivedrab3")) +
  stat_summary(fun.y = mean, geom = "crossbar", width=0.3, color="black") +
  theme_pubr() 
#

########
r<-ggplot(vmr %>% 
         filter(tot_reads >=1000000  & bac_depth >= 5) %>% 
         select(-c(vir_depth, vir_cov, vir_reads, bac_depth, bac_cov, bac_reads, tot_reads, 
                   bracken_cacnes, ApproxAge, Region, Sex, StudyName, `Acne?`, 
                   BioProject, NumberIndividuals, ReadLength, `L.1`)) %>% 
         group_by(SampleName, vmr) %>% 
         reshape2::melt(id.vars = c("SampleName", "vmr")) %>% 
         group_by(variable) %>% 
         reframe(freq=mean(value, na.rm=T),
                 myfill=ifelse(variable %in% c("D.1","H.1","H.2","H.3","K.1.1","K.1.2","K.2","K.3"), "PDE", "noPDE")) %>% 
         unique(),
       aes(x=factor(myfill), y=freq, color=myfill, fill=myfill)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions="all", alpha=0.5 ) +
  stat_compare_means(label.y=0.25) +
  scale_fill_manual(values = c("grey50","olivedrab3"))+   
  scale_color_manual(values = c("grey50","olivedrab3")) +
  stat_summary(fun.y = mean, geom = "crossbar", width=0.3, color="black") +
  theme_pubr()
#


grid.arrange(h,r, ncol=2)


# Save plot to a PDF file
ggsave(
  filename = "phalme_subphylo_rel_abund_pvalues.pdf",
  plot = grid.arrange(h,r, ncol=2),
  device = "pdf",
  width = 8,    # adjust width (in inches) as needed
  height = 4,    # adjust height (in inches) as needed
  units = "in"
)

#########

