########
########
# Set-up environment
########
########

# library
library(ggplot2)
library(readr)
library(tidyverse)
library(data.table)
library(tools)
library(ggpubr)
library(gridExtra)

options(scipen=999)
`%notin%` <- Negate(`%in%`)

#Define custom color palette
mycolors <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',  '#008080', '#e6194b', '#ffe119', 
              '#aaffc3', '#808000', '#3cb44b', '#e6beff', '#9a6324', '#fffac8', '#bcf60c','#800000',
              '#fabebe','#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')


########
########
# Import and wrangling data
########
########

#Load in metadata
global_skin_metadata <- read_csv("global_skin_metadata.csv")



####
#Load in bracken concatenated output and extract total reads and fraction C. acnes
bracken <- read_csv("merged_bracken_output.csv",
                    col_types = cols(taxonomy_id = col_skip(),
                                     taxonomy_lvl = col_skip(),
                                     kraken_assigned_reads = col_skip(),
                                     added_reads = col_skip(),
                                     fraction_total_reads = col_skip()) )
head(bracken)
names(bracken)[1] = "taxa"
bracken = bracken %>% 
  filter(!grepl("Homo sapiens", taxa)) %>% 
  group_by(filenames) %>%
  mutate(SampleName = gsub(".bracken","",filenames),
         tot_reads = sum(new_est_reads)) %>%
  ungroup() %>%
  select(SampleName, tot_reads, taxa, new_est_reads) %>%
  filter(!is.na(taxa))
head(bracken)


#Simplified dataset
bracken_melt = left_join(
  bracken %>%
  mutate(ColorCode = ifelse(grepl("Cutibacterium", taxa), "Cuti",
                            ifelse(grepl("Staphylococcus", taxa), "Staph",
                                   ifelse(grepl("Streptococcus", taxa), "Strep",
                                          ifelse(grepl("Corynebacterium", taxa), "Coryne",
                                          "Other"))))) %>%
  group_by(SampleName, tot_reads, ColorCode) %>%
  summarise(new_est_reads = sum(new_est_reads)) %>%
  mutate(fraction_total_reads = new_est_reads/tot_reads),
  global_skin_metadata %>%
    select(SampleName, BioProject, StudyName, Region) )


bracken_dcast = bracken_melt %>%
  reshape2::dcast(., StudyName + BioProject + Region + SampleName + tot_reads ~ ColorCode, fun.aggregate = sum, value.var = "new_est_reads") 

rm(bracken)

# #Plot

#Total reads v C. acnes boxplot
study_reads <- ggplot() +
  geom_violin(data = bracken_dcast, aes(y=StudyName,x=log10(tot_reads))) +
  geom_boxplot(data = bracken_dcast, aes(y=StudyName,x=log10(Cuti), alpha=0.1), width=0.1, color="blue") +
  theme_classic() + theme(legend.position = "none") +
  ggtitle("Total (black) Cutibacterium (blue)")  +
  xlab("log10(reads)") 
# Save plot to a PDF file
ggsave(
  filename = "reads_total_cacnes.pdf",
  plot = study_reads ,
  device = "pdf",
  width = 6,    # adjust width (in inches) as needed
  height = 6,    # adjust height (in inches) as needed
  units = "in"
)


#Save and export
write.csv(bracken_dcast, "bracken_dcast-taxa.csv", row.names = F, quote = F)


##Plot

bracken_taxa_plot <-
  ggplot(bracken_melt,
       aes(x=SampleName, y=fraction_total_reads, fill=ColorCode)) +
  geom_col(alpha=0.65, width = 0.7)  +
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "white"),
        panel.spacing = unit(0.8, "lines")) +
  facet_grid(BioProject ~ . , scales="free_y", space="free_y") +
  scale_fill_manual(values=c("orange","green4","grey","blue4","purple")) +
  coord_flip()


ggsave(
  filename = "bracken.pdf",
  plot = bracken_taxa_plot,
  device = "pdf",
  width = 4,    # adjust width (in inches) as needed
  height = 30,    # adjust height (in inches) as needed
  units = "in"
)

