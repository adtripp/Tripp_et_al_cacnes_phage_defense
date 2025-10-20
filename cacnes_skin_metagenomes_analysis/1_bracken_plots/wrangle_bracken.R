################################################################################
################################################################################
# ENVIRONMENT SETUP
################################################################################
################################################################################

# Load required libraries for data manipulation and visualization
library(ggplot2)      # Grammar of graphics plotting system
library(readr)        # Fast CSV reading/writing
library(tidyverse)    # Collection of data science packages (includes dplyr, tidyr, etc.)
library(data.table)   # High-performance data manipulation
library(tools)        # File and path utilities
library(ggpubr)       # Publication-ready plots
library(gridExtra)    # Arrange multiple grid-based plots

# Disable scientific notation for clearer numeric output
options(scipen=999)

# Create custom 'not in' operator as opposite of %in%
`%notin%` <- Negate(`%in%`)

# Define custom color palette for consistent visualization across plots
# 22 distinct colors for categorical data representation
mycolors <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',  '#008080', 
              '#e6194b', '#ffe119', '#aaffc3', '#808000', '#3cb44b', '#e6beff', 
              '#9a6324', '#fffac8', '#bcf60c','#800000', '#fabebe','#ffd8b1', 
              '#000075', '#808080', '#ffffff', '#000000')


################################################################################
################################################################################
# DATA IMPORT AND WRANGLING
################################################################################
################################################################################

# Load sample metadata containing study information and sample characteristics
global_skin_metadata <- read_csv("global_skin_metadata.csv")


################################################################################
# PROCESS BRACKEN OUTPUT DATA
################################################################################
# Bracken is a Bayesian reestimation tool for abundance estimation
# Load concatenated bracken output containing taxonomic classification results
# and total read counts

bracken <- read_csv("merged_bracken_output.csv")
head(bracken)  # Preview first few rows


################################################################################
# CREATE SIMPLIFIED DATASET WITH TAXONOMIC GROUPING
################################################################################

# Join bracken data with metadata and categorize taxa into major groups
bracken_melt = left_join(
  # Categorize taxa into color-coded groups based on genus
  bracken %>%
  mutate(ColorCode = ifelse(grepl("Cutibacterium", taxa), "Cuti",
                            ifelse(grepl("Staphylococcus", taxa), "Staph",
                                   ifelse(grepl("Streptococcus", taxa), "Strep",
                                          ifelse(grepl("Corynebacterium", taxa), "Coryne",
                                          "Other"))))) %>%
  # Aggregate read counts by sample and taxonomic group
  group_by(SampleName, tot_reads, ColorCode) %>%
  summarise(new_est_reads = sum(new_est_reads)) %>%
  # Calculate fraction of total reads for each taxonomic group
  mutate(fraction_total_reads = new_est_reads/tot_reads),
  
  # Join with metadata to add study context
  global_skin_metadata %>%
    select(SampleName, BioProject, StudyName, Region) )


################################################################################
# RESHAPE DATA TO WIDE FORMAT
################################################################################
# Convert from long format (one row per sample-taxa combination) to 
# wide format (one row per sample with taxa as columns)

bracken_dcast = bracken_melt %>%
  reshape2::dcast(., StudyName + BioProject + Region + SampleName + tot_reads ~ ColorCode, 
                  fun.aggregate = sum, 
                  value.var = "new_est_reads") 

# Remove original bracken object to free memory
rm(bracken)


################################################################################
################################################################################
# VISUALIZATION: TOTAL READS VS CUTIBACTERIUM READS
################################################################################
################################################################################

# Create combined violin and boxplot comparing total reads to Cutibacterium reads
# across different studies
study_reads <- ggplot() +
  # Violin plot showing distribution of total reads per study
  geom_violin(data = bracken_dcast, aes(y=StudyName, x=log10(tot_reads))) +
  # Overlaid boxplot showing Cutibacterium reads distribution (in blue)
  geom_boxplot(data = bracken_dcast, aes(y=StudyName, x=log10(Cuti), alpha=0.1), 
               width=0.1, color="blue") +
  theme_classic() + 
  theme(legend.position = "none") +  # Hide legend
  ggtitle("Total (black) Cutibacterium (blue)")  +
  xlab("log10(reads)") 

# Export plot to PDF
ggsave(
  filename = "reads_total_cacnes.pdf",
  plot = study_reads,
  device = "pdf",
  width = 6,    # Width in inches
  height = 6,   # Height in inches
  units = "in"
)


################################################################################
# EXPORT PROCESSED DATA
################################################################################

# Save wide-format bracken data with taxonomic groupings to CSV
write.csv(bracken_dcast, "bracken_dcast-taxa.csv", row.names = F, quote = F)


################################################################################
################################################################################
# VISUALIZATION: TAXONOMIC COMPOSITION STACKED BARPLOT
################################################################################
################################################################################

# Create stacked bar chart showing relative abundance of major taxa groups
# across all samples, faceted by BioProject
bracken_taxa_plot <-
  ggplot(bracken_melt,
       aes(x=SampleName, y=fraction_total_reads, fill=ColorCode)) +
  # Stacked columns with slight transparency
  geom_col(alpha=0.65, width = 0.7)  +
  theme_classic()+
  theme(
    axis.text.y = element_blank(),     # Hide y-axis labels (too many samples)
    axis.ticks.y = element_blank(),    # Hide y-axis ticks
    strip.background = element_rect(fill = "white"),  # White facet labels
    panel.spacing = unit(0.8, "lines")  # Space between faceted panels
  ) +
  # Separate panels for each BioProject with independent y-axes
  facet_grid(BioProject ~ . , scales="free_y", space="free_y") +
  # Apply custom colors to taxonomic groups
  # Cuti=blue4, Coryne=orange, Other=grey, Staph=green4, Strep=purple
  scale_fill_manual(values=c("orange","green4","grey","blue4","purple")) +
  coord_flip()  # Horizontal bars for better sample name readability


# Export taxonomic composition plot to PDF
ggsave(
  filename = "bracken.pdf",
  plot = bracken_taxa_plot,
  device = "pdf",
  width = 4,     # Width in inches
  height = 30,   # Tall height to accommodate many samples
  units = "in"
)
