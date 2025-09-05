# =============================================================================
# CRISPR Spacer Jaccard Distance Analysis
# 
# This script analyzes the output of SpacerExtractor, a tool that extracts 
# CRISPR-Cas spacers from raw short reads. It calculates Jaccard distances
# between samples based on shared spacer sequences and compares distances
# across different phylogenetic groupings.
#
# Input: SpacerExtractor-Filt_Output.csv (filtered spacer extraction results)
# Output: Jaccard distance boxplot comparing within-lineage vs between-lineage distances
# =============================================================================

# Libraries
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(vegan)
library(magrittr)
library(ade4)

# Import datasets
spacer_data <- read_csv("SpacerExtractor-Filt_Output.csv")


# Process data for Jaccard distance calculation
df_for_jaccard <- spacer_data %>%
  filter(Array == "Ac_14669") %>%
  select(SampleName, spacer_seq) %>%
  mutate(random = 1) %>%
  dcast(SampleName ~ spacer_seq, value.var = "random") 

# Replace NA values with 0
df_for_jaccard[is.na(df_for_jaccard)] <- 0

# Calculate Jaccard distances
distance_matrix <- dist.binary(df_for_jaccard[, -1], method = 1)
distance_df <- as.data.frame(as.matrix(distance_matrix)) 

# Format distance matrix
colnames(distance_df) <- df_for_jaccard$SampleName
distance_df <- cbind(df_for_jaccard$SampleName, distance_df)
colnames(distance_df)[1] <- "SampleName"

# Prepare final dataset for plotting
df_jaccard <- left_join(
  left_join(distance_df %>%
              pivot_longer(!SampleName) %>%
              filter(SampleName != name),
            iso_metadata) %>%
    set_colnames(c("SampleA", "SampleName", "Jaccard", "PhyA", "LinA")),
  iso_metadata) %>%
  set_colnames(c("SampleA", "SampleB", "Jaccard", "PhyA", "LinA", "PhyB", "LinB")) %>%
  mutate(shared = ifelse(PhyA == PhyB, "2-subphylo", "1-phylo"),
         shared = ifelse(LinA == LinB, "3-lineage", shared))

# Create boxplot
jaccard_plot <- ggplot(df_jaccard, aes(y = Jaccard, x = shared)) +
  geom_boxplot() +
  theme_classic() +
  stat_compare_means(
    aes(label = ..p.signif..), 
    comparisons = list(
      c("1-phylo", "3-lineage"), 
      c("3-lineage", "2-subphylo"), 
      c("1-phylo", "2-subphylo")
    )
  )

# Display plot
print(jaccard_plot)

# Optional: Save plot
# ggsave("jaccard_distance_plot.png", jaccard_plot, width = 8, height = 6, dpi = 300)

