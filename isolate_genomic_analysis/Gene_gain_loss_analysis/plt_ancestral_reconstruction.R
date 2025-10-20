################################################################################
################################################################################
# ENVIRONMENT SETUP
################################################################################
################################################################################

# Load required libraries
library(tidyverse)    # Collection of data science packages (dplyr, tidyr, ggplot2, etc.)
library(reshape2)     # Data reshaping (melt/dcast functions)
library(ggpubr)       # Publication-ready plots and statistical comparisons
library(gridExtra)    # Arrange multiple grid-based plots
library(patchwork)    # Compose multiple ggplot objects
library(readr)        # Fast CSV reading/writing

# Disable scientific notation for clearer numeric output
options(scipen=999)

# Create custom 'not in' operator as opposite of %in%
`%notin%` <- Negate(`%in%`)


################################################################################
################################################################################
# DATA IMPORT
################################################################################
################################################################################

################################################################################
# PASTML OUTPUT - ALL GENES
################################################################################
# PastML: Phylogenetic reconstruction of ancestral states
# Load PastML output containing gain/loss events for all C. acnes genes
# across the phylogenetic tree

cacnes_pastmlout <- read_csv("cacnes_pastmlout.csv", 
                             col_types = cols(...1 = col_skip()))  # Skip first unnamed column
View(cacnes_pastmlout)


################################################################################
# PASTML OUTPUT - PDE GENES (FILTERED)
################################################################################
# Load filtered PastML output specifically for prophage-deficient element (PDE) genes
# These are genes associated with prophage regions that have been lost

cacnes_pde_pastmlout_filt <- read_csv("cacnes_pde_pastmlout-filt.csv")
View(cacnes_pde_pastmlout_filt)


################################################################################
################################################################################
# GENE CATEGORIZATION BY PANGENOME CLASSIFICATION
################################################################################
################################################################################

################################################################################
# LOAD GENE PRESENCE/ABSENCE DATA
################################################################################
# Binary matrix of gene presence (1) or absence (0) across C. acnes lineages

cacnes_gene_presence_absence_lineage <- read_csv("cacnes_gene_presence_absence_lineage.csv")
View(cacnes_gene_presence_absence_lineage)


################################################################################
# CATEGORIZE GENES INTO PANGENOME CLASSES
################################################################################
# Classify genes based on prevalence across lineages:
#   - Core genes: present in ≥95% of lineages (not explicitly labeled here)
#   - Shell genes: present in 15-95% of lineages
#   - Cloud genes: present in <15% of lineages (rare/accessory genes)

gene_categories = cacnes_gene_presence_absence_lineage %>%
  group_by(Lineage) %>%
  # Convert from wide to long format (each gene becomes a row)
  reshape2::melt() %>%
  ungroup() %>%
  mutate(gene = variable) %>%  # Rename 'variable' column to 'gene'
  group_by(gene) %>%
  # Keep only lineages where the gene is present (value = 1)
  filter(value == 1) %>%
  summarise(
    gene_cnt = sum(value, na.rm = F),  # Count number of lineages with gene
    # Calculate gene carriage: proportion of lineages carrying the gene
    gene_carriage = sum(value, na.rm = F)/(nrow(cacnes_gene_presence_absence_lineage))
  ) %>%
  ungroup() %>%
  # Classify genes into pangenome categories based on carriage frequency
  mutate(pangenome = ifelse(gene_carriage < 0.15, "cloud_gene",
                            ifelse(gene_carriage < 0.95 & gene_carriage >= 0.15, "shell_gene", NA))) 
  # Note: Core genes (≥95%) are assigned NA and will be excluded from analysis


################################################################################
################################################################################
# COMBINE GAIN/LOSS DATA WITH PANGENOME CATEGORIES
################################################################################
################################################################################

# Merge pangenome classification with PastML gain/loss events for:
#   1. Shell and cloud genes (from general pangenome)
#   2. PDE genes (prophage-deficient elements)

gainloss = rbind(
  # Process shell/cloud genes
  left_join(gene_categories, cacnes_pastmlout) %>%
    # Calculate total evolutionary events (gains + losses)
    mutate(event_cnt = gains + losses) %>%
    select(gene, gains, losses, count, event_cnt, pangenome),
  
  # Add PDE genes as a separate category
  cacnes_pde_pastmlout_filt %>%
    mutate(event_cnt = gains + losses,
           pangenome = "pdes")  # Label all PDE genes with "pdes" category
)

View(gainloss)


################################################################################
################################################################################
# VISUALIZATION: GAIN/LOSS EVENTS BY GENE CATEGORY
################################################################################
################################################################################

################################################################################
# DEFINE GENE CATEGORY ORDER FOR PLOTTING
################################################################################
# Custom ordering: shell → cloud → PDEs
gene_order = c("shell_gene", "cloud_gene", "pdes")


################################################################################
# CREATE COMPARATIVE BOXPLOT WITH STATISTICAL TESTS
################################################################################
# Compare number of gain/loss events across different gene categories
# Hypothesis: PDE genes may show different evolutionary dynamics than
# typical accessory (cloud) or moderately conserved (shell) genes

ggplot(gainloss, aes(x = factor(pangenome, levels = gene_order), y = event_cnt)) + 
  # Individual gene points with transparency and horizontal jitter
  geom_jitter(alpha = 0.3, width = 0.05) +
  # Pairwise statistical comparisons (Wilcoxon rank-sum test by default)
  stat_compare_means(comparisons = list(c("pdes", "shell_gene"),
                                        c("pdes", "cloud_gene"),
                                        c("shell_gene", "cloud_gene"))) +
  theme_pubr() +
  # Add blue crossbar showing mean event count for each category
  stat_summary(fun.y = mean, geom = "crossbar", width = 0.4, color = "blue") +
  labs(title = "C.acnes gain/loss") +
  theme(axis.title.x = element_blank())  # Remove x-axis title for cleaner look
