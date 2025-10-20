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

# Disable scientific notation for clearer numeric output
options(scipen=999)

# Create custom 'not in' operator as opposite of %in%
`%notin%` <- Negate(`%in%`)


################################################################################
################################################################################
# DATA IMPORT AND PREPROCESSING
################################################################################
################################################################################

################################################################################
# PHAGE MAPPING DATA (VIRAL COMPONENT)
################################################################################
# Load samtools coverage statistics for Cutibacterium phage (NC_028967)
# This represents viral reads mapped to the phage reference genome

samtools_coverage_stats_cutiphage_update <- read_csv("samtools_coverage_stats-cutiphage-update.csv")
names(samtools_coverage_stats_cutiphage_update)[1] = "ref"  # Rename first column to "ref"

# Extract phage-specific metrics and clean sample names
phage = samtools_coverage_stats_cutiphage_update %>%
  # Remove file suffix from filename to get clean SampleName
  mutate(SampleName = gsub("_ref_Pacnes_phage_NC028967_cov.tsv","", filename),
         vir_depth = meandepth,    # Mean sequencing depth for viral genome
         vir_cov = coverage,        # Genome coverage percentage for virus
         vir_reads = numreads) %>%  # Number of reads mapped to viral genome
  select(vir_depth, vir_cov, vir_reads, SampleName) 


################################################################################
# BACTERIAL MAPPING DATA (C. ACNES COMPONENT)
################################################################################
# Load samtools coverage statistics for C. acnes reference genome (strain C1)
# This represents bacterial reads mapped to the host genome

samtools_coverage_stats_cacnes <- read_csv("samtools_coverage_stats-cacnes.csv")
names(samtools_coverage_stats_cacnes)[1] = "ref"  # Rename first column to "ref"

# Extract C. acnes-specific metrics and clean sample names
cacnes = samtools_coverage_stats_cacnes %>%
  # Remove file suffix from filename to get clean SampleName
  mutate(SampleName = gsub("_ref_Pacnes_C1_aligned.sorted.bam","", filename),
         bac_depth = meandepth,    # Mean sequencing depth for bacterial genome
         bac_cov = coverage,        # Genome coverage percentage for bacteria
         bac_reads = numreads) %>%  # Number of reads mapped to bacterial genome
  select(bac_depth, bac_cov, bac_reads, SampleName)


################################################################################
# METADATA IMPORT
################################################################################

# Load global skin microbiome study metadata (sample information, study details)
global_skin_metadata <- read_csv("global_skin_metadata.csv")

# Load read length information for each study
global_skin_metadata_readlength <- read_csv("global_skin_metadata-readlength.csv")
View(global_skin_metadata_readlength)


################################################################################
# PHYLOGENETIC GROUPING DATA (PHLAME RESULTS)
################################################################################
# PHyLogenetic Analysis of Metagenomic Estimates (PhAME/PHLAME)
# Load C. acnes subphylogroup frequency data from PHLAME analysis
# Note: Phylogroup-level data commented out, using subphylogroup-level instead

# Cacnes_global_phylogroup_frequencies <- read_csv("Cacnes_global_phylogroup_frequencies.csv")
Cacnes_global_subphylogroup_frequencies <- read_csv("Cacnes_global_subphylogroup_frequencies.csv")


################################################################################
# TAXONOMIC ABUNDANCE DATA (BRACKEN RESULTS)
################################################################################
# Load Bracken abundance estimates specific to C. acnes
bracken_dcast_cacnes <- read_csv("bracken_dcast-cacnes.csv")


################################################################################
# MERGE ALL DATASETS AND CALCULATE VIRUS-TO-MICROBE RATIO (VMR)
################################################################################
# Create comprehensive dataset combining phage, bacterial, phylogenetic,
# metadata, and abundance information

vmr = left_join(
        left_join(
          left_join(
            left_join(
              left_join(phage, cacnes),  # Join phage and bacterial mapping data
              Cacnes_global_subphylogroup_frequencies),  # Add subphylogroup frequencies
            global_skin_metadata),  # Add sample metadata
          bracken_dcast_cacnes)  # Add bracken abundance data
      ) %>% 
      group_by(SampleName) %>%
      mutate(
        # Calculate virus-to-microbe ratio (VMR) as viral depth / bacterial depth
        vmr = vir_depth/bac_depth,
        # Set VMR to 0 if viral coverage is less than 10% (quality filter)
        vmr = ifelse(vir_cov >= 10, vmr, 0)
      ) %>%
      # Remove unnecessary metadata columns to simplify dataset
      select(-c(Instrument, LibraryLayout, LibrarySelection, Site, `Public data?`, 
                `Run Number`, AvgSpotLen, City, Country, 
                `United States, Europe or Canada`, Ethnicity, Age, BioProject))

# Add read length information
vmr = left_join(vmr, global_skin_metadata_readlength)

View(vmr)


################################################################################
################################################################################
# VISUALIZATION: SUBPHYLOGROUP PREVALENCE AND RELATIVE ABUNDANCE
################################################################################
################################################################################

################################################################################
# DEFINE SUBPHYLOGROUP ORDERING
################################################################################
# Custom order for subphylogroups (likely ordered by prevalence or phylogeny)
# subphylo_order = c("A.1","A.2","A.3","B.1","C.1","F.1","F.2","F.3","E.1","D.1","H.1","H.2","H.3","K.1.1","K.1.2","K.2","K.3")
subphylo_order = c("A.3","A.1","D.1","F.3","H.2","K.2","E.1","K.1.1","K.1.2",
                   "F.2","A.2","K.3","C.1","F.1","H.1","B.1","H.3")


################################################################################
# PLOT 1: PREVALENCE BAR CHART (g)
################################################################################
# Calculate and visualize prevalence of each subphylogroup across samples
# Prevalence = proportion of samples where subphylogroup is present (≥1% abundance)

g <- ggplot(
  vmr %>% 
    # Filter for high-quality samples: ≥1M total reads and ≥5x bacterial depth
    filter(tot_reads >= 1000000 & bac_depth >= 5) %>% 
    # Remove non-phylogenetic columns
    select(-c(vir_depth, vir_cov, vir_reads, bac_depth, bac_cov, bac_reads, 
              tot_reads, bracken_cacnes, ApproxAge, Region, Sex, StudyName, 
              `Acne?`, BioProject, NumberIndividuals, ReadLength, `L.1`)) %>% 
    group_by(SampleName, vmr) %>% 
    # Convert from wide to long format (each subphylogroup becomes a row)
    reshape2::melt(id.vars = c("SampleName", "vmr")) %>% 
    select(SampleName, variable, value) %>% 
    # Keep only samples where subphylogroup is ≥1% abundant
    filter(value >= 0.01) %>% 
    group_by(variable) %>% 
    # Calculate prevalence: number of samples with subphylogroup / 471 total samples
    summarise(prevalence = n()/471),
  aes(y = factor(variable, levels = rev(subphylo_order)), 
      x = prevalence, 
      # Color code by PDE status (prophage-deficient vs. non-PDE lineages)
      fill = ifelse(variable %in% c("D.1","H.1","H.2","H.3","K.1.1","K.1.2","K.2","K.3"), 
                    "PDE", "noPDE"))
) +
  geom_col(alpha = 0.8) +
  labs(y = "Subphylogroup", x = "prevalence") +
  theme_pubr() + 
  xlim(0, 0.75) +
  # PDE = olivedrab3 (green), non-PDE = grey50
  scale_fill_manual(values = c("grey50","olivedrab3")) +
  theme(legend.position = "none") 


################################################################################
# PLOT 2: RELATIVE ABUNDANCE SCATTER PLOT (q)
################################################################################
# Visualize distribution of relative abundance values for each subphylogroup

q <- ggplot(
  vmr %>%
    # Filter for high-quality samples
    filter(tot_reads >= 1000000 & bac_depth >= 5) %>%
    # Remove non-phylogenetic columns
    select(-c(vir_depth, vir_cov, vir_reads, bac_depth, bac_cov, bac_reads, 
              tot_reads, bracken_cacnes, ApproxAge, Region, Sex, StudyName, 
              `Acne?`, BioProject, NumberIndividuals, ReadLength, `L.1`)) %>%
    group_by(SampleName, vmr) %>%
    # Convert to long format
    reshape2::melt(id.vars = c("SampleName", "vmr")),
  aes(y = factor(variable, levels = rev(subphylo_order)), 
      x = value,
      # Color code by PDE status
      color = ifelse(variable %in% c("D.1","H.1","H.2","H.3","K.1.1","K.1.2","K.2","K.3"), 
                     "PDE", "noPDE"))
) +
  geom_point(alpha = 0.3, size = 3) +  # Individual sample points
  # Add black crossbar showing mean abundance for each subphylogroup
  stat_summary(fun.y = mean, width = 0.75, color = "black", size = 0.3, 
               geom = "crossbar", alpha = 0.5) +
  labs(y = "Subphylogroup", x = "phlame relative abundance") +
  theme_pubr() +
  scale_color_manual(values = c("grey50","olivedrab3")) +
  theme(legend.position = "none")


# Display both plots side by side
grid.arrange(g, q, ncol = 2, widths = c(1/3, 2/3))


################################################################################
# EXPORT COMBINED PLOT
################################################################################

# Save combined prevalence and abundance plot to PDF
ggsave(
  filename = "phalme_subphylo_rel_abund.pdf",
  plot = grid.arrange(g, q, ncol = 2, widths = c(1/3, 2/3)),
  device = "pdf",
  width = 11,    # Width in inches
  height = 6,    # Height in inches
  units = "in"
)


################################################################################
################################################################################
# STATISTICAL COMPARISON: PDE vs NON-PDE SUBPHYLOGROUPS
################################################################################
################################################################################

################################################################################
# PLOT 3: PREVALENCE COMPARISON (h)
################################################################################
# Compare prevalence between PDE and non-PDE subphylogroups using dotplot

h <- ggplot(
  vmr %>% 
    # Filter for high-quality samples
    filter(tot_reads >= 1000000 & bac_depth >= 5) %>% 
    # Remove non-phylogenetic columns
    select(-c(vir_depth, vir_cov, vir_reads, bac_depth, bac_cov, bac_reads, 
              tot_reads, bracken_cacnes, ApproxAge, Region, Sex, StudyName, 
              `Acne?`, BioProject, NumberIndividuals, ReadLength, `L.1`)) %>% 
    group_by(SampleName, vmr) %>% 
    reshape2::melt(id.vars = c("SampleName", "vmr")) %>% 
    select(SampleName, variable, value) %>% 
    # Keep only samples where subphylogroup is ≥1% abundant
    filter(value >= 0.01) %>% 
    group_by(variable) %>% 
    # Calculate prevalence for each subphylogroup
    summarise(prevalence = n()/471) %>%
    # Categorize as PDE or non-PDE
    mutate(myfill = ifelse(variable %in% c("D.1","H.1","H.2","H.3","K.1.1","K.1.2","K.2","K.3"), 
                           "PDE", "noPDE")),
  aes(x = factor(myfill), 
      y = prevalence, 
      fill = myfill, 
      color = myfill)
) +
  # Dotplot showing individual subphylogroup prevalence values
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions = "all", alpha = 0.5) +
  # Add statistical comparison (Wilcoxon test p-value)
  stat_compare_means(label.y = 0.8) +
  scale_fill_manual(values = c("grey50","olivedrab3")) +   
  scale_color_manual(values = c("grey50","olivedrab3")) +
  # Add black crossbar showing mean prevalence for each group
  stat_summary(fun.y = mean, geom = "crossbar", width = 0.3, color = "black") +
  theme_pubr() 


################################################################################
# PLOT 4: MEAN RELATIVE ABUNDANCE COMPARISON (r)
################################################################################
# Compare mean relative abundance between PDE and non-PDE subphylogroups

r <- ggplot(
  vmr %>% 
    # Filter for high-quality samples
    filter(tot_reads >= 1000000 & bac_depth >= 5) %>% 
    # Remove non-phylogenetic columns
    select(-c(vir_depth, vir_cov, vir_reads, bac_depth, bac_cov, bac_reads, 
              tot_reads, bracken_cacnes, ApproxAge, Region, Sex, StudyName, 
              `Acne?`, BioProject, NumberIndividuals, ReadLength, `L.1`)) %>% 
    group_by(SampleName, vmr) %>% 
    reshape2::melt(id.vars = c("SampleName", "vmr")) %>% 
    group_by(variable) %>% 
    # Calculate mean relative abundance for each subphylogroup across all samples
    reframe(freq = mean(value, na.rm = T),
            myfill = ifelse(variable %in% c("D.1","H.1","H.2","H.3","K.1.1","K.1.2","K.2","K.3"), 
                            "PDE", "noPDE")) %>% 
    unique(),
  aes(x = factor(myfill), 
      y = freq, 
      color = myfill, 
      fill = myfill)
) +
  # Dotplot showing mean abundance for each subphylogroup
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions = "all", alpha = 0.5) +
  # Add statistical comparison (Wilcoxon test p-value)
  stat_compare_means(label.y = 0.25) +
  scale_fill_manual(values = c("grey50","olivedrab3")) +   
  scale_color_manual(values = c("grey50","olivedrab3")) +
  # Add black crossbar showing overall mean for each group
  stat_summary(fun.y = mean, geom = "crossbar", width = 0.3, color = "black") +
  theme_pubr()


# Display both comparison plots side by side
grid.arrange(h, r, ncol = 2)


################################################################################
# EXPORT STATISTICAL COMPARISON PLOTS
################################################################################

# Save PDE vs non-PDE comparison plots to PDF
ggsave(
  filename = "phalme_subphylo_rel_abund_pvalues.pdf",
  plot = grid.arrange(h, r, ncol = 2),
  device = "pdf",
  width = 8,     # Width in inches
  height = 4,    # Height in inches
  units = "in"
)
