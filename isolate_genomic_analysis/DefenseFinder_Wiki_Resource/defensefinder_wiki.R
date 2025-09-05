# =============================================================================
# Analysis of Defense Systems in Skin-Associated Bacterial Genera
# =============================================================================
#
# STUDY OVERVIEW:
# This analysis examines the distribution of bacterial defense systems across
# four key skin-associated bacterial genera: Corynebacterium, Cutibacterium, 
# Staphylococcus, and Streptococcus. These genera represent important components 
# of the human skin microbiome and are relevant to skin health and disease.
#
# DATA SOURCE:
# Defense system annotations from the DefenseFinder database, derived from 
# RefSeq genomes spanning two major bacterial phyla: Actinomycetota and Bacillota.
#
# METHODOLOGY:
# - Quantifies defense system abundance per genome assembly
# - Performs statistical comparisons between genera using Wilcoxon rank-sum tests
# - Provides detailed species-level analysis for Cutibacterium
# - Generates visualizations
#
# STATISTICAL APPROACH:
# Non-parametric tests (Mann-Whitney U) account for non-normal distribution
# of defense system counts. All pairwise comparisons between genera are performed
# with appropriate visualization of statistical significance.
#
# =============================================================================

# Load required libraries for data processing, statistical analysis, and visualization
library(readr)      # Data import functions
library(dplyr)      # Data manipulation and processing
library(ggplot2)    # Core plotting functionality
library(ggpubr)     # Publication-ready plots and statistical testing
library(gridExtra)  # Multi-panel figure arrangement

# =============================================================================
# DATA IMPORT AND INITIAL PROCESSING
# =============================================================================

# Define file paths for reproducibility and easy modification
# Modify these paths to match your directory structure
data_dir <- "."
actinomycetota_file <- file.path(data_dir, "df-RefSeq_phylum_actinomycetota.csv")
bacillota_file <- file.path(data_dir, "df-RefSeq_phylum_bacillota.csv")

# Create reusable function to read and clean defense system data
# Removes redundant annotation fields as described in methods:
# - System coordinates (sys_beg, sys_end, sys_id)
# - Protein identifiers (protein_in_syst, name_of_profiles_in_sys, accession_in_sys)
# - Higher-level taxonomic classifications (Superkingdom, phylum)
read_defense_data <- function(file_path) {
  cat("Reading defense system data from:", basename(file_path), "\n")
  read_csv(file_path, show_col_types = FALSE) %>%
    select(-c(sys_beg, sys_end, sys_id, protein_in_syst, 
              name_of_profiles_in_sys, accession_in_sys, 
              Superkingdom, phylum))
}

# Import defense system data from both major bacterial phyla
cat("Importing defense system annotations from DefenseFinder database...\n")
df_actinomycetota <- read_defense_data(actinomycetota_file)
df_bacillota <- read_defense_data(bacillota_file)

# =============================================================================
# TAXONOMIC FILTERING AND DEFENSE SYSTEM CLASSIFICATION
# =============================================================================

# Define target genera: key skin-associated bacteria relevant to skin microbiome research
target_genera <- c("Corynebacterium", "Cutibacterium", "Staphylococcus", "Streptococcus")
target_pattern <- paste(target_genera, collapse = "|")

cat("Processing defense systems for skin-associated bacterial genera:", 
    paste(target_genera, collapse = ", "), "\n")

# Combine phyla data and process defense system annotations
# This implements the methodology described in the paper:
defense_systems <- rbind(df_actinomycetota, df_bacillota) %>%
  # Filter for skin-associated bacterial genera of interest
  filter(grepl(target_pattern, species)) %>%
  
  # Classify defense systems: present (subtype available) vs absent (NA subtype)
  # This binary classification handles mixed annotations in the database
  mutate(system = ifelse(is.na(subtype), "no sys", "sys")) %>%
  
  # Count defense systems per genome assembly
  # Groups by Assembly to ensure each genome is represented once
  group_by(Assembly, genus, species, system) %>%
  tally() %>%
  
  # Convert "no defense system" classifications to count of 0
  # This standardizes the quantification approach
  mutate(n = ifelse(system == "no sys", 0, n)) %>%
  ungroup() %>%
  
  # Handle assemblies with mixed annotations (both defense systems and no-system regions)
  # Retain maximum count per assembly as described in methods
  group_by(Assembly) %>%
  slice_max(order_by = n, n = 1) %>%
  ungroup()

cat("Processed", nrow(defense_systems), "genome assemblies\n")
cat("Genera represented:", paste(unique(defense_systems$genus), collapse = ", "), "\n")

# =============================================================================
# STATISTICAL SUMMARY AND DATA EXPORT
# =============================================================================

cat("Calculating summary statistics by genus and species...\n")

# Generate comprehensive summary statistics as described in methods
# Includes sample sizes, distribution measures for defense system counts
summary_stats <- defense_systems %>%
  group_by(genus, species) %>%
  summarize(
    n_isolates = n(),                    # Sample size per species
    min_defense_systems = min(n),        # Minimum defense system count
    max_defense_systems = max(n),        # Maximum defense system count  
    median_defense_systems = median(n),  # Median for non-parametric analysis
    mean_defense_systems = round(mean(n), 2),  # Mean for reference
    .groups = "drop"
  )

# Export summary statistics for supplementary materials
output_file <- "defense_system_summary_stats.csv"
write_csv(summary_stats, output_file)
cat("Summary statistics exported to:", output_file, "\n")

# Display summary for immediate review
print(summary_stats)

# =============================================================================
# STATISTICAL TESTING SETUP
# =============================================================================

# Prepare pairwise comparisons for statistical testing
# All possible genus pairs will be compared using Wilcoxon rank-sum tests
# This non-parametric approach accounts for non-normal distribution of counts
genera_levels <- levels(as.factor(defense_systems$genus))
pairwise_comparisons <- combn(genera_levels, 2, simplify = FALSE)

cat("Prepared", length(pairwise_comparisons), "pairwise comparisons between genera\n")

# =============================================================================
# VISUALIZATION 1: GENUS-LEVEL COMPARISON
# =============================================================================

cat("Generating genus-level comparison plot...\n")

# Create publication-ready violin plot comparing defense systems across genera
# Combines probability density (violin) with summary statistics (boxplot)
# Statistical significance displayed directly on plot
plot_genus_comparison <- ggplot(defense_systems, aes(x = genus, y = n)) +
  
  # Violin plot shows probability density distribution of defense system counts
  geom_violin(alpha = 0.6, trim = FALSE) +
  
  # Boxplot overlay highlights median, quartiles (outliers suppressed for clarity)
  geom_boxplot(width = 0.15, alpha = 0.7, outlier.shape = NA) +
  
  # Publication-ready labels and formatting
  labs(
    title = "Defense System Abundance Across Skin-Associated Bacterial Genera",
    x = "Bacterial Genus",
    y = "Number of Defense Systems per Genome"
  ) +
  
  # Clean, publication-appropriate theme
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1),
    legend.position = "right",
    plot.title = element_text(size = 12)
  ) +
  
  # Add statistical comparisons using Wilcoxon rank-sum test
  # Non-parametric test appropriate for defense system count distributions
  stat_compare_means(
    comparisons = pairwise_comparisons,
    method = "wilcox.test",
    label = "p.signif"  # Display significance symbols (*, **, ***)
  )

# =============================================================================
# VISUALIZATION 2: CUTIBACTERIUM SPECIES-LEVEL ANALYSIS
# =============================================================================

cat("Generating detailed Cutibacterium species analysis...\n")

# Filter for Cutibacterium genus for detailed species-level analysis
# This genus is of particular interest in skin microbiome research
cutibacterium_data <- defense_systems %>%
  filter(genus == "Cutibacterium")

# Prepare species-level comparisons within Cutibacterium
cutibacterium_species <- levels(as.factor(cutibacterium_data$species))
cutibacterium_comparisons <- combn(cutibacterium_species, 2, simplify = FALSE)

cat("Analyzing", length(unique(cutibacterium_data$species)), "Cutibacterium species\n")

# Create detailed species comparison plot
# Uses boxplot + dotplot combination suitable for smaller sample sizes
plot_cutibacterium_species <- ggplot(cutibacterium_data, aes(x = species, y = n)) +
  
  # Boxplot shows distribution summary statistics
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
  
  # Dot plot shows individual data points (appropriate for smaller sample sizes)
  geom_dotplot(binaxis = "y", stackdir = "center", alpha = 0.5) +
  
  # Species-specific labels and formatting
  labs(
    title = "Defense System Distribution in Cutibacterium Species",
    x = "Cutibacterium Species",
    y = "Number of Defense Systems per Genome"
  ) +
  
  # Clean formatting for species names
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1),
    legend.position = "none",
    plot.title = element_text(size = 12)
  ) +
  
  # Statistical comparisons between Cutibacterium species
  stat_compare_means(
    comparisons = cutibacterium_comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )

# =============================================================================
# FINAL VISUALIZATION AND RESULTS SUMMARY
# =============================================================================

cat("Creating combined publication figure...\n")

# Arrange plots in publication-ready format
# Left panel: genus comparison, Right panel: Cutibacterium species detail
combined_plot <- grid.arrange(
  plot_genus_comparison, 
  plot_cutibacterium_species, 
  ncol = 2,
  top = "Bacterial Defense Systems in Skin-Associated Microbiota"
)

# =============================================================================
# ANALYSIS COMPLETION SUMMARY
# =============================================================================

cat("\n" %+% "="*60 %+% "\n")
cat("ANALYSIS COMPLETE - SUMMARY REPORT\n")
cat("="*60 %+% "\n")
cat("Total genome assemblies analyzed:", nrow(defense_systems), "\n")
cat("Bacterial genera included:", paste(unique(defense_systems$genus), collapse = ", "), "\n")
cat("Species analyzed:", length(unique(defense_systems$species)), "\n")
cat("Defense systems range:", min(defense_systems$n), "-", max(defense_systems$n), "per genome\n")
cat("Summary statistics saved to:", output_file, "\n")
cat("Figures generated: Genus comparison + Cutibacterium species analysis\n")
cat("Statistical method: Wilcoxon rank-sum test (Mann-Whitney U)\n")
cat("\nAnalysis ready for publication submission.\n")
cat("="*60 %+% "\n")

# End of analysis script