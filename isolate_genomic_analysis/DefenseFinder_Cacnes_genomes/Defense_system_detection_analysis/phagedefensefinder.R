# =============================================================================
# DefenseFinder Analysis for Cutibacterium acnes Genomes
# =============================================================================
#
# STUDY OVERVIEW:
# This script processes DefenseFinder output to identify and quantify bacterial
# defense systems in Cutibacterium acnes genomes. DefenseFinder is a tool that
# detects defense systems against phages and other mobile genetic elements
# using profile Hidden Markov Models (HMMs).
#
# INPUT DATA:
# - filtered_blast.csv: DefenseFinder output containing defense system annotations
#   with columns for samples, defense system types, subtypes, and proteins
#
# OUTPUT:
# - Defense system inventory per C. acnes genome
# - Gene counts organized by defense system type and subtype
# - Summary statistics for downstream comparative analysis
#
# METHODOLOGY:
# DefenseFinder identifies defense systems by searching for characteristic
# protein profiles. This analysis processes those results to create a
# comprehensive catalog of defense systems across C. acnes genomes.
#
# =============================================================================

# Load required libraries
library(tidyverse)  # Data manipulation and visualization
library(readr)      # Enhanced data reading functions

# =============================================================================
# DATA IMPORT AND VALIDATION
# =============================================================================

cat("Loading DefenseFinder results for C. acnes genome analysis...\n")

# Define input/output paths for reproducibility
input_dir <- "DefenseFinder_Cacnes_genomes/Defense_system_detection_analysis"
input_file <- file.path(input_dir, "filtered_blast.csv")
output_dir <- "results"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Load DefenseFinder results with error handling
if (!file.exists(input_file)) {
  stop("Error: Input file not found at ", input_file, 
       "\nPlease ensure DefenseFinder output is available.")
}

cat("Reading DefenseFinder results from:", input_file, "\n")
filtered_blast <- read_csv(input_file, show_col_types = FALSE)

# Remove index column if present (common in CSV exports from R/Python)
if (ncol(filtered_blast) > 0 && names(filtered_blast)[1] %in% c("X1", "...1", "")) {
  filtered_blast <- filtered_blast[, -1]
  cat("Removed index column from imported data\n")
}

# Validate required columns are present
required_columns <- c("samples", "type", "subtype", "protein")
missing_columns <- setdiff(required_columns, names(filtered_blast))

if (length(missing_columns) > 0) {
  stop("Error: Missing required columns: ", paste(missing_columns, collapse = ", "),
       "\nExpected columns: ", paste(required_columns, collapse = ", "))
}

cat("Successfully loaded", nrow(filtered_blast), "defense system annotations\n")
cat("Columns available:", paste(names(filtered_blast), collapse = ", "), "\n")

# =============================================================================
# DATA PROCESSING AND FORMATTING
# =============================================================================

cat("Processing defense system annotations...\n")

# Create clean working dataframe with standardized column names
# Select essential columns for defense system analysis
defense_systems_raw <- filtered_blast %>%
  select(samples, type, subtype, protein) %>%
  # Standardize sample name column for consistency
  rename(SampleName = samples) %>%
  # Remove any rows with missing critical information
  filter(!is.na(SampleName), !is.na(type), !is.na(protein))

cat("Processed data contains:\n")
cat("- Samples analyzed:", length(unique(defense_systems_raw$SampleName)), "\n")
cat("- Defense system types:", length(unique(defense_systems_raw$type)), "\n")
cat("- Total protein annotations:", nrow(defense_systems_raw), "\n")

# Display sample of processed data for verification
cat("\nSample of processed data:\n")
print(head(defense_systems_raw, 10))

# Export processed defense system inventory
inventory_file <- file.path(output_dir, "cacnes_defense_systems_inventory.csv")
write_csv(defense_systems_raw, inventory_file)
cat("Defense system inventory exported to:", inventory_file, "\n")

# =============================================================================
# DEFENSE SYSTEM QUANTIFICATION
# =============================================================================

cat("Quantifying defense systems by type and subtype...\n")

# Count defense system genes per sample, organized by type and subtype
# This creates a comprehensive catalog of defense system abundance
gene_counts <- defense_systems_raw %>%
  group_by(SampleName, type, subtype) %>%
  tally(name = "gene_count") %>%
  ungroup() %>%
  # Sort for consistent output
  arrange(SampleName, type, subtype)

cat("Generated gene counts for", nrow(gene_counts), "defense system categories\n")

# Create summary of defense system types found
defense_type_summary <- gene_counts %>%
  group_by(type, subtype) %>%
  summarize(
    genomes_with_system = n(),
    total_genes = sum(gene_count),
    mean_genes_per_genome = round(mean(gene_count), 2),
    .groups = "drop"
  ) %>%
  arrange(desc(total_genes))

cat("\nDefense system types identified:\n")
print(defense_type_summary)

# Export quantified gene counts
gene_count_file <- file.path(output_dir, "cacnes_defense_system_gene_counts.csv")
write_csv(gene_counts, gene_count_file)
cat("Gene counts exported to:", gene_count_file, "\n")

# Export defense system summary
summary_file <- file.path(output_dir, "cacnes_defense_system_summary.csv")
write_csv(defense_type_summary, summary_file)
cat("Defense system summary exported to:", summary_file, "\n")

# =============================================================================
# GENOME-LEVEL DEFENSE SYSTEM ANALYSIS
# =============================================================================

cat("Analyzing defense system profiles per genome...\n")

# Calculate total defense systems per genome
genome_defense_profile <- gene_counts %>%
  group_by(SampleName) %>%
  summarize(
    total_defense_genes = sum(gene_count),
    defense_system_types = n(),
    unique_types = length(unique(type)),
    unique_subtypes = length(unique(subtype)),
    .groups = "drop"
  ) %>%
  arrange(desc(total_defense_genes))

cat("Defense system statistics per genome:\n")
cat("- Mean defense genes per genome:", round(mean(genome_defense_profile$total_defense_genes), 2), "\n")
cat("- Median defense genes per genome:", median(genome_defense_profile$total_defense_genes), "\n")
cat("- Range:", min(genome_defense_profile$total_defense_genes), "-", 
    max(genome_defense_profile$total_defense_genes), "genes\n")

# Export genome-level profiles
genome_profile_file <- file.path(output_dir, "cacnes_genome_defense_profiles.csv")
write_csv(genome_defense_profile, genome_profile_file)
cat("Genome defense profiles exported to:", genome_profile_file, "\n")


# End of DefenseFinder analysis script