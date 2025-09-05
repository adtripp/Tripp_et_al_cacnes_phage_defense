# =============================================================================
# Prophage and Mobile Genetic Element Analysis in Cutibacterium acnes
# =============================================================================
#
# STUDY OVERVIEW:
# This script analyzes the presence and distribution of mobile genetic elements
# (prophages, plasmids) across Cutibacterium acnes isolates from multiple datasets.
# The analysis examines carriage patterns by phylogenetic groups and individual
# subjects to understand the evolutionary and epidemiological dynamics of
# mobile genetic elements in C. acnes populations.
#
# INPUT DATA:
# - summaryDT_sample_coverage.csv: Coverage and presence/absence data for mobile elements
# - Cacnes_metatree_GTR_isonames_metadata.csv: Phylogenetic classification data
# - cacnes_isolates_metadata_updated.csv: Sample metadata including subject information
#
# MOBILE GENETIC ELEMENTS ANALYZED:
# - pCRESS plasmids: pCRESS DNA plasmids
# - dsDNA pseudolysogens: Double-stranded DNA prophages  
# - Linear plasmids: Homologous to CP003294 reference
#
# OUTPUT:
# - Heatmap showing mobile element carriage by subphylogroup
# - Subject-level analysis of prophage distribution patterns
# - Family-based comparative visualizations
#
# METHODOLOGY:
# Presence/absence is determined by coverage depth ≥1X and breadth ≥50%.
# Analysis includes phylogenetic context and family-level comparisons.
#
# =============================================================================

# Load required libraries
library(readr)        # Enhanced data reading functions
library(tidyverse)    # Data manipulation and visualization
library(data.table)   # High-performance data operations  
library(ggExtra)      # Additional ggplot2 functionality
library(cowplot)      # Publication-quality plot layouts
library(ggbreak)      # Axis breaks for large ranges
library(patchwork)    # Advanced plot composition
library(zoo)          # Time series and data manipulation
library(reshape2)     # Data reshaping functions

# =============================================================================
# CONFIGURATION AND FILE PATHS
# =============================================================================

# Define input files and output directory
input_dir <- "."
output_dir <- "results"
plot_dir <- "plots"

# Create output directories if they don't exist
for (dir in c(output_dir, plot_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Created directory:", dir, "\n")
  }
}

# Input file paths
coverage_file <- file.path(input_dir, "summaryDT_sample_coverage.csv")
phylo_file <- file.path(input_dir, "Cacnes_megatree_GTR_isonames_metadata.csv")
metadata_file <- file.path(input_dir, "cacnes_isolates_metadata_updated.csv")

# =============================================================================
# DATA IMPORT AND VALIDATION
# =============================================================================

cat("Loading mobile genetic element coverage data...\n")

# Import coverage data for mobile genetic elements
if (!file.exists(coverage_file)) {
  stop("Coverage file not found: ", coverage_file)
}
summaryDT_sample_coverage <- read_csv(coverage_file, show_col_types = FALSE)
cat("Loaded coverage data for", nrow(summaryDT_sample_coverage), "samples\n")

cat("Loading phylogenetic classification data...\n")

# Import phylogroup, subphylogroup, and lineage information
if (!file.exists(phylo_file)) {
  stop("Phylogenetic file not found: ", phylo_file)
}
cacnes_phylo_data <- read_csv(phylo_file, show_col_types = FALSE)
cat("Loaded phylogenetic data for", nrow(cacnes_phylo_data), "isolates\n")

cat("Loading sample metadata...\n")

# Import comprehensive metadata
if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}
sample_metadata <- read_csv(metadata_file, show_col_types = FALSE)
cat("Loaded metadata for", nrow(sample_metadata), "samples\n")

# =============================================================================
# MOBILE GENETIC ELEMENT PRESENCE/ABSENCE DETERMINATION
# =============================================================================

cat("Processing mobile genetic element presence/absence...\n")

# Define presence/absence criteria based on coverage metrics
# Criteria: mean_depth ≥ 1X AND breadth ≥ 50% coverage
DEPTH_THRESHOLD <- 1
BREADTH_THRESHOLD <- 0.5

# Process coverage data and determine presence/absence
mobile_element_summary <- summaryDT_sample_coverage %>%
  # Join with phylogenetic data (remove duplicate columns)
  left_join(cacnes_phylo_data %>% select(-2), by = "SampleName") %>%
  # Filter for C. acnes samples with phylogenetic classification
  filter(!is.na(Phylogroup), Species == "C_acnes") %>%
  
  # Apply presence/absence criteria
  mutate(Present = ifelse(mean_depth >= DEPTH_THRESHOLD & 
                            breadth >= BREADTH_THRESHOLD, 1, 0)) %>%
  
  # Select essential columns for analysis
  select(SampleName, Species, Dataset, Reference_short, 
         Phylogroup, Subphylogroup, Lineage, Present) %>%
  
  # Reshape data to wide format for mobile element analysis
  group_by(SampleName) %>%
  reshape2::dcast(SampleName + Dataset + Phylogroup + Subphylogroup + Lineage ~ 
                    Reference_short, value.var = "Present", fun.aggregate = sum) %>%
  
  # Replace missing values with 0 (absence)
  replace(is.na(.), 0) %>%
  
  # Handle unclassified phylogroups
  mutate(Phylogroup = ifelse(Phylogroup == 0, "Unclassified", Phylogroup))

cat("Processed mobile genetic element data for", nrow(mobile_element_summary), "C. acnes isolates\n")

# =============================================================================
# MOBILE GENETIC ELEMENT CLASSIFICATION
# =============================================================================

cat("Classifying mobile genetic elements...\n")

# Create standardized mobile genetic element categories
mobile_elements_classified <- mobile_element_summary %>%
  mutate(

    # Single-stranded DNA prophages
    pCRESS_plasmid = ifelse(`V: ssDNA_lysogen` == 1, 1, 0),
    
    # Double-stranded DNA prophages
    dsDNA_prophage = ifelse(`V: dsDNA_pseudolysogen` == 1, 1, 0),
    
    # Linear plasmids (CP003294 homologs)
    linear_plasmid = ifelse(`P: CP003294_homologous_region` == 1, 1, 0)
  ) %>%
  
  # Remove raw reference columns and bacterial elements
  select(-starts_with("B: "),           # Remove bacterial elements
         -starts_with("P: "),           # Remove plasmid raw data
         -c(`V: Cryptic_prophage_region`, `V: ssDNA_lysogen`, `V: dsDNA_pseudolysogen`))

cat("Classified mobile genetic elements into standardized categories\n")

# =============================================================================
# METADATA INTEGRATION AND FINAL DATASET PREPARATION
# =============================================================================

cat("Integrating sample metadata...\n")

# Combine mobile element data with comprehensive metadata
final_dataset <- sample_metadata %>%
  left_join(mobile_elements_classified, by = "SampleName") %>%
  
  # Filter for samples with complete information
  filter(!is.na(subject), !is.na(Subphylogroup)) %>%
  
  # Create family groupings based on dataset
  mutate(Family = ifelse(Dataset == "Baker2023", 
                         substr(subject, 1, 1),  # Extract family ID from Baker2023
                         "Conwill2022")) %>%
  
  # Remove control samples (0CC prefix)
  filter(!grepl("0CC", SampleName))

cat("Final dataset contains", nrow(final_dataset), "samples from", 
    length(unique(final_dataset$subject)), "subjects\n")

# =============================================================================
# SUBPHYLOGROUP ANALYSIS AND VISUALIZATION
# =============================================================================

cat("Analyzing mobile element carriage by subphylogroup...\n")

# Define subphylogroup order for consistent visualization
subphylogroup_order <- c("A.1", "A.2", "A.3", "B.1", "C.1", "F.1", "F.3", 
                         "E.1", "D.1", "H.2", "H.1", "H.3", "K.1.1", "K.1.2", "K.2", "K.3")

# Create supplementary data for complete visualization
# This ensures all subphylogroups show up even if they have 0% carriage
supplementary_data <- data.frame(
  Subphylogroup = c("H.3", "K.1.1", "H.3", "K.1.1", "H.3", "K.1.1", 
                    "A.2", "A.2", "B.1", "F.1", "H.2", "K.2", "K.2", "K.3"),
  n_subphylo = c(20, 24, 20, 24, 20, 24, 11, 11, 72, 10, 298, 172, 172, 31),
  value = rep(1, 14),
  variable = c("pCRESS_plasmid", "pCRESS_plasmid", "dsDNA_prophage", "dsDNA_prophage", "linear_plasmid", "linear_plasmid",
               "pCRESS_plasmid", "dsDNA_prophage", "dsDNA_prophage", "dsDNA_prophage", "pCRESS_plasmid", "pCRESS_plasmid", 
               "dsDNA_prophage", "dsDNA_prophage"),
  n = rep(0, 14),
  per_carriage = rep(0, 14)
)

# Calculate carriage frequencies by subphylogroup
subphylogroup_carriage <- final_dataset %>%
  group_by(Subphylogroup) %>%
  mutate(n_subphylo = n()) %>%
  select(n_subphylo, pCRESS_plasmid, dsDNA_prophage, linear_plasmid) %>%
  reshape2::melt(id.vars = c("Subphylogroup", "n_subphylo")) %>%
  group_by(Subphylogroup, n_subphylo, value, variable) %>%
  tally() %>%
  mutate(per_carriage = n / n_subphylo)

# Combine real data with supplementary data for complete visualization
combined_carriage_data <- rbind(subphylogroup_carriage, supplementary_data) %>%
  filter(value == 1)  # Only show presence data

# Create subphylogroup heatmap
subphylogroup_heatmap <- ggplot(combined_carriage_data, 
                                aes(x = variable, 
                                    y = factor(Subphylogroup, levels = rev(subphylogroup_order)), 
                                    fill = per_carriage)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient(low = "white", high = "midnightblue", 
                      name = "Carriage\nFrequency",
                      labels = scales::percent) +
  labs(title = "Mobile Genetic Element Carriage by C. acnes Subphylogroup",
       x = "Mobile Genetic Element Type",
       y = "Subphylogroup") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

# Save subphylogroup heatmap
ggsave(file.path(plot_dir, "subphylogroup_mobile_element_heatmap.png"), 
       subphylogroup_heatmap, width = 8, height = 6, dpi = 300)

cat("Generated subphylogroup carriage heatmap\n")

# =============================================================================
# SUBJECT-LEVEL PROPHAGE CARRIAGE ANALYSIS
# =============================================================================

cat("Analyzing prophage carriage patterns by subject...\n")

# Calculate subject-level prophage carriage by phylogroup
subject_prophage_data <- final_dataset %>%
  select(subject, Phylogroup, Family, dsDNA_prophage) %>%
  group_by(subject, Family) %>%
  mutate(n_subject_total = n()) %>%
  group_by(subject, Phylogroup, Family) %>%
  mutate(n_phylogroup = n()) %>%
  group_by(subject, Phylogroup, Family, dsDNA_prophage) %>%
  mutate(n_prophage = n()) %>%
  unique() %>% 
  ungroup() %>%
  
  # Calculate prophage counts
  mutate(prophage_count = ifelse(dsDNA_prophage == 0, 0, n_prophage)) %>%
  group_by(subject, Phylogroup, Family, n_subject_total, n_phylogroup) %>%
  summarise(total_prophage = sum(prophage_count), .groups = "drop") %>%
  unique()

# Create main heatmap showing prophage carriage by subject and phylogroup
subject_prophage_heatmap <- ggplot(subject_prophage_data, 
                                   aes(y = subject, x = Phylogroup, 
                                       fill = total_prophage / n_phylogroup)) +
  geom_tile(color = ifelse(subject_prophage_data$total_prophage == 0, "black", "red"),
            size = 0.3) +
  scale_fill_gradient(low = "white", high = "midnightblue",
                      name = "Prophage\nFrequency") +
  labs(title = "Prophage Carriage by Subject and Phylogroup",
       x = "C. acnes Phylogroup",
       y = "Subject ID") +
  facet_grid(Family ~ ., scales = "free_y", space = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0))

# Create supplementary barplot showing sample sizes per subject
subject_sample_sizes <- final_dataset %>%
  select(subject, Family) %>%
  group_by(subject, Family) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  unique()

subject_barplot <- ggplot(subject_sample_sizes, 
                          aes(x = n_samples, y = subject)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  labs(title = "Sample Size",
       x = "Number of Isolates",
       y = "") +
  facet_grid(Family ~ ., scales = "free_y", space = "free_y") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank())

# Combine plots using cowplot
combined_subject_plot <- plot_grid(subject_prophage_heatmap, subject_barplot, 
                                   nrow = 1, rel_widths = c(4, 1),
                                   align = "h")

# Save combined subject analysis plot
ggsave(file.path(plot_dir, "subject_prophage_carriage_analysis.png"), 
       combined_subject_plot, width = 12, height = 8, dpi = 300)

cat("Generated subject-level prophage carriage analysis\n")

# =============================================================================
# SUMMARY STATISTICS AND DATA EXPORT
# =============================================================================

cat("Generating summary statistics...\n")

# Calculate overall carriage statistics
carriage_summary <- final_dataset %>%
  summarise(
    total_samples = n(),
    total_subjects = length(unique(subject)),
    ssDNA_carriage = sum(pCRESS_plasmid) / n(),
    dsDNA_carriage = sum(dsDNA_prophage) / n(),
    linear_plasmid_carriage = sum(linear_plasmid) / n()
  )

# Export processed data and summary statistics
write_csv(final_dataset, file.path(output_dir, "cacnes_mobile_elements_final_dataset.csv"))
write_csv(carriage_summary, file.path(output_dir, "mobile_element_carriage_summary.csv"))
write_csv(subphylogroup_carriage, file.path(output_dir, "subphylogroup_carriage_frequencies.csv"))


# End of mobile genetic element analysis script