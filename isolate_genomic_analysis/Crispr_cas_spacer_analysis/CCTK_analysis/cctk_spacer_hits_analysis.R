# =============================================================================
# CRISPR Spacer Target Analysis in Cutibacterium acnes
# =============================================================================
#
# STUDY OVERVIEW:
# This script analyzes CRISPR spacer sequences from C. acnes genomes and their
# targeting patterns against viral and plasmid databases. CRISPR-Cas systems
# provide adaptive immunity against phages and other mobile genetic elements
# by storing "memory" sequences (spacers) that match foreign DNA.
#
# ANALYSIS COMPONENTS:
# 1. Spacer target identification across multiple databases
# 2. Contig analysis of spacer hits (assembly contigs vs reference databases)
# 3. IMG/VR and IMG/PR database comparisons (viral and prokaryotic targets)
# 4. Custom database analysis against known C. acnes mobile elements
# 5. Phylogenetic context and spacer frequency analysis
# 6. Clustering analysis of spacer carriage patterns
#
# INPUT DATABASES:
# - IMGVR: IMG/VR viral reference database
# - IMGPR: IMG/PR prokaryotic reference database  
# - Custom: C. acnes-specific prophage and plasmid sequences
#
# METHODOLOGY:
# Spacer sequences are BLASTed against target databases with 90% identity
# threshold. Results are analyzed for targeting patterns, frequency across
# lineages, and phylogenetic distribution.
#
# OUTPUT:
# - Comprehensive spacer targeting profiles
# - Phylogenetic heatmaps of spacer carriage
# - Clustering analysis based on spacer presence/absence
# - Target genome visualizations
#
# =============================================================================

# Load required libraries
library(readr)        # Enhanced data reading
library(tidyverse)    # Data manipulation and visualization
library(readxl)       # Excel file reading
library(reshape2)     # Data reshaping
library(data.table)   # High-performance data operations
library(RColorBrewer) # Color palettes
library(ggExtra)      # Additional ggplot2 functionality
library(magrittr)     # Pipe operators
library(ggtext)       # Enhanced text rendering
library(gridExtra)    # Plot arrangement

# =============================================================================
# CONFIGURATION AND FILE PATHS
# =============================================================================

# Define directory structure for reproducibility
input_dir <- "./in_files"
output_dir <- "./out_files"
plot_dir <- "./plots"

# Create output directories if they don't exist
for (dir in c(output_dir, plot_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Created directory:", dir, "\n")
  }
}

# Define analysis parameters
IDENTITY_THRESHOLD <- 90  # Minimum identity for spacer-target matches
MIN_SPACER_LENGTH <- 20   # Minimum spacer length for analysis

# =============================================================================
# SPACER DATA IMPORT AND INITIAL PROCESSING
# =============================================================================

cat("Loading CRISPR spacer sequence identifiers...\n")

# Import spacer identifiers
spacer_names_file <- file.path(input_dir, "spacer_names-filt.txt")
if (!file.exists(spacer_names_file)) {
  stop("Spacer names file not found: ", spacer_names_file)
}

spacer_names <- read.delim(spacer_names_file, header = FALSE) %>%
  set_colnames("Spacer_ID")

cat("Loaded", nrow(spacer_names), "spacer sequences for analysis\n")

# =============================================================================
# CONTIG TARGET ANALYSIS
# =============================================================================

cat("Analyzing spacer targets in assembled contigs...\n")

# Import and process contig hit data
contig_hits_file <- file.path(input_dir, "spacer_names-filt-90.tsv")
if (!file.exists(contig_hits_file)) {
  stop("Contig hits file not found: ", contig_hits_file)
}

# Process contig targeting data
# Extract contig length information from target names
contig_hits <- read.delim(contig_hits_file) %>%
  # Parse contig length from target contig names (format: name_length_XXXX_cov_XX)
  separate(Target_contig, sep = "length_", into = c("dump", "keep"), remove = FALSE) %>%
  separate(keep, sep = "_cov", into = c("contig_length", "dump1")) %>%
  select(-c(dump, dump1)) %>%
  
  # Summarize hits per spacer
  group_by(Spacer_ID) %>%
  summarise(
    n_contigs = n(),                                           # Number of target contigs
    median_contig_length = median(as.numeric(contig_length)),  # Median target length
    .groups = "drop"
  )

cat("Processed contig targets for", nrow(contig_hits), "spacers\n")
cat("Spacers targeting assembled contigs:", sum(contig_hits$n_contigs > 0), "\n")

# =============================================================================
# IMG DATABASE TARGET ANALYSIS
# =============================================================================

cat("Analyzing spacer targets in IMG databases...\n")

# Define standard BLAST output column names
blast_columns <- c("qcov", "qcovhsp", "Spacer_ID", "Target_contig", "PID", "Length", 
                   "Mismatch", "Gap", "qstart", "qend", "sstart", "send", 
                   "evalue", "bitscore", "qlen", "slen")

# Import IMGVR (viral) database hits
imgvr_file <- file.path(input_dir, "IMGVR_db-spacersfilt.tsv")
imgpr_file <- file.path(input_dir, "IMGPR_db-spacersfilt.tsv")

img_hits <- NULL

# Process IMGVR database hits (viral targets)
if (file.exists(imgvr_file)) {
  imgvr_hits <- read.delim(imgvr_file, header = FALSE) %>%
    set_colnames(blast_columns) %>%
    mutate(IMGdb = "IMGVR")
  img_hits <- imgvr_hits
  cat("Loaded IMGVR hits for", length(unique(imgvr_hits$Spacer_ID)), "spacers\n")
}

# Process IMGPR database hits (prokaryotic targets)
if (file.exists(imgpr_file)) {
  imgpr_hits <- read.delim(imgpr_file, header = FALSE) %>%
    set_colnames(blast_columns) %>%
    mutate(IMGdb = "IMGPR")
  
  if (is.null(img_hits)) {
    img_hits <- imgpr_hits
  } else {
    img_hits <- rbind(img_hits, imgpr_hits)
  }
  cat("Loaded IMGPR hits for", length(unique(imgpr_hits$Spacer_ID)), "spacers\n")
}

# Summarize IMG database targeting
if (!is.null(img_hits)) {
  img_summary <- img_hits %>%
    group_by(Spacer_ID) %>%
    summarise(IMGdb = paste(unique(IMGdb), collapse = ";"), .groups = "drop")
  
  cat("IMG database summary: ", nrow(img_summary), "spacers with targets\n")
} else {
  cat("Warning: No IMG database files found\n")
  img_summary <- data.frame(Spacer_ID = character(), IMGdb = character())
}

# =============================================================================
# COMPREHENSIVE SPACER TARGET SUMMARY
# =============================================================================

cat("Creating comprehensive spacer targeting summary...\n")

# Combine all targeting information
spacer_targeting_summary <- spacer_names %>%
  # Add IMG database targeting information
  left_join(img_summary, by = "Spacer_ID") %>%
  # Add contig targeting information
  left_join(contig_hits, by = "Spacer_ID") %>%
  # Fill missing values
  mutate(
    IMGdb = ifelse(is.na(IMGdb), "None", IMGdb),
    n_contigs = ifelse(is.na(n_contigs), 0, n_contigs),
    median_contig_length = ifelse(is.na(median_contig_length), 0, median_contig_length)
  )

# Export comprehensive targeting summary
summary_output_file <- file.path(output_dir, "spacer_targeting_comprehensive_summary.csv")
write_csv(spacer_targeting_summary, summary_output_file)
cat("Exported comprehensive summary to:", summary_output_file, "\n")

# =============================================================================
# CUSTOM DATABASE ANALYSIS - C. ACNES MOBILE ELEMENTS
# =============================================================================

cat("Analyzing targeting of known C. acnes mobile genetic elements...\n")

# Import custom database results (C. acnes prophages and plasmids)
custom_db_file <- file.path(input_dir, "spacer_names-Cacnes_prophage_plasmid_genomes-90.tsv")
if (file.exists(custom_db_file)) {
  custom_db_hits <- read.delim(custom_db_file)
  cat("Loaded custom database hits:", nrow(custom_db_hits), "total hits\n")
  
  # Identify key mobile genetic elements for visualization
  target_elements <- c("CP003294_homologous_region_length_31925",  # Linear plasmid
                       "NC_028967.1")                              # Known prophage
  
  # Create targeting visualizations for key mobile elements
  for (element in target_elements) {
    if (element %in% custom_db_hits$Target_contig) {
      cat("Creating targeting map for:", element, "\n")
      
      element_hits <- custom_db_hits %>%
        filter(Target_contig == element) %>%
        select(Spacer_ID, Protospacer_start, Protospacer_end) %>%
        reshape2::melt(id.vars = "Spacer_ID")
      
      # Create targeting visualization
      element_plot <- ggplot(element_hits, 
                             aes(x = value, y = 0, group = Spacer_ID)) +
        geom_line(size = 3, alpha = 0.7, color = "darkblue") +
        labs(title = paste("CRISPR Spacer Targeting Positions:", element),
             x = "Genomic Position (bp)",
             y = "") +
        theme_classic() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.title = element_text(hjust = 0.5))
      
      # Save individual targeting plots
      element_name <- gsub("[^A-Za-z0-9]", "_", element)
      ggsave(file.path(plot_dir, paste0("targeting_map_", element_name, ".png")), 
             element_plot, width = 10, height = 3, dpi = 300)
    }
  }
} else {
  cat("Warning: Custom database file not found:", custom_db_file, "\n")
}

# =============================================================================
# PHYLOGENETIC AND SPACER FREQUENCY ANALYSIS
# =============================================================================

cat("Loading phylogenetic context and spacer frequency data...\n")

# Import phylogenetic classification data
phylo_file <- file.path(input_dir, "Cacnes_ALL_summary_tree.csv")
if (file.exists(phylo_file)) {
  cacnes_phylo_summary <- read_csv(phylo_file)
  cat("Loaded phylogenetic data for", nrow(cacnes_phylo_summary), "isolates\n")
} else {
  cat("Warning: Phylogenetic summary file not found:", phylo_file, "\n")
}

# Import CRISPR array location data
array_locations_file <- file.path(input_dir, "Array_locations.bed")
if (file.exists(array_locations_file)) {
  array_locations <- read.delim(array_locations_file, header = TRUE, sep = "\t")
  cat("Loaded CRISPR array locations for", nrow(array_locations), "arrays\n")
} else {
  cat("Warning: Array locations file not found:", array_locations_file, "\n")
}

# Import spacer frequency analysis
spacer_freq_file <- file.path(output_dir, "spacer_lineage_permenance.csv")
if (file.exists(spacer_freq_file)) {
  spacer_frequency <- read_csv(spacer_freq_file)
  names(spacer_frequency)[2] <- "Spacer_ID"
  
  # Integrate with targeting data if available
  manual_annotation_file <- file.path(output_dir, "spacerblast_names-filt-manual.csv")
  if (file.exists(manual_annotation_file)) {
    manual_annotations <- read_csv(manual_annotation_file)
    spacer_frequency <- left_join(spacer_frequency, manual_annotations, by = "Spacer_ID")
  }
  
  cat("Loaded spacer frequency data for", length(unique(spacer_frequency$Spacer_ID)), "spacers\n")
} else {
  cat("Note: Spacer frequency file not found. Skipping frequency analysis.\n")
  spacer_frequency <- NULL
}

# =============================================================================
# SPACER CARRIAGE PATTERN VISUALIZATION AND CLUSTERING
# =============================================================================

if (!is.null(spacer_frequency)) {
  cat("Performing clustering analysis of spacer carriage patterns...\n")
  
  # Create matrix for clustering analysis
  spacer_matrix_data <- spacer_frequency %>%
    select(TreeTipLabel, Spacer_ID, percentage) %>%
    tidyr::pivot_wider(names_from = "Spacer_ID", values_from = "percentage")
  
  # Convert to matrix format
  spacer_matrix <- as.matrix(spacer_matrix_data[, -1])
  spacer_matrix[is.na(spacer_matrix)] <- 0
  
  # Perform hierarchical clustering based on Euclidean distance
  spacer_clustering <- hclust(dist(t(spacer_matrix)))
  spacer_order <- colnames(spacer_matrix)[spacer_clustering$order]
  
  cat("Completed clustering of", ncol(spacer_matrix), "spacers\n")
  
  # Create phylogenetic heatmap with clustered spacer order
  if ("Plot_Label" %in% names(spacer_frequency)) {
    phylo_heatmap <- ggplot(spacer_frequency, 
                            aes(x = Spacer_ID, y = as.factor(order), fill = Plot_Label)) +
      geom_tile(color = "white", size = 0.1) +
      scale_x_discrete(limits = spacer_order) +
      scale_fill_manual(values = c("lightgrey", "cadetblue3", "green4"),
                        name = "Carriage\nStatus") +
      labs(title = "CRISPR Spacer Carriage Patterns Across C. acnes Phylogeny",
           x = "CRISPR Spacer (Clustered by Similarity)",
           y = "Phylogenetic Position") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "top",
            plot.title = element_text(hjust = 0.5))
  } else {
    # Alternative visualization using percentage data
    phylo_heatmap <- ggplot(spacer_frequency,
                            aes(x = Spacer_ID, y = TreeTipLabel, fill = percentage)) +
      geom_tile() +
      scale_x_discrete(limits = spacer_order) +
      scale_fill_gradient(low = "grey90", high = "grey10", name = "Frequency") +
      labs(title = "CRISPR Spacer Frequency Across C. acnes Isolates",
           x = "CRISPR Spacer (Clustered)",
           y = "C. acnes Isolate") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(size = 6))
  }
  
  # Create spacer frequency distribution boxplot
  frequency_boxplot <- ggplot(spacer_frequency, 
                              aes(x = Spacer_ID, y = percentage)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    scale_x_discrete(limits = spacer_order) +
    labs(title = "Distribution of Spacer Carriage Frequencies",
         x = "CRISPR Spacer (Clustered by Similarity)",
         y = "Carriage Frequency (%)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
  
  # Combine plots
  combined_spacer_plot <- gridExtra::grid.arrange(phylo_heatmap, frequency_boxplot, 
                                                  nrow = 2, heights = c(2, 1))
  
  # Save combined visualization
  ggsave(file.path(plot_dir, "spacer_carriage_clustering_analysis.png"), 
         combined_spacer_plot, width = 14, height = 10, dpi = 300)
  
  cat("Generated spacer carriage pattern visualization\n")
}


# End of CRISPR spacer target analysis script