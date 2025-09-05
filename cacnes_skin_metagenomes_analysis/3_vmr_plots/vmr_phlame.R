################################################################################
# Title: C. acnes Phage Defense Analysis
# Author: [Your Name]
# Date: [Date]
# Description: 
#   This script processes sequencing coverage data for C. acnes bacteria and 
#   phages, merges metadata and taxonomic abundance data, calculates virus-to-
#   microbe ratios (VMR), and generates summary plots comparing phage presence,
#   bacterial abundance, and resistance frequencies across studies.
# 
# Phage mapping data: Coverage statistics for C. acnes phages from sequencing reads mapped to phage reference genomes.
# Bacteria mapping data: Coverage statistics for C. acnes bacterial genomes from sequencing reads.
# Metadata: Sample and study-level metadata including study names, geographic regions, sequencing read lengths, and other sample characteristics.
# Phylogroup frequencies: Frequencies of C. acnes subphylogroups or resistant phylogroups (e.g., groups labeled D, H, K), likely derived from phylogenetic or taxonomic profiling.
# Bracken abundance data: Taxonomic abundance estimates from Bracken software for C. acnes in samples.
# The analysis calculates virus-to-microbe ratios (VMR), evaluates the prevalence of phages and bacteria across samples and studies, and investigates resistance frequencies in C. acnes populations.
# 
################################################################################


# ---------------------------- #
# 1. Setup Environment          #
# ---------------------------- #

# Set working directory (modify as needed)
setwd(".")

# Load required libraries
library(tidyverse)    # Data manipulation and visualization
library(reshape2)     # Data reshaping
library(ggpubr)       # Publication ready plots
library(gridExtra)    # Arrange multiple plots
library(patchwork)    # Combine ggplots
library(ggridges)     # Ridgeline plots
library(broom) 

# Global options
options(scipen = 999) # Disable scientific notation for readability


# ---------------------------- #
# 2. Data Import and Wrangling #
# ---------------------------- #

## 2.1 Import Phage Mapping Data
phage_coverage <- read_csv("samtools_coverage_stats-cutiphage-update.csv")

# Rename first column for clarity
names(phage_coverage)[1] <- "ref"

# Extract relevant columns and create sample names
phage_data <- phage_coverage %>%
  mutate(
    SampleName = gsub("_ref_Pacnes_phage_NC028967_cov.tsv", "", filename),
    vir_depth = meandepth,    # Mean sequencing depth for phage
    vir_cov = coverage,       # Coverage percentage for phage
    vir_reads = numreads      # Number of reads mapped to phage
  ) %>%
  select(vir_depth, vir_cov, vir_reads, SampleName)


## 2.2 Import Bacteria Mapping Data
bacteria_coverage <- read_csv("samtools_coverage_stats-cacnes.csv")

# Rename first column for clarity
names(bacteria_coverage)[1] <- "ref"

# Extract relevant columns and create sample names
bacteria_data <- bacteria_coverage %>%
  mutate(
    SampleName = gsub("_ref_Pacnes_C1_aligned.sorted.bam", "", filename),
    bac_depth = meandepth,    # Mean sequencing depth for bacteria
    bac_cov = coverage,       # Coverage percentage for bacteria
    bac_reads = numreads      # Number of reads mapped to bacteria
  ) %>%
  select(bac_depth, bac_cov, bac_reads, SampleName)


## 2.3 Import Metadata and Additional Data

# Global skin sample metadata
global_skin_metadata <- read_csv("global_skin_metadata.csv")

# Read lengths by study
global_skin_readlength <- read_csv("global_skin_metadata-readlength.csv")

# Phylogroup frequencies from Phlame results
phylogroup_freq <- read_csv("Cacnes_global_phylogroup_frequencies.csv")

# Bracken taxonomic abundance results
bracken_abundance <- read_csv("bracken_dcast-cacnes.csv")


# ---------------------------- #
# 3. Data Integration and VMR Calculation #
# ---------------------------- #

# Merge all datasets by SampleName
vmr_data <- phage_data %>%
  left_join(bacteria_data, by = "SampleName") %>%
  left_join(phylogroup_freq, by = "SampleName") %>%
  left_join(global_skin_metadata, by = "SampleName") %>%
  left_join(bracken_abundance, by = "SampleName") %>%
  group_by(SampleName) %>%
  mutate(
    vmr = vir_depth / bac_depth,            # Virus-to-microbe ratio
    vmr = ifelse(vir_cov >= 25, vmr, 0)    # Set VMR to 0 if phage coverage < 25%
  ) %>%
  # Remove unneeded metadata columns to clean dataset
  select(
    -c(Instrument, LibraryLayout, LibrarySelection, Site, `Public data?`, 
       `Run Number`, AvgSpotLen, City, Country, `United States, Europe or Canada`, 
       Ethnicity, Age, BioProject)
  ) %>%
  ungroup() %>%
  # Join read length metadata for further analysis
  left_join(global_skin_readlength)

# Optional: View the final combined data
# View(vmr_data)


# ---------------------------- #
# 4. Data Visualization        #
# ---------------------------- #

## 4.1 Prepare data for plotting
vmr_plot <- vmr_data

# Filter samples with sufficient total reads for reliable analysis
min_reads_threshold <- 1e6



################################################################################
# 
################################################################################

# Define coverage thresholds
coverage_thresholds <- seq(0, 100, by = 5)

# Calculate phage prevalence and confidence intervals per study and threshold
phage_prevalence_ci <- map_df(coverage_thresholds, function(thresh) {
  vmr_plot %>%
    filter(tot_reads >= 1e6) %>%
    mutate(has_phage = ifelse(vir_cov >= thresh, 1, 0)) %>%
    group_by(StudyName, BioProject) %>%
    summarise(
      n_samples = n(),
      n_phage = sum(has_phage),
      .groups = "drop"
    ) %>%
    mutate(
      phage_pos_frac = n_phage / n_samples,
      coverage_threshold = thresh
    ) %>%
    rowwise() %>%
    mutate(
      # Calculate 95% binomial confidence interval using Wilson method
      conf_int = list(binom.test(n_phage, n_samples)$conf.int)
    ) %>%
    mutate(
      conf_low = conf_int[[1]],
      conf_high = conf_int[[2]]
    ) %>%
    select(-conf_int)
})

# Plot with confidence intervals
prevalence_thres = ggplot(phage_prevalence_ci, 
                          aes(x = coverage_threshold, y = phage_pos_frac, 
                              color = BioProject, fill = BioProject)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.2, color = NA) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Phage Prevalence vs. Coverage Threshold by Study",
    x = "Phage Horizontal Coverage Threshold (%)",
    y = "Fraction of Samples with Phage",
    color = "Study",
    fill = "Study"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Save plot to a PDF file
ggsave(
  filename = "4_strain_coexistence/final_plots/prevalence_cov_threshold.pdf",
  plot = prevalence_thres ,
  device = "pdf",
  width = 8,    # adjust width (in inches) as needed
  height = 6,    # adjust height (in inches) as needed
  units = "in"
)




################################################################################
# 
################################################################################

# Define coverage thresholds
coverage_thresholds <- seq(0, 10, by = 1)

# Calculate phage prevalence and confidence intervals per study and threshold
bac_prevalence_ci <- map_df(coverage_thresholds, function(thresh) {
  vmr_plot %>%
    filter(tot_reads >= 1e6) %>%
    mutate(has_bac = ifelse(bac_depth >= thresh, 1, 0)) %>%
    group_by(StudyName, BioProject) %>%
    summarise(
      n_samples = n(),
      n_bac = sum(has_bac),
      .groups = "drop"
    ) %>%
    mutate(
      bac_pos_frac = n_bac / n_samples,
      coverage_threshold = thresh
    ) %>%
    rowwise() %>%
    mutate(
      # Calculate 95% binomial confidence interval using Wilson method
      conf_int = list(binom.test(n_bac, n_samples)$conf.int)
    ) %>%
    mutate(
      conf_low = conf_int[[1]],
      conf_high = conf_int[[2]]
    ) %>%
    select(-conf_int)
})

# Plot with confidence intervals
prevalence_thres = ggplot(bac_prevalence_ci, 
                          aes(x = coverage_threshold, y = bac_pos_frac, 
                              color = BioProject, fill = BioProject)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.2, color = NA) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Bacteria Prevalence vs. Depth Threshold by Study",
    x = "Bacteria Depth Threshold",
    y = "Fraction of Samples with Bacteria",
    color = "Study",
    fill = "Study"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Save plot to a PDF file
ggsave(
  filename = "4_strain_coexistence/final_plots/bac_prevalence_cov_threshold.pdf",
  plot = prevalence_thres ,
  device = "pdf",
  width = 8,    # adjust width (in inches) as needed
  height = 6,    # adjust height (in inches) as needed
  units = "in"
)



## 4.2 Plot 1: Percent of individuals with bacteria vs. log10 total sequenced bases per study
plot_bac_pos_frac <- ggscatter(
  right_join(
    # Calculate fraction of samples with bacteria presence per study
    vmr_plot %>%
      filter(tot_reads >= min_reads_threshold) %>%
      mutate(has_bac = ifelse(bac_depth >= 5, "Y", "N")) %>%
      group_by(StudyName, BioProject, has_bac) %>%
      tally() %>%
      ungroup() %>%
      group_by(StudyName, BioProject) %>%
      mutate(sample_counts = sum(n)) %>%
      filter(has_bac == "Y") %>%
      mutate(bac_pos_frac = n / sample_counts) %>%
      select(bac_pos_frac),
    # Median sequencing metrics per study
    vmr_plot %>%
      filter(tot_reads >= min_reads_threshold) %>%
      group_by(StudyName, BioProject, Region, ReadLength) %>%
      summarise(
        median_tot_reads = median(tot_reads),
        median_ca = median(bracken_cacnes),
        median_vmr = median(vmr),
        sample_counts = n(),
        .groups = "drop"
      )
  ) %>%
    mutate(median_tot_bps = median_tot_reads * ReadLength),
  y = "bac_pos_frac",
  x = "median_tot_bps",
  size = "sample_counts",
  alpha = 0.7,
  color = "Region",
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  xlab = "Median # Sequenced Bases",
  ylab = "% C. acnes Bacteria"
) +
  scale_x_log10() +
  ylim(0, 1) +
  theme(legend.position = "right")


## 4.3 Plot 2: Percent of individuals with phage (≥25% coverage) vs. log10 total sequenced bases per study
plot_phage_pos_frac <- ggscatter(
  right_join(
    vmr_plot %>%
      filter(tot_reads >= min_reads_threshold) %>%
      mutate(has_phage = ifelse(vir_cov >= 25, "Y", "N")) %>%
      group_by(StudyName, BioProject, has_phage) %>%
      tally() %>%
      ungroup() %>%
      group_by(StudyName, BioProject) %>%
      mutate(sample_counts = sum(n)) %>%
      filter(has_phage == "Y") %>%
      mutate(phage_pos_frac = n / sample_counts) %>%
      select(phage_pos_frac),
    vmr_plot %>%
      filter(tot_reads >= min_reads_threshold) %>%
      group_by(StudyName, BioProject, Region, ReadLength) %>%
      summarise(
        median_tot_reads = median(tot_reads),
        median_ca = median(bracken_cacnes),
        median_vmr = median(vmr),
        sample_counts = n(),
        .groups = "drop"
      )
  ) %>%
    mutate(median_tot_bps = median_tot_reads * ReadLength),
  y = "phage_pos_frac",
  x = "median_tot_bps",
  size = "sample_counts",
  alpha = 0.7,
  color = "Region",
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  xlab = "Median # Sequenced Bases",
  ylab = "% C. acnes Phage"
) +
  scale_x_log10() +
  ylim(0, 1) +
  theme(legend.position = "right")


## 4.5 Arrange plots vertically
pdf("phage_bacteria_prevalence_plots.pdf", width = 7, height = 10)  # Adjust width/height as needed
grid.arrange(plot_bac_pos_frac, plot_phage_pos_frac,ncol = 1)
dev.off()




# ---------------------------- #
# 4.6. Additional Plot: VMR #
# ---------------------------- #

## 4.7 Scatter plot with VMR
plot_vmr <- vmr_plot %>%
  filter(tot_reads >= min_reads_threshold & bac_depth >= 5) %>%
  ggplot(aes(y = vir_depth, x = bac_depth, shape = BioProject, color = Region)) +
  geom_point(alpha = 0.75) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype="dashed", alpha=0.35) + 
  geom_abline(intercept = mean(vmr_plot_filtered$log_diff, na.rm = TRUE), slope = 1, color = "black", linetype="solid", alpha=0.35) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"), name = "Region") +
  scale_shape_manual(values = c(1, 15, 4, 6, 18, 17, 0, 15)) +
  scale_x_continuous(trans = "log10") + # Log-transform x-axis
  scale_y_continuous(trans = "log10") + # Log-transform y-axis
  theme_pubr() +
  theme(legend.position = "left") +
  labs(y = "phage coverage", x = "bacteria coverage")

ggsave(
  filename = "vmr_scatter.pdf",
  plot = plot_vmr,
  device = "pdf",
  width = 8,    # adjust width (in inches) as needed
  height = 6,    # adjust height (in inches) as needed
  units = "in"
)



# ---------------------------- #
# Resistance Frequency vs VMR #
# ---------------------------- #

## 6.1 Scatter plot with linear regression of resistance frequency vs. VMR
plot_resistance_vs_vmr <- vmr_plot %>%
  filter(tot_reads >= min_reads_threshold & bac_depth >= 5) %>%
  mutate(res_freq = D + H + K) %>%
  ggplot(aes(y = res_freq, x = vmr, shape = BioProject )) + #, color = ifelse(vmr >= 1, "Y", "N"))) +
  geom_point(alpha = 0.5) +
  # scale_color_manual(values = c("grey20", "tomato"), name = "VMR ≥ 1") +
  scale_shape_manual(values = c(1, 15, 4, 6, 18, 17, 0, 15)) +
  scale_x_continuous(trans = "log10") + # Log-transform x-axis
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  # Uncomment below for regression equation and correlation stats
  # stat_regline_equation(aes(group = 1, label = ..rr.label..), label.x.npc = "left", label.y = 1.2) +
  # stat_cor(aes(group = 1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
  #          method = "pearson", label.x.npc = "left", label.y = 1.25, p.accuracy = 0.001,
  #          r.accuracy = 0.01, show.legend = FALSE) +
  theme_pubr() +
  theme(legend.position = "bottom") +
  labs(y = "Resistance Frequency", x = "VMR")


## 6.2 Density plots for VMR and resistance frequency, colored by VMR threshold
density_vmr <- ggplot(vmr_plot %>% filter(tot_reads >= min_reads_threshold & bac_depth >= 5),
                      aes(x = vmr, color = "white" )) + #, fill = ifelse(vmr >= 1, "Y", "N"))) +
  geom_histogram(alpha = 0.9) +
  theme_void() +
  theme(legend.position = "none") +
  scale_x_continuous(trans = "log10") +
  # scale_fill_manual(values = c("grey20", "tomato")) +
  scale_color_manual(values = "white")

density_res_freq <- ggplot(vmr_plot %>% filter(tot_reads >= min_reads_threshold & bac_depth >= 5) %>% 
                             mutate(res_freq = D + H + K),
                           aes(x = res_freq, color = "white" )) + #, fill = ifelse(vmr >= 1, "Y", "N"))) +
  geom_histogram(alpha = 0.9) +
  theme_void() +
  theme(legend.position = "none") +
  coord_flip() +
  # scale_fill_manual(values = c("grey20", "tomato")) +
  scale_color_manual(values = "white")


## 6.3 Combine density plots and scatter plot with patchwork
combined_plot <- density_vmr + plot_spacer() + plot_resistance_vs_vmr + density_res_freq +
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
combined_plot
# Save the combined plot to a PDF file
ggsave(
  filename = "resistance_vs_vmr_density_plot-nocolor.pdf",
  plot = combined_plot,
  device = "pdf",
  width = 8,    # adjust width (in inches) as needed
  height = 8,    # adjust height (in inches) as needed
  units = "in"
)




################################################################################
# Plotting VMR by human metadata
################################################################################

### SEX

# First, filter your data and summarize means
vmr_summary <- vmr_plot %>%
  filter(tot_reads >= 1e6 & bac_depth >= 5) %>%
  group_by(Sex) %>%
  summarize(mean_vmr = mean(vmr, na.rm=TRUE),
            n = n())

# Perform statistical test (Wilcoxon test on vmr)
stat_test <- vmr_plot %>%
  filter(tot_reads >= 1e6 & bac_depth >= 5) %>%
  wilcox.test(vmr ~ Sex, data = .)

p_val <- stat_test$p.value

# Create the plot
ggplot(vmr_summary, aes(x = as.character(Sex), y = mean_vmr, group = 1)) +
  
  # Line connecting means
  geom_line(size = 1) +
  
  # Points at means
  geom_point(size = 5, shape = 21, fill = "white", stroke = 1.5) +
  
  # Add bracket and p-value annotation
  geom_signif(comparisons = list(unique(as.character(vmr_summary$Sex))),
              annotations = paste0("P = ", signif(p_val, 2)),
              y_position = max(vmr_summary$mean_vmr)*1.1,
              tip_length = 0.02, textsize = 5) +
  
  # Adjust scales and theme
  scale_y_continuous(trans = "log10", expand = expansion(mult = c(0.1, 0.2)), limits = c(1e-1, NA),
                     breaks = c(0.1, 0.3, 1.0, 3.0),
                     labels = c("0.1", "0.3", "1.0", "3.0")) +
  labs(x = "Sex Status", y = "Virus-to-Microbe Ratio (VMR)") +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        plot.title = element_text(hjust=0.5, size=15, face="bold")) +
  ggtitle("VMR by Sex Status")


### AGE

# First, filter your data and summarize means
vmr_summary <- vmr_plot %>%
  filter(tot_reads >= 1e6 & bac_depth >= 5) %>%
  mutate(Age_Category = ifelse(ApproxAge == -1, NA,
                               ifelse(ApproxAge != -1 & ApproxAge >= 40, "40+", "Under 40"))) %>%
  filter(!is.na(Age_Category)) %>%
  group_by(Age_Category) %>%
  summarize(mean_vmr = mean(vmr, na.rm=TRUE),
            n = n())

# Perform statistical test (Wilcoxon test on vmr)
stat_test <- vmr_plot %>%
  filter(tot_reads >= 1e6 & bac_depth >= 5) %>%
  mutate(Age_Category = ifelse(ApproxAge == -1, NA,
                               ifelse(ApproxAge != -1 & ApproxAge >= 40, "40+", "Under 40"))) %>%
  filter(!is.na(Age_Category)) %>%
  wilcox.test(vmr ~ Age_Category, data = .)

p_val <- stat_test$p.value

# Create the plot
ggplot(vmr_summary, aes(x = as.character(Age_Category), y = mean_vmr, group = 1)) +
  
  # Line connecting means
  geom_line(size = 1) +
  
  # Points at means
  geom_point(size = 5, shape = 21, fill = "white", stroke = 1.5) +
  
  # Add bracket and p-value annotation
  geom_signif(comparisons = list(unique(as.character(vmr_summary$Age_Category))),
              annotations = paste0("P = ", signif(p_val, 2)),
              y_position = max(vmr_summary$mean_vmr)*1.1,
              tip_length = 0.02, textsize = 5) +
  
  # Adjust scales and theme
  scale_y_continuous(trans = "log10", expand = expansion(mult = c(0.1, 0.2)), limits = c(1e-1, NA),
                     breaks = c(0.1, 0.3, 1.0, 3.0),
                     labels = c("0.1", "0.3", "1.0", "3.0")) +
  labs(x = "Age Status", y = "Virus-to-Microbe Ratio (VMR)") +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        plot.title = element_text(hjust=0.5, size=15, face="bold")) +
  ggtitle("VMR by Age Status")

