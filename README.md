
## Phage-mediated lysis does not determine *Cutibacterium acnes* colonization on human skin

---

## About
This repository contains code, scripts, and supplementary materials necessary to reproduce the analyses described in:

> Tripp, A. Delphine, et al. "Phage-mediated lysis does not determine Cutibacterium acnes colonization on human skin." bioRxiv (2025): 2025-09.
> DOI: [https://doi.org/10.1101/2025.09.09.675206](https://doi.org/10.1101/2025.09.09.675206)

This project is part of my doctoral research in the [Lieberman Lab](https://lieberman.science/research/).


---

## Usage
To reproduce the results from the paper:

1. **Prepare datasets**  
   Download raw genomic and metagenomic data from the sources listed in the paper (Tables S1 and S5).

2. **Run genomic analyses**  

   These scripts performs genome assembly, annotation, defense system detection, and phylogenetic reconstruction.

```bash
isolate_genome_analysis/Crispr_cas_spacer_analysis/ # CRISPR array and spacer analyses
isolate_genome_analysis/DefenseFinder_Cacnes_genomes/  # Anti-phage defense systems detection in C. acnes isolate genomes
isolate_genome_analysis/DefenseFinder_Wiki_Resource/  # Anti-phage defense systems in other skin-associated genera
isolate_genome_analysis/Gene_gain_loss_analysis/  # Ancestral trait reconstruction of C. acnes accessory genome
isolate_genome_analysis/MGE_presence_analysis/  # Carriage of prophages and plasmids across C. acnes isolate genomes
```

3. **Run metagenomic analyses**  

   These scripts processes metagenomic reads, calculates virus-to-microbe ratios, and estimates strain-level diversity using PHLAME.

```bash
cacnes_skin_metagenomes_analysis/1_bracken_plots  # Bracken results
cacnes_skin_metagenomes_analysis/2_phlame_plots  # PHLAME results
cacnes_skin_metagenomes_analysis/3_vmr_plots  # Virus-to-mircobe ratio results
```

4. **Run growth curve and phage infectivity assay analyses**  

   This script processes OD600 time series data and calculates growth rate.

```bash
liquid_phage_infection_assay/growth_curve.py # OD600 and growth data
```
   
---

## Data Availability
All original code is available in this repository:  
[https://github.com/adtripp/Tripp_et_al_cacnes_phage_defense](https://github.com/adtripp/Tripp_et_al_cacnes_phage_defense)

Datasets used in this study are available from:
- Natural isolate genomes (Lieberman Lab collections) — see Tables S1–S2.
- Public metagenomic datasets from NCBI SRA — see Table S5.
- Defense system annotations — see Table S2.
- CRISPR spacers — see Table S4.

Any additional information required to reanalyze the data is available from the lead contact upon request.

---

## Citation
If you use this code or data in your research, please cite:

```
Tripp, A. Delphine, et al. "Phage-mediated lysis does not determine Cutibacterium acnes colonization on human skin." bioRxiv (2025): 2025-09.
```

---

## License
This project is licensed under the MIT License.

---


