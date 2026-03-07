# Data Manifest

This document maps the manuscript-facing scripts to the data files they read and the outputs they write.

## 1. Diversity derivation

Script:
- `scripts/1. Diversity_index_formation_notthingham.R`

Primary inputs:
- `results/Raw/metaphlan_results_Notthingham_v1.csv`
- `results/Raw/HUMAnN_merged_pathabundance_cpm.tsv`

Primary outputs:
- `results/Processed/02_taxonomic_profiling/species_profile_unfiltered.csv`
- `results/Processed/02_taxonomic_profiling/species_profile_filtered.csv`
- `results/Processed/02_taxonomic_profiling/species_prevalence_breakdown.csv`
- `results/Processed/03_functional_profiling/functional_profile_filtered.csv`
- `results/Processed/03_functional_profiling/functional_prevalence_breakdown.csv`
- `results/Processed/04_diversity_analysis/species_alpha_diversity.csv`
- `results/Processed/04_diversity_analysis/functional_alpha_diversity.csv`
- `results/Processed/04_diversity_analysis/species_beta_diversity_distance_matrix.csv`
- `results/Processed/04_diversity_analysis/functional_beta_diversity_distance_matrix.csv`

## 2. Bayesian alpha-diversity modelling

Script:
- `scripts/2. Turmeric & Gut Microbiome Diversity –  Bayesian Mixed Model Results.R`

Primary inputs:
- `results/Processed/04_diversity_analysis/species_alpha_diversity.csv`
- `results/Processed/04_diversity_analysis/functional_alpha_diversity.csv`

Primary outputs:
- `results/Processed/05_statistical_analysis/combined_alpha_diversity_results.csv`
- figure outputs under `figures/05_statistical_analysis/`

## 3. Beta diversity

Script:
- `scripts/3. Beta diversity.R`

Primary inputs:
- `results/Processed/02_taxonomic_profiling/species_profile_filtered.csv`
- `results/Processed/03_functional_profiling/functional_profile_filtered.csv`
- `results/Processed/04_diversity_analysis/species_beta_diversity_distance_matrix.csv`
- `results/Processed/04_diversity_analysis/functional_beta_diversity_distance_matrix.csv`

Primary outputs:
- statistical summaries in `results/Processed/05_statistical_analysis/`
- ordination and dispersion figures in `figures/05_statistical_analysis/`

## 4. MaAsLin3 small-sample validation

Script:
- `scripts/4. maaslin3_results_SMALL_SAMPLE.R`

Primary inputs:
- `results/Processed/02_taxonomic_profiling/species_profile_filtered.csv`
- `results/Processed/03_functional_profiling/functional_profile_filtered.csv`
- `results/Processed/05_statistical_analysis/differential_abundance/species_maaslin3_all_results.csv`
- `results/Processed/05_statistical_analysis/differential_abundance/functional_maaslin3_all_results.csv`

Primary outputs:
- `results/Processed/05_statistical_analysis/differential_abundance/`
- supporting figures under `figures/05_statistical_analysis/differential_abundance/`

## 5. Exploratory paired effect-size pipeline

Script:
- `scripts/4e_honest_smallN_effectsize_taxa_function_pipeline.R`

Primary inputs:
- `results/Processed/02_taxonomic_profiling/species_profile_filtered.csv`
- `results/Reports/species_profile_filtered.csv` (legacy fallback)
- `results/05_short_read_functional_profiling/HUMAnN_merged_pathabundance_cpm.tsv`
- `results/Reports/HUMAnN_merged_pathabundance_cpm.tsv` (legacy fallback)

Primary outputs:
- `results/Processed/05_statistical_analysis/honest_small_n_strain/strain_effect_sizes.csv`
- `results/Processed/05_statistical_analysis/honest_small_n_strain/genus_effect_sizes.csv`
- `results/Processed/05_statistical_analysis/honest_small_n_strain/pathway_effect_sizes.csv`
- `results/Processed/05_statistical_analysis/honest_small_n_strain/taxonomy_function_candidate_pairs.csv`
- full-panel supplementary figures in the same directory

## 6. GPS / injury association layer

Script:
- `scripts/5_microbiome_gps_injury_associations.R`

Primary microbiome inputs:
- `results/Reports/alpha_diversity.csv`
- `results/Reports/species_profile_filtered.csv`
- `results/Reports/HUMAnN_merged_pathabundance_cpm.tsv`
- `results/Processed/05_statistical_analysis/honest_small_n_strain/strain_effect_sizes.csv`
- `results/Processed/05_statistical_analysis/honest_small_n_strain/genus_effect_sizes.csv`
- `results/Processed/05_statistical_analysis/honest_small_n_strain/pathway_effect_sizes.csv`

GPS/injury inputs:
- preferred raw workbook via `GPS_INJURY_XLSX`
- repo-local fallback:
  - `results/Processed/05_statistical_analysis/gps_injury_associations/gps_injury_clean_all_participants.csv`

Primary outputs:
- `results/Processed/05_statistical_analysis/gps_injury_associations/gps_injury_clean_all_participants.csv`
- `results/Processed/05_statistical_analysis/gps_injury_associations/gps_injury_analysis_cohort_n11.csv`
- `results/Processed/05_statistical_analysis/gps_injury_associations/participant_alpha_deltas_n11.csv`
- `results/Processed/05_statistical_analysis/gps_injury_associations/participant_signal_deltas_n11.csv`
- `results/Processed/05_statistical_analysis/gps_injury_associations/association_summary.txt`
- frequentist and Bayesian association result tables in the same directory

## 7. Alignment audit workbook

Script:
- `scripts/export_gps_microbiome_alignment_workbook.py`

Primary inputs:
- GPS/injury association outputs in `results/Processed/05_statistical_analysis/gps_injury_associations/`

Primary outputs:
- `results/Processed/05_statistical_analysis/gps_injury_associations/gps_microbiome_alignment_bundle_n11_v2.xlsx`

## Notes

- Sample metadata are encoded in the microbiome sample IDs and reconstructed inside the scripts.
- The processed tables in `results/Processed/` are the main GitHub-available analysis inputs.
- `results/Reports/` remains in use as a legacy export location for some scripts and is therefore kept visible in the repository.
