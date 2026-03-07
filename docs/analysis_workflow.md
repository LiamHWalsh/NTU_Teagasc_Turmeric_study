# Analysis Workflow Guide

This document explains how the main manuscript-facing analysis scripts fit together.

For exact input and output file paths, see [data_manifest.md](data_manifest.md).

## 1. Diversity metric generation

Script:
- `scripts/1. Diversity_index_formation_notthingham.R`

Purpose:
- derives species- and function-level alpha diversity
- derives beta-diversity summaries and distance matrices
- writes the processed taxonomic and functional tables used downstream

Primary output directories:
- `results/Processed/02_taxonomic_profiling/`
- `results/Processed/03_functional_profiling/`
- `results/Processed/04_diversity_analysis/`

## 2. Bayesian alpha-diversity modelling

Script:
- `scripts/2. Turmeric & Gut Microbiome Diversity –  Bayesian Mixed Model Results.R`

Purpose:
- fits Bayesian mixed-effects models for richness, Shannon, and Simpson diversity
- produces the alpha-diversity results used in the manuscript

Depends on:
- diversity outputs from script 1

## 3. Beta-diversity analysis

Script:
- `scripts/3. Beta diversity.R`

Purpose:
- evaluates species- and pathway-level beta diversity
- current manuscript-facing logic uses:
  - exact permutation testing on participant-level CLR deltas
  - PERMANOVA as a secondary sensitivity analysis
  - PERMDISP to qualify interpretation
  - ordination for visualization only

Depends on:
- processed species and pathway profile tables
- beta-distance matrices from script 1

## 4. MaAsLin3 small-sample differential modelling

Script:
- `scripts/4. maaslin3_results_SMALL_SAMPLE.R`

Purpose:
- runs the model-based differential layer alongside the exploratory effect-size pipeline
- validates and visualizes small-sample MaAsLin3 outputs

Depends on:
- processed species and functional profiles
- MaAsLin3 result tables

## 5. Exploratory paired effect-size pipeline

Script:
- `scripts/4e_honest_smallN_effectsize_taxa_function_pipeline.R`

Purpose:
- computes participant-level baseline-to-6-month deltas
- compares treatment vs control deltas using mean difference, bootstrap CIs, Cohen's d, probability of superiority, LOO sensitivity, and CLR sensitivity
- produces the current strain/genus/pathway exploratory results used in the manuscript

Current configuration of note:
- taxonomic resolution: strain (`t__` entries)
- genus panel: literature-anchored
  - `Bacteroides`
  - `Bifidobacterium`
  - `Alistipes`
  - `Parabacteroides`
  - `Blautia`
- inference carried forward on:
  - all retained strains within those genera
  - all retained pathways after prevalence/abundance filtering

Current output directory:
- `results/Processed/05_statistical_analysis/honest_small_n_strain/`

Most important files there:
- `strain_effect_sizes.csv`
- `genus_effect_sizes.csv`
- `pathway_effect_sizes.csv`
- `taxonomy_function_candidate_pairs.csv`

## 6. GPS / injury association analysis

Script:
- `scripts/5_microbiome_gps_injury_associations.R`

Purpose:
- tests whether participant-level microbiome deltas co-vary with training load and injury/illness burden
- uses both frequentist rank correlation and Bayesian rank-model sensitivity analysis

Important implementation details:
- restricted to the same paired `n = 11` treatment/control subset used in the exploratory microbiome pipeline
- microbiome outcome delta is defined as `month6 - mean(Baseline, 2_weeks)`
- GPS metadata uses two explicit windows:
  - pre turmeric = 14-week pre-season up to Visit 2
  - post turmeric = 24-week season up to Visit 3
- if the raw GPS workbook is unavailable, the script can fall back to the committed cleaned cohort table in `results/Processed/05_statistical_analysis/gps_injury_associations/`
- workbook-header validation is retained when the raw Excel file is used
- all retained `4e` signals are carried forward, not a top-N subset

Current output directory:
- `results/Processed/05_statistical_analysis/gps_injury_associations/`

Useful outputs:
- `association_summary.txt`
- `gps_injury_clean_all_participants.csv`
- `gps_injury_analysis_cohort_n11.csv`
- `participant_alpha_deltas_n11.csv`
- `participant_signal_deltas_n11.csv`
- `gps_microbiome_alignment_bundle_n11_v2.xlsx`

## 7. Manuscript-facing text outputs

Main manuscript-support file:
- `results/Processed/06_reports/manuscript_revision_v5_exploratory_sections.md`

This file is currently the most useful single location for:
- exploratory Methods text
- exploratory Results text
- Discussion language aligned to the latest reruns

## 8. Practical rerun order

If you want to rerun the current manuscript stack, the relevant order is:

1. `scripts/1. Diversity_index_formation_notthingham.R`
2. `scripts/2. Turmeric & Gut Microbiome Diversity –  Bayesian Mixed Model Results.R`
3. `scripts/3. Beta diversity.R`
4. `scripts/4. maaslin3_results_SMALL_SAMPLE.R`
5. `scripts/4e_honest_smallN_effectsize_taxa_function_pipeline.R`
6. `scripts/5_microbiome_gps_injury_associations.R`
7. `scripts/export_gps_microbiome_alignment_workbook.py`

## 9. Repository status

This repository is usable and now exposes the processed inputs required for the manuscript-facing analyses, but it is still an analysis repository rather than a polished software package. Some older scripts remain legacy or duplicated, and several filenames still contain spaces and punctuation.
