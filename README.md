# NTU / Teagasc Turmeric Study

This repository contains the analysis code, processed input tables, derived results, and manuscript-support material for the Nottingham / Teagasc turmeric study in elite soccer players.

The repository is organized as an analysis workbench rather than an R package. The useful unit is the manuscript workflow: diversity derivation, Bayesian alpha-diversity modelling, beta-diversity testing, MaAsLin3 screening, exploratory paired effect-size analysis, and exploratory GPS/injury associations.

## Main Questions

1. Does turmeric exposure associate with alpha-diversity change?
2. Does turmeric exposure associate with global beta-diversity separation?
3. Which taxa and pathways show the clearest small-sample change-score signals?
4. Do those microbiome signals co-vary with external load or injury/illness burden?

## Repository Layout

```text
scripts/                                  Core analysis scripts
results/Processed/                        Processed analysis inputs and derived outputs
results/Reports/                          Legacy export tables still used by some scripts
results/Raw/                              Raw merged microbiome tables used as upstream inputs
figures/                                  Manuscript-facing and analysis figures
docs/                                     Workflow and data-dependency documentation
README.md                                 GitHub landing page
```

## Core Scripts

| Script | Role | Primary inputs | Primary outputs |
|---|---|---|---|
| `scripts/1. Diversity_index_formation_notthingham.R` | Derives alpha/beta diversity from microbiome profile tables | `results/Raw/paired_without_unclassified_abundance_table_species.tsv`, `results/Raw/HUMAnN_merged_pathabundance_cpm.tsv`, `results/Raw/analysis_metadata_current.csv` | `results/Processed/02_taxonomic_profiling/`, `results/Processed/03_functional_profiling/`, `results/Processed/04_diversity_analysis/` |
| `scripts/2. Turmeric & Gut Microbiome Diversity –  Bayesian Mixed Model Results.R` | Bayesian mixed-effects models for richness, Shannon, Simpson | Alpha-diversity CSVs in `results/Processed/04_diversity_analysis/`, `results/Raw/analysis_metadata_current.csv` | `results/Processed/05_statistical_analysis/`, `figures/05_statistical_analysis/` |
| `scripts/3. Beta diversity.R` | Exact CLR-delta test, PERMANOVA/PERMDISP sensitivity analyses, ordinations | Filtered and unfiltered profiles, beta-distance tables in `results/Processed/`, `results/Raw/analysis_metadata_current.csv` | `results/Processed/05_statistical_analysis/`, `figures/05_statistical_analysis/` |
| `scripts/4. maaslin3_results_SMALL_SAMPLE.R` | Small-sample MaAsLin3 validation/screening layer | Processed species and functional profiles plus MaAsLin3 result tables | `results/Processed/05_statistical_analysis/differential_abundance/` |
| `scripts/4e_honest_smallN_effectsize_taxa_function_pipeline.R` | Exploratory paired strain/genus/pathway effect-size pipeline | `results/Processed/02_taxonomic_profiling/species_profile_filtered.csv`, pathway abundance tables, participant metadata encoded in sample IDs | `results/Processed/05_statistical_analysis/honest_small_n_strain/` |
| `scripts/5_microbiome_gps_injury_associations.R` | Exploratory microbiome vs GPS/injury association analysis | Alpha diversity, strain/genus/pathway outputs from `4e`, GPS/injury workbook or committed prepared cohort CSV | `results/Processed/05_statistical_analysis/gps_injury_associations/` |
| `scripts/export_gps_microbiome_alignment_workbook.py` | Exports an audit workbook showing GPS + microbiome alignment | GPS association outputs + participant-level microbiome deltas | `results/Processed/05_statistical_analysis/gps_injury_associations/` |

## Data Availability in This Repo

The processed inputs needed for the manuscript-facing scripts are committed in the repository.

Key directories:

- `results/Processed/02_taxonomic_profiling/`
  - filtered and unfiltered species-profile tables
- `results/Processed/03_functional_profiling/`
  - filtered pathway profile table
- `results/Processed/04_diversity_analysis/`
  - alpha-diversity tables and beta-distance matrices
- `results/Processed/05_statistical_analysis/honest_small_n_strain/`
  - exploratory strain/genus/pathway outputs
- `results/Processed/05_statistical_analysis/gps_injury_associations/`
  - cleaned GPS cohort tables and association outputs

The only non-repo raw dependency used by one script is the original GPS/injury Excel workbook. The current GPS script can also run from the committed cleaned cohort CSV in `results/Processed/05_statistical_analysis/gps_injury_associations/` when the workbook is not available.

For a file-by-file dependency map, see [docs/data_manifest.md](docs/data_manifest.md).

## Workflow Guide

See [docs/analysis_workflow.md](docs/analysis_workflow.md) for script ordering and high-level workflow.

Current rerun order:

1. `scripts/1. Diversity_index_formation_notthingham.R`
2. `scripts/2. Turmeric & Gut Microbiome Diversity –  Bayesian Mixed Model Results.R`
3. `scripts/3. Beta diversity.R`
4. `scripts/4. maaslin3_results_SMALL_SAMPLE.R`
5. `scripts/4e_honest_smallN_effectsize_taxa_function_pipeline.R`
6. `scripts/5_microbiome_gps_injury_associations.R`
7. `scripts/export_gps_microbiome_alignment_workbook.py`

## Environment Notes

Before rerunning locally, create a `.Renviron` from `.Renviron.example` and set values as needed:

- `GITHUB_PAT` or `GITHUB_TOKEN`
  - only needed for authenticated GitHub API reads or optional upload helpers
- `ENABLE_GITHUB_UPLOAD=1`
  - opt-in flag; default workflow is local-first
- `NTU_R_LIB`
  - optional override for user R library path
- `GPS_INJURY_XLSX`
  - optional path to the original GPS/injury workbook

Most manuscript-facing scripts now save locally under `results/Processed/` and `figures/` unless explicit GitHub upload is enabled.

