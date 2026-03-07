# Set library path safely
.libPaths("F:/SDXC/STORE N GO/R/R-4.0.2/win-library/4.0")

# Load packages using pacman for reliability
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readxl, readr, reshape2, dplyr, gplots, Heatplus, vegan, RColorBrewer,
  tidyr, gtools, stringr, tidyverse, ComplexHeatmap, magick, viridis, 
  Hotelling, remotes, ggnewscale, ggtreeExtra, jmvcore, gtsummary, 
  officer, fpc, data.table, ggpubr, ggstatsplot, ggtree, grid
)

# Ensure reproducibility
set.seed(123)

# Define data source paths
data_source_details <- list(
  NottinghamJuly2025 = data.frame(
    metaphlan = "//mkvx341/Projects2/CoDiet/liam walsh/NottinghamJuly2025/results/04_short_read_taxonomic_profiling/metaphlan_results_Notthingham_v1.csv",
    preprocessing = "//mkvx341/Projects2/CoDiet/liam walsh/NottinghamJuly2025/results/01_preprocessing/total_fascinar_results.tsv",
    output_path="//mkvx341/Projects2/CoDiet/liam walsh/NottinghamJuly2025/officer_output/Reports"
   # officer_path = "//mkvx341/Projects2/CoDiet/liam walsh/NottinghamJuly2025/officer_output/",
    #figure_path = "//mkvx341/Projects2/CoDiet/liam walsh/NottinghamJuly2025/plots/"
  )
)

# Load input data
data <- list()
dataset <- "NottinghamJuly2025"
paths <- data_source_details[[dataset]]

for (file_key in names(paths)) {
  if (file_key %in% c("figure_path", "officer_path","output_path")) next
  file_path <- paths[[file_key]]
  
  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    data[[file_key]] <- read_csv(file_path)
  } else if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
    data[[file_key]] <- read_delim(file_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  } else {
    warning(paste("Unknown file format for:", file_path))
  }
}

# Clean and harmonise preprocessing metadata
data[["preprocessing"]] <- data[["preprocessing"]] %>%
  separate(`Filename,#samplename`, into = c("X1", "X2", "X3", "sample_id_and_contig"), sep = "--", extra = "merge") %>%
  separate(sample_id_and_contig, into = c("sample_id", "file_name"), sep = ",", extra = "merge") %>%
  filter(sample_id == file_name) %>%
  select(-c("X1", "X2", "X3","sample_id"  ))


# Standardise and filter species profile
data$metaphlan <- data$metaphlan %>% rename_with(~ gsub("_metaphlan", "", .))
species_profile <- data$metaphlan %>%
  column_to_rownames("clade_name") %>%
  t() %>%
  as.data.frame() %>%
  select(which(grepl("\\|t__", names(.)))) %>%
  t() %>%
  as.data.frame()



# Filter: Keep species with max relative abundance > 0.1
maxab_species <- apply(species_profile, 1, max, na.rm = TRUE)
species_profile <- species_profile[names(maxab_species[maxab_species > 0.1]), ]
species_profile <- species_profile[, colSums(species_profile) > 0]

# Calculate prevalence
species_prevalence <- data.frame(
  sum_species = rowSums(species_profile > 0.1)
) %>%
  rownames_to_column("specie") %>%
  filter(sum_species > ncol(species_profile) * 0.1) %>%
  cbind(str_split_fixed(.$specie, "\\|", 8)) %>%
  rename(
    kingdom = `1`, phylum = `2`, class = `3`, order = `4`,
    family = `5`, genus = `6`, species = `7`, strain = `8`
  )

# Alpha diversity
alpha_summary <- data.frame(
  Richness = hillR::hill_taxa(t(species_profile), q = 0),
  Shannon = hillR::hill_taxa(t(species_profile), q = 1),
  Simpson = hillR::hill_taxa(t(species_profile), q = 2)
)

# Beta diversity
data_matrix <- as.matrix(t(species_profile))
beta_summary <- hill_taxa_parti_pairwise(
  data_matrix, q = 2, show_warning = TRUE, .progress = TRUE, rel_then_pool = FALSE
)

dist <- bind_rows(
  beta_summary %>% select(site1, site2, TD_beta),
  beta_summary %>% select(site2, site1, TD_beta) %>% rename(site1 = 1, site2 = 2)
) %>%
  pivot_wider(names_from = site2, values_from = TD_beta) %>%
  column_to_rownames("site1") %>%
  replace(is.na(.), 0)


# Ensure sample IDs match across species_profile and alpha_summary
common_samples <- intersect(colnames(species_profile), rownames(alpha_summary))

# Reorder species_profile and alpha_summary to match sample order
species_profile <- species_profile[, common_samples]
alpha_summary <- alpha_summary[common_samples, ]

# Confirm dimensions match
stopifnot(all(colnames(species_profile) == rownames(alpha_summary)))

# reporting time 
sample_summary <- data.frame(
  SampleID = colnames(species_profile),
  TotalReads = colSums(species_profile),
  Richness = alpha_summary$Richness,
  Shannon = alpha_summary$Shannon,
  Simpson = alpha_summary$Simpson
)

#  sanity check: do sums match?
message("✅ Sample summary generated with ", nrow(sample_summary), " samples.")
message("Total read sum: ", sum(sample_summary$TotalReads))

# Define path
report_dir <- data_source_details$NottinghamJuly2025$output_path

# Create directory if it doesn't exist
if (!dir.exists(report_dir)) {
  dir.create(report_dir, recursive = TRUE)
}

# 1. Save cleaned preprocessing metadata
preprocessing_out_path <- file.path(report_dir, "sequencing_run_summary_statistics.csv")
write_csv(data[["preprocessing"]], preprocessing_out_path)

# 2. Save filtered species profile
species_out_path <- file.path(report_dir, "species_profile_filtered.csv")
write_csv(as.data.frame(species_profile) %>% rownames_to_column("Taxa"), species_out_path)

# 3. Save unfiltered species profile
unfiltered_species_out_path <- file.path(report_dir, "species_profile_unfiltered.csv")
write_csv(as.data.frame(data$metaphlan ), unfiltered_species_out_path)

# 4. Save species prevalence data
prevalence_out_path <- file.path(report_dir, "species_prevalence_breakdown.csv")
write_csv(species_prevalence, prevalence_out_path)

# 5. Save alpha diversity summary
alpha_out_path <- file.path(report_dir, "alpha_diversity.csv")
write_csv(sample_summary %>% rownames_to_column("sample_id"), alpha_out_path)

# 6. Save beta diversity matrix
beta_out_path <- file.path(report_dir, "beta_diversity_TD_q2_matrix.csv")
write_csv(beta_summary, beta_out_path)

# 7. Summary log file of exports
summary_file <- file.path(report_dir, "summary_of_exports.txt")
cat(
  "Report contents generated on:", Sys.time(), "\n\n",
  "- preprocessing_cleaned_metadata.csv: Cleaned preprocessing metadata\n",
  "- species_profile_filtered.csv: Filtered species-level relative abundance matrix\n",
  "- species_profile_unfiltered.csv: Unfiltered species-level matrix\n",
  "- species_prevalence_breakdown.csv: Prevalence of species across samples\n",
  "- alpha_diversity_summary.csv: Richness, Shannon, and Simpson indices\n",
  "- beta_diversity_TD_q2_matrix.csv: Pairwise total dissimilarity (q = 2)\n",
  file = summary_file
)

message("✅ All reports exported to: ", report_dir)






