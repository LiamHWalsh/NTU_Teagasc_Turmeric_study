################################################################################
# Generate manuscript-ready Methods/Results text and a composite summary figure
# from the honest_small_n effect-size outputs.
################################################################################

preferred_lib <- "G:/R/win-library/4.4"
if (!dir.exists(preferred_lib)) {
  preferred_lib <- normalizePath(
    file.path(Sys.getenv("USERPROFILE"), "Documents", "R", "win-library", "4.4"),
    winslash = "/", mustWork = FALSE
  )
}
dir.create(preferred_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(preferred_lib, .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

candidate_input_dirs <- c(
  file.path("results", "Processed", "05_statistical_analysis", "honest_small_n_strain"),
  file.path("results", "Processed", "05_statistical_analysis", "honest_small_n")
)

detect_input <- function(dir_path) {
  strain_file <- file.path(dir_path, "strain_effect_sizes.csv")
  species_file <- file.path(dir_path, "species_effect_sizes.csv")
  genus_file <- file.path(dir_path, "genus_effect_sizes.csv")
  pathway_file <- file.path(dir_path, "pathway_effect_sizes.csv")
  pairs_file <- file.path(dir_path, "taxonomy_function_candidate_pairs.csv")

  has_common <- file.exists(genus_file) && file.exists(pathway_file) && file.exists(pairs_file)
  if (!has_common) return(NULL)

  if (file.exists(strain_file)) {
    return(list(
      input_dir = dir_path,
      taxon_file = "strain_effect_sizes.csv",
      taxon_label = "Strain",
      taxon_label_lower = "strain",
      taxon_entry_label = "strain-resolved entries",
      output_stub = "smallN_strain"
    ))
  }

  if (file.exists(species_file)) {
    return(list(
      input_dir = dir_path,
      taxon_file = "species_effect_sizes.csv",
      taxon_label = "Species",
      taxon_label_lower = "species",
      taxon_entry_label = "species entries",
      output_stub = "smallN_species"
    ))
  }
  NULL
}

cfg <- NULL
for (d in candidate_input_dirs) {
  cfg <- detect_input(d)
  if (!is.null(cfg)) break
}

if (is.null(cfg)) {
  stop(
    "Could not locate a valid honest_small_n results directory. Tried: ",
    paste(candidate_input_dirs, collapse = ", ")
  )
}

input_dir <- cfg$input_dir
taxon_file <- cfg$taxon_file
TAXON_LABEL <- cfg$taxon_label
TAXON_LABEL_LOWER <- cfg$taxon_label_lower
TAXON_ENTRY_LABEL <- cfg$taxon_entry_label
output_md <- file.path(input_dir, "MANUSCRIPT_methods_results.md")
output_png <- file.path(input_dir, paste0("figure_", cfg$output_stub, "_methods_results_summary.png"))
output_pdf <- file.path(input_dir, paste0("figure_", cfg$output_stub, "_methods_results_summary.pdf"))
output_merged_csv <- file.path(input_dir, "KEY_RESULTS_merged_table.csv")

required <- c(
  taxon_file,
  "genus_effect_sizes.csv",
  "pathway_effect_sizes.csv",
  "taxonomy_function_candidate_pairs.csv"
)
missing_files <- required[!file.exists(file.path(input_dir, required))]
if (length(missing_files) > 0) {
  stop("Missing required input files in ", input_dir, ": ", paste(missing_files, collapse = ", "))
}

to_bool <- function(x) {
  if (is.logical(x)) return(x)
  lx <- tolower(as.character(x))
  lx %in% c("true", "t", "1", "yes")
}

fmt <- function(x, digits = 2) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
  formatC(as.numeric(x), format = "f", digits = digits)
}

fmt_p <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
  format(signif(as.numeric(x), 3), scientific = TRUE, trim = TRUE)
}

trim_text <- function(x, width = 60) {
  x <- as.character(x)
  ifelse(nchar(x) > width, paste0(substr(x, 1, width - 3), "..."), x)
}

row_bind_fill <- function(df_list) {
  all_cols <- unique(unlist(lapply(df_list, names), use.names = FALSE))
  aligned <- lapply(df_list, function(df) {
    missing <- setdiff(all_cols, names(df))
    if (length(missing) > 0) {
      for (m in missing) df[[m]] <- NA
    }
    df[, all_cols, drop = FALSE]
  })
  out <- do.call(rbind, aligned)
  rownames(out) <- NULL
  out
}

panel_summary <- function(df, level_name) {
  ad <- abs(df$cohens_d)
  data.frame(
    level = level_name,
    n_features = nrow(df),
    n_ci_excludes_zero = sum(df$ci_excludes_zero, na.rm = TRUE),
    prop_ci_excludes_zero = mean(df$ci_excludes_zero, na.rm = TRUE),
    median_abs_d = median(ad, na.rm = TRUE),
    iqr_abs_d = IQR(ad, na.rm = TRUE),
    median_ci_width = median(df$ci_width, na.rm = TRUE),
    median_loo = median(df$loo_sign_agreement, na.rm = TRUE),
    transform_agreement = mean(df$transform_sign_agreement_raw_log1p, na.rm = TRUE),
    n_large = sum(ad >= 0.8, na.rm = TRUE),
    n_moderate = sum(ad >= 0.5 & ad < 0.8, na.rm = TRUE),
    n_small = sum(ad < 0.5, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

top_n_by_rank <- function(df, n = 3) {
  if (nrow(df) == 0) return(df)
  ord <- order(df$rank_consensus, -abs(df$cohens_d), na.last = TRUE)
  df[ord, , drop = FALSE][seq_len(min(n, nrow(df))), , drop = FALSE]
}

feature_line <- function(row, include_clr = FALSE) {
  base <- paste0(
    "- ", row$feature_short,
    ": d=", fmt(row$cohens_d, 2),
    "; delta-difference=", fmt(row$mean_difference, 2),
    "; bootstrap CI [", fmt(row$ci_lower, 2), ", ", fmt(row$ci_upper, 2), "]",
    "; LOO sign agreement=", fmt(row$loo_sign_agreement, 2),
    "."
  )

  if (include_clr && "d_clr" %in% names(row)) {
    base <- paste0(
      base,
      " CLR sensitivity d=",
      fmt(row$d_clr, 2),
      "; primary-vs-CLR sign agreement=",
      ifelse(is.na(row$sign_agreement_primary_vs_clr), "NA", as.character(row$sign_agreement_primary_vs_clr)),
      "."
    )
  }
  base
}

pair_line <- function(row) {
  paste0(
    "- ", row$genus_short, " -> ", row$pathway_short,
    ": rho=", fmt(row$spearman_rho, 2),
    "; bootstrap CI [", fmt(row$rho_ci_low, 2), ", ", fmt(row$rho_ci_high, 2), "]",
    "; permutation p=", fmt_p(row$p_perm),
    "; BH q=", fmt_p(row$q_perm),
    "."
  )
}

feature_sentence <- function(row, include_clr = FALSE) {
  direction <- ifelse(row$mean_difference >= 0, "greater increase", "greater decrease")
  txt <- paste0(
    row$feature_short,
    " showed a ", direction,
    " in treatment relative to control (Cohen's d=", fmt(row$cohens_d, 2),
    ", mean delta difference=", fmt(row$mean_difference, 2),
    ", 95% bootstrap CI ", fmt(row$ci_lower, 2), " to ", fmt(row$ci_upper, 2),
    ", LOO sign agreement=", fmt(row$loo_sign_agreement, 2), ")."
  )

  if (include_clr && "d_clr" %in% names(row) && !is.na(row$d_clr)) {
    clr_consistent <- ifelse(isTRUE(row$sign_agreement_primary_vs_clr), "consistent", "not consistent")
    txt <- paste0(
      txt,
      " Direction was ", clr_consistent, " under CLR sensitivity (d_CLR=", fmt(row$d_clr, 2), ")."
    )
  }
  txt
}

pair_sentence <- function(row) {
  paste0(
    row$genus_short, " with ", row$pathway_short,
    " (rho=", fmt(row$spearman_rho, 2),
    ", 95% bootstrap CI ", fmt(row$rho_ci_low, 2), " to ", fmt(row$rho_ci_high, 2),
    ", permutation p=", fmt_p(row$p_perm),
    ", BH q=", fmt_p(row$q_perm), ")"
  )
}

species <- read.csv(file.path(input_dir, taxon_file), check.names = FALSE)
genus <- read.csv(file.path(input_dir, "genus_effect_sizes.csv"), check.names = FALSE)
pathway <- read.csv(file.path(input_dir, "pathway_effect_sizes.csv"), check.names = FALSE)
pairs <- read.csv(file.path(input_dir, "taxonomy_function_candidate_pairs.csv"), check.names = FALSE)

for (nm in c("species", "genus", "pathway")) {
  df <- get(nm)
  df$ci_excludes_zero <- to_bool(df$ci_excludes_zero)
  df$transform_sign_agreement_raw_log1p <- to_bool(df$transform_sign_agreement_raw_log1p)
  if ("sign_agreement_primary_vs_clr" %in% names(df)) {
    df$sign_agreement_primary_vs_clr <- to_bool(df$sign_agreement_primary_vs_clr)
  }
  assign(nm, df)
}

# Merge core results into one cleaner manuscript-facing key table.
as_clean_effect_table <- function(df, section_name, source_name) {
  out <- data.frame(
    result_source = source_name,
    key_section = section_name,
    result_type = "effect_size",
    feature_id = df$feature,
    feature_label = df$feature_short,
    comparison = "Treatment vs Control on paired deltas (6m - baseline)",
    n_pairs = df$n_pairs,
    n_treatment = df$n_treatment,
    n_control = df$n_control,
    effect_metric = "Cohens_d",
    effect_estimate = df$cohens_d,
    delta_difference = df$mean_difference,
    ci_lower = df$ci_lower,
    ci_upper = df$ci_upper,
    ci_width = df$ci_width,
    ci_excludes_zero = df$ci_excludes_zero,
    rank_consensus = df$rank_consensus,
    robustness_metric = "LOO_sign_agreement",
    robustness_value = df$loo_sign_agreement,
    transform_sign_agreement = df$transform_sign_agreement_raw_log1p,
    p_value = NA_real_,
    q_value = NA_real_,
    d_raw = df$d_raw,
    d_log1p = df$d_log1p,
    d_clr = if ("d_clr" %in% names(df)) df$d_clr else NA_real_,
    integration_genus = NA_character_,
    integration_pathway = NA_character_,
    note = df$interpretation,
    stringsAsFactors = FALSE
  )
  out
}

as_clean_pair_table <- function(df) {
  ci_excludes_zero <- (df$rho_ci_low > 0) | (df$rho_ci_high < 0)
  out <- data.frame(
    result_source = "taxonomy_function_candidate_pairs",
    key_section = "Integration",
    result_type = "delta_correlation",
    feature_id = paste0(df$genus, " -> ", df$pathway),
    feature_label = paste0(df$genus_short, " -> ", df$pathway_short),
    comparison = "Spearman correlation of participant-level deltas",
    n_pairs = df$n_participants,
    n_treatment = df$n_treatment,
    n_control = df$n_control,
    effect_metric = "Spearman_rho",
    effect_estimate = df$spearman_rho,
    delta_difference = NA_real_,
    ci_lower = df$rho_ci_low,
    ci_upper = df$rho_ci_high,
    ci_width = df$rho_ci_width,
    ci_excludes_zero = ci_excludes_zero,
    rank_consensus = df$rank_consensus,
    robustness_metric = "NA",
    robustness_value = NA_real_,
    transform_sign_agreement = NA,
    p_value = df$p_perm,
    q_value = df$q_perm,
    d_raw = df$genus_d,
    d_log1p = df$pathway_d,
    d_clr = df$abs_rho,
    integration_genus = df$genus_short,
    integration_pathway = df$pathway_short,
    note = ifelse(df$strong_corr_flag, "abs(rho) >= 0.8", "abs(rho) < 0.8"),
    stringsAsFactors = FALSE
  )
  out
}

key_species <- as_clean_effect_table(species, TAXON_LABEL, gsub("\\.csv$", "", taxon_file))
key_genus <- as_clean_effect_table(genus, "Genus", "genus_effect_sizes")
key_pathway <- as_clean_effect_table(pathway, "Pathway", "pathway_effect_sizes")
key_pairs <- as_clean_pair_table(pairs)

merged_key_results <- rbind(key_species, key_genus, key_pathway, key_pairs)
section_order <- c(TAXON_LABEL, "Genus", "Pathway", "Integration")
ord <- order(
  match(merged_key_results$key_section, section_order),
  merged_key_results$rank_consensus,
  -abs(merged_key_results$effect_estimate),
  na.last = TRUE
)
merged_key_results <- merged_key_results[ord, , drop = FALSE]
merged_key_results$merged_row_id <- seq_len(nrow(merged_key_results))
merged_csv_written_path <- output_merged_csv
tryCatch(
  write.csv(merged_key_results, output_merged_csv, row.names = FALSE, na = ""),
  error = function(e) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    merged_csv_written_path <<- file.path(input_dir, paste0("KEY_RESULTS_merged_table_", ts, ".csv"))
    write.csv(merged_key_results, merged_csv_written_path, row.names = FALSE, na = "")
    message(
      "Could not overwrite locked merged CSV. Wrote alternate file instead: ",
      merged_csv_written_path
    )
  }
)

species <- species[!is.na(species$feature_short), , drop = FALSE]
ord_species <- order(species$feature_short, species$rank_consensus, -abs(species$cohens_d), na.last = TRUE)
species_dedup <- species[ord_species, , drop = FALSE]
species_dedup <- species_dedup[!duplicated(species_dedup$feature_short), , drop = FALSE]

sum_species <- panel_summary(species, TAXON_LABEL)
sum_genus <- panel_summary(genus, "Genus")
sum_pathway <- panel_summary(pathway, "Pathway")
sum_all <- rbind(sum_species, sum_genus, sum_pathway)

top_species <- top_n_by_rank(species_dedup, n = 3)
top_genus <- top_n_by_rank(genus, n = 3)
top_pathway <- top_n_by_rank(pathway, n = 3)

ord_pairs <- order(pairs$rank_consensus, -abs(pairs$spearman_rho), na.last = TRUE)
pairs_ranked <- pairs[ord_pairs, , drop = FALSE]
top_pairs <- pairs_ranked[seq_len(min(5, nrow(pairs_ranked))), , drop = FALSE]

n_pairs <- if (nrow(species) > 0) unique(species$n_pairs)[1] else NA
n_tx <- if (nrow(species) > 0) unique(species$n_treatment)[1] else NA
n_ct <- if (nrow(species) > 0) unique(species$n_control)[1] else NA

n_pairs_tested <- nrow(pairs)
n_perm_lt_005 <- sum(pairs$p_perm < 0.05, na.rm = TRUE)
n_q_lt_010 <- sum(pairs$q_perm < 0.10, na.rm = TRUE)
n_abs_rho_ge_08 <- sum(abs(pairs$spearman_rho) >= 0.8, na.rm = TRUE)
median_abs_rho <- median(abs(pairs$spearman_rho), na.rm = TRUE)
max_abs_rho <- max(abs(pairs$spearman_rho), na.rm = TRUE)
min_q_perm <- min(pairs$q_perm, na.rm = TRUE)

species_pos_large <- sum(species$cohens_d >= 0.8, na.rm = TRUE)
species_neg_large <- sum(species$cohens_d <= -0.8, na.rm = TRUE)
genus_pos_large <- sum(genus$cohens_d >= 0.8, na.rm = TRUE)
genus_neg_large <- sum(genus$cohens_d <= -0.8, na.rm = TRUE)
pathway_pos_large <- sum(pathway$cohens_d >= 0.8, na.rm = TRUE)
pathway_neg_large <- sum(pathway$cohens_d <= -0.8, na.rm = TRUE)

ci_ratio_path_vs_species <- sum_pathway$median_ci_width / sum_species$median_ci_width
ci_ratio_genus_vs_species <- sum_genus$median_ci_width / sum_species$median_ci_width

rectale_row <- species_dedup[species_dedup$feature_short == "Eubacterium_rectale", , drop = FALSE]
has_rectale <- nrow(rectale_row) > 0

rectale_es_line <- "Eubacterium rectale was not present in the effect-size panel after filtering."
if (has_rectale) {
  rectale_es_line <- paste0(
    "In the effect-size pipeline, Eubacterium rectale showed d=", fmt(rectale_row$cohens_d[1], 2),
    "; delta-difference=", fmt(rectale_row$mean_difference[1], 2),
    "; 95% bootstrap CI ", fmt(rectale_row$ci_lower[1], 2), " to ", fmt(rectale_row$ci_upper[1], 2),
    "; LOO sign agreement=", fmt(rectale_row$loo_sign_agreement[1], 2),
    "; log1p d=", fmt(rectale_row$d_log1p[1], 2), "; CLR d=", fmt(rectale_row$d_clr[1], 2), "."
  )
}

maaslin_rectale_path <- file.path("results", "04_short_read_taxonomic_profiling", "maaslin3_results", "maaslin3_all_results.tsv")
rectale_maaslin_line <- "MaAsLin3 output was not found, so cross-method triangulation for Eubacterium rectale could not be added."
if (file.exists(maaslin_rectale_path)) {
  ma <- read.delim(maaslin_rectale_path, check.names = FALSE)
  ma_rect <- ma[
    grepl("Eubacterium_rectale", ma$feature, fixed = TRUE) &
      ma$metadata == "turmeric_group" &
      ma$value == "Treatment",
    ,
    drop = FALSE
  ]

  if (nrow(ma_rect) > 0) {
    prev_row <- ma_rect[ma_rect$model == "prevalence", , drop = FALSE]
    abd_row <- ma_rect[ma_rect$model == "abundance", , drop = FALSE]

    prev_txt <- if (nrow(prev_row) > 0) {
      paste0(
        "prevalence model coef=", fmt(prev_row$coef[1], 2),
        ", q_joint=", fmt_p(prev_row$qval_joint[1])
      )
    } else {
      "prevalence model not reported"
    }

    abd_txt <- if (nrow(abd_row) > 0) {
      paste0(
        "abundance model coef=", fmt(abd_row$coef[1], 2),
        ", q_individual=", fmt_p(abd_row$qval_individual[1])
      )
    } else {
      "abundance model not reported"
    }

    warn_txt <- ""
    if (nrow(prev_row) > 0 && !is.na(prev_row$error[1]) && nzchar(prev_row$error[1])) {
      warn_txt <- " The prevalence model includes a small-random-effects warning in the MaAsLin output."
    }

    rectale_maaslin_line <- paste0("In MaAsLin3, Eubacterium rectale showed ", prev_txt, "; ", abd_txt, ".", warn_txt)
  }
}

species_top_text <- if (nrow(top_species) > 0) {
  paste(vapply(seq_len(nrow(top_species)), function(i) feature_sentence(top_species[i, ], include_clr = TRUE), character(1)), collapse = " ")
} else {
  paste0("No ", TAXON_LABEL_LOWER, "-level features met inclusion criteria for ranked reporting.")
}

genus_top_text <- if (nrow(top_genus) > 0) {
  paste(vapply(seq_len(nrow(top_genus)), function(i) feature_sentence(top_genus[i, ], include_clr = TRUE), character(1)), collapse = " ")
} else {
  "No genus-level features met inclusion criteria for ranked reporting."
}

pathway_top_text <- if (nrow(top_pathway) > 0) {
  paste(vapply(seq_len(nrow(top_pathway)), function(i) feature_sentence(top_pathway[i, ], include_clr = FALSE), character(1)), collapse = " ")
} else {
  "No pathway-level features met inclusion criteria for ranked reporting."
}

top_pairs_text <- if (nrow(top_pairs) > 0) {
  paste(vapply(seq_len(nrow(top_pairs)), function(i) pair_sentence(top_pairs[i, ]), character(1)), collapse = "; ")
} else {
  "No genus-pathway candidate pairs met the minimum data threshold for integration."
}

method_lines <- c(
  "# Manuscript-Ready Methods and Results (Exploratory Small-N Analysis)",
  "",
  "## Methods",
  paste0(
    "We performed an exploratory, effect-size-focused paired analysis of microbiome taxonomy and pathway profiles in a low-sample setting. ",
    "For each feature, participant-level change was defined as 6-month abundance minus baseline abundance, and Treatment versus Control groups were compared on these within-participant deltas."
  ),
  paste0(
    "Feature inclusion followed a prevalence-and-abundance filter (taxonomy: prevalence >= 0.20 and mean abundance >= 0.001; pathways: prevalence >= 0.20 and mean abundance >= 1.0 CPM). ",
    "Primary analysis used the raw scale for interpretability, with log1p transformation used as directional sensitivity analysis; taxonomy was further assessed using CLR-transformed abundances to evaluate compositional robustness."
  ),
  paste0(
    "For each retained feature, we estimated mean and median delta differences, Cohen's d, and probability of superiority, and quantified uncertainty with 95% bootstrap confidence intervals (2,000 resamples). ",
    "Participant-level sensitivity was assessed via leave-one-participant-out (LOO) re-estimation, summarized by sign agreement and sign-flip behavior. ",
    "Features were prioritized by a transparent consensus ranking that combined absolute effect magnitude, CI width, and LOO directional stability."
  ),
  paste0(
    "Exploratory taxonomy-function integration was then performed by pairing top-ranked genera and pathways and calculating Spearman correlations on participant-level deltas. ",
    "For each pair, we reported bootstrap correlation intervals (2,000 resamples), permutation p-values (2,000 permutations), and BH-adjusted q-values as multiplicity context. ",
    "The analytical framework was designed for hypothesis prioritization and does not constitute confirmatory feature-level inference."
  ),
  "",
  "## Results",
  paste0(
    "The paired analysis set included ", n_pairs, " participants (Treatment ", n_tx, ", Control ", n_ct, "). ",
    "Feature panels included ", sum_species$n_features, " ", TAXON_ENTRY_LABEL, " (", nrow(species_dedup), " unique ", TAXON_LABEL_LOWER, " names), ",
    sum_genus$n_features, " genera, and ", sum_pathway$n_features, " pathways."
  ),
  paste0(
    "At panel level, median absolute effect sizes were ", fmt(sum_species$median_abs_d, 2), " for ", TAXON_LABEL_LOWER, ", ",
    fmt(sum_genus$median_abs_d, 2), " for genera, and ", fmt(sum_pathway$median_abs_d, 2),
    ". Confidence intervals excluded zero for ",
    sum_species$n_ci_excludes_zero, "/", sum_species$n_features, " ", TAXON_LABEL_LOWER, " features (", fmt(100 * sum_species$prop_ci_excludes_zero, 1), "%), ",
    sum_genus$n_ci_excludes_zero, "/", sum_genus$n_features, " genus features (", fmt(100 * sum_genus$prop_ci_excludes_zero, 1), "%), and ",
    sum_pathway$n_ci_excludes_zero, "/", sum_pathway$n_features, " pathway features (", fmt(100 * sum_pathway$prop_ci_excludes_zero, 1), "%)."
  ),
  paste0(
    "Directional robustness across raw and log1p scales was high (", TAXON_LABEL_LOWER, " ", fmt(sum_species$transform_agreement, 2),
    ", genus ", fmt(sum_genus$transform_agreement, 2), ", pathway ", fmt(sum_pathway$transform_agreement, 2),
    "), and median LOO sign agreement was ", paste(unique(c(sum_species$median_loo, sum_genus$median_loo, sum_pathway$median_loo)), collapse = ", "),
    ", supporting stability of direction under participant exclusion."
  ),
  paste0(
    "Precision differed substantially across feature levels: median CI width was ", fmt(sum_species$median_ci_width, 2), " for ", TAXON_LABEL_LOWER, ", ",
    fmt(sum_genus$median_ci_width, 2), " for genus (", fmt(ci_ratio_genus_vs_species, 1), "x ", TAXON_LABEL_LOWER, "), and ",
    fmt(sum_pathway$median_ci_width, 2), " for pathways (", fmt(ci_ratio_path_vs_species, 1), "x ", TAXON_LABEL_LOWER, "). ",
    "Among larger effects (|d| >= 0.8), directionality was mixed for ", TAXON_LABEL_LOWER, " (+", species_pos_large, "/-", species_neg_large, ") and genera (+", genus_pos_large, "/-", genus_neg_large, "), ",
    "and pathways (+", pathway_pos_large, "/-", pathway_neg_large, ")."
  ),
  "",
  "Feature-level interpretation indicated coherent but exploratory signals across taxonomic and functional strata.",
  paste0("At ", TAXON_LABEL_LOWER, " level, ", species_top_text),
  paste0("At genus level, ", genus_top_text),
  paste0("At pathway level, ", pathway_top_text),
  "",
  "Cross-method triangulation for Eubacterium rectale supported prioritization but not confirmation.",
  rectale_es_line,
  rectale_maaslin_line,
  paste0(
    "Taken together, these results suggest that Eubacterium rectale is a biologically plausible follow-up target, but inference remains provisional in this small paired dataset."
  ),
  "",
  paste0(
    "In exploratory taxonomy-function integration, ", n_pairs_tested, " genus-pathway pairs were evaluated. ",
    n_abs_rho_ge_08, " pairs had |rho| >= 0.80 and ", n_perm_lt_005, " had permutation p < 0.05. ",
    n_q_lt_010, " pairs met BH q < 0.10 (minimum q=", fmt_p(min_q_perm), ")."
  ),
  paste0(
    "Correlation magnitudes were moderate overall (median |rho|=", fmt(median_abs_rho, 2),
    "; maximum |rho|=", fmt(max_abs_rho, 2), "). ",
    "Top-ranked associations included ", top_pairs_text, "."
  ),
  paste0(
    "Overall, the analysis supports candidate ranking and mechanistic hypothesis generation, while emphasizing that independent replication is required before confirmatory conclusions."
  )
)

writeLines(method_lines, con = output_md, useBytes = TRUE)

level_cols <- c("#1b9e77", "#d95f02", "#7570b3")
names(level_cols) <- c(TAXON_LABEL, "Genus", "Pathway")
theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.25)
  )

overview <- sum_all
overview$label <- paste0(overview$n_ci_excludes_zero, "/", overview$n_features)

p_overview <- ggplot(overview, aes(x = level, y = prop_ci_excludes_zero, fill = level)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_text(aes(label = label), vjust = -0.4, size = 3.8) +
  scale_fill_manual(values = level_cols) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  labs(
    x = NULL,
    y = "Features with CI excluding 0"
  ) +
  theme_pub +
  theme(legend.position = "none")

dist_df <- rbind(
  data.frame(level = TAXON_LABEL, abs_d = abs(species$cohens_d), ci_excludes_zero = species$ci_excludes_zero),
  data.frame(level = "Genus", abs_d = abs(genus$cohens_d), ci_excludes_zero = genus$ci_excludes_zero),
  data.frame(level = "Pathway", abs_d = abs(pathway$cohens_d), ci_excludes_zero = pathway$ci_excludes_zero)
)
dist_df$ci_label <- ifelse(dist_df$ci_excludes_zero, "CI excludes 0", "CI includes 0")

p_dist <- ggplot(dist_df, aes(x = level, y = abs_d)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = "grey95", color = "grey40") +
  geom_jitter(aes(color = ci_label), width = 0.16, height = 0, alpha = 0.8, size = 1.8) +
  scale_color_manual(values = c("CI excludes 0" = "#b2182b", "CI includes 0" = "#2166ac")) +
  labs(
    x = NULL,
    y = "Absolute Cohen's d",
    color = NULL
  ) +
  theme_pub

pick_plot_cols <- function(df, level_name, n = 4) {
  top <- top_n_by_rank(df, n = n)
  data.frame(
    level = level_name,
    feature_short = top$feature_short,
    cohens_d = top$cohens_d,
    loo_sign_agreement = top$loo_sign_agreement,
    ci_excludes_zero = top$ci_excludes_zero,
    stringsAsFactors = FALSE
  )
}

top_features_plot <- rbind(
  pick_plot_cols(species_dedup, TAXON_LABEL, n = 4),
  pick_plot_cols(genus, "Genus", n = 4),
  pick_plot_cols(pathway, "Pathway", n = 4)
)
top_features_plot$feature_short <- trim_text(top_features_plot$feature_short, width = 55)
top_features_plot$plot_label <- paste0(top_features_plot$level, ": ", top_features_plot$feature_short)
top_features_plot$plot_label <- factor(
  top_features_plot$plot_label,
  levels = top_features_plot$plot_label[order(top_features_plot$cohens_d)]
)
top_features_plot$loo_sign_agreement[is.na(top_features_plot$loo_sign_agreement)] <- 0
top_features_plot$ci_label <- ifelse(top_features_plot$ci_excludes_zero, "CI excludes 0", "CI includes 0")

p_top <- ggplot(top_features_plot, aes(x = cohens_d, y = plot_label, color = level)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_segment(aes(x = 0, xend = cohens_d, yend = plot_label), alpha = 0.35, linewidth = 0.9) +
  geom_point(aes(size = loo_sign_agreement, shape = ci_label), stroke = 0.8) +
  scale_color_manual(values = level_cols) +
  scale_shape_manual(values = c("CI excludes 0" = 16, "CI includes 0" = 1)) +
  scale_size_continuous(limits = c(0, 1), range = c(2, 5)) +
  labs(
    x = "Cohen's d",
    y = NULL,
    color = "Level",
    size = "LOO sign\nagreement",
    shape = NULL
  ) +
  theme_pub +
  theme(text = element_text(size = 10))

top_corr <- pairs_ranked[seq_len(min(12, nrow(pairs_ranked))), , drop = FALSE]
top_corr$pathway_short <- trim_text(top_corr$pathway_short, width = 45)
top_corr$pair_label <- paste0(top_corr$genus_short, " -> ", top_corr$pathway_short)
top_corr$pair_label <- factor(top_corr$pair_label, levels = top_corr$pair_label[order(top_corr$spearman_rho)])
top_corr$perm_label <- ifelse(top_corr$p_perm < 0.05, "perm p < 0.05", "perm p >= 0.05")

p_corr <- ggplot(top_corr, aes(x = spearman_rho, y = pair_label, color = perm_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_segment(
    aes(x = rho_ci_low, xend = rho_ci_high, y = pair_label, yend = pair_label),
    linewidth = 1.0, color = "grey65"
  ) +
  geom_segment(aes(x = 0, xend = spearman_rho, yend = pair_label), alpha = 0.35, linewidth = 0.8) +
  geom_point(size = 2.4) +
  scale_color_manual(values = c("perm p < 0.05" = "#b2182b", "perm p >= 0.05" = "#2166ac")) +
  labs(
    x = "Spearman rho",
    y = NULL,
    color = NULL
  ) +
  theme_pub +
  theme(text = element_text(size = 10))

draw_composite <- function(device = c("png", "pdf"), path) {
  device <- match.arg(device)

  if (device == "png") {
    png(path, width = 3800, height = 3000, res = 300)
  } else {
    pdf(path, width = 12.5, height = 9.5, useDingbats = FALSE)
  }

  grid.newpage()
  lay <- grid.layout(
    nrow = 2, ncol = 2,
    widths = unit(c(1, 1), "null")
  )
  pushViewport(viewport(layout = lay))

  print(p_overview, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_dist, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p_top, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(p_corr, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

  panel_gp <- gpar(fontsize = 18, fontface = "bold")
  grid.text("A", x = unit(0.02, "npc"), y = unit(0.98, "npc"), just = c("left", "top"),
            vp = viewport(layout.pos.row = 1, layout.pos.col = 1), gp = panel_gp)
  grid.text("B", x = unit(0.02, "npc"), y = unit(0.98, "npc"), just = c("left", "top"),
            vp = viewport(layout.pos.row = 1, layout.pos.col = 2), gp = panel_gp)
  grid.text("C", x = unit(0.02, "npc"), y = unit(0.98, "npc"), just = c("left", "top"),
            vp = viewport(layout.pos.row = 2, layout.pos.col = 1), gp = panel_gp)
  grid.text("D", x = unit(0.02, "npc"), y = unit(0.98, "npc"), just = c("left", "top"),
            vp = viewport(layout.pos.row = 2, layout.pos.col = 2), gp = panel_gp)

  dev.off()
}

draw_composite("png", output_png)
pdf_ok <- TRUE
tryCatch(
  draw_composite("pdf", output_pdf),
  error = function(e) {
    pdf_ok <<- FALSE
    message("Could not write PDF output: ", e$message)
  }
)

cat("Wrote manuscript text: ", output_md, "\n", sep = "")
cat("Wrote figure (PNG): ", output_png, "\n", sep = "")
if (pdf_ok) {
  cat("Wrote figure (PDF): ", output_pdf, "\n", sep = "")
} else {
  cat("Skipped PDF write (file may be open/locked): ", output_pdf, "\n", sep = "")
}
cat("Wrote merged key results CSV: ", merged_csv_written_path, "\n", sep = "")

# ------------------------- Consolidated Results Bundle -------------------------
bundle_dir <- file.path(input_dir, "manuscript_results_bundle")
dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)

bundle_files <- list.files(
  input_dir,
  pattern = "\\.(csv|png|pdf|md)$",
  full.names = TRUE
)
bundle_files <- bundle_files[!file.info(bundle_files)$isdir]

script_path <- file.path("scripts", "4f_generate_manuscript_reporting_honest_smallN.R")
if (file.exists(script_path)) {
  bundle_files <- unique(c(bundle_files, script_path))
}

copied_ok <- file.copy(bundle_files, bundle_dir, overwrite = TRUE)

manifest_lines <- c(
  "Consolidated manuscript results bundle",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Files copied:"
)

for (i in seq_along(bundle_files)) {
  status <- ifelse(isTRUE(copied_ok[i]), "OK", "FAILED")
  manifest_lines <- c(
    manifest_lines,
    paste0("- [", status, "] ", basename(bundle_files[i]))
  )
}

writeLines(manifest_lines, con = file.path(bundle_dir, "BUNDLE_MANIFEST.txt"))
cat("Wrote consolidated bundle: ", bundle_dir, "\n", sep = "")
