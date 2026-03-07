################################################################################
# Small-Sample Paired Taxonomy + Function Analysis (Exploratory)
# Objective: quantify treatment-associated change from baseline to 6 months
# using participant-level deltas (6_months - baseline) in Treatment vs Control.
# Inference is effect-size focused (Cohen's d, probability of superiority,
# bootstrap CI) with robustness diagnostics (leave-one-participant-out,
# transform agreement, and CLR compositional sensitivity).
# Outputs are intended for hypothesis prioritization in low-N settings,
# not confirmatory feature-level significance claims.
################################################################################

set.seed(123)

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
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(glue)
  library(data.table)
  library(gh)
  library(base64enc)
})

owner <- "LiamHWalsh"
repo <- "NTU_Teagasc_Turmeric_study"
branch <- "main"
token <- Sys.getenv("GITHUB_PAT")
if (token == "") token <- Sys.getenv("GITHUB_TOKEN")

# ------------------------- Configuration -------------------------
SUPPORTED_TRANSFORMS <- c("raw", "log1p")
PRIMARY_TRANSFORM <- "raw"
SECONDARY_TRANSFORM <- "log1p"
RUN_CLR_SENSITIVITY <- TRUE
TAXON_RESOLUTION <- "strain" # "strain" (t__) or "species_tagged" (contains s__)

MIN_PREVALENCE_TAXA <- 0.20
MIN_PREVALENCE_PATHWAY <- 0.20
MIN_MEAN_ABUND_TAXA <- 0.001
MIN_MEAN_ABUND_PATHWAY <- 1.0

N_BOOTSTRAP <- 2000
N_BOOT_RHO <- 2000
N_PERM_RHO <- 2000
CI_LEVEL <- 0.95
MIN_GROUP_N <- 3

MAX_ABUNDANT_GENUS <- 15
MAX_SPECIES_PANEL <- 60
MAX_PATHWAY_PANEL <- 40
TOP_FOR_PLOTS <- 20
TOP_FOR_INTEGRATION_GENUS <- 8
TOP_FOR_INTEGRATION_PATHWAY <- 20
USE_ALL_RETAINED_FOR_INFERENCE <- TRUE
ENABLE_GITHUB_UPLOAD <- FALSE
GENUS_PANEL_MODE <- "pilot_2018_curcumin" # "pilot_2018_curcumin" or "pilot_plus_abundant"
FILTER_PATHWAYS_BY_KEYWORD <- FALSE
PATHWAY_KEYWORDS <- c("CARBOHYDRATE", "LIPID", "AMINO", "BIOSYN", "DEGRAD", "FERMENT", "GLYCOL", "TCA", "COA")

if (!TAXON_RESOLUTION %in% c("strain", "species_tagged")) {
  stop("TAXON_RESOLUTION must be 'strain' or 'species_tagged'")
}

TAXON_FEATURE_PATTERN <- if (TAXON_RESOLUTION == "strain") "\\|t__" else "\\|s__"
TAXON_LEVEL_LABEL <- if (TAXON_RESOLUTION == "strain") "Strain" else "Species"
TAXON_FILE_PREFIX <- if (TAXON_RESOLUTION == "strain") "strain" else "species"

output_subdir <- if (TAXON_RESOLUTION == "strain") "honest_small_n_strain" else "honest_small_n"
output_dir <- file.path("results", "Processed", "05_statistical_analysis", output_subdir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!PRIMARY_TRANSFORM %in% SUPPORTED_TRANSFORMS) stop("Unsupported PRIMARY_TRANSFORM")
if (!SECONDARY_TRANSFORM %in% SUPPORTED_TRANSFORMS) stop("Unsupported SECONDARY_TRANSFORM")

# ------------------------- IO Helpers -------------------------
transform_value <- function(x, mode = PRIMARY_TRANSFORM) {
  if (mode == "log1p") return(log1p(pmax(x, 0)))
  x
}

read_text_table <- function(txt, ext_hint = "tsv") {
  ext <- tolower(ext_hint)
  if (ext == "csv") return(readr::read_csv(I(txt), show_col_types = FALSE, name_repair = "minimal"))

  first <- strsplit(txt, "\n", fixed = TRUE)[[1]][1]
  if (!is.na(first) && startsWith(first, "#mpa_")) {
    return(readr::read_tsv(I(txt), show_col_types = FALSE, comment = "#", name_repair = "minimal"))
  }
  readr::read_tsv(I(txt), show_col_types = FALSE, comment = "", name_repair = "minimal")
}

load_github_file <- function(repo_path, token, owner, repo, branch = "main") {
  if (!nzchar(token)) stop("No GitHub token found.")
  info <- gh::gh(
    "GET /repos/:owner/:repo/contents/:path",
    owner = owner, repo = repo, path = repo_path, ref = branch, .token = token
  )
  txt <- rawToChar(base64enc::base64decode(gsub("\n", "", info$content)))
  read_text_table(txt, ext_hint = tolower(tools::file_ext(repo_path)))
}

load_local_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(readr::read_csv(path, show_col_types = FALSE, name_repair = "minimal"))

  first <- readLines(path, n = 1, warn = FALSE)
  if (length(first) == 0) stop("Empty file: ", path)
  if (startsWith(first, "#mpa_")) {
    return(readr::read_tsv(path, show_col_types = FALSE, comment = "#", name_repair = "minimal"))
  }
  readr::read_tsv(path, show_col_types = FALSE, comment = "", name_repair = "minimal")
}

load_first_available <- function(github_paths, local_paths, label) {
  if (nzchar(token)) {
    for (p in github_paths) {
      out <- tryCatch(load_github_file(p, token, owner, repo, branch), error = function(e) NULL)
      if (!is.null(out)) {
        cat("Loaded", label, "from GitHub:", p, "\n")
        return(out)
      }
    }
  }

  for (p in local_paths) {
    if (file.exists(p)) {
      cat("Loaded", label, "from local:", p, "\n")
      return(load_local_table(p))
    }
  }
  stop("Could not load ", label)
}

upload_to_github <- function(local_path, github_path, message) {
  if (!nzchar(token)) return(invisible(NULL))
  encoded <- base64enc::base64encode(local_path)
  info <- tryCatch(
    gh::gh("GET /repos/:owner/:repo/contents/:path", owner = owner, repo = repo, path = github_path, ref = branch, .token = token),
    error = function(e) NULL
  )
  sha <- if (!is.null(info)) info$sha else NULL
  gh::gh(
    "PUT /repos/:owner/:repo/contents/:path",
    owner = owner, repo = repo, path = github_path,
    message = message, content = encoded, branch = branch, sha = sha, .token = token
  )
  cat("Uploaded:", github_path, "\n")
}

to_github_path <- function(path) gsub("\\\\", "/", path)

# ------------------------- Metadata + Parsing -------------------------
normalize_sample_id <- function(x) {
  x %>%
    str_replace("_metaphlan$", "") %>%
    str_replace("_L001_concat_Abundance-CPM$", "")
}

extract_rank <- function(feature_string, prefix) {
  parts <- str_split(feature_string, "\\|")[[1]]
  hit <- parts[str_detect(parts, paste0("^", prefix, "__"))]
  if (length(hit) == 0) return(NA_character_)
  str_remove(hit[1], paste0("^", prefix, "__"))
}

extract_short_name <- function(feature_string) {
  if (is.na(feature_string) || feature_string == "") return(NA_character_)
  if (TAXON_RESOLUTION == "strain" && str_detect(feature_string, "\\|t__")) return(extract_rank(feature_string, "t"))
  if (str_detect(feature_string, "\\|s__")) return(extract_rank(feature_string, "s"))
  if (str_detect(feature_string, "\\|t__")) return(extract_rank(feature_string, "t"))
  if (str_detect(feature_string, "\\|")) return(tail(str_split(feature_string, "\\|")[[1]], 1))
  feature_string
}

prepare_metadata <- function(sample_ids) {
  split_ids <- str_split_fixed(sample_ids, "_", 3)
  md <- data.frame(
    sample_id = sample_ids,
    block_code = split_ids[, 2],
    block_number = as.integer(gsub("[^0-9]", "", split_ids[, 2])),
    block_letter = str_to_upper(gsub("[0-9]", "", split_ids[, 2])),
    stringsAsFactors = FALSE
  )

  turmeric_status <- data.frame(
    ID = c(3, 4, 5, 6, 7, 9, 11, 14, 15, 16, 18, 19, 20, 21, 22, 23),
    status = c(
      "Didn't take turmeric", "Turmeric 6 months", "Turmeric 6 months",
      "Didn't take turmeric", "Turmeric 6 months", "Turmeric 6 months",
      "Turmeric 3 months", "Baseline sample only", "Turmeric 6 months",
      "Turmeric 6 months", "Didn't take turmeric", "Baseline sample only",
      "Turmeric 3 months", "Turmeric 6 months", "Turmeric 6 months",
      "Baseline sample only"
    ),
    stringsAsFactors = FALSE
  )

  md <- merge(md, turmeric_status, by.x = "block_number", by.y = "ID", all.x = TRUE) %>%
    filter(!is.na(block_number), !is.na(status))

  md$timepoint <- case_when(
    md$block_letter == "A" ~ "Baseline",
    md$block_letter == "B" ~ "2_weeks",
    md$block_letter == "C" ~ "6_months",
    TRUE ~ NA_character_
  )

  md$turmeric_group <- "Control"
  md$turmeric_group[grepl("Turmeric 6 months", md$status) & md$block_letter == "C"] <- "Treatment"
  md$turmeric_group[md$status == "Turmeric 3 months" & md$block_letter == "C"] <- "Stopped"

  for (i in seq_len(nrow(md))) {
    if (!is.na(md$timepoint[i]) && md$timepoint[i] %in% c("Baseline", "2_weeks")) {
      pid <- md$block_number[i]
      has_tx <- any(md$block_number == pid & md$turmeric_group == "Treatment", na.rm = TRUE)
      if (has_tx) md$turmeric_group[i] <- "Treatment"
    }
  }

  md %>%
    filter(turmeric_group != "Stopped") %>%
    mutate(
      participant = as.character(block_number),
      timepoint = factor(timepoint, levels = c("Baseline", "2_weeks", "6_months")),
      turmeric_group = factor(turmeric_group, levels = c("Control", "Treatment"))
    )
}

coerce_numeric_samples <- function(df, sample_cols) {
  out <- df
  for (sc in sample_cols) out[[sc]] <- suppressWarnings(as.numeric(out[[sc]]))
  out
}

filter_feature_table <- function(df, sample_cols, min_prev, min_mean_abund) {
  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0

  prev <- rowMeans(mat > 0)
  mean_abund <- rowMeans(mat)
  keep <- prev >= min_prev & mean_abund >= min_mean_abund
  df[keep, , drop = FALSE]
}

# ------------------------- Effect + Sensitivity Functions -------------------------
get_change_table <- function(feature_values, metadata, transform_mode = PRIMARY_TRANSFORM) {
  long <- data.frame(
    sample_id = names(feature_values),
    abundance = as.numeric(feature_values),
    stringsAsFactors = FALSE
  )

  dat <- metadata %>%
    inner_join(long, by = "sample_id") %>%
    mutate(abundance = transform_value(abundance, mode = transform_mode))

  dat %>%
    filter(timepoint %in% c("Baseline", "6_months")) %>%
    select(participant, turmeric_group, timepoint, abundance) %>%
    pivot_wider(names_from = timepoint, values_from = abundance) %>%
    filter(!is.na(Baseline), !is.na(`6_months`)) %>%
    mutate(delta = `6_months` - Baseline)
}

compute_effect_no_boot <- function(feature_values, metadata, transform_mode = PRIMARY_TRANSFORM) {
  change <- get_change_table(feature_values, metadata, transform_mode = transform_mode)

  tx <- change %>% filter(turmeric_group == "Treatment") %>% pull(delta)
  ct <- change %>% filter(turmeric_group == "Control") %>% pull(delta)
  if (length(tx) < MIN_GROUP_N || length(ct) < MIN_GROUP_N) return(NULL)

  mean_diff <- mean(tx) - mean(ct)
  median_diff <- median(tx) - median(ct)
  pooled_sd <- sqrt((var(tx) + var(ct)) / 2)
  cohens_d <- ifelse(is.finite(pooled_sd) && pooled_sd > 0, mean_diff / pooled_sd, NA_real_)

  cmp <- outer(tx, ct, "-")
  prob_superiority <- mean(cmp > 0) + 0.5 * mean(cmp == 0)

  data.frame(
    n_pairs = nrow(change),
    n_treatment = length(tx),
    n_control = length(ct),
    mean_change_treatment = mean(tx),
    sd_change_treatment = sd(tx),
    mean_change_control = mean(ct),
    sd_change_control = sd(ct),
    mean_difference = mean_diff,
    median_difference = median_diff,
    cohens_d = cohens_d,
    prob_superiority = prob_superiority,
    stringsAsFactors = FALSE
  )
}

compute_effect_with_ci <- function(feature_values, metadata, transform_mode = PRIMARY_TRANSFORM,
                                   n_boot = N_BOOTSTRAP, ci_level = CI_LEVEL) {
  base <- compute_effect_no_boot(feature_values, metadata, transform_mode = transform_mode)
  if (is.null(base)) return(NULL)

  change <- get_change_table(feature_values, metadata, transform_mode = transform_mode)
  tx <- change %>% filter(turmeric_group == "Treatment") %>% pull(delta)
  ct <- change %>% filter(turmeric_group == "Control") %>% pull(delta)

  alpha <- (1 - ci_level) / 2
  boot_diffs <- replicate(n_boot, {
    tx_boot <- sample(tx, replace = TRUE)
    ct_boot <- sample(ct, replace = TRUE)
    mean(tx_boot) - mean(ct_boot)
  })
  ci <- quantile(boot_diffs, c(alpha, 1 - alpha), na.rm = TRUE)

  cbind(
    base,
    data.frame(
      ci_lower = as.numeric(ci[1]),
      ci_upper = as.numeric(ci[2]),
      ci_width = as.numeric(ci[2] - ci[1]),
      ci_excludes_zero = (ci[1] > 0) || (ci[2] < 0),
      stringsAsFactors = FALSE
    )
  )
}

assess_loo_sensitivity <- function(feature_values, metadata, transform_mode = PRIMARY_TRANSFORM,
                                   reference_mean_diff = NA_real_) {
  participants <- unique(metadata$participant)

  loo_stats <- lapply(participants, function(pid) {
    md_sub <- metadata %>% filter(participant != pid)
    sample_sub <- intersect(names(feature_values), md_sub$sample_id)
    if (length(sample_sub) < 6) return(NULL)
    compute_effect_no_boot(feature_values[sample_sub], md_sub, transform_mode = transform_mode)
  })

  loo_stats <- bind_rows(loo_stats)
  if (nrow(loo_stats) == 0) return(NULL)

  ref_sign <- sign(reference_mean_diff)
  loo_sign <- sign(loo_stats$mean_difference)

  sign_agreement <- NA_real_
  sign_flip_frac <- NA_real_
  if (is.finite(ref_sign) && ref_sign != 0) {
    sign_agreement <- mean(loo_sign == ref_sign, na.rm = TRUE)
    sign_flip_frac <- mean(loo_sign == (-ref_sign), na.rm = TRUE)
  }

  data.frame(
    loo_n_folds = nrow(loo_stats),
    loo_sign_agreement = sign_agreement,
    loo_sign_flip_frac = sign_flip_frac,
    loo_median_d = median(loo_stats$cohens_d, na.rm = TRUE),
    loo_range_d = diff(range(loo_stats$cohens_d, na.rm = TRUE)),
    loo_iqr_d = IQR(loo_stats$cohens_d, na.rm = TRUE),
    loo_median_mean_diff = median(loo_stats$mean_difference, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

compute_clr_table <- function(df, feature_col, sample_cols) {
  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0

  col_sum <- colSums(mat)
  for (j in seq_len(ncol(mat))) {
    if (col_sum[j] > 0) {
      mat[, j] <- mat[, j] / col_sum[j]
    } else {
      mat[, j] <- 0
    }
  }

  pos <- mat[mat > 0]
  pseudo <- if (length(pos) > 0) min(pos) / 2 else 1e-6
  clr <- sweep(log(mat + pseudo), 2, colMeans(log(mat + pseudo)), FUN = "-")

  out <- as.data.frame(clr)
  colnames(out) <- sample_cols
  out[[feature_col]] <- df[[feature_col]]
  out[, c(feature_col, sample_cols), drop = FALSE]
}

compute_panel_clr_effects <- function(clr_df, feature_col, feature_list, metadata) {
  rows <- lapply(feature_list, function(fid) {
    row <- clr_df[clr_df[[feature_col]] == fid, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)

    x <- unlist(row[1, setdiff(colnames(row), feature_col), drop = TRUE])
    names(x) <- setdiff(colnames(row), feature_col)
    eff <- compute_effect_no_boot(x, metadata, transform_mode = "raw")
    if (is.null(eff)) return(NULL)

    data.frame(feature = fid, d_clr = eff$cohens_d, mean_diff_clr = eff$mean_difference, stringsAsFactors = FALSE)
  })
  bind_rows(rows)
}

analyze_feature_panel <- function(abundance_df, feature_col, feature_list, metadata, level_name) {
  cat("\nAnalyzing", level_name, "panel:", length(feature_list), "features\n")

  rows <- lapply(feature_list, function(feat) {
    row <- abundance_df[abundance_df[[feature_col]] == feat, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)

    x <- unlist(row[1, setdiff(colnames(row), feature_col), drop = TRUE])
    names(x) <- setdiff(colnames(row), feature_col)

    main_eff <- compute_effect_with_ci(x, metadata, transform_mode = PRIMARY_TRANSFORM)
    if (is.null(main_eff)) return(NULL)

    raw_eff <- compute_effect_no_boot(x, metadata, transform_mode = "raw")
    log_eff <- compute_effect_no_boot(x, metadata, transform_mode = "log1p")

    loo <- assess_loo_sensitivity(
      feature_values = x,
      metadata = metadata,
      transform_mode = PRIMARY_TRANSFORM,
      reference_mean_diff = main_eff$mean_difference
    )

    cbind(
      data.frame(
        feature = feat,
        feature_short = extract_short_name(feat),
        level = level_name,
        transform_primary = PRIMARY_TRANSFORM,
        stringsAsFactors = FALSE
      ),
      main_eff,
      data.frame(
        d_raw = if (!is.null(raw_eff)) raw_eff$cohens_d else NA_real_,
        d_log1p = if (!is.null(log_eff)) log_eff$cohens_d else NA_real_,
        stringsAsFactors = FALSE
      ),
      loo
    )
  })

  out <- bind_rows(rows)
  if (nrow(out) == 0) return(out)

  out <- out %>%
    mutate(
      abs_d = abs(cohens_d),
      transform_sign_agreement_raw_log1p = case_when(
        is.na(d_raw) | is.na(d_log1p) ~ NA,
        sign(d_raw) == sign(d_log1p) ~ TRUE,
        TRUE ~ FALSE
      ),
      ci_width_rank = dense_rank(ci_width),
      abs_d_rank = dense_rank(desc(abs_d)),
      loo_rank = dense_rank(desc(ifelse(is.na(loo_sign_agreement), 0, loo_sign_agreement))),
      rank_consensus = abs_d_rank + ci_width_rank + loo_rank
    )

  ci_ref <- median(out$ci_width, na.rm = TRUE)
  out %>%
    mutate(
      effect_class = case_when(abs_d >= 0.8 ~ "Large", abs_d >= 0.5 ~ "Moderate", TRUE ~ "Small"),
      precision_class = ifelse(ci_width <= ci_ref, "Narrower CI", "Wider CI"),
      interpretation = paste(effect_class, "effect,", precision_class)
    ) %>%
    arrange(rank_consensus, abs_d_rank, ci_width_rank)
}

extract_feature_delta <- function(abundance_df, feature_col, feature_id, metadata,
                                  transform_mode = PRIMARY_TRANSFORM) {
  row <- abundance_df[abundance_df[[feature_col]] == feature_id, , drop = FALSE]
  if (nrow(row) == 0) return(NULL)

  x <- unlist(row[1, setdiff(colnames(row), feature_col), drop = TRUE])
  names(x) <- setdiff(colnames(row), feature_col)

  change <- get_change_table(x, metadata, transform_mode = transform_mode)
  if (nrow(change) == 0) return(NULL)
  out <- change %>% select(participant, turmeric_group, delta)
  out$feature <- feature_id
  out
}

# ------------------------- Correlation + Plot Helpers -------------------------
bootstrap_spearman_ci <- function(x, y, n_boot = N_BOOT_RHO, ci_level = CI_LEVEL) {
  if (length(x) != length(y) || length(x) < 6) return(c(NA_real_, NA_real_))
  n <- length(x)
  vals <- replicate(n_boot, {
    idx <- sample(seq_len(n), replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman", use = "complete.obs"))
  })
  vals <- vals[is.finite(vals)]
  if (length(vals) < 50) return(c(NA_real_, NA_real_))
  a <- (1 - ci_level) / 2
  as.numeric(quantile(vals, c(a, 1 - a), na.rm = TRUE))
}

perm_spearman_p <- function(x, y, n_perm = N_PERM_RHO) {
  if (length(x) != length(y) || length(x) < 6) return(NA_real_)
  obs <- suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
  if (!is.finite(obs)) return(NA_real_)

  perm_vals <- replicate(n_perm, {
    suppressWarnings(cor(x, sample(y), method = "spearman", use = "complete.obs"))
  })
  (sum(abs(perm_vals) >= abs(obs), na.rm = TRUE) + 1) / (n_perm + 1)
}

make_forest_plot <- function(results_df, title, top_n = TOP_FOR_PLOTS) {
  if (is.null(results_df) || nrow(results_df) == 0) return(NULL)

  plot_df <- results_df %>%
    filter(!is.na(cohens_d), !is.na(ci_lower), !is.na(ci_upper)) %>%
    arrange(rank_consensus, abs_d_rank) %>%
    slice_head(n = top_n) %>%
    mutate(
      feature_label = paste0(feature_short, " (n=", n_pairs, ")"),
      feature_label = reorder(feature_label, cohens_d)
    )

  ggplot(plot_df, aes(y = feature_label, x = cohens_d, color = interpretation)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0.18, linewidth = 0.8) +
    geom_point(size = 2.8, alpha = 0.9) +
    scale_color_manual(
      values = c(
        "Large effect, Narrower CI" = "#b22222",
        "Large effect, Wider CI" = "#d95f02",
        "Moderate effect, Narrower CI" = "#1f77b4",
        "Moderate effect, Wider CI" = "#6baed6",
        "Small effect, Narrower CI" = "#6a3d9a",
        "Small effect, Wider CI" = "#7f7f7f"
      ),
      drop = FALSE
    ) +
    labs(
      title = title,
      subtitle = glue("Primary scale: {PRIMARY_TRANSFORM}; {round(100 * CI_LEVEL)}% bootstrap CI"),
      x = "Cohen's d (Treatment delta vs Control delta)",
      y = NULL,
      color = "Effect/precision"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.major.y = element_line(color = "gray92"))
}

make_all_pairs_plot <- function(pairs_df, title = "All Genus-Pathway Delta Correlations (Supplementary)") {
  if (is.null(pairs_df) || nrow(pairs_df) == 0) return(NULL)

  plot_df <- pairs_df %>%
    mutate(
      pair_label = paste0(genus_short, " -> ", pathway_short),
      pair_label = reorder(pair_label, spearman_rho),
      perm_sig = ifelse(p_perm < 0.05, "perm p < 0.05", "perm p >= 0.05")
    )

  ggplot(plot_df, aes(y = pair_label, x = spearman_rho, color = perm_sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-0.8, 0.8), linetype = "dotted", color = "gray65") +
    geom_errorbar(aes(xmin = rho_ci_low, xmax = rho_ci_high), orientation = "y", width = 0.16, alpha = 0.8) +
    geom_point(size = 1.9, alpha = 0.9) +
    scale_color_manual(values = c("perm p < 0.05" = "#b22222", "perm p >= 0.05" = "#7f7f7f")) +
    labs(
      title = title,
      subtitle = glue("All {nrow(plot_df)} tested pairs; BH q-values in table output"),
      x = "Spearman rho",
      y = NULL,
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.y = element_line(color = "gray92"), legend.position = "bottom")
}

summarize_panel <- function(df, label, top_n = 6, include_clr = FALSE) {
  if (is.null(df) || nrow(df) == 0) {
    return(c(glue("### {label}"), "No analyzable features after filtering.", ""))
  }

  top <- df %>% arrange(rank_consensus, abs_d_rank) %>% slice_head(n = top_n)
  out <- c(
    glue("### {label}"),
    glue("- Features analyzed: {nrow(df)}"),
    glue("- Median |d|: {round(median(abs(df$cohens_d), na.rm = TRUE), 2)}"),
    glue("- Median CI width: {round(median(df$ci_width, na.rm = TRUE), 3)}"),
    glue("- Median LOO sign agreement: {round(median(df$loo_sign_agreement, na.rm = TRUE), 2)}"),
    glue("- Raw-vs-log1p sign agreement: {round(mean(df$transform_sign_agreement_raw_log1p, na.rm = TRUE), 2)}"),
    "- Top features (effect size + precision + LOO sensitivity):"
  )

  for (i in seq_len(nrow(top))) {
    line <- glue(
      "  - {top$feature_short[i]}: d={round(top$cohens_d[i], 2)}, CI [{round(top$ci_lower[i], 2)}, {round(top$ci_upper[i], 2)}], width={round(top$ci_width[i], 2)}, LOO sign agreement={round(top$loo_sign_agreement[i], 2)}, d_raw={round(top$d_raw[i], 2)}, d_log1p={round(top$d_log1p[i], 2)}"
    )
    if (include_clr && "d_clr" %in% colnames(top)) {
      line <- paste0(line, glue(", d_clr={round(top$d_clr[i], 2)}, primary-vs-CLR sign agreement={top$sign_agreement_primary_vs_clr[i]}"))
    }
    out <- c(out, line)
  }
  c(out, "")
}

# ------------------------- Data Loading -------------------------
cat("\n============================\n")
cat("Loading input tables\n")
cat("============================\n")

species_abundance <- load_first_available(
  github_paths = c(
    "results/Processed/02_taxonomic_profiling/species_profile_filtered.csv",
    "results/Reports/species_profile_filtered.csv"
  ),
  local_paths = c(
    file.path("results", "04_short_read_taxonomic_profiling", "merged_abundance_table.txt"),
    file.path("officer_output", "Reports", "species_profile_filtered.csv")
  ),
  label = "species abundance"
)

pathway_abundance <- tryCatch(
  load_first_available(
    github_paths = c(
      "results/Processed/03_functional_profiling/functional_profile_filtered.csv",
      "results/Reports/pathway_abundance.tsv"
    ),
    local_paths = c(
      file.path("results", "05_short_read_functional_profiling", "HUMAnN_merged_pathabundance_cpm.tsv"),
      file.path("officer_output", "Reports", "pathway_abundance.tsv")
    ),
    label = "pathway abundance"
  ),
  error = function(e) {
    cat("Pathway table not found; functional analysis skipped.\n")
    NULL
  }
)

feature_col <- colnames(species_abundance)[1]
colnames(species_abundance)[1] <- "feature"
feature_col <- "feature"

sample_cols <- colnames(species_abundance)[-1]
sample_cols_norm <- normalize_sample_id(sample_cols)
colnames(species_abundance) <- c(feature_col, sample_cols_norm)

metadata <- prepare_metadata(sample_cols_norm)
keep_samples <- intersect(metadata$sample_id, colnames(species_abundance))
metadata <- metadata %>% filter(sample_id %in% keep_samples)
species_abundance <- species_abundance %>% select(all_of(c(feature_col, keep_samples)))
species_abundance <- coerce_numeric_samples(species_abundance, keep_samples)

paired_design <- metadata %>%
  filter(timepoint %in% c("Baseline", "6_months")) %>%
  select(participant, turmeric_group, timepoint, sample_id) %>%
  pivot_wider(names_from = timepoint, values_from = sample_id) %>%
  filter(!is.na(Baseline), !is.na(`6_months`))

cat("\nStudy design summary\n")
cat("  Samples retained:", nrow(metadata), "\n")
cat("  Participants retained:", length(unique(metadata$participant)), "\n")
cat("  Paired baseline->6m:", nrow(paired_design), "\n")
cat("  Treatment pairs:", sum(paired_design$turmeric_group == "Treatment"), "\n")
cat("  Control pairs:", sum(paired_design$turmeric_group == "Control"), "\n")

# ------------------------- Taxonomy + Pathway Processing -------------------------
species_only <- species_abundance %>% filter(str_detect(.data[[feature_col]], TAXON_FEATURE_PATTERN))
species_only <- filter_feature_table(species_only, keep_samples, MIN_PREVALENCE_TAXA, MIN_MEAN_ABUND_TAXA)
species_only$genus <- sapply(species_only[[feature_col]], extract_rank, prefix = "g")

genus_df <- species_only %>%
  filter(!is.na(genus), genus != "", genus != "unclassified") %>%
  select(genus, all_of(keep_samples)) %>%
  group_by(genus) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

genus_df <- filter_feature_table(genus_df, keep_samples, MIN_PREVALENCE_TAXA, MIN_MEAN_ABUND_TAXA)

cat("\nTaxonomy filtering\n")
cat(" ", TAXON_LEVEL_LABEL, " entries retained:", nrow(species_only), "\n", sep = "")
cat("  Genera retained:", nrow(genus_df), "\n")

pathway_df <- NULL
if (!is.null(pathway_abundance)) {
  pathway_col <- colnames(pathway_abundance)[1]
  colnames(pathway_abundance)[1] <- "pathway"

  p_samples <- normalize_sample_id(colnames(pathway_abundance)[-1])
  colnames(pathway_abundance) <- c("pathway", p_samples)

  keep_p <- intersect(keep_samples, colnames(pathway_abundance))
  pathway_abundance <- pathway_abundance %>% select(all_of(c("pathway", keep_p)))
  pathway_abundance <- coerce_numeric_samples(pathway_abundance, keep_p)

  pathway_df <- pathway_abundance %>%
    filter(!str_detect(pathway, "^UNMAPPED$|^UNINTEGRATED$")) %>%
    filter(!str_detect(pathway, "\\|"))

  pathway_df <- filter_feature_table(pathway_df, keep_p, MIN_PREVALENCE_PATHWAY, MIN_MEAN_ABUND_PATHWAY)
  cat("\nFunctional filtering\n")
  cat("  Pathways retained:", nrow(pathway_df), "\n")
}

# ------------------------- Panel Definition -------------------------
cat("\n============================\n")
cat("Defining objective panels\n")
cat("============================\n")

literature_genera <- c(
  "Bacteroides", "Bifidobacterium", "Alistipes", "Parabacteroides", "Blautia"
)
available_lit_genera <- intersect(literature_genera, genus_df$genus)

top_abundant_genera <- genus_df %>%
  mutate(mean_abundance = rowMeans(across(all_of(keep_samples)), na.rm = TRUE)) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = MAX_ABUNDANT_GENUS) %>%
  pull(genus)

panel_genera <- if (GENUS_PANEL_MODE == "pilot_plus_abundant") {
  unique(c(available_lit_genera, top_abundant_genera))
} else {
  available_lit_genera
}
if (length(panel_genera) == 0) {
  stop("No pilot-prioritized genera were retained after filtering.")
}

species_with_stats <- species_only %>%
  mutate(
    mean_abundance = rowMeans(across(all_of(keep_samples)), na.rm = TRUE),
    genus = sapply(.data[[feature_col]], extract_rank, prefix = "g")
  ) %>%
  filter(genus %in% panel_genera) %>%
  arrange(desc(mean_abundance))

if (USE_ALL_RETAINED_FOR_INFERENCE) {
  panel_species <- species_with_stats %>% pull(.data[[feature_col]])
} else {
  panel_species <- species_with_stats %>%
    slice_head(n = MAX_SPECIES_PANEL) %>%
    pull(.data[[feature_col]])
}

panel_pathways <- character(0)
if (!is.null(pathway_df) && nrow(pathway_df) > 0) {
  pathway_ranked <- pathway_df %>%
    mutate(mean_abundance = rowMeans(across(where(is.numeric)), na.rm = TRUE))

  if (FILTER_PATHWAYS_BY_KEYWORD) {
    pathway_focus <- pathway_ranked %>%
      filter(str_detect(str_to_upper(pathway), paste(PATHWAY_KEYWORDS, collapse = "|"))) %>%
      arrange(desc(mean_abundance))
    if (nrow(pathway_focus) == 0) {
      warning("No keyword-matched pathways found after filtering; falling back to all retained pathways.")
      pathway_focus <- pathway_ranked %>% arrange(desc(mean_abundance))
    }
  } else {
    pathway_focus <- pathway_ranked %>% arrange(desc(mean_abundance))
  }

  if (USE_ALL_RETAINED_FOR_INFERENCE) {
    panel_pathways <- pathway_focus %>% pull(pathway)
  } else {
    panel_pathways <- pathway_focus %>% slice_head(n = MAX_PATHWAY_PANEL) %>% pull(pathway)
  }
}

cat("Panel sizes\n")
cat("  Genus panel mode:", GENUS_PANEL_MODE, "\n")
cat("  Pilot genera available:", length(available_lit_genera), "of", length(literature_genera), "\n")
cat("  Genus panel:", length(panel_genera), "\n")
cat(" ", TAXON_LEVEL_LABEL, " panel:", length(panel_species), "\n", sep = "")
cat("  Pathway panel:", length(panel_pathways), "\n")
cat("  Pathway selection:", ifelse(FILTER_PATHWAYS_BY_KEYWORD, "keyword-matched retained pathways", "all retained pathways"), "\n")
cat("  Inference mode:", ifelse(USE_ALL_RETAINED_FOR_INFERENCE, "all retained features", "capped feature panels"), "\n")

# ------------------------- Main Analyses -------------------------
cat("\n============================\n")
cat("Running effect-size analyses\n")
cat("============================\n")

species_results <- analyze_feature_panel(species_only, feature_col, panel_species, metadata, TAXON_LEVEL_LABEL)
genus_results <- analyze_feature_panel(genus_df, "genus", panel_genera, metadata, "Genus")

pathway_results <- NULL
if (!is.null(pathway_df) && length(panel_pathways) > 0) {
  pathway_results <- analyze_feature_panel(pathway_df, "pathway", panel_pathways, metadata, "Pathway")
}

if (RUN_CLR_SENSITIVITY) {
  cat("\nRunning CLR sensitivity for taxonomy\n")
  clr_species <- compute_clr_table(species_only, feature_col, keep_samples)
  clr_genus <- compute_clr_table(genus_df, "genus", keep_samples)

  species_clr <- compute_panel_clr_effects(clr_species, feature_col, panel_species, metadata)
  genus_clr <- compute_panel_clr_effects(clr_genus, "genus", panel_genera, metadata)

  species_results <- species_results %>%
    left_join(species_clr, by = "feature") %>%
    mutate(sign_agreement_primary_vs_clr = case_when(
      is.na(d_clr) ~ NA,
      sign(cohens_d) == sign(d_clr) ~ TRUE,
      TRUE ~ FALSE
    ))

  genus_results <- genus_results %>%
    left_join(genus_clr, by = "feature") %>%
    mutate(sign_agreement_primary_vs_clr = case_when(
      is.na(d_clr) ~ NA,
      sign(cohens_d) == sign(d_clr) ~ TRUE,
      TRUE ~ FALSE
    ))
}

# ------------------------- Taxonomy-Function Integration -------------------------
candidate_pairs <- NULL
if (!is.null(pathway_results) && nrow(pathway_results) > 0 && nrow(genus_results) > 0) {
  cat("\nRunning taxonomy-function integration (exploratory)\n")

  top_genus <- genus_results %>% arrange(rank_consensus, abs_d_rank) %>% slice_head(n = TOP_FOR_INTEGRATION_GENUS) %>% pull(feature)
  top_path <- pathway_results %>% arrange(rank_consensus, abs_d_rank) %>% slice_head(n = TOP_FOR_INTEGRATION_PATHWAY) %>% pull(feature)

  pairs <- expand.grid(genus = top_genus, pathway = top_path, stringsAsFactors = FALSE)

  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    g <- pairs$genus[i]
    p <- pairs$pathway[i]

    g_delta <- extract_feature_delta(genus_df, "genus", g, metadata, transform_mode = PRIMARY_TRANSFORM)
    p_delta <- extract_feature_delta(pathway_df, "pathway", p, metadata, transform_mode = PRIMARY_TRANSFORM)
    if (is.null(g_delta) || is.null(p_delta)) return(NULL)

    merged <- g_delta %>%
      select(participant, turmeric_group, g_delta = delta) %>%
      inner_join(p_delta %>% select(participant, p_delta = delta), by = "participant")

    if (nrow(merged) < 6) return(NULL)

    rho <- suppressWarnings(cor(merged$g_delta, merged$p_delta, method = "spearman", use = "complete.obs"))
    if (!is.finite(rho)) return(NULL)

    ci <- bootstrap_spearman_ci(merged$g_delta, merged$p_delta)
    p_perm <- perm_spearman_p(merged$g_delta, merged$p_delta)

    n_tx <- sum(merged$turmeric_group == "Treatment")
    n_ct <- sum(merged$turmeric_group == "Control")

    rho_tx <- if (n_tx >= 4) suppressWarnings(cor(merged$g_delta[merged$turmeric_group == "Treatment"], merged$p_delta[merged$turmeric_group == "Treatment"], method = "spearman")) else NA_real_
    rho_ct <- if (n_ct >= 4) suppressWarnings(cor(merged$g_delta[merged$turmeric_group == "Control"], merged$p_delta[merged$turmeric_group == "Control"], method = "spearman")) else NA_real_

    g_d <- genus_results$cohens_d[genus_results$feature == g][1]
    p_d <- pathway_results$cohens_d[pathway_results$feature == p][1]

    data.frame(
      genus = g,
      pathway = p,
      genus_short = extract_short_name(g),
      pathway_short = str_trunc(gsub("^.*: ", "", p), 80),
      n_participants = nrow(merged),
      n_treatment = n_tx,
      n_control = n_ct,
      spearman_rho = rho,
      rho_ci_low = ci[1],
      rho_ci_high = ci[2],
      rho_ci_width = ci[2] - ci[1],
      p_perm = p_perm,
      rho_treatment = rho_tx,
      rho_control = rho_ct,
      genus_d = g_d,
      pathway_d = p_d,
      same_effect_direction = sign(g_d) == sign(p_d),
      stringsAsFactors = FALSE
    )
  })

  candidate_pairs <- bind_rows(rows)
  if (nrow(candidate_pairs) > 0) {
    candidate_pairs <- candidate_pairs %>%
      mutate(
        abs_rho = abs(spearman_rho),
        q_perm = p.adjust(p_perm, method = "BH"),
        rank_abs_rho = dense_rank(desc(abs_rho)),
        rank_ci = dense_rank(rho_ci_width),
        rank_effect_support = dense_rank(desc(pmin(abs(genus_d), abs(pathway_d)))),
        rank_consensus = rank_abs_rho + rank_ci + rank_effect_support,
        strong_corr_flag = abs_rho >= 0.8
      ) %>%
      arrange(rank_consensus, rank_abs_rho)
  }
}

# ------------------------- Save Outputs -------------------------
species_path <- file.path(output_dir, paste0(TAXON_FILE_PREFIX, "_effect_sizes.csv"))
genus_path <- file.path(output_dir, "genus_effect_sizes.csv")
fwrite(species_results, species_path)
fwrite(genus_results, genus_path)

pathway_path <- NULL
if (!is.null(pathway_results)) {
  pathway_path <- file.path(output_dir, "pathway_effect_sizes.csv")
  fwrite(pathway_results, pathway_path)
}

pairs_path <- NULL
if (!is.null(candidate_pairs) && nrow(candidate_pairs) > 0) {
  pairs_path <- file.path(output_dir, "taxonomy_function_candidate_pairs.csv")
  fwrite(candidate_pairs, pairs_path)
}

species_plot <- make_forest_plot(species_results, glue("{TAXON_LEVEL_LABEL}-Level Effect Estimates"))
genus_plot <- make_forest_plot(genus_results, "Genus-Level Effect Estimates")
pathway_plot <- make_forest_plot(pathway_results, "Pathway-Level Effect Estimates")

species_plot_path <- file.path(output_dir, paste0(TAXON_FILE_PREFIX, "_forest_plot.png"))
genus_plot_path <- file.path(output_dir, "genus_forest_plot.png")
pathway_plot_path <- file.path(output_dir, "pathway_forest_plot.png")

if (!is.null(species_plot)) ggsave(species_plot_path, species_plot, width = 10, height = 8, dpi = 300, bg = "white")
if (!is.null(genus_plot)) ggsave(genus_plot_path, genus_plot, width = 10, height = 8, dpi = 300, bg = "white")
if (!is.null(pathway_plot)) ggsave(pathway_plot_path, pathway_plot, width = 11, height = 8, dpi = 300, bg = "white")

# Supplementary plots: all analyzed features/pairs (not top-N only)
species_plot_all <- make_forest_plot(
  species_results,
  glue("{TAXON_LEVEL_LABEL}-Level Effect Estimates (All Analyzed Features)"),
  top_n = ifelse(is.null(species_results), TOP_FOR_PLOTS, nrow(species_results))
)
genus_plot_all <- make_forest_plot(
  genus_results,
  "Genus-Level Effect Estimates (All Analyzed Features)",
  top_n = ifelse(is.null(genus_results), TOP_FOR_PLOTS, nrow(genus_results))
)
pathway_plot_all <- make_forest_plot(
  pathway_results,
  "Pathway-Level Effect Estimates (All Analyzed Features)",
  top_n = ifelse(is.null(pathway_results), TOP_FOR_PLOTS, nrow(pathway_results))
)
pairs_plot_all <- make_all_pairs_plot(candidate_pairs)

species_plot_all_png <- file.path(output_dir, paste0("supplementary_", TAXON_FILE_PREFIX, "_forest_all_features.png"))
species_plot_all_pdf <- file.path(output_dir, paste0("supplementary_", TAXON_FILE_PREFIX, "_forest_all_features.pdf"))
genus_plot_all_png <- file.path(output_dir, "supplementary_genus_forest_all_features.png")
genus_plot_all_pdf <- file.path(output_dir, "supplementary_genus_forest_all_features.pdf")
pathway_plot_all_png <- file.path(output_dir, "supplementary_pathway_forest_all_features.png")
pathway_plot_all_pdf <- file.path(output_dir, "supplementary_pathway_forest_all_features.pdf")
pairs_plot_all_png <- file.path(output_dir, "supplementary_taxonomy_function_all_pairs.png")
pairs_plot_all_pdf <- file.path(output_dir, "supplementary_taxonomy_function_all_pairs.pdf")

if (!is.null(species_plot_all)) {
  h <- max(8, min(28, 0.24 * nrow(species_results)))
  ggsave(species_plot_all_png, species_plot_all, width = 11, height = h, dpi = 300, bg = "white")
  ggsave(species_plot_all_pdf, species_plot_all, width = 11, height = h, bg = "white")
}
if (!is.null(genus_plot_all)) {
  h <- max(8, min(18, 0.28 * nrow(genus_results)))
  ggsave(genus_plot_all_png, genus_plot_all, width = 11, height = h, dpi = 300, bg = "white")
  ggsave(genus_plot_all_pdf, genus_plot_all, width = 11, height = h, bg = "white")
}
if (!is.null(pathway_plot_all)) {
  h <- max(10, min(24, 0.24 * nrow(pathway_results)))
  ggsave(pathway_plot_all_png, pathway_plot_all, width = 12, height = h, dpi = 300, bg = "white")
  ggsave(pathway_plot_all_pdf, pathway_plot_all, width = 12, height = h, bg = "white")
}
if (!is.null(pairs_plot_all)) {
  h <- max(14, min(60, 0.16 * nrow(candidate_pairs)))
  ggsave(pairs_plot_all_png, pairs_plot_all, width = 12, height = h, dpi = 300, bg = "white")
  ggsave(pairs_plot_all_pdf, pairs_plot_all, width = 12, height = h, bg = "white")
}

robust_df <- bind_rows(species_results, genus_results, pathway_results) %>%
  filter(!is.na(cohens_d), !is.na(loo_sign_agreement), !is.na(ci_width))

robust_plot_path <- file.path(output_dir, "robustness_summary.png")
if (nrow(robust_df) > 0) {
  p_robust <- ggplot(robust_df, aes(x = abs(cohens_d), y = loo_sign_agreement, color = level)) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = 1 / (ci_width + 1e-6)), alpha = 0.65) +
    labs(
      title = "Effect Magnitude vs LOO Sensitivity",
      subtitle = "LOO sign agreement is sensitivity, not external validation",
      x = "Absolute Cohen's d",
      y = "LOO sign agreement",
      size = "Precision\n(1 / CI width)",
      color = "Level"
    ) +
    theme_bw(base_size = 11)

  ggsave(robust_plot_path, p_robust, width = 10, height = 7, dpi = 300, bg = "white")
}

# ------------------------- Interpretation Document -------------------------
interpretation <- c(
  "# Honest Small-Sample Microbiome and Functional Analysis",
  "",
  "## Study Design",
  glue("- Retained participants: {length(unique(metadata$participant))}"),
  glue("- Paired baseline->6 month participants: {nrow(paired_design)}"),
  glue("- Treatment pairs: {sum(paired_design$turmeric_group == 'Treatment')}, Control pairs: {sum(paired_design$turmeric_group == 'Control')}"),
  "",
  "## Statistical Position",
  "- This is an effect-estimation analysis for a low-N study.",
  "- No feature-level significance threshold claims are made.",
  "- Uncertainty and sensitivity diagnostics are reported explicitly.",
  "",
  "## Transformation and Compositional Strategy",
  glue("- Primary effect scale: {PRIMARY_TRANSFORM} (for interpretability)."),
  glue("- Secondary transform sensitivity: {SECONDARY_TRANSFORM} (directional concordance reported)."),
  glue("- Taxonomy CLR sensitivity enabled: {RUN_CLR_SENSITIVITY}."),
  "",
  "## Methods (Concise)",
  glue("- Taxonomy filtering: prevalence >= {MIN_PREVALENCE_TAXA}, mean abundance >= {MIN_MEAN_ABUND_TAXA}"),
  glue("- Pathway filtering: prevalence >= {MIN_PREVALENCE_PATHWAY}, mean abundance >= {MIN_MEAN_ABUND_PATHWAY}"),
  "- Paired delta per participant: value_6m - value_baseline.",
  glue("- Effect metrics: Cohen's d, probability of superiority, {round(100 * CI_LEVEL)}% bootstrap CI ({N_BOOTSTRAP} resamples)."),
  "- LOO summary: sign agreement with full-sample effect (participant-sensitivity, not validation).",
  "- Feature ranking: transparent multi-rank (effect size, CI width, LOO sign agreement).",
  ""
)

interpretation <- c(interpretation, summarize_panel(species_results, glue("{TAXON_LEVEL_LABEL} Panel"), top_n = 6, include_clr = RUN_CLR_SENSITIVITY))
interpretation <- c(interpretation, summarize_panel(genus_results, "Genus Panel", top_n = 6, include_clr = RUN_CLR_SENSITIVITY))
interpretation <- c(interpretation, summarize_panel(pathway_results, "Pathway Panel", top_n = 6, include_clr = FALSE))

if (!is.null(candidate_pairs) && nrow(candidate_pairs) > 0) {
  top_pairs <- candidate_pairs %>% slice_head(n = 10)
  interpretation <- c(
    interpretation,
    "## Taxonomy-Function Integration (Exploratory)",
    "- Delta-based Spearman correlations with bootstrap CI and permutation p-values are provided.",
    "- BH-adjusted q-values are included as multiplicity context across tested pairs.",
    "- Low N means these are hypothesis-prioritization signals, not confirmatory findings.",
    "- Top candidate pairs:"
  )

  for (i in seq_len(nrow(top_pairs))) {
    interpretation <- c(
      interpretation,
      glue("  - {top_pairs$genus_short[i]} -> {top_pairs$pathway_short[i]}: rho={round(top_pairs$spearman_rho[i], 2)}, CI [{round(top_pairs$rho_ci_low[i], 2)}, {round(top_pairs$rho_ci_high[i], 2)}], p_perm={signif(top_pairs$p_perm[i], 3)}, q_perm={signif(top_pairs$q_perm[i], 3)}, n={top_pairs$n_participants[i]}, rho_tx={round(top_pairs$rho_treatment[i], 2)}, rho_ct={round(top_pairs$rho_control[i], 2)}")
    )
  }
  interpretation <- c(interpretation, "")
}

interpretation <- c(
  interpretation,
  "## What This Supports",
  "- Ranked hypotheses for follow-up based on effect size, precision, and participant-level sensitivity.",
  "- Directional signal stability checks across transform and compositional sensitivity analyses.",
  "- Supplementary Figures S1-S4 map to complete panels: S1 (all strain-resolved features), S2 (all genera), S3 (all pathways), S4 (all tested genus-pathway pairs).",
  "",
  "## What This Does Not Support",
  "- Confirmatory causal claims from this dataset alone.",
  "- Population-level generalization without independent replication.",
  "",
  "## Recommended Next Step",
  "- Replicate top-ranked taxonomy and pathway candidates in a larger independent cohort.",
  "",
  glue("Report generated: {Sys.Date()}")
)

interpretation_path <- file.path(output_dir, "INTERPRETATION_GUIDE.md")
writeLines(interpretation, interpretation_path)

# ------------------------- Optional Upload -------------------------
if (nzchar(token) && isTRUE(ENABLE_GITHUB_UPLOAD)) {
  upload_to_github(species_path, to_github_path(file.path(output_dir, basename(species_path))), glue("Update honest small-N {tolower(TAXON_LEVEL_LABEL)} effect sizes"))
  upload_to_github(genus_path, to_github_path(file.path(output_dir, basename(genus_path))), "Update honest small-N genus effect sizes")

  if (!is.null(pathway_path) && file.exists(pathway_path)) {
    upload_to_github(pathway_path, to_github_path(file.path(output_dir, basename(pathway_path))), "Update honest small-N pathway effect sizes")
  }
  if (!is.null(pairs_path) && file.exists(pairs_path)) {
    upload_to_github(pairs_path, to_github_path(file.path(output_dir, basename(pairs_path))), "Update taxonomy-function candidate pairs")
  }

  if (file.exists(species_plot_path)) upload_to_github(species_plot_path, to_github_path(file.path(output_dir, basename(species_plot_path))), glue("Update {tolower(TAXON_LEVEL_LABEL)} forest plot"))
  if (file.exists(genus_plot_path)) upload_to_github(genus_plot_path, to_github_path(file.path(output_dir, basename(genus_plot_path))), "Update genus forest plot")
  if (file.exists(pathway_plot_path)) upload_to_github(pathway_plot_path, to_github_path(file.path(output_dir, basename(pathway_plot_path))), "Update pathway forest plot")
  if (file.exists(robust_plot_path)) upload_to_github(robust_plot_path, to_github_path(file.path(output_dir, basename(robust_plot_path))), "Update LOO sensitivity summary plot")
  upload_to_github(interpretation_path, to_github_path(file.path(output_dir, basename(interpretation_path))), "Update honest small-N interpretation guide")
}

cat("\nAnalysis complete. Outputs in:", output_dir, "\n")

