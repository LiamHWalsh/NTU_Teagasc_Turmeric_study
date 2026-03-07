# ────────────────────────────────────────────────────────────────────────────────
# Script: 09_visualize_maaslin3_results_SMALL_SAMPLE.R
# Description: Small-sample optimized validation of MaAsLin3 results
# Author: Liam Walsh
# Date: 2025
# ────────────────────────────────────────────────────────────────────────────────

# ────────────────────────────────────────────────────────────────────────────────
# Section: Setup
# ────────────────────────────────────────────────────────────────────────────────

# Use a writable user library if the preferred path is unavailable.
preferred_lib <- "G:/R/win-library/4.4"
if (!dir.exists(preferred_lib)) {
  preferred_lib <- normalizePath(file.path(Sys.getenv("USERPROFILE"), "Documents", "R", "win-library", "4.4"),
                                 winslash = "/", mustWork = FALSE)
}
dir.create(preferred_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(preferred_lib, .libPaths()))
set.seed(123)

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cran.rstudio.com/", lib = .libPaths()[1])
}

suppressPackageStartupMessages({
  pacman::p_load(
    readxl, readr, dplyr, tidyr, stringr, tidyverse, httr,
    data.table, ggplot2, viridis, gh, base64enc,
    patchwork, glue, utf8, scales,
    ggbeeswarm, ggsignif, coin, rstatix, ggpubr,
    exactRankTests, boot
  )
})

# ────────────────────────────────────────────────────────────────────────────────
# Section: GitHub Configuration
# ────────────────────────────────────────────────────────────────────────────────

owner <- "LiamHWalsh"
repo <- "NTU_Teagasc_Turmeric_study"
branch <- "main"

pipeline_paths <- list(
  taxonomic = "results/Processed/02_taxonomic_profiling",
  functional = "results/Processed/03_functional_profiling",
  statistics = "results/Processed/05_statistical_analysis"
)

figure_paths <- list(
  statistics = "figures/05_statistical_analysis"
)

output_path_diff <- file.path(pipeline_paths$statistics, "differential_abundance")
figure_output_diff <- file.path(figure_paths$statistics, "differential_abundance")

# Small-sample analysis defaults (explicitly split confirmatory vs exploratory reporting).
PRIMARY_ALPHA <- 0.05
EXPLORATORY_Q_ALPHA <- 0.25
SCREEN_P_THRESHOLD <- 0.20
MIN_PREVALENCE_PROP <- 0.10
MIN_PAIRED_PARTICIPANTS <- 4
MIN_GROUP_N <- 2
MAX_FEATURES_TO_VALIDATE <- 25
BOOTSTRAP_RESAMPLES <- 10000
PERMUTATION_RESAMPLES <- 10000

token <- Sys.getenv("GITHUB_PAT")
if (token == "") token <- Sys.getenv("GITHUB_TOKEN")
allow_github_write <- nzchar(token)

if (!allow_github_write) {
  cat("⚠️  No GitHub token found in GITHUB_PAT or GITHUB_TOKEN.\n")
  cat("    Running in read-only mode (uploads disabled).\n")
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: GitHub Functions
# ────────────────────────────────────────────────────────────────────────────────

load_github_file <- function(repo_path, token = "", owner, repo, branch = "main") {
  file_text <- NULL

  if (nzchar(token)) {
    file_info <- gh::gh(
      "GET /repos/:owner/:repo/contents/:path",
      owner = owner, repo = repo, path = repo_path, ref = branch, .token = token
    )
    if (is.null(file_info$content)) {
      stop(glue::glue("GitHub API returned no content for: {repo_path}"))
    }
    decoded <- base64enc::base64decode(gsub("\\n", "", file_info$content))
    file_text <- rawToChar(decoded)
  } else {
    file_url <- paste0("https://raw.githubusercontent.com/", owner, "/", repo, "/", branch, "/", repo_path)
    response <- GET(file_url)
    if (status_code(response) != 200) {
      stop(glue::glue("Failed to load file: {file_url}"))
    }
    file_text <- rawToChar(response$content)
  }

  file_ext <- tolower(tools::file_ext(repo_path))
  if (file_ext == "csv") {
    return(readr::read_csv(I(file_text), show_col_types = FALSE, name_repair = "minimal"))
  }
  if (file_ext %in% c("tsv", "txt")) {
    return(readr::read_tsv(I(file_text), show_col_types = FALSE, comment = "#", name_repair = "minimal"))
  }
  readr::read_tsv(I(file_text), show_col_types = FALSE, comment = "#", name_repair = "minimal")
}

load_github_first_available <- function(repo_paths, token, owner, repo, branch, label) {
  for (repo_path in repo_paths) {
    out <- tryCatch(
      load_github_file(repo_path, token = token, owner = owner, repo = repo, branch = branch),
      error = function(e) NULL
    )
    if (!is.null(out)) {
      cat("✅", label, "loaded from GitHub:", repo_path, "\n")
      return(out)
    }
  }
  NULL
}

save_and_upload_to_github <- function(object = NULL, object_type = c("data", "plot", "file"), 
                                      local_path = NULL, github_path, commit_message,
                                      owner, repo, branch = "main") {
  object_type <- match.arg(object_type)
  gh_token <- Sys.getenv("GITHUB_PAT")
  if (gh_token == "") gh_token <- Sys.getenv("GITHUB_TOKEN")
  if (gh_token == "") {
    message("⚠️  Upload skipped (missing GitHub token): ", github_path)
    return(invisible(NULL))
  }
  
  if (is.null(local_path)) {
    file_ext <- switch(object_type, 
                       data = ".csv",
                       plot = ".png",
                       stop("For 'file' type, 'local_path' must be provided."))
    local_path <- tempfile(fileext = file_ext)
  }
  
  if (object_type == "data") {
    if (is.null(object)) stop("❌ No data object provided.")
    data.table::fwrite(object, local_path, row.names = FALSE)
  } else if (object_type == "plot") {
    if (is.null(object)) stop("❌ No plot object provided.")
    ggplot2::ggsave(local_path, plot = object, width = 12, height = 8, dpi = 300, bg = "white")
  } else if (object_type == "file") {
    if (!file.exists(local_path)) stop(glue::glue("❌ File does not exist: {local_path}"))
  }
  
  encoded_content <- base64enc::base64encode(local_path)
  file_info <- tryCatch(
    gh::gh("GET /repos/:owner/:repo/contents/:path",
           owner = owner, repo = repo, path = github_path, ref = branch, .token = gh_token),
    error = function(e) {
      if (grepl("404", e$message)) return(NULL)
      stop(e)
    }
  )
  
  sha <- if (!is.null(file_info)) file_info$sha else NULL
  
  gh::gh("PUT /repos/:owner/:repo/contents/:path",
         owner = owner, repo = repo, path = github_path,
         message = commit_message, content = encoded_content,
         branch = branch, sha = sha, .token = gh_token)
  
  message(glue::glue("🚀 File uploaded: https://github.com/{owner}/{repo}/blob/{branch}/{github_path}"))
  
  if (object_type %in% c("data", "plot")) unlink(local_path)
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: Helper Functions
# ────────────────────────────────────────────────────────────────────────────────

extract_species_name <- function(feature_string) {
  parts <- str_split(feature_string, "\\|")[[1]]
  species_part <- parts[str_detect(parts, "^s__")]
  if (length(species_part) > 0) {
    return(str_remove(species_part, "^s__"))
  }
  return(feature_string)
}

simplify_feature_names <- function(df) {
  df %>%
    mutate(
      feature_simple = sapply(feature, extract_species_name),
      feature_short = str_trunc(feature_simple, width = 50, ellipsis = "...")
    )
}

load_local_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    return(readr::read_csv(path, show_col_types = FALSE, name_repair = "minimal"))
  }
  if (ext %in% c("tsv", "txt")) {
    return(readr::read_tsv(path, show_col_types = FALSE, comment = "#", name_repair = "minimal"))
  }
  readr::read_tsv(path, show_col_types = FALSE, comment = "#", name_repair = "minimal")
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: Load Data
# ────────────────────────────────────────────────────────────────────────────────

github_base <- paste0("https://raw.githubusercontent.com/", owner, "/", repo, "/", branch)

cat("\n═══════════════════════════════════════\n")
cat("Loading Data (GitHub with local fallback)\n")
cat("═══════════════════════════════════════\n")

profiles <- c("species", "functional")
maaslin_data <- list()
abundance_data <- list()

for (profile in profiles) {
  if (profile == "species") {
    maaslin_candidates <- c(
      "results/Processed/05_statistical_analysis/differential_abundance/species_maaslin3_all_results.csv",
      "results/04_short_read_taxonomic_profiling/maaslin3_results/maaslin3_all_results.tsv"
    )
    abundance_candidates <- c(
      "results/Processed/02_taxonomic_profiling/species_profile_filtered.csv",
      "results/04_short_read_taxonomic_profiling/merged_abundance_table.txt"
    )
  } else {
    maaslin_candidates <- c(
      "results/Processed/05_statistical_analysis/differential_abundance/functional_maaslin3_all_results.csv",
      "results/05_short_read_functional_profiling/maaslin3_results/maaslin3_all_results.tsv"
    )
    abundance_candidates <- c(
      "results/Processed/03_functional_profiling/functional_profile_filtered.csv",
      "results/05_short_read_functional_profiling/HUMAnN_merged_pathabundance_cpm.tsv"
    )
  }

  maaslin_data[[profile]] <- load_github_first_available(
    maaslin_candidates, token = token, owner = owner, repo = repo, branch = branch,
    label = paste(profile, "MaAsLin3")
  )
  if (is.null(maaslin_data[[profile]])) {
    cat("WARNING:", profile, "MaAsLin3 - no candidate GitHub paths were available\n")
  }

  abundance_data[[profile]] <- load_github_first_available(
    abundance_candidates, token = token, owner = owner, repo = repo, branch = branch,
    label = paste(profile, "abundance")
  )
  if (is.null(abundance_data[[profile]])) {
    cat("WARNING:", profile, "abundance - no candidate GitHub paths were available\n")
  }
}
# Local fallbacks when GitHub paths are unavailable in this environment.
if (is.null(maaslin_data[["species"]])) {
  local_species_maaslin <- file.path("results", "04_short_read_taxonomic_profiling", "maaslin3_results", "maaslin3_all_results.tsv")
  if (file.exists(local_species_maaslin)) {
    maaslin_data[["species"]] <- load_local_table(local_species_maaslin)
    cat("✅ species MaAsLin3 loaded from local fallback\n")
  }
}

if (is.null(abundance_data[["species"]])) {
  local_species_abundance <- file.path("results", "04_short_read_taxonomic_profiling", "merged_abundance_table.txt")
  if (file.exists(local_species_abundance)) {
    abundance_data[["species"]] <- load_local_table(local_species_abundance)
    cat("✅ species abundance loaded from local fallback\n")
  }
}

if (is.null(maaslin_data[["functional"]])) {
  local_functional_maaslin <- file.path("results", "05_short_read_functional_profiling", "maaslin3_results", "maaslin3_all_results.tsv")
  if (file.exists(local_functional_maaslin)) {
    maaslin_data[["functional"]] <- load_local_table(local_functional_maaslin)
    cat("✅ functional MaAsLin3 loaded from local fallback\n")
  }
}

if (is.null(abundance_data[["functional"]])) {
  local_functional_abundance <- file.path("results", "05_short_read_functional_profiling", "HUMAnN_merged_pathabundance_cpm.tsv")
  if (file.exists(local_functional_abundance)) {
    abundance_data[["functional"]] <- load_local_table(local_functional_abundance)
    cat("✅ functional abundance loaded from local fallback\n")
  }
}

# Load metadata
turmeric_status <- data.frame(
  ID = c(3,4,5,6,7,9,11,14,15,16,18,19,20,21,22,23),
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

prepare_metadata <- function(sample_ids) {
  split_ids <- str_split_fixed(sample_ids, "_", 3)
  
  metadata <- data.frame(
    sample_id = sample_ids,
    block_code = split_ids[,2],
    block_number = as.integer(gsub("[^0-9]", "", split_ids[,2])),
    block_letter = gsub("[0-9]", "", split_ids[,2]),
    stringsAsFactors = FALSE
  )
  
  metadata <- merge(metadata, turmeric_status,
                    by.x = "block_number", by.y = "ID", all.x = TRUE)

  # Drop controls/negatives or unmatched sample IDs before deriving groups.
  metadata <- metadata %>% filter(!is.na(block_number), !is.na(status))
  
  metadata$timepoint <- "Baseline"
  metadata$timepoint[metadata$block_letter %in% c("a","A")] <- "Baseline"
  metadata$timepoint[metadata$block_letter %in% c("b","B")] <- "2_weeks"
  metadata$timepoint[metadata$block_letter %in% c("c","C")] <- "6_months"
  
  metadata$turmeric_group <- "Control"
  metadata$turmeric_group[
    grepl("Turmeric 6 months", metadata$status) & 
      metadata$block_letter %in% c("c", "C")
  ] <- "Treatment"
  
  metadata$turmeric_group[
    metadata$status == "Turmeric 3 months" & 
      metadata$block_letter == "C"
  ] <- "Stopped"
  
  for (i in 1:nrow(metadata)) {
    if (!is.na(metadata$timepoint[i]) && metadata$timepoint[i] %in% c("Baseline", "2_weeks")) {
      participant_id <- metadata$block_number[i]
      has_treatment <- any(
        metadata$block_number == participant_id & 
          metadata$turmeric_group == "Treatment",
        na.rm = TRUE
      )
      if (has_treatment) metadata$turmeric_group[i] <- "Treatment"
    }
  }
  
  metadata$block_letter <- str_to_upper(metadata$block_letter)
  
  metadata$turmeric_status <- "Not taken"
  metadata$turmeric_status[
    metadata$block_letter == "C" & 
      metadata$status == "Turmeric 6 months"
  ] <- "Taken"
  
  metadata <- metadata %>%
    filter(turmeric_group != "Stopped") %>%
    mutate(
      participant = as.character(block_number),
      timepoint = factor(timepoint, levels = c("Baseline", "2_weeks", "6_months")),
      turmeric_group = factor(turmeric_group, levels = c("Control", "Treatment")),
      turmeric_status = factor(turmeric_status, levels = c("Not taken", "Taken"))
    )
  
  return(metadata)
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: Identify Trends
# ────────────────────────────────────────────────────────────────────────────────

identify_trends <- function(results_df, profile_type, screen_p_threshold = SCREEN_P_THRESHOLD,
                            min_non_zero = 0) {
  
  cat("\n═══════════════════════════════════════\n")
  cat(glue::glue("IDENTIFYING TRENDS: {toupper(profile_type)}\n"))
  cat("═══════════════════════════════════════\n")
  
  turmeric_metadata_name <- if ("metadata" %in% colnames(results_df) &&
                                any(results_df$metadata == "turmeric_status", na.rm = TRUE)) {
    "turmeric_status"
  } else {
    "turmeric_group"
  }

  results_df <- simplify_feature_names(results_df) %>%
    filter(metadata == turmeric_metadata_name, model == "abundance")

  if ("N_not_zero" %in% colnames(results_df) && min_non_zero > 0) {
    results_df <- results_df %>% filter(N_not_zero >= min_non_zero)
  }
  
  trends <- results_df %>%
    filter(pval_individual < screen_p_threshold) %>%
    mutate(
      pval_individual = as.numeric(pval_individual),
      qval_individual = as.numeric(qval_individual)
    ) %>%
    arrange(qval_individual, pval_individual, desc(abs(coef)))

  confirmatory_hits <- trends %>%
    filter(!is.na(qval_individual), qval_individual < PRIMARY_ALPHA)

  exploratory_hits <- trends %>%
    filter(!is.na(qval_individual), qval_individual < EXPLORATORY_Q_ALPHA)
  
  cat(glue::glue("Screening set (p < {screen_p_threshold}): {nrow(trends)}\n"))
  cat(glue::glue("Confirmatory hits (q < {PRIMARY_ALPHA}): {nrow(confirmatory_hits)}\n"))
  cat(glue::glue("Exploratory hits (q < {EXPLORATORY_Q_ALPHA}): {nrow(exploratory_hits)}\n\n"))
  
  if (nrow(trends) > 0) {
    cat("Top 15 screening candidates (not confirmatory):\n")
    print(trends %>%
            select(feature_short, coef, pval_individual, qval_individual, N_not_zero) %>%
            head(15))
  }
  
  return(list(
    all_trends = trends,
    turmeric_trends = trends,
    confirmatory_hits = confirmatory_hits,
    exploratory_hits = exploratory_hits
  ))
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: Paired Within-Subject Analysis
# ────────────────────────────────────────────────────────────────────────────────

comprehensive_small_sample_validation <- function(feature_name, abundance_df, metadata) {
  
  feature_col <- colnames(abundance_df)[1]
  
  if (!feature_name %in% abundance_df[[feature_col]]) {
    return(NULL)
  }
  
  feature_data <- abundance_df %>%
    filter(.data[[feature_col]] == feature_name) %>%
    select(-1) %>%
    pivot_longer(everything(), names_to = "sample_id", values_to = "abundance")
  
  plot_data <- metadata %>%
    inner_join(feature_data, by = "sample_id") %>%
    filter(!is.na(abundance))
  
  if (nrow(plot_data) < 4) {
    return(list(plot_data = plot_data, change_data = NULL, tests = NULL))
  }
  
  # Calculate baseline-to-6month changes
  change_data <- plot_data %>%
    filter(timepoint %in% c("Baseline", "6_months")) %>%
    select(participant, timepoint, turmeric_group, abundance) %>%
    pivot_wider(names_from = timepoint, values_from = abundance) %>%
    filter(!is.na(Baseline) & !is.na(`6_months`)) %>%
    mutate(
      change = `6_months` - Baseline,
      fold_change = `6_months` / (Baseline + 0.001),
      log_fold_change = log2(fold_change)
    )
  
  cat("    Paired analysis: n(participants with both timepoints)=", nrow(change_data), "\n")
  
  if (nrow(change_data) < MIN_PAIRED_PARTICIPANTS) {
    cat("    ⚠️  Insufficient paired data (need ≥4 participants)\n")
    return(list(plot_data = plot_data, change_data = change_data, tests = NULL))
  }
  
  treatment_changes <- change_data %>% filter(turmeric_group == "Treatment") %>% pull(change)
  control_changes <- change_data %>% filter(turmeric_group == "Control") %>% pull(change)
  
  cat("    n(Treatment)=", length(treatment_changes), 
      ", n(Control)=", length(control_changes), "\n")
  
  if (length(treatment_changes) < MIN_GROUP_N || length(control_changes) < MIN_GROUP_N) {
    cat("    ⚠️  Need ≥2 participants per group\n")
    return(list(plot_data = plot_data, change_data = change_data, tests = NULL))
  }
  
  test_results <- list()
  
  # TEST 1A: Exact Wilcoxon on CHANGES
  tryCatch({
    exact_wilcox_change <- wilcox.exact(treatment_changes, control_changes, exact = TRUE)
    
    test_results$exact_wilcoxon_paired <- data.frame(
      method = "Exact Wilcoxon (paired changes)",
      p_value = exact_wilcox_change$p.value,
      statistic = exact_wilcox_change$statistic,
      best_for = "Small n, most conservative",
      comparison = "Treatment Δ vs Control Δ"
    )
  }, error = function(e) {
    test_results$exact_wilcoxon_paired <<- NULL
  })
  
  # TEST 1B: Paired t-test on CHANGES
  tryCatch({
    if (length(treatment_changes) >= 3) {
      shapiro_treatment <- shapiro.test(treatment_changes)
      shapiro_control <- shapiro.test(control_changes)
      
      if (shapiro_treatment$p.value > 0.05 && shapiro_control$p.value > 0.05) {
        t_test_change <- t.test(treatment_changes, control_changes)
        
        test_results$t_test_paired <- data.frame(
          method = "Welch t-test (paired changes)",
          p_value = t_test_change$p.value,
          statistic = t_test_change$statistic,
          ci_lower = t_test_change$conf.int[1],
          ci_upper = t_test_change$conf.int[2],
          best_for = "Normal data, parametric power",
          comparison = "Treatment Δ vs Control Δ"
        )
      }
    }
  }, error = function(e) {
    test_results$t_test_paired <<- NULL
  })
  
  # TEST 2: Permutation test
  tryCatch({
    perm_data <- data.frame(
      change = c(treatment_changes, control_changes),
      group = factor(c(
        rep("Treatment", length(treatment_changes)),
        rep("Control", length(control_changes))
      ), levels = c("Control", "Treatment"))
    )
    
    perm_test <- coin::wilcox_test(
      change ~ group,
      data = perm_data,
      distribution = approximate(nresample = PERMUTATION_RESAMPLES)
    )
    
    test_results$permutation_paired <- data.frame(
      method = "Permutation (10k, paired)",
      p_value = coin::pvalue(perm_test),
      statistic = coin::statistic(perm_test),
      best_for = "Small samples, exact distribution",
      comparison = "Treatment Δ vs Control Δ"
    )
  }, error = function(e) {
    test_results$permutation_paired <<- NULL
  })
  
  # TEST 3: Bootstrap
  n_boot <- BOOTSTRAP_RESAMPLES
  
  boot_diffs <- replicate(n_boot, {
    mean(sample(treatment_changes, replace = TRUE)) - 
      mean(sample(control_changes, replace = TRUE))
  })
  
  obs_diff <- mean(treatment_changes) - mean(control_changes)
  boot_p <- 2 * min(mean(boot_diffs <= 0), mean(boot_diffs >= 0))
  boot_p <- min(1, boot_p)
  boot_ci <- quantile(boot_diffs, probs = c(0.025, 0.975), na.rm = TRUE)
  
  test_results$bootstrap_paired <- data.frame(
    method = "Bootstrap (10k, paired)",
    p_value = boot_p,
    statistic = obs_diff,
    ci_lower = boot_ci[1],
    ci_upper = boot_ci[2],
    best_for = "Effect size CI, distribution-free",
    comparison = "Treatment Δ vs Control Δ"
  )
  
  # TEST 4: Within-group paired tests
  treatment_data <- change_data %>% filter(turmeric_group == "Treatment")
  if (nrow(treatment_data) >= 3) {
    tryCatch({
      wilcox_treatment <- wilcox.test(
        treatment_data$Baseline,
        treatment_data$`6_months`,
        paired = TRUE,
        exact = TRUE
      )
      
      test_results$within_treatment <- data.frame(
        method = "Paired Wilcoxon (Treatment only)",
        p_value = wilcox_treatment$p.value,
        statistic = wilcox_treatment$statistic,
        best_for = "Within-subject change",
        comparison = "Treatment: Baseline vs 6m"
      )
    }, error = function(e) {
      test_results$within_treatment <<- NULL
    })
  }
  
  control_data <- change_data %>% filter(turmeric_group == "Control")
  if (nrow(control_data) >= 3) {
    tryCatch({
      wilcox_control <- wilcox.test(
        control_data$Baseline,
        control_data$`6_months`,
        paired = TRUE,
        exact = TRUE
      )
      
      test_results$within_control <- data.frame(
        method = "Paired Wilcoxon (Control only)",
        p_value = wilcox_control$p.value,
        statistic = wilcox_control$statistic,
        best_for = "Within-subject change",
        comparison = "Control: Baseline vs 6m"
      )
    }, error = function(e) {
      test_results$within_control <<- NULL
    })
  }
  
  # Cross-sectional at 6 months (SECONDARY)
  data_6m <- plot_data %>% filter(timepoint == "6_months")
  
  if (nrow(data_6m) >= 4) {
    treatment_6m <- data_6m %>% filter(turmeric_group == "Treatment") %>% pull(abundance)
    control_6m <- data_6m %>% filter(turmeric_group == "Control") %>% pull(abundance)
    
    if (length(treatment_6m) >= 2 && length(control_6m) >= 2) {
      tryCatch({
        wilcox_6m <- wilcox.exact(treatment_6m, control_6m, exact = TRUE)
        
        test_results$cross_sectional_6m <- data.frame(
          method = "Exact Wilcoxon (6m only)",
          p_value = wilcox_6m$p.value,
          statistic = wilcox_6m$statistic,
          best_for = "Cross-sectional comparison",
          comparison = "6m Treatment vs 6m Control"
        )
      }, error = function(e) {
        test_results$cross_sectional_6m <<- NULL
      })
    }
  }
  
  # EFFECT SIZES
  pooled_sd <- sqrt((var(treatment_changes, na.rm = TRUE) + var(control_changes, na.rm = TRUE)) / 2)
  cohens_d_change <- ifelse(
    is.finite(pooled_sd) && pooled_sd > 0,
    (mean(treatment_changes) - mean(control_changes)) / pooled_sd,
    NA_real_
  )
  
  if (!is.null(test_results$exact_wilcoxon_paired)) {
    U_stat <- test_results$exact_wilcoxon_paired$statistic
  } else {
    U_stat <- wilcox.test(treatment_changes, control_changes, exact = FALSE)$statistic
  }
  
  r_rb <- 1 - (2 * U_stat) / (length(treatment_changes) * length(control_changes))
  prob_superiority <- U_stat / (length(treatment_changes) * length(control_changes))
  
  test_results$effect_sizes <- data.frame(
    cohens_d_change = cohens_d_change,
    d_interpretation = case_when(
      abs(cohens_d_change) < 0.2 ~ "negligible",
      abs(cohens_d_change) < 0.5 ~ "small",
      abs(cohens_d_change) < 0.8 ~ "medium",
      TRUE ~ "large"
    ),
    rank_biserial = r_rb,
    rb_interpretation = case_when(
      abs(r_rb) < 0.1 ~ "negligible",
      abs(r_rb) < 0.3 ~ "small",
      abs(r_rb) < 0.5 ~ "medium",
      TRUE ~ "large"
    ),
    prob_superiority = prob_superiority,
    mean_change_treatment = mean(treatment_changes),
    mean_change_control = mean(control_changes),
    mean_change_diff = mean(treatment_changes) - mean(control_changes),
    median_change_treatment = median(treatment_changes),
    median_change_control = median(control_changes),
    median_change_diff = median(treatment_changes) - median(control_changes),
    n_treatment = length(treatment_changes),
    n_control = length(control_changes),
    n_total_paired = nrow(change_data),
    treatment_pct_increased = sum(treatment_changes > 0) / length(treatment_changes) * 100,
    control_pct_increased = sum(control_changes > 0) / length(control_changes) * 100,
    comparison = "Baseline-to-6m change: Treatment vs Control"
  )
  
  return(list(
    plot_data = plot_data,
    change_data = change_data,
    tests = test_results
  ))
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: Visualization
# ────────────────────────────────────────────────────────────────────────────────

create_enhanced_validation_plot <- function(validation_result, feature_name, profile_type) {
  
  if (is.null(validation_result) || is.null(validation_result$tests)) return(NULL)
  
  plot_data <- validation_result$plot_data
  change_data <- if (!is.null(validation_result$change_data)) validation_result$change_data else NULL
  feature_short <- extract_species_name(feature_name)
  tests <- validation_result$tests
  
  # Pre-specify exact Wilcoxon on paired changes as the primary p-value.
  primary_p <- if (!is.null(tests$exact_wilcoxon_paired)) {
    tests$exact_wilcoxon_paired$p_value
  } else {
    NA_real_
  }
  effect_info <- tests$effect_sizes
  
  # Create subtitle
  subtitle_parts <- c()
  if (!is.null(tests$exact_wilcoxon_paired)) {
    subtitle_parts <- c(subtitle_parts, 
                        glue("Exact Wilcoxon (Δ): p={format.pval(tests$exact_wilcoxon_paired$p_value, digits=3)}"))
  }
  if (!is.null(tests$permutation_paired)) {
    subtitle_parts <- c(subtitle_parts,
                        glue("Permutation (Δ): p={format.pval(tests$permutation_paired$p_value, digits=3)}"))
  }
  if (!is.null(effect_info)) {
    subtitle_parts <- c(subtitle_parts,
                        glue("Cohen's d={round(effect_info$cohens_d_change, 2)} ({effect_info$d_interpretation})"))
    subtitle_parts <- c(subtitle_parts,
                        glue("n(pairs)={effect_info$n_total_paired}"))
  }
  
  subtitle_text <- paste(subtitle_parts, collapse=" | ")
  
  # PLOT 1: Longitudinal trajectories
  p1 <- ggplot(plot_data, aes(x = timepoint, y = abundance, 
                              color = turmeric_group, group = participant)) +
    geom_line(alpha = 0.4, linewidth = 0.8) +
    geom_point(size = 2.5, alpha = 0.7) +
    stat_summary(aes(group = turmeric_group), fun = mean, geom = "line", 
                 linewidth = 1.5, linetype = "dashed") +
    stat_summary(aes(group = turmeric_group), fun = mean, geom = "point", 
                 size = 4, shape = 18) +
    scale_color_manual(values = c("Control" = "#3498DB", "Treatment" = "#E74C3C")) +
    labs(
      title = feature_short,
      subtitle = subtitle_text,
      x = "Timepoint",
      y = "Relative Abundance",
      color = "Group"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 7.5, lineheight = 1.1),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  # Add significance annotation
  if (!is.na(primary_p) && primary_p < PRIMARY_ALPHA) {
    sig_label <- case_when(
      primary_p < 0.001 ~ "***",
      primary_p < 0.01 ~ "**",
      primary_p < 0.05 ~ "*",
      TRUE ~ "ns"
    )
    
    p1 <- p1 +
      annotate("text",
               x = 2, y = max(plot_data$abundance, na.rm = TRUE) * 1.1,
               label = paste("Δ comparison:", sig_label),
               size = 4.5, fontface = "bold", color = "#E74C3C")
  }
  
  # PLOT 2: Change distributions (inset)
  if (!is.null(change_data) && nrow(change_data) >= 4) {
    p2 <- ggplot(change_data, aes(x = turmeric_group, y = change, fill = turmeric_group)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
      geom_jitter(width = 0.15, size = 2, alpha = 0.7, shape = 21, color = "black") +
      scale_fill_manual(values = c("Control" = "#3498DB", "Treatment" = "#E74C3C")) +
      labs(
        title = "Baseline → 6m Change",
        x = NULL,
        y = "Δ Abundance"
      ) +
      theme_bw(base_size = 9) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    # Combine plots
    combined_plot <- p1 + inset_element(p2, left = 0.6, bottom = 0.6, right = 0.98, top = 0.98)
    return(combined_plot)
  }
  
  return(p1)
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: Main Analysis Loop
# ────────────────────────────────────────────────────────────────────────────────

all_validation_results <- list()

for (profile in names(maaslin_data)) {
  
  if (is.null(maaslin_data[[profile]]) || is.null(abundance_data[[profile]])) {
    cat("\n⚠️  Skipping", profile, "- missing data\n")
    next
  }
  
  cat("\n", strrep("=", 80), "\n")
  cat(glue::glue("PROCESSING: {toupper(profile)}\n"))
  cat(strrep("=", 80), "\n")
  
  # Prepare metadata
  sample_ids <- colnames(abundance_data[[profile]])[-1]
  metadata <- prepare_metadata(sample_ids)

  min_non_zero <- ceiling(MIN_PREVALENCE_PROP * nrow(metadata))
  cat(glue::glue("Minimum non-zero prevalence count for screening: {min_non_zero}\n"))

  # Identify trends (screening only; inferential claims remain q-value based)
  trends <- identify_trends(
    maaslin_data[[profile]],
    profile,
    screen_p_threshold = SCREEN_P_THRESHOLD,
    min_non_zero = min_non_zero
  )
  
  # Validate top trends
  cat("\n═══════════════════════════════════════\n")
  cat("PAIRED ANALYSIS VALIDATION\n")
  cat("═══════════════════════════════════════\n")
  
  top_features <- trends$turmeric_trends %>%
    head(MAX_FEATURES_TO_VALIDATE) %>%
    pull(feature)
  
  cat("Validating top", length(top_features), "turmeric trends\n\n")
  
  validation_results <- list()
  
  for (i in seq_along(top_features)) {
    feature <- top_features[i]
    feature_short <- extract_species_name(feature)
    
    cat(glue::glue("  [{i}/{length(top_features)}] {feature_short}\n"))
    
    # Comprehensive validation
    validation <- comprehensive_small_sample_validation(
      feature,
      abundance_data[[profile]],
      metadata
    )
    
    if (!is.null(validation) && !is.null(validation$tests)) {
      validation_results[[feature]] <- validation
      
      # Create plot
      val_plot <- create_enhanced_validation_plot(validation, feature, profile)
      
      if (!is.null(val_plot)) {
        safe_name <- str_replace_all(feature_short, "[^[:alnum:]_]", "_")
        safe_name <- str_trunc(safe_name, width = 50)
        
        save_and_upload_to_github(
          object = val_plot,
          object_type = "plot",
          github_path = file.path(figure_output_diff, "validation_paired_analysis",
                                  glue::glue("{profile}_val_{i}_{safe_name}.png")),
          commit_message = glue::glue("Add paired analysis validation plot for {feature_short}"),
          owner = owner, repo = repo, branch = branch
        )
      }
    }
  }
  
  # Compile validation summary
  cat("\n═══════════════════════════════════════\n")
  cat("VALIDATION SUMMARY\n")
  cat("═══════════════════════════════════════\n")
  
  if (length(validation_results) > 0) {
    validation_summary <- bind_rows(lapply(names(validation_results), function(feat) {
      val <- validation_results[[feat]]
      tests <- val$tests
      
      # Extract p-values
      exact_wilcox_p <- if (!is.null(tests$exact_wilcoxon_paired)) {
        as.numeric(tests$exact_wilcoxon_paired$p_value)
      } else NA
      
      permutation_p <- if (!is.null(tests$permutation_paired)) {
        as.numeric(tests$permutation_paired$p_value)
      } else NA
      
      bootstrap_p <- if (!is.null(tests$bootstrap_paired)) {
        as.numeric(tests$bootstrap_paired$p_value)
      } else NA
      
      within_treatment_p <- if (!is.null(tests$within_treatment)) {
        as.numeric(tests$within_treatment$p_value)
      } else NA
      
      within_control_p <- if (!is.null(tests$within_control)) {
        as.numeric(tests$within_control$p_value)
      } else NA
      
      cross_sectional_p <- if (!is.null(tests$cross_sectional_6m)) {
        as.numeric(tests$cross_sectional_6m$p_value)
      } else NA
      
      # Effect sizes
      effect_sizes <- tests$effect_sizes
      cohens_d <- if (!is.null(effect_sizes)) effect_sizes$cohens_d_change else NA
      rank_biserial <- if (!is.null(effect_sizes)) effect_sizes$rank_biserial else NA
      effect_interp <- if (!is.null(effect_sizes)) effect_sizes$d_interpretation else NA
      mean_change_diff <- if (!is.null(effect_sizes)) effect_sizes$mean_change_diff else NA
      mean_change_treatment <- if (!is.null(effect_sizes)) effect_sizes$mean_change_treatment else NA
      mean_change_control <- if (!is.null(effect_sizes)) effect_sizes$mean_change_control else NA
      n_treatment <- if (!is.null(effect_sizes)) effect_sizes$n_treatment else NA
      n_control <- if (!is.null(effect_sizes)) effect_sizes$n_control else NA
      pct_increased_treatment <- if (!is.null(effect_sizes)) effect_sizes$treatment_pct_increased else NA
      pct_increased_control <- if (!is.null(effect_sizes)) effect_sizes$control_pct_increased else NA

      data.frame(
        feature = feat,
        feature_short = extract_species_name(feat),
        primary_paired_p = exact_wilcox_p,
        permutation_paired_p = permutation_p,
        bootstrap_paired_p = bootstrap_p,
        within_treatment_p = within_treatment_p,
        within_control_p = within_control_p,
        cross_sectional_6m_p = cross_sectional_p,
        cohens_d = cohens_d,
        rank_biserial = rank_biserial,
        effect_interpretation = effect_interp,
        mean_change_diff = mean_change_diff,
        mean_change_treatment = mean_change_treatment,
        mean_change_control = mean_change_control,
        n_treatment = n_treatment,
        n_control = n_control,
        pct_increased_treatment = pct_increased_treatment,
        pct_increased_control = pct_increased_control,
        stringsAsFactors = FALSE
      )
    }))
    
    # Add MaAsLin3 info
    validation_summary <- validation_summary %>%
      left_join(
        trends$turmeric_trends %>%
          select(feature, maaslin_p = pval_individual, maaslin_coef = coef, maaslin_qval = qval_individual),
        by = "feature"
      ) %>%
      mutate(
        primary_paired_q = p.adjust(primary_paired_p, method = "BH"),
        exploratory_support = ifelse(
          !is.na(permutation_paired_p) & !is.na(bootstrap_paired_p),
          permutation_paired_p < PRIMARY_ALPHA & bootstrap_paired_p < PRIMARY_ALPHA,
          FALSE
        )
      ) %>%
      arrange(primary_paired_q, primary_paired_p)
    
    cat("\nTop 20 validation results (sorted by primary paired q-value):\n")
    print(validation_summary %>%
            select(feature_short, primary_paired_p, primary_paired_q, permutation_paired_p,
                   bootstrap_paired_p, cohens_d, effect_interpretation,
                   within_treatment_p, within_control_p) %>%
            head(20))
    
    # Save results
    save_and_upload_to_github(
      object = validation_summary,
      object_type = "data",
      github_path = file.path(output_path_diff,
                              glue::glue("{profile}_validation_paired_analysis.csv")),
      commit_message = glue::glue("Add {profile} paired analysis validation results"),
      owner = owner, repo = repo, branch = branch
    )
    
    # Report findings
    sig_q_005 <- sum(validation_summary$primary_paired_q < PRIMARY_ALPHA, na.rm = TRUE)
    sig_q_025 <- sum(validation_summary$primary_paired_q < EXPLORATORY_Q_ALPHA, na.rm = TRUE)
    sig_p_005 <- sum(validation_summary$primary_paired_p < PRIMARY_ALPHA, na.rm = TRUE)
    
    cat(glue::glue("\nPAIRED ANALYSIS FINDINGS:\n"))
    cat(glue::glue("  Primary test nominal p < {PRIMARY_ALPHA}: {sig_p_005}\n"))
    cat(glue::glue("  Primary test BH q < {PRIMARY_ALPHA}: {sig_q_005}\n"))
    cat(glue::glue("  Primary test BH q < {EXPLORATORY_Q_ALPHA} (exploratory): {sig_q_025}\n\n"))
    
    # Show strongest trends
    if (sig_q_025 > 0) {
      cat("Strongest validated changes (baseline→6m):\n")
      strong_trends <- validation_summary %>%
        filter(primary_paired_q < EXPLORATORY_Q_ALPHA) %>%
        arrange(primary_paired_q, primary_paired_p) %>%
        mutate(
          interpretation = case_when(
            within_treatment_p < 0.05 & within_control_p >= 0.10 ~ "Treatment-specific",
            within_treatment_p < 0.05 & within_control_p < 0.05 ~ "Both changed",
            within_treatment_p >= 0.10 & within_control_p < 0.05 ~ "Control-specific",
            TRUE ~ "Between-group diff"
          )
        ) %>%
        select(feature_short, primary_paired_p, primary_paired_q, cohens_d, effect_interpretation, 
               mean_change_treatment, mean_change_control, interpretation)
      print(strong_trends)
    }
    
    # Treatment-specific effects
    treatment_specific <- validation_summary %>%
      filter(primary_paired_q < EXPLORATORY_Q_ALPHA,
             within_treatment_p < 0.10,
             within_control_p >= 0.10)
    
    if (nrow(treatment_specific) > 0) {
      cat(glue::glue("\n🎯 TREATMENT-SPECIFIC EFFECTS (n={nrow(treatment_specific)}):\n"))
      cat("(Treatment group changed significantly, Control did not)\n")
      print(treatment_specific %>%
              select(feature_short, primary_paired_p, primary_paired_q, within_treatment_p, 
                     within_control_p, mean_change_treatment))
    }
    
  } else {
    cat("⚠️  No validation results to summarize\n")
    validation_summary <- NULL
  }
  
  all_validation_results[[profile]] <- list(
    trends = trends,
    validation = validation_results,
    summary = validation_summary
  )
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: Final Summary
# ────────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 80), "\n")
cat("PAIRED ANALYSIS VALIDATION COMPLETE!\n")
cat(strrep("=", 80), "\n")

for (profile in names(all_validation_results)) {
  res <- all_validation_results[[profile]]
  cat(glue::glue("\n{toupper(profile)}:\n"))
  cat(glue::glue("  MaAsLin3 screened trends (p<{SCREEN_P_THRESHOLD}): {nrow(res$trends$turmeric_trends)}\n"))
  
  if (!is.null(res$summary)) {
    n_validated_p <- sum(res$summary$primary_paired_p < PRIMARY_ALPHA, na.rm = TRUE)
    n_validated_q_confirm <- sum(res$summary$primary_paired_q < PRIMARY_ALPHA, na.rm = TRUE)
    n_validated_q_explore <- sum(res$summary$primary_paired_q < EXPLORATORY_Q_ALPHA, na.rm = TRUE)
    cat(glue::glue("  Primary paired p<{PRIMARY_ALPHA}: {n_validated_p}\n"))
    cat(glue::glue("  Primary paired BH q<{PRIMARY_ALPHA}: {n_validated_q_confirm}\n"))
    cat(glue::glue("  Primary paired BH q<{EXPLORATORY_Q_ALPHA} (exploratory): {n_validated_q_explore}\n"))
  }
}

cat("\nMETHODS USED (paired within-subject analysis):\n")
cat("  1. Paired baseline→6m changes calculated for each participant\n")
cat("  2. Exact Wilcoxon test on changes (primary endpoint)\n")
cat("  3. Permutation test (secondary sensitivity analysis)\n")
cat("  4. Bootstrap null test + 95% CI on mean change difference\n")
cat("  5. Within-group paired tests (Treatment & Control separately)\n")
cat("  6. Effect sizes: Cohen's d on changes, rank-biserial\n")
cat("  7. BH correction applied to primary paired p-values across validated features\n")
cat("\nFILES GENERATED:\n")
cat("  • Paired analysis validation results CSV\n")
cat("  • Longitudinal trajectory plots with change insets\n")
cat("  • Effect size estimates and interpretations\n")
cat("\n")







