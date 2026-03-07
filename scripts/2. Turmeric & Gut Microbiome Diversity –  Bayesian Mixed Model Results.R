# ────────────────────────────────────────────────────────────────────────────────
# Script: 06_statistical_analysis_alpha_diversity_BAYESIAN.R
# Description: Bayesian statistical analysis of alpha diversity metrics
# Author: Liam Walsh
# Date: 2025
# ────────────────────────────────────────────────────────────────────────────────

# ────────────────────────────────────────────────────────────────────────────────
# Section: Setup and Configuration
# ────────────────────────────────────────────────────────────────────────────────

preferred_lib <- Sys.getenv("NTU_R_LIB")
if (!nzchar(preferred_lib)) preferred_lib <- "G:/R/win-library/4.4"
if (!dir.exists(preferred_lib)) {
  preferred_lib <- normalizePath(
    file.path(Sys.getenv("USERPROFILE"), "Documents", "R", "win-library", "4.4"),
    winslash = "/", mustWork = FALSE
  )
}
dir.create(preferred_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(preferred_lib, .libPaths()))
set.seed(123)

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cran.rstudio.com/")
}

suppressPackageStartupMessages({
  pacman::p_load(
    readxl, readr, dplyr, tidyr, stringr, tidyverse, httr,
    data.table, ggplot2, viridis, gh, base64enc,
    brms, bayesplot, tidybayes, posterior, glue, broom.mixed
  )
})

# GitHub configuration (same as original)
owner <- "LiamHWalsh"
repo <- "NTU_Teagasc_Turmeric_study"
branch <- "main"

pipeline_paths <- list(
  raw_data = "results/Raw",
  preprocessing = "results/Processed/01_preprocessing",
  taxonomic = "results/Processed/02_taxonomic_profiling",
  functional = "results/Processed/03_functional_profiling",
  diversity = "results/Processed/04_diversity_analysis",
  statistics = "results/Processed/05_statistical_analysis",
  reports = "results/Processed/06_reports"
)

figure_paths <- list(
  preprocessing = "figures/01_preprocessing",
  taxonomic = "figures/02_taxonomic_profiling",
  functional = "figures/03_functional_profiling",
  diversity = "figures/04_diversity_analysis",
  statistics = "figures/05_statistical_analysis",
  reports = "figures/06_reports"
)

input_path <- pipeline_paths$diversity
output_path <- pipeline_paths$statistics
figure_output_path <- figure_paths$statistics

token <- Sys.getenv("GITHUB_PAT")
if (token == "") token <- Sys.getenv("GITHUB_TOKEN")
enable_github_upload <- identical(Sys.getenv("ENABLE_GITHUB_UPLOAD"), "1")
allow_github_write <- nzchar(token) && enable_github_upload
if (nzchar(token)) {
  message(glue::glue("GitHub token detected (length: {nchar(token)} characters)"))
} else {
  message("No GitHub token detected. Running in read-only / local-first mode.")
}
if (!allow_github_write) {
  message("GitHub upload disabled. Set ENABLE_GITHUB_UPLOAD=1 with a valid token to upload outputs.")
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: GitHub Functions (same as original)
# ────────────────────────────────────────────────────────────────────────────────

load_github_file <- function(file_url, token) {
  response <- if (nzchar(token)) {
    GET(file_url, authenticate(token, ""))
  } else {
    GET(file_url)
  }
  if (status_code(response) != 200) {
    stop(glue::glue("Failed to load file: {file_url}"))
  }
  return(read_csv(I(rawToChar(response$content)), show_col_types = FALSE))
}

save_and_upload_to_github <- function(object = NULL, object_type = c("data", "plot", "file"), 
                                      local_path = NULL, github_path, commit_message,
                                      owner, repo, branch = "main") {
  object_type <- match.arg(object_type)
  
  if (is.null(local_path)) {
    file_ext <- switch(object_type, data = ".csv", plot = ".png",
                       stop("For 'file' type, 'local_path' must be provided."))
    local_path <- tempfile(fileext = file_ext)
  }
  
  if (object_type == "data") {
    if (is.null(object)) stop("❌ No data object provided.")
    data.table::fwrite(object, local_path, row.names = FALSE)
  } else if (object_type == "plot") {
    if (is.null(object)) stop("❌ No plot object provided.")
    ggplot2::ggsave(local_path, plot = object, width = 10, height = 6, dpi = 300)
  } else if (object_type == "file") {
    if (!file.exists(local_path)) stop(glue::glue("❌ File does not exist: {local_path}"))
  }
  
  if (!allow_github_write) {
    local_output_path <- normalizePath(github_path, winslash = "/", mustWork = FALSE)
    dir.create(dirname(local_output_path), recursive = TRUE, showWarnings = FALSE)
    if (normalizePath(local_path, winslash = "/", mustWork = FALSE) != local_output_path) {
      file.copy(local_path, local_output_path, overwrite = TRUE)
      if (object_type %in% c("data", "plot")) unlink(local_path)
    }
    message(glue::glue("Saved locally to: {local_output_path}"))
    return(invisible(local_output_path))
  }
  
  encoded_content <- base64enc::base64encode(local_path)
  gh_token <- Sys.getenv("GITHUB_PAT")
  if (gh_token == "") gh_token <- Sys.getenv("GITHUB_TOKEN")
  
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
# Section: Bayesian LMM Analysis Function
# ────────────────────────────────────────────────────────────────────────────────

fit_bayesian_lmm <- function(data, response_var, group_var = "turmeric_status") {
  
  # Data structure assessment
  design_matrix <- table(data$timepoint, data[[group_var]])
  n_timepoints <- length(unique(data$timepoint))
  has_full_design <- all(design_matrix > 0)
  
  participant_counts <- data %>%
    group_by(block_number) %>%
    summarise(n_obs = n(), .groups = "drop")
  
  n_participants <- nrow(participant_counts)
  n_observations <- nrow(data)
  participants_with_repeats <- sum(participant_counts$n_obs > 1)
  
  cat("\n═══════════════════════════════════════\n")
  cat("Bayesian Mixed Model for:", response_var, "\n")
  cat("═══════════════════════════════════════\n")
  cat("\nData Structure:\n")
  cat("  Total observations:", n_observations, "\n")
  cat("  Unique participants:", n_participants, "\n")
  cat("  Participants with 2+ samples:", participants_with_repeats, "\n")
  
  # Ensure proper factor levels
  data <- data %>%
    mutate(
      !!sym(group_var) := factor(.data[[group_var]]),
      timepoint = factor(timepoint, levels = c("Baseline", "2_weeks", "6_months"))
    )
  
  # Define priors
  priors <- c(
    prior(normal(0, 5), class = "b"),           # Fixed effects
    prior(student_t(3, 0, 5), class = "sd"),    # Random intercept SD
    prior(student_t(3, 0, 10), class = "sigma") # Residual SD
  )
  
  # Build formula
  can_use_random_effects <- (n_observations > n_participants + 5) && 
    (participants_with_repeats >= 3)
  
  if (!can_use_random_effects) {
    cat("⚠️  Using BAYESIAN LINEAR MODEL (insufficient data for random effects)\n\n")
    
    if (!has_full_design || n_timepoints == 1) {
      formula_str <- paste0(response_var, " ~ ", group_var, 
                            ifelse(n_timepoints > 1, " + timepoint", ""))
    } else {
      formula_str <- paste0(response_var, " ~ ", group_var, " * timepoint")
    }
    
    model_type <- "bayesian_lm"
    priors <- priors[priors$class != "sd"]  # Remove random effect prior
    
  } else {
    cat("✓ Using BAYESIAN MIXED EFFECTS MODEL\n\n")
    
    if (!has_full_design) {
      formula_str <- paste0(response_var, " ~ ", group_var, 
                            " + timepoint + (1|block_number)")
    } else {
      formula_str <- paste0(response_var, " ~ ", group_var, 
                            " * timepoint + (1|block_number)")
    }
    
    model_type <- "bayesian_lmm"
  }
  
  cat("Formula:", formula_str, "\n")
  cat("Fitting model with brms...\n\n")
  
  # Fit Bayesian model
  b_model <- brm(
    formula = as.formula(formula_str),
    data = data,
    family = gaussian(),
    prior = priors,
    chains = 4,
    iter = 4000,
    warmup = 2000,
    cores = 4,
    seed = 123,
    control = list(adapt_delta = 0.95, max_treedepth = 15),
    silent = 2,
    refresh = 0
  )
  
  # Extract posterior summaries
  model_summary <- summary(b_model)
  
  # Test specific hypotheses
  hypotheses <- list()
  
  # Get coefficient names
  coef_names <- rownames(model_summary$fixed)
  
  # Treatment effect hypothesis
  treatment_coef <- grep(paste0(group_var, "Taken"), coef_names, value = TRUE)
  if (length(treatment_coef) > 0) {
    hypotheses$treatment_decrease <- hypothesis(b_model, 
                                                paste0(treatment_coef, " < 0"))
    hypotheses$treatment_increase <- hypothesis(b_model, 
                                                paste0(treatment_coef, " > 0"))
  }
  
  # Interaction hypotheses if present
  interaction_coefs <- grep(":", coef_names, value = TRUE)
  if (length(interaction_coefs) > 0) {
    for (coef in interaction_coefs) {
      hyp_name <- gsub(":", "_x_", coef)
      hypotheses[[paste0(hyp_name, "_nonzero")]] <- hypothesis(b_model, 
                                                               paste0(coef, " != 0"))
    }
  }
  
  # Convergence diagnostics
  convergence <- list(
    rhat_max = max(rhat(b_model)),
    rhat_ok = all(rhat(b_model) < 1.01),
    neff_min = min(neff_ratio(b_model)),
    neff_ok = all(neff_ratio(b_model) > 0.1)
  )
  
  cat("\nConvergence Diagnostics:\n")
  cat("  Max Rhat:", round(convergence$rhat_max, 4), 
      ifelse(convergence$rhat_ok, "✓", "✗"), "\n")
  cat("  Min ESS ratio:", round(convergence$neff_min, 4), 
      ifelse(convergence$neff_ok, "✓", "✗"), "\n\n")
  
  return(list(
    model = b_model,
    model_type = model_type,
    formula = formula_str,
    summary = model_summary,
    hypotheses = hypotheses,
    convergence = convergence
  ))
}



# ────────────────────────────────────────────────────────────────────────────────
# Section: Visualization Functions
# ────────────────────────────────────────────────────────────────────────────────

make_alpha_plot <- function(data, y_var, profile_type, y_label) {
  p <- ggplot(data, aes(x = timepoint, y = .data[[y_var]], 
                        fill = turmeric_status, color = turmeric_status)) +
    geom_line(aes(group = block_number), alpha = 0.3, linewidth = 0.5) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(0.8)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                size = 3, alpha = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +
    labs(
      title = paste(profile_type, "-", y_label),#, "across timepoints"
      x = "Timepoint",
      y = y_label,
      fill = "Turmeric Status",
      color = "Turmeric Status"
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.position = "top"
    )
  
  return(p)
}

make_posterior_plot <- function(bayesian_result, response_var, profile_type) {
  # Extract posterior draws
  posterior_draws <- as_draws_df(bayesian_result$model)
  
  # Create trace plots
  trace_plot <- mcmc_trace(posterior_draws, 
                           regex_pars = "^b_",
                           facet_args = list(ncol = 2)) +
    ggtitle(paste(profile_type, response_var, "- MCMC Trace Plots")) +
    theme_minimal()
  
  return(trace_plot)
}

make_posterior_intervals_plot <- function(bayesian_result, response_var, profile_type) {
  # Posterior intervals
  intervals_plot <- mcmc_intervals(bayesian_result$model, 
                                   regex_pars = "^b_",
                                   prob = 0.5, prob_outer = 0.95) +
    ggtitle(paste(profile_type, response_var, "- Posterior Intervals (50% & 95%)")) +
    theme_minimal()
  
  return(intervals_plot)
}




# ────────────────────────────────────────────────────────────────────────────────
# Section: Load Data from GitHub
# ────────────────────────────────────────────────────────────────────────────────

github_base <- paste0("https://raw.githubusercontent.com/", owner, "/", repo, "/", branch)

cat("\n═══════════════════════════════════════\n")
cat("Loading Data from GitHub\n")
cat("═══════════════════════════════════════\n")

data <- list()
files_to_load <- c(
  "04_diversity_analysis/species_alpha_diversity.csv",
  "04_diversity_analysis/functional_alpha_diversity.csv",
  "03_functional_profiling/functional_profile_filtered.csv",
  "02_taxonomic_profiling/species_profile_filtered.csv"
)

for (file_name in files_to_load) {
  file_url <- file.path(github_base, "results","Processed", file_name)
  clean_name <- tools::file_path_sans_ext(basename(file_name))
  data[[clean_name]] <- load_github_file(file_url, token)
  cat("✅", file_name, "\n")
}

# Turmeric status lookup (same as original)
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

# ────────────────────────────────────────────────────────────────────────────────
# Section: Main Bayesian Analysis Loop
# ────────────────────────────────────────────────────────────────────────────────

cat("\nSTARTING BAYESIAN STATISTICAL ANALYSIS\n")

all_results <- list()
plot=list()
for (alpha_type in c("species_alpha_diversity", "functional_alpha_diversity")) {
  
  profile_type <- ifelse(grepl("species", alpha_type), "Species", "Functional")
  
  cat(glue::glue("\nANALYZING: {profile_type} Alpha Diversity\n"))
  
  # Prepare metadata (same as original)
  ids <- unique(data[[alpha_type]]$sample_id)
  split_ids <- str_split_fixed(ids, "_", 3)
  
  sample_ids <- data.frame(
    full_id = ids,
    block_code = split_ids[,2],
    block_number = as.integer(gsub("[^0-9]", "", split_ids[,2])),
    block_letter = gsub("[0-9]", "", split_ids[,2]),
    stringsAsFactors = FALSE
  )
  
  sample_ids <- merge(sample_ids, turmeric_status,
                      by.x = "block_number", by.y = "ID", all.x = TRUE)
  
  sample_ids$timepoint <- "Baseline"
  sample_ids$timepoint[sample_ids$block_letter %in% c("a","A")] <- "Baseline"
  sample_ids$timepoint[sample_ids$block_letter %in% c("b","B")] <- "2_weeks"
  sample_ids$timepoint[sample_ids$block_letter %in% c("c","C")] <- "6_months"
  
  sample_ids$turmeric_group <- "Control"
  sample_ids$turmeric_group[
    grepl("Turmeric 6 months", sample_ids$status) & 
      sample_ids$block_letter %in% c("c", "C")
  ] <- "Treatment"
  
  sample_ids$turmeric_group[
    sample_ids$status == "Turmeric 3 months" & 
      sample_ids$block_letter == "C"
  ] <- "Stopped"
  
  for (i in 1:nrow(sample_ids)) {
    if (sample_ids$timepoint[i] %in% c("Baseline", "2_weeks")) {
      participant_id <- sample_ids$block_number[i]
      has_treatment <- any(
        sample_ids$block_number == participant_id & 
          sample_ids$turmeric_group == "Treatment"
      )
      if (has_treatment) sample_ids$turmeric_group[i] <- "Treatment"
    }
  }
  
  sample_ids$block_letter <- str_to_upper(sample_ids$block_letter)
  
  analysis_data <- merge(sample_ids, data[[alpha_type]],
                         by.x = "full_id", by.y = "sample_id", all.x = TRUE)
  
  filtered_data <- analysis_data %>%
    filter(turmeric_group != "Stopped") %>%
    mutate(
      timepoint = factor(timepoint, levels = c("Baseline", "2_weeks", "6_months")),
      turmeric_group = factor(turmeric_group, levels = c("Control", "Treatment"))
    )
  
  filtered_data$turmeric_status <- "Not taken"
  filtered_data$turmeric_status[
    filtered_data$block_letter == "C" & 
      filtered_data$status == "Turmeric 6 months"
  ] <- "Taken"
  filtered_data$turmeric_status <- factor(filtered_data$turmeric_status, 
                                          levels = c("Not taken", "Taken"))
  
  # Run Bayesian analyses
  metrics <- c("Richness", "Shannon", "Simpson")
  results <- list()
  
  for (metric in metrics) {
    cat("\n", strrep("-", 80), "\n")
    cat(glue::glue("Analyzing: {metric}\n"))
    cat(strrep("-", 80), "\n")
    
    # Fit Bayesian model
    bayesian_result <- fit_bayesian_lmm(filtered_data, metric, 
                                        group_var = "turmeric_status")
    
    results[[metric]] <- bayesian_result
    
    # Create and upload main plot
    plot[[alpha_type]][[metric]] <- make_alpha_plot(filtered_data, metric, profile_type, 
                            paste(metric, "Diversity"))
    
    plot_filename <- tolower(glue::glue("{profile_type}_{metric}_alpha_diversity_bayesian.png"))
    
    save_and_upload_to_github(
      object =  plot[[alpha_type]][[metric]],
      object_type = "plot",
      github_path = file.path(figure_output_path, plot_filename),
      commit_message = glue::glue("Add Bayesian {profile_type} {metric} diversity plot"),
      owner = owner, repo = repo, branch = branch
    )
    
    # Create and upload posterior plots
    trace_plot <- make_posterior_plot(bayesian_result, metric, profile_type)
    trace_filename <- tolower(glue::glue("{profile_type}_{metric}_trace_plots.png"))
    
    # save_and_upload_to_github(
    #   object = trace_plot,
    #   object_type = "plot",
    #   github_path = file.path(figure_output_path, trace_filename),
    #   commit_message = glue::glue("Add {profile_type} {metric} trace plots"),
    #   owner = owner, repo = repo, branch = branch
    # )
    
    # intervals_plot <- make_posterior_intervals_plot(bayesian_result, metric, profile_type)
    # intervals_filename <- tolower(glue::glue("{profile_type}_{metric}_posterior_intervals.png"))
    # 
    # save_and_upload_to_github(
    #   object = intervals_plot,
    #   object_type = "plot",
    #   github_path = file.path(figure_output_path, intervals_filename),
    #   commit_message = glue::glue("Add {profile_type} {metric} posterior intervals"),
    #   owner = owner, repo = repo, branch = branch
    # )
    
    # Save results
    # Fixed effects summary
    fixed_effects <- as.data.frame(bayesian_result$summary$fixed)
    fixed_effects$parameter <- rownames(fixed_effects)
    fixed_effects$profile <- profile_type
    fixed_effects$metric <- metric
    fixed_effects$model_type <- bayesian_result$model_type
    
    save_and_upload_to_github(
      object = fixed_effects,
      object_type = "data",
      github_path = file.path(output_path, 
                              tolower(glue::glue("{profile_type}_{metric}_bayesian_fixed_effects.csv"))),
      commit_message = glue::glue("Add Bayesian {profile_type} {metric} fixed effects"),
      owner = owner, repo = repo, branch = branch
    )
    
    # Hypothesis tests
    if (length(bayesian_result$hypotheses) > 0) {
      hyp_df <- bind_rows(lapply(names(bayesian_result$hypotheses), function(h_name) {
        h <- bayesian_result$hypotheses[[h_name]]
        data.frame(
          profile = profile_type,
          metric = metric,
          hypothesis = h_name,
          estimate = h$hypothesis$Estimate,
          est_error = h$hypothesis$Est.Error,
          CI_lower = h$hypothesis$CI.Lower,
          CI_upper = h$hypothesis$CI.Upper,
          evid_ratio = h$hypothesis$Evid.Ratio,
          post_prob = h$hypothesis$Post.Prob
        )
      }))
      
      save_and_upload_to_github(
        object = hyp_df,
        object_type = "data",
        github_path = file.path(output_path,
                                tolower(glue::glue("{profile_type}_{metric}_bayesian_hypotheses.csv"))),
        commit_message = glue::glue("Add Bayesian {profile_type} {metric} hypothesis tests"),
        owner = owner, repo = repo, branch = branch
      )
    }
  }
  
  all_results[[alpha_type]] <- results
}

library(ggpubr)

alpha_diversity_figure <- ggarrange(
  plot[["species_alpha_diversity"]][["Richness"]],
  plot[["species_alpha_diversity"]][["Shannon"]],
  plot[["species_alpha_diversity"]][["Simpson"]],
  plot[["functional_alpha_diversity"]][["Richness"]],
  plot[["functional_alpha_diversity"]][["Shannon"]],
  plot[["functional_alpha_diversity"]][["Simpson"]],
  ncol = 3,
  nrow = 2,
  labels = c("A", "B", "C", "D", "E", "F"),
  common.legend = TRUE,
  legend = "bottom"
)


# Force white background everywhere (including legend)
alpha_diversity_figure <- ggpar(
  alpha_diversity_figure,
  background = "white"
)
library(ggpubr)

# Create a temporary file
fig_path <- tempfile(fileext = ".jpeg")

# Save the figure
ggsave(
  filename = fig_path,
  plot = alpha_diversity_figure,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

save_and_upload_to_github(
  object_type = "file",
  local_path = fig_path,
  github_path = file.path("figures", "manuscript_figures", "figure_1.jpeg"),
  commit_message = "Add Figure 1 for the manuscript",
  owner = owner,
  repo = repo,
  branch = branch
)

# ────────────────────────────────────────────────────────────────────────────────
# Section: Summary Report
# ────────────────────────────────────────────────────────────────────────────────

summary_text <- glue::glue("
═══════════════════════════════════════════════════════════════════════════════
BAYESIAN STATISTICAL ANALYSIS REPORT: ALPHA DIVERSITY
═══════════════════════════════════════════════════════════════════════════════

Script: 06_statistical_analysis_alpha_diversity_BAYESIAN.R
Generated: {Sys.time()}

ANALYSIS OVERVIEW:
  Method: Bayesian Linear Mixed Models (brms)
  Chains: 4
  Iterations: 4000 (2000 warmup)
  Priors: Weakly informative

═══════════════════════════════════════════════════════════════════════════════
BAYESIAN APPROACH ADVANTAGES
═══════════════════════════════════════════════════════════════════════════════

1. Full posterior distributions for all parameters
2. Natural quantification of uncertainty
3. Direct probability statements about effects
4. Robust to small sample sizes
5. Hypothesis testing via posterior probabilities

═══════════════════════════════════════════════════════════════════════════════
INTERPRETATION GUIDE
═══════════════════════════════════════════════════════════════════════════════

Fixed Effects:
  • Estimate: Posterior mean
  • Est.Error: Posterior standard deviation
  • CI.Lower/Upper: 95% credible intervals
  • Rhat: Convergence diagnostic (should be < 1.01)
  • ESS: Effective sample size (should be > 400)

Hypothesis Tests:
  • Post.Prob: Posterior probability hypothesis is true
  • Evid.Ratio: Evidence ratio (odds in favor)
  • > 0.95: Strong evidence
  • 0.90-0.95: Moderate evidence
  • < 0.90: Weak evidence

═══════════════════════════════════════════════════════════════════════════════
CONVERGENCE STATUS
═══════════════════════════════════════════════════════════════════════════════
")

for (alpha_type in names(all_results)) {
  profile <- ifelse(grepl("species", alpha_type), "Species", "Functional")
  summary_text <- paste0(summary_text, "\n", profile, " Profile:\n")
  
  for (metric in names(all_results[[alpha_type]])) {
    conv <- all_results[[alpha_type]][[metric]]$convergence
    summary_text <- paste0(summary_text, 
                           glue::glue("  {metric}: Rhat={round(conv$rhat_max, 4)} ",
                                      "{ifelse(conv$rhat_ok, '✓', '✗')}, ",
                                      "ESS ratio={round(conv$neff_min, 4)} ",
                                      "{ifelse(conv$neff_ok, '✓', '✗')}\n"))
  }
}

summary_text <- paste0(summary_text, "\n\n✅ BAYESIAN ANALYSIS COMPLETE\n")

temp_summary <- tempfile(fileext = ".txt")
cat(summary_text, file = temp_summary)








save_and_upload_to_github(
  object_type = "file",
  local_path = temp_summary,
  github_path = file.path(pipeline_paths$reports, 
                          "06_bayesian_analysis_alpha_diversity_REPORT.txt"),
  commit_message = "Add Bayesian analysis alpha diversity report",
  owner = owner, repo = repo, branch = branch
)

unlink(temp_summary)

cat("\n✅ BAYESIAN ANALYSIS COMPLETE\n")
cat(glue::glue("\n🔗 View results: https://github.com/{owner}/{repo}/tree/{branch}/{output_path}\n"))
cat(glue::glue("🔗 View figures: https://github.com/{owner}/{repo}/tree/{branch}/{figure_output_path}\n"))
