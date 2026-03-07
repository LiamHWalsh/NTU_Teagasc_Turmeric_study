# ────────────────────────────────────────────────────────────────────────────────
# Script: 08_differential_abundance_analysis.R
# Description: Differential abundance analysis using DESeq2 and other methods
#              comparing turmeric treatment effects across timepoints
# Input: Species and functional profiles from earlier steps
# Output: Statistical results to 05_statistical_analysis/differential_abundance
#         Figures to figures/05_statistical_analysis/differential_abundance
# Author: Liam Walsh
# Date: 2025
# ────────────────────────────────────────────────────────────────────────────────

# ────────────────────────────────────────────────────────────────────────────────
# Section: Setup and Configuration
# ────────────────────────────────────────────────────────────────────────────────

set.seed(123)

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cran.rstudio.com/")
}

suppressPackageStartupMessages({
  pacman::p_load(
    readxl, readr, dplyr, tidyr, stringr, tidyverse, httr,
    data.table, ggplot2, viridis, gh, base64enc,
    DESeq2, edgeR, limma, ggrepel, patchwork, glue, pheatmap
  )
})

# ────────────────────────────────────────────────────────────────────────────────
# Section: GitHub Configuration
# ────────────────────────────────────────────────────────────────────────────────

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

# Organized output paths
output_path_beta <- file.path(pipeline_paths$statistics, "beta_diversity")
output_path_diff <- file.path(pipeline_paths$statistics, "differential_abundance")
output_path_alpha <- file.path(pipeline_paths$statistics, "alpha_diversity")

figure_output_beta <- file.path(figure_paths$statistics, "beta_diversity")
figure_output_diff <- file.path(figure_paths$statistics, "differential_abundance")
figure_output_alpha <- file.path(figure_paths$statistics, "alpha_diversity")

input_path_taxonomic <- pipeline_paths$taxonomic
input_path_functional <- pipeline_paths$functional

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
message(glue::glue("📥 Input path (taxonomic): {input_path_taxonomic}"))
message(glue::glue("📥 Input path (functional): {input_path_functional}"))
message(glue::glue("📤 Output path (differential abundance): {output_path_diff}"))
message(glue::glue("📊 Output path (figures): {figure_output_diff}"))

# ────────────────────────────────────────────────────────────────────────────────
# Section: GitHub Functions
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
    ggplot2::ggsave(local_path, plot = object, width = 12, height = 8, dpi = 300)
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
# Section: Load Data from GitHub
# ────────────────────────────────────────────────────────────────────────────────

github_base <- paste0("https://raw.githubusercontent.com/", owner, "/", repo, "/", branch)

cat("\n═══════════════════════════════════════\n")
cat("Loading Data from GitHub\n")
cat("═══════════════════════════════════════\n")

data <- list()
files_to_load <- c(
  "02_taxonomic_profiling/species_profile_filtered.csv",
  "03_functional_profiling/functional_profile_filtered.csv"
)

for (file_name in files_to_load) {
  file_url <- file.path(github_base, "results", "Processed", file_name)
  clean_name <- tools::file_path_sans_ext(basename(file_name))
  data[[clean_name]] <- load_github_file(file_url, token)
  cat("✅", file_name, "\n")
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: Metadata Preparation
# ────────────────────────────────────────────────────────────────────────────────

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
    if (metadata$timepoint[i] %in% c("Baseline", "2_weeks")) {
      participant_id <- metadata$block_number[i]
      has_treatment <- any(
        metadata$block_number == participant_id & 
          metadata$turmeric_group == "Treatment"
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
      timepoint = factor(timepoint, levels = c("Baseline", "2_weeks", "6_months")),
      turmeric_group = factor(turmeric_group, levels = c("Control", "Treatment")),
      turmeric_status = factor(turmeric_status, levels = c("Not taken", "Taken"))
    )
  
  return(metadata)
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: Alternative Methods for Sparse Data
# ────────────────────────────────────────────────────────────────────────────────

run_wilcoxon_tests <- function(abundance_data, metadata, group_var = "turmeric_status", 
                               group1 = "Taken", group2 = "Not taken") {
  
  cat("\n═══════════════════════════════════════\n")
  cat("Running Wilcoxon Rank-Sum Tests\n")
  cat("═══════════════════════════════════════\n")
  cat(glue::glue("Comparing: {group1} vs {group2}\n\n"))
  
  # Prepare abundance matrix
  feature_col <- colnames(abundance_data)[1]
  abundance_matrix <- abundance_data %>%
    column_to_rownames(feature_col) %>%
    select(all_of(metadata$sample_id)) %>%
    as.matrix()
  
  # Get group indices
  group1_samples <- metadata %>% filter(!!sym(group_var) == group1) %>% pull(sample_id)
  group2_samples <- metadata %>% filter(!!sym(group_var) == group2) %>% pull(sample_id)
  
  cat(glue::glue("Group 1 ({group1}): {length(group1_samples)} samples\n"))
  cat(glue::glue("Group 2 ({group2}): {length(group2_samples)} samples\n\n"))
  
  # Run Wilcoxon test for each feature
  results_list <- list()
  
  for (i in 1:nrow(abundance_matrix)) {
    feature_name <- rownames(abundance_matrix)[i]
    
    group1_values <- abundance_matrix[i, group1_samples]
    group2_values <- abundance_matrix[i, group2_samples]
    
    # Skip if all zeros in both groups
    if (sum(group1_values) == 0 && sum(group2_values) == 0) {
      next
    }
    
    # Wilcoxon rank-sum test
    test_result <- wilcox.test(group1_values, group2_values, exact = FALSE)
    
    # Calculate mean abundances and fold change
    mean_group1 <- mean(group1_values)
    mean_group2 <- mean(group2_values)
    
    # Log2 fold change (add pseudocount to avoid log(0))
    log2fc <- log2((mean_group1 + 1) / (mean_group2 + 1))
    
    results_list[[feature_name]] <- data.frame(
      feature = feature_name,
      mean_group1 = mean_group1,
      mean_group2 = mean_group2,
      log2FoldChange = log2fc,
      pvalue = test_result$p.value,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results
  results_df <- bind_rows(results_list)
  
  # Adjust p-values
  results_df$padj <- p.adjust(results_df$pvalue, method = "BH")
  
  # Add significance
  results_df <- results_df %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    mutate(
      significant = ifelse(padj < 0.05, "Significant", "Not significant"),
      direction = case_when(
        padj >= 0.05 ~ "Not significant",
        log2FoldChange >= 1 ~ "Upregulated",
        log2FoldChange <= -1 ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  # Summary
  cat("\nResults Summary:\n")
  cat(glue::glue("  Total features tested: {nrow(results_df)}\n"))
  cat(glue::glue("  Significant (padj < 0.05): {sum(results_df$padj < 0.05, na.rm = TRUE)}\n"))
  cat(glue::glue("  Upregulated (LFC >= 1): {sum(results_df$direction == 'Upregulated', na.rm = TRUE)}\n"))
  cat(glue::glue("  Downregulated (LFC <= -1): {sum(results_df$direction == 'Downregulated', na.rm = TRUE)}\n\n"))
  
  return(results_df)
}

run_paired_wilcoxon_tests <- function(abundance_data, metadata, 
                                     timepoint1 = "Baseline", timepoint2 = "2_weeks") {
  
  cat("\n═══════════════════════════════════════\n")
  cat("Running Paired Wilcoxon Signed-Rank Tests\n")
  cat("═══════════════════════════════════════\n")
  cat(glue::glue("Comparing: {timepoint2} vs {timepoint1} (paired)\n\n"))
  
  # Get samples for each timepoint
  tp1_metadata <- metadata %>% filter(timepoint == timepoint1)
  tp2_metadata <- metadata %>% filter(timepoint == timepoint2)
  
  # Find paired samples (same block_number)
  paired_blocks <- intersect(tp1_metadata$block_number, tp2_metadata$block_number)
  
  cat(glue::glue("Paired individuals: {length(paired_blocks)}\n\n"))
  
  if (length(paired_blocks) < 3) {
    cat("⚠️  Too few paired samples. Cannot run paired test.\n")
    return(NULL)
  }
  
  # Get sample IDs for paired samples
  tp1_samples <- tp1_metadata %>% 
    filter(block_number %in% paired_blocks) %>%
    arrange(block_number) %>%
    pull(sample_id)
  
  tp2_samples <- tp2_metadata %>% 
    filter(block_number %in% paired_blocks) %>%
    arrange(block_number) %>%
    pull(sample_id)
  
  # Prepare abundance matrix
  feature_col <- colnames(abundance_data)[1]
  abundance_matrix <- abundance_data %>%
    column_to_rownames(feature_col) %>%
    as.matrix()
  
  # Run paired Wilcoxon test for each feature
  results_list <- list()
  
  for (i in 1:nrow(abundance_matrix)) {
    feature_name <- rownames(abundance_matrix)[i]
    
    tp1_values <- abundance_matrix[i, tp1_samples]
    tp2_values <- abundance_matrix[i, tp2_samples]
    
    # Skip if all zeros at both timepoints
    if (sum(tp1_values) == 0 && sum(tp2_values) == 0) {
      next
    }
    
    # Paired Wilcoxon signed-rank test
    test_result <- wilcox.test(tp2_values, tp1_values, paired = TRUE, exact = FALSE)
    
    # Calculate mean abundances and fold change
    mean_tp1 <- mean(tp1_values)
    mean_tp2 <- mean(tp2_values)
    
    # Log2 fold change
    log2fc <- log2((mean_tp2 + 1) / (mean_tp1 + 1))
    
    results_list[[feature_name]] <- data.frame(
      feature = feature_name,
      mean_timepoint1 = mean_tp1,
      mean_timepoint2 = mean_tp2,
      log2FoldChange = log2fc,
      pvalue = test_result$p.value,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results
  results_df <- bind_rows(results_list)
  
  # Adjust p-values
  results_df$padj <- p.adjust(results_df$pvalue, method = "BH")
  
  # Add significance
  results_df <- results_df %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    mutate(
      baseMean = (mean_timepoint1 + mean_timepoint2) / 2,
      significant = ifelse(padj < 0.05, "Significant", "Not significant"),
      direction = case_when(
        padj >= 0.05 ~ "Not significant",
        log2FoldChange >= 1 ~ "Upregulated",
        log2FoldChange <= -1 ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  # Summary
  cat("\nResults Summary:\n")
  cat(glue::glue("  Total features tested: {nrow(results_df)}\n"))
  cat(glue::glue("  Significant (padj < 0.05): {sum(results_df$padj < 0.05, na.rm = TRUE)}\n"))
  cat(glue::glue("  Upregulated: {sum(results_df$direction == 'Upregulated', na.rm = TRUE)}\n"))
  cat(glue::glue("  Downregulated: {sum(results_df$direction == 'Downregulated', na.rm = TRUE)}\n\n"))
  
  return(results_df)
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: Differential Abundance Functions - DESeq2
# ────────────────────────────────────────────────────────────────────────────────

check_design_matrix <- function(metadata, design_formula) {
  
  cat("\nChecking design matrix...\n")
  
  # Create contingency table
  cat("\nSample distribution:\n")
  table_result <- table(metadata$turmeric_status, metadata$timepoint)
  print(table_result)
  cat("\n")
  
  # Check if design matrix is full rank
  model_matrix <- model.matrix(as.formula(design_formula), data = metadata)
  rank_mm <- qr(model_matrix)$rank
  ncol_mm <- ncol(model_matrix)
  
  cat(glue::glue("Model matrix rank: {rank_mm}\n"))
  cat(glue::glue("Model matrix columns: {ncol_mm}\n"))
  
  if (rank_mm < ncol_mm) {
    cat("⚠️  Warning: Model matrix is not full rank!\n")
    cat("    Interaction term cannot be estimated.\n")
    cat("    Will use simpler model without interaction.\n\n")
    return(FALSE)
  } else {
    cat("✅ Model matrix is full rank.\n\n")
    return(TRUE)
  }
}

run_deseq2_analysis <- function(abundance_data, metadata, design_formula = "~ turmeric_status + timepoint + turmeric_status:timepoint") {
  
  cat("\n═══════════════════════════════════════\n")
  cat("Running DESeq2 Analysis\n")
  cat("═══════════════════════════════════════\n")
  cat("Design formula:", design_formula, "\n")
  
  # Check if design matrix is full rank
  is_full_rank <- check_design_matrix(metadata, design_formula)
  
  if (!is_full_rank) {
    # Use simpler model without interaction
    design_formula <- "~ turmeric_status + timepoint"
    cat("Updated design formula:", design_formula, "\n\n")
  }
  
  # Prepare count matrix
  feature_col <- colnames(abundance_data)[1]
  count_matrix <- abundance_data %>%
    column_to_rownames(feature_col) %>%
    select(all_of(metadata$sample_id)) %>%
    as.matrix()
  
  # Convert to integer counts (required for DESeq2)
  # Round and ensure non-negative
  count_matrix <- round(count_matrix)
  count_matrix[count_matrix < 0] <- 0
  
  # Remove features with zero counts across all samples
  count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]
  
  # Filter features: keep only those present in at least 10% of samples
  min_samples <- ceiling(ncol(count_matrix) * 0.1)
  count_matrix <- count_matrix[rowSums(count_matrix > 0) >= min_samples, ]
  
  cat("Features after filtering:", nrow(count_matrix), "\n")
  cat("Samples:", ncol(count_matrix), "\n\n")
  
  # Check if we have enough features
  if (nrow(count_matrix) < 10) {
    stop("❌ Too few features remaining after filtering. Cannot run DESeq2.")
  }
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = as.formula(design_formula)
  )
  
  # Run DESeq2 with appropriate size factor estimation
  cat("Running DESeq2...\n")
  cat("Estimating size factors using 'poscounts' method (for sparse data)...\n")
  
  # Use type="poscounts" for data with many zeros
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  
  cat("✅ DESeq2 analysis complete\n\n")
  
  return(dds)
}

extract_deseq2_results <- function(dds, contrast, alpha = 0.05, lfc_threshold = 1) {
  
  cat(glue::glue("Extracting results for contrast: {paste(contrast, collapse = ' vs ')}\n"))
  
  # Get results
  res <- results(dds, contrast = contrast, alpha = alpha)
  
  # Convert to data frame
  res_df <- as.data.frame(res) %>%
    rownames_to_column("feature") %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  # Add significance categories
  res_df <- res_df %>%
    mutate(
      significant = ifelse(padj < alpha & abs(log2FoldChange) >= lfc_threshold, 
                          "Significant", "Not significant"),
      direction = case_when(
        padj >= alpha ~ "Not significant",
        log2FoldChange >= lfc_threshold ~ "Upregulated",
        log2FoldChange <= -lfc_threshold ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  # Summary
  cat("\nResults Summary:\n")
  cat(glue::glue("  Total features: {nrow(res_df)}\n"))
  cat(glue::glue("  Significant (padj < {alpha}): {sum(res_df$padj < alpha, na.rm = TRUE)}\n"))
  cat(glue::glue("  Upregulated (LFC >= {lfc_threshold}): {sum(res_df$direction == 'Upregulated', na.rm = TRUE)}\n"))
  cat(glue::glue("  Downregulated (LFC <= -{lfc_threshold}): {sum(res_df$direction == 'Downregulated', na.rm = TRUE)}\n\n"))
  
  return(res_df)
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: Visualization Functions
# ────────────────────────────────────────────────────────────────────────────────

plot_volcano <- function(results_df, title = "Volcano Plot", 
                        alpha = 0.05, lfc_threshold = 1, 
                        label_top = 10) {
  
  # Prepare data
  plot_data <- results_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      neg_log10_padj = -log10(padj),
      label = ifelse(row_number() <= label_top, feature, "")
    )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_padj, 
                             color = direction, label = label)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), 
               linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(alpha), 
               linetype = "dashed", color = "gray50") +
    geom_text_repel(size = 3, max.overlaps = 20, 
                    box.padding = 0.5, point.padding = 0.3) +
    scale_color_manual(values = c(
      "Upregulated" = "#E74C3C",
      "Downregulated" = "#3498DB",
      "Not significant" = "gray70"
    )) +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Regulation"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

plot_ma <- function(results_df, title = "MA Plot", alpha = 0.05) {
  
  # Prepare data
  plot_data <- results_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange), !is.na(baseMean)) %>%
    mutate(
      log10_baseMean = log10(baseMean + 1),
      significant = ifelse(padj < alpha, "Significant", "Not significant")
    )
  
  p <- ggplot(plot_data, aes(x = log10_baseMean, y = log2FoldChange, 
                             color = significant)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = c(
      "Significant" = "#E74C3C",
      "Not significant" = "gray70"
    )) +
    labs(
      title = title,
      x = "Log10 Mean Expression",
      y = "Log2 Fold Change",
      color = "Status"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

plot_top_features_heatmap <- function(abundance_data, metadata, results_df, 
                                     top_n = 30, profile_type = "Species") {
  
  # Get top features
  top_features <- results_df %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    head(top_n) %>%
    pull(feature)
  
  # Prepare abundance matrix
  feature_col <- colnames(abundance_data)[1]
  abundance_matrix <- abundance_data %>%
    filter(!!sym(feature_col) %in% top_features) %>%
    column_to_rownames(feature_col) %>%
    select(all_of(metadata$sample_id)) %>%
    as.matrix()
  
  # Log transform and scale
  abundance_matrix_log <- log10(abundance_matrix + 1)
  abundance_matrix_scaled <- t(scale(t(abundance_matrix_log)))
  
  # Prepare annotation
  annotation_col <- metadata %>%
    column_to_rownames("sample_id") %>%
    select(turmeric_status, timepoint)
  
  # Create heatmap
  pheatmap(
    abundance_matrix_scaled,
    annotation_col = annotation_col,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 8,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = glue::glue("{profile_type} - Top {top_n} Differential Features"),
    border_color = NA
  )
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: Main Analysis Loop
# ────────────────────────────────────────────────────────────────────────────────

cat("\n═══════════════════════════════════════\n")
cat("STARTING DIFFERENTIAL ABUNDANCE ANALYSIS\n")
cat("═══════════════════════════════════════\n")

all_diff_results <- list()

# Define profile mappings
profile_mappings <- list(
  Species = "species_profile_filtered",
  Functional = "functional_profile_filtered"
)

for (profile_type in names(profile_mappings)) {
  
  cat("\n", strrep("=", 80), "\n")
  cat(glue::glue("ANALYZING: {profile_type} Profile\n"))
  cat(strrep("=", 80), "\n")
  
  # Get profile data
  profile_name <- profile_mappings[[profile_type]]
  abundance_data <- data[[profile_name]]
  
  # Prepare metadata
  sample_ids <- colnames(abundance_data)[-1]
  metadata <- prepare_metadata(sample_ids)
  
  # Keep only samples in metadata
  abundance_data <- abundance_data %>%
    select(1, all_of(metadata$sample_id))
  
  # Try DESeq2 first
  deseq2_success <- TRUE
  dds <- NULL
  
  tryCatch({
    # For repeated measures, include individual as blocking factor
    dds <- run_deseq2_analysis(abundance_data, metadata, 
                               design_formula = "~ block_number + timepoint")
  }, error = function(e) {
    cat("\n⚠️  DESeq2 with blocking failed. Trying simpler model...\n")
    tryCatch({
      dds <<- run_deseq2_analysis(abundance_data, metadata, 
                                 design_formula = "~ timepoint")
    }, error = function(e2) {
      cat("\n⚠️  DESeq2 failed with error:\n")
      cat(as.character(e2), "\n")
      cat("\nFalling back to Wilcoxon rank-sum tests...\n")
      deseq2_success <<- FALSE
    })
  })
  
  results_list <- list()
  
  if (deseq2_success && !is.null(dds)) {
    # ─────────────────────────────────────────────────────────────────────────────
    # DESeq2 Analysis - Timepoint Effects (Repeated Measures)
    # ─────────────────────────────────────────────────────────────────────────────
    
    cat("\n", strrep("=", 60), "\n")
    cat("TIMEPOINT EFFECTS (Repeated Measures Design)\n")
    cat(strrep("=", 60), "\n\n")
    
    # Extract results for timepoint contrasts
    timepoint_contrasts <- list(
      timepoint_2weeks_vs_baseline = c("timepoint", "2_weeks", "Baseline"),
      timepoint_6months_vs_baseline = c("timepoint", "6_months", "Baseline"),
      timepoint_6months_vs_2weeks = c("timepoint", "6_months", "2_weeks")
    )
    
    for (contrast_name in names(timepoint_contrasts)) {
      cat("\n", strrep("-", 60), "\n")
      cat(glue::glue("Contrast: {contrast_name}\n"))
      cat(strrep("-", 60), "\n")
      
      results_df <- extract_deseq2_results(dds, timepoint_contrasts[[contrast_name]])
      results_list[[contrast_name]] <- results_df
      
      # Save results
      save_and_upload_to_github(
        object = results_df,
        object_type = "data",
        github_path = file.path(output_path_diff,
                                tolower(glue::glue("{profile_type}_{contrast_name}_deseq2_results.csv"))),
        commit_message = glue::glue("Add {profile_type} {contrast_name} DESeq2 results"),
        owner = owner, repo = repo, branch = branch
      )
      
      # Create volcano plot
      volcano_plot <- plot_volcano(results_df, 
                                   title = glue::glue("{profile_type} - {contrast_name}"))
      
      save_and_upload_to_github(
        object = volcano_plot,
        object_type = "plot",
        github_path = file.path(figure_output_diff,
                                tolower(glue::glue("{profile_type}_{contrast_name}_volcano.png"))),
        commit_message = glue::glue("Add {profile_type} {contrast_name} volcano plot"),
        owner = owner, repo = repo, branch = branch
      )
      
      # Create MA plot
      ma_plot <- plot_ma(results_df, 
                        title = glue::glue("{profile_type} - {contrast_name}"))
      
      save_and_upload_to_github(
        object = ma_plot,
        object_type = "plot",
        github_path = file.path(figure_output_diff,
                                tolower(glue::glue("{profile_type}_{contrast_name}_ma.png"))),
        commit_message = glue::glue("Add {profile_type} {contrast_name} MA plot"),
        owner = owner, repo = repo, branch = branch
      )
    }
  } else {
    # ─────────────────────────────────────────────────────────────────────────────
    # Paired Wilcoxon Tests (fallback for sparse data)
    # ─────────────────────────────────────────────────────────────────────────────
    
    cat("\n═══════════════════════════════════════\n")
    cat("Using Paired Wilcoxon Tests for Timepoint Effects\n")
    cat("═══════════════════════════════════════\n\n")
    
    # Run paired tests for timepoint comparisons
    timepoint_pairs <- list(
      "2_weeks_vs_Baseline" = c("Baseline", "2_weeks"),
      "6_months_vs_Baseline" = c("Baseline", "6_months"),
      "6_months_vs_2_weeks" = c("2_weeks", "6_months")
    )
    
    for (comparison_name in names(timepoint_pairs)) {
      tp_pair <- timepoint_pairs[[comparison_name]]
      
      # Run paired Wilcoxon test
      paired_results <- run_paired_wilcoxon_tests(
        abundance_data, metadata, 
        timepoint1 = tp_pair[1], 
        timepoint2 = tp_pair[2]
      )
      
      results_list[[comparison_name]] <- paired_results
      
      # Save results
      save_and_upload_to_github(
        object = paired_results,
        object_type = "data",
        github_path = file.path(output_path_diff,
                                tolower(glue::glue("{profile_type}_{comparison_name}_paired_wilcoxon_results.csv"))),
        commit_message = glue::glue("Add {profile_type} {comparison_name} paired Wilcoxon results"),
        owner = owner, repo = repo, branch = branch
      )
      
      # Create volcano plot
      volcano_plot <- plot_volcano(paired_results, 
                                   title = glue::glue("{profile_type} - {comparison_name} (Paired)"))
      
      save_and_upload_to_github(
        object = volcano_plot,
        object_type = "plot",
        github_path = file.path(figure_output_diff,
                                tolower(glue::glue("{profile_type}_{comparison_name}_paired_volcano.png"))),
        commit_message = glue::glue("Add {profile_type} {comparison_name} paired volcano plot"),
        owner = owner, repo = repo, branch = branch
      )
    }
  }
  
  # ─────────────────────────────────────────────────────────────────────────────
  # Turmeric Effect at 6 Months (Taken vs Not Taken)
  # ─────────────────────────────────────────────────────────────────────────────
  
  cat("\n", strrep("=", 60), "\n")
  cat("TURMERIC EFFECT AT 6 MONTHS (Taken vs Not Taken)\n")
  cat(strrep("=", 60), "\n\n")
  
  # Subset to 6 months only
  metadata_6m <- metadata %>% filter(timepoint == "6_months")
  
  if (nrow(metadata_6m) >= 4 && length(unique(metadata_6m$turmeric_status)) == 2) {
    
    # Subset abundance data
    abundance_6m <- abundance_data %>%
      select(1, all_of(metadata_6m$sample_id))
    
    # Try DESeq2 first
    turmeric_success <- TRUE
    
    tryCatch({
      dds_6m <- run_deseq2_analysis(abundance_6m, metadata_6m, 
                                    design_formula = "~ turmeric_status")
      
      turmeric_results <- extract_deseq2_results(dds_6m, 
                                                 c("turmeric_status", "Taken", "Not taken"))
      results_list$turmeric_effect_6months <- turmeric_results
      
      # Save results
      save_and_upload_to_github(
        object = turmeric_results,
        object_type = "data",
        github_path = file.path(output_path_diff,
                                tolower(glue::glue("{profile_type}_turmeric_effect_6months_deseq2_results.csv"))),
        commit_message = glue::glue("Add {profile_type} turmeric effect at 6 months DESeq2 results"),
        owner = owner, repo = repo, branch = branch
      )
      
      # Create volcano plot
      volcano_plot <- plot_volcano(turmeric_results, 
                                   title = glue::glue("{profile_type} - Turmeric Effect at 6 Months"))
      
      save_and_upload_to_github(
        object = volcano_plot,
        object_type = "plot",
        github_path = file.path(figure_output_diff,
                                tolower(glue::glue("{profile_type}_turmeric_effect_6months_volcano.png"))),
        commit_message = glue::glue("Add {profile_type} turmeric effect at 6 months volcano plot"),
        owner = owner, repo = repo, branch = branch
      )
      
    }, error = function(e) {
      cat("\n⚠️  DESeq2 failed for turmeric effect. Using Wilcoxon test...\n")
      
      turmeric_results <- run_wilcoxon_tests(abundance_6m, metadata_6m, 
                                            "turmeric_status", "Taken", "Not taken")
      results_list$turmeric_effect_6months <<- turmeric_results
      
      # Save results
      save_and_upload_to_github(
        object = turmeric_results,
        object_type = "data",
        github_path = file.path(output_path_diff,
                                tolower(glue::glue("{profile_type}_turmeric_effect_6months_wilcoxon_results.csv"))),
        commit_message = glue::glue("Add {profile_type} turmeric effect at 6 months Wilcoxon results"),
        owner = owner, repo = repo, branch = branch
      )
      
      # Create volcano plot
      volcano_plot <- plot_volcano(turmeric_results, 
                                   title = glue::glue("{profile_type} - Turmeric Effect at 6 Months (Wilcoxon)"))
      
      save_and_upload_to_github(
        object = volcano_plot,
        object_type = "plot",
        github_path = file.path(figure_output_diff,
                                tolower(glue::glue("{profile_type}_turmeric_effect_6months_wilcoxon_volcano.png"))),
        commit_message = glue::glue("Add {profile_type} turmeric effect at 6 months Wilcoxon volcano plot"),
        owner = owner, repo = repo, branch = branch
      )
    })
    
  } else {
    cat("⚠️  Insufficient samples for turmeric effect analysis at 6 months.\n")
  }
  
  # ─────────────────────────────────────────────────────────────────────────────
  # Create heatmap for top features
  # ─────────────────────────────────────────────────────────────────────────────
  heatmap_file <- tempfile(fileext = ".png")
  png(heatmap_file, width = 12, height = 10, units = "in", res = 300)
  plot_top_features_heatmap(abundance_data, metadata, 
                           results_list$turmeric_effect, 
                           top_n = 30, profile_type = profile_type)
  dev.off()
  
  save_and_upload_to_github(
    object = NULL,
    object_type = "file",
    local_path = heatmap_file,
    github_path = file.path(figure_output_diff,
                            tolower(glue::glue("{profile_type}_top_features_heatmap.png"))),
    commit_message = glue::glue("Add {profile_type} top features heatmap"),
    owner = owner, repo = repo, branch = branch
  )
  
  unlink(heatmap_file)
  
  # Store results
  all_diff_results[[profile_type]] <- list(
    dds = dds,
    results = results_list
  )
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: Summary
# ────────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 80), "\n")
cat("DIFFERENTIAL ABUNDANCE ANALYSIS COMPLETE\n")
cat(strrep("=", 80), "\n")
cat("\nResults saved to:", output_path_diff, "\n")
cat("Figures saved to:", figure_output_diff, "\n")
cat("\nAnalyses completed:\n")
cat("  ✅ DESeq2 differential abundance\n")
cat("  ✅ Volcano plots\n")
cat("  ✅ MA plots\n")
cat("  ✅ Heatmaps\n")
cat("\nProfiles analyzed:\n")
for (profile in names(all_diff_results)) {
  cat(glue::glue("  ✅ {profile}\n"))
}
cat("\n")