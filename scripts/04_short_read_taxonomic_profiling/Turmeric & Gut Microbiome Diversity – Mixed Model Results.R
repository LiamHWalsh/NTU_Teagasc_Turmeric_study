
# ────────────────────────────────────────────────────────────────────────────────
# Section: Reproducibility
# ────────────────────────────────────────────────────────────────────────────────

set.seed(123)

# ────────────────────────────────────────────────────────────────────────────────
# Section: Load Metadata Files
# ────────────────────────────────────────────────────────────────────────────────

metadata_dir <- "Q:/CoDiet/liam walsh/NottinghamJuly2025/officer_output/Reports"

if (!dir.exists(metadata_dir)) stop("❌ Metadata directory not found.")

file_list <- list.files(metadata_dir, full.names = TRUE)
valid_files <- file_list[!grepl("\\.docx$", file_list, ignore.case = TRUE)]

data <- list()

for (file_path in valid_files) {
  tryCatch({
    if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
      data[[file_path]] <- read_csv(file_path, show_col_types = FALSE)
    } else if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
      data[[file_path]] <- read_delim(file_path, delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE, show_col_types = FALSE)
    }
  }, error = function(e) {
    warning(glue::glue("❌ Error reading file {file_path}: {e$message}"))
  })
}

names(data) <- gsub(paste0(metadata_dir, "/|\\.csv$|\\.tsv$"), "", names(data))
names(data)[grep("metaphlan", names(data))] <- "metaphlan"

# ────────────────────────────────────────────────────────────────────────────────
# Section: Create Turmeric Status Lookup
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

# ────────────────────────────────────────────────────────────────────────────────
# Section: Extract Sample IDs and Create Metadata
# ────────────────────────────────────────────────────────────────────────────────

ids <- unique(data$alpha_diversity$sample_id)
split_ids <- str_split_fixed(ids, "_", 3)

sample_ids <- data.frame(
  full_id = ids,
  block_code = split_ids[,2],
  block_number = as.integer(gsub("[^0-9]", "", split_ids[,2])),
  block_letter = gsub("[0-9]", "", split_ids[,2]),
  stringsAsFactors = FALSE
)

# Merge with turmeric status
sample_ids <- merge(sample_ids, turmeric_status,
                    by.x = "block_number", by.y = "ID", all.x = TRUE)

# Assign timepoints
sample_ids$timepoint <- "Baseline"
sample_ids$timepoint[sample_ids$block_letter %in% c("a","A")] <- "Baseline"
sample_ids$timepoint[sample_ids$block_letter %in% c("b","B")] <- "2_weeks"
sample_ids$timepoint[sample_ids$block_letter %in% c("c","C")] <- "6_months"

# Create turmeric group indicator
sample_ids$turmeric_group <- "Control"

# Participants who took turmeric at 6 months visit
sample_ids$turmeric_group[
  grepl("Turmeric 6 months", sample_ids$status) & 
    sample_ids$block_letter %in% c("c", "C")
] <- "Treatment"

# Mark stopped group (3-month turmeric at 6-month visit)
sample_ids$turmeric_group[
  sample_ids$status == "Turmeric 3 months" & 
    sample_ids$block_letter == "C"
] <- "Stopped"

# For baseline and 2-week visits, assign based on eventual status
# This ensures we track the same individuals over time
for (i in 1:nrow(sample_ids)) {
  if (sample_ids$timepoint[i] %in% c("Baseline", "2_weeks")) {
    participant_id <- sample_ids$block_number[i]
    
    # Check if this participant has a treatment sample at 6 months
    has_treatment <- any(
      sample_ids$block_number == participant_id & 
        sample_ids$turmeric_group == "Treatment"
    )
    
    if (has_treatment) {
      sample_ids$turmeric_group[i] <- "Treatment"
    }
  }
}

# Standardize letter case
sample_ids$block_letter <- str_to_upper(sample_ids$block_letter)

# Merge with alpha diversity
data$alpha_diversity <- merge(sample_ids, data$alpha_diversity,
                              by.x = "full_id", by.y = "sample_id", all.x = TRUE)
data$alpha_diversity$block_letter <- str_to_upper(data$alpha_diversity$block_letter)
# ────────────────────────────────────────────────────────────────────────────────
# Section: Alpha Diversity Analysis with Mixed Models - FIXED VERSION
# ────────────────────────────────────────────────────────────────────────────────

outdir <- "Q:/CoDiet/liam walsh/NottinghamJuly2025/plots"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Filter data (exclude stopped participants)
filtered_data <- data$alpha_diversity %>%
  filter(turmeric_group != "Stopped") %>%
  mutate(
    timepoint = factor(timepoint, levels = c("Baseline", "2_weeks", "6_months")),
    turmeric_group = factor(turmeric_group, levels = c("Control", "Treatment"))
  )

# Create turmeric_status variable
filtered_data$turmeric_status <- "Not taken"
filtered_data$turmeric_status[
  filtered_data$block_letter == "C" & 
    filtered_data$status == "Turmeric 6 months"
] <- "Taken"

# Convert to factor
filtered_data$turmeric_status <- factor(filtered_data$turmeric_status, 
                                        levels = c("Not taken", "Taken"))

# Summary of data structure
cat("\n═══════════════════════════════════════\n")
cat("Data Structure Summary\n")
cat("═══════════════════════════════════════\n")
cat("Total samples:", nrow(filtered_data), "\n")
cat("Unique participants:", length(unique(filtered_data$block_number)), "\n")

# Check repeated measures structure
participant_counts <- filtered_data %>%
  group_by(block_number) %>%
  summarise(n_samples = n(), .groups = "drop")

cat("\nSamples per participant:\n")
print(table(participant_counts$n_samples))
cat("\nParticipants with only 1 sample:", sum(participant_counts$n_samples == 1), "\n")
cat("Participants with 2+ samples:", sum(participant_counts$n_samples > 1), "\n")

cat("\nSamples per status and timepoint:\n")
table(filtered_data$timepoint, filtered_data$turmeric_status) %>% print()

# Check for design issues
design_check <- filtered_data %>%
  group_by(timepoint, turmeric_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = turmeric_status, values_from = n, values_fill = 0)

cat("\nDesign matrix:\n")
print(design_check)

# Identify the issue
if (any(design_check[,-1] == 0)) {
  cat("\n⚠️  WARNING: Some group×timepoint combinations have ZERO samples!\n")
  cat("This causes rank deficiency in the interaction model.\n")
  cat("Solution: Using appropriate model based on data structure.\n\n")
}


library(brms)
library(dplyr)

# ──────────────────────────────────────────────────────
# Function: Fit Bayesian LMM for Shannon Diversity
# ──────────────────────────────────────────────────────
fit_shannon_bayesian <- function(data) {
  
  # Ensure factors are correctly specified
  data <- data %>%
    mutate(
      turmeric_status = factor(turmeric_status, levels = c("Not taken", "Taken")),
      timepoint = factor(timepoint, levels = c("Baseline", "2_weeks", "6_months"))
    )
  
  # Priors
  priors <- c(
    prior(normal(0, 5), class = "b"),          # Fixed effects
    prior(student_t(3, 0, 5), class = "sd"),   # Random intercept SD
    prior(student_t(3, 0, 10), class = "sigma") # Residual SD
  )
  
  # Fit model
  b_model <- brm(
    formula = Shannon ~ turmeric_status + timepoint + (1 | block_number),
    data = data,
    family = gaussian(),
    prior = priors,
    chains = 4,
    iter = 4000,
    warmup = 2000,
    cores = 4,
    seed = 123,
    control = list(adapt_delta = 0.95, max_treedepth = 15)
  )
  
  # Posterior check: probability that turmeric decreases Shannon diversity
  turmeric_post <- hypothesis(b_model, "turmeric_statusTaken < 0")
  
  # Summarize posterior distributions
  summary_model <- summary(b_model)
  
  list(
    model = b_model,
    summary = summary_model,
    turmeric_hypothesis = turmeric_post
  )
}

# ──────────────────────────────────────────────────────
# Run analysis
# ──────────────────────────────────────────────────────
results <- fit_shannon_bayesian(filtered_data)

# ──────────────────────────────────────────────────────
# View results
# ──────────────────────────────────────────────────────
results$summary
results$turmeric_hypothesis







# ────────────────────────────────────────────────────────────────────────────────
# Function for robust linear mixed model analysis
# ────────────────────────────────────────────────────────────────────────────────

run_lmm_robust <- function(data, response_var, group_var = "turmeric_status") {
  # Check data structure
  design_matrix <- table(data$timepoint, data[[group_var]])
  n_timepoints <- length(unique(data$timepoint))
  n_groups <- length(unique(data[[group_var]]))
  has_full_design <- all(design_matrix > 0)
  
  # Check repeated measures structure
  participant_counts <- data %>%
    group_by(block_number) %>%
    summarise(n_obs = n(), .groups = "drop")
  
  n_participants <- nrow(participant_counts)
  n_observations <- nrow(data)
  participants_with_repeats <- sum(participant_counts$n_obs > 1)
  
  cat("\n═══════════════════════════════════════\n")
  cat("Linear Mixed Model for:", response_var, "\n")
  cat("═══════════════════════════════════════\n")
  cat("\nData Structure:\n")
  cat("  Total observations:", n_observations, "\n")
  cat("  Unique participants:", n_participants, "\n")
  cat("  Participants with 2+ samples:", participants_with_repeats, "\n")
  cat("  Ratio (obs/participants):", round(n_observations / n_participants, 2), "\n")
  
  cat("\nDesign Matrix (samples per cell):\n")
  print(design_matrix)
  cat("\n")
  
  # Decision tree for model selection
  can_use_random_effects <- (n_observations > n_participants + 5) && (participants_with_repeats >= 3)
  
  if (!can_use_random_effects) {
    cat("⚠️  INSUFFICIENT DATA for random effects (need obs >> participants)\n")
    cat("Solution: Using LINEAR MODEL (no random effects)\n\n")
    
    # Use regular linear model
    if (!has_full_design) {
      formula_str <- paste0(response_var, " ~ ", group_var, " + timepoint")
      cat("Model formula:", formula_str, "\n\n")
    } else if (n_timepoints > 1) {
      formula_str <- paste0(response_var, " ~ ", group_var, " * timepoint")
      cat("Model formula:", formula_str, "\n\n")
    } else {
      formula_str <- paste0(response_var, " ~ ", group_var)
      cat("Model formula:", formula_str, "\n\n")
    }
    
    model <- lm(as.formula(formula_str), data = data)
    print(summary(model))
    
    # ANOVA
    cat("\n--- ANOVA ---\n")
    print(anova(model))
    
    # Emmeans
    cat("\n--- Estimated Marginal Means ---\n")
    
    if (n_timepoints > 1) {
      cat("\n1. Main Effect of Treatment:\n")
      emm_group <- emmeans(model, specs = group_var)
      print(emm_group)
      cat("\nPairwise comparisons:\n")
      print(pairs(emm_group, adjust = "fdr"))
      
      cat("\n2. Main Effect of Timepoint:\n")
      emm_time <- emmeans(model, specs = "timepoint")
      print(emm_time)
      cat("\nPairwise comparisons:\n")
      print(pairs(emm_time, adjust = "fdr"))
      
      if (has_full_design) {
        cat("\n3. Treatment Effect at Each Timepoint:\n")
        emm_simple <- emmeans(model, specs = as.formula(paste0("~ ", group_var, " | timepoint")))
        print(pairs(emm_simple, adjust = "fdr"))
        
        return(list(
          model = model,
          model_type = "lm",
          emm_group = emm_group,
          emm_time = emm_time,
          emm_simple = emm_simple,
          balanced = has_full_design
        ))
      } else {
        return(list(
          model = model,
          model_type = "lm",
          emm_group = emm_group,
          emm_time = emm_time,
          balanced = has_full_design
        ))
      }
    } else {
      emm <- emmeans(model, specs = group_var)
      print(emm)
      cat("\n")
      print(pairs(emm, adjust = "fdr"))
      
      return(list(
        model = model,
        model_type = "lm",
        emmeans = emm,
        balanced = has_full_design
      ))
    }
    
  } else {
    # Can use mixed model
    cat("✓ SUFFICIENT DATA for mixed effects model\n")
    
    if (!has_full_design) {
      cat("⚠️  UNBALANCED DESIGN: Some cells have zero samples\n")
      cat("Strategy: Using MAIN EFFECTS with random intercept\n\n")
      
      empty_cells <- which(design_matrix == 0, arr.ind = TRUE)
      if (nrow(empty_cells) > 0) {
        cat("Empty cells:\n")
        for (i in 1:nrow(empty_cells)) {
          cat(sprintf("  - %s × %s\n", 
                      rownames(design_matrix)[empty_cells[i, "row"]],
                      colnames(design_matrix)[empty_cells[i, "col"]]))
        }
        cat("\n")
      }
      
      if (n_timepoints > 1) {
        formula_str <- paste0(response_var, " ~ ", group_var, " + timepoint + (1|block_number)")
        cat("Model formula:", formula_str, "\n\n")
        
        model <- lmer(as.formula(formula_str), data = data)
        print(summary(model))
        
        cat("\n--- Type III ANOVA ---\n")
        print(anova(model, type = 3))
        
        cat("\n--- Estimated Marginal Means ---\n")
        
        cat("\n1. Main Effect of Treatment:\n")
        emm_group <- emmeans(model, specs = group_var)
        print(emm_group)
        cat("\nPairwise comparisons:\n")
        print(pairs(emm_group, adjust = "fdr"))
        
        cat("\n2. Main Effect of Timepoint:\n")
        emm_time <- emmeans(model, specs = "timepoint")
        print(emm_time)
        cat("\nPairwise comparisons:\n")
        print(pairs(emm_time, adjust = "fdr"))
        
        return(list(
          model = model,
          model_type = "lmer",
          emm_group = emm_group,
          emm_time = emm_time,
          balanced = FALSE
        ))
      } else {
        formula_str <- paste0(response_var, " ~ ", group_var, " + (1|block_number)")
        cat("Model formula:", formula_str, "\n\n")
        
        model <- lmer(as.formula(formula_str), data = data)
        print(summary(model))
        
        emm <- emmeans(model, specs = group_var)
        print(emm)
        cat("\n")
        print(pairs(emm, adjust = "fdr"))
        
        return(list(
          model = model,
          model_type = "lmer",
          emmeans = emm,
          balanced = FALSE
        ))
      }
      
    } else {
      # Balanced design with sufficient data
      cat("✓ BALANCED DESIGN with sufficient repeated measures\n")
      cat("Strategy: Using FULL INTERACTION with random intercept\n\n")
      
      if (n_timepoints > 1) {
        formula_str <- paste0(response_var, " ~ ", group_var, " * timepoint + (1|block_number)")
        cat("Model formula:", formula_str, "\n\n")
        
        model <- lmer(as.formula(formula_str), data = data)
        print(summary(model))
        
        cat("\n--- Type III ANOVA ---\n")
        print(anova(model, type = 3))
        
        cat("\n--- Estimated Marginal Means ---\n")
        
        cat("\n1. Main Effect of Treatment:\n")
        emm_group <- emmeans(model, specs = group_var)
        print(emm_group)
        cat("\n")
        print(pairs(emm_group, adjust = "fdr"))
        
        cat("\n2. Main Effect of Timepoint:\n")
        emm_time <- emmeans(model, specs = "timepoint")
        print(emm_time)
        cat("\n")
        print(pairs(emm_time, adjust = "fdr"))
        
        cat("\n3. Treatment Effect at Each Timepoint:\n")
        emm_simple <- emmeans(model, specs = as.formula(paste0("~ ", group_var, " | timepoint")))
        print(pairs(emm_simple, adjust = "fdr"))
        
        return(list(
          model = model,
          model_type = "lmer",
          emm_group = emm_group,
          emm_time = emm_time,
          emm_simple = emm_simple,
          balanced = TRUE
        ))
      } else {
        formula_str <- paste0(response_var, " ~ ", group_var, " + (1|block_number)")
        cat("Model formula:", formula_str, "\n\n")
        
        model <- lmer(as.formula(formula_str), data = data)
        print(summary(model))
        
        emm <- emmeans(model, specs = group_var)
        print(emm)
        cat("\n")
        print(pairs(emm, adjust = "fdr"))
        
        return(list(
          model = model,
          model_type = "lmer",
          emmeans = emm,
          balanced = TRUE
        ))
      }
    }
  }
}

# ────────────────────────────────────────────────────────────────────────────────
# ALTERNATIVE: Analysis by Timepoint (always works)
# ────────────────────────────────────────────────────────────────────────────────

analyze_by_timepoint <- function(data, response_var, group_var = "turmeric_status") {
  cat("\n\n═══════════════════════════════════════\n")
  cat("TIMEPOINT-SPECIFIC ANALYSIS:", response_var, "\n")
  cat("═══════════════════════════════════════\n")
  cat("Analyzing each timepoint separately\n")
  cat("(This approach doesn't require repeated measures)\n\n")
  
  results_list <- list()
  
  for (tp in unique(data$timepoint)) {
    tp_data <- data %>% filter(timepoint == tp)
    
    # Get group counts dynamically
    group_counts <- table(tp_data[[group_var]])
    
    cat("\n──────────────────────────────────────\n")
    cat("Timepoint:", tp, "\n")
    for (grp in names(group_counts)) {
      cat(sprintf("  %s: n = %d\n", grp, group_counts[grp]))
    }
    
    # Check if we have at least 2 groups with sufficient data
    if (length(group_counts) >= 2 && all(group_counts >= 3)) {
      # Summary statistics
      stats <- tp_data %>%
        group_by(.data[[group_var]]) %>%
        summarise(
          n = n(),
          mean = mean(.data[[response_var]], na.rm = TRUE),
          sd = sd(.data[[response_var]], na.rm = TRUE),
          se = sd / sqrt(n),
          .groups = "drop"
        )
      print(stats)
      
      # Test normality
      shapiro_p <- shapiro.test(tp_data[[response_var]])$p.value
      cat(sprintf("\nNormality test (Shapiro-Wilk): p = %.4f %s\n", 
                  shapiro_p,
                  ifelse(shapiro_p > 0.05, "(normal)", "(non-normal)")))
      
      # Choose appropriate test
      if (shapiro_p > 0.05 && nrow(tp_data) >= 6) {
        # Parametric test
        test <- t.test(
          as.formula(paste(response_var, "~", group_var)),
          data = tp_data,
          var.equal = FALSE
        )
        cat(sprintf("Independent t-test: t(%.1f) = %.2f, p = %.4f %s\n",
                    test$parameter, test$statistic, test$p.value,
                    ifelse(test$p.value < 0.05, "***", "")))
        
        results_list[[tp]] <- list(
          timepoint = tp,
          test_type = "t-test",
          statistic = test$statistic,
          p_value = test$p.value,
          group_counts = as.list(group_counts)
        )
      } else {
        # Non-parametric test
        test <- wilcox.test(
          as.formula(paste(response_var, "~", group_var)),
          data = tp_data,
          exact = FALSE
        )
        cat(sprintf("Wilcoxon test: W = %.1f, p = %.4f %s\n",
                    test$statistic, test$p.value,
                    ifelse(test$p.value < 0.05, "***", "")))
        
        results_list[[tp]] <- list(
          timepoint = tp,
          test_type = "Wilcoxon",
          statistic = test$statistic,
          p_value = test$p.value,
          group_counts = as.list(group_counts)
        )
      }
      
      # Effect size
      if (!requireNamespace("effsize", quietly = TRUE)) {
        install.packages("effsize")
      }
      library(effsize)
      
      es <- cohen.d(as.formula(paste(response_var, "~", group_var)), 
                    data = tp_data, na.rm = TRUE)
      cat(sprintf("Effect size (Cohen's d): %.2f (%s)\n", 
                  es$estimate, es$magnitude))
      
      results_list[[tp]]$cohens_d <- es$estimate
      results_list[[tp]]$effect_magnitude <- es$magnitude
      
    } else {
      cat("⚠️  Insufficient data for comparison (need n≥3 per group)\n")
      results_list[[tp]] <- list(
        timepoint = tp,
        test_type = "Insufficient data",
        group_counts = as.list(group_counts)
      )
    }
  }
  
  # Multiple testing correction
  p_values <- sapply(results_list, function(x) x$p_value)
  p_values <- p_values[!is.na(p_values)]
  
  if (length(p_values) > 1) {
    cat("\n──────────────────────────────────────\n")
    cat("Multiple Testing Correction (FDR)\n")
    cat("──────────────────────────────────────\n")
    
    adjusted_p <- p.adjust(p_values, method = "fdr")
    
    correction_df <- data.frame(
      Timepoint = names(p_values),
      p_value = p_values,
      q_value = adjusted_p,
      Significant = ifelse(adjusted_p < 0.05, "Yes", "No")
    )
    print(correction_df)
    
    # Add adjusted p-values to results
    for (tp in names(adjusted_p)) {
      results_list[[tp]]$q_value <- adjusted_p[tp]
    }
  }
  
  return(results_list)
}

# ────────────────────────────────────────────────────────────────────────────────
# Run all analyses
# ────────────────────────────────────────────────────────────────────────────────

# Try mixed model approach first (will fall back to lm if needed)
cat("\n")
cat("════════════════════════════════════════════════════════\n")
cat("APPROACH 1: MIXED MODEL (or LM if insufficient data)\n")
cat("════════════════════════════════════════════════════════\n")

lmm_shannon <- run_lmm_robust(filtered_data, "Shannon", group_var = "turmeric_status")
lmm_richness <- run_lmm_robust(filtered_data, "Richness", group_var = "turmeric_status")
lmm_simpson <- run_lmm_robust(filtered_data, "Simpson", group_var = "turmeric_status")

# Always run timepoint-specific analysis as supplementary
cat("\n\n════════════════════════════════════════════════════════\n")
cat("APPROACH 2: TIMEPOINT-SPECIFIC ANALYSES\n")
cat("════════════════════════════════════════════════════════\n")
cat("(Independent tests at each timepoint - no repeated measures needed)\n")

shannon_by_tp <- analyze_by_timepoint(filtered_data, "Shannon", group_var = "turmeric_status")
richness_by_tp <- analyze_by_timepoint(filtered_data, "Richness", group_var = "turmeric_status")
simpson_by_tp <- analyze_by_timepoint(filtered_data, "Simpson", group_var = "turmeric_status")

# ────────────────────────────────────────────────────────────────────────────────
# Visualization
# ────────────────────────────────────────────────────────────────────────────────

make_alpha_plot <- function(data, y_var, filename, y_label) {
  p <- ggplot(data, aes(x = timepoint, y = .data[[y_var]], 
                     fill = turmeric_status, color = turmeric_status)) +
    geom_line(aes(group = block_number), alpha = 0.3, linewidth = 0.5) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(0.8)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                size = 3, alpha = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    scale_colour_brewer(palette = "Set2") +
    labs(
      title = paste(y_label, "across timepoints"),
      x = "Timepoint",
      y = y_label,
      fill = "Turmeric Status",
      color = "Turmeric Status"
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top"
    )
  
  ggsave(file.path(outdir, filename), plot = p, width = 10, height = 6, dpi = 300)
  return(p)
}

# Generate plots
make_alpha_plot(filtered_data, "Shannon", "alpha_Shannon_lmm.png", "Shannon Diversity")
make_alpha_plot(filtered_data, "Richness", "alpha_Richness_lmm.png", "Species Richness")
make_alpha_plot(filtered_data, "Simpson", "alpha_Simpson_lmm.png", "Simpson Diversity")

cat("\n✅ Robust analysis complete!\n")
cat("\nModel types used:\n")
cat("  Shannon:", lmm_shannon$model_type, "\n")
cat("  Richness:", lmm_richness$model_type, "\n")
cat("  Simpson:", lmm_simpson$model_type, "\n")