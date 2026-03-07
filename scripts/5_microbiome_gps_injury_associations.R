suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(readxl)
  library(purrr)
})

# -----------------------------------------------------------------------------
# Exploratory association analysis:
# microbiome diversity/signals vs GPS + injury/illness metadata
# -----------------------------------------------------------------------------

GPS_FILE <- Sys.getenv("GPS_INJURY_XLSX", unset = "")
GPS_PREPARED_FILE <- "results/Processed/05_statistical_analysis/gps_injury_associations/gps_injury_clean_all_participants.csv"
ALPHA_FILE <- "results/Reports/alpha_diversity.csv"
SPECIES_FILE <- "results/Reports/species_profile_filtered.csv"
PATHWAY_FILE <- "results/Reports/HUMAnN_merged_pathabundance_cpm.tsv"

STRAIN_EFF_FILE <- "results/Processed/05_statistical_analysis/honest_small_n_strain/strain_effect_sizes.csv"
GENUS_EFF_FILE <- "results/Processed/05_statistical_analysis/honest_small_n_strain/genus_effect_sizes.csv"
PATHWAY_EFF_FILE <- "results/Processed/05_statistical_analysis/honest_small_n_strain/pathway_effect_sizes.csv"

OUT_DIR <- "results/Processed/05_statistical_analysis/gps_injury_associations"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_N_ASSOC <- 8

# Cohort coding used across existing scripts/manuscript
turmeric_status <- tibble::tibble(
  participant_id = c(3, 4, 5, 6, 7, 9, 11, 14, 15, 16, 18, 19, 20, 21, 22, 23),
  status = c(
    "Didn't take turmeric", "Turmeric 6 months", "Turmeric 6 months",
    "Didn't take turmeric", "Turmeric 6 months", "Turmeric 6 months",
    "Turmeric 3 months", "Baseline sample only", "Turmeric 6 months",
    "Turmeric 6 months", "Didn't take turmeric", "Baseline sample only",
    "Turmeric 3 months", "Turmeric 6 months", "Turmeric 6 months",
    "Baseline sample only"
  )
) %>%
  mutate(
    turmeric_group = case_when(
      status == "Turmeric 6 months" ~ "Treatment",
      status == "Didn't take turmeric" ~ "Control",
      TRUE ~ "Excluded"
    )
  )

analysis_ids <- turmeric_status %>%
  filter(turmeric_group %in% c("Treatment", "Control")) %>%
  pull(participant_id)

parse_num <- function(x) {
  suppressWarnings(
    readr::parse_number(as.character(x), locale = readr::locale(decimal_mark = "."))
  )
}

validate_gps_layout <- function(gps_raw) {
  normalize_header <- function(x) {
    if (is.na(x)) return(NA_character_)
    stringr::str_squish(as.character(x))
  }

  expect_cell <- function(row_idx, col_idx, expected) {
    observed <- gps_raw[[col_idx]][row_idx]
    observed_chr <- normalize_header(observed)
    expected_chr <- normalize_header(expected)
    if (!identical(observed_chr, expected_chr)) {
      stop(
        sprintf(
          "GPS workbook layout mismatch at row %d, col %d: expected '%s' but found '%s'",
          row_idx, col_idx, expected_chr, observed_chr %||% "NA"
        ),
        call. = FALSE
      )
    }
  }

  # Explicitly validate the two-row workbook header so column-position ingest
  # cannot silently drift if the Excel layout changes.
  expect_cell(1, 2, "Pre turmeric (14-week pre-season up to Visit 2)")
  expect_cell(1, 7, "Post turmeric period (visit 3) (24-week season up to visit 3)")
  expect_cell(1, 25, "Days absent injury")
  expect_cell(1, 29, "Days absent illness")

  expect_cell(2, 2, "Training ")
  expect_cell(2, 3, "Training minutes per week")
  expect_cell(2, 4, "Match")
  expect_cell(2, 5, "Match minutes per week")
  expect_cell(2, 7, "Training Total")
  expect_cell(2, 8, "Training minutes per week")
  expect_cell(2, 9, "Match Total")
  expect_cell(2, 10, "Match minutes per week")
  expect_cell(2, 25, "pre turmeric")
  expect_cell(2, 26, "% days absent")
  expect_cell(2, 27, "post turmeric")
  expect_cell(2, 28, "% days absent")
  expect_cell(2, 29, "pre turmeric")
  expect_cell(2, 30, "% days absent")
  expect_cell(2, 31, "post turmeric")
  expect_cell(2, 32, "% days absent")

  invisible(TRUE)
}

extract_sample_metadata <- function(sample_cols) {
  base <- sub("_L001.*$", "", sample_cols)
  m <- stringr::str_match(base, "^Nottingham_([0-9]+)([A-Za-z])_S\\d+$")
  tibble::tibble(
    sample_col = sample_cols,
    sample_id = base,
    participant_id = suppressWarnings(as.integer(m[, 2])),
    block_letter = toupper(m[, 3])
  ) %>%
    mutate(
      timepoint = case_when(
        block_letter == "A" ~ "Baseline",
        block_letter == "B" ~ "2_weeks",
        block_letter == "C" ~ "6_months",
        TRUE ~ NA_character_
      )
    )
}

compute_feature_delta <- function(long_df) {
  long_df %>%
    group_by(participant_id, outcome, level) %>%
    summarise(
      pre_value = ifelse(
        any(timepoint %in% c("Baseline", "2_weeks") & !is.na(value)),
        mean(value[timepoint %in% c("Baseline", "2_weeks")], na.rm = TRUE),
        NA_real_
      ),
      post_value = ifelse(
        any(timepoint == "6_months" & !is.na(value)),
        value[which(timepoint == "6_months")[1]],
        NA_real_
      ),
      delta = post_value - pre_value,
      .groups = "drop"
    )
}

run_association_grid <- function(outcome_df, predictor_df, predictor_cols, family_label) {
  all_outcomes <- unique(outcome_df$outcome)
  res <- purrr::map_dfr(all_outcomes, function(outc) {
    purrr::map_dfr(predictor_cols, function(pred) {
      d <- outcome_df %>%
        filter(outcome == outc) %>%
        select(participant_id, level, outcome, delta) %>%
        left_join(predictor_df %>% select(participant_id, all_of(pred)), by = "participant_id") %>%
        rename(predictor_value = all_of(pred)) %>%
        filter(is.finite(delta), is.finite(predictor_value))

      n <- nrow(d)
      if (n < MIN_N_ASSOC) {
        return(tibble::tibble(
          family = family_label,
          level = unique(d$level)[1] %||% NA_character_,
          outcome = outc,
          predictor = pred,
          n = n,
          rho = NA_real_,
          p_value = NA_real_
        ))
      }

      if (sd(d$delta, na.rm = TRUE) == 0 || sd(d$predictor_value, na.rm = TRUE) == 0) {
        return(tibble::tibble(
          family = family_label,
          level = unique(d$level)[1] %||% NA_character_,
          outcome = outc,
          predictor = pred,
          n = n,
          rho = NA_real_,
          p_value = NA_real_
        ))
      }

      ct <- suppressWarnings(cor.test(d$delta, d$predictor_value, method = "spearman", exact = FALSE))
      tibble::tibble(
        family = family_label,
        level = unique(d$level)[1] %||% NA_character_,
        outcome = outc,
        predictor = pred,
        n = n,
        rho = as.numeric(ct$estimate),
        p_value = ct$p.value
      )
    })
  })

  res %>%
    mutate(
      q_value_bh = p.adjust(p_value, method = "BH"),
      abs_rho = abs(rho)
    ) %>%
    arrange(p_value, desc(abs_rho))
}

run_association_grid_adjusted <- function(outcome_df, predictor_df, predictor_cols, family_label) {
  all_outcomes <- unique(outcome_df$outcome)
  res <- purrr::map_dfr(all_outcomes, function(outc) {
    purrr::map_dfr(predictor_cols, function(pred) {
      d <- outcome_df %>%
        filter(outcome == outc) %>%
        select(participant_id, level, outcome, delta) %>%
        left_join(predictor_df %>% select(participant_id, turmeric_group, all_of(pred)), by = "participant_id") %>%
        rename(predictor_value = all_of(pred)) %>%
        filter(is.finite(delta), is.finite(predictor_value), !is.na(turmeric_group))

      n <- nrow(d)
      if (n < MIN_N_ASSOC) {
        return(tibble::tibble(
          family = family_label,
          level = unique(d$level)[1] %||% NA_character_,
          outcome = outc,
          predictor = pred,
          n = n,
          rho = NA_real_,
          p_value = NA_real_
        ))
      }

      if (length(unique(d$turmeric_group)) < 2) {
        return(tibble::tibble(
          family = family_label,
          level = unique(d$level)[1] %||% NA_character_,
          outcome = outc,
          predictor = pred,
          n = n,
          rho = NA_real_,
          p_value = NA_real_
        ))
      }

      # Partial Spearman approach:
      # rank-transform both variables, regress each on turmeric_group,
      # then correlate residuals (Pearson on residualized ranks).
      ry <- rank(d$delta, ties.method = "average")
      rx <- rank(d$predictor_value, ties.method = "average")

      fit_y <- lm(ry ~ turmeric_group, data = d)
      fit_x <- lm(rx ~ turmeric_group, data = d)
      ey <- residuals(fit_y)
      ex <- residuals(fit_x)

      if (sd(ey, na.rm = TRUE) == 0 || sd(ex, na.rm = TRUE) == 0) {
        return(tibble::tibble(
          family = family_label,
          level = unique(d$level)[1] %||% NA_character_,
          outcome = outc,
          predictor = pred,
          n = n,
          rho = NA_real_,
          p_value = NA_real_
        ))
      }

      ct <- suppressWarnings(cor.test(ey, ex, method = "pearson"))
      tibble::tibble(
        family = family_label,
        level = unique(d$level)[1] %||% NA_character_,
        outcome = outc,
        predictor = pred,
        n = n,
        rho = as.numeric(ct$estimate),
        p_value = ct$p.value
      )
    })
  })

  res %>%
    mutate(
      q_value_bh = p.adjust(p_value, method = "BH"),
      abs_rho = abs(rho)
    ) %>%
    arrange(p_value, desc(abs_rho))
}

fit_bayes_rank_regression <- function(delta, predictor, turmeric_group = NULL, adjusted = FALSE,
                                      prior_var = 100, a0 = 0.5, b0 = 0.5) {
  y <- as.numeric(scale(rank(delta, ties.method = "average")))
  x <- as.numeric(scale(rank(predictor, ties.method = "average")))

  if (!all(is.finite(y)) || !all(is.finite(x)) || sd(y) == 0 || sd(x) == 0) {
    return(NULL)
  }

  if (adjusted) {
    d <- tibble::tibble(y = y, x = x, turmeric_group = as.factor(turmeric_group))
    if (length(unique(d$turmeric_group)) < 2) return(NULL)
    X <- model.matrix(~ x + turmeric_group, data = d)
  } else {
    d <- tibble::tibble(y = y, x = x)
    X <- model.matrix(~ x, data = d)
  }

  p <- ncol(X)
  n <- nrow(X)
  if (n <= p) return(NULL)

  coef_idx <- which(colnames(X) == "x")
  if (length(coef_idx) != 1) return(NULL)

  # Conjugate Bayesian linear model:
  # beta | sigma^2 ~ N(0, sigma^2 * prior_var * I), sigma^2 ~ IG(a0, b0)
  V0_inv <- diag(1 / prior_var, p)
  XtX <- crossprod(X)
  Vn_inv <- V0_inv + XtX
  Vn <- tryCatch(solve(Vn_inv), error = function(e) NULL)
  if (is.null(Vn)) return(NULL)

  m_n <- Vn %*% crossprod(X, y)
  a_n <- a0 + n / 2
  b_n <- b0 + 0.5 * as.numeric(crossprod(y) - crossprod(m_n, Vn_inv %*% m_n))
  if (!is.finite(b_n) || b_n <= 0) return(NULL)

  mean_beta <- as.numeric(m_n[coef_idx, 1])
  var_beta <- as.numeric((b_n / a_n) * Vn[coef_idx, coef_idx])
  if (!is.finite(var_beta) || var_beta <= 0) return(NULL)

  sd_beta <- sqrt(var_beta)
  df <- 2 * a_n

  ci_low <- mean_beta + qt(0.025, df = df) * sd_beta
  ci_high <- mean_beta + qt(0.975, df = df) * sd_beta
  p_beta_gt0 <- 1 - pt((0 - mean_beta) / sd_beta, df = df)
  p_beta_lt0 <- 1 - p_beta_gt0
  p_direction <- max(p_beta_gt0, p_beta_lt0)
  p_tail_two_sided <- 2 * min(p_beta_gt0, p_beta_lt0)

  tibble::tibble(
    beta_rank = mean_beta,
    ci_low = ci_low,
    ci_high = ci_high,
    p_beta_gt0 = p_beta_gt0,
    p_beta_lt0 = p_beta_lt0,
    p_direction = p_direction,
    p_tail_two_sided = p_tail_two_sided
  )
}

run_bayesian_grid <- function(outcome_df, predictor_df, predictor_cols, family_label, adjusted = FALSE) {
  all_outcomes <- unique(outcome_df$outcome)
  res <- purrr::map_dfr(all_outcomes, function(outc) {
    purrr::map_dfr(predictor_cols, function(pred) {
      join_cols <- c("participant_id", pred)
      if (adjusted) join_cols <- c("participant_id", "turmeric_group", pred)

      d <- outcome_df %>%
        filter(outcome == outc) %>%
        select(participant_id, level, outcome, delta) %>%
        left_join(predictor_df %>% select(all_of(join_cols)), by = "participant_id") %>%
        rename(predictor_value = all_of(pred)) %>%
        filter(is.finite(delta), is.finite(predictor_value))

      if (adjusted) d <- d %>% filter(!is.na(turmeric_group))

      n <- nrow(d)
      if (n < MIN_N_ASSOC) {
        return(tibble::tibble(
          family = family_label,
          level = unique(d$level)[1] %||% NA_character_,
          outcome = outc,
          predictor = pred,
          n = n,
          beta_rank = NA_real_,
          ci_low = NA_real_,
          ci_high = NA_real_,
          p_beta_gt0 = NA_real_,
          p_beta_lt0 = NA_real_,
          p_direction = NA_real_,
          p_tail_two_sided = NA_real_
        ))
      }

      fit <- fit_bayes_rank_regression(
        delta = d$delta,
        predictor = d$predictor_value,
        turmeric_group = if (adjusted) d$turmeric_group else NULL,
        adjusted = adjusted
      )

      if (is.null(fit)) {
        return(tibble::tibble(
          family = family_label,
          level = unique(d$level)[1] %||% NA_character_,
          outcome = outc,
          predictor = pred,
          n = n,
          beta_rank = NA_real_,
          ci_low = NA_real_,
          ci_high = NA_real_,
          p_beta_gt0 = NA_real_,
          p_beta_lt0 = NA_real_,
          p_direction = NA_real_,
          p_tail_two_sided = NA_real_
        ))
      }

      tibble::tibble(
        family = family_label,
        level = unique(d$level)[1] %||% NA_character_,
        outcome = outc,
        predictor = pred,
        n = n
      ) %>%
        bind_cols(fit)
    })
  })

  res %>%
    mutate(
      q_tail_bh = p.adjust(p_tail_two_sided, method = "BH"),
      abs_beta_rank = abs(beta_rank)
    ) %>%
    arrange(p_tail_two_sided, desc(abs_beta_rank))
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

load_gps_predictors <- function() {
  if (nzchar(GPS_FILE) && file.exists(GPS_FILE)) {
    gps_raw <- readxl::read_excel(GPS_FILE, sheet = 1, col_names = FALSE)
    validate_gps_layout(gps_raw)

    return(
      gps_raw %>%
        transmute(
          participant_id = parse_num(...1),
          pre_training_total = parse_num(...2),
          pre_training_min_week = parse_num(...3),
          pre_match_total = parse_num(...4),
          pre_match_min_week = parse_num(...5),
          post_training_total = parse_num(...7),
          post_training_min_week = parse_num(...8),
          post_match_total = parse_num(...9),
          post_match_min_week = parse_num(...10),
          injury_days_pre = parse_num(...25),
          injury_pct_pre = parse_num(...26),
          injury_days_post = parse_num(...27),
          injury_pct_post = parse_num(...28),
          illness_days_pre = parse_num(...29),
          illness_pct_pre = parse_num(...30),
          illness_days_post = parse_num(...31),
          illness_pct_post = parse_num(...32)
        )
    )
  }

  if (file.exists(GPS_PREPARED_FILE)) {
    message("GPS workbook not found; using committed prepared cohort table: ", GPS_PREPARED_FILE)
    return(readr::read_csv(GPS_PREPARED_FILE, show_col_types = FALSE))
  }

  stop(
    "GPS/injury input not found. Set GPS_INJURY_XLSX to the source workbook or provide ",
    GPS_PREPARED_FILE,
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 1) Clean GPS + injury workbook (participant-level)
# -----------------------------------------------------------------------------
gps <- load_gps_predictors() %>%
  filter(!is.na(participant_id)) %>%
  mutate(participant_id = as.integer(participant_id)) %>%
  left_join(turmeric_status, by = "participant_id", suffix = c("", ".status")) %>%
  mutate(
    status = coalesce(status, status.status),
    turmeric_group = coalesce(turmeric_group, turmeric_group.status)
  ) %>%
  select(-matches("\\.status$")) %>%
  mutate(
    delta_training_total = post_training_total - pre_training_total,
    delta_training_min_week = post_training_min_week - pre_training_min_week,
    delta_match_total = post_match_total - pre_match_total,
    delta_match_min_week = post_match_min_week - pre_match_min_week,
    injury_days_total = injury_days_pre + injury_days_post,
    illness_days_total = illness_days_pre + illness_days_post
  ) %>%
  arrange(participant_id)

gps_analysis <- gps %>% filter(participant_id %in% analysis_ids)

predictor_cols <- c(
  "pre_training_total", "pre_training_min_week", "pre_match_total", "pre_match_min_week",
  "post_training_total", "post_training_min_week", "post_match_total", "post_match_min_week",
  "delta_training_total", "delta_training_min_week", "delta_match_total", "delta_match_min_week",
  "injury_days_post", "illness_days_post", "injury_days_total", "illness_days_total",
  "injury_pct_post", "illness_pct_post"
)

predictor_cols <- predictor_cols[predictor_cols %in% names(gps_analysis)]

# -----------------------------------------------------------------------------
# 2) Alpha-diversity deltas (participant-level)
# -----------------------------------------------------------------------------
alpha <- read.csv(ALPHA_FILE, check.names = FALSE)

alpha_long <- alpha %>%
  select(sample_id, Richness, Shannon, Simpson) %>%
  pivot_longer(cols = c(Richness, Shannon, Simpson), names_to = "outcome", values_to = "value") %>%
  left_join(extract_sample_metadata(unique(alpha$sample_id)), by = c("sample_id")) %>%
  filter(participant_id %in% analysis_ids, !is.na(timepoint)) %>%
  mutate(level = "alpha_diversity")

alpha_delta <- compute_feature_delta(alpha_long)

# -----------------------------------------------------------------------------
# 3) Strain/genus/pathway signal deltas carried forward from 4e
# -----------------------------------------------------------------------------
strain_eff <- read.csv(STRAIN_EFF_FILE, check.names = FALSE)
genus_eff <- read.csv(GENUS_EFF_FILE, check.names = FALSE)
pathway_eff <- read.csv(PATHWAY_EFF_FILE, check.names = FALSE)

all_strains <- strain_eff %>%
  arrange(rank_consensus, desc(abs(cohens_d))) %>%
  select(feature, feature_short)

all_genera <- genus_eff %>%
  arrange(rank_consensus, desc(abs(cohens_d))) %>%
  select(feature, feature_short)

all_pathways <- pathway_eff %>%
  arrange(rank_consensus, desc(abs(cohens_d))) %>%
  select(feature, feature_short)

# Strains from species table
species <- read.csv(SPECIES_FILE, check.names = FALSE)
species_sample_cols <- setdiff(names(species), "Taxa")
species_meta <- extract_sample_metadata(species_sample_cols)

strain_long <- species %>%
  filter(Taxa %in% all_strains$feature) %>%
  rename(feature = Taxa) %>%
  pivot_longer(cols = all_of(species_sample_cols), names_to = "sample_col", values_to = "value") %>%
  left_join(species_meta, by = "sample_col") %>%
  filter(participant_id %in% analysis_ids, !is.na(timepoint)) %>%
  left_join(all_strains, by = c("feature")) %>%
  mutate(
    outcome = paste0("strain_", feature_short),
    level = "strain"
  )

strain_delta <- compute_feature_delta(strain_long)

# Genera aggregated from species table
genus_abund <- species %>%
  mutate(genus = str_extract(Taxa, "(?<=\\|g__)[^|]+")) %>%
  filter(!is.na(genus), genus %in% all_genera$feature) %>%
  select(genus, all_of(species_sample_cols)) %>%
  group_by(genus) %>%
  summarise(across(all_of(species_sample_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  rename(feature = genus)

genus_long <- genus_abund %>%
  pivot_longer(cols = all_of(species_sample_cols), names_to = "sample_col", values_to = "value") %>%
  left_join(species_meta, by = "sample_col") %>%
  filter(participant_id %in% analysis_ids, !is.na(timepoint)) %>%
  left_join(all_genera, by = "feature") %>%
  mutate(
    outcome = paste0("genus_", feature_short),
    level = "genus"
  )

genus_delta <- compute_feature_delta(genus_long)

# Pathways from HUMAnN table
pathways <- read.delim(PATHWAY_FILE, check.names = FALSE)
path_feature_col <- names(pathways)[1]
path_sample_cols <- setdiff(names(pathways), path_feature_col)
path_meta <- extract_sample_metadata(path_sample_cols)

pathway_long <- pathways %>%
  rename(feature = !!path_feature_col) %>%
  filter(feature %in% all_pathways$feature) %>%
  pivot_longer(cols = all_of(path_sample_cols), names_to = "sample_col", values_to = "value") %>%
  left_join(path_meta, by = "sample_col") %>%
  filter(participant_id %in% analysis_ids, !is.na(timepoint)) %>%
  left_join(all_pathways, by = "feature") %>%
  mutate(
    outcome = paste0("pathway_", feature_short),
    level = "pathway"
  )

pathway_delta <- compute_feature_delta(pathway_long)

signal_delta <- bind_rows(strain_delta, genus_delta, pathway_delta)

# -----------------------------------------------------------------------------
# 4) Association testing
# -----------------------------------------------------------------------------
alpha_assoc <- run_association_grid(
  outcome_df = alpha_delta,
  predictor_df = gps_analysis,
  predictor_cols = predictor_cols,
  family_label = "alpha"
)

signal_assoc <- run_association_grid(
  outcome_df = signal_delta,
  predictor_df = gps_analysis,
  predictor_cols = predictor_cols,
  family_label = "signals"
)

alpha_assoc_adj <- run_association_grid_adjusted(
  outcome_df = alpha_delta,
  predictor_df = gps_analysis,
  predictor_cols = predictor_cols,
  family_label = "alpha"
)

signal_assoc_adj <- run_association_grid_adjusted(
  outcome_df = signal_delta,
  predictor_df = gps_analysis,
  predictor_cols = predictor_cols,
  family_label = "signals"
)

alpha_bayes <- run_bayesian_grid(
  outcome_df = alpha_delta,
  predictor_df = gps_analysis,
  predictor_cols = predictor_cols,
  family_label = "alpha",
  adjusted = FALSE
)

signal_bayes <- run_bayesian_grid(
  outcome_df = signal_delta,
  predictor_df = gps_analysis,
  predictor_cols = predictor_cols,
  family_label = "signals",
  adjusted = FALSE
)

alpha_bayes_adj <- run_bayesian_grid(
  outcome_df = alpha_delta,
  predictor_df = gps_analysis,
  predictor_cols = predictor_cols,
  family_label = "alpha",
  adjusted = TRUE
)

signal_bayes_adj <- run_bayesian_grid(
  outcome_df = signal_delta,
  predictor_df = gps_analysis,
  predictor_cols = predictor_cols,
  family_label = "signals",
  adjusted = TRUE
)

# -----------------------------------------------------------------------------
# 5) Save outputs
# -----------------------------------------------------------------------------
write.csv(gps, file.path(OUT_DIR, "gps_injury_clean_all_participants.csv"), row.names = FALSE)
write.csv(gps_analysis, file.path(OUT_DIR, "gps_injury_analysis_cohort_n11.csv"), row.names = FALSE)

write.csv(alpha_delta, file.path(OUT_DIR, "participant_alpha_deltas_n11.csv"), row.names = FALSE)
write.csv(signal_delta, file.path(OUT_DIR, "participant_signal_deltas_n11.csv"), row.names = FALSE)

write.csv(alpha_assoc, file.path(OUT_DIR, "associations_alpha_vs_gps_injury_spearman.csv"), row.names = FALSE)
write.csv(signal_assoc, file.path(OUT_DIR, "associations_signals_vs_gps_injury_spearman.csv"), row.names = FALSE)
write.csv(alpha_assoc_adj, file.path(OUT_DIR, "associations_alpha_vs_gps_injury_partial_spearman_adjusted_turmeric_group.csv"), row.names = FALSE)
write.csv(signal_assoc_adj, file.path(OUT_DIR, "associations_signals_vs_gps_injury_partial_spearman_adjusted_turmeric_group.csv"), row.names = FALSE)
write.csv(alpha_bayes, file.path(OUT_DIR, "bayes_associations_alpha_vs_gps_injury_rank_regression.csv"), row.names = FALSE)
write.csv(signal_bayes, file.path(OUT_DIR, "bayes_associations_signals_vs_gps_injury_rank_regression.csv"), row.names = FALSE)
write.csv(alpha_bayes_adj, file.path(OUT_DIR, "bayes_associations_alpha_vs_gps_injury_rank_regression_adjusted_turmeric_group.csv"), row.names = FALSE)
write.csv(signal_bayes_adj, file.path(OUT_DIR, "bayes_associations_signals_vs_gps_injury_rank_regression_adjusted_turmeric_group.csv"), row.names = FALSE)

top_alpha <- alpha_assoc %>%
  filter(!is.na(p_value)) %>%
  arrange(p_value) %>%
  slice_head(n = 20)

top_signals <- signal_assoc %>%
  filter(!is.na(p_value)) %>%
  arrange(p_value) %>%
  slice_head(n = 30)

top_alpha_adj <- alpha_assoc_adj %>%
  filter(!is.na(p_value)) %>%
  arrange(p_value) %>%
  slice_head(n = 20)

top_signals_adj <- signal_assoc_adj %>%
  filter(!is.na(p_value)) %>%
  arrange(p_value) %>%
  slice_head(n = 30)

top_alpha_bayes <- alpha_bayes %>%
  filter(!is.na(p_tail_two_sided)) %>%
  arrange(p_tail_two_sided) %>%
  slice_head(n = 20)

top_signals_bayes <- signal_bayes %>%
  filter(!is.na(p_tail_two_sided)) %>%
  arrange(p_tail_two_sided) %>%
  slice_head(n = 30)

top_alpha_bayes_adj <- alpha_bayes_adj %>%
  filter(!is.na(p_tail_two_sided)) %>%
  arrange(p_tail_two_sided) %>%
  slice_head(n = 20)

top_signals_bayes_adj <- signal_bayes_adj %>%
  filter(!is.na(p_tail_two_sided)) %>%
  arrange(p_tail_two_sided) %>%
  slice_head(n = 30)

summary_lines <- c(
  "# GPS/Injury vs Microbiome Association Summary",
  "",
  glue::glue("Analysis date: {Sys.Date()}"),
  glue::glue("GPS workbook: {GPS_FILE}"),
  glue::glue("Analysis cohort (aligned with microbiome treatment/control pipeline): n = {nrow(gps_analysis)} participants"),
  "",
  "## Predictors tested",
  paste0("- ", predictor_cols),
  "",
  "## Alpha-diversity outcomes",
  glue::glue("- Outcomes: {paste(unique(alpha_delta$outcome), collapse=', ')}"),
  glue::glue("- Tests run: {sum(!is.na(alpha_assoc$p_value))}"),
  glue::glue("- Nominal p < 0.05: {sum(alpha_assoc$p_value < 0.05, na.rm = TRUE)}"),
  glue::glue("- BH q < 0.10: {sum(alpha_assoc$q_value_bh < 0.10, na.rm = TRUE)}"),
  glue::glue("- Adjusted (partial Spearman controlling for turmeric_group): p < 0.05 = {sum(alpha_assoc_adj$p_value < 0.05, na.rm = TRUE)}, q < 0.10 = {sum(alpha_assoc_adj$q_value_bh < 0.10, na.rm = TRUE)}"),
  glue::glue("- Bayesian rank-regression (unadjusted): posterior direction >= 0.95 = {sum(alpha_bayes$p_direction >= 0.95, na.rm = TRUE)}, BH q_tail < 0.10 = {sum(alpha_bayes$q_tail_bh < 0.10, na.rm = TRUE)}"),
  glue::glue("- Bayesian rank-regression (adjusted for turmeric_group): posterior direction >= 0.95 = {sum(alpha_bayes_adj$p_direction >= 0.95, na.rm = TRUE)}, BH q_tail < 0.10 = {sum(alpha_bayes_adj$q_tail_bh < 0.10, na.rm = TRUE)}"),
  "",
  "## Signal outcomes (all retained strain/genus/pathway features from 4e)",
  glue::glue("- Strains: {nrow(all_strains)} | Genera: {nrow(all_genera)} | Pathways: {nrow(all_pathways)}"),
  glue::glue("- Tests run: {sum(!is.na(signal_assoc$p_value))}"),
  glue::glue("- Nominal p < 0.05: {sum(signal_assoc$p_value < 0.05, na.rm = TRUE)}"),
  glue::glue("- BH q < 0.10: {sum(signal_assoc$q_value_bh < 0.10, na.rm = TRUE)}"),
  glue::glue("- Adjusted (partial Spearman controlling for turmeric_group): p < 0.05 = {sum(signal_assoc_adj$p_value < 0.05, na.rm = TRUE)}, q < 0.10 = {sum(signal_assoc_adj$q_value_bh < 0.10, na.rm = TRUE)}"),
  glue::glue("- Bayesian rank-regression (unadjusted): posterior direction >= 0.95 = {sum(signal_bayes$p_direction >= 0.95, na.rm = TRUE)}, BH q_tail < 0.10 = {sum(signal_bayes$q_tail_bh < 0.10, na.rm = TRUE)}"),
  glue::glue("- Bayesian rank-regression (adjusted for turmeric_group): posterior direction >= 0.95 = {sum(signal_bayes_adj$p_direction >= 0.95, na.rm = TRUE)}, BH q_tail < 0.10 = {sum(signal_bayes_adj$q_tail_bh < 0.10, na.rm = TRUE)}"),
  "",
  "## Top alpha associations (by p-value)",
  paste(capture.output(print(top_alpha, n = 20, width = 140)), collapse = "\n"),
  "",
  "## Top alpha associations adjusted for turmeric_group (partial Spearman)",
  paste(capture.output(print(top_alpha_adj, n = 20, width = 140)), collapse = "\n"),
  "",
  "## Top signal associations (by p-value)",
  paste(capture.output(print(top_signals, n = 30, width = 160)), collapse = "\n"),
  "",
  "## Top signal associations adjusted for turmeric_group (partial Spearman)",
  paste(capture.output(print(top_signals_adj, n = 30, width = 160)), collapse = "\n"),
  "",
  "## Top alpha Bayesian rank-regression associations (by posterior tail probability)",
  paste(capture.output(print(top_alpha_bayes, n = 20, width = 160)), collapse = "\n"),
  "",
  "## Top alpha Bayesian rank-regression associations adjusted for turmeric_group",
  paste(capture.output(print(top_alpha_bayes_adj, n = 20, width = 160)), collapse = "\n"),
  "",
  "## Top signal Bayesian rank-regression associations (by posterior tail probability)",
  paste(capture.output(print(top_signals_bayes, n = 30, width = 180)), collapse = "\n"),
  "",
  "## Top signal Bayesian rank-regression associations adjusted for turmeric_group",
  paste(capture.output(print(top_signals_bayes_adj, n = 30, width = 180)), collapse = "\n")
)

writeLines(summary_lines, file.path(OUT_DIR, "association_summary.txt"))

cat("Association analysis complete.\n")
cat("Output directory:", OUT_DIR, "\n")
