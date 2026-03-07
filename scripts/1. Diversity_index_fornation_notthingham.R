

# ────────────────────────────────────────────────────────────────────────────────
# Script: 05_diversity_indices_calculation.R
# Description: Calculate alpha and beta diversity indices for species and 
#              functional profiles from preprocessed microbiome data
# Input: Raw data from 00_raw_data
# Output: Diversity metrics to 04_diversity_analysis
#         Figures to figures/04_diversity_analysis
# Author: Liam Walsh
# Date: 2025
# ────────────────────────────────────────────────────────────────────────────────

# ────────────────────────────────────────────────────────────────────────────────
# Section: Reproducibility
# ────────────────────────────────────────────────────────────────────────────────
# Set custom library path and seed
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
# ────────────────────────────────────────────────────────────────────────────────
# Section: Load Required Packages
# ────────────────────────────────────────────────────────────────────────────────
pacman::p_load(gh, base64enc, readr, readxl, httr, glue, jsonlite, hillR, dplyr, 
               tidyr, textshape, tibble, stringr, data.table)

# ────────────────────────────────────────────────────────────────────────────────
# Section: GitHub Configuration
# ────────────────────────────────────────────────────────────────────────────────
owner <- "LiamHWalsh"
repo <- "NTU_Teagasc_Turmeric_study"
branch <- "main"
token <- Sys.getenv("GITHUB_PAT")
if (token == "") token <- Sys.getenv("GITHUB_TOKEN")
enable_github_upload <- identical(Sys.getenv("ENABLE_GITHUB_UPLOAD"), "1")
allow_github_write <- nzchar(token) && enable_github_upload
# Define pipeline structure - each step's output becomes next step's input
pipeline_paths <- list(
  # Step 0: Raw data (input for this script)
  raw_data = "results/Raw",
  
  # Step 1: Quality control and preprocessing
  preprocessing = "results/Processed/01_preprocessing",
  
  # Step 2: Taxonomic profiling (species-level)
  taxonomic = "results/Processed/02_taxonomic_profiling",
  
  # Step 3: Functional profiling (pathways)
  functional = "results/Processed/03_functional_profiling",
  
  # Step 4: Diversity analysis (alpha/beta)
  diversity = "results/Processed/04_diversity_analysis",
  
  # Step 5: Statistical analysis
  statistics = "results/Processed/05_statistical_analysis",
  
  # Step 6: Final reports
  reports = "results/Processed/06_reports"
)

# Define parallel figure structure - mirrors results structure
figure_paths <- list(
  # Figures organized by analysis step
  preprocessing = "figures/01_preprocessing",
  taxonomic = "figures/02_taxonomic_profiling",
  functional = "figures/03_functional_profiling",
  diversity = "figures/04_diversity_analysis",
  statistics = "figures/05_statistical_analysis",
  reports = "figures/06_reports"
)

# Set current script's input and output paths
input_path <- pipeline_paths$raw_data
output_path <- pipeline_paths$diversity  # This script outputs diversity analysis
figure_output_path <- figure_paths$diversity  # Figures go here

if (nzchar(token)) {
  message(glue::glue("GitHub token detected (length: {nchar(token)} characters)"))
} else {
  message("No GitHub token detected. Running in read-only / local-first mode.")
}
if (!allow_github_write) {
  message("GitHub upload disabled. Set ENABLE_GITHUB_UPLOAD=1 with a valid token to upload outputs.")
}
message(glue::glue("Input path: {input_path}"))

# Check if input directory exists, if not try alternative paths
message("🔍 Checking for available data directories...")

# List of possible input directories to check
possible_paths <- c(
  "results/00_raw_data",
  "results/Raw",
  "results/raw",
  "results/raw_data",
  "Raw",
  "raw"
)

# Function to check if directory exists on GitHub
check_github_directory <- function(path) {
  api_url <- glue::glue("https://api.github.com/repos/{owner}/{repo}/contents/{path}")
  response <- if (nzchar(token)) {
    GET(api_url, authenticate(token, ""), add_headers("Accept" = "application/vnd.github.v3+json"))
  } else {
    GET(api_url, add_headers("Accept" = "application/vnd.github.v3+json"))
  }
  return(status_code(response) == 200)
}

# Try to find existing directory
found_path <- NULL
for (path in possible_paths) {
  if (check_github_directory(path)) {
    found_path <- path
    message(glue::glue("✅ Found data directory: {path}"))
    break
  }
}

if (is.null(found_path)) {
  stop(glue::glue("
  ❌ No raw data directory found in the repository.
  
  Tried the following paths:
  {paste('  -', possible_paths, collapse = '\n')}
  
  Please either:
  1. Upload your raw data files to one of these directories, OR
  2. Update the 'input_path' variable to match your actual directory name
  
  Current repository structure can be viewed at:
  https://github.com/{owner}/{repo}
  "))
}

# Use the found path
input_path <- found_path
message(glue::glue("📤 Output path (data): {output_path}"))
message(glue::glue("📊 Output path (figures): {figure_output_path}"))

# ────────────────────────────────────────────────────────────────────────────────
# Section: GitHub Upload Function (Defined Early)
# ────────────────────────────────────────────────────────────────────────────────
save_and_upload_to_github <- function(
    object = NULL,                        
    object_type = c("data", "plot", "file"), 
    local_path = NULL,                   
    github_path,                         
    commit_message,                      
    owner, repo, branch = "main"         
) {
  
  # Validate object type
  object_type <- match.arg(object_type)
  
  # Use a temporary file path if none provided (for data/plot)
  if (is.null(local_path)) {
    file_ext <- switch(object_type,
                       data = ".csv",
                       plot = ".png",
                       stop("For 'file' type, 'local_path' must be provided.")
    )
    local_path <- tempfile(fileext = file_ext)
    message(glue::glue("ℹ️  No local_path provided. Using temporary file: {local_path}"))
  }
  
  # Save the object locally (temporary)
  if (object_type == "data") {
    if (is.null(object)) stop("❌ No data object provided.")
    data.table::fwrite(object, local_path, row.names = FALSE)
    message(glue::glue("✅ Data saved temporarily to: {local_path}"))
    
  } else if (object_type == "plot") {
    if (is.null(object)) stop("❌ No plot object provided.")
    ggplot2::ggsave(local_path, plot = object)
    message(glue::glue("✅ Plot saved temporarily to: {local_path}"))
    
  } else if (object_type == "file") {
    if (!file.exists(local_path)) stop(glue::glue("❌ File does not exist: {local_path}"))
    message(glue::glue("ℹ️  Using provided file: {local_path}"))
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

  # Encode file content to Base64 for GitHub API
  encoded_content <- base64enc::base64encode(local_path)
  
  # Get token from environment
  gh_token <- Sys.getenv("GITHUB_PAT")
  if (gh_token == "") gh_token <- Sys.getenv("GITHUB_TOKEN")
  
  # Check if the file already exists on GitHub
  file_info <- tryCatch(
    gh::gh("GET /repos/:owner/:repo/contents/:path",
           owner = owner, repo = repo, path = github_path, ref = branch,
           .token = gh_token),
    error = function(e) {
      if (grepl("404", e$message)) return(NULL)  # File doesn't exist yet
      stop(e)  # Other errors should be raised
    }
  )
  
  sha <- if (!is.null(file_info)) file_info$sha else NULL
  
  # Upload or update the file on GitHub
  gh::gh("PUT /repos/:owner/:repo/contents/:path",
         owner = owner,
         repo = repo,
         path = github_path,
         message = commit_message,
         content = encoded_content,
         branch = branch,
         sha = sha,
         .token = gh_token)
  
  # Confirmation message
  message(glue::glue("🚀 File uploaded to GitHub: https://github.com/{owner}/{repo}/blob/{branch}/{github_path}"))
  
  # Clean up temporary file if we created it
  if (object_type %in% c("data", "plot")) {
    unlink(local_path)
  }
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: List All Files in GitHub Directory
# ────────────────────────────────────────────────────────────────────────────────
message(glue::glue("📂 Fetching file list from: {owner}/{repo}/{input_path}"))

# Use GitHub API to list directory contents
api_url <- glue::glue("https://api.github.com/repos/{owner}/{repo}/contents/{input_path}")
response <- GET(api_url, 
                authenticate(token, ""),
                add_headers("Accept" = "application/vnd.github.v3+json"))

if (status_code(response) != 200) {
  stop(glue::glue("❌ Failed to access directory. Status code: {status_code(response)}
                   Check: 
                   1. Repository exists and is accessible
                   2. Directory path is correct: {input_path}
                   3. GitHub token has proper permissions"))
}

# Parse the response
content_list <- content(response, as = "text", encoding = "UTF-8")
files_info <- fromJSON(content_list)

# Filter for actual files (not directories) and exclude .docx files
file_data <- files_info[files_info$type == "file", ]
file_data <- file_data[!grepl("\\.docx$", file_data$name, ignore.case = TRUE), ]

message(glue::glue("✅ Found {nrow(file_data)} files to import (excluding .docx files)"))

# ────────────────────────────────────────────────────────────────────────────────
# Section: Download and Read Data Files from GitHub
# ────────────────────────────────────────────────────────────────────────────────
data <- list()  # Initialize an empty list to store data

for (i in 1:nrow(file_data)) {
  
  file_name <- file_data$name[i]
  file_url <- file_data$download_url[i]
  
  message(glue::glue("📥 Loading: {file_name}"))
  
  # Attempt to retrieve the file from GitHub
  response <- GET(file_url, authenticate(token, ""))
  
  if (status_code(response) == 200) {
    
    # Extract file extension
    file_ext <- tools::file_ext(file_name)
    
    # Create a clean name for the list (remove extension)
    clean_name <- gsub(paste0("\\.", file_ext, "$"), "", file_name)
    
    # Process based on file type
    tryCatch({
      if (file_ext == "csv") {
        file_content <- content(response, as = "text", encoding = "UTF-8")
        data[[clean_name]] <- read_csv(file_content, show_col_types = FALSE)
        message(glue::glue("   ✅ Loaded {file_name} ({nrow(data[[clean_name]])} rows)"))
        
      } else if (file_ext == "tsv") {
        file_content <- content(response, as = "text", encoding = "UTF-8")
        data[[clean_name]] <- read_delim(file_content, delim = "\t", 
                                         escape_double = FALSE, 
                                         trim_ws = TRUE, 
                                         show_col_types = FALSE)
        message(glue::glue("   ✅ Loaded {file_name} ({nrow(data[[clean_name]])} rows)"))
        
      } else if (file_ext %in% c("xls", "xlsx")) {
        # Download the file to a temporary location
        temp_file <- tempfile(fileext = paste0(".", file_ext))
        writeBin(content(response, as = "raw"), temp_file)
        data[[clean_name]] <- read_excel(temp_file)
        unlink(temp_file)  # Remove temporary file after reading
        message(glue::glue("   ✅ Loaded {file_name} ({nrow(data[[clean_name]])} rows)"))
        
      } else {
        warning(glue::glue("   ⚠️  Skipping {file_name} - unsupported file type"))
      }
      
    }, error = function(e) {
      warning(glue::glue("   ❌ Error reading {file_name}: {e$message}"))
    })
    
  } else {
    warning(glue::glue("   ❌ Unable to download {file_name} - Status: {status_code(response)}"))
  }
}

# Clean up metadata names for consistency
names(data)[grep("metaphlan", names(data), ignore.case = TRUE)] <- "metaphlan"

# ────────────────────────────────────────────────────────────────────────────────
# Section: Summary
# ────────────────────────────────────────────────────────────────────────────────
message("\n" , strrep("─", 80))
message(glue::glue("📊 Import Summary:"))
message(glue::glue("   Total files loaded: {length(data)}"))
message(glue::glue("   Available datasets: {paste(names(data), collapse = ', ')}"))
message(strrep("─", 80))

# ────────────────────────────────────────────────────────────────────────────────
# Section: Data Cleaning and Preparation
# ────────────────────────────────────────────────────────────────────────────────

# Rename preprocessing data
names(data)[which(names(data) == "total_fascinar_results")] <- "preprocessing"

# Rename pathway data
names(data)[which(names(data) == "HUMAnN_merged_pathabundance_cpm")] <- "pathways"

# Clean and harmonise preprocessing metadata
data[["preprocessing"]] <- data[["preprocessing"]] %>%
  separate(`Filename,#samplename`, into = c("X1", "X2", "X3", "sample_id_and_contig"), sep = "--", extra = "merge") %>%
  separate(sample_id_and_contig, into = c("sample_id", "file_name"), sep = ",", extra = "merge") %>%
  filter(sample_id == file_name) %>%
  select(-c("X1", "X2", "X3", "sample_id"))

# Prepare functional profile
data[["functional_profile"]] <- data$pathways %>%
  rename_with(~ gsub("_L001_concat_Abundance-CPM", "", .x))

data[["functional_profile"]] <- data[["functional_profile"]][-c(grep("\\|", data[["functional_profile"]]$`# Pathway`)), ]

# Convert to base data.frame so rownames are preserved
fp <- as.data.frame(data[["functional_profile"]])
rownames(fp) <- fp$`# Pathway`
fp$`# Pathway` <- NULL
data[["functional_profile"]] <- fp

# Standardise and filter species profile
data$metaphlan <- data$metaphlan %>% rename_with(~ gsub("_metaphlan", "", .))
data[["species_profile"]] <- data$metaphlan %>%
  column_to_rownames("clade_name") %>%
  t() %>%
  as.data.frame() %>%
  select(which(grepl("\\|t__", names(.)))) %>%
  t() %>%
  as.data.frame()

# ────────────────────────────────────────────────────────────────────────────────
# Section: Diversity Analysis Loop (Species & Functional)
# ────────────────────────────────────────────────────────────────────────────────

# Initialize list to store results
diversity_results <- list()

for (profile in c("species_profile", "functional_profile")) {
  
  # Filter: Keep features with max relative abundance > 0.1
  maxab <- apply(data[[profile]], 1, max, na.rm = TRUE)
  data[[profile]] <- data[[profile]][names(maxab[maxab > 0.1]), ]
  data[[profile]] <- data[[profile]][, colSums(data[[profile]]) > 0]
  
  # Calculate prevalence
  prevalence_df <- data.frame(
    sum_features = rowSums(data[[profile]] > 0.1)
  ) %>%
    rownames_to_column("feature") %>%
    filter(sum_features > ncol(data[[profile]]) * 0.1)
  
  # Add taxonomic breakdown for species profile
  if (profile == "species_profile") {
    prevalence_df <- prevalence_df %>%
      cbind(str_split_fixed(.$feature, "\\|", 8)) %>%
      rename(
        kingdom = `1`, phylum = `2`, class = `3`, order = `4`,
        family = `5`, genus = `6`, species = `7`, strain = `8`
      )
  }
  
  # Alpha diversity
  alpha_summary <- data.frame(
    Richness = hillR::hill_taxa(t(data[[profile]]), q = 0),
    Shannon = hillR::hill_taxa(t(data[[profile]]), q = 1),
    Simpson = hillR::hill_taxa(t(data[[profile]]), q = 2)
  ) %>%
    rownames_to_column("sample_id")
  
  # Beta diversity
  data_matrix <- as.matrix(t(data[[profile]]))
  beta_summary <- hill_taxa_parti_pairwise(
    data_matrix, q = 2, show_warning = FALSE, .progress = TRUE, rel_then_pool = FALSE
  )
  
  # Create distance matrix
  dist_matrix <- bind_rows(
    beta_summary %>% select(site1, site2, TD_beta),
    beta_summary %>% select(site2, site1, TD_beta) %>% rename(site1 = 1, site2 = 2)
  ) %>%
    pivot_wider(names_from = site2, values_from = TD_beta) %>%
    column_to_rownames("site1") %>%
    replace(is.na(.), 0)
  
  # Sample summary
  sample_summary <- data.frame(
    sample_id = colnames(data[[profile]]),
    total_abundance = colSums(data[[profile]]),
    richness = alpha_summary$Richness,
    shannon = alpha_summary$Shannon,
    simpson = alpha_summary$Simpson
  )
  
  # Store results
  diversity_results[[profile]] <- list(
    filtered_profile = data[[profile]],
    prevalence = prevalence_df,
    alpha_diversity = alpha_summary,
    beta_diversity = beta_summary,
    beta_distance_matrix = dist_matrix,
    sample_summary = sample_summary
  )
  
  message(glue::glue("✅ {profile} analysis complete:"))
  message(glue::glue("   - Features retained: {nrow(data[[profile]])}"))
  message(glue::glue("   - Samples: {ncol(data[[profile]])}"))
  message(glue::glue("   - Total abundance: {sum(sample_summary$total_abundance)}"))
}

# ────────────────────────────────────────────────────────────────────────────────
# Section: Upload Results to GitHub
# ────────────────────────────────────────────────────────────────────────────────

message("\n" , strrep("─", 80))
message("📤 Creating directory structure and uploading results to GitHub...")
message(strrep("─", 80))

# ────────────────────────────────────────────────────────────────────────────────
# Create directory structure with .gitkeep files
# ────────────────────────────────────────────────────────────────────────────────
message("📁 Ensuring all directories exist...")

# Combine all directory paths that need to be created
all_directories <- c(
  unlist(pipeline_paths),
  unlist(figure_paths)
)

# Create a .gitkeep file in each directory to ensure it exists
for (dir_path in all_directories) {
  tryCatch({
    # Create a temporary .gitkeep file
    temp_gitkeep <- tempfile()
    cat("# This file ensures the directory exists in Git\n", file = temp_gitkeep)
    cat("# Created by: 05_diversity_indices_calculation.R\n", file = temp_gitkeep, append = TRUE)
    cat(paste0("# Date: ", Sys.time(), "\n"), file = temp_gitkeep, append = TRUE)
    
    # Upload .gitkeep to create the directory
    save_and_upload_to_github(
      object_type = "file",
      local_path = temp_gitkeep,
      github_path = file.path(dir_path, ".gitkeep"),
      commit_message = glue::glue("Create directory: {dir_path}"),
      owner = owner, repo = repo, branch = branch
    )
    
    unlink(temp_gitkeep)
    message(glue::glue("   ✅ Created: {dir_path}"))
    
  }, error = function(e) {
    # Directory might already exist, which is fine
    if (grepl("already exists", e$message, ignore.case = TRUE)) {
      message(glue::glue("   ℹ️  Already exists: {dir_path}"))
    } else {
      warning(glue::glue("   ⚠️  Issue with {dir_path}: {e$message}"))
    }
  })
}

message("\n" , strrep("─", 80))
message("📤 Uploading analysis results...")
message(strrep("─", 80))

# 1. Upload preprocessing summary
save_and_upload_to_github(
  object = data[["preprocessing"]],
  object_type = "data",
  github_path = file.path(output_path, "sequencing_run_summary_statistics.csv"),
  commit_message = "Update sequencing run summary statistics",
  owner = owner, repo = repo, branch = branch
)

# 2. Upload unfiltered species profile
save_and_upload_to_github(
  object = data$metaphlan,
  object_type = "data",
  github_path = file.path(pipeline_paths$taxonomic, "species_profile_unfiltered.csv"),
  commit_message = "Update unfiltered species profile",
  owner = owner, repo = repo, branch = branch
)

# 3. Upload species-specific results
save_and_upload_to_github(
  object = diversity_results$species_profile$filtered_profile %>% rownames_to_column("Taxa"),
  object_type = "data",
  github_path = file.path(pipeline_paths$taxonomic, "species_profile_filtered.csv"),
  commit_message = "Update filtered species profile",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$species_profile$prevalence,
  object_type = "data",
  github_path = file.path(pipeline_paths$taxonomic, "species_prevalence_breakdown.csv"),
  commit_message = "Update species prevalence breakdown",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$species_profile$alpha_diversity,
  object_type = "data",
  github_path = file.path(output_path, "species_alpha_diversity.csv"),
  commit_message = "Update species alpha diversity",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$species_profile$beta_diversity,
  object_type = "data",
  github_path = file.path(output_path, "species_beta_diversity_TD_q2.csv"),
  commit_message = "Update species beta diversity",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$species_profile$beta_distance_matrix %>% rownames_to_column("sample_id"),
  object_type = "data",
  github_path = file.path(output_path, "species_beta_diversity_distance_matrix.csv"),
  commit_message = "Update species beta diversity distance matrix",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$species_profile$sample_summary,
  object_type = "data",
  github_path = file.path(output_path, "species_sample_summary.csv"),
  commit_message = "Update species sample summary",
  owner = owner, repo = repo, branch = branch
)

# 4. Upload functional profile-specific results
save_and_upload_to_github(
  object = diversity_results$functional_profile$filtered_profile %>% rownames_to_column("Pathway"),
  object_type = "data",
  github_path = file.path(pipeline_paths$functional, "functional_profile_filtered.csv"),
  commit_message = "Update filtered functional profile",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$functional_profile$prevalence,
  object_type = "data",
  github_path = file.path(pipeline_paths$functional, "functional_prevalence_breakdown.csv"),
  commit_message = "Update functional prevalence breakdown",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$functional_profile$alpha_diversity,
  object_type = "data",
  github_path = file.path(output_path, "functional_alpha_diversity.csv"),
  commit_message = "Update functional alpha diversity",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$functional_profile$beta_diversity,
  object_type = "data",
  github_path = file.path(output_path, "functional_beta_diversity_TD_q2.csv"),
  commit_message = "Update functional beta diversity",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$functional_profile$beta_distance_matrix %>% rownames_to_column("sample_id"),
  object_type = "data",
  github_path = file.path(output_path, "functional_beta_diversity_distance_matrix.csv"),
  commit_message = "Update functional beta diversity distance matrix",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$functional_profile$sample_summary,
  object_type = "data",
  github_path = file.path(output_path, "functional_sample_summary.csv"),
  commit_message = "Update functional sample summary",
  owner = owner, repo = repo, branch = branch
)

# 5. Create and upload summary manifest
summary_text <- glue::glue("
═══════════════════════════════════════════════════════════════════════════════
ANALYSIS REPORT: DIVERSITY INDICES CALCULATION
═══════════════════════════════════════════════════════════════════════════════

Script: 05_diversity_indices_calculation.R
Generated: {Sys.time()}
Repository: {owner}/{repo}
Branch: {branch}

PIPELINE WORKFLOW:
  Input:  {input_path}
  Output (Data): {output_path}
  Output (Figures): {figure_output_path}

═══════════════════════════════════════════════════════════════════════════════
REPOSITORY STRUCTURE
═══════════════════════════════════════════════════════════════════════════════

{repo}/
├── results/                          (Data outputs, analysis-ready files)
│   ├── 00_raw_data/                 Input: Raw sequencing data
│   ├── 01_preprocessing/            QC metrics, filtered reads
│   ├── 02_taxonomic_profiling/      Species abundance matrices
│   ├── 03_functional_profiling/     Pathway abundance matrices
│   ├── 04_diversity_analysis/       Alpha/beta diversity ← CURRENT
│   ├── 05_statistical_analysis/     Stats, models, tests
│   └── 06_reports/                  Final summaries
│
└── figures/                          (All visualizations)
    ├── 01_preprocessing/            QC plots, read quality
    ├── 02_taxonomic_profiling/      Taxonomy barplots, heatmaps
    ├── 03_functional_profiling/     Pathway plots
    ├── 04_diversity_analysis/       Alpha/beta plots ← FIGURES HERE
    ├── 05_statistical_analysis/     Statistical plots
    └── 06_reports/                  Publication figures

═══════════════════════════════════════════════════════════════════════════════
FILES GENERATED
═══════════════════════════════════════════════════════════════════════════════

PREPROCESSING OUTPUTS → {pipeline_paths$preprocessing}
  ├── sequencing_run_summary_statistics.csv
  │   └── Cleaned preprocessing metadata with QC metrics

TAXONOMIC PROFILING OUTPUTS → {pipeline_paths$taxonomic}
  ├── species_profile_unfiltered.csv
  │   └── Complete unfiltered species abundance matrix
  ├── species_profile_filtered.csv
  │   └── Filtered species abundance matrix (max abundance > 0.1, prevalence > 10%)
  └── species_prevalence_breakdown.csv
      └── Species prevalence with full taxonomic breakdown (kingdom to strain)

FUNCTIONAL PROFILING OUTPUTS → {pipeline_paths$functional}
  ├── functional_profile_filtered.csv
  │   └── Filtered pathway abundance matrix (max abundance > 0.1, prevalence > 10%)
  └── functional_prevalence_breakdown.csv
      └── Pathway prevalence across all samples

DIVERSITY ANALYSIS OUTPUTS → {output_path}
  
  Species Diversity Metrics:
  ├── species_alpha_diversity.csv
  │   └── Richness (q=0), Shannon (q=1), Simpson (q=2) indices per sample
  ├── species_beta_diversity_TD_q2.csv
  │   └── Pairwise total dissimilarity (Hill numbers, q=2)
  ├── species_beta_diversity_distance_matrix.csv
  │   └── Square distance matrix for ordination analyses
  └── species_sample_summary.csv
      └── Per-sample summary: total abundance + all diversity metrics
  
  Functional Diversity Metrics:
  ├── functional_alpha_diversity.csv
  │   └── Richness (q=0), Shannon (q=1), Simpson (q=2) indices per sample
  ├── functional_beta_diversity_TD_q2.csv
  │   └── Pairwise total dissimilarity (Hill numbers, q=2)
  ├── functional_beta_diversity_distance_matrix.csv
  │   └── Square distance matrix for ordination analyses
  └── functional_sample_summary.csv
      └── Per-sample summary: total abundance + all diversity metrics

FUTURE FIGURES (to be generated) → {figure_output_path}
  
  Suggested visualizations for this analysis step:
  ├── species_alpha_diversity_boxplot.png
  ├── species_beta_diversity_pcoa.png
  ├── species_rarefaction_curves.png
  ├── functional_alpha_diversity_boxplot.png
  ├── functional_beta_diversity_pcoa.png
  └── diversity_comparison_heatmap.png

═══════════════════════════════════════════════════════════════════════════════
ANALYSIS PARAMETERS
═══════════════════════════════════════════════════════════════════════════════

Filtering Criteria:
  • Maximum abundance threshold: > 0.1
  • Prevalence threshold: Present in > 10% of samples
  • Samples with zero abundance removed

Diversity Metrics:
  • Alpha diversity: Hill numbers
    - q = 0 (Richness): Number of features
    - q = 1 (Shannon): Exponential of Shannon entropy
    - q = 2 (Simpson): Inverse Simpson concentration
  
  • Beta diversity: Hill-based total dissimilarity
    - q = 2 (Simpson-type): Emphasizes dominant features
    - Pairwise comparisons between all samples
    - Distance matrices suitable for PCoA, NMDS, etc.

═══════════════════════════════════════════════════════════════════════════════
RESULTS SUMMARY
═══════════════════════════════════════════════════════════════════════════════

SPECIES ANALYSIS:
  Features retained: {nrow(diversity_results$species_profile$filtered_profile)}
  Samples analyzed:  {ncol(diversity_results$species_profile$filtered_profile)}
  Total abundance:   {format(sum(diversity_results$species_profile$sample_summary$total_abundance), big.mark=',')}
  
  Alpha Diversity Range:
    Richness: {min(diversity_results$species_profile$alpha_diversity$Richness)} - {max(diversity_results$species_profile$alpha_diversity$Richness)}
    Shannon:  {round(min(diversity_results$species_profile$alpha_diversity$Shannon), 2)} - {round(max(diversity_results$species_profile$alpha_diversity$Shannon), 2)}
    Simpson:  {round(min(diversity_results$species_profile$alpha_diversity$Simpson), 2)} - {round(max(diversity_results$species_profile$alpha_diversity$Simpson), 2)}

FUNCTIONAL ANALYSIS:
  Pathways retained: {nrow(diversity_results$functional_profile$filtered_profile)}
  Samples analyzed:  {ncol(diversity_results$functional_profile$filtered_profile)}
  Total abundance:   {format(sum(diversity_results$functional_profile$sample_summary$total_abundance), big.mark=',')}
  
  Alpha Diversity Range:
    Richness: {min(diversity_results$functional_profile$alpha_diversity$Richness)} - {max(diversity_results$functional_profile$alpha_diversity$Richness)}
    Shannon:  {round(min(diversity_results$functional_profile$alpha_diversity$Shannon), 2)} - {round(max(diversity_results$functional_profile$alpha_diversity$Shannon), 2)}
    Simpson:  {round(min(diversity_results$functional_profile$alpha_diversity$Simpson), 2)} - {round(max(diversity_results$functional_profile$alpha_diversity$Simpson), 2)}

═══════════════════════════════════════════════════════════════════════════════
NEXT STEPS
═══════════════════════════════════════════════════════════════════════════════

The following analyses can now be performed using these outputs:

1. Create Visualizations (→ {figure_output_path})
   • Alpha diversity boxplots/violin plots by group
   • PCoA/NMDS ordinations using beta diversity matrices
   • Rarefaction curves
   • Diversity heatmaps

2. Statistical Analysis (→ {pipeline_paths$statistics})
   • Compare diversity between groups (diet, treatment, timepoint)
   • PERMANOVA on beta diversity matrices
   • Differential abundance testing
   • Correlation analyses
   
   Figures for statistics (→ {figure_paths$statistics})

3. Final Reports (→ {pipeline_paths$reports})
   • Integrated analysis summaries
   • Publication-ready figures
   • Supplementary tables

═══════════════════════════════════════════════════════════════════════════════
")

# Save summary to temporary file and upload
temp_summary <- tempfile(fileext = ".txt")
cat(summary_text, file = temp_summary)

save_and_upload_to_github(
  object_type = "file",
  local_path = temp_summary,
  github_path = file.path("results/Processed/06_reports", "ANALYSIS_SUMMARY_r_script_1.Diversity_index_formation_notthingham.txt"),
  commit_message = "Update analysis summary",
  owner = owner, repo = repo, branch = branch
)

unlink(temp_summary)

# ────────────────────────────────────────────────────────────────────────────────
# Section: Final Summary
# ────────────────────────────────────────────────────────────────────────────────

message("✅ ALL ANALYSES COMPLETE AND UPLOADED TO GITHUB")

message(glue::glue("📊 Species Profile: {nrow(diversity_results$species_profile$filtered_profile)} features, {ncol(diversity_results$species_profile$filtered_profile)} samples"))
message(glue::glue("📊 Functional Profile: {nrow(diversity_results$functional_profile$filtered_profile)} pathways, {ncol(diversity_results$functional_profile$filtered_profile)} samples"))
message(glue::glue("\n📁 Complete Repository Structure Created:"))
message(glue::glue("\n   Results Pipeline:"))
message(glue::glue("   ├── {pipeline_paths$raw_data} (raw input data)"))
message(glue::glue("   ├── {pipeline_paths$preprocessing} (QC & preprocessing)"))
message(glue::glue("   ├── {pipeline_paths$taxonomic} (species profiles)"))
message(glue::glue("   ├── {pipeline_paths$functional} (pathway profiles)"))
message(glue::glue("   ├── {output_path} (diversity metrics) ← CURRENT OUTPUT"))
message(glue::glue("   ├── {pipeline_paths$statistics} (stats - ready for next step)"))
message(glue::glue("   └── {pipeline_paths$reports} (final reports)"))
message(glue::glue("\n   Figures Pipeline:"))
message(glue::glue("   ├── {figure_paths$preprocessing}"))
message(glue::glue("   ├── {figure_paths$taxonomic}"))
message(glue::glue("   ├── {figure_paths$functional}"))
message(glue::glue("   ├── {figure_output_path} ← FIGURES GO HERE"))
message(glue::glue("   ├── {figure_paths$statistics}"))
message(glue::glue("   └── {figure_paths$reports}"))
message(glue::glue("\n🔗 View results at: https://github.com/{owner}/{repo}/tree/{branch}/{output_path}"))
message(glue::glue("🔗 View figures at: https://github.com/{owner}/{repo}/tree/{branch}/{figure_output_path}"))
message(glue::glue("🔗 View repository: https://github.com/{owner}/{repo}"))
