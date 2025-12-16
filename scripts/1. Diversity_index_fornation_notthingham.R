# ────────────────────────────────────────────────────────────────────────────────
# Section: Reproducibility
# ────────────────────────────────────────────────────────────────────────────────
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
directory_path <- "results/Raw"
output_directory_path <- "results/processed"  # Where to upload processed results
token <- "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"  # Replace or use a .Renviron variable
Sys.setenv(GITHUB_PAT = "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")
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
  
  # Encode file content to Base64 for GitHub API
  encoded_content <- base64enc::base64encode(local_path)
  
  # Check if the file already exists on GitHub
  file_info <- tryCatch(
    gh::gh("GET /repos/:owner/:repo/contents/:path",
           owner = owner, repo = repo, path = github_path, ref = branch),
    error = function(e) NULL
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
         sha = sha)
  
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
message(glue::glue("📂 Fetching file list from: {owner}/{repo}/{directory_path}"))

# Use GitHub API to list directory contents
api_url <- glue::glue("https://api.github.com/repos/{owner}/{repo}/contents/{directory_path}")
response <- GET(api_url, 
                authenticate(token, ""),
                add_headers("Accept" = "application/vnd.github.v3+json"))

if (status_code(response) != 200) {
  stop(glue::glue("❌ Failed to access directory. Status code: {status_code(response)}
                   Check: 
                   1. Repository exists and is accessible
                   2. Directory path is correct: {directory_path}
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
message("📤 Uploading results to GitHub...")
message(strrep("─", 80))

# 1. Upload preprocessing summary
save_and_upload_to_github(
  object = data[["preprocessing"]],
  object_type = "data",
  github_path = file.path(output_directory_path, "sequencing_run_summary_statistics.csv"),
  commit_message = "Update sequencing run summary statistics",
  owner = owner, repo = repo, branch = branch
)

# 2. Upload unfiltered species profile
save_and_upload_to_github(
  object = data$metaphlan,
  object_type = "data",
  github_path = file.path(output_directory_path, "species_profile_unfiltered.csv"),
  commit_message = "Update unfiltered species profile",
  owner = owner, repo = repo, branch = branch
)

# 3. Upload species-specific results
save_and_upload_to_github(
  object = diversity_results$species_profile$filtered_profile %>% rownames_to_column("Taxa"),
  object_type = "data",
  github_path = file.path(output_directory_path, "species_profile_filtered.csv"),
  commit_message = "Update filtered species profile",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$species_profile$prevalence,
  object_type = "data",
  github_path = file.path(output_directory_path, "species_prevalence_breakdown.csv"),
  commit_message = "Update species prevalence breakdown",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$species_profile$alpha_diversity,
  object_type = "data",
  github_path = file.path(output_directory_path, "species_alpha_diversity.csv"),
  commit_message = "Update species alpha diversity",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$species_profile$beta_diversity,
  object_type = "data",
  github_path = file.path(output_directory_path, "species_beta_diversity_TD_q2.csv"),
  commit_message = "Update species beta diversity",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$species_profile$beta_distance_matrix %>% rownames_to_column("sample_id"),
  object_type = "data",
  github_path = file.path(output_directory_path, "species_beta_diversity_distance_matrix.csv"),
  commit_message = "Update species beta diversity distance matrix",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$species_profile$sample_summary,
  object_type = "data",
  github_path = file.path(output_directory_path, "species_sample_summary.csv"),
  commit_message = "Update species sample summary",
  owner = owner, repo = repo, branch = branch
)

# 4. Upload functional profile-specific results
save_and_upload_to_github(
  object = diversity_results$functional_profile$filtered_profile %>% rownames_to_column("Pathway"),
  object_type = "data",
  github_path = file.path(output_directory_path, "functional_profile_filtered.csv"),
  commit_message = "Update filtered functional profile",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$functional_profile$prevalence,
  object_type = "data",
  github_path = file.path(output_directory_path, "functional_prevalence_breakdown.csv"),
  commit_message = "Update functional prevalence breakdown",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$functional_profile$alpha_diversity,
  object_type = "data",
  github_path = file.path(output_directory_path, "functional_alpha_diversity.csv"),
  commit_message = "Update functional alpha diversity",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$functional_profile$beta_diversity,
  object_type = "data",
  github_path = file.path(output_directory_path, "functional_beta_diversity_TD_q2.csv"),
  commit_message = "Update functional beta diversity",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$functional_profile$beta_distance_matrix %>% rownames_to_column("sample_id"),
  object_type = "data",
  github_path = file.path(output_directory_path, "functional_beta_diversity_distance_matrix.csv"),
  commit_message = "Update functional beta diversity distance matrix",
  owner = owner, repo = repo, branch = branch
)

save_and_upload_to_github(
  object = diversity_results$functional_profile$sample_summary,
  object_type = "data",
  github_path = file.path(output_directory_path, "functional_sample_summary.csv"),
  commit_message = "Update functional sample summary",
  owner = owner, repo = repo, branch = branch
)

# 5. Create and upload summary manifest
summary_text <- glue::glue("
Analysis Report Generated: {Sys.time()}
Repository: {owner}/{repo}
Branch: {branch}

═══════════════════════════════════════════════════════════════════════════════

FILES GENERATED:

PREPROCESSING:
  - sequencing_run_summary_statistics.csv
    Description: Cleaned preprocessing metadata

SPECIES PROFILE:
  - species_profile_unfiltered.csv
    Description: Complete unfiltered species abundance matrix
  - species_profile_filtered.csv
    Description: Filtered species abundance matrix (max abundance > 0.1, prevalence > 10%)
  - species_prevalence_breakdown.csv
    Description: Species prevalence with taxonomic breakdown
  - species_alpha_diversity.csv
    Description: Alpha diversity metrics (Richness, Shannon, Simpson)
  - species_beta_diversity_TD_q2.csv
    Description: Pairwise beta diversity (Total Dissimilarity, q=2)
  - species_beta_diversity_distance_matrix.csv
    Description: Beta diversity distance matrix
  - species_sample_summary.csv
    Description: Per-sample summary statistics

FUNCTIONAL PROFILE:
  - functional_profile_filtered.csv
    Description: Filtered pathway abundance matrix (max abundance > 0.1, prevalence > 10%)
  - functional_prevalence_breakdown.csv
    Description: Pathway prevalence across samples
  - functional_alpha_diversity.csv
    Description: Alpha diversity metrics (Richness, Shannon, Simpson)
  - functional_beta_diversity_TD_q2.csv
    Description: Pairwise beta diversity (Total Dissimilarity, q=2)
  - functional_beta_diversity_distance_matrix.csv
    Description: Beta diversity distance matrix
  - functional_sample_summary.csv
    Description: Per-sample summary statistics

═══════════════════════════════════════════════════════════════════════════════

ANALYSIS PARAMETERS:
  - Filtering threshold: Maximum abundance > 0.1
  - Prevalence threshold: Present in > 10% of samples
  - Beta diversity: Hill numbers with q = 2 (Simpson-type)
  - Alpha diversity: Hill numbers q = 0 (Richness), q = 1 (Shannon), q = 2 (Simpson)

SPECIES ANALYSIS SUMMARY:
  - Features retained: {nrow(diversity_results$species_profile$filtered_profile)}
  - Samples analyzed: {ncol(diversity_results$species_profile$filtered_profile)}
  - Total abundance: {sum(diversity_results$species_profile$sample_summary$total_abundance)}

FUNCTIONAL ANALYSIS SUMMARY:
  - Pathways retained: {nrow(diversity_results$functional_profile$filtered_profile)}
  - Samples analyzed: {ncol(diversity_results$functional_profile$filtered_profile)}
  - Total abundance: {sum(diversity_results$functional_profile$sample_summary$total_abundance)}
")

# Save summary to temporary file and upload
temp_summary <- tempfile(fileext = ".txt")
cat(summary_text, file = temp_summary)

save_and_upload_to_github(
  object_type = "file",
  local_path = temp_summary,
  github_path = file.path(output_directory_path, "ANALYSIS_SUMMARY.txt"),
  commit_message = "Update analysis summary",
  owner = owner, repo = repo, branch = branch
)

unlink(temp_summary)

# ────────────────────────────────────────────────────────────────────────────────
# Section: Final Summary
# ────────────────────────────────────────────────────────────────────────────────
message("\n" , strrep("═", 80))
message("✅ ALL ANALYSES COMPLETE AND UPLOADED TO GITHUB")
message(strrep("═", 80))
message(glue::glue("📊 Species Profile: {nrow(diversity_results$species_profile$filtered_profile)} features, {ncol(diversity_results$species_profile$filtered_profile)} samples"))
message(glue::glue("📊 Functional Profile: {nrow(diversity_results$functional_profile$filtered_profile)} pathways, {ncol(diversity_results$functional_profile$filtered_profile)} samples"))
message(glue::glue("🔗 View results at: https://github.com/{owner}/{repo}/tree/{branch}/{output_directory_path}"))
message(strrep("═", 80))