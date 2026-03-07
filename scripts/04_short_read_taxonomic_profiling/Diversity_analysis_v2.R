



# ────────────────────────────────────────────────────────────────────────────────
# Pipeline: Microbiome–Metabolome Association Tests
# Purpose:
#   - Process MetaPhlAn species profiles
#   - Align microbiome & metabolomics datasets
#   - Apply CLR transformation (compositionality correction)
#   - Compute three association tests: Mantel, Procrustes, RV coefficient
#   - Save results in a single structured object
# Author: [Your Name]
# Date: [YYYY-MM-DD]
# ────────────────────────────────────────────────────────────────────────────────

# ────────────────────────────────────────────────────────────────────────────────
# Section: Set Custom R Library Path
# ────────────────────────────────────────────────────────────────────────────────

# Candidate paths (modify as needed for your system)

.libPaths("G:/R/win-library/4.4")

# ────────────────────────────────────────────────────────────────────────────────
# Section: Load Required Packages Safely
# ────────────────────────────────────────────────────────────────────────────────

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

#BiocManager::install("mixOmics")
#library(mixOmics) 
# Load packages using pacman
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
suppressPackageStartupMessages(
  pacman::p_load(
    readxl, readr, reshape2, dplyr, gplots, Heatplus, vegan, RColorBrewer,
    tidyr, gtools, stringr, tidyverse, ComplexHeatmap, magick, viridis, 
    Hotelling, remotes, ggnewscale, ggtreeExtra, jmvcore, gtsummary, 
    officer, fpc, data.table, ggpubr, grid, cluster,dplyr, tibble,ggstatsplot
  )
)

# ────────────────────────────────────────────────────────────────────────────────
# Section: Reproducibility
# ────────────────────────────────────────────────────────────────────────────────

set.seed(123)

# ────────────────────────────────────────────────────────────────────────────────
# Section: Load Metadata Files
# ────────────────────────────────────────────────────────────────────────────────

metadata_dir <- "Q:/CoDiet/liam walsh/NottinghamJuly2025/officer_output/Reports"



if (!dir.exists(metadata_dir)) stop("❌ Metadata directory not found.")

# Load all files except .docx
file_list <- list.files(metadata_dir, full.names = TRUE)
valid_files <- file_list[!grepl("\\.docx$", file_list, ignore.case = TRUE)]

data <- list()

for (file_path in valid_files) {
  tryCatch({
    if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
      data[[file_path]] <- read_csv(file_path, show_col_types = FALSE)
    } else if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
      data[[file_path]] <- read_delim(file_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
    } else {
      warning(glue::glue("⚠️ Skipping unrecognised file format: {file_path}"))
    }
  }, error = function(e) {
    warning(glue::glue("❌ Error reading file {file_path}: {e$message}"))
  })
}

# Clean file names
names(data) <- gsub(paste0(metadata_dir, "/|\\.csv$|\\.tsv$"), "", names(data))



names(data )[grep("metaphlan", names(data ))]="metaphlan"

# ------------------------------------------------------------------------------
# Create turmeric-status lookup table
# ID = participant number (from sample ID, e.g. "CD_003" → 3)
# Status = duration of turmeric supplementation or "Baseline only"
# ------------------------------------------------------------------------------

turmeric_status <- data.frame(
  ID  = c(3,4,5,6,7,9,11,14,15,16,18,19,20,21,22,23),
  status = c(
    "Didn't take turmeric",
    "Turmeric 6 months",
    "Turmeric 6 months",
    "Didn't take turmeric",
    "Turmeric 6 months",
    "Turmeric 6 months",
    "Turmeric 3 months",
    "Baseline sample only",
    "Turmeric 6 months",
    "Turmeric 6 months",
    "Didn't take turmeric",
    "Baseline sample only",
    "Turmeric 3 months",
    "Turmeric 6 months",
    "Turmeric 6 months",
    "Baseline sample only"
  ),
  stringsAsFactors = FALSE
)

# Quick check of the distribution of statuses
table(turmeric_status$status)


# ------------------------------------------------------------------------------
# Extract participant ID and centre code from sample IDs
# (ids look like "Nottingham_11A_S19").
# We keep:
#   block_code   → e.g. "11A"
#   block_number → numeric part e.g. 11
#   block_letter → letter part  e.g. "A"
# ------------------------------------------------------------------------------

ids <- unique(data$alpha_diversity$sample_id)
split_ids <- str_split_fixed(ids, "_", 3)

sample_ids <- data.frame(
  full_id      = ids,
  block_code   = split_ids[,2],
  block_number = as.integer(gsub("[^0-9]", "", split_ids[,2])),
  block_letter = gsub("[0-9]", "", split_ids[,2]),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# Add turmeric-status to sample table (join on block_number = ID)
# ------------------------------------------------------------------------------

sample_ids <- merge(sample_ids,
                    turmeric_status,
                    by.x = "block_number",
                    by.y = "ID",
                    all.x = TRUE)
sample_ids$timepoint="Baseline"

sample_ids$timepoint[which(sample_ids$block_letter%in% c("a","A"))]="Baseline"

sample_ids$timepoint[which(sample_ids$block_letter%in% c("b","B"))]="2 weeks"


sample_ids$timepoint[which(sample_ids$block_letter%in% c("C","C"))]="6 months"


# ------------------------------------------------------------------------------
# Create indicator flag for turmeric use
#   TRUE  → participant was given turmeric
#   FALSE → no turmeric
# ------------------------------------------------------------------------------
# 
 sample_ids$turmeric_group <- "FALSE"

sample_ids$turmeric_group[
  grepl("Turmeric", sample_ids$status) & sample_ids$block_letter %in% c("c", "C")
]  <- "TRUE"
# 
# # Optionally: mark "stopped" for 3-month turmeric in visit C
 sample_ids$turmeric_group[
   which(sample_ids$status == "Turmeric 3 months" & sample_ids$block_letter == "C")
 ] <- "Stopped"
# 

 sample_ids$block_letter[which( sample_ids$turmeric_group=="TRUE")]


data$alpha_diversity=
  merge(sample_ids, 
        data$alpha_diversity,
        by.x = "full_id",
        by.y="sample_id",
        all.x="TRUE")

data$alpha_diversity$block_letter=
  str_to_upper(data$alpha_diversity$block_letter)






# Load required package
library(ggstatsplot)
library(ggplot2)

# Define output directory
outdir <- "Q:/CoDiet/liam walsh/NottinghamJuly2025/plots"

# Create directory if it doesn’t exist
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Custom plotting function for reusability
make_plot <- function(data, x_var, y_var, filename) {
  p <- ggbetweenstats(
    data = data,
    x = !!sym(x_var),
    y = !!sym(y_var),
    fill = turmeric_group,
    color = turmeric_group,
    type = "nonparametric",
    plot.type = "box",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    centrality.plotting = FALSE,
    ggsignif.args = list(textsize = 1, tip_length = 0.01),
    bf.message = FALSE,
    xlab = "Group",
    ylab = "Alpha diversity values",
    ggtheme = ggplot2::theme_bw(),
    ggplot.component = list(
      theme(
        text = element_text(size = 7),
        panel.border = element_blank()
      )
    )
  )
  
  # Save plot
  ggsave(
    filename = file.path(outdir, filename),
    plot = p,
    width = 6,
    height = 4,
    units = "in",
    dpi = 300
  )
}

#------------------------------
# 1. Block-based comparisons
#------------------------------
make_plot(data$alpha_diversity, "block_letter", "Shannon", "alpha_Shannon_block.png")
make_plot(data$alpha_diversity, "block_letter", "Richness", "alpha_Richness_block.png")
make_plot(data$alpha_diversity, "block_letter", "Simpson", "alpha_Simpson_block.png")

#------------------------------
# 2. Turmeric group comparisons (excluding "Stopped")
#------------------------------



ggbetweenstats(
  data =   
    filtered_data ,
  x = block_letter,
  y = Shannon,
  fill = turmeric_group,
  color = turmeric_group,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  ggsignif.args = list(textsize = 1, tip_length = 0.01),
  bf.message = FALSE,
  xlab = "Group",
  ylab = "Alpha diversity values",
  ggtheme = ggplot2::theme_bw(),
  ggplot.component = list(
    theme(
      text = element_text(size = 7),
      panel.border = element_blank()
    )
  )
)


  
filtered_data <- data$alpha_diversity[
  which(data$alpha_diversity$turmeric_group!="Stopped"),]


filtered_data <- data$alpha_diversity[
  which(data$alpha_diversity$status=="Turmeric 6 months"),]



make_plot(filtered_data, "turmeric_group", "Shannon", "alpha_Shannon_turmeric.png")
make_plot(filtered_data, "turmeric_group", "Richness", "alpha_Richness_turmeric.png")
make_plot(filtered_data, "turmeric_group", "Simpson", "alpha_Simpson_turmeric.png")


library(ggplot2)
library(dplyr)

data$alpha_diversity %>%
  ggplot(aes(x = block_letter, y = Richness, fill = turmeric_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, colour = "black") +
  geom_jitter(aes(colour = turmeric_group),
              width = 0.15, size = 4.5, alpha = 0.6, show.legend = FALSE) +
  scale_fill_brewer(palette = "Set2") +
  scale_colour_brewer(palette = "Set2") +
  labs(
    title = "Richness across blocks by turmeric group",
    x = "Block letter",
    y = "Richness (observed species)",
    fill = "Turmeric group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 0.6),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)
  )



library(dplyr)
library(purrr)

# Filter out stopped group
df_filtered <- data$alpha_diversity %>% filter(turmeric_group != "Stopped")

# Split by block_letter
groups <- df_filtered %>% group_by(block_letter) %>% group_split()

# Run Kruskal-Wallis test per block
kruskal_results <- map(groups, ~ {
  tibble(
    block_letter = unique(.x$block_letter),
    kruskal_p = kruskal.test(Richness ~ turmeric_group, data = .x)$p.value
  )
})

# Combine results
kruskal_results <- bind_rows(kruskal_results)
kruskal_results



# Confirmation message
message("✅ All ggstatsplot figures saved to: ", outdir)

beta_summary=dplyr::select( data$beta_diversity_TD_q2_matrix, site1,site2,TD_beta)


dist <- rbind(
  beta_summary %>% dplyr::select(site1, site2, TD_beta),
  beta_summary %>%
    dplyr::select(site2, site1, TD_beta) %>%
    dplyr::rename(site1 = 1, site2 = 2)
) %>%
  tidyr::pivot_wider(names_from = site2, values_from = TD_beta) %>%
  tibble::column_to_rownames("site1")

# Fill NA entries with 0 to ensure complete square distance matrix
dist[is.na(dist)] <- 0


# ────────────────────────────────────────────────────────────────────────────────
# Perform PCoA (Principal Coordinates Analysis)
# ────────────────────────────────────────────────────────────────────────────────
pcoa <- cmdscale(dist, k = 2, eig = TRUE, add = TRUE)
percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)  # % variance explained

# Convert PCoA results to dataframe
pcoa_df <- data.frame(PC1 = pcoa$points[, 1], PC2 = pcoa$points[, 2])

# ────────────────────────────────────────────────────────────────────────────────
# Merge PCoA Output with Sample Metadata
# ────────────────────────────────────────────────────────────────────────────────


library(ggplot2)

# Compute % variance explained for PC1 and PC2
percent_explained <- round(100 * pcoa$eig / sum(pcoa$eig), 1)



pcoa_df <- merge(
  pcoa_df,
  sample_ids,
  by.x = 0,
  by.y = "full_id",
  all.x = TRUE
)

# Ensure turmeric_group is a factor
pcoa_df$turmeric_group <- factor(pcoa_df$turmeric_group, levels = c("TRUE", "FALSE", "Stopped"))



# Plot
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = turmeric_group)) +
  geom_point(size = 5, alpha = 0.8) +                                      # Points
  stat_ellipse(aes(fill = turmeric_group), type = "t", alpha = 0.2, 
               geom = "polygon", color = NA) +                             # Ellipses
  labs(
    x = paste0("PC1 (", percent_explained[1], "%)"),
    y = paste0("PC2 (", percent_explained[2], "%)"),
    color = "Turmeric Group",
    fill = "Turmeric Group"
  ) +
  scale_color_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02", "Stopped" = "#7570b3")) +
  scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02", "Stopped" = "#7570b3")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold")
  )





ggplot(pcoa_df, aes(PC1, PC2)) +
  geom_point(size = 5, aes(colour = as.factor(`turmeric_group`)))

pcoa_df$block_letter=
  str_to_upper(pcoa_df$block_letter)

ggplot(pcoa_df, aes(PC1, PC2)) +
  geom_point(size = 5, aes(colour = as.factor(block_letter)))




#────────────────────────────────────────────────────────────────────────────────
# Calculate Centroids for Each Inflammatory Group
# ────────────────────────────────────────────────────────────────────────────────
centroid <- data.frame(
  data_source = as.character(levels(as.factor(pcoa_df$`turmeric_group`))),
  PC1 = 0, PC2 = 0
)

for (i in centroid$data_source) {
  centroid$PC1[centroid$data_source == i] <- mean(pcoa_df$PC1[pcoa_df$`turmeric_group` == i])
  centroid$PC2[centroid$data_source == i] <- mean(pcoa_df$PC2[pcoa_df$`turmeric_group` == i])
}

# Reset row names for compatibility
pcoa_df <- pcoa_df %>%
  remove_rownames() %>%
  column_to_rownames("Row.names")

pcoa_df <- pcoa_df[
  which(pcoa_df$turmeric_group!="Stopped"),]

# ────────────────────────────────────────────────────────────────────────────────
# Step 1: Subset distance matrix and metadata to common sample IDs
# ────────────────────────────────────────────────────────────────────────────────

# Identify common samples between the distance matrix and metadata
common_samples <- intersect(rownames(dist), rownames(pcoa_df))

# Stop if too few samples are shared
if (length(common_samples) < 2) {
  stop("Too few common samples between distance matrix and metadata to perform Adonis.")
}

# Subset both the distance matrix and metadata to the shared sample set
dist_subset <- dist[common_samples, common_samples]
pcoa_df_subset <- pcoa_df[common_samples, , drop = FALSE]

# Check that the sample order matches between the distance matrix and metadata
if (!identical(rownames(dist_subset), rownames(pcoa_df_subset))) {
  stop("Sample names do not align after subsetting — check your inputs.")
}

# ────────────────────────────────────────────────────────────────────────────────
# Step 2: Run PERMANOVA (Adonis2) — Single Factor
# ────────────────────────────────────────────────────────────────────────────────

# Test for overall effect of inflammatory status
adonis_result <- adonis2(
  dist_subset ~ turmeric_group,
  data = pcoa_df_subset,
  permutations = 10000,
  na.rm = TRUE
)


# Test for overall effect of inflammatory status
adonis_result <- adonis2(
  dist_subset ~ block_letter ,
  data = pcoa_df_subset,
  permutations = 10000,
  na.rm = TRUE
)






# ──────────────────────────────────────────────────
# 6. Metadata Processing
# ──────────────────────────────────────────────────

library(dplyr)
library(stringr)
library(tibble)

# Clean participant ID and block
sample_ids <- sample_ids %>%
  mutate(participant = gsub("[0-9]+","", block_code),
         block_letter = str_to_upper(block_letter))

# Filter out stopped participants and select relevant columns
s_meta <- sample_ids %>%
  filter(turmeric_group != "Stopped") %>%
  select(full_id, turmeric_group, block_letter, participant) %>%
  mutate(sampleNames_row = full_id) %>%
  column_to_rownames("sampleNames_row") %>%
  select(turmeric_group, block_letter, participant)

# ──────────────────────────────────────────────────
# 7. OTU Table Setup
# ──────────────────────────────────────────────────

data$species_profile_filtered <- data$species_profile_filtered %>%
  select(c("Taxa", intersect(colnames(.), rownames(s_meta))))

s_otu_tab <- data$species_profile_filtered %>%
  column_to_rownames("Taxa")

# Remove NAs and zero-sum taxa/samples
s_otu_tab[is.na(s_otu_tab)] <- 0
s_otu_tab <- s_otu_tab[rowSums(s_otu_tab) > 0, , drop = FALSE]
s_otu_tab <- s_otu_tab[, colSums(s_otu_tab) > 0, drop = FALSE]

# Align metadata with OTU table
s_meta <- s_meta[rownames(s_meta) %in% colnames(s_otu_tab), , drop = FALSE]

# Filter low-variance and low-mean taxa
s_otu_tab <- s_otu_tab[apply(s_otu_tab, 1, sd) > 0, , drop = FALSE]
s_otu_tab <- s_otu_tab[rowMeans(s_otu_tab) > 0, , drop = FALSE]

# Taxonomy placeholder
taxonomy_df <- data.frame(Species = rownames(s_otu_tab), row.names = rownames(s_otu_tab))

# ──────────────────────────────────────────────────
# 8. Build Phyloseq Object
# ──────────────────────────────────────────────────

library(phyloseq)

ps_codiet <- phyloseq(
  sample_data(s_meta),
  otu_table(as.matrix(s_otu_tab), taxa_are_rows = TRUE),
  tax_table(as.matrix(taxonomy_df))
)

sample_data(ps_codiet)$turmeric_group <- factor(sample_data(ps_codiet)$turmeric_group)
sample_data(ps_codiet)$block_letter   <- factor(sample_data(ps_codiet)$block_letter)
sample_data(ps_codiet)$participant    <- factor(sample_data(ps_codiet)$participant)

# ──────────────────────────────────────────────────
# 9. Run ANCOM-BC2 with repeated measures and block effects
# ──────────────────────────────────────────────────

library(ANCOMBC)
library(microbiome)
ancombc2_out <- ancombc2(
  data         = ps_codiet,
  assay.type   = "counts",
  rank         = NULL,
  fix_formula  = "turmeric_group + block_letter",   # fixed effects only
  # random effect for repeated measures
  p_adj_method = "fdr",
  prv_cut      = 0.70,
  lib_cut      = 0,
  group        = "turmeric_group",
  struc_zero   = TRUE,
  neg_lb       = TRUE,
  alpha        = 0.05,
  global       = TRUE,
  verbose      = TRUE
)


res2 <- ancombc2_out$res


#
ancom_summary <- res2 %>%
  select(taxon, lfc_turmeric_groupTRUE, se_turmeric_groupTRUE,
         W_turmeric_groupTRUE, p_turmeric_groupTRUE, q_turmeric_groupTRUE,
         passed_ss_turmeric_groupTRUE) %>%
  rename(
    Taxon = taxon,
    Log2FC = lfc_turmeric_groupTRUE,
    SE = se_turmeric_groupTRUE,
    W_stat = W_turmeric_groupTRUE,
    p_value = p_turmeric_groupTRUE,
    q_value = q_turmeric_groupTRUE,
    Significant = passed_ss_turmeric_groupTRUE
  )

# View table
ancom_summary




# ──────────────────────────────────────────────────
# 6. Metadata Processing
# ──────────────────────────────────────────────────


library(dplyr)

sample_ids <- sample_ids %>%
  mutate(participant = gsub("[0-9]+","", block_code) )%>% 
  mutate(block_letter=str_to_upper(block_letter))





s_meta <- sample_ids %>%
  filter(turmeric_group !="Stopped") %>% 
  dplyr::select(full_id, turmeric_group, block_letter ,participant ) %>%
  #filter(!is.na(HealthStatus)) %>%
  mutate(sampleNames_row =full_id) %>%
  column_to_rownames("sampleNames_row") %>%
  dplyr::select(turmeric_group)



# ──────────────────────────────────────────────────
# 7. OTU Table Setup
# ──────────────────────────────────────────────────


data$species_profile_filtered<- data$species_profile_filtered %>%
  #filter(species != "UNCLASSIFIED") %>%
  dplyr::select(c("Taxa", intersect(colnames(.), rownames(s_meta))))

s_otu_tab <- data$species_profile_filtered%>%
  column_to_rownames("Taxa")

s_otu_tab[is.na(s_otu_tab)] <- 0
s_otu_tab <- s_otu_tab[rowSums(s_otu_tab) > 0, , drop = FALSE]
s_otu_tab <- s_otu_tab[, colSums(s_otu_tab) > 0, drop = FALSE]

nonzero_samples <- colSums(s_otu_tab) > 0
s_otu_tab <- s_otu_tab[, nonzero_samples, drop = FALSE]
s_meta <- s_meta[rownames(s_meta) %in% colnames(s_otu_tab), , drop = FALSE]

s_otu_tab <- s_otu_tab[apply(s_otu_tab, 1, sd) > 0, , drop = FALSE]
s_otu_tab <- s_otu_tab[rowMeans(s_otu_tab) > 0, , drop = FALSE]

taxonomy_df <- data.frame(Species = rownames(s_otu_tab), row.names = rownames(s_otu_tab))



# ──────────────────────────────────────────────────
# 8. Build Phyloseq Object
# ──────────────────────────────────────────────────

library( phyloseq)
ps_codiet <- phyloseq(
  sample_data(s_meta),
  otu_table(as.matrix(s_otu_tab), taxa_are_rows = TRUE),
  tax_table(as.matrix(taxonomy_df))
)

sample_data(ps_codiet)$turmeric_group <- factor(sample_data(ps_codiet)$turmeric_group)



library(ANCOMBC)
library(microbiome)
ancombc2_out <- ancombc2(
  data         = ps_codiet,
  assay.type   = "counts",
  rank         = NULL,
  fix_formula  = "turmeric_group",
  p_adj_method = "fdr",
  prv_cut      = 0.70,  # relax filter to avoid 10-taxa warning
  lib_cut      = 0,
  group        = "turmeric_group",
  struc_zero   = TRUE,
  neg_lb       = TRUE,
  alpha        = 0.05,
  global       = TRUE,
  verbose      = TRUE  # place here, not inside control lists
)



res2 <- ancombc2_out$res










