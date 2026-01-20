###############################################################################
# Step1: T Cell Clonotype Expansion Analysis
# Compares pre vs post treatment clonotype frequencies using Poisson tests
# Analyzes both TEC (T cell Engaging) and CAR-T (Cilta-cel) therapies
###############################################################################

# -----------------------------------------------------------------------------
# 1. Load Dependencies
# -----------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(qs)
library(ggplot2)
library(readxl)
library(patchwork)
library(data.table)
library(ggpubr)
library(pheatmap)
library(vegan)        # For PERMANOVA

# -----------------------------------------------------------------------------
# 2. Configuration
# -----------------------------------------------------------------------------
output_dir <- "~/TCR_project/"
plot_dir <- paste0(output_dir, "plots/")

# Statistical thresholds
FDR_THRESHOLD <- 0.1
LARGE_CLONE_THRESHOLD <- 0.01

# Plot aesthetics
CLONE_COLORS <- c(

"large clone contraction" = "darkblue",
"small clone contraction" = "#00A6D6",
"stable" = "grey60",
"small clone expansion" = "#7E57C2",
"large clone expansion" = "#D81B60"
)

AXIS_BREAKS <- c(1e-4, 1e-3, 1e-2, 1e-1)
AXIS_LABELS <- c("0.01%", "0.1%", "1%", "10%")

# -----------------------------------------------------------------------------
# 3. Helper Functions
# -----------------------------------------------------------------------------

#' Perform Poisson test comparing two counts with different sample sizes
#' @param count1 Count in condition 1
#' @param count2 Count in condition 2
#' @param size1 Total sample size for condition 1
#' @param size2 Total sample size for condition 2
#' @return p-value from Poisson test
poisson_test_single <- function(count1, count2, size1, size2) {
  result <- poisson.test(c(count1, count2), T = c(size1, size2))
  return(result$p.value)
}

#' Perform Poisson tests for all clonotypes between two conditions
#' @param counts1 Vector of counts for condition 1
#' @param counts2 Vector of counts for condition 2
#' @return Tibble with counts, proportions, p-values, FDR, and significance
poisson_test_clonotypes <- function(counts1, counts2) {
  size1 <- sum(counts1)
  size2 <- sum(counts2)

  p_values <- mapply(
    poisson_test_single, counts1, counts2,
    MoreArgs = list(size1 = size1, size2 = size2)
  )

  tibble(
    n1 = counts1,
    n2 = counts2,
    p1 = counts1 / size1,
    p2 = counts2 / size2,
    p = p_values
  ) %>%
    mutate(
      fdr = p.adjust(p, method = "BH"),
      sig = fdr < 0.05
    )
}

#' Prepare Seurat object for clonotype analysis
#' @param seurat_obj Seurat object with TCR data
#' @return Filtered Seurat object with combined clonotype column
prepare_clonotype_data <- function(seurat_obj) {
  # Remove cells without clonotype information
 remove_cells <- rownames(seurat_obj@meta.data[is.na(seurat_obj@meta.data$clonotypeID), ])
  keep_cells <- setdiff(rownames(seurat_obj@meta.data), remove_cells)

  obj_filtered <- subset(seurat_obj, cells = keep_cells)

 # Create combined clonotype identifier (TRA + TRB CDR3 sequences)
  obj_filtered@meta.data <- obj_filtered@meta.data %>%
    mutate(
      combined_clonotype = paste(TRA_cdr3_aa, TRB_cdr3_aa, sep = "_"),
      new_clones = paste0(PT_ID, "_", combined_clonotype)
    )

  return(obj_filtered)
}

#' Run Poisson analysis for all patients
#' @param meta Metadata with PT_ID, combined_clonotype, Timepoint columns
#' @param pre_label Label for pre-treatment timepoint
#' @param post_label Label for post-treatment timepoint
#' @return Data frame with Poisson test results for all patients
run_patient_poisson_analysis <- function(meta, pre_label = "Pre", post_label = "Post") {
  # Summarize clonotype counts per patient and timepoint
  clono_counts <- meta %>%
    group_by(PT_ID, Timepoint, combined_clonotype) %>%
    summarise(n = n(), .groups = "drop")

  # Keep only patients with paired samples
  paired_patients <- clono_counts %>%
    group_by(PT_ID) %>%
    summarise(n_timepoints = n_distinct(Timepoint), .groups = "drop") %>%
    filter(n_timepoints > 1) %>%
    pull(PT_ID)

  clono_counts <- clono_counts %>% filter(PT_ID %in% paired_patients)

  # Run Poisson test for each patient
  results_list <- list()

  for (pt in unique(clono_counts$PT_ID)) {
    pt_data <- clono_counts %>%
      filter(PT_ID == pt) %>%
      pivot_wider(
        names_from = Timepoint,
        values_from = n,
        values_fill = 0
      ) %>%
      arrange(combined_clonotype)

    # Extract pre and post counts
    pre_counts <- pt_data[[pre_label]]
    post_counts <- pt_data[[post_label]]

    # Run test and add clonotype info
    res <- poisson_test_clonotypes(pre_counts, post_counts) %>%
      mutate(combined_clonotype = pt_data$combined_clonotype) %>%
      arrange(fdr)

    results_list[[pt]] <- res
  }

  # Combine all patient results
  bind_rows(results_list, .id = "PT_ID")
}

#' Prepare data frame for plotting clonotype changes
#' @param res_df Results data frame from Poisson analysis
#' @param fdr_thr FDR threshold for significance
#' @param large_thr Frequency threshold for large clones
#' @return Data frame ready for plotting
prepare_plot_data <- function(res_df, fdr_thr = FDR_THRESHOLD, large_thr = LARGE_CLONE_THRESHOLD) {
  res_df %>%
    mutate(
      freq_pre = as.numeric(freq_pre),
      freq_post = as.numeric(freq_post),

      # Add pseudocount for log scale visualization
      eps = 1e-6,
      freq_pre_pc = pmax(freq_pre, eps),
      freq_post_pc = pmax(freq_post, eps),

      # Facet label: Patient ID + Disease + Best Overall Response
      pt_label = paste0(PT_ID, "\n", Disease, " | ", BOR_FC),

      # Determine direction of change
      change_dir = case_when(
        fdr < fdr_thr & freq_post > freq_pre ~ "expansion",
        fdr < fdr_thr & freq_post < freq_pre ~ "contraction",
        TRUE ~ "stable"
      ),

      # Classify clone size
      size_class = ifelse(pmax(freq_pre, freq_post) >= large_thr, "large", "small"),

      # Combined category
      clonal_change = case_when(
        change_dir == "stable" ~ "stable",
        change_dir == "expansion" & size_class == "small" ~ "small clone expansion",
        change_dir == "expansion" & size_class == "large" ~ "large clone expansion",
        change_dir == "contraction" & size_class == "small" ~ "small clone contraction",
        change_dir == "contraction" & size_class == "large" ~ "large clone contraction"
      )
    ) %>%
    mutate(
      clonal_change = factor(
        clonal_change,
        levels = c(
          "large clone contraction",
          "small clone contraction",
          "stable",
          "small clone expansion",
          "large clone expansion"
        )
      )
    )
}

#' Create clonotype tracking scatter plot
#' @param plot_df Prepared plot data frame
#' @param disease_filter Disease to filter for (e.g., "HRSMM", "RRMM")
#' @param x_label X-axis label
#' @param y_label Y-axis label
#' @return ggplot object
create_tracking_plot <- function(plot_df, disease_filter, x_label, y_label) {
  plot_df %>%
    filter(Disease == disease_filter) %>%
    ggplot(aes(x = freq_pre_pc, y = freq_post_pc)) +
    # Reference lines at large clone threshold
    geom_vline(xintercept = LARGE_CLONE_THRESHOLD, linewidth = 0.4) +
    geom_hline(yintercept = LARGE_CLONE_THRESHOLD, linewidth = 0.4) +
    # Diagonal reference (no change)
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50", linewidth = 0.4) +
    # Data points with jitter
    geom_jitter(
      aes(color = clonal_change),
      width = 0.08,
      height = 0.08,
      size = 1.6,
      alpha = 0.5
    ) +
    scale_x_log10(breaks = AXIS_BREAKS, labels = AXIS_LABELS) +
    scale_y_log10(breaks = AXIS_BREAKS, labels = AXIS_LABELS) +
    scale_color_manual(values = CLONE_COLORS, name = "Clonal change category:") +
    facet_wrap(~pt_label, ncol = 4) +
    labs(x = x_label, y = y_label) +
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank()
    )
}

# -----------------------------------------------------------------------------
# 4. TEC Analysis
# -----------------------------------------------------------------------------
message("Loading TEC T cell data...")
T_cells_raw <- qread("~/Immune_project/data_to_share_CART_TEC/T_cells_nov25_deid.qs")

# Prepare and filter data
T_cells_tec <- prepare_clonotype_data(T_cells_raw)
T_cells_tec <- subset(T_cells_tec, Therapy == "TEC" & Tissue == "PB")

# Extract metadata for analysis
meta_tec <- T_cells_tec@meta.data %>%
  select(PT_ID, combined_clonotype, Timepoint)

# Run Poisson analysis
message("Running Poisson tests for TEC patients...")
res_tec <- run_patient_poisson_analysis(meta_tec, pre_label = "Pre", post_label = "Post")

# Load clinical metadata for annotation
output_dir_car <- "~/Immune_project/Paper_Figures_June/"
date_of_sampling <- read_xlsx(paste0(output_dir_car, "Swimmer_Plot_TEC.xlsx"))
deidentified <- read_csv("~/Immune_project/data_to_share_CART_TEC/Dictionary_for_deidentifying_seurat_object.csv") %>%
  filter(Therapy %in% c("Healthy", "TEC", "TALQ")) %>%
  select(ID, PT_ID) %>%
  distinct(PT_ID, .keep_all = TRUE)

date_of_sampling <- date_of_sampling %>%
  left_join(deidentified, by = "ID")

# Get paired patient IDs
paired_tec <- meta_tec %>%
  group_by(PT_ID) %>%
  summarise(n = n_distinct(Timepoint), .groups = "drop") %>%
  filter(n > 1) %>%
  pull(PT_ID)

paired_meta_tec <- date_of_sampling %>%
  filter(PT_ID %in% paired_tec)

# Add clinical annotations to results
res_tec <- res_tec %>%
  left_join(paired_meta_tec[c("PT_ID", "Disease", "BOR_FC")], by = "PT_ID") %>%
  rename(freq_pre = p1, freq_post = p2) %>%
  mutate(
    Disease = replace_na(Disease, "NA"),
    BOR_FC = replace_na(BOR_FC, "NA")
  )

# Save TEC results
write_csv(res_tec, paste0(output_dir, "Poisson_TEC.csv"))

# Generate and save TEC plots
message("Creating TEC plots...")
plot_tec <- prepare_plot_data(res_tec)

p_tec_smm <- create_tracking_plot(
  plot_tec, "HRSMM",
  "Relative abundance pre TEC (log10)",
  "Relative abundance post TEC (log10)"
)
ggsave(p_tec_smm, filename = paste0(plot_dir, "clonotype_tracking_smm_tec.pdf"), width = 8, height = 6)

p_tec_rrmm <- create_tracking_plot(
  plot_tec, "RRMM",
  "Relative abundance pre TEC (log10)",
  "Relative abundance post TEC (log10)"
)
ggsave(p_tec_rrmm, filename = paste0(plot_dir, "clonotype_tracking_rrmm_tec.pdf"), width = 8, height = 6)

# -----------------------------------------------------------------------------
# 5. CAR-T (Cilta-cel) Analysis
# -----------------------------------------------------------------------------
message("Loading CAR-T endogenous T cell data...")
endogenous_t_cells <- qread("~/Immune_project/data_to_share_CART_TEC/endogenous_t_cells_deid.qs")

# Prepare and filter data
T_cells_cart <- prepare_clonotype_data(endogenous_t_cells)
T_cells_cart <- subset(T_cells_cart, Product == "Cilta-cel" & Tissue == "PB" & Timepoint %in% c("Pre", "Day30"))

# Extract metadata for analysis
meta_cart <- T_cells_cart@meta.data %>%
  select(PT_ID, combined_clonotype, Timepoint)

# Run Poisson analysis
message("Running Poisson tests for CAR-T patients...")
res_cart <- run_patient_poisson_analysis(meta_cart, pre_label = "Pre", post_label = "Day30")

# Add clinical annotations from Seurat metadata
cart_clinical <- T_cells_cart@meta.data %>%
  distinct(PT_ID, .keep_all = TRUE) %>%
  select(PT_ID, Disease, BOR_FC)

res_cart <- res_cart %>%
  left_join(cart_clinical, by = "PT_ID") %>%
  rename(freq_pre = p1, freq_post = p2) %>%
  mutate(
    Disease = replace_na(Disease, "NA"),
    BOR_FC = replace_na(BOR_FC, "NA")
  )

# Save CAR-T results
write_csv(res_cart, paste0(output_dir, "Poisson_CART.csv"))

# Generate and save CAR-T plots
message("Creating CAR-T plots...")
plot_cart <- prepare_plot_data(res_cart)

p_cart_smm <- create_tracking_plot(
  plot_cart, "HRSMM",
  "Relative abundance pre Cilta (log10)",
  "Relative abundance post Cilta (log10)"
)
ggsave(p_cart_smm, filename = paste0(plot_dir, "clonotype_tracking_smm_CAR.pdf"), width = 8, height = 6)

p_cart_rrmm <- create_tracking_plot(
  plot_cart, "RRMM",
  "Relative abundance pre Cilta (log10)",
  "Relative abundance post Cilta (log10)"
)
ggsave(p_cart_rrmm, filename = paste0(plot_dir, "clonotype_tracking_rrmm_car.pdf"), width = 8, height = 6)

# -----------------------------------------------------------------------------
# 6. Combined Analysis: Expanding Clone Statistics
# -----------------------------------------------------------------------------
message("Performing combined expansion analysis...")

# Combine TEC and CAR-T results
res_combined <- bind_rows(res_tec, res_cart) %>%
  mutate(
    category = case_when(
      sig == TRUE & freq_post > freq_pre ~ "expanding",
      sig == TRUE & freq_post < freq_pre ~ "contracting",
      sig == FALSE ~ "stable",
      TRUE ~ NA_character_
    )
  )

# Reload TEC data for cell-level analysis
T_cells_tec_full <- prepare_clonotype_data(T_cells_raw)
T_cells_tec_full <- subset(T_cells_tec_full, Therapy == "TEC" & Tissue == "PB")

# Filter for paired patients only
paired_ids_tec <- T_cells_tec_full@meta.data %>%
  group_by(PT_ID) %>%
  summarise(n = n_distinct(Timepoint), .groups = "drop") %>%
  filter(n > 1) %>%
  pull(PT_ID)

T_cells_tec_full <- subset(T_cells_tec_full, PT_ID %in% paired_ids_tec)

# Annotate cells with expansion category
res_for_join <- res_combined %>% select(-Disease)
T_cells_tec_full@meta.data <- T_cells_tec_full@meta.data %>%
  left_join(res_for_join, by = c("PT_ID", "combined_clonotype"))

# Analyze post-treatment samples only
post_tec <- T_cells_tec_full@meta.data %>%
  filter(Timepoint == "Post")

# -----------------------------------------------------------------------------
# 7. Patient-Level Expansion Analysis
# -----------------------------------------------------------------------------
message("Analyzing patient-level expansion patterns...")

# Calculate percentage of patients with at least one expanding clone
pt_expansion_status <- post_tec %>%
  group_by(PT_ID, Disease) %>%
  summarise(
    has_expanding = any(category == "expanding" & n() > 0),
    .groups = "drop"
  )

# Create contingency table
expansion_table <- pt_expansion_status %>%
  count(Disease, has_expanding) %>%
  pivot_wider(names_from = has_expanding, values_from = n, values_fill = 0) %>%
  arrange(Disease)

# Fisher's exact test
expansion_matrix <- as.matrix(expansion_table[, c("FALSE", "TRUE")])
rownames(expansion_matrix) <- expansion_table$Disease
colnames(expansion_matrix) <- c("no_expanding", "has_expanding")

fisher_result <- fisher.test(expansion_matrix)
message(sprintf("Fisher's exact test p-value: %.4f", fisher_result$p.value))

# Create bar plot of expansion status by disease
expansion_plot_data <- as.data.frame(expansion_matrix) %>%
  rownames_to_column("Disease") %>%
  pivot_longer(
    cols = c(no_expanding, has_expanding),
    names_to = "status",
    values_to = "n"
  ) %>%
  group_by(Disease) %>%
  mutate(perc = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    status = factor(
      status,
      levels = c("no_expanding", "has_expanding"),
      labels = c("Not Expanding", "Expanding")
    )
  )

p_expansion_bar <- ggplot(expansion_plot_data, aes(x = Disease, y = perc, fill = status)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  scale_fill_manual(values = c("Not Expanding" = "lightgrey", "Expanding" = "steelblue")) +
  labs(y = "% of Patients", x = NULL) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme_classic(base_size = 22) +
  theme(
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

ggsave(p_expansion_bar, filename = paste0(plot_dir, "n_pts_expanding_TEC.pdf"), width = 6, height = 6)

# -----------------------------------------------------------------------------
# 8. CD8+ T Cell Expansion Analysis
# -----------------------------------------------------------------------------
message("Analyzing CD8+ T cell expansion...")

# Calculate fraction of cells in expanding clones per patient and lineage
lineage_expansion <- post_tec %>%
  group_by(PT_ID, Disease, lineage, category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PT_ID, lineage) %>%
  mutate(
    tot = sum(n),
    fraction = n / tot
  ) %>%
  ungroup()

# Filter to expanding CD8 cells for boxplot
cd8_expanding <- lineage_expansion %>%
  filter(category == "expanding" & lineage == "CD8") %>%
  mutate(
    Disease = factor(Disease, levels = c("HRSMM", "RRMM")),
    jitter_x = as.numeric(Disease) + runif(n(), -0.12, 0.12)
  )

p_cd8_expansion <- ggplot(cd8_expanding, aes(x = Disease, y = fraction * 100, fill = Disease)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.7) +
  geom_point(
    aes(x = jitter_x),
    fill = "darkgrey",
    alpha = 0.5,
    shape = 21,
    size = 4,
    stroke = 0.6,
    color = "black"
  ) +
  geom_pwc(method = "wilcox.test", method.args = list(alternative = "greater"), label.size = 5) +
  scale_fill_manual(values = c("HRSMM" = "lightgrey", "RRMM" = "gold")) +
  labs(y = "% of CD8+ T cells", x = NULL, title = "Cells in an expanding clone") +
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 13),
    legend.position = "none"
  )

ggsave(p_cd8_expansion, filename = paste0(plot_dir, "CD8_expanding.pdf"), width = 6, height = 4)

# Save lineage expansion data
write_csv(lineage_expansion, "~/Immune_project/data_to_share_CART_TEC/clonotype_tracking_TEC.csv")

# -----------------------------------------------------------------------------
# 9. Phenotype Composition of Expanding Clones
# -----------------------------------------------------------------------------
message("Analyzing phenotype composition of expanding clones...")

# Identify patients with at least one expanding clone
patients_with_expansion <- post_tec %>%
  filter(category == "expanding") %>%
  distinct(PT_ID) %>%
  pull(PT_ID)

# Count cells per phenotype in expanding clones
expanding_phenotype <- post_tec %>%
  filter(category == "expanding", PT_ID %in% patients_with_expansion) %>%
  count(PT_ID, final_annotation, name = "n")

# Get all phenotype annotations present in expanding cells
all_annotations <- sort(unique(expanding_phenotype$final_annotation))

# Complete with 0s ONLY for patients who have expanding clones
# Rationale: if a patient has expanding cells but none of phenotype X,
# the percentage of phenotype X in their expanding population is truly 0%
expanding_phenotype <- expanding_phenotype %>%
  complete(
    PT_ID = patients_with_expansion,
    final_annotation = all_annotations,
    fill = list(n = 0)
  ) %>%
  left_join(post_tec %>% distinct(PT_ID, Disease), by = "PT_ID") %>%
  group_by(PT_ID) %>%
  mutate(
    total = sum(n),
    percentage = (n / total) * 100
  ) %>%
  ungroup() %>%
  mutate(Disease = factor(Disease, levels = c("HRSMM", "RRMM")))

# Create heatmap matrix
heatmap_matrix <- expanding_phenotype %>%
  select(PT_ID, final_annotation, percentage) %>%
  pivot_wider(names_from = final_annotation, values_from = percentage) %>%
  column_to_rownames("PT_ID") %>%
  as.matrix()

# Row annotation with disease
row_annotation <- expanding_phenotype %>%
  distinct(PT_ID, Disease) %>%
  column_to_rownames("PT_ID")
row_annotation <- row_annotation[rownames(heatmap_matrix), , drop = FALSE]

# Sort rows by disease
row_order <- order(row_annotation$Disease)
heatmap_matrix <- heatmap_matrix[row_order, , drop = FALSE]
row_annotation <- row_annotation[row_order, , drop = FALSE]

# Annotation colors
annotation_colors <- list(
  Disease = c(HRSMM = "grey70", RRMM = "gold")
)

# Generate heatmap
pdf(paste0(plot_dir, "expanding_phenotype_heatmap.pdf"), width = 10, height = 8)
pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  border_color = NA,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Phenotype composition of expanding clones (%)"
)
dev.off()

# Boxplot comparing phenotype percentages between disease groups
# with Wilcoxon test for each phenotype
p_phenotype_box <- ggplot(
  expanding_phenotype,
  aes(x = final_annotation, y = percentage, fill = Disease)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_point(
    aes(group = Disease),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 2,
    alpha = 0.6,
    shape = 21,
    color = "black"
  ) +
  stat_compare_means(
    aes(group = Disease),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE,
    label.y.npc = 0.95,
    size = 4
  ) +
  scale_fill_manual(values = c("HRSMM" = "grey70", "RRMM" = "gold")) +
  labs(
    x = "Cell Phenotype",
    y = "% of Expanding Cells",
    title = "Phenotype composition of expanding clones by disease"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  p_phenotype_box,
  filename = paste0(plot_dir, "expanding_phenotype_boxplot_by_disease.pdf"),
  width = 12,
  height = 6
)

# -----------------------------------------------------------------------------
# 10. Multivariate Analysis: Overall Compositional Difference
# -----------------------------------------------------------------------------
message("Performing multivariate compositional analysis...")

# PERMANOVA: Test if overall phenotype composition differs by disease
# This tests the hypothesis that centroids and/or dispersion differ between groups
permanova_res <- adonis2(
  heatmap_matrix ~ Disease,
  data = row_annotation,
  method = "euclidean",
  permutations = 999
)

message("PERMANOVA results:")
print(permanova_res)

# Save PERMANOVA results
permanova_df <- as.data.frame(permanova_res)
write_csv(permanova_df, paste0(output_dir, "permanova_phenotype_composition.csv"))

# Test for homogeneity of dispersion (assumption check for PERMANOVA)
# If significant, PERMANOVA result may be driven by dispersion rather than location
dist_matrix <- dist(heatmap_matrix, method = "euclidean")
betadisper_res <- betadisper(dist_matrix, row_annotation$Disease)
permutest_disp <- permutest(betadisper_res, permutations = 999)

message("Homogeneity of dispersion test (betadisper):")
print(permutest_disp)

# -----------------------------------------------------------------------------
# 11. PCA Visualization of Phenotype Composition
# -----------------------------------------------------------------------------
message("Performing PCA on phenotype composition...")

# Remove zero-variance columns (phenotypes absent in all patients)
col_vars <- apply(heatmap_matrix, 2, var)
zero_var_cols <- names(col_vars[col_vars == 0])
if (length(zero_var_cols) > 0) {
  message(sprintf("Removing %d zero-variance phenotypes: %s",
                  length(zero_var_cols), paste(zero_var_cols, collapse = ", ")))
}
heatmap_matrix_pca <- heatmap_matrix[, col_vars > 0, drop = FALSE]

# PCA on the filtered phenotype composition matrix
pca_res <- prcomp(heatmap_matrix_pca, scale. = TRUE, center = TRUE)

# Extract variance explained
var_explained <- summary(pca_res)$importance[2, ] * 100

# Create PCA data frame for plotting
pca_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("PT_ID") %>%
  left_join(row_annotation %>% rownames_to_column("PT_ID"), by = "PT_ID")

# PCA scatter plot
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = Disease)) +
  geom_point(size = 4, shape = 21, color = "black", alpha = 0.8) +
  stat_ellipse(aes(color = Disease), level = 0.95, linewidth = 1, linetype = "dashed") +
  scale_fill_manual(values = c("HRSMM" = "grey70", "RRMM" = "gold")) +
  scale_color_manual(values = c("HRSMM" = "grey50", "RRMM" = "darkgoldenrod")) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    title = "PCA of expanding clone phenotype composition"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(p_pca, filename = paste0(plot_dir, "pca_phenotype_composition.pdf"), width = 7, height = 6)

# PCA biplot showing which phenotypes drive separation
# Extract loadings
loadings_df <- as.data.frame(pca_res$rotation[, 1:2]) %>%
  rownames_to_column("phenotype") %>%
  mutate(
    # Scale loadings for visualization
    PC1_scaled = PC1 * max(abs(pca_df$PC1)) * 0.8,
    PC2_scaled = PC2 * max(abs(pca_df$PC2)) * 0.8
  )

p_biplot <- ggplot() +
  # Sample points
  geom_point(
    data = pca_df,
    aes(x = PC1, y = PC2, fill = Disease),
    size = 4, shape = 21, color = "black", alpha = 0.7
  ) +
  stat_ellipse(
    data = pca_df,
    aes(x = PC1, y = PC2, color = Disease),
    level = 0.95, linewidth = 1, linetype = "dashed"
  ) +
  # Loading arrows
  geom_segment(
    data = loadings_df,
    aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "red", alpha = 0.7
  ) +
  # Loading labels
  geom_text(
    data = loadings_df,
    aes(x = PC1_scaled * 1.1, y = PC2_scaled * 1.1, label = phenotype),
    size = 3, color = "red"
  ) +
  scale_fill_manual(values = c("HRSMM" = "grey70", "RRMM" = "gold")) +
  scale_color_manual(values = c("HRSMM" = "grey50", "RRMM" = "darkgoldenrod")) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    title = "PCA biplot: phenotypes driving disease separation"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(p_biplot, filename = paste0(plot_dir, "pca_biplot_phenotype.pdf"), width = 9, height = 7)

# -----------------------------------------------------------------------------
# 12. Per-Phenotype Wilcoxon Tests with FDR Correction
# -----------------------------------------------------------------------------
message("Performing per-phenotype Wilcoxon tests with FDR correction...")

# Run Wilcoxon test for each phenotype
wilcox_results <- expanding_phenotype %>%
  group_by(final_annotation) %>%
  summarise(
    median_HRSMM = median(percentage[Disease == "HRSMM"], na.rm = TRUE),
    median_RRMM = median(percentage[Disease == "RRMM"], na.rm = TRUE),
    n_HRSMM = sum(Disease == "HRSMM"),
    n_RRMM = sum(Disease == "RRMM"),
    p_value = tryCatch(
      wilcox.test(
        percentage[Disease == "HRSMM"],
        percentage[Disease == "RRMM"]
      )$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    significant = fdr < 0.05,
    direction = case_when(
      median_RRMM > median_HRSMM ~ "Higher in RRMM",
      median_RRMM < median_HRSMM ~ "Higher in HRSMM",
      TRUE ~ "No difference"
    )
  ) %>%
  arrange(p_value)

message("Per-phenotype Wilcoxon test results (FDR-corrected):")
print(wilcox_results)

# Save results
write_csv(wilcox_results, paste0(output_dir, "wilcox_phenotype_by_disease.csv"))

# Create summary plot with FDR-corrected significance
wilcox_for_plot <- wilcox_results %>%
  mutate(
    sig_label = case_when(
      fdr < 0.001 ~ "***",
      fdr < 0.01 ~ "**",
      fdr < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Boxplot with FDR-corrected p-values annotated
p_phenotype_fdr <- ggplot(
  expanding_phenotype,
  aes(x = final_annotation, y = percentage, fill = Disease)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_point(
    aes(group = Disease),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 2,
    alpha = 0.6,
    shape = 21,
    color = "black"
  ) +
  # Add FDR-corrected significance labels
  geom_text(
    data = wilcox_for_plot,
    aes(x = final_annotation, y = max(expanding_phenotype$percentage) * 1.05, label = sig_label),
    inherit.aes = FALSE,
    size = 5
  ) +
  scale_fill_manual(values = c("HRSMM" = "grey70", "RRMM" = "gold")) +
  labs(
    x = "Cell Phenotype",
    y = "% of Expanding Cells",
    title = "Phenotype composition by disease (FDR-corrected)",
    caption = "Significance: * FDR<0.05, ** FDR<0.01, *** FDR<0.001"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(hjust = 0, face = "italic")
  )

ggsave(
  p_phenotype_fdr,
  filename = paste0(plot_dir, "expanding_phenotype_boxplot_FDR.pdf"),
  width = 12,
  height = 6
)

# -----------------------------------------------------------------------------
# 13. Summary Statistics
# -----------------------------------------------------------------------------
message("\n========== SUMMARY ==========")
message(sprintf("Patients with expanding clones: %d", length(patients_with_expansion)))
message(sprintf("  HRSMM: %d", sum(row_annotation$Disease == "HRSMM")))
message(sprintf("  RRMM: %d", sum(row_annotation$Disease == "RRMM")))
message(sprintf("\nPERMANOVA p-value: %.4f", permanova_res$`Pr(>F)`[1]))
message(sprintf("Dispersion test p-value: %.4f", permutest_disp$tab$`Pr(>F)`[1]))
message(sprintf("\nSignificant phenotypes (FDR < 0.05): %d", sum(wilcox_results$significant, na.rm = TRUE)))

if (any(wilcox_results$significant, na.rm = TRUE)) {
  sig_phenos <- wilcox_results %>% filter(significant) %>% pull(final_annotation)
  message(sprintf("  %s", paste(sig_phenos, collapse = ", ")))
}

message("=============================\n")

# -----------------------------------------------------------------------------
# 14. Clonotype-Level Analysis: Dominant Phenotype per Expanding Clone
# -----------------------------------------------------------------------------
message("Performing clonotype-level analysis (dominant phenotype)...")

# For each expanding clonotype, determine its dominant (most common) phenotype
# This gives each clonotype equal weight, regardless of size
clonotype_dominant <- post_tec %>%
  filter(category == "expanding") %>%
  count(PT_ID, Disease, combined_clonotype, final_annotation, name = "n_cells") %>%
  group_by(PT_ID, combined_clonotype) %>%
  slice_max(n_cells, n = 1, with_ties = FALSE) %>%  # Keep most common phenotype
  ungroup() %>%
  rename(dominant_phenotype = final_annotation)

message(sprintf("Total expanding clonotypes: %d", nrow(clonotype_dominant)))
message(sprintf("  HRSMM: %d clonotypes", sum(clonotype_dominant$Disease == "HRSMM")))
message(sprintf("  RRMM: %d clonotypes", sum(clonotype_dominant$Disease == "RRMM")))

# Count clonotypes per dominant phenotype per disease
clonotype_by_phenotype <- clonotype_dominant %>%
  count(Disease, dominant_phenotype, name = "n_clonotypes") %>%
  group_by(Disease) %>%
  mutate(
    total_clonotypes = sum(n_clonotypes),
    percentage = (n_clonotypes / total_clonotypes) * 100
  ) %>%
  ungroup()

# Ensure all phenotypes are represented for both diseases
all_dom_phenotypes <- unique(clonotype_by_phenotype$dominant_phenotype)
clonotype_by_phenotype <- clonotype_by_phenotype %>%
  complete(
    Disease = c("HRSMM", "RRMM"),
    dominant_phenotype = all_dom_phenotypes,
    fill = list(n_clonotypes = 0, percentage = 0)
  ) %>%
  group_by(Disease) %>%
  mutate(total_clonotypes = sum(n_clonotypes)) %>%
  ungroup() %>%
  mutate(Disease = factor(Disease, levels = c("HRSMM", "RRMM")))

# Chi-square test: overall association between disease and dominant phenotype
contingency_table <- clonotype_dominant %>%
  count(Disease, dominant_phenotype) %>%
  pivot_wider(names_from = dominant_phenotype, values_from = n, values_fill = 0) %>%
  column_to_rownames("Disease") %>%
  as.matrix()

# Use Fisher's exact test if any expected count < 5, otherwise chi-square
expected_counts <- chisq.test(contingency_table)$expected
if (any(expected_counts < 5)) {
  message("Using Fisher's exact test (some expected counts < 5)...")
  fisher_clonotype <- fisher.test(contingency_table)
  overall_test <- list(
    method = "Fisher's exact test",
    p.value = fisher_clonotype$p.value
  )
} else {
  chisq_clonotype <- chisq.test(contingency_table)
  overall_test <- list(
    method = "Chi-square test",
    p.value = chisq_clonotype$p.value
  )
}

message(sprintf("%s p-value: %.4f", overall_test$method, overall_test$p.value))

# Per-phenotype comparison: Fisher's exact test for each phenotype
# Is phenotype X more likely to be dominant in RRMM vs HRSMM?
phenotype_fisher <- clonotype_dominant %>%
  mutate(Disease = factor(Disease, levels = c("HRSMM", "RRMM"))) %>%
  group_by(dominant_phenotype) %>%
  summarise(
    n_HRSMM = sum(Disease == "HRSMM"),
    n_RRMM = sum(Disease == "RRMM"),
    .groups = "drop"
  ) %>%
  mutate(
    total_HRSMM = sum(n_HRSMM),
    total_RRMM = sum(n_RRMM),
    pct_HRSMM = (n_HRSMM / total_HRSMM) * 100,
    pct_RRMM = (n_RRMM / total_RRMM) * 100
  ) %>%
  rowwise() %>%
  mutate(
    # 2x2 Fisher test: this phenotype vs all others, HRSMM vs RRMM
    p_value = fisher.test(matrix(
      c(n_HRSMM, total_HRSMM - n_HRSMM,
        n_RRMM, total_RRMM - n_RRMM),
      nrow = 2
    ))$p.value
  ) %>%
  ungroup() %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    significant = fdr < 0.05,
    enriched_in = case_when(
      pct_RRMM > pct_HRSMM ~ "RRMM",
      pct_RRMM < pct_HRSMM ~ "HRSMM",
      TRUE ~ "Neither"
    )
  ) %>%
  arrange(p_value)

message("\nPer-phenotype Fisher's exact test (FDR-corrected):")
print(phenotype_fisher %>% select(dominant_phenotype, pct_HRSMM, pct_RRMM, p_value, fdr, enriched_in))

# Save results
write_csv(phenotype_fisher, paste0(output_dir, "clonotype_dominant_phenotype_fisher.csv"))
write_csv(clonotype_by_phenotype, paste0(output_dir, "clonotype_dominant_phenotype_counts.csv"))

# Stacked bar plot: distribution of dominant phenotypes by disease
p_clonotype_bar <- ggplot(
  clonotype_by_phenotype,
  aes(x = Disease, y = percentage, fill = dominant_phenotype)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  labs(
    x = NULL,
    y = "% of Expanding Clonotypes",
    fill = "Dominant Phenotype",
    title = "Dominant phenotype of expanding clonotypes",
    subtitle = sprintf("%s p = %.3f", overall_test$method, overall_test$p.value)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(color = "black", size = 12)
  )

ggsave(
  p_clonotype_bar,
  filename = paste0(plot_dir, "clonotype_dominant_phenotype_stacked.pdf"),
  width = 8,
  height = 6
)

# Side-by-side bar plot for easier comparison
p_clonotype_dodge <- ggplot(
  clonotype_by_phenotype,
  aes(x = dominant_phenotype, y = percentage, fill = Disease)
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  # Add significance stars from Fisher tests
  geom_text(
    data = phenotype_fisher %>%
      mutate(
        sig_label = case_when(
          fdr < 0.001 ~ "***",
          fdr < 0.01 ~ "**",
          fdr < 0.05 ~ "*",
          TRUE ~ ""
        ),
        y_pos = pmax(pct_HRSMM, pct_RRMM) + 2
      ),
    aes(x = dominant_phenotype, y = y_pos, label = sig_label),
    inherit.aes = FALSE,
    size = 5
  ) +
  scale_fill_manual(values = c("HRSMM" = "grey70", "RRMM" = "gold")) +
  labs(
    x = "Dominant Phenotype",
    y = "% of Expanding Clonotypes",
    title = "Dominant phenotype distribution by disease",
    caption = "* FDR<0.05, ** FDR<0.01, *** FDR<0.001"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(hjust = 0, face = "italic")
  )

ggsave(
  p_clonotype_dodge,
  filename = paste0(plot_dir, "clonotype_dominant_phenotype_dodge.pdf"),
  width = 10,
  height = 6
)

# Patient-level analysis: percentage of each dominant phenotype per patient
clonotype_patient_level <- clonotype_dominant %>%
  count(PT_ID, dominant_phenotype, name = "n_clonotypes") %>%
  group_by(PT_ID) %>%
  mutate(
    total = sum(n_clonotypes),
    percentage = (n_clonotypes / total) * 100
  ) %>%
  ungroup()

# Complete with zeros for missing phenotypes per patient
patients_in_analysis <- unique(clonotype_patient_level$PT_ID)
clonotype_patient_level <- clonotype_patient_level %>%
  complete(
    PT_ID = patients_in_analysis,
    dominant_phenotype = all_dom_phenotypes,
    fill = list(n_clonotypes = 0, percentage = 0)
  ) %>%
  left_join(
    clonotype_dominant %>% distinct(PT_ID, Disease),
    by = "PT_ID"
  ) %>%
  mutate(Disease = factor(Disease, levels = c("HRSMM", "RRMM")))

# Boxplot at patient level with Wilcoxon tests
p_clonotype_patient <- ggplot(
  clonotype_patient_level,
  aes(x = dominant_phenotype, y = percentage, fill = Disease)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_point(
    aes(group = Disease),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 2,
    alpha = 0.6,
    shape = 21,
    color = "black"
  ) +
  stat_compare_means(
    aes(group = Disease),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE,
    label.y.npc = 0.95,
    size = 4
  ) +
  scale_fill_manual(values = c("HRSMM" = "grey70", "RRMM" = "gold")) +
  labs(
    x = "Dominant Phenotype",
    y = "% of Patient's Expanding Clonotypes",
    title = "Clonotype-level: dominant phenotype by disease (patient-level)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  p_clonotype_patient,
  filename = paste0(plot_dir, "clonotype_dominant_phenotype_patient_boxplot.pdf"),
  width = 12,
  height = 6
)

message("\nClonotype-level analysis complete!")

# -----------------------------------------------------------------------------
# 15. UMAP Visualization of Expanding/Contracting/Stable Clonotypes
# -----------------------------------------------------------------------------
message("Generating UMAP plots for clonotype expansion categories...")

# Colors matching the clonotype tracking plot theme
CATEGORY_COLORS <- c(
  "expanding" = "#D81B60",
  "contracting" = "darkblue",
  "stable" = "grey80"
)

# Extract UMAP coordinates from the Seurat object
# T_cells_tec_full already has the 'category' column from the join with res_combined
umap_coords <- as.data.frame(T_cells_tec_full@reductions$harmony_umap_t@cell.embeddings)
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$cell_id <- rownames(umap_coords)

# Add metadata
umap_df <- umap_coords %>%
  left_join(
    T_cells_tec_full@meta.data %>%
      rownames_to_column("cell_id") %>%
      select(cell_id, PT_ID, Disease, Timepoint, category, final_annotation, lineage),
    by = "cell_id"
  ) %>%
  filter(!is.na(category)) %>%
  mutate(
    category = factor(category, levels = c("expanding", "stable", "contracting")),
    Disease = factor(Disease, levels = c("HRSMM", "RRMM"))
  )

# Function to create UMAP plot for a specific disease
create_category_umap <- function(data, disease_name, split_by_timepoint = TRUE) {

  df <- data %>% filter(Disease == disease_name)

  if (split_by_timepoint) {
    # Order: stable first (background), then contracting, then expanding (on top)
    df <- df %>%
      arrange(factor(category, levels = c("stable", "contracting", "expanding")))

    p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = category)) +
      geom_point(size = 0.5, alpha = 0.6) +
      scale_color_manual(
        values = CATEGORY_COLORS,
        name = "Clonal dynamics:",
        drop = FALSE
      ) +
      facet_wrap(~Timepoint, ncol = 2) +
      labs(
        x = "UMAP-1",
        y = "UMAP-2",
        title = disease_name
      ) +
      theme_classic(base_size = 14) +
      theme(
        legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.text = element_text(size = 10),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5)
      ) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  } else {
    # Single UMAP without timepoint split
    df <- df %>%
      arrange(factor(category, levels = c("stable", "contracting", "expanding")))

    p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = category)) +
      geom_point(size = 0.5, alpha = 0.6) +
      scale_color_manual(
        values = CATEGORY_COLORS,
        name = "Clonal dynamics:",
        drop = FALSE
      ) +
      labs(
        x = "UMAP-1",
        y = "UMAP-2",
        title = disease_name
      ) +
      theme_classic(base_size = 14) +
      theme(
        legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.text = element_text(size = 10),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5)
      ) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  }

  return(p)
}

# Create UMAPs for HRSMM and RRMM, split by timepoint
p_umap_hrsmm <- create_category_umap(umap_df, "HRSMM", split_by_timepoint = TRUE)
p_umap_rrmm <- create_category_umap(umap_df, "RRMM", split_by_timepoint = TRUE)

# Save individual plots
ggsave(
  p_umap_hrsmm,
  filename = paste0(plot_dir, "umap_clonotype_category_HRSMM.pdf"),
  width = 8,
  height = 5
)

ggsave(
  p_umap_rrmm,
  filename = paste0(plot_dir, "umap_clonotype_category_RRMM.pdf"),
  width = 8,
  height = 5
)

# Combined plot with both diseases
p_umap_combined <- (p_umap_hrsmm + theme(legend.position = "none")) /
  (p_umap_rrmm + theme(legend.position = "bottom"))

ggsave(
  p_umap_combined,
  filename = paste0(plot_dir, "umap_clonotype_category_combined.pdf"),
  width = 8,
  height = 10
)

# Also create a version highlighting only expanding cells
umap_df_highlight <- umap_df %>%
  mutate(
    highlight = ifelse(category == "expanding", "expanding", "other"),
    highlight = factor(highlight, levels = c("other", "expanding"))
  ) %>%
  arrange(highlight)  # Plot "other" first, expanding on top

HIGHLIGHT_COLORS <- c("other" = "grey85", "expanding" = "#D81B60")

create_highlight_umap <- function(data, disease_name) {
  df <- data %>%
    filter(Disease == disease_name) %>%
    arrange(highlight)

  ggplot(df, aes(x = UMAP1, y = UMAP2, color = highlight)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_manual(
      values = HIGHLIGHT_COLORS,
      name = "",
      labels = c("other" = "Other", "expanding" = "Expanding")
    ) +
    facet_wrap(~Timepoint, ncol = 2) +
    labs(
      x = "UMAP-1",
      y = "UMAP-2",
      title = disease_name
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "bottom",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

p_highlight_hrsmm <- create_highlight_umap(umap_df_highlight, "HRSMM")
p_highlight_rrmm <- create_highlight_umap(umap_df_highlight, "RRMM")

p_highlight_combined <- (p_highlight_hrsmm + theme(legend.position = "none")) /
  (p_highlight_rrmm + theme(legend.position = "bottom"))

ggsave(
  p_highlight_combined,
  filename = paste0(plot_dir, "umap_expanding_highlighted_combined.pdf"),
  width = 8,
  height = 10
)

# Print summary of cells per category per disease/timepoint
umap_summary <- umap_df %>%
  group_by(Disease, Timepoint, category) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = category, values_from = n, values_fill = 0)

message("\nUMAP cell counts by category:")
print(umap_summary)

message("\nUMAP plots complete!")

# -----------------------------------------------------------------------------
# 16. Balanced DEG Analysis: Expanding Cells RRMM vs HRSMM (Post-treatment)
# -----------------------------------------------------------------------------
message("Performing balanced DEG analysis: Expanding cells RRMM vs HRSMM...")

# Load required library for DEG
library(SingleCellExperiment)
library(ggrepel)

# Get post-treatment expanding cells only
expanding_post <- post_tec %>%
  filter(category == "expanding") %>%
  filter(Disease %in% c("HRSMM", "RRMM"))

# Check we have cells from both groups
n_hrsmm <- sum(expanding_post$Disease == "HRSMM")
n_rrmm <- sum(expanding_post$Disease == "RRMM")
message(sprintf("Expanding cells: HRSMM = %d, RRMM = %d", n_hrsmm, n_rrmm))

if (n_hrsmm > 0 && n_rrmm > 0) {

  # Subset the Seurat object to expanding post-treatment cells
  expanding_barcodes <- rownames(expanding_post)
  T_cells_expanding <- subset(T_cells_tec_full, cells = expanding_barcodes)

  # Create SingleCellExperiment object
  meta_exp <- T_cells_expanding@meta.data
  meta_exp$Condition <- factor(meta_exp$Disease, levels = c("HRSMM", "RRMM"))

  sce_exp <- SingleCellExperiment(
    assays = list(
      counts = T_cells_expanding@assays$RNA@counts,
      data = T_cells_expanding@assays$RNA@data
    ),
    colData = meta_exp
  )

  # Filter genes: keep genes with >3 counts in at least 10 cells
  sce_exp <- sce_exp[rowSums(counts(sce_exp) > 3) >= 10, ]
  keep_genes <- rownames(sce_exp)

  # Remove non-coding genes, mitochondrial, ribosomal, IG/TCR
  noncoding_idx <- grep(
    "^RP[SL]|^MT-|^IG[HKL][VDJCGAMDE]|^TR[ABGD][VDJC]|-AS1$|orf|^ENSG",
    keep_genes,
    ignore.case = TRUE
  )
  if (length(noncoding_idx) > 0) {
    sce_exp <- sce_exp[-noncoding_idx, ]
  }

  message(sprintf("Genes after filtering: %d", nrow(sce_exp)))

  # -----------------------------------------------------------------------------
  # Balanced sampling
  # -----------------------------------------------------------------------------
  set.seed(42)
  min_n <- 50  # Minimum cells per patient (lower threshold for expanding cells)

  stmp <- meta_exp

  # Count cells per patient
  pt_n <- data.frame(table(stmp$PT_ID))
  colnames(pt_n) <- c("PT_ID", "Freq")

  # Remove patients with fewer than min_n cells
  remove_pts <- as.character(pt_n$PT_ID[pt_n$Freq < min_n])
  stmp <- stmp[!stmp$PT_ID %in% remove_pts, ]

  # Check valid patients
  valid_pts <- pt_n$Freq[pt_n$Freq >= min_n]

  if (length(valid_pts) >= 2 && length(unique(stmp$Condition)) == 2) {

    new_min <- min(valid_pts)

    # Count cells per disease group
    group_df <- data.frame(table(stmp$Condition))
    colnames(group_df) <- c("Condition", "Freq")

    # Count patients per disease group
    stmp_un <- stmp[!duplicated(stmp$PT_ID), ]
    pt_df <- data.frame(table(stmp_un$Condition))
    colnames(pt_df) <- c("Condition", "Freq")

    n_pts_hrsmm <- pt_df$Freq[pt_df$Condition == "HRSMM"]
    n_pts_rrmm <- pt_df$Freq[pt_df$Condition == "RRMM"]

    message(sprintf("Patients after filtering: HRSMM = %d, RRMM = %d", n_pts_hrsmm, n_pts_rrmm))

    # Calculate feasible balanced sample sizes
    ceil_cells <- min(group_df$Freq)

    prog_n <- lapply(seq_len(ceil_cells), function(x) x / n_pts_hrsmm)
    non_n <- lapply(seq_len(ceil_cells), function(x) x / n_pts_rrmm)

    is_int <- function(x) is.integer(x) || (is.numeric(x) && identical(round(x), x))

    prog_cand <- which(sapply(prog_n, is_int))
    non_cand <- which(sapply(non_n, is_int))

    prog_cand <- prog_cand[prog_cand / n_pts_hrsmm <= new_min]
    non_cand <- non_cand[non_cand / n_pts_rrmm <= new_min]

    common_samplesize <- intersect(prog_cand, non_cand)

    # Perform balanced sampling
    stmp$CellID <- rownames(stmp)

    if (length(common_samplesize) == 0) {
      # Fallback: sample 30 cells per patient
      fallback_n <- min(30, new_min)
      message(sprintf("Using fallback sampling: %d cells per patient", fallback_n))

      group1_s <- stmp %>%
        filter(Condition == "HRSMM") %>%
        group_by(PT_ID) %>%
        sample_n(size = min(fallback_n, n()))

      group2_s <- stmp %>%
        filter(Condition == "RRMM") %>%
        group_by(PT_ID) %>%
        sample_n(size = min(fallback_n, n()))

      full <- rbind(group1_s, group2_s)
      total_n <- paste0("HRSMM: ", nrow(group1_s), "; RRMM: ", nrow(group2_s))

    } else {
      total_n <- max(common_samplesize)
      message(sprintf("Using balanced sampling: %d cells per group", total_n))

      tmp_n_hrsmm <- total_n / n_pts_hrsmm
      tmp_n_rrmm <- total_n / n_pts_rrmm

      group1_s <- stmp %>%
        filter(Condition == "HRSMM") %>%
        group_by(PT_ID) %>%
        sample_n(size = tmp_n_hrsmm)

      group2_s <- stmp %>%
        filter(Condition == "RRMM") %>%
        group_by(PT_ID) %>%
        sample_n(size = tmp_n_rrmm)

      full <- rbind(group1_s, group2_s)
    }

    rownames(full) <- full$CellID

    message(sprintf("Sampled cells: HRSMM = %d, RRMM = %d",
                    sum(full$Condition == "HRSMM"),
                    sum(full$Condition == "RRMM")))

    # -----------------------------------------------------------------------------
    # Wilcoxon test per gene
    # -----------------------------------------------------------------------------
    cellids <- rownames(full)
    levels <- c("HRSMM", "RRMM")

    # Extract expression matrix
    tmp <- data.frame(t(assay(sce_exp[, cellids], "data")))
    tmp$Condition <- factor(
      as.character(meta_exp$Condition[match(rownames(tmp), rownames(meta_exp))]),
      levels = levels
    )

    message("Running Wilcoxon tests...")

    # Wilcoxon test per gene
    p_values <- apply(tmp[, 1:nrow(sce_exp)], 2, function(x) {
      tryCatch(
        wilcox.test(x ~ tmp$Condition)$p.value,
        error = function(e) NA_real_
      )
    })

    # Mean log2 expression per group
    mean_hrsmm <- apply(
      tmp[tmp$Condition == "HRSMM", 1:nrow(sce_exp)], 2,
      function(x) mean(log2(x + 1), na.rm = TRUE)
    )
    mean_rrmm <- apply(
      tmp[tmp$Condition == "RRMM", 1:nrow(sce_exp)], 2,
      function(x) mean(log2(x + 1), na.rm = TRUE)
    )

    # Build results table
    res_deg <- data.frame(
      Wilcox_p = p_values,
      Mean_Log2_HRSMM = mean_hrsmm,
      Mean_Log2_RRMM = mean_rrmm
    )

    res_deg$Log2FC <- res_deg$Mean_Log2_RRMM - res_deg$Mean_Log2_HRSMM
    res_deg$padj <- p.adjust(res_deg$Wilcox_p, method = "BH")
    res_deg$Gene <- rownames(res_deg)
    res_deg$CellsPerGroup <- total_n
    res_deg$PatientsPerGroup <- paste0("HRSMM: ", n_pts_hrsmm, "; RRMM: ", n_pts_rrmm)

    # Save DEG results
    write_csv(
      res_deg,
      paste0(output_dir, "DEG_expanding_RRMMvsHRSMM_Balanced.csv")
    )

    message(sprintf("DEG results saved: %d genes tested", nrow(res_deg)))

    # -----------------------------------------------------------------------------
    # Volcano plot
    # -----------------------------------------------------------------------------
    message("Creating volcano plot...")

    logfc_thresh <- 0.25
    pval_thresh <- 0.05

    # Annotate significance
    res_deg <- res_deg %>%
      mutate(
        sig = case_when(
          padj < pval_thresh & Log2FC > logfc_thresh ~ "RRMM_up",
          padj < pval_thresh & Log2FC < -logfc_thresh ~ "HRSMM_up",
          TRUE ~ "NS"
        )
      )

    # Count significant genes
    n_rrmm_up <- sum(res_deg$sig == "RRMM_up", na.rm = TRUE)
    n_hrsmm_up <- sum(res_deg$sig == "HRSMM_up", na.rm = TRUE)
    message(sprintf("Significant genes: RRMM_up = %d, HRSMM_up = %d", n_rrmm_up, n_hrsmm_up))

    # Top genes for labeling
    top_up <- res_deg %>%
      filter(sig == "RRMM_up") %>%
      arrange(desc(Log2FC)) %>%
      slice_head(n = 20)

    top_down <- res_deg %>%
      filter(sig == "HRSMM_up") %>%
      arrange(Log2FC) %>%
      slice_head(n = 20)

    top_genes <- bind_rows(top_up, top_down) %>%
      distinct(Gene, .keep_all = TRUE)

    # Volcano plot
    p_volcano <- ggplot(res_deg, aes(x = Log2FC, y = -log10(padj), color = sig)) +
      geom_point(alpha = 0.7, size = 2) +
      scale_color_manual(
        values = c(
          "RRMM_up" = "gold",
          "HRSMM_up" = "grey50",
          "NS" = "grey85"
        ),
        name = "Significance"
      ) +
      geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed", color = "black") +
      geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed", color = "black") +
      geom_text_repel(
        data = top_genes,
        aes(label = Gene),
        size = 3.5,
        max.overlaps = 20,
        segment.size = 0.3,
        color = "black",
        force = 2
      ) +
      labs(
        title = "DEG: Expanding cells RRMM vs HRSMM (Post-treatment)",
        subtitle = sprintf(
          "Cells: HRSMM = %d, RRMM = %d | Patients: HRSMM = %d, RRMM = %d",
          sum(full$Condition == "HRSMM"),
          sum(full$Condition == "RRMM"),
          n_pts_hrsmm,
          n_pts_rrmm
        ),
        x = "log2FC (RRMM vs HRSMM)",
        y = expression(-log[10] * "(adj. p-value)")
      ) +
      theme_classic(base_size = 16) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14),
        legend.position = "top"
      ) +
      guides(color = guide_legend(override.aes = list(size = 4)))

    ggsave(
      p_volcano,
      filename = paste0(plot_dir, "volcano_expanding_RRMMvsHRSMM.pdf"),
      width = 10,
      height = 8
    )

    message("Volcano plot saved!")

  } else {
    message("Not enough patients with sufficient cells for balanced DEG analysis")
  }

} else {
  message("Not enough expanding cells in both disease groups for DEG analysis")
}

message("\nDEG analysis complete!")

# -----------------------------------------------------------------------------
# 17. Signature Analysis: Expanding Cells RRMM vs HRSMM (Sample-Level)
# -----------------------------------------------------------------------------
# NOTE: Comparison is done at the SAMPLE level (mean signature score per sample)
# to account for patient-level variation and avoid pseudo-replication.
# This follows the approach in Step11.BsAbs_workflow.R
# -----------------------------------------------------------------------------
message("Performing signature analysis on expanding cells (sample-level)...")

# Load pre-computed signature scores from SignatureAnalyzer
sig_file <- "~/Immune_project/data_to_share_CART_TEC/dt_cleaned_TEC_deid.csv"

if (file.exists(sig_file)) {

  sig_dt <- read_csv(sig_file, show_col_types = FALSE)

  # Get signature columns (S3, S5, S6, etc.)
  sig_cols <- grep("^S[0-9]+$", colnames(sig_dt), value = TRUE)
  message(sprintf("Found %d signatures: %s", length(sig_cols), paste(sig_cols, collapse = ", ")))

  # Filter to expanding post-treatment cells
  # Match by cell barcode
  expanding_cells <- rownames(expanding_post)

  sig_expanding <- sig_dt %>%
    filter(cell %in% expanding_cells) %>%
    filter(Timepoint == "Post") %>%
    filter(Disease %in% c("HRSMM", "RRMM")) %>%
    mutate(Disease = factor(Disease, levels = c("HRSMM", "RRMM")))

  n_sig_cells <- nrow(sig_expanding)
  message(sprintf("Expanding cells with signature data: %d", n_sig_cells))

  if (n_sig_cells > 0) {

    # -------------------------------------------------------------------------
    # Aggregate to SAMPLE level: mean signature score per sample
    # This is the correct approach for comparing between disease groups
    # -------------------------------------------------------------------------
    message("Aggregating signature scores to sample level...")

    sig_by_sample <- sig_expanding %>%
      group_by(ID, Disease) %>%
      summarise(
        n_cells = n(),
        across(all_of(sig_cols), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      mutate(Disease = factor(Disease, levels = c("HRSMM", "RRMM")))

    n_hrsmm_samples <- sum(sig_by_sample$Disease == "HRSMM")
    n_rrmm_samples <- sum(sig_by_sample$Disease == "RRMM")
    message(sprintf("Samples with expanding cells: HRSMM = %d, RRMM = %d",
                    n_hrsmm_samples, n_rrmm_samples))

    # -------------------------------------------------------------------------
    # Compare signatures between HRSMM and RRMM at SAMPLE level
    # -------------------------------------------------------------------------
    message("Comparing signatures between disease groups (sample-level)...")

    # Wilcoxon test for each signature using sample-level means
    sig_comparison <- map_dfr(sig_cols, function(sig) {
      hrsmm_vals <- sig_by_sample[[sig]][sig_by_sample$Disease == "HRSMM"]
      rrmm_vals <- sig_by_sample[[sig]][sig_by_sample$Disease == "RRMM"]

      # Skip if not enough samples
      if (length(hrsmm_vals) < 3 || length(rrmm_vals) < 3) {
        return(tibble(
          signature = sig,
          median_HRSMM = NA_real_,
          median_RRMM = NA_real_,
          p_value = NA_real_
        ))
      }

      test_res <- tryCatch(
        wilcox.test(hrsmm_vals, rrmm_vals),
        error = function(e) list(p.value = NA_real_)
      )

      tibble(
        signature = sig,
        median_HRSMM = median(hrsmm_vals, na.rm = TRUE),
        median_RRMM = median(rrmm_vals, na.rm = TRUE),
        mean_HRSMM = mean(hrsmm_vals, na.rm = TRUE),
        mean_RRMM = mean(rrmm_vals, na.rm = TRUE),
        n_samples_HRSMM = length(hrsmm_vals),
        n_samples_RRMM = length(rrmm_vals),
        p_value = test_res$p.value
      )
    })

    # FDR correction
    sig_comparison <- sig_comparison %>%
      mutate(
        fdr = p.adjust(p_value, method = "BH"),
        significant = fdr < 0.05,
        log2FC = log2((mean_RRMM + 1) / (mean_HRSMM + 1)),
        enriched_in = case_when(
          significant & mean_RRMM > mean_HRSMM ~ "RRMM",
          significant & mean_RRMM < mean_HRSMM ~ "HRSMM",
          TRUE ~ "NS"
        )
      ) %>%
      arrange(p_value)

    message("\nSignature comparison results (sample-level):")
    print(sig_comparison %>% select(signature, median_HRSMM, median_RRMM, p_value, fdr, enriched_in))

    # Save results
    write_csv(
      sig_comparison,
      paste0(output_dir, "signature_comparison_expanding_RRMMvsHRSMM.csv")
    )

    # -------------------------------------------------------------------------
    # Boxplots for all signatures (sample-level)
    # -------------------------------------------------------------------------
    message("Creating signature boxplots (sample-level)...")

    # Pivot sample-level data to long format for plotting
    sig_long <- sig_by_sample %>%
      pivot_longer(
        cols = all_of(sig_cols),
        names_to = "signature",
        values_to = "score"
      ) %>%
      mutate(signature = factor(signature, levels = sig_cols))

    # Add significance labels
    sig_labels <- sig_comparison %>%
      mutate(
        sig_label = case_when(
          fdr < 0.001 ~ "***",
          fdr < 0.01 ~ "**",
          fdr < 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      select(signature, sig_label, fdr)

    # Overall boxplot of all signatures (sample-level means)
    p_sig_all <- ggplot(sig_long, aes(x = signature, y = score, fill = Disease)) +
      geom_boxplot(
        position = position_dodge(width = 0.8),
        width = 0.7,
        outlier.shape = NA,
        alpha = 0.7
      ) +
      geom_point(
        aes(group = Disease),
        position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
        size = 2,
        alpha = 0.6,
        shape = 21,
        color = "black"
      ) +
      stat_compare_means(
        aes(group = Disease),
        method = "wilcox.test",
        label = "p.signif",
        hide.ns = TRUE,
        label.y.npc = 0.95,
        size = 4
      ) +
      scale_fill_manual(values = c("HRSMM" = "grey70", "RRMM" = "gold")) +
      labs(
        x = "Signature",
        y = "Mean Signature Score (per sample)",
        title = "Signature scores in expanding cells: RRMM vs HRSMM",
        subtitle = sprintf("Sample-level comparison (HRSMM: n=%d, RRMM: n=%d)",
                           n_hrsmm_samples, n_rrmm_samples)
      ) +
      theme_classic(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10)
      ) +
      coord_cartesian(ylim = c(0, quantile(sig_long$score, 0.99, na.rm = TRUE)))

    ggsave(
      p_sig_all,
      filename = paste0(plot_dir, "signature_boxplot_all_expanding.pdf"),
      width = 14,
      height = 6
    )

    # -------------------------------------------------------------------------
    # Focused boxplots for significant signatures
    # -------------------------------------------------------------------------
    sig_significant <- sig_comparison %>%
      filter(significant) %>%
      pull(signature)

    if (length(sig_significant) > 0) {
      message(sprintf("Significant signatures (FDR < 0.05): %s", paste(sig_significant, collapse = ", ")))

      sig_long_sig <- sig_long %>%
        filter(signature %in% sig_significant)

      p_sig_significant <- ggplot(sig_long_sig, aes(x = signature, y = score, fill = Disease)) +
        geom_boxplot(
          position = position_dodge(width = 0.8),
          width = 0.7,
          outlier.shape = NA,
          alpha = 0.7
        ) +
        geom_point(
          aes(group = Disease),
          position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
          size = 2.5,
          alpha = 0.7,
          shape = 21,
          color = "black"
        ) +
        stat_compare_means(
          aes(group = Disease),
          method = "wilcox.test",
          label = "p.format",
          label.y.npc = 0.95,
          size = 3.5
        ) +
        scale_fill_manual(values = c("HRSMM" = "grey70", "RRMM" = "gold")) +
        labs(
          x = "Signature",
          y = "Mean Signature Score (per sample)",
          title = "Significantly different signatures in expanding cells",
          subtitle = "Sample-level comparison"
        ) +
        theme_classic(base_size = 14) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10)
        )

      ggsave(
        p_sig_significant,
        filename = paste0(plot_dir, "signature_boxplot_significant_expanding.pdf"),
        width = max(6, length(sig_significant) * 1.2),
        height = 6
      )
    } else {
      message("No signatures reached significance (FDR < 0.05)")
    }

    # -------------------------------------------------------------------------
    # Heatmap of signature scores by patient (already at sample level)
    # -------------------------------------------------------------------------
    message("Creating signature heatmap...")

    # Create matrix for heatmap (already aggregated by patient)
    sig_matrix <- sig_by_sample %>%
      select(ID, all_of(sig_cols)) %>%
      column_to_rownames("ID") %>%
      as.matrix()

    # Scale signatures for visualization
    sig_matrix_scaled <- scale(sig_matrix)

    # Row annotation
    row_annot <- sig_by_sample %>%
      select(ID, Disease) %>%
      column_to_rownames("ID")

    # Sort by disease
    row_order <- order(row_annot$Disease)
    sig_matrix_scaled <- sig_matrix_scaled[row_order, , drop = FALSE]
    row_annot <- row_annot[row_order, , drop = FALSE]

    # Annotation colors
    annot_colors <- list(
      Disease = c(HRSMM = "grey70", RRMM = "gold")
    )

    pdf(paste0(plot_dir, "signature_heatmap_expanding.pdf"), width = 12, height = 8)
    pheatmap(
      sig_matrix_scaled,
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      annotation_row = row_annot,
      annotation_colors = annot_colors,
      border_color = NA,
      fontsize_row = 8,
      fontsize_col = 9,
      main = "Signature scores in expanding cells (z-scaled, sample means)",
      color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
    )
    dev.off()

    # -------------------------------------------------------------------------
    # Summary barplot: number of signatures enriched in each group
    # -------------------------------------------------------------------------
    enrichment_summary <- sig_comparison %>%
      filter(significant) %>%
      group_by(enriched_in) %>%
      summarise(n_signatures = n(), .groups = "drop")

    if (nrow(enrichment_summary) > 0) {
      p_enrichment <- ggplot(enrichment_summary, aes(x = enriched_in, y = n_signatures, fill = enriched_in)) +
        geom_bar(stat = "identity", width = 0.6, color = "black") +
        geom_text(aes(label = n_signatures), vjust = -0.5, size = 5) +
        scale_fill_manual(values = c("HRSMM" = "grey70", "RRMM" = "gold")) +
        labs(
          x = "Enriched in",
          y = "Number of Signatures",
          title = "Signatures enriched in expanding cells",
          subtitle = "Sample-level comparison"
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_classic(base_size = 14) +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10)
        )

      ggsave(
        p_enrichment,
        filename = paste0(plot_dir, "signature_enrichment_summary.pdf"),
        width = 5,
        height = 5
      )
    }

    message("\nSignature analysis complete (sample-level)!")

  } else {
    message("No expanding cells found with signature data")
  }

} else {
  message(sprintf("Signature file not found: %s", sig_file))
}

message("Analysis complete!")

# -----------------------------------------------------------------------------
# 18. Pseudobulk DEG Analysis: Expanding Cells RRMM vs HRSMM
# -----------------------------------------------------------------------------
# More robust than cell-level analysis - aggregates counts per patient first
# then performs DEG at the patient level (no pseudo-replication)
# -----------------------------------------------------------------------------
message("\nPerforming pseudobulk DEG analysis: Expanding cells RRMM vs HRSMM...")

library(DESeq2)
library(jsonlite)

if (exists("T_cells_expanding") && ncol(T_cells_expanding) > 0) {

  # Get expanding cells metadata
  meta_pb <- T_cells_expanding@meta.data %>%
    filter(Disease %in% c("HRSMM", "RRMM")) %>%
    mutate(Disease = factor(Disease, levels = c("HRSMM", "RRMM")))

  # Check minimum patients per group
  n_pts_per_disease <- meta_pb %>%
    distinct(PT_ID, Disease) %>%
    count(Disease)

  message("Patients per disease group:")
  print(n_pts_per_disease)

  if (all(n_pts_per_disease$n >= 2)) {

    # -------------------------------------------------------------------------
    # Aggregate counts to patient level (pseudobulk)
    # -------------------------------------------------------------------------
    message("Aggregating counts to patient level...")

    # Get raw counts for expanding cells
    counts_mat <- T_cells_expanding@assays$RNA@counts[, rownames(meta_pb)]

    # Create pseudobulk by summing counts per patient
    patients <- unique(meta_pb$PT_ID)
    pseudobulk_counts <- matrix(0, nrow = nrow(counts_mat), ncol = length(patients))
    rownames(pseudobulk_counts) <- rownames(counts_mat)
    colnames(pseudobulk_counts) <- patients

    for (pt in patients) {
      pt_cells <- rownames(meta_pb)[meta_pb$PT_ID == pt]
      if (length(pt_cells) > 1) {
        pseudobulk_counts[, pt] <- rowSums(counts_mat[, pt_cells, drop = FALSE])
      } else if (length(pt_cells) == 1) {
        pseudobulk_counts[, pt] <- as.numeric(counts_mat[, pt_cells])
      }
    }

    # Convert to integer matrix (DESeq2 requirement)
    pseudobulk_counts <- round(pseudobulk_counts)
    storage.mode(pseudobulk_counts) <- "integer"

    # Create sample metadata for DESeq2
    sample_meta <- meta_pb %>%
      distinct(PT_ID, Disease) %>%
      column_to_rownames("PT_ID")
    sample_meta <- sample_meta[colnames(pseudobulk_counts), , drop = FALSE]

    # Count cells per patient (for filtering and info)
    cells_per_pt <- meta_pb %>%
      count(PT_ID, name = "n_cells")
    message("Cells per patient:")
    print(cells_per_pt)

    # Filter patients with at least 20 cells
    min_cells_pb <- 20
    valid_pts <- cells_per_pt %>%
      filter(n_cells >= min_cells_pb) %>%
      pull(PT_ID)

    if (length(valid_pts) >= 4) {  # Need at least 2 per group

      pseudobulk_counts <- pseudobulk_counts[, valid_pts, drop = FALSE]
      sample_meta <- sample_meta[valid_pts, , drop = FALSE]

      # Filter genes: keep genes with at least 10 counts in at least 2 samples
      keep_genes <- rowSums(pseudobulk_counts >= 10) >= 2
      pseudobulk_counts <- pseudobulk_counts[keep_genes, , drop = FALSE]

      # Remove non-coding genes
      gene_names <- rownames(pseudobulk_counts)
      noncoding_idx <- grep(
        "^RP[SL]|^MT-|^IG[HKL][VDJCGAMDE]|^TR[ABGD][VDJC]|-AS1$|orf|^ENSG",
        gene_names,
        ignore.case = TRUE
      )
      if (length(noncoding_idx) > 0) {
        pseudobulk_counts <- pseudobulk_counts[-noncoding_idx, , drop = FALSE]
      }

      message(sprintf("Genes after filtering: %d", nrow(pseudobulk_counts)))
      message(sprintf("Patients after filtering: %d", ncol(pseudobulk_counts)))

      # -----------------------------------------------------------------------
      # DESeq2 analysis
      # -----------------------------------------------------------------------
      message("Running DESeq2...")

      dds <- DESeqDataSetFromMatrix(
        countData = pseudobulk_counts,
        colData = sample_meta,
        design = ~ Disease
      )

      # Run DESeq2
      dds <- DESeq(dds)

      # Get results: RRMM vs HRSMM
      res_pb <- results(dds, contrast = c("Disease", "RRMM", "HRSMM"))
      res_pb <- as.data.frame(res_pb) %>%
        rownames_to_column("Gene") %>%
        arrange(pvalue) %>%
        mutate(
          significant = padj < 0.05 & !is.na(padj),
          direction = case_when(
            significant & log2FoldChange > 0 ~ "RRMM_up",
            significant & log2FoldChange < 0 ~ "HRSMM_up",
            TRUE ~ "NS"
          )
        )

      # Save results
      write_csv(
        res_pb,
        paste0(output_dir, "DEG_expanding_RRMMvsHRSMM_Pseudobulk.csv")
      )

      # Summary
      n_sig <- sum(res_pb$significant, na.rm = TRUE)
      n_rrmm_up <- sum(res_pb$direction == "RRMM_up", na.rm = TRUE)
      n_hrsmm_up <- sum(res_pb$direction == "HRSMM_up", na.rm = TRUE)

      message(sprintf("\nPseudobulk DEG results:"))
      message(sprintf("  Total significant (FDR < 0.05): %d", n_sig))
      message(sprintf("  Upregulated in RRMM: %d", n_rrmm_up))
      message(sprintf("  Upregulated in HRSMM: %d", n_hrsmm_up))

      # -----------------------------------------------------------------------
      # Volcano plot for pseudobulk
      # -----------------------------------------------------------------------
      message("Creating pseudobulk volcano plot...")

      # Top genes for labeling
      top_up_pb <- res_pb %>%
        filter(direction == "RRMM_up") %>%
        slice_head(n = 15)

      top_down_pb <- res_pb %>%
        filter(direction == "HRSMM_up") %>%
        slice_head(n = 15)

      top_genes_pb <- bind_rows(top_up_pb, top_down_pb)

      p_volcano_pb <- ggplot(res_pb, aes(x = log2FoldChange, y = -log10(padj), color = direction)) +
        geom_point(alpha = 0.7, size = 2) +
        scale_color_manual(
          values = c("RRMM_up" = "gold", "HRSMM_up" = "grey50", "NS" = "grey85"),
          name = "Significance"
        ) +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        geom_text_repel(
          data = top_genes_pb,
          aes(label = Gene),
          size = 3.5,
          max.overlaps = 20,
          segment.size = 0.3,
          color = "black"
        ) +
        labs(
          title = "Pseudobulk DEG: Expanding cells RRMM vs HRSMM",
          subtitle = sprintf("Patients: HRSMM = %d, RRMM = %d",
                             sum(sample_meta$Disease == "HRSMM"),
                             sum(sample_meta$Disease == "RRMM")),
          x = "log2 Fold Change (RRMM vs HRSMM)",
          y = expression(-log[10] * "(adj. p-value)")
        ) +
        theme_classic(base_size = 16) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "top"
        )

      ggsave(
        p_volcano_pb,
        filename = paste0(plot_dir, "volcano_expanding_RRMMvsHRSMM_pseudobulk.pdf"),
        width = 10,
        height = 8
      )

      message("Pseudobulk DEG analysis complete!")

    } else {
      message("Not enough patients with sufficient cells for pseudobulk analysis")
    }
  } else {
    message("Need at least 2 patients per disease group for pseudobulk analysis")
  }
} else {
  message("T_cells_expanding object not available for pseudobulk analysis")
}

# -----------------------------------------------------------------------------
# 19. Viral Specificity Matching: CDR3B vs VDJdb/SearchTable (with HLA matching)
# -----------------------------------------------------------------------------
message("\n========================================")
message("Matching CDR3B sequences to viral epitopes...")
message("========================================\n")

# Load SearchTable (VDJdb-like database of known TCR-epitope pairs)
search_table_file <- "~/TCR_project/SearchTable-2026-01-16 00_26_28.282.tsv"
hla_json_dir <- "~/TCR_project/genotype_jsons/"

# -------------------------------------------------------------------------
# First, load HLA genotypes from JSON files
# -------------------------------------------------------------------------
hla_patient_df <- NULL

if (dir.exists(hla_json_dir)) {

  json_files <- list.files(hla_json_dir, pattern = "\\.json$", full.names = TRUE)
  message(sprintf("Found %d HLA genotype files", length(json_files)))

  if (length(json_files) > 0) {

    hla_list <- list()

    for (f in json_files) {
      sample_name <- gsub("\\.genotype\\.json$", "", basename(f))
      hla_data <- fromJSON(f)

      # Extract and simplify to 2-digit resolution (e.g., A*02)
      simplify_hla <- function(allele) {
        if (is.null(allele) || is.na(allele)) return(NA_character_)
        # Extract gene and first field: A*02:01:01 -> A*02
        match <- regmatches(allele, regexpr("^[A-Z]+\\*[0-9]+", allele))
        if (length(match) > 0) return(match) else return(NA_character_)
      }

      hla_list[[sample_name]] <- tibble(
        Sample = sample_name,
        HLA_A = list(unique(na.omit(c(simplify_hla(hla_data$A[1]), simplify_hla(hla_data$A[2]))))),
        HLA_B = list(unique(na.omit(c(simplify_hla(hla_data$B[1]), simplify_hla(hla_data$B[2]))))),
        HLA_C = list(unique(na.omit(c(simplify_hla(hla_data$C[1]), simplify_hla(hla_data$C[2])))))
      )
    }

    hla_patient_df <- bind_rows(hla_list)

    # Extract patient ID from sample name (P10, P11, etc.)
    hla_patient_df <- hla_patient_df %>%
      mutate(PT_ID_raw = gsub("^(P[0-9]+)_.*", "\\1", Sample)) %>%
      group_by(PT_ID_raw) %>%
      summarise(
        HLA_A = list(unique(unlist(HLA_A))),
        HLA_B = list(unique(unlist(HLA_B))),
        HLA_C = list(unique(unlist(HLA_C))),
        .groups = "drop"
      )

    message("HLA genotypes loaded and simplified to 2-digit resolution")
    message("Example HLA alleles:")
    print(head(hla_patient_df))
  }
} else {
  message("HLA genotype directory not found - will proceed without HLA matching")
}

if (file.exists(search_table_file)) {

  # Read the database
  vdjdb <- read_tsv(search_table_file, show_col_types = FALSE) %>%
    filter(Gene == "TRB"& !`Epitope species`=="HomoSapiens") %>%
    select(CDR3, V, J, `MHC A`, `MHC class`, Epitope, `Epitope gene`, `Epitope species`, Score) %>%
    rename(
      CDR3B_db = CDR3,
      V_gene = V,
      J_gene = J,
      HLA_restriction = `MHC A`,
      MHC_class = `MHC class`,
      epitope = Epitope,
      epitope_gene = `Epitope gene`,
      epitope_species = `Epitope species`,
      confidence_score = Score
    ) %>%
    # Simplify HLA restriction to 2-digit (e.g., HLA-A*02:01 -> A*02)
    mutate(
      HLA_restriction_simple = gsub("^HLA-", "", HLA_restriction),
      HLA_restriction_simple = gsub("^([A-Z]+\\*[0-9]+).*", "\\1", HLA_restriction_simple)
    ) %>%
    distinct(CDR3B_db, epitope, epitope_species, HLA_restriction_simple, .keep_all = TRUE)

  message(sprintf("Loaded %d unique CDR3B-epitope pairs from database", nrow(vdjdb)))
  message(sprintf("Epitope species: %s", paste(unique(vdjdb$epitope_species), collapse = ", ")))

  # -------------------------------------------------------------------------
  # Get unique CDR3B sequences from our dataset (TEC cohort)
  # NOTE: Do NOT include 'category' here - it causes duplicates since the same
  # CDR3B can appear in multiple categories for the same patient
  # -------------------------------------------------------------------------
  tcr_data <- T_cells_tec_full@meta.data %>%
    filter(!is.na(TRB_cdr3_aa)) %>%
    select(PT_ID, Disease, TRB_cdr3_aa) %>%
    rename(CDR3B = TRB_cdr3_aa) %>%
    distinct()

  # Keep category information separately for category-specific analyses
  tcr_category <- T_cells_tec_full@meta.data %>%
    filter(!is.na(TRB_cdr3_aa)) %>%
    select(PT_ID, TRB_cdr3_aa, category) %>%
    rename(CDR3B = TRB_cdr3_aa) %>%
    distinct()

  message(sprintf("Unique CDR3B-patient combinations: %d", nrow(tcr_data)))
  message(sprintf("Unique CDR3B sequences: %d", n_distinct(tcr_data$CDR3B)))

  # -------------------------------------------------------------------------
  # Match CDR3B sequences (exact match)
  # -------------------------------------------------------------------------
  message("\nPerforming exact CDR3B matching...")

  tcr_with_epitope <- tcr_data %>%
    inner_join(vdjdb, by = c("CDR3B" = "CDR3B_db"))

  n_matches <- nrow(tcr_with_epitope)
  n_unique_cdr3_matched <- n_distinct(tcr_with_epitope$CDR3B)

  message(sprintf("Cells with viral specificity match (CDR3B only): %d", n_matches))
  message(sprintf("Unique CDR3B sequences matched: %d", n_unique_cdr3_matched))

  if (n_matches > 0) {

    # -------------------------------------------------------------------------
    # Check HLA concordance if HLA data available
    # -------------------------------------------------------------------------
    if (!is.null(hla_patient_df) && nrow(hla_patient_df) > 0) {

      message("\nChecking HLA concordance...")

      # Function to check if patient has the required HLA allele
      check_hla_match <- function(pt_id, hla_required, hla_df) {
        # Try to match patient ID
        pt_match <- hla_df %>%
          filter(grepl(pt_id, PT_ID_raw, ignore.case = TRUE) | PT_ID_raw == pt_id)

        if (nrow(pt_match) == 0) return(NA)  # No HLA data for this patient

        # Get the HLA gene (A, B, or C) from the restriction
        hla_gene <- gsub("^([A-Z]+)\\*.*", "\\1", hla_required)

        # Get patient's alleles for that gene
        patient_alleles <- switch(hla_gene,
          "A" = unlist(pt_match$HLA_A),
          "B" = unlist(pt_match$HLA_B),
          "C" = unlist(pt_match$HLA_C),
          NULL
        )

        if (is.null(patient_alleles) || length(patient_alleles) == 0) return(NA)

        # Check if patient has the required allele
        return(hla_required %in% patient_alleles)
      }

      # Apply HLA matching
      tcr_with_epitope <- tcr_with_epitope %>%
        rowwise() %>%
        mutate(
          HLA_concordant = check_hla_match(PT_ID, HLA_restriction_simple, hla_patient_df)
        ) %>%
        ungroup()

      # Summarize HLA concordance
      hla_summary <- tcr_with_epitope %>%
        filter(!is.na(HLA_concordant)) %>%
        count(HLA_concordant, name = "n_cells")

      message("\nHLA concordance summary:")
      print(hla_summary)

      n_hla_concordant <- sum(tcr_with_epitope$HLA_concordant == TRUE, na.rm = TRUE)
      n_hla_discordant <- sum(tcr_with_epitope$HLA_concordant == FALSE, na.rm = TRUE)
      n_hla_unknown <- sum(is.na(tcr_with_epitope$HLA_concordant))

      message(sprintf("  HLA-concordant matches: %d", n_hla_concordant))
      message(sprintf("  HLA-discordant matches: %d", n_hla_discordant))
      message(sprintf("  HLA unknown (no genotype data): %d", n_hla_unknown))

    } else {
      tcr_with_epitope$HLA_concordant <- NA
      message("No HLA data available - skipping HLA concordance check")
    }

    # -------------------------------------------------------------------------
    # Summarize matches by epitope species
    # -------------------------------------------------------------------------
    epitope_summary <- tcr_with_epitope %>%
      count(epitope_species, epitope, epitope_gene, name = "n_cells") %>%
      arrange(desc(n_cells))

    message("\nTop epitope matches (all CDR3B matches):")
    print(head(epitope_summary, 20))

    # HLA-concordant only summary
    if ("HLA_concordant" %in% colnames(tcr_with_epitope)) {
      epitope_summary_hla <- tcr_with_epitope %>%
        filter(HLA_concordant == TRUE) %>%
        count(epitope_species, epitope, epitope_gene, HLA_restriction_simple, name = "n_cells") %>%
        arrange(desc(n_cells))

      if (nrow(epitope_summary_hla) > 0) {
        message("\nTop epitope matches (HLA-concordant only):")
        print(head(epitope_summary_hla, 20))
      }
    }

    # -------------------------------------------------------------------------
    # Summarize by patient and disease
    # -------------------------------------------------------------------------
    patient_epitope_summary <- tcr_with_epitope %>%
      group_by(PT_ID, Disease, epitope_species, HLA_concordant) %>%
      summarise(
        n_cells = n(),
        n_clonotypes = n_distinct(CDR3B),
        .groups = "drop"
      ) %>%
      arrange(PT_ID, desc(n_cells))

    # -------------------------------------------------------------------------
    # Focus on expanding clones with viral specificity
    # Join category information from tcr_category
    # -------------------------------------------------------------------------
    tcr_with_epitope_category <- tcr_with_epitope %>%
      left_join(tcr_category, by = c("PT_ID", "CDR3B"))

    expanding_viral <- tcr_with_epitope_category %>%
      filter(category == "expanding")

    message(sprintf("\nExpanding cells with viral specificity: %d", nrow(expanding_viral)))

    if ("HLA_concordant" %in% colnames(expanding_viral)) {
      n_expanding_hla_concordant <- sum(expanding_viral$HLA_concordant == TRUE, na.rm = TRUE)
      message(sprintf("Expanding cells with HLA-concordant viral specificity: %d", n_expanding_hla_concordant))
    }

    if (nrow(expanding_viral) > 0) {
      expanding_epitope_summary <- expanding_viral %>%
        count(Disease, epitope_species, epitope, HLA_concordant, name = "n_cells") %>%
        arrange(Disease, desc(n_cells))

      message("\nExpanding clones - epitope matches by disease:")
      print(expanding_epitope_summary)
    }

    # -------------------------------------------------------------------------
    # Save results
    # -------------------------------------------------------------------------
    write_csv(tcr_with_epitope, paste0(output_dir, "CDR3B_viral_specificity_matches.csv"))
    write_csv(epitope_summary, paste0(output_dir, "CDR3B_epitope_summary.csv"))
    write_csv(patient_epitope_summary, paste0(output_dir, "CDR3B_patient_epitope_summary.csv"))

    # Save HLA-concordant only
    if ("HLA_concordant" %in% colnames(tcr_with_epitope)) {
      tcr_hla_concordant <- tcr_with_epitope %>% filter(HLA_concordant == TRUE)
      if (nrow(tcr_hla_concordant) > 0) {
        write_csv(tcr_hla_concordant, paste0(output_dir, "CDR3B_viral_specificity_HLA_concordant.csv"))
      }
    }

    # -------------------------------------------------------------------------
    # Visualization: Comprehensive CDR3B and VDJdb matching analysis
    # -------------------------------------------------------------------------
    message("Creating viral specificity plots...")

    # =========================================================================
    # PLOT 1: Total CDR3B per patient (all unique CDR3B sequences)
    # =========================================================================
    total_cdr3b_per_patient <- tcr_data %>%
      group_by(PT_ID, Disease) %>%
      summarise(
        n_unique_cdr3b = n_distinct(CDR3B),
        .groups = "drop"
      ) %>%
      arrange(desc(n_unique_cdr3b))

    # Calculate median and range
    median_cdr3b <- median(total_cdr3b_per_patient$n_unique_cdr3b)
    range_cdr3b <- range(total_cdr3b_per_patient$n_unique_cdr3b)
    stats_label_cdr3b <- sprintf("Median: %.0f\nRange: %.0f - %.0f",
                                  median_cdr3b, range_cdr3b[1], range_cdr3b[2])

    p_total_cdr3b <- ggplot(total_cdr3b_per_patient,
                            aes(x = reorder(PT_ID, n_unique_cdr3b), y = n_unique_cdr3b, fill = Disease)) +
      geom_col(color = "black", linewidth = 0.3) +
      geom_text(aes(label = n_unique_cdr3b), hjust = -0.1, size = 3) +
      coord_flip() +
      labs(
        x = "Patient",
        y = "Number of Unique CDR3B Sequences",
        title = "Total CDR3B sequences per patient",
        fill = "Disease"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      annotate("text", x = Inf, y = Inf, label = stats_label_cdr3b,
               hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )

    ggsave(
      p_total_cdr3b,
      filename = paste0(plot_dir, "total_CDR3B_per_patient.pdf"),
      width = 10,
      height = 8
    )

    # =========================================================================
    # PLOT 2: Total clonotypes per patient matching VDJdb (HLA concordant)
    # =========================================================================
    # Filter for HLA-concordant matches only
    tcr_vdjdb_matches <- tcr_with_epitope %>%
      filter(HLA_concordant == TRUE)

    vdjdb_per_patient <- tcr_vdjdb_matches %>%
      group_by(PT_ID, Disease) %>%
      summarise(
        n_matching_clonotypes = n_distinct(CDR3B),
        .groups = "drop"
      ) %>%
      arrange(desc(n_matching_clonotypes))

    # Calculate median and range for VDJdb matches
    median_vdjdb <- median(vdjdb_per_patient$n_matching_clonotypes)
    range_vdjdb <- range(vdjdb_per_patient$n_matching_clonotypes)
    stats_label_vdjdb <- sprintf("Median: %.0f\nRange: %.0f - %.0f",
                                  median_vdjdb, range_vdjdb[1], range_vdjdb[2])

    p_vdjdb_per_patient <- ggplot(vdjdb_per_patient,
                                   aes(x = reorder(PT_ID, n_matching_clonotypes),
                                       y = n_matching_clonotypes, fill = Disease)) +
      geom_col(color = "black", linewidth = 0.3) +
      geom_text(aes(label = n_matching_clonotypes), hjust = -0.1, size = 3) +
      coord_flip() +
      labs(
        x = "Patient",
        y = "Number of Unique Clonotypes",
        title = "Clonotypes matching VDJdb per patient\n(HLA concordant)",
        fill = "Disease"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      annotate("text", x = Inf, y = Inf, label = stats_label_vdjdb,
               hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )

    ggsave(
      p_vdjdb_per_patient,
      filename = paste0(plot_dir, "VDJdb_matching_clonotypes_per_patient.pdf"),
      width = 10,
      height = 8
    )

    # =========================================================================
    # PLOT 3: Number of unique clonotypes matching different species
    # =========================================================================
    clonotypes_by_species <- tcr_vdjdb_matches %>%
      group_by(epitope_species) %>%
      summarise(
        n_unique_clonotypes = n_distinct(CDR3B),
        .groups = "drop"
      ) %>%
      arrange(desc(n_unique_clonotypes)) %>%
      mutate(epitope_species = fct_reorder(epitope_species, n_unique_clonotypes))

    # Calculate median and range
    median_species <- median(clonotypes_by_species$n_unique_clonotypes)
    range_species <- range(clonotypes_by_species$n_unique_clonotypes)
    stats_label_species <- sprintf("Median: %.0f\nRange: %.0f - %.0f",
                                    median_species, range_species[1], range_species[2])

    p_clonotypes_by_species <- ggplot(clonotypes_by_species,
                                       aes(x = epitope_species, y = n_unique_clonotypes,
                                           fill = epitope_species)) +
      geom_col(color = "black", linewidth = 0.3) +
      geom_text(aes(label = n_unique_clonotypes), hjust = -0.1, size = 3.5) +
      coord_flip() +
      labs(
        x = "Epitope Species",
        y = "Number of Unique Clonotypes",
        title = "Unique clonotypes matching different species"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      annotate("text", x = Inf, y = Inf, label = stats_label_species,
               hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
      theme_classic(base_size = 14) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )

    ggsave(
      p_clonotypes_by_species,
      filename = paste0(plot_dir, "unique_clonotypes_by_species.pdf"),
      width = 10,
      height = 6
    )

    # =========================================================================
    # PLOT 4: Unique clonotypes matching species per category (clonal dynamics)
    # =========================================================================
    # Join category information from tcr_category
    tcr_vdjdb_with_category <- tcr_vdjdb_matches %>%
      left_join(tcr_category, by = c("PT_ID", "CDR3B"))

    if ("category" %in% colnames(tcr_vdjdb_with_category)) {

      clonotypes_species_category <- tcr_vdjdb_with_category %>%
        filter(!is.na(category)) %>%
        group_by(category, epitope_species) %>%
        summarise(
          n_unique_clonotypes = n_distinct(CDR3B),
          .groups = "drop"
        ) %>%
        mutate(category = factor(category, levels = c("expanding", "stable", "contracting")))

      # Create summary for totals per category
      category_totals <- clonotypes_species_category %>%
        group_by(category) %>%
        summarise(total = sum(n_unique_clonotypes), .groups = "drop")

      # Stacked bar with absolute numbers
      p_species_category <- ggplot(clonotypes_species_category,
                                    aes(x = category, y = n_unique_clonotypes, fill = epitope_species)) +
        geom_col(color = "black", linewidth = 0.3, position = "stack") +
        geom_text(aes(label = ifelse(n_unique_clonotypes > 0, n_unique_clonotypes, "")),
                  position = position_stack(vjust = 0.5), size = 3, color = "white", fontface = "bold") +
        geom_text(data = category_totals, aes(x = category, y = total, label = total, fill = NULL),
                  vjust = -0.3, size = 4, fontface = "bold") +
        labs(
          x = "Clonal Dynamics Category",
          y = "Number of Unique Clonotypes",
          fill = "Epitope Species",
          title = "Unique clonotypes by species per clonal category"
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme_classic(base_size = 14) +
        theme(
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold")
        )

      ggsave(
        p_species_category,
        filename = paste0(plot_dir, "unique_clonotypes_species_by_category.pdf"),
        width = 10,
        height = 7
      )

      # Also create faceted version for clearer per-species view
      p_species_category_facet <- ggplot(clonotypes_species_category,
                                          aes(x = category, y = n_unique_clonotypes, fill = category)) +
        geom_col(color = "black", linewidth = 0.3) +
        geom_text(aes(label = n_unique_clonotypes), vjust = -0.3, size = 3.5) +
        facet_wrap(~ epitope_species, scales = "free_y", ncol = 3) +
        labs(
          x = "Clonal Dynamics Category",
          y = "Number of Unique Clonotypes",
          title = "Unique clonotypes per species by clonal category"
        ) +
        scale_fill_manual(values = c("expanding" = "#D81B60", "stable" = "grey60", "contracting" = "#00A6D6")) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_classic(base_size = 12) +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold", size = 10)
        )

      ggsave(
        p_species_category_facet,
        filename = paste0(plot_dir, "unique_clonotypes_species_by_category_faceted.pdf"),
        width = 12,
        height = 10
      )
    }

    # =========================================================================
    # PLOT 5: Clonotypes matching VDJdb by Disease (with absolute numbers)
    # =========================================================================
    clonotypes_by_disease <- tcr_vdjdb_matches %>%
      group_by(Disease, epitope_species) %>%
      summarise(
        n_unique_clonotypes = n_distinct(CDR3B),
        .groups = "drop"
      )

    disease_totals <- clonotypes_by_disease %>%
      group_by(Disease) %>%
      summarise(total = sum(n_unique_clonotypes), .groups = "drop")

    p_clonotypes_disease <- ggplot(clonotypes_by_disease,
                                    aes(x = Disease, y = n_unique_clonotypes, fill = epitope_species)) +
      geom_col(color = "black", linewidth = 0.3, position = "stack") +
      geom_text(aes(label = ifelse(n_unique_clonotypes > 0, n_unique_clonotypes, "")),
                position = position_stack(vjust = 0.5), size = 3, color = "white", fontface = "bold") +
      geom_text(data = disease_totals, aes(x = Disease, y = total, label = total, fill = NULL),
                vjust = -0.3, size = 4, fontface = "bold") +
      labs(
        x = "Disease",
        y = "Number of Unique Clonotypes",
        fill = "Epitope Species",
        title = "Unique clonotypes matching VDJdb by disease"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      theme_classic(base_size = 14) +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )

    ggsave(
      p_clonotypes_disease,
      filename = paste0(plot_dir, "unique_clonotypes_by_disease.pdf"),
      width = 10,
      height = 7
    )

    # =========================================================================
    # Summary table: per patient breakdown
    # =========================================================================
    patient_species_summary <- tcr_vdjdb_matches %>%
      group_by(PT_ID, Disease, epitope_species) %>%
      summarise(
        n_unique_clonotypes = n_distinct(CDR3B),
        .groups = "drop"
      ) %>%
      pivot_wider(names_from = epitope_species, values_from = n_unique_clonotypes, values_fill = 0)

    write_csv(patient_species_summary, paste0(output_dir, "CDR3B_patient_species_clonotype_summary.csv"))

    message("Viral specificity analysis complete!")

  } else {
    message("No CDR3B matches found in the database")
  }

} else {
  message(sprintf("SearchTable file not found: %s", search_table_file))
}









# -------------------------------------------------------------------------
#20. pTRT Signature Analysis
# -------------------------------------------------------------------------
message("\n========================================")
message("pTRT Signature Analysis")
message("========================================\n")

# Define pTRT gene signature
pTRT <- c("GNLY", "ZNF683", "GZMH", "FGFBP2", "GZMB", "NKG7", "CCL5", "HOPX",
           "KLRD1", "EFHD2", "CD8A", "CTSW", "CST7", "ITGB1", "BHLHE40")

# Score cells using AddModuleScore
message("Scoring cells with pTRT signature...")
T_cells_tec_full <- AddModuleScore(
  T_cells_tec_full,
  features = list(pTRT = pTRT),
  name = "pTRT_score"
)

# Rename the score column (AddModuleScore appends "1" to the name)
T_cells_tec_full$pTRT_score <- T_cells_tec_full$pTRT_score1
T_cells_tec_full$pTRT_score1 <- NULL

message(sprintf("pTRT signature scored for %d cells", ncol(T_cells_tec_full)))

# -------------------------------------------------------------------------
# 20a. VlnPlot per final_annotation
# -------------------------------------------------------------------------
message("\nCreating VlnPlot of pTRT score by final annotation...")

T_cells_tec_full@meta.data <- T_cells_tec_full@meta.data %>%
  mutate(
    final_annotation = fct_reorder(
      final_annotation,
      pTRT_score,
      .fun = median,
      .desc = TRUE   # highest → lowest
    )
  )
p_pTRT_vln <- VlnPlot(
  T_cells_tec_full,
  features = "pTRT_score",
  group.by = "final_annotation",
  pt.size = 0
) &
  labs(
    title = "pTRT Signature Score by Cell Type",
    x = "Cell Annotation",
    y = "pTRT Score"
  ) &
  theme_classic(base_size = 12) &
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) &
  NoLegend()

ggsave(p_pTRT_vln, filename = paste0(plot_dir, "pTRT_score_by_annotation.pdf"),
       width = 10, height = 6)
message("Saved: pTRT_score_by_annotation.pdf")

# -------------------------------------------------------------------------
# 20b. Mean pTRT score per sample by Timepoint and category (Post vs Pre)
# -------------------------------------------------------------------------
message("\nCalculating mean pTRT scores per sample...")

# Calculate mean pTRT score per sample, Timepoint, and category
pTRT_means <- T_cells_tec_full@meta.data %>%
  filter(final_annotation%in%c("CD8+GZMB+TEM",
                               "CD8+GZMK+TEM",
                               "CD8+KIR+TEM",
                               "CD8+IFN+",
                               "CD8+TNFRSF9+",
                               "CD4+CTL")&!category=="contracting")%>%
  filter(!is.na(category), Timepoint %in% c("Pre", "Post")) %>%
  group_by(PT_ID, Disease, Timepoint, category) %>%
  summarise(
    mean_pTRT = mean(pTRT_score, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )%>%filter(n_cells>10)

# Create boxplot comparing Pre vs Post for each category
p_pTRT_timepoint <- ggplot(
  pTRT_means,
  aes(x = category, y = mean_pTRT, fill = category)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(color = Disease), position = position_jitter(width = 0.1), size = 2) +
  facet_wrap(~Timepoint, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     comparisons = list(c("expanding", "stable"))) +
  scale_fill_manual(values = c("Pre" = "#4DBEEE", "Post" = "#D95319")) +
  scale_color_manual(values = c("HRSMM" = "#77AC30", "RRMM" = "#7E2F8E")) +
  labs(
    title = "pTRT Signature: Pre vs Post Treatment",
    y = "Mean pTRT Score",
    x = "Timepoint"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(p_pTRT_timepoint, filename = paste0(plot_dir, "pTRT_score_pre_vs_post_by_category.pdf"),
       width = 10, height = 5)
message("Saved: pTRT_score_pre_vs_post_by_category.pdf")

# Perform statistical tests for each category
message("\nStatistical tests (Wilcoxon) for Pre vs Post by category:")
pTRT_stats <- pTRT_means %>%
  group_by(category) %>%
  summarise(
    n_pre = sum(Timepoint == "Pre"),
    n_post = sum(Timepoint == "Post"),
    median_pre = median(mean_pTRT[Timepoint == "Pre"], na.rm = TRUE),
    median_post = median(mean_pTRT[Timepoint == "Post"], na.rm = TRUE),
    p_value = tryCatch(
      wilcox.test(mean_pTRT[Timepoint == "Pre"],
                  mean_pTRT[Timepoint == "Post"])$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    significant = fdr < 0.05
  )

print(pTRT_stats)
write_csv(pTRT_stats, paste0(output_dir, "pTRT_score_stats_by_category.csv"))

# -------------------------------------------------------------------------
# 20c. Repeat analysis EXCLUDING viral-specific CDR3B sequences
# -------------------------------------------------------------------------
message("\n----------------------------------------")
message("Repeating pTRT analysis excluding viral-specific CDR3B...")
message("----------------------------------------\n")

# Load viral specificity matches
viral_cdr3b_file <- paste0(output_dir, "CDR3B_viral_specificity_matches.csv")

if (file.exists(viral_cdr3b_file)) {
  viral_cdr3b <- read_csv(viral_cdr3b_file, show_col_types = FALSE)

  # Get unique viral-specific CDR3B sequences
  viral_cdr3b_list <- unique(viral_cdr3b$CDR3B)
  message(sprintf("Loaded %d unique viral-specific CDR3B sequences", length(viral_cdr3b_list)))

  # Identify cells with viral-specific CDR3B
  cells_viral <- T_cells_tec_full@meta.data %>%
    filter(TRB_cdr3_aa %in% viral_cdr3b_list) %>%
    rownames()

  message(sprintf("Cells with viral-specific CDR3B: %d (%.1f%%)",
                  length(cells_viral),
                  100 * length(cells_viral) / ncol(T_cells_tec_full)))

  # Subset to exclude viral-specific cells
  cells_nonviral <- setdiff(colnames(T_cells_tec_full), cells_viral)
  T_cells_nonviral <- subset(T_cells_tec_full, cells = cells_nonviral)

  message(sprintf("Cells remaining after exclusion: %d", ncol(T_cells_nonviral)))

  # -------------------------------------------------------------------------
  # 20c.1 VlnPlot per final_annotation (non-viral)
  # -------------------------------------------------------------------------
  message("\nCreating VlnPlot (excluding viral-specific CDR3B)...")

  p_pTRT_vln_nonviral <- VlnPlot(
    T_cells_nonviral,
    features = "pTRT_score",
    group.by = "final_annotation",
    pt.size = 0
  ) &
    labs(
      title = "pTRT Signature Score by Cell Type\n(Excluding Viral-Specific CDR3B)",
      x = "Cell Annotation",
      y = "pTRT Score"
    ) &
    theme_classic(base_size = 12) &
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) &
    NoLegend()

  ggsave(p_pTRT_vln_nonviral,
         filename = paste0(plot_dir, "pTRT_score_by_annotation_nonviral.pdf"),
         width = 10, height = 6)
  message("Saved: pTRT_score_by_annotation_nonviral.pdf")

  # -------------------------------------------------------------------------
  # 20c.2 Mean pTRT score per sample by Timepoint and category (non-viral)
  # -------------------------------------------------------------------------
  message("\nCalculating mean pTRT scores (excluding viral-specific CDR3B)...")

  pTRT_means_nonviral <- T_cells_nonviral@meta.data %>%
    filter(final_annotation%in%c("CD8+GZMB+TEM",
                                 "CD8+GZMK+TEM",
                                 "CD8+KIR+TEM",
                                 "CD8+IFN+",
                                 "CD8+TNFRSF9+",
                                 "CD4+CTL")&!category=="contracting")%>%
    
    filter(!is.na(category), Timepoint %in% c("Pre", "Post")) %>%
    group_by(PT_ID, Disease, Timepoint, category) %>%
    summarise(
      mean_pTRT = mean(pTRT_score, na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    )%>%filter(n_cells>10)

  # Create boxplot comparing Pre vs Post for each category
  p_pTRT_timepoint_nonviral <- ggplot(
    pTRT_means_nonviral,
    aes(x = category, y = mean_pTRT, fill = category)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_point(aes(color = Disease), position = position_jitter(width = 0.1), size = 2) +
    facet_wrap(~Timepoint,  scales = "free_y") +
    stat_compare_means(method = "wilcox.test", label = "p.format",
                       comparisons = list(c("expanding", "stable"))) +
    scale_fill_manual(values = c("Pre" = "#4DBEEE", "Post" = "#D95319")) +
    scale_color_manual(values = c("HRSMM" = "#77AC30", "RRMM" = "#7E2F8E")) +
    labs(
      title = "pTRT Signature: Pre vs Post Treatment\n(Excluding Viral-Specific CDR3B)",
      y = "Mean pTRT Score",
      x = "Timepoint"
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  ggsave(p_pTRT_timepoint_nonviral,
         filename = paste0(plot_dir, "pTRT_score_pre_vs_post_by_category_nonviral.pdf"),
         width = 10, height = 5)
  message("Saved: pTRT_score_pre_vs_post_by_category_nonviral.pdf")

  # Perform statistical tests for each category (non-viral)
  message("\nStatistical tests (Wilcoxon) for Pre vs Post by category (non-viral):")
  pTRT_stats_nonviral <- pTRT_means_nonviral %>%
    group_by(category) %>%
    summarise(
      n_pre = sum(Timepoint == "Pre"),
      n_post = sum(Timepoint == "Post"),
      median_pre = median(mean_pTRT[Timepoint == "Pre"], na.rm = TRUE),
      median_post = median(mean_pTRT[Timepoint == "Post"], na.rm = TRUE),
      p_value = tryCatch(
        wilcox.test(mean_pTRT[Timepoint == "Pre"],
                    mean_pTRT[Timepoint == "Post"])$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      fdr = p.adjust(p_value, method = "BH"),
      significant = fdr < 0.05
    )

  print(pTRT_stats_nonviral)
  write_csv(pTRT_stats_nonviral,
            paste0(output_dir, "pTRT_score_stats_by_category_nonviral.csv"))

  # Compare viral vs non-viral analysis
  message("\n----------------------------------------")
  message("Comparison: All cells vs Non-viral cells")
  message("----------------------------------------")
  comparison_df <- pTRT_stats %>%
    select(category, p_value_all = p_value, fdr_all = fdr) %>%
    left_join(
      pTRT_stats_nonviral %>%
        select(category, p_value_nonviral = p_value, fdr_nonviral = fdr),
      by = "category"
    )
  print(comparison_df)
  write_csv(comparison_df, paste0(output_dir, "pTRT_score_viral_vs_nonviral_comparison.csv"))

} else {
  message("WARNING: Viral specificity file not found. Skipping non-viral analysis.")
  message(sprintf("Expected file: %s", viral_cdr3b_file))
}

message("\npTRT signature analysis complete!")

