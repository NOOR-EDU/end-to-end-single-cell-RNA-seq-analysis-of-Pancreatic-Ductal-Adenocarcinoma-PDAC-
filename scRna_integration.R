# Core packages (already installed in Part 2)
# Seurat, ggplot2, dplyr, patchwork

# Install Bioconductor packages for advanced integration
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "batchelor",  # FastMNN integration (used by Seurat's IntegrateLayers)
  "scran"       # Additional normalization methods
))
# Install Seurat 5 (the star of this tutorial!)
install.packages("Seurat")

# Install SeuratObject (Seurat's data structure package)
install.packages("SeuratObject")
# Install Harmony for integration
install.packages("harmony")

# Install additional visualization and utility packages
install.packages(c(
  "remotes",         # Install R packages stored in GitHub
  "ggrepel",         # Non-overlapping text labels
  "RColorBrewer",    # Color palettes
  "viridis"          # Perceptually uniform color scales
))
# Install visualization packages
install.packages(c(
  "ggplot2",        # For custom plots when needed
  "patchwork",      # For combining Seurat plots
  "dplyr"           # Data manipulation
))
# Wrapper functions for integration methods
remotes::install_github('satijalab/seurat-wrappers')

# Install future for parallel processing (used by RPCA integration)
install.packages("future")

# Install FNN for k-nearest neighbor calculations (used in quality metrics)
install.packages("FNN")

# Install cluster for clustering quality metrics (silhouette scores)
install.packages("cluster")

# Install reshape2 for data reshaping (used in cluster quality assessment)
install.packages("reshape2")



# STEP 1: Load required libraries
#-----------------------------------------------

# Core single-cell analysis
library(Seurat)
library(SeuratObject)

# Integration methods
library(harmony)
library(batchelor)
library(SeuratWrappers)

# Visualization and data manipulation
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(viridis)

# Quality metrics and utilities
library(FNN)              # K-nearest neighbor calculations for mixing metrics
library(cluster)          # Silhouette scores for clustering quality
library(reshape2)         # Data reshaping for visualization

# Parallel processing
library(future)
options(future.globals.maxSize = 20 * 1024^3)  # Increase to 20GB for large datasets


# Set working directory (adjust to your path)
setwd("~/GSE212966_scRNA/integration_analysis")

# Create output directories
dir.create("plots", showWarnings = FALSE)
dir.create("plots/integration_comparison", showWarnings = FALSE)
dir.create("plots/clustering", showWarnings = FALSE)
dir.create("integrated_data", showWarnings = FALSE)
dir.create("metadata", showWarnings = FALSE)

# Set random seed for reproducibility
set.seed(42)

# Configure plotting defaults
theme_set(theme_classic(base_size = 12))

cat("Environment setup complete\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
# Seurat version: 5.3.1


#--------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
# file correcton: 
list.files("~/galaxy_inputs")
files <- list.files(
  "~/galaxy_inputs",
  pattern = "cell_metadata_cancerous.*\\.csv$",
  full.names = TRUE
)

for (f in files) {
  df <- read.csv(f)
  cat("\n", basename(f), "\n")
  print(unique(df$patient_id))
}
mapping <- c(
  "cancerous_1" = "Pancreatic_cancer_1",
  "cancerous_2" = "Pancreatic_cancer_2",
  "cancerous_3" = "Pancreatic_cancer_3",
  "cancerous_4" = "Pancreatic_cancer_4",
  "cancerous_5" = "Pancreatic_cancer_5",
  "cancerous_6" = "Pancreatic_cancer_6"
)


dir.create("~/corrected_metadata", showWarnings = FALSE)
for (f in files) {
  
  sample_num <- sub(".*cancerous_([0-9]+)\\.csv", "\\1", basename(f))
  
  df <- read.csv(f, stringsAsFactors = FALSE)
  
  df$patient_id <- paste0("Pancreatic_cancer_", sample_num)
  
  out_file <- file.path(
    "~/corrected_metadata",
    basename(f)
  )
  
  write.csv(df, out_file, row.names = FALSE)
  
  cat("Created:", out_file, "\n")
}

df <- read.csv("~/corrected_metadata/1_cell_metadata_cancerous_6.csv")
unique(df$patient_id)

files <- list.files(
  "~/galaxy_inputs",
  pattern = "cell_metadata_noncancerous.*\\.csv$",
  full.names = TRUE
)

for (f in files) {
  
  sample_num <- sub(".*noncancerous_([0-9]+)\\.csv", "\\1", basename(f))
  
  df <- read.csv(f, stringsAsFactors = FALSE)
  
  df$patient_id <- paste0("Adjacent_noncancerous_", sample_num)
  
  out_file <- file.path(
    "~/corrected_metadata",
    basename(f)
  )
  
  write.csv(df, out_file, row.names = FALSE)
  
  cat("Created:", basename(f),
      "-> Adjacent_noncancerous_", sample_num, "\n", sep = "")
}
files2 <- list.files(
  "~/corrected_metadata",
  pattern = "noncancerous.*\\.csv$",
  full.names = TRUE
)

for (f in files2) {
  cat("\n", basename(f), "\n")
  print(unique(read.csv(f)$patient_id))
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define sample metadata
sample_metadata <- data.frame(
  sample_id = c(
    "Adjacent_noncancerous_1", "Adjacent_noncancerous_2", "Adjacent_noncancerous_3", "Adjacent_noncancerous_4", "Adjacent_noncancerous_5", "Adjacent_noncancerous_6",
    "Pancreatic_cancer_1", "Pancreatic_cancer_2", "Pancreatic_cancer_3", "Pancreatic_cancer_4", "Pancreatic_cancer_5", "Pancreatic_cancer_6"
  ),
  condition = c(
    rep("noncancerous", 6),
    rep("cancerous", 6)
  ),
  patient_id = c(
    "Adjacent_noncancerous_1", "Adjacent_noncancerous_2", "Adjacent_noncancerous_3", "Adjacent_noncancerous_4", "Adjacent_noncancerous_5", "Adjacent_noncancerous_6",
    "Pancreatic_cancer_1", "Pancreatic_cancer_2", "Pancreatic_cancer_3", "Pancreatic_cancer_4", "Pancreatic_cancer_5", "Pancreatic_cancer_6"
  ),
  stringsAsFactors = FALSE
)


dir.create("~/GSE212966_scRNA/QC_analysis/Filtered_data", recursive = TRUE, showWarnings = FALSE)
source_dir <- "~/galaxy_inputs"
dest_dir <- "~/GSE212966_scRNA/QC_analysis/Filtered_data"
list.files("~/galaxy_inputs")

# 1. Define the directory path
dir_path <- "~/galaxy_inputs"

# 2. List all the .rds files in that directory
file_list <- list.files(dir_path, pattern = "\\.rds$", full.names = TRUE)

for (f in file_list) {
  
  clean_name <- sub("^[0-9]+_", "", basename(f))
  
  # Remove duplicated .rds extension (.rds.rds -> .rds)
  clean_name <- sub("\\.rds\\.rds$", ".rds", clean_name, ignore.case = TRUE)
  
  dest_file <- file.path(dest_dir, clean_name)
  
  file.copy(f, dest_file, overwrite = TRUE)
  
  cat("Moved:", basename(f), "→", clean_name, "\n")
}
qc_data_path <- "~/GSE212966_scRNA/QC_analysis/Filtered_data"

#trouble shooting: 
qc_dir <- "~/GSE212966_scRNA/QC_analysis/Filtered_data"

file_list <- list.files(qc_dir, full.names = TRUE)

for (f in file_list) {
  
  fname <- basename(f)
  
  # Extract the sample name
  sample_name <- sub(
    "^cellranger_output_(.+?)_QC_analysis_filtered_data_.*$",
    "\\1",
    fname
  )
  
  # New filename
  new_name <- paste0(sample_name, "_qc_filtered.rds")
  
  # Rename in place
  file.rename(f, file.path(qc_dir, new_name))
  
  cat("Renamed:", fname, "→", new_name, "\n")
}


all(file.exists(file.path(qc_data_path,
                          paste0(sample_metadata$sample_id, "_qc_filtered.rds"))))

list.files("~/GSE212966_scRNA/QC_analysis/Filtered_data")
print(sample_metadata$sample_id)
# Load QC-filtered Seurat objects
qc_data_path <- "~/GSE212966_scRNA/QC_analysis/Filtered_data"

seurat_list <- lapply(sample_metadata$sample_id, function(sample_id) {
  obj <- readRDS(file.path(qc_data_path, paste0(sample_id, "_qc_filtered.rds")))
  
  # Ensure metadata is consistent
  obj$sample_id <- sample_id
  obj$condition <- sample_metadata$condition[sample_metadata$sample_id == sample_id]
  obj$patient_id <- sample_metadata$patient_id[sample_metadata$sample_id == sample_id]
  
  return(obj)
})

names(seurat_list) <- sample_metadata$sample_id

# Ensure all samples have the same genes (critical for FastMNN integration)
all_genes <- lapply(seurat_list, rownames)
common_genes <- Reduce(intersect, all_genes)

# Subset all samples to common genes
seurat_list <- lapply(seurat_list, function(obj) {
  obj[common_genes, ]
})


# Report dimensions
cat("\nLoaded", length(seurat_list), "QC-filtered samples:\n")
for (sample_id in names(seurat_list)) {
  cat(sprintf("  %s: %d cells × %d genes\n", 
              sample_id, 
              ncol(seurat_list[[sample_id]]), 
              nrow(seurat_list[[sample_id]])))
}

cat("\nTotal cells across all samples:", 
    sum(sapply(seurat_list, ncol)), "\n")



# STEP 3: Naive merge without integration
#-----------------------------------------------

# Merge all samples
merged_naive <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE212966_Naive_Merge"
)

cat("Merged dataset:", ncol(merged_naive), "cells ×", nrow(merged_naive), "genes\n")
# Merged dataset: 37131 cells × 13874 genes

# Standard Seurat workflow
merged_naive <- NormalizeData(merged_naive, verbose = FALSE)
merged_naive <- FindVariableFeatures(merged_naive, nfeatures = 2000, verbose = FALSE)
merged_naive <- ScaleData(merged_naive, verbose = FALSE)
merged_naive <- RunPCA(merged_naive, npcs = 50, verbose = FALSE)
merged_naive <- RunUMAP(merged_naive, dims = 1:30, reduction = "pca", verbose = FALSE)
merged_naive <- FindNeighbors(merged_naive, dims = 1:30, verbose = FALSE)
merged_naive <- FindClusters(merged_naive, resolution = 0.6, verbose = FALSE)

cat("Clustering complete:", length(unique(merged_naive$seurat_clusters)), "clusters identified\n")
# Clustering complete: 28 clusters identified


# STEP 4: Visualize batch effects in naive merge
#-----------------------------------------------

# Define color palettes for visualization
# 8 samples need 8 distinct, visible colors
sample_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", "#00CED1", # noncancerous 1-6: Red, Blue, Green, Purple
  "#FF7F00", "#A65628", "#F781BF", "#F0E442", "#009E73", "#88CCEE"   # cancerous 1-6: Orange, Brown, Pink, Gray
)
names(sample_colors) <- sample_metadata$sample_id

# 2 conditions need 2 distinct colors
condition_colors <- c(
  "noncancerous" = "#2E86AB",       # Blue
  "cancerous" = "#F18F01"  # Orange
)
# UMAP colored by sample - shows batch effects
p1_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "sample_id", 
                    pt.size = 0.05, cols = sample_colors) +
  ggtitle("Naive Merge: Colored by Sample") +
  theme(legend.position = "right", legend.text = element_text(size = 8))
p1_naive

# UMAP colored by condition
p2_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "condition",
                    pt.size = 0.05, cols = condition_colors) +
  ggtitle("Naive Merge: Colored by Condition")
p2_naive

# UMAP colored by clusters
p3_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "seurat_clusters",
                    pt.size = 0.05, label = TRUE, label.size = 5) +
  ggtitle("Naive Merge: Clusters") +
  NoLegend()
p3_naive
# Split by sample to see separation
p4_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "condition",
                    split.by = "sample_id", pt.size = 0.05, ncol = 6, 
                    cols = condition_colors) +
  ggtitle("Naive Merge: Split by Sample") +
  theme(strip.text = element_text(size = 9, face = "bold"))
p4_naive
getwd()
# Combine plots
combined_naive <- (p1_naive | p2_naive) / (p3_naive | p4_naive)
ggsave("plots/01_naive_merge_batch_effects.png", combined_naive, 
       width = 16, height = 12, dpi = 300)


# STEP 5: Quantify batch effects with mixing metrics
#-----------------------------------------------

# Calculate per-cluster sample composition
cluster_composition <- table(merged_naive$seurat_clusters, merged_naive$sample_id)
cluster_composition_pct <- prop.table(cluster_composition, margin = 1) * 100

# Find clusters dominated by single samples (>50% from one sample = batch-driven)
dominant_sample_clusters <- apply(cluster_composition_pct, 1, max) > 50
n_batch_clusters <- sum(dominant_sample_clusters)

cat("\nBatch effect assessment:\n")
cat("  Clusters dominated by single sample (>50%):", n_batch_clusters, 
    "/", nrow(cluster_composition), "\n")
# Clusters dominated by single sample (>50%): 8 / 28

if (n_batch_clusters > 0) {
  cat("⚠ Naive merge shows significant batch effects\n")
}
# → Integration is necessary




# STEP 6: Prepare data for integration
#-----------------------------------------------

# Merge all samples with split layers (required for Seurat 5 integration)
merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE212966_Integration"
)

# Add all metadata
# Extract sample_id from cell names (added via add.cell.ids)
# Cell names format: "SampleID_BARCODE-1"
cell_names <- colnames(merged_seurat)
merged_seurat$sample_id <- gsub("_[ACGT].*$", "", cell_names)  # Remove barcode, keep sample_id

# Add other metadata by matching sample_id
merged_seurat$condition <- sample_metadata$condition[match(merged_seurat$sample_id, sample_metadata$sample_id)]
merged_seurat$patient_id <- sample_metadata$patient_id[match(merged_seurat$sample_id, sample_metadata$sample_id)]

# Check if layers are already split
current_layers <- Layers(merged_seurat[["RNA"]])
if (length(current_layers) > 1) {
  cat("Layers already split by merge operation (", length(current_layers), " layers)\n", sep = "")
} else {
  # Split layers by sample if not already split
  cat("Splitting layers by sample\n")
  merged_seurat[["RNA"]] <- split(merged_seurat[["RNA"]], f = merged_seurat$sample_id)
}

# Normalize and find variable features on split layers
merged_seurat <- NormalizeData(merged_seurat, verbose = FALSE)
merged_seurat <- FindVariableFeatures(merged_seurat, nfeatures = 2000, verbose = FALSE)

# Scale data and run PCA on split layers
merged_seurat <- ScaleData(merged_seurat, verbose = FALSE)
merged_seurat <- RunPCA(merged_seurat, npcs = 50, verbose = FALSE)

cat("Data prepared for integration\n")
cat("Layers:", length(Layers(merged_seurat, search = "data")), "samples\n")
# Layers: 12 samples

# STEP 7: CCA Integration
#-----------------------------------------------

# Integrate using CCA
integrated_cca <- IntegrateLayers(
  object = merged_seurat,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  dims = 1:30,
  verbose = FALSE
)

# Standard downstream workflow
integrated_cca <- FindNeighbors(integrated_cca, reduction = "integrated.cca", dims = 1:30, verbose = FALSE)
integrated_cca <- FindClusters(integrated_cca, resolution = 0.6, verbose = FALSE)
integrated_cca <- RunUMAP(integrated_cca, reduction = "integrated.cca", dims = 1:30, 
                          reduction.name = "umap.cca", verbose = FALSE)

cat("CCA integration complete:", length(unique(integrated_cca$seurat_clusters)), "clusters\n")
# CCA integration complete: 16 clusters

# STEP 8: RPCA Integration
#-----------------------------------------------

# Integrate using RPCA
integrated_rpca <- IntegrateLayers(
  object = merged_seurat,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  dims = 1:30,
  verbose = FALSE
)

# Downstream workflow
integrated_rpca <- FindNeighbors(integrated_rpca, reduction = "integrated.rpca", dims = 1:30, verbose = FALSE)
integrated_rpca <- FindClusters(integrated_rpca, resolution = 0.6, verbose = FALSE)
integrated_rpca <- RunUMAP(integrated_rpca, reduction = "integrated.rpca", dims = 1:30,
                           reduction.name = "umap.rpca", verbose = FALSE)

cat("RPCA integration complete:", length(unique(integrated_rpca$seurat_clusters)), "clusters\n")
# RPCA integration complete: 25 clusters

# STEP 9: Harmony Integration
#-----------------------------------------------

# Harmony works directly on PCA embedding
# First, join layers and run standard workflow
merged_harmony <- JoinLayers(merged_seurat)
merged_harmony <- RunPCA(merged_harmony, npcs = 50, verbose = FALSE)

# Run Harmony (corrects batch effects on PCA embedding)
integrated_harmony <- RunHarmony(merged_harmony, "sample_id")

# Downstream workflow
integrated_harmony <- FindNeighbors(integrated_harmony, reduction = "harmony", dims = 1:30, verbose = FALSE)
integrated_harmony <- FindClusters(integrated_harmony, resolution = 0.6, verbose = FALSE)
integrated_harmony <- RunUMAP(integrated_harmony, reduction = "harmony", dims = 1:30,
                              reduction.name = "umap.harmony", verbose = FALSE)

cat("Harmony integration complete:", length(unique(integrated_harmony$seurat_clusters)), "clusters\n")
# Harmony integration complete: 24 clusters


# STEP 10: FastMNN Integration
#-----------------------------------------------
# 2. Install SeuratWrappers from GitHub
# 1. Install the missing dependency from CRAN
install.packages("R.utils")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("satijalab/seurat-wrappers")
library(SeuratWrappers)
# Integrate using FastMNN (works on split layers)
integrated_fastmnn <- IntegrateLayers(
  object = merged_seurat,
  method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

cat("FastMNN integration complete\n")
#sanity check:
# 1. Check if the integrated MNN reduction exists
Reductions(merged_seurat)

# 2. Check if the newly integrated MNN dimensional reduction has data
print(integrated_fastmnn[["integrated.mnn"]])
# Downstream workflow
integrated_fastmnn <- FindNeighbors(integrated_fastmnn, reduction = "integrated.mnn", dims = 1:30, verbose = FALSE)
integrated_fastmnn <- FindClusters(integrated_fastmnn, resolution = 0.6, verbose = FALSE)
integrated_fastmnn <- RunUMAP(integrated_fastmnn, reduction = "integrated.mnn", dims = 1:30,
                              reduction.name = "umap.mnn", verbose = FALSE)

cat("Clustering complete:", length(unique(integrated_fastmnn$seurat_clusters)), "clusters\n")
# FastMNN integration complete: 24 clusters

# STEP 11: Generate comparison UMAPs
#-----------------------------------------------

cat("\n=== Generating Comparison Plots ===\n")

# Create a function to make standardized UMAP plots
make_comparison_plot <- function(seurat_obj, reduction_name, title, group_by = "sample_id", colors = NULL) {
  p <- DimPlot(seurat_obj, reduction = reduction_name, group.by = group_by, pt.size = 0.05) +
    ggtitle(title) +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.3, "cm"))
  
  # Apply custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }
  
  return(p)
}

# UMAPs colored by sample (assesses mixing)
p_sample_naive <- make_comparison_plot(merged_naive, "umap", "Naive Merge", 
                                       colors = sample_colors)
p_sample_cca <- make_comparison_plot(integrated_cca, "umap.cca", "CCA", 
                                     colors = sample_colors)
p_sample_rpca <- make_comparison_plot(integrated_rpca, "umap.rpca", "RPCA", 
                                      colors = sample_colors)
p_sample_harmony <- make_comparison_plot(integrated_harmony, "umap.harmony", "Harmony", 
                                         colors = sample_colors)
p_sample_mnn <- make_comparison_plot(integrated_fastmnn, "umap.mnn", "FastMNN", 
                                     colors = sample_colors)

# Combine sample-colored UMAPs
combined_samples <- (p_sample_naive | p_sample_cca | p_sample_rpca) / 
  (p_sample_harmony | p_sample_mnn | plot_spacer())


combined_samples
ggsave("plots/integration_comparison/02_integration_by_sample.png", 
       combined_samples, width = 16, height = 12, dpi = 300)

# UMAPs colored by condition (assesses biological preservation)
p_cond_naive <- make_comparison_plot(merged_naive, "umap", "Naive Merge", "condition", 
                                     colors = condition_colors)
p_cond_cca <- make_comparison_plot(integrated_cca, "umap.cca", "CCA", "condition", 
                                   colors = condition_colors)
p_cond_rpca <- make_comparison_plot(integrated_rpca, "umap.rpca", "RPCA", "condition", 
                                    colors = condition_colors)
p_cond_harmony <- make_comparison_plot(integrated_harmony, "umap.harmony", "Harmony", "condition", 
                                       colors = condition_colors)
p_cond_mnn <- make_comparison_plot(integrated_fastmnn, "umap.mnn", "FastMNN", "condition", 
                                   colors = condition_colors)

# Combine condition-colored UMAPs
combined_conditions <- (p_cond_naive | p_cond_cca | p_cond_rpca) /
  (p_cond_harmony | p_cond_mnn | plot_spacer())

combined_conditions

ggsave("plots/integration_comparison/03_integration_by_condition.png",
       combined_conditions, width = 16, height = 12, dpi = 300)


# STEP 12: Calculate integration quality metrics
#-----------------------------------------------

# Function to calculate mixing metric (local inverse Simpson's Index)
# This measures sample diversity in each cell's k-nearest neighborhood
calculate_mixing_metric <- function(seurat_obj, group_by = "sample_id", reduction = "umap", k = 50) {
  # Get embedding
  embedding <- Embeddings(seurat_obj, reduction = reduction)
  
  # Calculate k-nearest neighbors using FNN package
  nn_result <- FNN::get.knn(embedding[, 1:2], k = k)
  
  # For each cell, calculate diversity of samples in its neighborhood
  group_labels <- seurat_obj@meta.data[[group_by]]
  mixing_scores <- sapply(1:nrow(nn_result$nn.index), function(i) {
    neighbors <- nn_result$nn.index[i, ]
    neighbor_groups <- group_labels[neighbors]
    props <- table(neighbor_groups) / length(neighbor_groups)
    # Inverse Simpson's Index (higher = more mixed)
    1 / sum(props^2)
  })
  
  return(mean(mixing_scores))
}

# Calculate for all methods
methods_list <- list(
  "Naive" = list(obj = merged_naive, reduction = "umap"),
  "CCA" = list(obj = integrated_cca, reduction = "umap.cca"),
  "RPCA" = list(obj = integrated_rpca, reduction = "umap.rpca"),
  "Harmony" = list(obj = integrated_harmony, reduction = "umap.harmony"),
  "FastMNN" = list(obj = integrated_fastmnn, reduction = "umap.mnn")
)

# Calculate mixing metrics
mixing_results <- data.frame(
  method = names(methods_list),
  mixing_score = sapply(methods_list, function(x) {
    calculate_mixing_metric(x$obj, group_by = "sample_id", reduction = x$reduction)
  }),
  n_clusters = sapply(methods_list, function(x) {
    length(unique(x$obj$seurat_clusters))
  })
)

# Display results
cat("\nIntegration Quality Metrics:\n")
cat("(Higher mixing score = better sample mixing)\n\n")
print(mixing_results)

# Save metrics
write.csv(mixing_results, "metadata/integration_comparison_metrics.csv", row.names = FALSE)

# Visualize mixing scores
p_mixing <- ggplot(mixing_results, aes(x = reorder(method, -mixing_score), y = mixing_score, fill = method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(mixing_score, 2)), vjust = -0.5) +
  labs(title = "Integration Method Comparison: Sample Mixing",
       subtitle = "Higher score = better integration (samples mix within cell types)",
       x = "Integration Method",
       y = "Mixing Score (Local Inverse Simpson's Index)") +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set2")

ggsave("plots/integration_comparison/04_mixing_scores.png", p_mixing, 
       width = 10, height = 6, dpi = 300)

# STEP 13: Assess biological signal preservation
#-----------------------------------------------

# Function to test if conditions are still distinguishable after integration
# Uses UMAP distances between conditions
assess_condition_separation <- function(seurat_obj, reduction) {
  umap_coords <- seurat_obj[[reduction]]@cell.embeddings[, 1:2]
  conditions <- seurat_obj$condition
  
  # Calculate mean UMAP position for each condition
  condition_centers <- aggregate(umap_coords, by = list(conditions), FUN = mean)
  
  # Calculate pairwise distances between condition centers
  dist_matrix <- dist(condition_centers[, -1])
  mean_dist <- mean(dist_matrix)
  
  return(mean_dist)
}

# Calculate for all methods
separation_results <- data.frame(
  method = names(methods_list),
  condition_separation = sapply(methods_list, function(x) {
    assess_condition_separation(x$obj, x$reduction)
  })
)

mixing_results$condition_separation <- separation_results$condition_separation

cat("\nBiological Signal Preservation:\n")
cat("(Higher separation = conditions remain distinguishable)\n\n")
print(mixing_results[, c("method", "mixing_score", "condition_separation", "n_clusters")])

# Ideal integration: High mixing + Moderate separation
# Over-integration: High mixing + Low separation (conditions lost)
# Under-integration: Low mixing + High separation (batches not corrected)

# Visualize the trade-off
p_tradeoff <- ggplot(mixing_results, aes(x = mixing_score, y = condition_separation, label = method)) +
  geom_point(size = 4, aes(color = method)) +
  geom_text_repel(size = 4, box.padding = 0.5) +
  labs(title = "Integration Quality: Mixing vs Biological Preservation",
       subtitle = "Ideal: Upper right (high mixing + preserved biology)",
       x = "Sample Mixing Score (higher = better integration)",
       y = "Condition Separation (higher = preserved biology)") +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set2")

p_tradeoff
ggsave("plots/integration_comparison/05_mixing_vs_preservation.png", p_tradeoff,
       width = 10, height = 7, dpi = 300)


# STEP 14: Select optimal integration method
#-----------------------------------------------

# Rank methods by mixing score (primary criterion)
mixing_results <- mixing_results[order(-mixing_results$mixing_score), ]

cat("\nMethod ranking (by sample mixing):\n")
for (i in 1:nrow(mixing_results)) {
  cat(sprintf("  %d. %s (mixing: %.2f, separation: %.2f)\n",
              i, 
              mixing_results$method[i],
              mixing_results$mixing_score[i],
              mixing_results$condition_separation[i]))
}

# Select top method
best_method <- mixing_results$method[1]
cat("\nRecommended method:", best_method, "\n")

# For subsequent analysis, we'll use the best method
# In most cases, this will be Harmony or CCA
if (best_method == "Harmony") {
  integrated_final <- integrated_harmony
  reduction_final <- "umap.harmony"
} else if (best_method == "CCA") {
  integrated_final <- integrated_cca
  reduction_final <- "umap.cca"
} else if (best_method == "RPCA") {
  integrated_final <- integrated_rpca
  reduction_final <- "umap.rpca"
} else if (best_method == "FastMNN") {
  integrated_final <- integrated_fastmnn
  reduction_final <- "umap.mnn"
} else {
  # Default to CCA if method not found
  integrated_final <- integrated_cca
  reduction_final <- "umap.cca"
  cat("Warning: Method not recognized, defaulting to CCA\n")
}

cat("\nUsing", best_method, "for downstream clustering analysis\n")
# Using CCA for downstream clustering analysis




# STEP 15: Perform clustering at multiple resolutions
#-----------------------------------------------

# Test multiple resolutions
resolutions <- c(0.4, 0.6, 0.8, 1.0, 1.2)

for (res in resolutions) {
  integrated_final <- FindClusters(
    integrated_final,
    resolution = res,
    verbose = FALSE
  )
  
  n_clusters <- length(unique(integrated_final@meta.data[[paste0("RNA_snn_res.", res)]]))
  cat(sprintf("Resolution %.1f: %d clusters\n", res, n_clusters))
}

# Rename cluster columns for clarity
colnames(integrated_final@meta.data) <- gsub("RNA_snn_res\\.", "clusters_res_", 
                                             colnames(integrated_final@meta.data))
# STEP 16: Visualize clustering results across resolutions
#-----------------------------------------------

# Create UMAP plots for each resolution
plot_list <- lapply(resolutions, function(res) {
  cluster_col <- paste0("clusters_res_", res)
  n_clusters <- length(unique(integrated_final@meta.data[[cluster_col]]))
  
  DimPlot(integrated_final, 
          reduction = reduction_final,
          group.by = cluster_col,
          label = TRUE,
          label.size = 4,
          pt.size = 0.3) +
    ggtitle(paste0("Resolution ", res, " (", n_clusters, " clusters)")) +
    NoLegend()
})

# Combine plots
combined_resolutions <- wrap_plots(plot_list, ncol = 2)
ggsave("plots/clustering/06_multi_resolution_clustering.png", 
       combined_resolutions, width = 14, height = 14, dpi = 300)


# STEP 17: Evaluate optimal clustering resolution
#-----------------------------------------------

# We'll use silhouette scores to assess cluster quality at each resolution
# Silhouette score measures how similar cells are to their own cluster vs other clusters
# Values range from -1 to 1 (higher is better)
library(cluster)

# Map UMAP reduction name to corresponding integrated reduction name
# We use integrated space (not UMAP) for silhouette calculation because:
# - Integrated space preserves more information (30+ dimensions vs UMAP's 2)
# - Silhouette scores are more meaningful in the space where clustering was performed
integrated_reduction <- switch(reduction_final,
                               "umap.cca" = "integrated.cca",
                               "umap.rpca" = "integrated.rpca",
                               "umap.harmony" = "harmony",
                               "umap.mnn" = "integrated.mnn",
                               "umap" = "pca"  # Fallback for naive merge
)

cat("\nEvaluating clustering resolutions using", integrated_reduction, "space\n")

# For large datasets (>10,000 cells), subsample for silhouette calculation (use 5000 cells in our case)
# Computing distance matrices for all cells is computationally prohibitive
n_cells <- ncol(integrated_final)
max_cells_for_silhouette <- 5000

if (n_cells > max_cells_for_silhouette) {
  cat("Dataset has", n_cells, "cells - subsampling", max_cells_for_silhouette, "cells for silhouette calculation\n")
  set.seed(42)
  subsample_idx <- sample(1:n_cells, max_cells_for_silhouette)
} else {
  subsample_idx <- 1:n_cells
}

silhouette_scores <- sapply(resolutions, function(res) {
  cluster_col <- paste0("clusters_res_", res)
  
  # Get clusters for subsampled cells
  clusters <- as.numeric(integrated_final@meta.data[[cluster_col]][subsample_idx])
  
  # Get integrated coordinates for subsampled cells (use first 30 dimensions)
  coords <- Embeddings(integrated_final, reduction = integrated_reduction)[subsample_idx, ]
  coords <- coords[, 1:min(30, ncol(coords))]
  
  # Calculate silhouette (only if we have at least 2 clusters)
  if (length(unique(clusters)) > 1) {
    dist_matrix <- dist(coords)
    sil <- silhouette(clusters, dist_matrix)
    mean(sil[, 3])
  } else {
    NA  # Return NA if only 1 cluster
  }
})

# Create resolution comparison table
resolution_comparison <- data.frame(
  resolution = resolutions,
  n_clusters = sapply(resolutions, function(res) {
    cluster_col <- paste0("clusters_res_", res)
    length(unique(integrated_final@meta.data[[cluster_col]]))
  }),
  silhouette_score = silhouette_scores
)

cat("\nClustering resolution comparison:\n")
print(resolution_comparison)

##resolution n_clusters silhouette_score
##1        0.4         16        0.2192736
##2        0.6         16        0.1989474
##3        0.8         21        0.1756458
##4        1.0         21        0.1709990
##5        1.2         23        0.1618541

# Recommend optimal resolution
# Choose resolution with good silhouette score and interpretable cluster number
optimal_idx <- which.max(resolution_comparison$silhouette_score)
optimal_resolution <- resolution_comparison$resolution[optimal_idx]

cat("\nRecommended resolution:", optimal_resolution, 
    "(", resolution_comparison$n_clusters[optimal_idx], "clusters )\n")

# Recommended resolution: 0.4 ( 17 clusters )

# Set the optimal clustering as the default
integrated_final$seurat_clusters <- integrated_final@meta.data[[paste0("clusters_res_", optimal_resolution)]]
Idents(integrated_final) <- "seurat_clusters"

# STEP 18: Assess cluster quality and stability
#-----------------------------------------------

# Check cluster sizes
cluster_sizes <- table(integrated_final$seurat_clusters)
cat("\nCluster sizes:\n")
print(cluster_sizes)

# Identify very small clusters (<1% of cells)
min_cluster_size <- 0.01 * ncol(integrated_final)
small_clusters <- which(cluster_sizes < min_cluster_size)

if (length(small_clusters) > 0) {
  cat("\nWarning: Small clusters detected (<1% of cells):\n")
  print(cluster_sizes[small_clusters])
  cat("   These may represent rare cell types or over-clustering\n")
}

# Check sample distribution across clusters
sample_cluster_table <- table(integrated_final$seurat_clusters, integrated_final$sample_id)
sample_cluster_pct <- prop.table(sample_cluster_table, margin = 1) * 100

# Identify sample-specific clusters (>70% from one sample)
sample_specific <- apply(sample_cluster_pct, 1, max) > 70
if (any(sample_specific)) {
  cat("\nWarning: Sample-dominated clusters detected (>70% from single sample):\n")
  print(names(which(sample_specific)))
  cat("   These may indicate incomplete batch correction\n")
}

# Visualize sample distribution in clusters
sample_dist_data <- as.data.frame.matrix(sample_cluster_pct)
sample_dist_data$cluster <- rownames(sample_dist_data)
sample_dist_long <- reshape2::melt(sample_dist_data, id.vars = "cluster",
                                   variable.name = "sample", value.name = "percentage")

p_sample_dist <- ggplot(sample_dist_long, aes(x = cluster, y = percentage, fill = sample)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Sample Distribution Across Clusters",
       subtitle = "Check for sample-dominated clusters (poor integration)",
       x = "Cluster", y = "Percentage of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

p_sample_dist
ggsave("plots/clustering/08_sample_distribution_clusters.png", p_sample_dist,
       width = 12, height = 6, dpi = 300)

# Check condition distribution
condition_dist <- table(integrated_final$seurat_clusters, integrated_final$condition)
condition_dist_pct <- prop.table(condition_dist, margin = 1) * 100

cat("\nCondition distribution across clusters:\n")
print(round(condition_dist_pct, 1))

# Visualize condition distribution
condition_dist_data <- as.data.frame.matrix(condition_dist_pct)
condition_dist_data$cluster <- rownames(condition_dist_data)
condition_dist_long <- reshape2::melt(condition_dist_data, id.vars = "cluster",
                                      variable.name = "condition", value.name = "percentage")

p_condition_dist <- ggplot(condition_dist_long, aes(x = cluster, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Condition Distribution Across Clusters",
       subtitle = "Check if biological conditions are preserved",
       x = "Cluster", y = "Percentage of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("noncancerous" = "#2E86AB", 
                               "cancerous" = "#F18F01"))
p_condition_dist
ggsave("plots/clustering/09_condition_distribution_clusters.png", p_condition_dist,
       width = 12, height = 6, dpi = 300)


# STEP 19: Create comprehensive visualization of final results
#-----------------------------------------------

# Main UMAP with clusters
p_final_clusters <- DimPlot(integrated_final, 
                            reduction = reduction_final,
                            group.by = "seurat_clusters",
                            label = TRUE,
                            label.size = 5,
                            pt.size = 0.05) +
  ggtitle(paste0("Final Clustering (", best_method, " Integration)")) +
  theme(plot.title = element_text(face = "bold", size = 14))

# Split by condition to show preservation
p_by_condition <- DimPlot(integrated_final,
                          reduction = reduction_final,
                          group.by = "seurat_clusters",
                          split.by = "condition",
                          label = TRUE,
                          label.size = 4,
                          pt.size = 0.05,
                          ncol = 2) +
  ggtitle("Clusters by Treatment Condition") +
  theme(strip.text = element_text(face = "bold"))

# Colored by sample (integration quality check)
p_by_sample <- DimPlot(integrated_final,
                       reduction = reduction_final,
                       group.by = "sample_id",
                       pt.size = 0.05,
                       cols = sample_colors) +
  ggtitle("Sample Mixing (Integration Quality)") +
  theme(legend.text = element_text(size = 7))

# Colored by condition
p_by_condition_single <- DimPlot(integrated_final,
                                 reduction = reduction_final,
                                 group.by = "condition",
                                 pt.size = 0.05,
                                 cols = condition_colors) +
  ggtitle("Treatment Conditions")

# Combine final visualizations
combined_final <- (p_final_clusters | p_by_condition_single) /
  (p_by_sample | plot_spacer())

ggsave("plots/clustering/10_final_integrated_clustered.png", combined_final,
       width = 16, height = 12, dpi = 300)

# Split by condition detailed view
ggsave("plots/clustering/11_clusters_by_condition_split.png", p_by_condition,
       width = 18, height = 6, dpi = 300)

# STEP 20: Save integrated and clustered data
#-----------------------------------------------

cat("\n=== Saving Final Results ===\n")

# Join layers for downstream analysis
integrated_final <- JoinLayers(integrated_final)

# Save the final integrated object
saveRDS(integrated_final, "integrated_data/integrated_clustered_seurat.rds")

# Save metadata
write.csv(integrated_final@meta.data, 
          "metadata/integrated_cell_metadata.csv",
          row.names = TRUE)

# Create comprehensive summary
integration_summary <- data.frame(
  dataset = "GSE212966",
  n_samples = length(unique(integrated_final$sample_id)),
  n_cells_total = ncol(integrated_final),
  n_genes = nrow(integrated_final),
  integration_method = best_method,
  optimal_resolution = optimal_resolution,
  n_clusters = length(unique(integrated_final$seurat_clusters)),
  mixing_score = mixing_results$mixing_score[mixing_results$method == best_method],
  condition_separation = mixing_results$condition_separation[mixing_results$method == best_method]
)

write.csv(integration_summary, "metadata/integration_summary.csv", row.names = FALSE)

list.dirs("~/GSE212966_scRNA")

# 1. Get a list of all files inside the GSE212966_scRNA folder
files_to_zip <- list.files("~/GSE212966_scRNA", full.names = TRUE, recursive = TRUE)

# 2. Create the zip file using that list
zip(zipfile = "GSE212966_scRNA_integration.zip", files = files_to_zip)

# 3. Move it to the folder that syncs with Galaxy history
gx_put('GSE212966_scRNA_integration.zip') 
