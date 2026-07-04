
# Install Bioconductor manager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# SingleR and Reference Datasets
BiocManager::install(c(
  "SingleR",
  "celldex",
  "SingleCellExperiment"
), update = FALSE, ask = FALSE)

# Install Seurat 5 (the star of this tutorial!)
install.packages("Seurat")

# Install SeuratObject (Seurat's data structure package)
install.packages("SeuratObject")

# scCATCH (Tissue-Specific Markers)
install.packages("scCATCH")

# scType Dependencies
install.packages(c(
  "HGNChelper",
  "openxlsx"
))
# Install visualization packages
install.packages(c(
  "ggplot2",        # For custom plots when needed
  "patchwork",      # For combining Seurat plots
  "dplyr"           # Data manipulation
))

# Visualization and Utilities
install.packages(c(
  "ggalluvial",
  "scales",
  "viridis"
))

# STEP 1: Load required libraries
#-----------------------------------------------

# Core single-cell analysis
library(Seurat)
library(SeuratObject)

# Automated annotation methods
library(SingleR)
library(celldex)
library(scCATCH)

# Data structures
library(SingleCellExperiment)


# Marker processing
library(HGNChelper)

# Visualization and data manipulation
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(patchwork)
library(scales)
library(viridis)


# Set working directory
setwd("~/GSE212966_scRNA/cell_type_annotation")

# Create output directories
dir.create("plots", showWarnings = FALSE)
dir.create("plots/manual_annotation", showWarnings = FALSE)
dir.create("plots/automated_annotation", showWarnings = FALSE)
dir.create("plots/method_comparison", showWarnings = FALSE)
dir.create("annotations", showWarnings = FALSE)

# Set random seed for reproducibility
set.seed(42)

# Configure plotting defaults
theme_set(theme_classic(base_size = 12))


# STEP 2: Load integrated data from Part 3
#-----------------------------------------------
# 1. Define your directories
source_dir <- "~/galaxy_inputs"
dest_dir <- "~/GSE212966_scRNA" # Ensure this folder exists

# 2. Get a list of all .rds.rds files
files <- list.files(source_dir, pattern = "\\.rds\\.rds$", full.names = TRUE)

# 3. Process the files
for (f in files) {
  # Define the new name for the destination (removing one .rds)
  new_name <- gsub("\\.rds\\.rds$", ".rds", basename(f))
  dest_path <- file.path(dest_dir, new_name)
  
  # Check if the file already exists in destination to avoid overwriting
  if (file.exists(dest_path)) {
    message(paste("Skipping:", new_name, "- already exists in destination."))
  } else {
    # Copy the file to the new location with the new name
    file.copy(from = f, to = dest_path)
    message(paste("Copied and renamed:", basename(f), "->", new_name))
  }
}

# Path to integrated Seurat object from Part 3
integrated_path <- "~/GSE212966_scRNA/integrated_clustered_seurat.rds"

# Load the object
integrated_seurat <- readRDS(integrated_path)


# Determine which UMAP to use (from best integration method in Part 3)
available_reductions <- names(integrated_seurat@reductions)
if ("umap.cca" %in% available_reductions) {
  umap_reduction <- "umap.cca"
} else if ("umap.rpca" %in% available_reductions) {
  umap_reduction <- "umap.rpca"
} else if ("umap.mnn" %in% available_reductions) {
  umap_reduction <- "umap.mnn"
} else {
  umap_reduction <- "umap"
}


# STEP 3: Visualize clusters before annotation
#-----------------------------------------------

# Create multi-panel overview
p1 <- DimPlot(integrated_seurat, reduction = umap_reduction,
              group.by = "seurat_clusters", label = TRUE, label.size = 5,
              pt.size = 0.1) +
  ggtitle("Clusters (Pre-Annotation)") +
  theme(plot.title = element_text(face = "bold", size = 14))

p2 <- DimPlot(integrated_seurat, reduction = umap_reduction,
              group.by = "sample_id", pt.size = 0.1) +
  ggtitle("Samples") +
  theme(legend.text = element_text(size = 7))

p3 <- DimPlot(integrated_seurat, reduction = umap_reduction,
              group.by = "condition", pt.size = 0.1) +
  ggtitle("Condition") +
  scale_color_manual(values = c("noncancerous" = "#2E86AB", 
                                "cancerous" = "#F18F01"))

# Combine
p_overview <- (p1 | p2) / (p3 | plot_spacer())

p_overview
ggsave("plots/00_starting_clusters.png", p_overview, 
       width = 14, height = 10, dpi = 300)

# STEP 4: Normalize and scale data
#-----------------------------------------------

# Normalize data (log-normalization)
integrated_seurat <- NormalizeData(integrated_seurat, 
                                   normalization.method = "LogNormalize",
                                   scale.factor = 10000, 
                                   verbose = FALSE)

# Scale all genes (required for DoHeatmap and marker analysis)
integrated_seurat <- ScaleData(integrated_seurat, 
                               features = rownames(integrated_seurat), 
                               verbose = FALSE)

# Define comprehensive marker panel for PDAC cell types
canonical_markers <- list(
  "Malignant_Ductal"    = c("KRT19", "EPCAM", "S100A4", "SOX9", "S100P", "AGR2", "CEACAM6", "FA2H", "HK2", "IL1RN", "OSBPL3", "GCNT3", "DNMBP"),
  "Normal_Ductal"       = c("KRT19", "CFTR", "AMBP"),
  "Normal_Acinar"       = c("PRSS1", "REG1A", "AMY2A", "CTRB1", "CTRB2"),
  "Islet_Endocrine"     = c("INS", "GCG", "SST"),
  "Pan_Fibroblast_CAF"  = c("COL1A1", "DCN", "ACTA2", "FAP", "LUM", "PDGFRB", "TAGLN", "COL11A1", "POSTN"),
  "iCAF"                = c("IL6", "CXCL12", "PDGFRA", "COMP"),
  "apCAF"               = c("CD74", "HLA-DRA", "HLA-DRB1"),
  "Endothelial"         = c("PECAM1", "VWF", "ENG", "CD34"),
  "T_NK_cells"          = c("CD3D", "CD3E", "TRAC", "NKG7", "CD8A", "CD4", "GNLY", "NCAM1"),
  "Myeloid_Cells"       = c("CD14", "LYZ", "CD68", "CSF1R", "SPP1"),
  "B_Plasma_Cells"      = c("CD79A", "MS4A1", "MZB1", "CD19")
)

canonical_markers <- list(
  # Tumor and Ductal
  "Malignant_Ductal"    = c("KRT19", "EPCAM", "S100P", "AGR2", "CEACAM6", "SOX9", "S100A4", "FA2H", "HK2", "IL1RN", "OSBPL3", "GCNT3", "DNMBP"),
  "Normal_Ductal"       = c("KRT19", "CFTR", "AMBP", "CA2", "TFF1", "KRT7"),
  
  # Exocrine & Endocrine (Often subject to dropout)
  "Normal_Acinar"       = c("PRSS1", "REG1A", "AMY2A", "CTRB1", "CTRB2", "PNLIP", "CEL"),
  "Islet_Endocrine"     = c("INS", "GCG", "SST", "PDX1"),
  
  # Fibroblasts / CAFs
  "Pan_Fibroblast_CAF"  = c("COL1A1", "DCN", "ACTA2", "FAP", "LUM", "PDGFRB", "TAGLN", "COL11A1", "POSTN"),
  "iCAF"                = c("IL6", "CXCL12", "PDGFRA", "COMP", "CCL2"),
  "apCAF"               = c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1"),
  
  # Endothelial & Immune
  "Endothelial"         = c("PECAM1", "VWF", "ENG", "CD34", "ACKR1"),
  "T_NK_cells"          = c("CD3D", "CD3E", "TRAC", "NKG7", "CD8A", "CD4", "GNLY", "NCAM1"),
  "Myeloid_Cells"       = c("CD14", "LYZ", "CD68", "CSF1R", "SPP1", "CD163", "ITGAX"),
  "B_Plasma_Cells"      = c("MS4A1", "CD79A", "CD19", "MZB1", "IGHG1", "IGKC", "IGLC2")
)
# Check which markers are present in dataset
all_genes <- rownames(integrated_seurat)

for (cell_type in names(canonical_markers)) {
  markers <- canonical_markers[[cell_type]]
  present <- markers %in% all_genes
  
  cat(sprintf("%-20s: %d/%d markers present\n", 
              cell_type, sum(present), length(markers)))
  
  if (!all(present)) {
    missing <- markers[!present]
    cat(sprintf("  Missing: %s\n", paste(missing, collapse = ", ")))
  }
}

# Platelets: Missing GP9


# STEP 6: Visualize canonical markers
#-----------------------------------------------

# Function to create violin plots for a marker set
plot_violin_markers <- function(markers, cell_type_name, filename_suffix) {
  present_markers <- markers[markers %in% rownames(integrated_seurat)]
  
  if (length(present_markers) == 0) {
    cat("!No markers present for", cell_type_name, "\n")
    return(NULL)
  }
  
  cat("  •", cell_type_name, ":", length(present_markers), "markers\n")
  
  p <- VlnPlot(integrated_seurat,
               features = present_markers,
               group.by = "seurat_clusters",
               pt.size = 0,
               ncol = 3) &
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 11, face = "bold"))
  
  n_rows <- ceiling(length(present_markers) / 3)
  height <- max(4, n_rows * 2.5)
  
  ggsave(paste0("plots/manual_annotation/violin_", filename_suffix, ".png"),
         p, width = 14, height = height, dpi = 300)
  
  return(p)
}

# Function to create UMAP feature plots for a marker set
plot_umap_markers <- function(markers, cell_type_name, filename_suffix) {
  present_markers <- markers[markers %in% rownames(integrated_seurat)]
  
  if (length(present_markers) == 0) {
    return(NULL)
  }
  
  p <- FeaturePlot(integrated_seurat,
                   features = present_markers,
                   reduction = umap_reduction,
                   ncol = 4,
                   pt.size = 0.1) &
    theme(plot.title = element_text(size = 10),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right")
  
  n_rows <- ceiling(length(present_markers) / 4)
  height <- max(4, n_rows * 3)
  
  ggsave(paste0("plots/manual_annotation/umap_", filename_suffix, ".png"),
         p, width = 16, height = height, dpi = 300)
  
  return(p)
}

# Generate plots for ALL cell types
for (cell_type_key in names(canonical_markers)) {
  markers <- canonical_markers[[cell_type_key]]
  cell_type_name <- gsub("_", " ", cell_type_key)
  
  plot_violin_markers(markers, cell_type_name, cell_type_key)
  plot_umap_markers(markers, cell_type_name, cell_type_key)
}

# Create comprehensive DotPlot showing all key markers
all_markers <- unique(unlist(canonical_markers))
all_markers_present <- all_markers[all_markers %in% rownames(integrated_seurat)]

p_dotplot <- DotPlot(integrated_seurat,
                     features = all_markers_present,
                     group.by = "seurat_clusters") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 8)) +
  labs(title = "All Canonical Markers by Cluster",
       subtitle = paste(length(all_markers_present), "markers across all PBMC cell types"))

p_dotplot
ggsave("plots/manual_annotation/dotplot_all_canonical_markers.png", p_dotplot,
       width = 12, height = max(8, length(all_markers_present) * 0.15), dpi = 300)


# STEP 7: Find cluster-specific marker genes
#-----------------------------------------------

# Find markers for each cluster vs all other cells
cluster_markers <- FindAllMarkers(
  integrated_seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

# Remove ribosomal and mitochondrial genes
cluster_markers_filtered <- cluster_markers %>%
  filter(!grepl("^RP[SL]", gene)) %>%
  filter(!grepl("^MT-", gene))

# Get top 5 markers per cluster
top_markers <- cluster_markers_filtered %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

cat("\nTop 5 markers per cluster:\n")
print(top_markers %>% select(cluster, gene, avg_log2FC, pct.1, pct.2))

# Save all markers
write.csv(cluster_markers_filtered, 
          "annotations/cluster_markers_all.csv",
          row.names = FALSE)

write.csv(top_markers,
          "annotations/top5_markers_per_cluster.csv",
          row.names = FALSE)

# Create heatmap of top markers
top_genes <- top_markers$gene

p_heatmap <- DoHeatmap(
  integrated_seurat,
  features = top_genes,
  group.by = "seurat_clusters",
  size = 3
) +
  theme(axis.text.y = element_text(size = 6)) +
  labs(title = "Top 5 Markers per Cluster")

p_heatmap
ggsave("plots/manual_annotation/06_heatmap_top_markers.png", p_heatmap,
       width = 12, height = 14, dpi = 300)

# STEP 8: Manual cell type assignment
#-----------------------------------------------

# EXAMPLE MAPPING - CUSTOMIZE FOR YOUR DATA
cluster_to_celltype <- c(
  "0" = "NK cells",
  "1" = "NK cells",
  "2" = "apCAF",
  "3" = "Pan_Fibroblast_CAF",
  "4" = "Malignant",
  "5" = "B_Plasma_Cells",
  "6" = "Myeloid_Cells",
  "7" = "endothelial",
  "8" = "Pan_Fibroblast_CAF",
  "9" = "NK cells",
  "10" = "Malignant",
  "11" = "Malignant",
  "12" = "Pan_Fibroblast_CAF",
  "13" = "apCAF",#
  "14" = "apCAF",
  "15" = "apCAF"
  )

# Apply to all cells
cluster_ids <- as.character(integrated_seurat$seurat_clusters)
integrated_seurat$manual_annotation <- unname(cluster_to_celltype[cluster_ids])

# Visualize manual annotation
p_manual <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "manual_annotation",
  label = TRUE,
  label.size = 4,
  pt.size = 0.1,
  repel = TRUE
) +
  ggtitle("Manual Cell Type Annotation") +
  theme(plot.title = element_text(face = "bold", size = 14))

p_manual
ggsave("plots/manual_annotation/07_manual_annotation_umap.png", p_manual,
       width = 10, height = 8, dpi = 300)

# Split by condition
p_manual_split <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "manual_annotation",
  split.by = "condition",
  pt.size = 0.1,
  ncol = 2
) +
  ggtitle("Manual Annotation by Condition")

ggsave("plots/manual_annotation/08_manual_annotation_by_condition.png",
       p_manual_split, width = 14, height = 6, dpi = 300)


# STEP 9: SingleR with multiple reference datasets
#-----------------------------------------------

# Convert Seurat to SingleCellExperiment for SingleR
sce <- as.SingleCellExperiment(integrated_seurat)

# Reference 1: HumanPrimaryCellAtlas (Broad)
hpca_ref <- celldex::HumanPrimaryCellAtlasData()

singler_hpca <- SingleR(
  test = sce,
  ref = hpca_ref,
  labels = hpca_ref$label.main,
  assay.type.test = "logcounts"
)

integrated_seurat$singler_hpca <- singler_hpca$labels

# Reference 2: BlueprintEncodeData (cell type-Specific)
blue_ref <- celldex::BlueprintEncodeData()

singler_blueprint <- SingleR(
  test = sce,
  ref = blue_ref,
  labels = blue_ref$label.main,
  assay.type.test = "logcounts"
)

integrated_seurat$singler_blueprint <- singler_blueprint$labels

# Reference 3: DatabaseImmuneCellExpression
dice_ref <- celldex::DatabaseImmuneCellExpressionData()

singler_dice <- SingleR(
  test = sce,
  ref = dice_ref,
  labels = dice_ref$label.main,
  assay.type.test = "logcounts"
)

integrated_seurat$singler_dice <- singler_dice$labels

# STEP 10: Compare SingleR references
#-----------------------------------------------

# Visualize all three references
p_hpca <- DimPlot(integrated_seurat, reduction = umap_reduction,
                  group.by = "singler_hpca", pt.size = 0.1) +
  labs(title = "HPCA Reference") +
  theme(legend.text = element_text(size = 7))

p_blueprint <- DimPlot(integrated_seurat, reduction = umap_reduction,
                    group.by = "singler_blueprint", pt.size = 0.1) +
  labs(title = "Blueprint Reference") +
  theme(legend.text = element_text(size = 7))

p_dice <- DimPlot(integrated_seurat, reduction = umap_reduction,
                  group.by = "singler_dice", pt.size = 0.1) +
  labs(title = "DICE Reference") +
  theme(legend.text = element_text(size = 7))

p_manual_ref <- DimPlot(integrated_seurat, reduction = umap_reduction,
                        group.by = "manual_annotation", pt.size = 0.1) +
  labs(title = "Manual") +
  theme(legend.text = element_text(size = 7))

p_ref_comparison <- (p_manual_ref | p_hpca) / (p_blueprint | p_dice)

ggsave("plots/automated_annotation/10_singler_reference_comparison.png",
       p_ref_comparison, width = 16, height = 12, dpi = 300)


# For this tutorial, we'll use HPCA as the primary SingleR annotation
integrated_seurat$singler_annotation <- integrated_seurat$singler_hpca

# Visualize chosen reference
p_singler_final <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "singler_annotation",
  label = TRUE,
  label.size = 3,
  pt.size = 0.1,
  repel = TRUE
) +
  ggtitle("SingleR Annotation (HPCA Reference)") +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("plots/automated_annotation/11_singler_final_annotation.png",
       p_singler_final, width = 11, height = 8, dpi = 300)


# STEP 11: Load scType functions
#-----------------------------------------------

# Load scType functions from GitHub
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# STEP 12: Load scType marker database
#-----------------------------------------------

# Download the full database
db_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
temp_file <- tempfile(fileext = ".xlsx")
download.file(db_url, destfile = temp_file, mode = "wb")

# Read the database for immune system
library(openxlsx)
db <- read.xlsx(temp_file)
print(db$tissueType)
db_pancreas <- db[db$tissueType == "Pancreas", ]

# Prepare gene sets
gs_list <- gene_sets_prepare(temp_file, "Pancreas")

cat("Available cell types in database:\n")
print(names(gs_list$gs_positive))


# STEP 13: Calculate scType scores
#-----------------------------------------------

# Get scaled data
scaled_data <- as.matrix(LayerData(integrated_seurat, layer = "scale.data"))

# Calculate scores
es_max <- sctype_score(
  scRNAseqData = scaled_data,
  scaled = TRUE,
  gs = gs_list$gs_positive,
  gs2 = gs_list$gs_negative
)

# Assign cell types to clusters
clusters <- integrated_seurat$seurat_clusters
cL_results <- do.call("rbind", lapply(unique(clusters), function(cl) {
  es_max_cl <- sort(rowSums(es_max[, clusters == cl]), decreasing = TRUE)
  
  top_score <- es_max_cl[1]
  top_type <- names(es_max_cl)[1]
  
  data.frame(
    cluster = cl,
    type = top_type,
    scores = top_score,
    ncells = sum(clusters == cl)
  )
}))

# Filter low-confidence assignments
score_threshold <- quantile(cL_results$scores, 0.6)
cL_results$type[cL_results$scores < score_threshold] <- "Unknown"

# Create mapping and apply to cells
cluster_to_sctype <- setNames(cL_results$type, cL_results$cluster)
integrated_seurat$sctype_annotation <- unname(cluster_to_sctype[
  as.character(integrated_seurat$seurat_clusters)
])

# Visualize
p_sctype <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "sctype_annotation",
  label = TRUE,
  label.size = 4,
  pt.size = 0.1,
  repel = TRUE
) +
  ggtitle("scType: Marker-Based Annotation") +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("plots/automated_annotation/12_sctype_annotation.png", p_sctype,
       width = 11, height = 8, dpi = 300)

# Save scType results
write.csv(cL_results, "annotations/sctype_cluster_scores.csv", row.names = FALSE)
head(cL_results)
print(cL_results)



#-----------------------------------------------
# STEP 14: scCATCH annotation
#-----------------------------------------------

# Prepare data for scCATCH
obj_sccatch <- createscCATCH(
  data = LayerData(integrated_seurat, layer = "data"),
  cluster = as.character(integrated_seurat$seurat_clusters)
)
# Access the internal cellmatch data
data("cellmatch")

# See unique cancer types available
unique(cellmatch$cancer)
unique(cellmatch$cancer[cellmatch$tissue == "Pancreas"])
# Find marker genes for each cluster
obj_sccatch <- findmarkergene(
  object = obj_sccatch,
  species = "Human",
  marker = cellmatch,
  tissue = "Pancreas",
  cancer = "Pancreatic Cancer",
  cell_min_pct = 0.25,
  logfc = 0.25,
  pvalue = 0.05
)

# Find cell types based on markers
obj_sccatch <- findcelltype(obj_sccatch)

# Extract results
sccatch_celltype <- obj_sccatch@celltype

# Create cluster to cell type mapping
cluster_to_sccatch <- setNames(
  sccatch_celltype$cell_type,
  sccatch_celltype$cluster
)

# Apply to all cells
integrated_seurat$sccatch_annotation <- unname(cluster_to_sccatch[
  as.character(integrated_seurat$seurat_clusters)
])

# Handle unmapped clusters
integrated_seurat$sccatch_annotation[
  is.na(integrated_seurat$sccatch_annotation)
] <- "Unknown"

# Visualize
p_sccatch <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "sccatch_annotation",
  label = TRUE,
  label.size = 4,
  pt.size = 0.1,
  repel = TRUE
) +
  ggtitle("scCATCH: Tissue-Specific Database") +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("plots/automated_annotation/13_sccatch_annotation.png", p_sccatch,
       width = 11, height = 8, dpi = 300)

# Save detailed results
write.csv(sccatch_celltype, 
          "annotations/sccatch_cluster_assignments.csv",
          row.names = FALSE)

# Save marker evidence
marker_evidence <- obj_sccatch@markergene
write.csv(marker_evidence,
          "annotations/sccatch_marker_evidence.csv",
          row.names = FALSE)


# STEP 15: Generate confusion matrices
#-----------------------------------------------

# Function to create confusion matrix
create_confusion_matrix <- function(method1_labels, method2_labels, 
                                    method1_name, method2_name) {
  conf_table <- table(
    Method1 = method1_labels,
    Method2 = method2_labels
  )
  
  # Optional: Dynamically rename the table dimensions using your arguments
  names(dimnames(conf_table)) <- c(method1_name, method2_name)
  
  # Crucial: You must explicitly return the table, and close the function with }
  return(as.data.frame.matrix(conf_table)) 
}


# Compare each automated method to manual
conf_manual_vs_singler <- create_confusion_matrix(
  integrated_seurat$manual_annotation,
  integrated_seurat$singler_annotation,
  "Manual", "SingleR"
)

conf_manual_vs_sctype <- create_confusion_matrix(
  integrated_seurat$manual_annotation,
  integrated_seurat$sctype_annotation,
  "Manual", "scType"
)

conf_manual_vs_sccatch <- create_confusion_matrix(
  integrated_seurat$manual_annotation,
  integrated_seurat$sccatch_annotation,
  "Manual", "scCATCH"
)

# Save confusion matrices
write.csv(conf_manual_vs_singler, 
          "annotations/confusion_matrix_manual_vs_singler.csv")
write.csv(conf_manual_vs_sctype,
          "annotations/confusion_matrix_manual_vs_sctype.csv")
write.csv(conf_manual_vs_sccatch,
          "annotations/confusion_matrix_manual_vs_sccatch.csv")  


# STEP 16: Create Sankey diagrams
#-----------------------------------------------

# Function to prepare Sankey data
prepare_sankey_data <- function(method1_labels, method2_labels, 
                                method1_name, method2_name) {
  sankey_df <- data.frame(
    Method1 = method1_labels,
    Method2 = method2_labels
  ) %>%
    group_by(Method1, Method2) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(Comparison = paste(method1_name, "vs", method2_name))
  
  return(sankey_df)
}

# Prepare data for each comparison
sankey_manual_singler <- prepare_sankey_data(
  integrated_seurat$manual_annotation,
  integrated_seurat$singler_annotation,
  "Manual", "SingleR"
)

sankey_manual_sctype <- prepare_sankey_data(
  integrated_seurat$manual_annotation,
  integrated_seurat$sctype_annotation,
  "Manual", "scType"
)

sankey_manual_sccatch <- prepare_sankey_data(
  integrated_seurat$manual_annotation,
  integrated_seurat$sccatch_annotation,
  "Manual", "scCATCH"
)

# Function to create Sankey plot
create_sankey_plot <- function(sankey_data, method1_name, method2_name) {
  p <- ggplot(sankey_data,
              aes(axis1 = Method1, axis2 = Method2, y = Count)) +
    geom_alluvium(aes(fill = Method1), width = 1/12, alpha = 0.7) +
    geom_stratum(width = 1/12, fill = "white", color = "grey") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
              size = 3, fontface = "bold") +
    scale_x_discrete(limits = c(method1_name, method2_name), 
                     expand = c(0.05, 0.05)) +
    labs(title = paste(method1_name, "vs", method2_name),
         subtitle = "Straight flows = agreement | Crossed flows = disagreement",
         y = "Number of Cells") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(face = "bold"))
  
  return(p)
}

# Create Sankey plots
p_sankey_singler <- create_sankey_plot(sankey_manual_singler, "Manual", "SingleR")
ggsave("plots/method_comparison/14_sankey_manual_vs_singler.png",
       p_sankey_singler, width = 12, height = 10, dpi = 300)

p_sankey_sctype <- create_sankey_plot(sankey_manual_sctype, "Manual", "scType")
ggsave("plots/method_comparison/15_sankey_manual_vs_sctype.png",
       p_sankey_sctype, width = 12, height = 10, dpi = 300)

p_sankey_sccatch <- create_sankey_plot(sankey_manual_sccatch, "Manual", "scCATCH")
ggsave("plots/method_comparison/16_sankey_manual_vs_sccatch.png",
       p_sankey_sccatch, width = 12, height = 10, dpi = 300)


# STEP 17: Create comprehensive UMAP comparison
#-----------------------------------------------

# Create individual plots with consistent styling
p_manual_comp <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "manual_annotation",
  pt.size = 0.1
) +
  labs(title = "Manual") +
  theme(legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold"))

p_singler_comp <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "singler_annotation",
  pt.size = 0.1
) +
  labs(title = "SingleR (HPCA)") +
  theme(legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold"))

p_sctype_comp <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "sctype_annotation",
  pt.size = 0.1
) +
  labs(title = "scType (Markers)") +
  theme(legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold"))

p_sccatch_comp <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "sccatch_annotation",
  pt.size = 0.1
) +
  labs(title = "scCATCH (Tissue DB)") +
  theme(legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold"))

# Combine all four methods
p_all_methods <- (p_manual_comp | p_singler_comp) /
  (p_sctype_comp | p_sccatch_comp)

ggsave("plots/method_comparison/17_all_methods_comparison.png",
       p_all_methods, width = 18, height = 14, dpi = 300)


# STEP 18: Export annotations for manual comparison
#-----------------------------------------------

# Create comprehensive annotation table
annotation_comparison <- data.frame(
  cell_barcode = colnames(integrated_seurat),
  cluster = integrated_seurat$seurat_clusters,
  manual = integrated_seurat$manual_annotation,
  singler = integrated_seurat$singler_annotation,
  sctype = integrated_seurat$sctype_annotation,
  sccatch = integrated_seurat$sccatch_annotation,
  sample_id = integrated_seurat$sample_id,
  condition = integrated_seurat$condition,
  stringsAsFactors = FALSE
)

# Export to CSV for manual review
write.csv(annotation_comparison,
          "annotations/all_methods_comparison.csv",
          row.names = FALSE)

# Create cluster-level summary (easier to review)
cluster_summary <- annotation_comparison %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    manual_type = names(sort(table(manual), decreasing = TRUE))[1],
    singler_type = names(sort(table(singler), decreasing = TRUE))[1],
    sctype_type = names(sort(table(sctype), decreasing = TRUE))[1],
    sccatch_type = names(sort(table(sccatch), decreasing = TRUE))[1],
    .groups = "drop"
  )

write.csv(cluster_summary,
          "annotations/cluster_level_comparison.csv",
          row.names = FALSE)
# STEP 19: Import consensus annotations
#-----------------------------------------------

# Read your manually edited cluster consensus file
# (Save your Excel/Google Sheets file as CSV first)
cluster_consensus <- read.csv("annotations/cluster_level_comparison.csv",
                              stringsAsFactors = FALSE)


colnames(cluster_consensus)

library(dplyr)

# 1. Create the consensus column (taking the most frequent label per row)
cluster_consensus <- cluster_consensus %>%
  rowwise() %>%
  mutate(consensus_annotation = {
    # Gather all labels for this cluster
    labels <- c(manual_type, singler_type, sctype_type, sccatch_type)
    # Remove any NA values
    labels <- na.omit(labels)
    # Return the most frequent label
    if(length(labels) == 0) NA else names(sort(table(labels), decreasing = TRUE))[1]
  }) %>%
  ungroup()

# 2. Now your original code will work perfectly!
cluster_to_consensus <- setNames(
  cluster_consensus$consensus_annotation,
  cluster_consensus$cluster
)
# Create mapping from cluster to consensus annotation
cluster_to_consensus <- setNames(
  cluster_consensus$consensus_annotation,  # Column you created
  cluster_consensus$cluster
)

# Apply consensus to all cells
integrated_seurat$final_annotation <- unname(cluster_to_consensus[
  as.character(integrated_seurat$seurat_clusters)
])

# Visualize final consensus
p_consensus <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "final_annotation",
  label = TRUE,
  label.size = 4,
  pt.size = 0.1,
  repel = TRUE
) +
  ggtitle("Final Consensus Annotation") +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("plots/method_comparison/18_final_consensus_annotation.png",
       p_consensus, width = 11, height = 8, dpi = 300)

# Split by condition to check biological validity
p_consensus_split <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "final_annotation",
  split.by = "condition",
  pt.size = 0.1,
  ncol = 2
) +
  ggtitle("Final Consensus by Condition")

ggsave("plots/method_comparison/19_consensus_by_condition.png",
       p_consensus_split, width = 14, height = 6, dpi = 300)


# STEP 20: Save final annotated Seurat object
#-----------------------------------------------

# Save the complete annotated Seurat object
saveRDS(integrated_seurat, "annotations/integrated_annotated_seurat.rds")

# Create metadata export
metadata_export <- integrated_seurat@meta.data %>%
  select(
    sample_id, condition, patient_id,
    seurat_clusters,
    manual_annotation,
    singler_annotation,
    sctype_annotation,
    sccatch_annotation,
    final_annotation
  )

write.csv(metadata_export,
          "annotations/cell_metadata_final.csv",
          row.names = TRUE)
list.dirs("~/GSE212966_scRNA")

# 1. Get a list of all files inside the GSE212966_scRNA folder
files_to_zip <- list.files("~/GSE212966_scRNA", full.names = TRUE, recursive = TRUE)

# 2. Create the zip file using that list
zip(zipfile = "GSE212966_scRNA_cellannotation.zip", files = files_to_zip)

# 3. Move it to the folder that syncs with Galaxy history
gx_put('GSE212966_scRNA_cellannotation.zip') 
