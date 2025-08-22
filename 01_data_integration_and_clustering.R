# ==============================================================================
# Single-Cell RNA-Seq Analysis of Rat Spinal Cord Injury (SCI)
# 
# This script performs the following steps:
# 1. Loads and preprocesses four individual scRNA-seq datasets (Intact, C2W, C4W, C8W).
# 2. Performs quality control (QC) filtering on each dataset.
# 3. Integrates the four datasets to correct for batch effects.
# 4. Performs dimensionality reduction, clustering, and visualization on the integrated data.
# 5. Identifies cluster-specific marker genes.
# 6. Annotates cell types based on canonical markers.
# ==============================================================================


# ---
# 1. SETUP & LIBRARY LOADING
# ---

# Load required libraries. The following code will install any missing packages from CRAN.
packages <- c("Seurat", "dplyr", "tidyverse", "patchwork")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Increase the maximum size for global variables to handle large Seurat objects
options(future.globals.maxSize = 128 * 1024^3) # Set limit to 128 GB

# Define gene sets for Quality Control (QC)
# NOTE: Ensure these gene symbols match the annotation in your reference genome.
# Mitochondrial genes
mt_genes <- c('ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 
              'ND3', 'ND4', 'ND4L', 'ND5', 'ND6')
# Hemoglobin genes
hb_genes <- c('Hbb', 'Hbg1', 'Hbe2', 'Hbe1', 'Hba-a3', 'Hba-a2', 'Hbq1b', 'Hbz')


# ---
# 2. DATA LOADING, QC, AND PREPROCESSING OF INDIVIDUAL SAMPLES
# ---

# This section processes each of the four samples individually.
# For each sample, we load the data, calculate QC metrics, filter cells,
# and perform a standard pre-processing workflow.

# -- Sample 1: Intact --
# IMPORTANT: Replace the placeholder path with the actual location of your data directory.
intact_data_path <- "path/to/your/data/Intact/filtered_feature_bc_matrix"
intact <- CreateSeuratObject(counts = Read10X(data.dir = intact_data_path), 
                             project = "Intact", min.cells = 3, min.features = 200)
intact[["percent.mt"]] <- PercentageFeatureSet(intact, features = mt_genes)
intact[["percent.hb"]] <- PercentageFeatureSet(intact, features = hb_genes)
VlnPlot(intact, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 4)
intact <- subset(intact, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
intact <- NormalizeData(intact)
intact <- FindVariableFeatures(intact, selection.method = "vst", nfeatures = 2000)

# -- Sample 2: C2W --
c2w_data_path <- "path/to/your/data/C2W/filtered_feature_bc_matrix"
c2w <- CreateSeuratObject(counts = Read10X(data.dir = c2w_data_path), 
                          project = "C2W", min.cells = 3, min.features = 200)
c2w[["percent.mt"]] <- PercentageFeatureSet(c2w, features = mt_genes)
c2w[["percent.hb"]] <- PercentageFeatureSet(c2w, features = hb_genes)
VlnPlot(c2w, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 4)
c2w <- subset(c2w, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
c2w <- NormalizeData(c2w)
c2w <- FindVariableFeatures(c2w, selection.method = "vst", nfeatures = 2000)

# -- Sample 3: C4W --
c4w_data_path <- "path/to/your/data/C4W/filtered_feature_bc_matrix"
c4w <- CreateSeuratObject(counts = Read10X(data.dir = c4w_data_path), 
                          project = "C4W", min.cells = 3, min.features = 200)
c4w[["percent.mt"]] <- PercentageFeatureSet(c4w, features = mt_genes)
c4w[["percent.hb"]] <- PercentageFeatureSet(c4w, features = hb_genes)
VlnPlot(c4w, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 4)
c4w <- subset(c4w, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
c4w <- NormalizeData(c4w)
c4w <- FindVariableFeatures(c4w, selection.method = "vst", nfeatures = 2000)

# -- Sample 4: C8W --
c8w_data_path <- "path/to/your/data/C8W/filtered_feature_bc_matrix"
c8w <- CreateSeuratObject(counts = Read10X(data.dir = c8w_data_path), 
                          project = "C8W", min.cells = 3, min.features = 200)
c8w[["percent.mt"]] <- PercentageFeatureSet(c8w, features = mt_genes)
c8w[["percent.hb"]] <- PercentageFeatureSet(c8w, features = hb_genes)
VlnPlot(c8w, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 4)
c8w <- subset(c8w, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
c8w <- NormalizeData(c8w)
c8w <- FindVariableFeatures(c8w, selection.method = "vst", nfeatures = 2000)


# ---
# 3. DATA INTEGRATION
# ---

# Create a list of the pre-processed Seurat objects for integration
object_list <- list(Intact = intact, C2W = c2w, C4W = c4w, C8W = c8w)

# Find integration anchors using the standard Seurat workflow.
# We use the previously identified variable features for anchoring.
integration_anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:20, anchor.features = 2000)

# Integrate the datasets to create a single, batch-corrected object
sci_integrated <- IntegrateData(anchorset = integration_anchors, dims = 1:20)

# Set the default assay to the new 'integrated' assay for downstream analysis
DefaultAssay(sci_integrated) <- "integrated"


# ---
# 4. INTEGRATED ANALYSIS: CLUSTERING & VISUALIZATION
# ---

# Scale the integrated data and perform Principal Component Analysis (PCA)
sci_integrated <- ScaleData(sci_integrated, verbose = FALSE)
sci_integrated <- RunPCA(sci_integrated, npcs = 30, verbose = FALSE)

# Visualize the "elbow plot" to help determine the number of PCs to use for clustering
ElbowPlot(sci_integrated)

# Perform clustering and non-linear dimensionality reduction (UMAP and t-SNE)
# We select 20 PCs based on the elbow plot.
sci_integrated <- FindNeighbors(sci_integrated, reduction = "pca", dims = 1:20)
sci_integrated <- FindClusters(sci_integrated, resolution = 0.6)
sci_integrated <- RunUMAP(sci_integrated, reduction = "pca", dims = 1:20)
sci_integrated <- RunTSNE(sci_integrated, reduction = "pca", dims = 1:20)

# Visualize the integrated clusters
# First, reorder the sample levels for consistent plotting
sci_integrated$orig.ident <- factor(sci_integrated$orig.ident, levels = c('Intact', 'C2W', 'C4W', 'C8W'))
# Plot UMAP split by sample to check integration quality
DimPlot(sci_integrated, label = TRUE, split.by = 'orig.ident', ncol = 2)


# ---
# 5. IDENTIFYING CELL TYPE MARKERS
# ---

# Switch the default assay to "RNA" to find markers from the original, non-integrated counts
DefaultAssay(sci_integrated) <- "RNA"

# Find markers for every cluster compared to all remaining cells
# `only.pos = TRUE` returns only the positive markers for each cluster
all_markers <- FindAllMarkers(sci_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Extract and save the top 30 markers for each cluster, ranked by log2 fold change
top30_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC)

# Save the marker lists to CSV files
# IMPORTANT: Replace placeholder with your desired output directory path
output_dir <- "path/to/your/output/directory/"
write.csv(all_markers, file = file.path(output_dir, "all_cluster_markers.csv"))
write.csv(top30_markers, file = file.path(output_dir, "top30_cluster_markers.csv"))


# ---
# 6. CELL TYPE ANNOTATION
# ---

# Assign cell type names to clusters based on the marker gene analysis
new_cluster_ids <- c(
  '0' = 'Microglia', '1' = 'OL', '2' = 'OL', '3' = 'Microglia', '4' = 'Macrophages',
  '5' = 'Fibroblast', '6' = 'OL', '7' = 'OPC', '8' = 'Fibroblast', '9' = 'Fibroblast',
  '10' = 'Neuron', '11' = 'Macrophages', '12' = 'Microglia', '13' = 'Schwann',
  '14' = 'Neuron', '15' = 'Astrocyte', '16' = 'OPC', '17' = 'Neuron', '18' = 'Neuron',
  '19' = 'T cell', '20' = 'Neuron', '21' = 'NPCs', '22' = 'Pericyte',
  '23' = 'Endothelial', '24' = 'Ependymal', '25' = 'OPC', '26' = 'OL'
)
sci_integrated <- RenameIdents(sci_integrated, new_cluster_ids)

# Store the final cell type annotations and sample condition in the metadata
sci_integrated$celltype <- Idents(sci_integrated)
sci_integrated$condition <- sci_integrated$orig.ident

# Visualize the final annotated cell types on a UMAP plot
DimPlot(sci_integrated, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(sci_integrated, split.by = 'condition', label = TRUE, repel = TRUE, ncol = 2)


# ---
# 7. SAVE FINAL OBJECT
# ---

# Save the final, annotated Seurat object for future analyses and sharing
saveRDS(sci_integrated, file = file.path(output_dir, "sci_integrated_annotated.rds"))

# --- END OF SCRIPT ---