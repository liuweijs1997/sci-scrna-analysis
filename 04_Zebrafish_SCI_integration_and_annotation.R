# ==================================================================================
# Zebrafish Spinal Cord Injury (SCI) scRNA-seq Analysis
#
# This script performs a complete workflow for multiple zebrafish scRNA-seq samples:
# 1. Loads raw 10x Genomics data for each of the 6 samples.
# 2. Applies sample-specific Quality Control (QC) filtering.
# 3. Integrates all samples using the Seurat v4 workflow to correct for batch effects.
# 4. Performs dimensionality reduction, clustering, and visualization.
# 5. Annotates cell types and performs final filtering based on annotations.
# ==================================================================================


# ---
# 1. SETUP & LIBRARY LOADING
# ---

# Load required libraries, installing any that are missing
packages <- c("Seurat", "dplyr", "tidyverse", "patchwork")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Set options for parallel processing to speed up computation
options(future.globals.maxSize = 128 * 1024^3) # Set limit to 128 GB
future::plan("multicore", workers = 8) # Adjust worker number based on your machine's cores


# ---
# 2. DATA LOADING & INDIVIDUAL SAMPLE QC
# ---

# Create a data frame to store sample information and QC parameters.
# This makes the code clean and easy to manage.
# IMPORTANT: Update 'data_path' for each sample to point to the correct
# 'filtered_feature_bc_matrix' directory.
sample_info <- data.frame(
  # A unique name for each sample run
  sample_id = c("Intact_1", "Intact_2", "C2W", "C4W", "C8W_1", "C8W_2"),
  
  # The project name, used to group biological replicates
  project_name = c("Intact", "Intact", "C2W", "C4W", "C8W", "C8W"),
  
  # Placeholder for the path to the 10X data directory
  data_path = c(
    "path/to/your/data/Intact/",
    "path/to/your/data/Intact/",
    "path/to/your/data/C2W/",
    "path/to/your/data/C4W/",
    "path/to/your/data/C8W/",
    "path/to/your/data/C8W/"
  ),
  
  # Sample-specific QC thresholds based on your notes
  max_features = c(4000, 3000, 4000, 3500, 2000, 3500),
  max_counts = c(15000, 10000, 20000, 15000, 5000, 15000)
)

# Initialize a list to store the processed Seurat objects
object_list <- list()

# Loop through each sample to perform loading and QC
for (i in 1:nrow(sample_info)) {
  sample_name <- sample_info$sample_id[i]
  project <- sample_info$project_name[i]
  path <- sample_info$data_path[i]
  
  message(paste("--- Processing sample:", sample_name, "---"))
  
  # Load the 10X Genomics data
  counts <- Read10X(data.dir = path)
  
  # Create a Seurat object
  # min.cells and min.features are initial broad filters
  seurat_obj <- CreateSeuratObject(counts = counts, project = project, 
                                   min.cells = 3, min.features = 200)
  
  # Calculate the percentage of mitochondrial genes
  # For zebrafish, mitochondrial genes are typically prefixed with "mt:"
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt:")
  
  # Apply sample-specific QC filtering based on the criteria in sample_info
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > 200 & 
      nFeature_RNA < sample_info$max_features[i] &
      nCount_RNA < sample_info$max_counts[i] &
      percent.mt < 5
  )
  
  # Standard pre-processing for integration
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # Add the processed object to the list
  object_list[[sample_name]] <- seurat_obj
}


# ---
# 3. DATA INTEGRATION
# ---

# Find integration anchors. This step identifies shared cell states across different samples.
# We use 30 dimensions, a common choice for complex datasets.
integration_anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:30, 
                                              anchor.features = 2000)

# Integrate the data to create a single, batch-corrected object
zebrafish_integrated <- IntegrateData(anchorset = integration_anchors, dims = 1:30)

# Set the default assay to the new 'integrated' assay
DefaultAssay(zebrafish_integrated) <- "integrated"


# ---
# 4. DOWNSTREAM ANALYSIS: CLUSTERING & VISUALIZATION
# ---

# Standard downstream analysis on the integrated data
zebrafish_integrated <- ScaleData(zebrafish_integrated, verbose = FALSE)
zebrafish_integrated <- RunPCA(zebrafish_integrated, npcs = 30, verbose = FALSE)

# Visualize the Elbow plot to confirm the number of PCs to use
ElbowPlot(zebrafish_integrated)

# Perform clustering and non-linear dimensionality reduction (UMAP and t-SNE)
# We use the first 20 PCs based on the Elbow plot.
zebrafish_integrated <- FindNeighbors(zebrafish_integrated, reduction = "pca", dims = 1:20)
zebrafish_integrated <- FindClusters(zebrafish_integrated, resolution = 0.6)
zebrafish_integrated <- RunUMAP(zebrafish_integrated, reduction = "pca", dims = 1:20)
zebrafish_integrated <- RunTSNE(zebrafish_integrated, reduction = "pca", dims = 1:20)

# Visualize the clusters, split by condition to check integration quality
# Re-order the sample levels for consistent plotting
zebrafish_integrated$orig.ident <- factor(zebrafish_integrated$orig.ident, levels = c('Intact', 'C2W', 'C4W', 'C8W'))
p1 <- DimPlot(zebrafish_integrated, label = TRUE, repel = TRUE)
p2 <- DimPlot(zebrafish_integrated, label = TRUE, split.by = 'orig.ident', ncol = 2)
print(p1 + p2)


# ---
# 5. CELL TYPE ANNOTATION
# ---

# Assign cell type names to clusters based on known marker genes
new_cluster_ids <- c("OL", "Microglia", "Neuron", "Erythrocyte", 'Monocyte', 
                     'Endothelial cell', 'Astrocyte', 'OPC', 'Pericyte', 
                     'ERG', 'Macrophage', 'Fibroblast', 'Undefined', 'Pigment cell')
names(new_cluster_ids) <- levels(zebrafish_integrated)
zebrafish_integrated <- RenameIdents(zebrafish_integrated, new_cluster_ids)

# Remove the 'Undefined' cluster from the dataset
zebrafish_integrated <- subset(zebrafish_integrated, idents = "Undefined", invert = TRUE)

# Re-order cell types for consistent plotting and analysis
cell_type_levels <- c('OL', 'Microglia', 'Astrocyte', 'OPC', 'Neuron', 'ERG',
                      'Monocyte', 'Macrophage', 'Pigment cell', 'Pericyte', 
                      'Fibroblast', 'Erythrocyte', 'Endothelial cell')
Idents(zebrafish_integrated) <- factor(Idents(zebrafish_integrated), levels = cell_type_levels)

# Store the final annotations and condition in the metadata
zebrafish_integrated$celltype <- Idents(zebrafish_integrated)
zebrafish_integrated$group <- zebrafish_integrated$orig.ident

# Final visualization of the annotated dataset
p_final <- DimPlot(zebrafish_integrated, group.by = 'celltype', label = TRUE, repel = TRUE) + NoLegend()
p_final_split <- DimPlot(zebrafish_integrated, group.by = 'celltype', split.by = 'group', label = TRUE, ncol = 2)
print(p_final)
print(p_final_split)


# ---
# 6. SAVE FINAL OBJECT
# ---

# IMPORTANT: Define your final output directory
output_dir <- "path/to/your/output_directory/"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Save the final, annotated Seurat object for downstream analyses
saveRDS(zebrafish_integrated, file = file.path(output_dir, "zebrafish_sci_integrated_annotated.rds"))

# --- END OF SCRIPT ---