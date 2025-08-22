# ==================================================================================
# Zebrafish SCI scRNA-seq: Downstream & CellChat Analysis
#
# This script performs downstream analysis on the integrated zebrafish dataset.
# The workflow includes:
# 1. Loading the annotated Seurat object and performing final visualizations.
# 2. Converting zebrafish gene names to mouse homologs for cross-species analysis.
# 3. A streamlined, function-based workflow for running CellChat on multiple
#    samples (Intact, C2W, C4W, C8W) and a combined SCI group.
# 4. Comparative analysis and visualization of cell-cell communication networks.
# ==================================================================================


# ---
# 1. SETUP & LIBRARY LOADING
# ---

# Load required libraries, installing any that are missing
packages <- c("Seurat", "dplyr", "tidyverse", "patchwork", "homologene", "CellChat", "NMF", "ggalluvial")
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
# 2. CUSTOM FUNCTIONS
# ---


RenameGenesSeurat <- function(obj, newnames, gene.use=NULL, de.assay="RNA") {
  print("Running RenameGenesSeurat to rename genes in the Seurat object...")
  lassays <- Assays(obj); assay.use <- tryCatch(obj@reductions$pca@assay.used, error=function(e) de.assay); DefaultAssay(obj) <- de.assay
  if (is.null(gene.use)) { all_genenames <- rownames(obj) } else { all_genenames <- gene.use; obj <- subset(obj, features=gene.use) }
  order_name <- function(v1, v2, ref){ v2 <- make.names(v2, unique=T); df1 <- data.frame(v1, v2); rownames(df1) <- df1$v1; df1 <- df1[ref, , drop=FALSE]; return(df1) }
  df1 <- order_name(v1=all_genenames, v2=newnames, ref=rownames(obj)); all_genenames <- df1$v1; newnames <- df1$v2
  if ('SCT' %in% lassays) { if ('SCTModel.list' %in%  slotNames(obj@assays$SCT)) { obj@assays$SCT@SCTModel.list[[1]]@feature.attributes <- obj@assays$SCT@SCTModel.list[[1]]@feature.attributes[all_genenames,]; rownames(obj@assays$SCT@SCTModel.list[[1]]@feature.attributes) <- newnames } }
  change_assay <- function(a1=de.assay, obj, newnames=NULL, all_genenames=NULL){
    RNA <- obj@assays[[a1]]
    if (nrow(RNA) == length(newnames)) {
      if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
      if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
      if (length(VariableFeatures(RNA))) { df1 <- order_name(v1=all_genenames, v2=newnames, ref=VariableFeatures(RNA)); VariableFeatures(RNA) <- df1$v2 }
      if (length(RNA@scale.data)){ df1 <- order_name(v1=all_genenames, v2=newnames, ref=rownames(RNA@scale.data)); rownames(RNA@scale.data) <- df1$v2 }
    } else { message("Unequal gene sets: nrow(RNA) != nrow(newnames)") }; obj@assays[[a1]] <- RNA; return(obj)
  }
  for (a in lassays) { df1 <- order_name(v1=all_genenames, v2=newnames, ref=rownames(obj@assays[[a]])); obj <- change_assay(obj=obj, a1=a, newnames=df1$v2, all_genenames=df1$v1) }
  if (length(obj@reductions$pca)){ df1 <- order_name(v1=all_genenames, v2=newnames, ref=rownames(obj@reductions$pca@feature.loadings)); rownames(obj@reductions$pca@feature.loadings) <- df1$v2 }
  try(obj[[de.assay]]@meta.features <- data.frame(row.names = rownames(obj[[de.assay]]))); return(obj)
}

convert_zebrafish_to_mouse_genes <- function(zebrafish_genes_list) {
  mouse_genes <- sapply(zebrafish_genes_list, function(gene) {
    mapped_gene <- gene_mapping$mouse_gene[gene_mapping$zebrafish_gene == gene]
    if (length(mapped_gene) > 0) { return(mapped_gene[1]) } else { return(NA) }
  })
  return(mouse_genes)
}

# ---
# 3. DATA LOADING AND PREPARATION
# ---

# IMPORTANT: Define your main input and output directories
input_dir <- "path/to/your/input_data/"
output_dir <- "path/to/your/output_plots_and_data/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load the previously saved integrated and annotated Seurat object
DRSCI <- readRDS(file.path(input_dir, ""))

# Finalize annotations and factor levels for plotting
DRSCI$orig.ident <- factor(DRSCI$orig.ident, levels = c('Intact','C2W','C4W','C8W'))
cell_type_levels <- c('OL','Microglia','Astrocyte','OPC','Neuron','ERG','Monocyte','Macrophage',
                      'Pigment cell','Pericyte','Fibroblast','Erythrocyte','Endothelial cell')
Idents(DRSCI) <- factor(DRSCI$celltype, levels = cell_type_levels)
DRSCI$celltype <- Idents(DRSCI)

# Define a consistent color palette for cell types
cell_colors <- c('OL'="#894FBE", 'Microglia'="#BCAC39", 'Astrocyte'="#51DD24", 'OPC'="#C3BAD1",
                 'Neuron'="#00ACAE", 'ERG'="#D3B5AA", 'Monocyte'="#ABABAB", 'Macrophage'='#ef7fdf',
                 'Pigment cell'="#889EC0", 'Pericyte'="#3B6591", 'Fibroblast'="#E8948D",
                 'Erythrocyte'="#d03800", 'Endothelial cell'="#F2a151")

# Generate and save final UMAP plots
p_final_umap <- DimPlot(DRSCI, pt.size = 0.1, cols = cell_colors, label = TRUE)
ggsave(file.path(output_dir, "UMAP_final_annotated.pdf"), plot = p_final_umap, width = 9, height = 8)

p_final_split <- DimPlot(DRSCI, pt.size = 1.5, split.by = 'orig.ident', cols = cell_colors, label = F) + NoLegend()
ggsave(file.path(output_dir, "UMAP_final_split.pdf"), plot = p_final_split, width = 16, height = 5)


# ---
# 4. CELLCHAT ANALYSIS: A STREAMLINED WORKFLOW
# ---

## 4A. Helper Functions for Automation ##

# This function automates the zebrafish-to-mouse gene conversion for a Seurat object
prepare_for_cellchat <- function(seurat_obj_zebrafish) {
  message("--- Converting zebrafish genes to mouse homologs ---")
  
  # Get mouse homologs using homologene. Zebrafish TaxID: 7955, Mouse TaxID: 10090
  zebrafish_genes <- rownames(seurat_obj_zebrafish)
  homologs <- homologene(zebrafish_genes, inTax = 7955, outTax = 10090)
  
  # Create the global 'gene_mapping' object required by the conversion function
  gene_mapping <<- as.data.frame(homologs)
  colnames(gene_mapping) <<- c("zebrafish_gene", "mouse_gene")
  
  # Perform the conversion
  mouse_genes_mapped <- convert_zebrafish_to_mouse_genes(zebrafish_genes)
  
  # Create and filter a clean gene conversion table
  gene_conversion_df <- data.frame(
    zebrafish_gene = zebrafish_genes,
    mouse_gene = mouse_genes_mapped
  ) %>% na.omit() %>% distinct(zebrafish_gene, .keep_all = TRUE) %>% distinct(mouse_gene, .keep_all = TRUE)
  
  # Subset and rename genes in the Seurat object
  seurat_subset <- subset(seurat_obj_zebrafish, features = gene_conversion_df$zebrafish_gene)
  seurat_obj_mouse <- RenameGenesSeurat(
    obj = seurat_subset,
    newnames = gene_conversion_df$mouse_gene
  )
  
  # Final prep for CellChat
  seurat_obj_mouse <- ScaleData(seurat_obj_mouse, features = rownames(seurat_obj_mouse))
  return(seurat_obj_mouse)
}

# This function runs the core CellChat pipeline
run_cellchat_pipeline <- function(seurat_obj_mouse, sample_name) {
  message(paste("--- Running CellChat pipeline for", sample_name, "---"))
  cellchat_obj <- createCellChat(object = seurat_obj_mouse, group.by = "celltype", assay = "RNA")
  cellchat_obj@DB <- CellChatDB.mouse
  cellchat_obj <- subsetData(cellchat_obj)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  cellchat_obj <- projectData(cellchat_obj, PPI.mouse) # Project data for mouse
  cellchat_obj <- computeCommunProb(cellchat_obj, raw.use = TRUE, population.size = TRUE)
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 10)
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)
  cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name = "netP")
  return(cellchat_obj)
}

## 4B. Running the Analysis ##

# Subset to cell types of interest for communication analysis
DRSCI_glia <- subset(DRSCI, idents = c("OL", "Microglia", "Astrocyte", 'OPC', 'Neuron'))

# Split the object by sample/timepoint
sp_zebrafish <- SplitObject(DRSCI_glia, split.by = 'orig.ident')

# Apply the automated workflow to each sample
object.list <- lapply(names(sp_zebrafish), function(name) {
  seurat_obj_zebrafish <- sp_zebrafish[[name]]
  seurat_obj_mouse <- prepare_for_cellchat(seurat_obj_zebrafish)
  run_cellchat_pipeline(seurat_obj_mouse, name)
})
names(object.list) <- names(sp_zebrafish)

# Merge all CellChat objects for comparative analysis
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


# ---
# 5. CELLCHAT COMPARATIVE VISUALIZATION
# ---

cellchat_output_dir <- file.path(output_dir, "CellChat_Analysis")
if (!dir.exists(cellchat_output_dir)) dir.create(cellchat_output_dir, recursive=TRUE)

# Define colors for the subset of cells in the chat
chat_colors <- cell_colors[names(cell_colors) %in% levels(DRSCI_glia)]

# Compare total interaction counts and strengths
p_compare <- compareInteractions(cellchat, show.legend = FALSE, group = 1:4)
ggsave(file.path(cellchat_output_dir, "compare_interactions_across_time.pdf"), plot = p_compare, width = 8, height = 5)

# Visualize differential interactions (e.g., C8W vs Intact)
pdf(file.path(cellchat_output_dir, "diff_interaction_C8W_vs_Intact.pdf"), width = 10, height = 10)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, comparison = c(1, 4), label.edge = FALSE, color.use = chat_colors)
dev.off()

# Visualize signaling roles (influencers, senders, etc.)
# Calculate comparable dot sizes across all plots
num_links <- sapply(object.list, function(x) { rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count) })
weight_min_max <- c(min(num_links), max(num_links))

plot_list <- lapply(names(object.list), function(name) {
  netAnalysis_signalingRole_scatter(object.list[[name]], title = name, weight.MinMax = weight_min_max, color.use = chat_colors) + NoLegend()
})
p_roles <- patchwork::wrap_plots(plots = plot_list, ncol = 4)
ggsave(file.path(cellchat_output_dir, "signaling_role_scatter_plots.pdf"), plot = p_roles, width = 16, height = 4)

# Save the data underlying the signaling role plots
for (i in 1:length(plot_list)) {
  write.csv(plot_list[[i]]$data, file = file.path(cellchat_output_dir, paste0(names(object.list)[i], "_signaling_roles.csv")))
}

# ---
# 6. SAVE FINAL CELLCHAT OBJECT
# ---

saveRDS(cellchat, file = file.path(cellchat_output_dir, "cellchat_merged_object.rds"))

# --- END OF SCRIPT ---