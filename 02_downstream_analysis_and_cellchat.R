# ==============================================================================
# Downstream Analysis of Rat SCI scRNA-seq Data
#
# This script performs:
# 1.  Loading of the integrated Seurat object and final cell type annotation.
# 2.  Generation of publication-quality visualizations (UMAP, marker plots, QC).
# 3.  Analysis of cell type proportions across different time points.
# 4.  Conversion of rat gene names to mouse homologs for cross-species analysis.
# 5.  In-depth cell-cell communication analysis using CellChat.
# ==============================================================================


# ---
# 1. SETUP & LIBRARY LOADING
# ---

# Load required libraries. This will also install any missing packages.
packages <- c("Seurat", "dplyr", "tidyverse", "patchwork", "RColorBrewer", 
              "ggrepel", "ggalluvial", "homologene", "CellChat", "NMF")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Set future options for parallel processing
options(future.globals.maxSize = 128 * 1024^3) # 128 GB limit
future::plan("multicore", workers = 4) # Use 4 cores for parallel tasks


# ---
# 2. CUSTOM FUNCTIONS
# ---

# This function renames genes across multiple slots in a Seurat object.
RenameGenesSeurat <- function(obj, newnames, gene.use=NULL, de.assay="RNA") {
  print("Running RenameGenesSeurat...")
  lassays <- Assays(obj)
  assay.use <- obj@reductions$pca@assay.used
  DefaultAssay(obj) <- de.assay
  if (is.null(gene.use)) {
    all_genenames <- rownames(obj)
  } else {
    all_genenames <- gene.use
    obj <- subset(obj, features=gene.use)
  }
  order_name <- function(v1, v2, ref){
    v2 <- make.names(v2, unique=T)
    df1 <- data.frame(v1, v2)
    rownames(df1) <- df1$v1
    df1 <- df1[ref, ]
    return(df1)
  }
  df1 <- order_name(v1=all_genenames, v2=newnames, ref=rownames(obj))
  all_genenames <- df1$v1
  newnames <- df1$v2
  if ('SCT' %in% lassays) {
    if ('SCTModel.list' %in%  slotNames(obj@assays$SCT)) {
      obj@assays$SCT@SCTModel.list$model1@feature.attributes <- obj@assays$SCT@SCTModel.list$model1@feature.attributes[all_genenames,]
      rownames(obj@assays$SCT@SCTModel.list$model1@feature.attributes) <- newnames
    }
  }
  change_assay <- function(a1=de.assay, obj, newnames=NULL, all_genenames=NULL){
    RNA <- obj@assays[a1][[1]]
    if (nrow(RNA) == length(newnames)) {
      if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
      if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
      if (length(RNA@var.features)) {
        df1 <- order_name(v1=all_genenames, v2=newnames, ref=RNA@var.features)
        RNA@var.features <- df1$v2
      }
      if (length(RNA@scale.data)){
        df1 <- order_name(v1=all_genenames, v2=newnames, ref=rownames(RNA@scale.data))
        rownames(RNA@scale.data) <- df1$v2
      }
    } else { "Unequal gene sets: nrow(RNA) != nrow(newnames)" }
    obj@assays[a1][[1]] <- RNA
    return(obj)
  }
  for (a in lassays) {
    DefaultAssay(obj) <- a
    df1 <- order_name(v1=all_genenames, v2=newnames, ref=rownames(obj))
    obj <- change_assay(obj=obj, a1=a, newnames=df1$v2, all_genenames=df1$v1)
  }
  hvg <- VariableFeatures(obj, assay=assay.use)
  if (length(obj@reductions$pca)){
    df1 <- order_name(v1=all_genenames, v2=newnames, ref=hvg)
    df1 <- df1[rownames(obj@reductions$pca@feature.loadings),]
    rownames(obj@reductions$pca@feature.loadings) <- df1$v2
  }
  try(obj[[de.assay]]@meta.features <- data.frame(row.names = rownames(obj[[de.assay]])))
  return(obj)
}

# This function converts a list of rat gene symbols to their mouse homologs.
# It requires a global data.frame named 'gene_mapping' with 'rat_gene' and 'mouse_gene' columns.
convert_rat_to_mouse_genes <- function(rat_genes_list) {
  mouse_genes <- sapply(rat_genes_list, function(gene) {
    # Find the corresponding mouse gene in the mapping table
    mapped_gene <- gene_mapping$mouse_gene[gene_mapping$rat_gene == gene]
    if (length(mapped_gene) > 0) {
      return(mapped_gene[1]) # Return the first match if multiple exist
    } else {
      return(NA) # Return NA if no mapping is found
    }
  })
  return(mouse_genes)
}


# ---
# 3. DATA LOADING AND ANNOTATION
# ---

# Load the previously saved integrated Seurat object
RATSCI <- readRDS(sci_integrated_path)

# Ensure proper ordering of samples for plotting
RATSCI$orig.ident <- factor(RATSCI$orig.ident, levels = c('Intact','C2W','C4W','C8W'))

# Define new cell type identities
new_cluster_ids <- c('Microglia','OL','OL','Microglia','Macrophages','Fibroblast',
                     'OL','OPC','Fibroblast','Fibroblast','Neuron',
                     'Macrophages','Microglia','Schwann','Neuron','Astrocyte',
                     'OPC','Neuron','Neuron','T cell','Neuron',
                     'NPCs','Pericyte','Endothelial cells','Ependymal cells','OPC',
                     'OL')
names(new_cluster_ids) <- levels(RATSCI)
RATSCI <- RenameIdents(RATSCI, new_cluster_ids)

# Re-order cell types for consistent plotting
cell_type_levels <- c('OL', 'Microglia', 'Astrocyte', 'OPC', 'Neuron', 'NPCs', 'T cell',
                      'Macrophages', 'Schwann', 'Pericyte', 'Fibroblast', 
                      'Ependymal cells', 'Endothelial cells')
Idents(RATSCI) <- factor(Idents(RATSCI), levels = cell_type_levels)
RATSCI$celltype <- Idents(RATSCI)

# Rescale data in the RNA assay
DefaultAssay(RATSCI) <- 'RNA'
RATSCI <- ScaleData(RATSCI, features = rownames(RATSCI))


# ---
# 4. VISUALIZATION AND PLOTTING
# ---

# Define a consistent color palette for cell types
cell_colors <- c(
  'OL' = "#894FBE" , 'Microglia' = "#BCAC39", 'Astrocyte' = "#51DD24",
  'OPC' = "#C3BAD1", 'Neuron' = "#00ACAE", 'NPCs' = "#D3B5AA",
  'T cell' = "#ABABAB", 'Macrophages' = '#ef7fdf', 'Schwann' = "#889EC0",
  'Pericyte' = "#3B6591", 'Fibroblast' = "#E8948D", 
  'Ependymal cells' = "#d03800", 'Endothelial cells' = "#F2a151"
)

# IMPORTANT: Define your output directory for plots
output_dir <- "path/to/your/output/plots/"
if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }

# Generate and save UMAP plots
umap_plot <- DimPlot(RATSCI, pt.size = 0.1, cols = cell_colors, label = TRUE) + NoLegend()
ggsave(file.path(output_dir, "UMAP_annotated_cells.pdf"), plot = umap_plot, width = 8, height = 8)

umap_split_plot <- DimPlot(RATSCI, pt.size = 1.5, split.by = 'orig.ident', cols = cell_colors, label = FALSE)
ggsave(file.path(output_dir, "UMAP_split_by_timepoint.pdf"), plot = umap_split_plot, width = 16, height = 5)

# Generate and save cell proportion bar plot
cell_prop_data <- as.data.frame(prop.table(table(RATSCI$celltype, RATSCI$orig.ident), margin = 2))
colnames(cell_prop_data) <- c("Celltype", "Group", "Proportion")

prop_plot <- ggplot(cell_prop_data, aes(x = Group, y = Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = cell_colors) +
  labs(y = "Cell Proportion", x = "Timepoint") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(file.path(output_dir, "cell_proportions_bar.svg"), plot = prop_plot, width = 8, height = 5)


# ---
# 5. CELLCHAT PREPARATION: CONVERT RAT GENES TO MOUSE HOMOLOGS
# ---


# Subset the object to the cell types of interest for communication analysis
RATSCI_glial <- subset(RATSCI, idents = c("OL", "Microglia", "Astrocyte", 'OPC', 'Neuron'))

# Get the list of all rat genes in the subset
rat_genes_list <- rownames(RATSCI_glial)

# Use homologene to find mouse homologs (Rat TaxID: 10116, Mouse TaxID: 10090)
homologs <- homologene(rat_genes_list, inTax = 10116, outTax = 10090)
# Create the 'gene_mapping' dataframe required by the conversion function
gene_mapping <- as.data.frame(homologs)
colnames(gene_mapping) <- c("rat_gene", "mouse_gene")

# Convert the original rat genes to mouse gene symbols
mouse_genes_mapped <- convert_rat_to_mouse_genes(rat_genes_list)

# Create a data frame to manage the gene conversion
gene_conversion_df <- data.frame(
  rat_gene = rat_genes_list,
  mouse_gene = mouse_genes_mapped
)

# Filter for genes that were successfully mapped and are unique
gene_conversion_df <- gene_conversion_df %>%
  na.omit() %>%
  distinct(rat_gene, .keep_all = TRUE) %>%
  distinct(mouse_gene, .keep_all = TRUE)

# Subset the Seurat object to only include genes we can map
RATSCI_for_cellchat <- subset(RATSCI_glial, features = gene_conversion_df$rat_gene)

# Rename the genes in the new object to their mouse homologs
RATSCI_mouse_genes <- RenameGenesSeurat(
  obj = RATSCI_for_cellchat,
  newnames = gene_conversion_df$mouse_gene
)


# ---
# 6. CELL-CELL COMMUNICATION ANALYSIS WITH CELLCHAT
# ---

# Streamline the analysis by creating a function to run the CellChat workflow
run_cellchat_analysis <- function(seurat_obj) {
  # Create CellChat object
  cellchat_obj <- createCellChat(
    object = seurat_obj,
    group.by = "celltype",
    assay = "RNA"
  )
  # Set the ligand-receptor interaction database
  cellchat_obj@DB <- CellChatDB.mouse
  # Pre-processing
  cellchat_obj <- subsetData(cellchat_obj)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  # Compute communication probabilities
  cellchat_obj <- computeCommunProb(cellchat_obj, raw.use = TRUE, population.size = TRUE)
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 10)
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)
  # Compute network centrality
  cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name = "netP")
  
  return(cellchat_obj)
}

# Split the Seurat object with mouse genes by timepoint
sp_mouse <- SplitObject(RATSCI_mouse_genes, split.by = 'orig.ident')

# Apply the CellChat workflow to each timepoint
object.list <- lapply(sp_mouse, run_cellchat_analysis)
names(object.list) <- c('Intact', 'C2W', 'C4W', 'C8W')

# Merge the CellChat objects for comparative analysis
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


# ---
# 7. CELLCHAT VISUALIZATION & COMPARATIVE ANALYSIS
# ---

# Define a color palette for the cell types used in CellChat
cellchat_colors <- cell_colors[names(cell_colors) %in% levels(RATSCI_glial)]
cellchat_output_dir <- file.path(output_dir, "CellChat/")
if (!dir.exists(cellchat_output_dir)) { dir.create(cellchat_output_dir) }


# Compare the total number and strength of interactions across timepoints
interaction_plot <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2,3,4))
ggsave(file.path(cellchat_output_dir, "compare_interactions.pdf"), plot = interaction_plot, width = 8, height = 4)

# **[USER REQUESTED BLOCK]** Visualize and export signaling roles
# This identifies which cell types are dominant senders, receivers, influencers, etc.
# ---
# Calculate the total interaction strength for each cell type to make dot sizes comparable across plots
num_links <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight_min_max <- c(min(num_links), max(num_links))

# Create a signaling role scatter plot for each timepoint
plot_list <- list()
for (i in 1:length(object.list)) {
  plot_list[[i]] <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    title = names(object.list)[i],
    weight.MinMax = weight_min_max,
    color.use = cellchat_colors,
    dot.size = 8, 
    label.size = 8, 
    font.size = 15
  ) + NoLegend()
}

# Arrange all scatter plots into a single figure and save
signaling_role_plot <- patchwork::wrap_plots(plots = plot_list, ncol = 4)
ggsave(file.path(cellchat_output_dir, "signaling_role_scatter_plots.pdf"), plot = signaling_role_plot, width = 16, height = 4)

# Export the data underlying the scatter plots to CSV files
for (i in 1:length(object.list)) {
  sample_name <- names(object.list)[i]
  data_to_save <- plot_list[[i]]$data
  write.csv(
    data_to_save,
    file = file.path(cellchat_output_dir, paste0(sample_name, "_signaling_roles.csv")),
    row.names = FALSE
  )
}
# ---

# Differential interaction analysis (e.g., C8W vs Intact)
diff_plot <- netVisual_diffInteraction(cellchat, weight.scale = TRUE, 
                                       comparison = c(1, 4), # Compare Intact (1) and C8W (4)
                                       label.edge = FALSE, color.use = cellchat_colors)
# Note: netVisual functions often plot directly. To save, wrap in pdf().
pdf(file.path(cellchat_output_dir, "diff_interaction_C8W_vs_Intact.pdf"), width=8, height=8)
print(diff_plot)
dev.off()

# Identify and visualize communication patterns
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 4)
river_plot_out <- netAnalysis_river(cellchat, pattern = "outgoing")
ggsave(file.path(cellchat_output_dir, "outgoing_patterns_river.pdf"), plot = river_plot_out, width = 10, height = 8)

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = 4)
river_plot_in <- netAnalysis_river(cellchat, pattern = "incoming")
ggsave(file.path(cellchat_output_dir, "incoming_patterns_river.pdf"), plot = river_plot_in, width = 10, height = 8)

# Signaling role analysis (e.g., for C8W sample)
role_heatmap <- netAnalysis_signalingRole_heatmap(object.list$C8W, pattern = "all")
pdf(file.path(cellchat_output_dir, "signaling_role_heatmap_C8W.pdf"), width=10, height=8)
print(role_heatmap)
dev.off()

# Save the final merged CellChat object
saveRDS(cellchat, file = file.path(output_dir, "cellchat_merged_object.rds"))

# --- END OF SCRIPT ---