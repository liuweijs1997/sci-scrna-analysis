# ==================================================================================
# Zebrafish SCI: Analysis of Oligodendrocyte (OL) Subtypes
#
# This script performs a detailed investigation of oligodendrocyte heterogeneity.
# The workflow includes:
# 1. Loading and annotating OL subtypes from a focused Seurat object.
# 2. Visualizing OL subtypes and their proportions across timepoints.
# 3. Performing functional enrichment analysis on subtype markers.
# 4. Integrating OL subtype annotations back into the main dataset.
# 5. A streamlined workflow for running CellChat on OL subtypes and Neurons.
# 6. Conducting pseudotime trajectory analysis on the OL lineage using Monocle 2.
# ==================================================================================


# ---
# 1. SETUP & LIBRARY LOADING
# ---

# Load required libraries, installing any that are missing
packages <- c(
  "Seurat", "dplyr", "tidyverse", "patchwork", "homologene", "CellChat",
  "clusterProfiler", "org.Rn.eg.db", "monocle", "gground", "ggprism"
)
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)) {
    # Install from CRAN or Bioconductor as appropriate
    if (pkg %in% c("clusterProfiler", "org.Rn.eg.db", "monocle")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# Set options for parallel processing
options(future.globals.maxSize = 128 * 1024^3)
future::plan("multicore", workers = 8) # Adjust based on your machine's cores


# ---
# 2. CUSTOM FUNCTIONS
# ---

RenameGenesSeurat <- function(obj, newnames, gene.use=NULL, de.assay="RNA") {
  # This function renames genes across multiple slots in a Seurat object.
  print("Running RenameGenesSeurat...")
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
    } else { message("Unequal gene sets.") }; obj@assays[[a1]] <- RNA; return(obj)
  }
  for (a in lassays) { df1 <- order_name(v1=all_genenames, v2=newnames, ref=rownames(obj@assays[[a]])); obj <- change_assay(obj=obj, a1=a, newnames=df1$v2, all_genenames=df1$v1) }
  if (!is.null(obj@reductions$pca)){ df1 <- order_name(v1=all_genenames, v2=newnames, ref=rownames(obj@reductions$pca@feature.loadings)); rownames(obj@reductions$pca@feature.loadings) <- df1$v2 }
  try(obj[[de.assay]]@meta.features <- data.frame(row.names = rownames(obj[[de.assay]]))); return(obj)
}

convert_zebrafish_to_mouse_genes <- function(zebrafish_genes_list) {
  # Converts a list of zebrafish gene symbols to mouse homologs.
  # Relies on a global 'gene_mapping' dataframe.
  mouse_genes <- sapply(zebrafish_genes_list, function(gene) {
    mapped_gene <- gene_mapping$mouse_gene[gene_mapping$zebrafish_gene == gene]
    if (length(mapped_gene) > 0) { return(mapped_gene[1]) } else { return(NA) }
  })
  return(mouse_genes)
}

erichplotGO <- function(data4plot){
  # Custom plotting function for GO enrichment results.
  data4plot <- data4plot[order(data4plot$qvalue, decreasing = FALSE), ][1:min(20, nrow(data4plot)), ]
  if (nrow(data4plot) == 0) return(NULL)
  p <- ggplot(data4plot, aes(x = Count, y = reorder(Description, Count))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(x = "Gene Count", y = "GO Term") + theme_bw()
  return(p)
}


# ---
# 3. ZEBRAFISH OL SUBTYPE ANNOTATION & VISUALIZATION
# ---

# IMPORTANT: Define your main input and output directories
input_dir <- "path/to/your/input_data/"
output_dir <- "path/to/your/output_plots_and_data/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load the pre-clustered zebrafish oligodendrocyte Seurat object
DROL <- readRDS(file.path(input_dir, ".rds"))

# Annotate OL subtypes based on clustering and set factor levels for plotting
new.cluster.ids <- c('MFOL','MFOL','MOL','DOL','MOL','MFOL','MOL','MFOL','NFOL','MOL','COL','COL','COL','COL','COL','COL')
names(new.cluster.ids) <- levels(DROL)
DROL <- RenameIdents(DROL, new.cluster.ids)
Idents(DROL) <- factor(Idents(DROL), levels = c('COL','MFOL','DOL','MOL','NFOL'))
DROL$celltype <- Idents(DROL)
DROL$group <- DROL$orig.ident

# Define a color palette for the OL subtypes
ol_colors <- c('COL'='#5C1CE2', 'MFOL'='#0083DF', 'DOL'='#DF00BC', 'MOL'='#FF6B1F', 'NFOL'='#00ACAE')

# Visualize and save UMAPs
p_umap <- DimPlot(DROL, label = TRUE, pt.size = 1.2, cols = ol_colors)
ggsave(file.path(output_dir, "UMAP_DROL_subtypes.pdf"), plot = p_umap, width = 8, height = 7)

p_umap_split <- DimPlot(DROL, split.by = 'orig.ident', label = TRUE, pt.size = 1.5, cols = ol_colors)
ggsave(file.path(output_dir, "UMAP_DROL_subtypes_split.pdf"), plot = p_umap_split, width = 12, height = 5)


# ---
# 4. FUNCTIONAL ENRICHMENT OF SUBTYPE MARKERS
# ---

# Find markers for each OL subtype
all_ol_markers <- FindAllMarkers(object = DROL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Create a dedicated directory for enrichment results
enrichment_dir <- file.path(output_dir, "Enrichment_Analysis")
if (!dir.exists(enrichment_dir)) dir.create(enrichment_dir, recursive = TRUE)

# Loop through each cluster to perform GO enrichment
clusters <- unique(all_ol_markers$cluster)
for (cluster_name in clusters) {
  message(paste("--- Running enrichment for OL subtype:", cluster_name, "---"))
  
  # Get marker genes for the current cluster
  marker_df <- all_ol_markers %>% filter(cluster == cluster_name)
  genes <- marker_df$gene
  
  gene_df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db, drop = TRUE)
  
  if (nrow(gene_df) > 0) {
    go_results <- enrichGO(
      gene = gene_df$ENTREZID,
      OrgDb = org.Rn.eg.db,
      ont = "BP", # Biological Process
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    
    if (!is.null(go_results) && nrow(go_results) > 0) {
      p_go <- erichplotGO(as.data.frame(go_results))
      if (!is.null(p_go)) {
        ggsave(file.path(enrichment_dir, paste0("GO_BP_", cluster_name, ".pdf")), plot = p_go, width = 10, height = 8)
      }
      write.csv(as.data.frame(go_results), file.path(enrichment_dir, paste0("GO_BP_", cluster_name, "_results.csv")))
    }
  }
}


# ---
# 5. CELLCHAT ANALYSIS: OL SUBTYPES & NEURONS
# ---

## 5A. Preparation: Integrate annotations and convert genes ##
# Load the main object and transfer the new, detailed OL subtype annotations
DRSCI <- readRDS(file.path(input_dir, "integrated_annotated.rds"))
Idents(DRSCI, cells = colnames(DROL)) <- Idents(DROL)

# Subset to the cell types of interest for communication analysis
celltypes_for_chat <- c('COL','MFOL','DOL','NFOL','MOL','Neuron')
DRSCI_for_chat <- subset(DRSCI, idents = celltypes_for_chat)
DRSCI_for_chat$celltype <- Idents(DRSCI_for_chat) # Ensure celltype metadata is updated

# Helper function to automate gene conversion and CellChat setup
prepare_and_run_cellchat <- function(seurat_obj_zebrafish, sample_name) {
  message(paste("--- Preparing and running CellChat for:", sample_name, "---"))
  
  # Convert zebrafish genes to mouse homologs
  zebrafish_genes <- rownames(seurat_obj_zebrafish)
  homologs <- homologene(zebrafish_genes, inTax = 7955, outTax = 10090) # Zebrafish -> Mouse
  gene_mapping <<- as.data.frame(homologs); colnames(gene_mapping) <<- c("zebrafish_gene", "mouse_gene")
  mouse_genes_mapped <- convert_zebrafish_to_mouse_genes(zebrafish_genes)
  gene_conversion_df <- data.frame(z_gene = zebrafish_genes, m_gene = mouse_genes_mapped) %>%
    na.omit() %>% distinct(z_gene, .keep_all = TRUE) %>% distinct(m_gene, .keep_all = TRUE)
  
  seurat_subset <- subset(seurat_obj_zebrafish, features = gene_conversion_df$z_gene)
  seurat_obj_mouse <- RenameGenesSeurat(obj = seurat_subset, newnames = gene_conversion_df$m_gene)
  seurat_obj_mouse <- ScaleData(seurat_obj_mouse)
  
  # Run CellChat pipeline
  cellchat_obj <- createCellChat(object = seurat_obj_mouse, group.by = "celltype", assay = "RNA")
  cellchat_obj@DB <- CellChatDB.mouse
  cellchat_obj <- subsetData(cellchat_obj)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  cellchat_obj <- projectData(cellchat_obj, PPI.mouse)
  cellchat_obj <- computeCommunProb(cellchat_obj, raw.use = TRUE, population.size = TRUE)
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 10)
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)
  
  return(cellchat_obj)
}

## 5B. Execute CellChat Analysis ##
# Split object by timepoint and apply the workflow
sp_zebrafish <- SplitObject(DRSCI_for_chat, split.by = 'orig.ident')
object.list <- lapply(names(sp_zebrafish), function(name) {
  prepare_and_run_cellchat(sp_zebrafish[[name]], name)
})
names(object.list) <- names(sp_zebrafish)

# Merge for comparative analysis
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


# ---
# 6. PSEUDOTIME TRAJECTORY ANALYSIS (MONOCLE 2)
# ---

trajectory_dir <- file.path(output_dir, "Trajectory_Analysis")
if (!dir.exists(trajectory_dir)) dir.create(trajectory_dir)

# Create Monocle object from the DROL Seurat object
data <- as(as.matrix(DROL@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = DROL@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

# Pre-process, find ordering genes, reduce dimensions, and order cells
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds)
mycds <- detectGenes(mycds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(mycds), num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(mycds[expressed_genes,], fullModelFormulaStr = "~celltype")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.05))

mycds <- setOrderingFilter(mycds, ordering_genes)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree', residualModelFormulaStr = '~orig.ident')
mycds <- orderCells(mycds, root_state = which.max(table(pData(mycds)$State[pData(mycds)$celltype == "MFOL"])))

# Visualize trajectory plots
p_traj_celltype <- plot_cell_trajectory(mycds, color_by = "celltype", cell_size = 0.7) + scale_color_manual(values = ol_colors)
p_traj_pseudotime <- plot_cell_trajectory(mycds, color_by = "Pseudotime", cell_size = 0.7)
p_traj_state <- plot_cell_trajectory(mycds, color_by = "State", cell_size = 0.7)

p_traj_combined <- p_traj_celltype | p_traj_pseudotime | p_traj_state
ggsave(file.path(trajectory_dir, "OL_trajectory_combined.pdf"), plot = p_traj_combined, width = 21, height = 6)


# --- END OF SCRIPT ---
