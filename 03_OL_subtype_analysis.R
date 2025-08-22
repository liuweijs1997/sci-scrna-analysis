# ==================================================================================
# Analysis of Oligodendrocyte (OL) Subtypes in Rat SCI
#
# This script performs a detailed investigation of oligodendrocyte heterogeneity.
# The workflow includes:
# 1. Loading and annotating OL subtypes from a focused Seurat object.
# 2. Visualizing OL subtypes and their proportions across timepoints.
# 3. Integrating OL subtype annotations back into the main dataset.
# 4. Performing in-depth CellChat analysis between OL subtypes and other cells,
#    followed by functional enrichment of interacting genes.
# 5. Conducting pseudotime trajectory analysis on OL subtypes using Monocle 2.
# ==================================================================================


# ---
# 1. SETUP & LIBRARY LOADING
# ---

# Load required libraries, installing any that are missing
packages <- c("Seurat", "dplyr", "tidyverse", "patchwork", "RColorBrewer", 
              "ggrepel", "ggalluvial", "homologene", "CellChat", "pheatmap",
              "Nebulosa", "monocle", "clusterProfiler", "org.Rn.eg.db")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)) {
    # Install from CRAN or Bioconductor as appropriate
    if (pkg %in% c("clusterProfiler", "org.Rn.eg.db")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# Set future options for parallel processing
options(future.globals.maxSize = 128 * 1024^3) # 128 GB limit
future::plan("multicore", workers = 4)


# ---
# 2. CUSTOM FUNCTIONS
# ---


RenameGenesSeurat <- function(obj, newnames, gene.use=NULL, de.assay="RNA") {
  print("Running RenameGenesSeurat...")
  lassays <- Assays(obj); assay.use <- obj@reductions$pca@assay.used; DefaultAssay(obj) <- de.assay
  if (is.null(gene.use)) { all_genenames <- rownames(obj) } else { all_genenames <- gene.use; obj <- subset(obj, features=gene.use) }
  order_name <- function(v1, v2, ref){ v2 <- make.names(v2, unique=T); df1 <- data.frame(v1, v2); rownames(df1) <- df1$v1; df1 <- df1[ref, ]; return(df1) }
  df1 <- order_name(v1=all_genenames, v2=newnames, ref=rownames(obj)); all_genenames <- df1$v1; newnames <- df1$v2
  if ('SCT' %in% lassays) { if ('SCTModel.list' %in%  slotNames(obj@assays$SCT)) { obj@assays$SCT@SCTModel.list$model1@feature.attributes <- obj@assays$SCT@SCTModel.list$model1@feature.attributes[all_genenames,]; rownames(obj@assays$SCT@SCTModel.list$model1@feature.attributes) <- newnames } }
  change_assay <- function(a1=de.assay, obj, newnames=NULL, all_genenames=NULL){
    RNA <- obj@assays[a1][[1]]
    if (nrow(RNA) == length(newnames)) {
      if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
      if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
      if (length(RNA@var.features)) { df1 <- order_name(v1=all_genenames, v2=newnames, ref=RNA@var.features); RNA@var.features <- df1$v2 }
      if (length(RNA@scale.data)){ df1 <- order_name(v1=all_genenames, v2=newnames, ref=rownames(RNA@scale.data)); rownames(RNA@scale.data) <- df1$v2 }
    } else { "Unequal gene sets: nrow(RNA) != nrow(newnames)" }; obj@assays[a1][[1]] <- RNA; return(obj)
  }
  for (a in lassays) { DefaultAssay(obj) <- a; df1 <- order_name(v1=all_genenames, v2=newnames, ref=rownames(obj)); obj <- change_assay(obj=obj, a1=a, newnames=df1$v2, all_genenames=df1$v1) }
  hvg <- VariableFeatures(obj, assay=assay.use)
  if (length(obj@reductions$pca)){ df1 <- order_name(v1=all_genenames, v2=newnames, ref=hvg); df1 <- df1[rownames(obj@reductions$pca@feature.loadings),]; rownames(obj@reductions$pca@feature.loadings) <- df1$v2 }
  try(obj[[de.assay]]@meta.features <- data.frame(row.names = rownames(obj[[de.assay]]))); return(obj)
}

convert_rat_to_mouse_genes <- function(rat_genes_list) {
  mouse_genes <- sapply(rat_genes_list, function(gene) {
    mapped_gene <- gene_mapping$mouse_gene[gene_mapping$rat_gene == gene]
    if (length(mapped_gene) > 0) { return(mapped_gene[1]) } else { return(NA) }
  })
  return(mouse_genes)
}

erichplotGO <- function(data4plot){
  data4plot <- data4plot[order(data4plot$qvalue, decreasing = FALSE)[1:20], ]
  data4plot$BgRatio <- apply(data4plot, 1, function(x) as.numeric(strsplit(x[4],'/')[[1]][1])) / apply(data4plot, 1, function(x) as.numeric(strsplit(x[5],'/')[[1]][1]))
  p <- ggplot(data4plot, aes(BgRatio, Description)) + geom_point(aes(size=Count, color=-1*log10(qvalue))) +
    scale_colour_gradient(low="#90EE90", high="red") + labs(color=expression(-log[10](qvalue)), size="Gene Count", x="Rich Factor", y="Term Description") + theme_bw()
  return(p)
}
erichplotKEGG <- function(data4plot){
  data4plot <- data4plot[order(data4plot$qvalue, decreasing = FALSE)[1:15], ]
  data4plot$Description <- factor(data4plot$Description, levels = rev(data4plot$Description))
  p <- ggplot(data4plot, aes(x = Description, y = Count, fill = -log10(qvalue))) + geom_bar(stat = "identity") + coord_flip() +
    scale_fill_gradientn(colors = c("blue", "purple", "red")) + labs(fill = "-log10(qvalue)", x = "Term Description", y = "Gene Count") + theme_bw()
  return(p)
}


# ---
# 3. OL SUBTYPE ANNOTATION & VISUALIZATION
# ---

# IMPORTANT: Define your main input and output directories
input_dir <- "path/to/your/input_data/"
output_dir <- "path/to/your/output_plots_and_data/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load the pre-clustered oligodendrocyte Seurat object
RATOL <- readRDS(file.path(input_dir, ""))

# Annotate OL subtypes based on clustering
new.cluster.ids <- c('MFOL','MOL','MFOL','MOL','COL','MFOL','MFOL','NF/DOL','COL','MOL')
names(new.cluster.ids) <- levels(RATOL)
RATOL <- RenameIdents(RATOL, new.cluster.ids)
Idents(RATOL) <- factor(Idents(RATOL), levels = c('COL','MFOL','NF/DOL','MOL'))
RATOL$celltype <- Idents(RATOL)

# Define color palette for OL subtypes
ol_colors <- c('COL'='#5C1CE2', 'MFOL'='#0083DF', 'NF/DOL'='#DF00BC', 'MOL'='#FF6B1F')

# Visualize OL subtypes on UMAP
p_ol_umap <- DimPlot(RATOL, label=TRUE, pt.size=1.5, cols=ol_colors)
ggsave(file.path(output_dir, "UMAP_OL_subtypes.pdf"), plot=p_ol_umap, width=8, height=7)

p_ol_umap_split <- DimPlot(RATOL, split.by='orig.ident', label=TRUE, pt.size=1.5, cols=ol_colors, ncol=2)
ggsave(file.path(output_dir, "UMAP_OL_subtypes_split.pdf"), plot=p_ol_umap_split, width=12, height=10)

# Plot cell proportions
RATOL$orig.ident <- factor(RATOL$orig.ident, levels = c('C8W','C4W','C2W','Intact'))
cell_prop_data <- as.data.frame(prop.table(table(RATOL$celltype, RATOL$orig.ident), 2))
colnames(cell_prop_data) <- c("Celltype", "Group", "Proportion")

p_prop <- ggplot(cell_prop_data, aes(Group, Proportion, fill=Celltype)) +
  geom_bar(stat="identity", position="fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=ol_colors) + coord_flip() +
  labs(y="Cell Proportion", x="") + theme_bw()
ggsave(file.path(output_dir, "OL_subtype_proportions.svg"), plot=p_prop, width=8, height=5)


# ---
# 4. INTEGRATE OL SUBTYPES BACK INTO THE MAIN OBJECT
# ---

# Load the main integrated object containing all cell types
RATSCI <- readRDS(file.path(input_dir, "sci_integrated_annotated.rds")) # Assuming this is the main object name

# Transfer the detailed OL subtype annotations to the main object
# This replaces the general 'OL' identity with the specific subtypes
sce.all <- RATSCI
sce.all <- RenameIdents(sce.all, new.cluster.ids) # First, rename all clusters based on the main annotation
Idents(sce.all, cells = colnames(RATOL)) <- Idents(RATOL) # Then, overwrite OL cells with subtype annotations

# Visualize the main UMAP with the new, detailed OL annotations
p_all_umap_with_ol <- DimPlot(sce.all, label=TRUE, repel=TRUE)
ggsave(file.path(output_dir, "UMAP_all_cells_with_OL_subtypes.pdf"), plot=p_all_umap_with_ol, width=10, height=9)


# ---
# 5. CELLCHAT ANALYSIS FOR OL SUBTYPES & OTHER GLIA/NEURONS
# ---

## 5A. Preparation: Gene Conversion ##
# Subset to the cell types of interest for communication analysis
celltypes_for_chat <- c("COL", 'MFOL', 'NF/DOL', 'MOL', "Microglia", "Astrocyte", 'OPC', 'Neuron')
RATSCI_for_chat <- subset(sce.all, idents = celltypes_for_chat)

# Convert rat genes to mouse homologs
rat_genes_list <- rownames(RATSCI_for_chat)
homologs <- homologene(rat_genes_list, inTax = 10116, outTax = 10090)
gene_mapping <- as.data.frame(homologs); colnames(gene_mapping) <- c("rat_gene", "mouse_gene")
mouse_genes_mapped <- convert_rat_to_mouse_genes(rat_genes_list)
gene_conversion_df <- data.frame(rat_gene = rat_genes_list, mouse_gene = mouse_genes_mapped) %>%
  na.omit() %>% distinct(rat_gene, .keep_all = TRUE) %>% distinct(mouse_gene, .keep_all = TRUE)

RATSCI_for_chat_subset <- subset(RATSCI_for_chat, features = gene_conversion_df$rat_gene)
RATSCI_mouse_genes <- RenameGenesSeurat(obj = RATSCI_for_chat_subset, newnames = gene_conversion_df$mouse_gene)

## 5B. Streamlined CellChat and Enrichment Analysis ##
sp_mouse <- SplitObject(RATSCI_mouse_genes, split.by = 'orig.ident')

# Function to run enrichment analysis on interaction data
run_enrichment <- function(interaction_df, gene_col, file_prefix, output_dir) {
  genes <- unique(interaction_df[[gene_col]])
  gene_df <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Rn.eg.db)
  
  if (nrow(gene_df) > 0) {
    # KEGG Enrichment
    ekegg <- enrichKEGG(unique(gene_df$ENTREZID), organism='rno', pvalueCutoff=0.05)
    if (!is.null(ekegg) && nrow(ekegg@result) > 0) {
      ekegg <- setReadable(ekegg, 'org.Rn.eg.db', 'ENTREZID')
      p_kegg <- erichplotKEGG(ekegg@result)
      ggsave(file.path(output_dir, paste0(file_prefix, "_KEGG.pdf")), plot=p_kegg, width=10, height=8)
      write.csv(ekegg@result, file.path(output_dir, paste0(file_prefix, "_KEGG_results.csv")))
    }
    # GO Enrichment
    go <- enrichGO(gene_df$SYMBOL, org.Rn.eg.db, keyType="SYMBOL", ont="all", pvalueCutoff=0.05)
    if (!is.null(go) && nrow(go@result) > 0) {
      p_go <- erichplotGO(go@result)
      ggsave(file.path(output_dir, paste0(file_prefix, "_GO.pdf")), plot=p_go, width=10, height=8)
      write.csv(go@result, file.path(output_dir, paste0(file_prefix, "_GO_results.csv")))
    }
  }
}

# Loop through each sample to run CellChat and Enrichment
cellchat_output_dir <- file.path(output_dir, "CellChat_OL_Subtypes")
if (!dir.exists(cellchat_output_dir)) dir.create(cellchat_output_dir)

object.list <- list()
for (sample_name in names(sp_mouse)) {
  message(paste("--- Running CellChat for", sample_name, "---"))
  
  # Run CellChat
  seurat_obj <- sp_mouse[[sample_name]]
  cellchat_obj <- createCellChat(object = seurat_obj, group.by = "celltype", assay = "RNA")
  cellchat_obj@DB <- CellChatDB.mouse
  cellchat_obj <- subsetData(cellchat_obj)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  cellchat_obj <- computeCommunProb(cellchat_obj, raw.use = TRUE, population.size = TRUE)
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 10)
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)
  object.list[[sample_name]] <- cellchat_obj
  
  # Subset interactions involving COL and run enrichment
  interactions_to_col <- subsetCommunication(cellchat_obj, targets.use = "COL")
  interactions_from_col <- subsetCommunication(cellchat_obj, sources.use = "COL")
  
  write.csv(interactions_to_col, file.path(cellchat_output_dir, paste0(sample_name, "_all_to_COL.csv")))
  write.csv(interactions_from_col, file.path(cellchat_output_dir, paste0(sample_name, "_COL_to_all.csv")))
  
  # Run enrichment for ligands and receptors
  run_enrichment(interactions_from_col, "ligand", paste0(sample_name, "_COL_ligands"), cellchat_output_dir)
  run_enrichment(interactions_to_col, "receptor", paste0(sample_name, "_COL_receptors"), cellchat_output_dir)
}

# Merge CellChat objects and perform comparative analysis
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


# ---
# 6. PSEUDOTIME TRAJECTORY ANALYSIS WITH MONOCLE 2
# ---

trajectory_output_dir <- file.path(output_dir, "Trajectory_Analysis")
if (!dir.exists(trajectory_output_dir)) dir.create(trajectory_output_dir)

# Subset to only OL subtypes for trajectory analysis
OL_for_traj <- subset(sce.all, idents = c("COL", 'MFOL', 'NF/DOL', 'MOL'))

# Create Monocle object from Seurat object
data <- as(as.matrix(OL_for_traj@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = OL_for_traj@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

# Pre-process the data
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds)

# Filter genes and identify ordering genes
mycds <- detectGenes(mycds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(mycds), num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(mycds[expressed_genes,], fullModelFormulaStr = "~celltype")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))

mycds <- setOrderingFilter(mycds, ordering_genes)

# Reduce dimensionality and order cells
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)

# Plot the trajectory colored by cell type, state, and pseudotime
p_traj_celltype <- plot_cell_trajectory(mycds, color_by = "celltype", cell_size = 0.5) + scale_color_manual(values = ol_colors)
p_traj_state <- plot_cell_trajectory(mycds, color_by = "State", cell_size = 0.5)
p_traj_pseudo <- plot_cell_trajectory(mycds, color_by = "Pseudotime", cell_size = 0.5)

combined_traj_plot <- p_traj_celltype | p_traj_state | p_traj_pseudo
ggsave(file.path(trajectory_output_dir, "OL_trajectory_plots.pdf"), plot=combined_traj_plot, width=21, height=6)

# Heatmap of gene expression changes along pseudotime
p_pseudo_heatmap <- plot_pseudotime_heatmap(mycds[ordering_genes[1:100],], num_clusters = 4, show_rownames = FALSE)
ggsave(file.path(trajectory_output_dir, "pseudotime_heatmap.pdf"), plot=p_pseudo_heatmap, width=8, height=10)


# --- END OF SCRIPT ---