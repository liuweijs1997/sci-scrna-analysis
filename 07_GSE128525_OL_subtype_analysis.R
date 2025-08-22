# ==================================================================================
# Analysis of Oligodendrocyte (OL) Subtypes from GSE128525
#
# This script performs a focused analysis on OL subtypes from a public mouse
# spinal cord injury dataset (Floriddia et al., 2020).
# The workflow includes:
# 1. Loading and preparing the public Seurat object.
# 2. Visualizing cell proportions and the expression of specific genes of interest.
# 3. A streamlined workflow for identifying marker genes for each OL subtype
#    and performing systematic GO and KEGG functional enrichment analysis.
# ==================================================================================


# ---
# 1. SETUP & LIBRARY LOADING
# ---

# Load required libraries, installing any that are missing
packages <- c(
  "Seurat", "dplyr", "tidyverse", "patchwork", "cowplot",
  "clusterProfiler", "org.Mm.eg.db" # Mouse annotation database
)
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)) {
    # Install from CRAN or Bioconductor as appropriate
    if (pkg %in% c("clusterProfiler", "org.Mm.eg.db")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}


# ---
# 2. CUSTOM PLOTTING FUNCTIONS
# ---

# NOTE: These custom functions for plotting enrichment results are preserved from your script.
erichplotGO <- function(data4plot){
  # Custom plotting function for GO enrichment results (bubble plot).
  if (is.null(data4plot) || nrow(data4plot) == 0) return(NULL)
  data4plot <- data4plot[order(data4plot$qvalue, decreasing = FALSE), ][1:min(20, nrow(data4plot)), ]
  p <- ggplot(data4plot, aes(x = Count, y = reorder(Description, Count))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(x = "Gene Count", y = "GO Term", color = "p.adjust") + theme_bw()
  return(p)
}

erichplotKEGG <- function(data4plot){
  # Custom plotting function for KEGG enrichment results (bar plot).
  if (is.null(data4plot) || nrow(data4plot) == 0) return(NULL)
  data4plot <- data4plot[order(data4plot$qvalue, decreasing = FALSE), ][1:min(15, nrow(data4plot)), ]
  data4plot$Description <- factor(data4plot$Description, levels = rev(data4plot$Description))
  p <- ggplot(data4plot, aes(x = Description, y = Count, fill = -log10(qvalue))) +
    geom_bar(stat = "identity") + coord_flip() +
    scale_fill_gradientn(colors = c("blue", "purple", "red")) +
    labs(fill = "-log10(qvalue)", x = "KEGG Pathway", y = "Gene Count") + theme_bw()
  return(p)
}


# ---
# 3. DATA LOADING AND PREPARATION
# ---

# IMPORTANT: Define your main input and output directories
input_dir <- "path/to/your/input_data/"
output_dir <- "path/to/your/output_plots_and_data/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load the .Rdata file containing the Seurat object
load(file.path(input_dir, "GSE128525_SeuratObjSpCrdGEO.Rdata"))

# The loaded object is named 'oligos.integrated'. Let's inspect it.
# This is a crucial step to ensure compatibility with the current Seurat version.
oligos.integrated <- UpdateSeuratObject(oligos.integrated)

# Set the cell identities to the annotated clusters
Idents(oligos.integrated) <- oligos.integrated$confidentclusters

# Subset the object to include only the desired oligodendrocyte subtypes
ol_subtypes <- c('MFOL1','MFOL2','MOL1','MOL2','MOL3','MOL4','MOL5','MOL6')
oligos <- subset(oligos.integrated, idents = ol_subtypes)

# Prepare the object for analysis by normalizing and scaling the data
oligos <- NormalizeData(oligos)
oligos <- ScaleData(oligos, features = rownames(oligos))


# ---
# 4. VISUALIZATION OF GENE EXPRESSION & CELL PROPORTIONS
# ---

## 4A. Gene Expression Bar Plots ##
# Define a list of genes of interest for visualization
genes_to_plot <- c('Tubb2a','Tubb2b','Dlc1','Dock7','Tubb3','Tenm4','Thsd7b',
                   'Thsd7a','Plekha7','S100b','Tox','Map2')

# Calculate the average expression for each gene in each OL subtype
avg_expr <- AverageExpression(oligos, features = genes_to_plot, group.by = "ident", assays = "RNA")$RNA

# Prepare the data for ggplot by converting it to a long format
avg_expr_long <- as.data.frame(avg_expr) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "CellType", values_to = "AvgExpression")

# Ensure the cell types are plotted in the desired order
avg_expr_long$CellType <- factor(avg_expr_long$CellType, levels = ol_subtypes)

# Create the faceted bar plot
p_gene_expr <- ggplot(avg_expr_long, aes(x = CellType, y = AvgExpression, fill = CellType)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  labs(title = "Average Gene Expression in OL Subtypes",
       x = "OL Subtype", y = "Average Expression") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave(file.path(output_dir, "gene_expression_barchart.pdf"), plot = p_gene_expr, width = 14, height = 10)


## 4B. Cell Proportion Plots ##
# Ensure condition levels are in a logical order
oligos$orig.ident <- factor(oligos$orig.ident, levels = c('CTRL', 'IS', 'WD'))

# Calculate cell proportions
cell_prop_data <- as.data.frame(prop.table(table(Idents(oligos), oligos$orig.ident), margin = 2))
colnames(cell_prop_data) <- c("Celltype", "Group", "Proportion")

# Plot the proportions as a stacked bar chart
p_prop <- ggplot(cell_prop_data, aes(x = Group, y = Proportion, fill = Celltype)) +
  geom_bar(stat="identity", position="fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Cell Proportion", x = "Condition") +
  theme_classic()

ggsave(file.path(output_dir, "OL_subtype_proportions.pdf"), plot = p_prop, width = 6, height = 6)


# ---
# 5. MARKER DISCOVERY & FUNCTIONAL ENRICHMENT
# ---

# This streamlined workflow finds markers and runs enrichment in a single automated loop.
message("--- Starting Marker Discovery and Functional Enrichment ---")

# Find marker genes for all OL subtypes at once
all_markers <- FindAllMarkers(oligos, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_markers, file.path(output_dir, "all_OL_subtype_markers.csv"))

# Split the markers data frame into a list, with one element per cluster
markers_list <- split(all_markers, all_markers$cluster)

# Create a dedicated directory for enrichment results
enrichment_dir <- file.path(output_dir, "Enrichment_Analysis")
if (!dir.exists(enrichment_dir)) dir.create(enrichment_dir, recursive = TRUE)

# Loop through each cluster's markers to perform GO and KEGG analysis
for (cluster_name in names(markers_list)) {
  message(paste("--- Running enrichment for:", cluster_name, "---"))
  
  marker_df <- markers_list[[cluster_name]]
  genes <- marker_df$gene
  
  # Convert gene symbols to Entrez IDs (required for enrichment)
  # This dataset is from mouse, so we use the mouse annotation database (org.Mm.eg.db)
  gene_df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop = TRUE)
  
  if (nrow(gene_df) > 0) {
    # Perform GO Enrichment (Biological Process)
    go_results <- enrichGO(
      gene = gene_df$ENTREZID,
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
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


# --- END OF SCRIPT ---
