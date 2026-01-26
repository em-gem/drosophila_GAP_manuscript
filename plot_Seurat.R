# Load required libraries
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)

# Set seed for reproducibility
set.seed(42)

# Define custom color palette
custom_colors <- c(
  '#E5A712',  # orange
  '#4F97C9',  # sky_blue
  '#007887',  # bluish_green
  '#F2DE61',  # yellow
  '#D37719',  # vermillion
  '#CC8D93',  # reddish_purple
  '#529098',  # blueish_green (saved for last)
  '#A7A9AC'   # gray (saved for last)
)

# Set paths and output directory

# Path to your saved Seurat object
seurat_object_path <- "PATH/TO/SEURAT/OBJECT"

# Output directory for figures
output_dir <- "PATH/TO/OUTPUT/DIRECTORY"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set working directory to output folder
setwd(output_dir)

# Load Seurat object
cat("Loading Seurat object...\n")
merged_seurat <- readRDS(seurat_object_path)

cat("Seurat object loaded successfully!\n")
cat("Number of cells:", ncol(merged_seurat), "\n")
cat("Number of genes:", nrow(merged_seurat), "\n")
cat("Number of clusters:", length(unique(merged_seurat$seurat_clusters)), "\n")

# Find marker genes
marker_genes <- c('slbo', 'fru', 'LifeActGFP', 'dsRedExpress')

# Join layers if needed
if (length(Layers(merged_seurat)) > 1) {
  cat("Joining data layers...\n")
  merged_seurat <- JoinLayers(merged_seurat)
}

# Find all markers
cluster_markers <- FindAllMarkers(merged_seurat, 
                                  only.pos = TRUE,
                                  verbose = FALSE)

cat("Marker genes found:", nrow(cluster_markers), "markers across", 
    length(unique(cluster_markers$cluster)), "clusters\n")


# 1. Heatmap of top 25 genes per cluster
cat("\nGenerating heatmap of top 25 marker genes...\n")

# Get top 25 markers per cluster
top25_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

# Get scaled data - works with both Assay and Assay5
scaled_data <- LayerData(merged_seurat, assay = "RNA", layer = "scale.data")
genes_in_scale_data <- rownames(scaled_data)

# Filter to only genes that are in scaled data
top25_markers_filtered <- top25_markers %>%
  filter(gene %in% genes_in_scale_data)

cat("Using", nrow(top25_markers_filtered), "genes present in scaled data (out of", 
    nrow(top25_markers), "total)\n")

# If no genes in scale.data, scale the top markers now
if (nrow(top25_markers_filtered) == 0) {
  cat("No genes in scale.data - scaling top markers now...\n")
  merged_seurat <- ScaleData(merged_seurat, features = unique(top25_markers$gene))
  scaled_data <- LayerData(merged_seurat, assay = "RNA", layer = "scale.data")
  genes_in_scale_data <- rownames(scaled_data)
  top25_markers_filtered <- top25_markers %>%
    filter(gene %in% genes_in_scale_data)
  cat("Now using", nrow(top25_markers_filtered), "genes\n")
}

# Create heatmap with fixed aspect ratio - calculate height based on number of genes
if (nrow(top25_markers_filtered) > 0) {
  # Calculate appropriate height: aim for ~30-40 genes per inch
  n_genes <- length(unique(top25_markers_filtered$gene))
  heatmap_height <- max(9, n_genes / 35)  # Minimum 9 inches, scale with gene count
  
  heatmap_plot <- DoHeatmap(merged_seurat, 
                            features = top25_markers_filtered$gene,
                            group.by = "seurat_clusters",
                            size = 2.5,
                            angle = 90,
                            draw.lines = TRUE) +  # Add lines between clusters
    scale_fill_gradientn(colors = c("#998ec3", "white", "#f1a340")) +
    ggtitle("Top 25 Marker Genes per Cluster") +
    theme(text = element_text(size = 8),
          axis.text = element_text(size = 6),
          axis.text.y = element_text(size = 5))  # Smaller gene names
  

  # Save as PDF with calculated height
  ggsave("cluster_heatmap_top25.pdf", 
         plot = heatmap_plot,
         width = 6.5, 
         height = heatmap_height,
         units = "in",
         dpi = 300,
         device = "pdf")
  
  # Backup as high-res PNG
  ggsave("cluster_heatmap_top25.png", 
         plot = heatmap_plot,
         width = 6.5, 
         height = heatmap_height,
         units = "in",
         dpi = 600)
  
  cat("Heatmap saved! (height:", heatmap_height, "inches)\n")
} else {
  cat("WARNING: Could not generate heatmap - no genes available\n")
}

# 2. Check which marker genes are present
cat("\nChecking for marker genes in dataset...\n")

present_genes <- marker_genes[marker_genes %in% rownames(merged_seurat)]
missing_genes <- marker_genes[!marker_genes %in% rownames(merged_seurat)]

cat("Present genes:", paste(present_genes, collapse = ", "), "\n")
if (length(missing_genes) > 0) {
  cat("WARNING - Missing genes:", paste(missing_genes, collapse = ", "), "\n")
}


# 3. Dot plot of marker genes
cat("\nGenerating dot plot...\n")

if (length(present_genes) > 0) {
  dotplot <- DotPlot(merged_seurat, 
                     features = present_genes,
                     group.by = "seurat_clusters",
                     dot.scale = 6,
                     cols=c("#A7A9AC", "#f1a340")) +
    coord_flip() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8, face = "italic"),
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          plot.title = element_text(size = 8)) +
    labs(x = "Gene", y = "Cluster") +
    ggtitle("Marker Gene Expression by Cluster")
  

  ggsave("marker_gene_dotplot.pdf", 
         plot = dotplot,
         width = 6, 
         height = 1.5,
         units = "in",
         dpi = 300,
         device = "pdf")
  
  ggsave("marker_gene_dotplot.png", 
         plot = dotplot,
         width = 6, 
         height = 1.5,
         units = "in",
         dpi = 600)
  
  cat("Dot plot saved!\n")
} else {
  cat("WARNING: No marker genes found in dataset, skipping dot plot\n")
}

# 4. UMAP feature plots
cat("\nGenerating feature plots...\n")

if (length(present_genes) > 0) {
  # Calculate grid dimensions
  n_genes <- length(present_genes)
  n_cols <- 3
  n_rows <- ceiling(n_genes / n_cols)
  
  # Make plots square by setting height = width * (n_rows / n_cols)
  plot_width <- 6.5
  plot_height <- plot_width * (n_rows / n_cols)
  
  feature_plot <- FeaturePlot(merged_seurat, 
                              features = present_genes,
                              ncol = n_cols,
                              pt.size = 0.1,
                              order = TRUE,
                              cols=c("gray", "#f1a340")) &
    theme_classic() &
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 8),
          plot.title = element_text(size = 8, face = "italic"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          aspect.ratio = 1)  # Force square aspect ratio for each subplot
  
  ggsave("marker_gene_featureplots.pdf", 
         plot = feature_plot,
         width = plot_width, 
         height = plot_height,
         units = "in",
         dpi = 300,
         device = "pdf")
  
  ggsave("marker_gene_featureplots.png", 
         plot = feature_plot,
         width = plot_width, 
         height = plot_height,
         units = "in",
         dpi = 600)
  
  cat("Feature plots saved! (", plot_width, "x", plot_height, "inches)\n")
}


# 5. Identify marker-defined border and polar cells
cat("\n", rep("=", 70), "\n", sep = "")
cat("MARKER-DEFINED CELL TYPE IDENTIFICATION\n")
cat(rep("=", 70), "\n\n", sep = "")

# Check required genes are present
border_genes <- c("LifeActGFP", "slbo", "fru")
polar_genes <- c("dsRedExpress")

border_genes_present <- border_genes[border_genes %in% rownames(merged_seurat)]
polar_genes_present <- polar_genes[polar_genes %in% rownames(merged_seurat)]

cat("Border cell markers available:", paste(border_genes_present, collapse = ", "), "\n")
cat("Polar cell markers available:", paste(polar_genes_present, collapse = ", "), "\n\n")

# Initialize marker-defined cell types
merged_seurat$marker_defined_type <- "Other"

# Identify POLAR CELLS: dsRedExpress+
if ("dsRedExpress" %in% polar_genes_present) {
  polar_cells_marker <- WhichCells(merged_seurat, expression = dsRedExpress > 0)
  merged_seurat$marker_defined_type[polar_cells_marker] <- "Polar_Cell"
  cat("Polar cells (dsRedExpress+):", length(polar_cells_marker), "\n")
} else {
  polar_cells_marker <- c()
}

# Identify BORDER CELLS: LifeActGFP+ AND slbo+ AND fru+
if (all(c("LifeActGFP", "slbo", "fru") %in% border_genes_present)) {
  border_cells_marker <- WhichCells(merged_seurat, 
                                    expression = LifeActGFP > 0 & slbo > 0 & fru > 0)
  
  # Exclude cells already classified as polar (in case of overlap)
  border_cells_marker <- setdiff(border_cells_marker, polar_cells_marker)
  
  merged_seurat$marker_defined_type[border_cells_marker] <- "Border_Cell"
  cat("Border cells (LifeActGFP+ & slbo+ & fru+):", length(border_cells_marker), "\n")
} else {
  border_cells_marker <- c()
  cat("WARNING: Not all border cell markers available\n")
}

# Check for overlap
if (length(polar_cells_marker) > 0 && length(border_cells_marker) > 0) {
  overlap <- intersect(
    WhichCells(merged_seurat, expression = dsRedExpress > 0),
    WhichCells(merged_seurat, expression = LifeActGFP > 0 & slbo > 0 & fru > 0)
  )
  if (length(overlap) > 0) {
    cat("WARNING:", length(overlap), "cells express both border and polar markers!\n")
  }
}

cat("\nMarker-defined cell type distribution:\n")
print(table(merged_seurat$marker_defined_type))

# 6. Analyze marker-defined border and polar cells
cat("\n", rep("=", 70), "\n", sep = "")
cat("MARKER-DEFINED CELL ANALYSIS\n")
cat(rep("=", 70), "\n\n", sep = "")

# Combine border and polar cells for analysis
marker_defined_cells <- c(border_cells_marker, polar_cells_marker)

if (length(marker_defined_cells) > 20) {  # Need at least some cells
  
  # Subset to marker-defined cells
  marker_subset <- subset(merged_seurat, cells = marker_defined_cells)
  
  cat("Subset created with", ncol(marker_subset), "cells\n")
  cat("  Border cells:", length(border_cells_marker), "\n")
  cat("  Polar cells:", length(polar_cells_marker), "\n\n")
  
  # Find markers comparing border vs polar vs other
  cat("Finding differentially expressed genes...\n")
  Idents(marker_subset) <- "marker_defined_type"
  
  marker_defined_markers <- FindAllMarkers(marker_subset, 
                                           only.pos = TRUE,
                                           verbose = FALSE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.5)
  
  cat("Found", nrow(marker_defined_markers), "differentially expressed genes\n")
  
  # Get top 25 markers per cell type
  top25_marker_defined <- marker_defined_markers %>%
    group_by(cluster) %>%
    top_n(n = 25, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))
  
  # Scale genes if needed
  genes_to_scale <- unique(top25_marker_defined$gene)
  marker_subset <- ScaleData(marker_subset, features = genes_to_scale)
  
  # Create heatmap
  cat("\nGenerating marker-defined cell type heatmap...\n")
  
  n_genes_marker <- length(unique(top25_marker_defined$gene))
  # Use more space per gene to prevent compression
  heatmap_height_marker <- max(12, n_genes_marker / 25)  # More spacing
  
  marker_heatmap <- DoHeatmap(marker_subset, 
                              features = top25_marker_defined$gene,
                              group.by = "marker_defined_type",
                              size = 3,  # Slightly larger text
                              angle = 90,
                              draw.lines = TRUE,
                              label=FALSE) +
    scale_fill_gradientn(colors = c("#998ec3", "white", "#f1a340")) +
    #ggtitle("Top 25 Genes: Border Cells vs Polar Cells") +
    theme(text = element_text(size = 9),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.text.x.top = element_blank())
  ggsave("marker_defined_cell_heatmap.pdf", 
         plot = marker_heatmap,
         width = 6.5, 
         height = 6.5,
         units = "in",
         dpi = 300,
         device = "pdf")
  
  ggsave("marker_defined_cell_heatmap.png", 
         plot = marker_heatmap,
         width = 6.5, 
         height = 6,
         units = "in",
         dpi = 600)
  
  cat("Heatmap saved! (height:", heatmap_height_marker, "inches)\n")
  
  # Creat expression matrix for export
  cat("\nCreating expression matrices for export...\n")
  
  # Get normalized expression data for all genes in marker-defined cells
  expression_matrix <- GetAssayData(marker_subset, assay = "RNA", layer = "data")
  
  # Convert to regular matrix and add cell type annotations
  expression_df <- as.data.frame(as.matrix(expression_matrix))
  expression_df <- t(expression_df)  # Transpose so cells are rows
  
  # Add metadata
  cell_metadata <- marker_subset@meta.data %>%
    select(marker_defined_type, seurat_clusters, sample_id, 
           nFeature_RNA, nCount_RNA, percent.mt)
  
  # Combine
  expression_with_metadata <- cbind(cell_metadata, expression_df)
  
  # Save full expression matrix
  write.csv(expression_with_metadata, 
            "marker_defined_cells_expression_matrix_FULL.csv",
            row.names = TRUE)
  
  cat("Full expression matrix saved:", nrow(expression_with_metadata), "cells x", 
      ncol(expression_df), "genes\n")
  
  # Save a focused matrix with just top variable genes and key markers
  focused_genes <- unique(c(
    top25_marker_defined$gene,
    present_genes,
    border_genes_present,
    polar_genes_present
  ))
  focused_genes <- focused_genes[focused_genes %in% colnames(expression_df)]
  
  expression_focused <- cbind(cell_metadata, expression_df[, focused_genes])
  
  write.csv(expression_focused, 
            "marker_defined_cells_expression_matrix_FOCUSED.csv",
            row.names = TRUE)
  
  cat("Focused expression matrix saved:", nrow(expression_focused), "cells x", 
      length(focused_genes), "genes\n")
  
  # Save as RDS for easy loading in R
  saveRDS(list(
    expression_matrix = expression_matrix,
    metadata = cell_metadata,
    border_cells = border_cells_marker,
    polar_cells = polar_cells_marker
  ), file = "marker_defined_cells_data.rds")
  
  cat("RDS object saved for easy R access\n")
  
  # Feature plots of marker-defined cell populations
  cat("\nGenerating feature plots for marker-defined cells...\n")
  
  present_in_marker <- present_genes[present_genes %in% rownames(marker_subset)]
  
  if (length(present_in_marker) > 0) {
    n_genes_feat <- length(present_in_marker)
    n_cols_feat <- 3
    n_rows_feat <- ceiling(n_genes_feat / n_cols_feat)
    plot_width_feat <- 6.5
    plot_height_feat <- plot_width_feat * (n_rows_feat / n_cols_feat)
    
    marker_featureplot <- FeaturePlot(marker_subset, 
                                      features = present_in_marker,
                                      ncol = n_cols_feat,
                                      pt.size = 0.3,
                                      order = TRUE,
                                      cols=c("#A7A9AC","#f1a340")) &
      theme_classic() &
      theme(axis.text = element_text(size = 6),
            axis.title = element_text(size = 7),
            plot.title = element_text(size = 8, face = "italic"),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 7),
            aspect.ratio = 1)
    
    ggsave("marker_defined_featureplots.pdf", 
           plot = marker_featureplot,
           width = plot_width_feat, 
           height = plot_height_feat,
           units = "in",
           dpi = 300,
           device = "pdf")
    
    ggsave("marker_defined_featureplots.png", 
           plot = marker_featureplot,
           width = plot_width_feat, 
           height = plot_height_feat,
           units = "in",
           dpi = 600)
  }
  
  # UMAP showing marker-defined cells
  marker_umap <- DimPlot(marker_subset, 
                         group.by = "marker_defined_type",
                         cols = c("Border Cell" = "#998ec3", 
                                  "Polar Cell" = "#f1a340"),
                         pt.size = 0.01) +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          plot.title = element_text(size = 10),
          aspect.ratio = 1) +
    theme(legend.position = c(0.85, 0.15)) +
    ggtitle("Marker-Defined Border and Polar Cells")
  
  ggsave("marker_defined_umap.pdf", 
         plot = marker_umap,
         width = 3.25, 
         height = 3.25,
         units = "in",
         dpi = 300,
         device = "pdf")
  
  ggsave("marker_defined_umap.png", 
         plot = marker_umap,
         width = 3.25, 
         height = 3.25,
         units = "in",
         dpi = 600)
  
  # Export marker gene lists
  write.csv(marker_defined_markers, 
            "marker_defined_cell_DEGs_all.csv", 
            row.names = FALSE)
  
  write.csv(top25_marker_defined, 
            "marker_defined_cell_DEGs_top25.csv", 
            row.names = FALSE)
  
  # Summary statistics
  marker_summary <- marker_subset@meta.data %>%
    group_by(marker_defined_type) %>%
    summarise(
      n_cells = n(),
      median_genes = median(nFeature_RNA),
      median_UMIs = median(nCount_RNA),
      median_mt_pct = median(percent.mt),
      clusters_represented = paste(sort(unique(seurat_clusters)), collapse = ",")
    )
  
  cat("\nMarker-defined cell summary:\n")
  print(marker_summary)
  
  write.csv(marker_summary, 
            "marker_defined_cell_summary.csv", 
            row.names = FALSE)
  
  cat("\nMarker-defined cell analysis complete!\n")
  
} else {
  cat("WARNING: Too few marker-defined cells (", length(marker_defined_cells), ") to analyze\n")
}

# 7. Transgene analysis

cat("\n", rep("=", 70), "\n", sep = "")
cat("TRANSGENE-POSITIVE CELL ANALYSIS (ALL TRANSGENE+ CELLS)\n")
cat(rep("=", 70), "\n\n", sep = "")

# Check if transgenes are present
has_LifeAct <- "LifeActGFP" %in% rownames(merged_seurat)
has_dsRed <- "dsRedExpress" %in% rownames(merged_seurat)

if (has_LifeAct || has_dsRed) {
  
  # Identify positive cells (expression > 0)
  if (has_LifeAct) {
    LifeAct_cells <- WhichCells(merged_seurat, expression = LifeActGFP > 0)
    cat("LifeActGFP+ cells:", length(LifeAct_cells), "\n")
  } else {
    LifeAct_cells <- c()
  }
  
  if (has_dsRed) {
    dsRed_cells <- WhichCells(merged_seurat, expression = dsRedExpress > 0)
    cat("dsRedExpress+ cells:", length(dsRed_cells), "\n")
  } else {
    dsRed_cells <- c()
  }
  
  # Combine transgene+ cells
  transgene_cells <- unique(c(LifeAct_cells, dsRed_cells))
  cat("Total transgene+ cells:", length(transgene_cells), "\n\n")
  
  if (length(transgene_cells) > 50) {  # Only proceed if enough cells
    
    # Subset to transgene+ cells
    transgene_subset <- subset(merged_seurat, cells = transgene_cells)
    
    cat("Subset created with", ncol(transgene_subset), "cells\n")
    cat("Clusters represented:", paste(sort(unique(transgene_subset$seurat_clusters)), collapse = ", "), "\n\n")
    
    # Find markers for transgene+ cells by cluster
    cat("Finding markers for transgene+ cells...\n")
    Idents(transgene_subset) <- "seurat_clusters"
    transgene_markers <- FindAllMarkers(transgene_subset, 
                                        only.pos = TRUE,
                                        verbose = FALSE,
                                        min.pct = 0.25)
    
    cat("Found", nrow(transgene_markers), "markers\n")
    
    # Get top 25 markers per cluster for transgene+ cells
    top25_transgene <- transgene_markers %>%
      group_by(cluster) %>%
      top_n(n = 25, wt = avg_log2FC) %>%
      arrange(cluster, desc(avg_log2FC))
    
    # Scale the top markers if needed
    if (!all(unique(top25_transgene$gene) %in% rownames(LayerData(transgene_subset, layer = "scale.data")))) {
      cat("Scaling transgene marker genes...\n")
      transgene_subset <- ScaleData(transgene_subset, features = unique(top25_transgene$gene))
    }
    
    # Heatmap for transgene+ cells
    cat("Generating transgene+ cell heatmap...\n")
    n_genes_transgene <- length(unique(top25_transgene$gene))
    heatmap_height_transgene <- max(12, n_genes_transgene / 25)  # More generous spacing
    
    transgene_heatmap <- DoHeatmap(transgene_subset, 
                                   features = top25_transgene$gene,
                                   group.by = "seurat_clusters",
                                   size = 3,
                                   angle = 90,
                                   draw.lines = TRUE) +
      scale_fill_gradientn(colors = c("#998ec3", "white", "#f1a340")) +
      ggtitle("Top 25 Markers: Transgene+ Cells Only") +
      theme(text = element_text(size = 9),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 6))
    
    ggsave("transgene_positive_heatmap.pdf", 
           plot = transgene_heatmap,
           width = 6.5, 
           height = heatmap_height_transgene,
           units = "in",
           dpi = 300,
           device = "pdf")
    
    ggsave("transgene_positive_heatmap.png", 
           plot = transgene_heatmap,
           width = 6.5, 
           height = heatmap_height_transgene,
           units = "in",
           dpi = 600)
    
    cat("Transgene+ heatmap saved! (height:", heatmap_height_transgene, "inches)\n")
    
    # Feature plots for transgene+ subset
    cat("Generating transgene+ feature plots...\n")
    
    # Use the same marker genes
    present_in_subset <- present_genes[present_genes %in% rownames(transgene_subset)]
    
    if (length(present_in_subset) > 0) {
      n_genes_sub <- length(present_in_subset)
      n_cols_sub <- 3
      n_rows_sub <- ceiling(n_genes_sub / n_cols_sub)
      plot_width_sub <- 6.5
      plot_height_sub <- plot_width_sub * (n_rows_sub / n_cols_sub)
      
      transgene_featureplot <- FeaturePlot(transgene_subset, 
                                           features = present_in_subset,
                                           ncol = n_cols_sub,
                                           pt.size = 0.2,
                                           order = TRUE,
                                           cols=c("#A7A9AC","#f1a340")) &
        theme_classic() &
        theme(axis.text = element_text(size = 6),
              axis.title = element_text(size = 7),
              plot.title = element_text(size = 8, face = "italic"),
              legend.text = element_text(size = 6),
              legend.title = element_text(size = 7),
              aspect.ratio = 1)
      
      ggsave("transgene_positive_featureplots.pdf", 
             plot = transgene_featureplot,
             width = plot_width_sub, 
             height = plot_height_sub,
             units = "in",
             dpi = 300,
             device = "pdf")
      
      ggsave("transgene_positive_featureplots.png", 
             plot = transgene_featureplot,
             width = plot_width_sub, 
             height = plot_height_sub,
             units = "in",
             dpi = 600)
      
      cat("Transgene+ feature plots saved!\n")
    }
    
    # UMAP of transgene+ cells
    cat("Generating transgene+ UMAP...\n")
    
    transgene_umap <- DimPlot(transgene_subset, 
                              group.by = "seurat_clusters",
                              label = TRUE,
                              label.size = 3,
                              pt.size = 0.3) +
      theme_classic() +
      theme(axis.text = element_text(size = 7),
            axis.title = element_text(size = 8),
            legend.text = element_text(size = 7),
            plot.title = element_text(size = 10),
            aspect.ratio = 1) +
      ggtitle("Transgene+ Cells by Cluster")
    
    ggsave("transgene_positive_umap.pdf", 
           plot = transgene_umap,
           width = 5, 
           height = 5,
           units = "in",
           dpi = 300,
           device = "pdf")
    
    ggsave("transgene_positive_umap.png", 
           plot = transgene_umap,
           width = 5, 
           height = 5,
           units = "in",
           dpi = 600)
    
    # Export transgene+ markers
    write.csv(transgene_markers, 
              "transgene_positive_markers_all.csv", 
              row.names = FALSE)
    
    write.csv(top25_transgene, 
              "transgene_positive_markers_top25.csv", 
              row.names = FALSE)
    
    # Summary statistics
    transgene_summary <- transgene_subset@meta.data %>%
      group_by(seurat_clusters) %>%
      summarise(
        n_cells = n(),
        median_genes = median(nFeature_RNA),
        median_UMIs = median(nCount_RNA)
      )
    
    cat("\nTransgene+ cells by cluster:\n")
    print(transgene_summary)
    
    write.csv(transgene_summary, 
              "transgene_positive_summary.csv", 
              row.names = FALSE)
    
    cat("\nTransgene+ analysis complete!\n")
    
  } else {
    cat("WARNING: Too few transgene+ cells (", length(transgene_cells), ") to analyze\n")
  }
  
} else {
  cat("WARNING: Neither LifeActGFP nor dsRedExpress found in dataset\n")
}

# 8. Cluster UMAP
cat("\nGenerating cluster UMAP...\n")

cluster_umap <- DimPlot(merged_seurat, 
                        group.by = "seurat_clusters",
                        label = TRUE,
                        label.size = 2,
                        pt.size = 0.01) +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10),
        aspect.ratio = 1) +
  theme(legend.position = "none") +
  ggtitle("Cell Clusters")

ggsave("cluster_umap.pdf", 
       plot = cluster_umap,
       width = 3.25, 
       height = 3.25,
       units = "in",
       dpi = 300,
       device = "pdf")

ggsave("cluster_umap.png", 
       plot = cluster_umap,
       width = 3.25, 
       height = 3.25,
       units = "in",
       dpi = 600)

cat("Cluster UMAP saved!\n")

# 9. Sample distribution UMAP to ensure integration
cat("\nGenerating sample distribution UMAP...\n")

sample_umap <- DimPlot(merged_seurat, 
                       reduction = "umap", 
                       group.by = "sample_id",
                       pt.size = 0.1) +
  theme_classic() +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        plot.title = element_text(size = 10),
        aspect.ratio = 1) +
  ggtitle("Cells by Sample")

ggsave("sample_distribution_umap.pdf", 
       plot = sample_umap,
       width = 5.5, 
       height = 5.5,
       units = "in",
       dpi = 300,
       device = "pdf")

ggsave("sample_distribution_umap.png", 
       plot = sample_umap,
       width = 5.5, 
       height = 5.5,
       units = "in",
       dpi = 600)

cat("Sample distribution UMAP saved!\n")

# 10. Marker-defined cells on full UMAP
cat("\nGenerating marker-defined cells overlay on full UMAP...\n")

# Show marker-defined cells highlighted on the full dataset
marker_overlay <- DimPlot(merged_seurat,
                          group.by = "marker_defined_type",
                          cols = c("Border_Cell" = "#542788",
                                   "Polar_Cell" = "#f1a340",
                                   "Other" = "gray"),
                          pt.size = 0.01,
                          order = c("Polar_Cell", "Border_Cell", "Other")) +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 10),
        aspect.ratio = 1) +
  theme(legend.position = c(0.85, 0.20)) +
  ggtitle("Border and Polar Cells on Full UMAP")

ggsave("marker_defined_overlay_umap.pdf", 
       plot = marker_overlay,
       width = 3.25, 
       height = 3.25,
       units = "in",
       dpi = 300,
       device = "pdf")

ggsave("marker_defined_overlay_umap.png", 
       plot = marker_overlay,
       width = 3.25, 
       height = 3.25,
       units = "in",
       dpi = 600)

# 11. Annotate cell types by cluster
cat("\n" , rep("=", 70), "\n", sep = "")
cat("CLUSTER-BASED CELL TYPE ANNOTATION\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Please examine the dot plot (marker_gene_dotplot.png/pdf) and identify:\n")
cat("- Which clusters express border cell markers (e.g., slbo, fru)?\n")
cat("- Which clusters express polar cell markers (e.g., upd1, dsRedExpress)?\n\n")

# Display cluster numbers for reference
cat("Available clusters:", paste(sort(unique(merged_seurat$seurat_clusters)), collapse = ", "), "\n\n")

# Show expression summary for key genes if present
if ("slbo" %in% present_genes) {
  slbo_by_cluster <- FetchData(merged_seurat, vars = c("seurat_clusters", "slbo")) %>%
    group_by(seurat_clusters) %>%
    summarise(mean_slbo = mean(slbo), pct_expressing = sum(slbo > 0) / n() * 100) %>%
    arrange(desc(mean_slbo))
  cat("\nslbo expression by cluster (top 5):\n")
  print(head(slbo_by_cluster, 5))
}

if ("upd1" %in% present_genes) {
  upd1_by_cluster <- FetchData(merged_seurat, vars = c("seurat_clusters", "upd1")) %>%
    group_by(seurat_clusters) %>%
    summarise(mean_upd1 = mean(upd1), pct_expressing = sum(upd1 > 0) / n() * 100) %>%
    arrange(desc(mean_upd1))
  cat("\nupd1 expression by cluster (top 5):\n")
  print(head(upd1_by_cluster, 5))
}

if ("dsRedExpress" %in% present_genes) {
  dsRed_by_cluster <- FetchData(merged_seurat, vars = c("seurat_clusters", "dsRedExpress")) %>%
    group_by(seurat_clusters) %>%
    summarise(mean_dsRed = mean(dsRedExpress), pct_expressing = sum(dsRedExpress > 0) / n() * 100) %>%
    arrange(desc(mean_dsRed))
  cat("\ndsRedExpress expression by cluster (top 5):\n")
  print(head(dsRed_by_cluster, 5))
}

# Prompt user for input
cat("\nEnter border cell cluster numbers (comma-separated, e.g., 0,3,5):\n")
border_input <- readline(prompt = "> ")

cat("\nEnter polar cell cluster numbers (comma-separated, e.g., 1,2):\n")
polar_input <- readline(prompt = "> ")

# Parse input
border_cell_clusters <- as.numeric(unlist(strsplit(border_input, ",")))
polar_cell_clusters <- as.numeric(unlist(strsplit(polar_input, ",")))

cat("\nBorder cell clusters:", paste(border_cell_clusters, collapse = ", "), "\n")
cat("Polar cell clusters:", paste(polar_cell_clusters, collapse = ", "), "\n\n")

# Create cell type annotation
merged_seurat$cell_type <- "Other"

# Assign border cells
merged_seurat$cell_type[merged_seurat$seurat_clusters %in% border_cell_clusters] <- "Border_Cell"

# Assign polar cells
merged_seurat$cell_type[merged_seurat$seurat_clusters %in% polar_cell_clusters] <- "Polar_Cell"

# Display summary
cat("\nCell type distribution:\n")
print(table(merged_seurat$cell_type))

cat("\nCell type by cluster:\n")
print(table(merged_seurat$cell_type, merged_seurat$seurat_clusters))

# Compare marker-defined vs cluster-defined
cat("\n=== COMPARISON: Marker-Defined vs Cluster-Defined ===\n")
comparison_table <- table(
  Marker_Defined = merged_seurat$marker_defined_type,
  Cluster_Defined = merged_seurat$cell_type
)
print(comparison_table)

write.csv(as.data.frame.matrix(comparison_table),
          "marker_vs_cluster_comparison.csv")

# 12. Cell type UMAP
cat("\nGenerating cell type UMAP...\n")

celltype_umap <- DimPlot(merged_seurat, 
                         group.by = "cell_type",
                         cols = c("Border_Cell" = "#998ec3",
                                  "Polar_Cell" = "#f1a340",
                                  "Other" = "lightgray"),
                         label = TRUE,
                         label.size = 3,
                         pt.size = 0.1,
                         repel = TRUE) +
  theme_classic() +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        plot.title = element_text(size = 10),
        aspect.ratio = 1) +
  ggtitle("Cell Type Annotations (Cluster-Based)")

ggsave("celltype_umap.pdf", 
       plot = celltype_umap,
       width = 5.5, 
       height = 5.5,
       units = "in",
       dpi = 300,
       device = "pdf")

ggsave("celltype_umap.png", 
       plot = celltype_umap,
       width = 5.5, 
       height = 5.5,
       units = "in",
       dpi = 600)

cat("Cell type UMAP saved!\n")

# 15. Export data tables
cat("\nExporting data tables...\n")

# Export all cluster markers
write.csv(cluster_markers, 
          "cluster_markers_all.csv", 
          row.names = FALSE)

# Export top 25 markers per cluster
write.csv(top25_markers, 
          "cluster_markers_top25.csv", 
          row.names = FALSE)

# Export cell type counts by sample
cell_counts <- as.data.frame(table(merged_seurat$cell_type, merged_seurat$sample_id))
colnames(cell_counts) <- c("Cell_Type", "Sample", "Count")
write.csv(cell_counts, 
          "cell_type_counts_by_sample.csv", 
          row.names = FALSE)

# Export cell type proportions
cell_proportions <- cell_counts %>%
  group_by(Sample) %>%
  mutate(Proportion = Count / sum(Count) * 100) %>%
  arrange(Sample, Cell_Type)
write.csv(cell_proportions, 
          "cell_type_proportions_by_sample.csv", 
          row.names = FALSE)

cat("Data tables exported!\n")

# 16. Save annotated Seurat object

cat("\nSaving annotated Seurat object...\n")

output_seurat_path <- file.path(dirname(seurat_object_path), "merged_seurat_annotated.rds")
saveRDS(merged_seurat, file = output_seurat_path)

cat("Annotated object saved to:", output_seurat_path, "\n")


# SUMMARY
cat("Output directory:", output_dir, "\n\n")

# 17. RhoGTPase and GAP feature plots (overlaid)
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("GENERATING RHOGTPASE FEATURE PLOTS (OVERLAID)\n")
cat(rep("=", 70), "\n\n", sep = "")

# Define gene groups (each list is one plot)
gene_groups <- list(
  Group1 = c("RhoGAP15B", "RhoGAP16F", "Cdc42"),
  Group2 = c("Ocrl", "Rab5", "Rab11"),
  Group3 = c("RacGAP84C", "RhoGAP19D", "RhoGAP5A", "RhoGAP92B", "Rlip", "tum", "Rac1", "Rac2"),
  Group4 = c("CG42684", "Nf1", "RasGAP1", "Ras64B", "Ras85D"),
  Group5 = c("cv-c", "Rho1"),
  Group6 = c("Graf", "RhoGAP102A", "RhoGAPp190", "Rho1", "Rac1", "Rac2"),
  Group7 = c("RhoGAP18B", "RhoGAP1A", "RhoGAP68F", "RhoGAP71E", "Rho1", "Rac1", "Rac2", "Cdc42"),
  Group8 = c("RhoGAP93B", "Rac1", "Rac2", "Cdc42"),
  Group9 = c("Zir", "Cdc42", "Rac1")
)

# Create RhoGTPase subdirectory
rho_output_dir <- file.path(output_dir, "RhoGTPase_FeaturePlots")
if (!dir.exists(rho_output_dir)) {
  dir.create(rho_output_dir, recursive = TRUE)
}

# Generate plots for each group
all_rho_plots <- list()

for (i in seq_along(gene_groups)) {
  group_name <- names(gene_groups)[i]
  genes <- gene_groups[[i]]
  
  cat("\n--- ", group_name, " ---\n", sep = "")
  cat("Genes:", paste(genes, collapse = ", "), "\n")
  
  # Check which genes are present
  genes_present <- genes[genes %in% rownames(merged_seurat)]
  genes_missing <- genes[!genes %in% rownames(merged_seurat)]
  
  if (length(genes_missing) > 0) {
    cat("WARNING - ", group_name, " missing genes:", 
        paste(genes_missing, collapse = ", "), "\n")
  }
  
  if (length(genes_present) == 0) {
    cat("ERROR - ", group_name, ": No genes found in dataset\n")
    next
  }
  
  # For single gene, create standard feature plot
  if (length(genes_present) == 1) {
    cat("Only 1 gene present, creating standard feature plot\n")
    p <- FeaturePlot(merged_seurat,
                     features = genes_present,
                     pt.size = 0.3,
                     order = TRUE) +
      theme_classic() +
      theme(axis.text = element_text(size = 8),
            axis.title = element_text(size = 9),
            plot.title = element_text(size = 10, face = "italic"),
            legend.text = element_text(size = 8),
            aspect.ratio = 1)
    
    plot_width <- 5
    plot_height <- 5
    
  } else {
    # For 2+ genes, create manual overlay plot with custom colors
    cat("Creating overlay plot for", length(genes_present), "genes\n")
    
    # Get UMAP coordinates
    umap_coords <- Embeddings(merged_seurat, reduction = "umap")
    umap_df <- as.data.frame(umap_coords)
    colnames(umap_df) <- c("UMAP_1", "UMAP_2")
    
    # Get expression data for all genes
    expr_data <- FetchData(merged_seurat, vars = genes_present)
    
    # Combine UMAP and expression
    plot_data <- cbind(umap_df, expr_data)
    
    # Create a categorical variable for which gene(s) each cell expresses
    plot_data$expressing <- "None"
    plot_data$n_expressing <- 0
    
    for (gene in genes_present) {
      expressing_cells <- plot_data[[gene]] > 0
      plot_data$n_expressing[expressing_cells] <- plot_data$n_expressing[expressing_cells] + 1
    }
    
    # Assign expression categories
    for (j in 1:nrow(plot_data)) {
      if (plot_data$n_expressing[j] == 0) {
        plot_data$expressing[j] <- "None"
      } else if (plot_data$n_expressing[j] == 1) {
        # Find which single gene is expressed
        for (gene in genes_present) {
          if (plot_data[[gene]][j] > 0) {
            plot_data$expressing[j] <- gene
            break
          }
        }
      } else {
        # Multiple genes expressed
        plot_data$expressing[j] <- "Combination"
      }
    }
    
    # Assign colors from custom palette
    n_genes <- length(genes_present)
    # Use first n colors from custom palette for genes
    gene_colors <- custom_colors[1:n_genes]
    names(gene_colors) <- genes_present
    
    # Add gray for non-expressing and black for combinations
    color_map <- c("None" = "#A7A9AC",      # gray
                   gene_colors, 
                   "Combination" = "#000000")  # black
    
    # Set factor levels to control legend order
    plot_data$expressing <- factor(plot_data$expressing, 
                                   levels = c("None", genes_present, "Combination"))
    
    # Reorder so expressing cells are plotted on top
    plot_data <- plot_data[order(plot_data$expressing == "None"), ]
    
    # Create plot
    p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = expressing)) +
      geom_point(size = 0.3, alpha = 0.8) +
      scale_color_manual(values = color_map, 
                         name = "Expression",
                         drop = FALSE) +
      theme_classic() +
      theme(axis.text = element_text(size = 8),
            axis.title = element_text(size = 9),
            legend.text = element_text(size = 8, face = "italic"),
            legend.title = element_text(size = 9),
            aspect.ratio = 1) +
      labs(title = paste(group_name, "Expression")) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    
    plot_width <- 6
    plot_height <- 5
  }
  
  # Save as pdf
  pdf_filename <- file.path(rho_output_dir, paste0(group_name, "_overlay_featureplot.pdf"))
  ggsave(pdf_filename,
         plot = p,
         width = plot_width,
         height = plot_height,
         units = "in",
         dpi = 300,
         device = "pdf")
  
  # Save as PNG
  png_filename <- file.path(rho_output_dir, paste0(group_name, "_overlay_featureplot.png"))
  ggsave(png_filename,
         plot = p,
         width = plot_width,
         height = plot_height,
         units = "in",
         dpi = 600)
  
  cat("Saved:", group_name, "(", plot_width, "x", plot_height, "inches)\n")
  
  all_rho_plots[[group_name]] <- list(plot = p, genes = genes_present)
}

cat("RHOGTPASE OVERLAY PLOTS COMPLETE\n")
cat("All plots saved in:", rho_output_dir, "\n")
cat("      Co-expressing cells shown as 'Combination' in black\n\n")

# Final summary
cat("Output directories:\n")
cat("  Main figures:", output_dir, "\n")
cat("  RhoGTPase analyses:", rho_output_dir, "\n\n")