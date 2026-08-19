library(Seurat)
library(SeuratObject)
library(ggplot2)
library(cowplot)
library(openxlsx)
library(stringr)
library(snakecase)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/immune_non_immune_cell_markers/")
dir.create(out_dir, showWarnings = F)

# read in integrated seurat
int_seu <- LoadSeuratRds(paste0(root_dir,
                                "results/celltype_naming/all_samples.celltype_named.seurat.RDS"))

# set assay to RNA and join layers
DefaultAssay(int_seu) <- "RNA"
int_seu <- JoinLayers(int_seu)

int_seu <- NormalizeData(int_seu)
int_seu <- FindVariableFeatures(int_seu)
int_seu <- ScaleData(int_seu)

Idents(int_seu) <- "harmony_clusters"

# read in markers
markers_raw <- read.xlsx(paste0(root_dir, "Immune and non Immune markers.xlsx"), sep.names = " ")

# convert to marker list

markers_list <- lapply(colnames(markers_raw), function(x) {
  
  markers <- markers_raw[,x]
  
  markers <- markers[!is.na(markers)]
  
  markers <- str_to_title(markers)
  
  return(markers)
  
})
names(markers_list) <- colnames(markers_raw)

# make sure all markers are in the data
genes <- rownames(int_seu)

lapply(markers_list, function(markers) {
  
  markers[!markers %in% genes]
  
})

# make some plots
height <- 8

for (cell_type in names(markers_list)) {
  
  markers <- markers_list[[cell_type]]
  
  width <- length(markers)-5
  
  if (width < 5) {
    width <- 5
  }
  
  DotPlot(int_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cell_type) 
  ggsave(paste0(out_dir, to_snake_case(cell_type), ".markers_dot_plot.png"),
         width=width, height=height, bg="white")
  
  
}

# further request, UMAP plots

# just do a plot for each marker, too complicated to figure out how to auto size

umap_dir <- paste0(out_dir, "gene_umaps/")

dir.create(umap_dir, showWarnings = F)

cell_type <- "Microglia"

markers <- markers_list[[cell_type]]

plot_list <- lapply(markers, function(gene) {
  
  p1 <- FeaturePlot(int_seu, features = gene, reduction="umap.harmony_filtered")
  p2 <- LabelClusters(plot = p1, id = "ident", repel=F)
  
  print(p2)
  ggsave(paste0(umap_dir, to_snake_case(cell_type), ".", gene, ".expr_umap.png"),
         width=6, height=5)
  
  return(p2)
  
})

plot_grid(plotlist = plot_list, ncol=5)
ggsave(paste0(out_dir, to_snake_case(cell_type), ".markers_expr_umap.png"),
       width=30, height=20)

# BAMS -------------------------------------------------------------------------

cell_type <- "BAMs"

markers <- markers_list[[cell_type]]

plot_list <- lapply(markers, function(gene) {
  
  p1 <- FeaturePlot(int_seu, features = gene, reduction="umap.harmony_filtered")
  p2 <- LabelClusters(plot = p1, id = "ident", repel=F)
  
  print(p2)
  ggsave(paste0(umap_dir, to_snake_case(cell_type), ".", gene, ".expr_umap.png"),
         width=6, height=5)
  
  return(p2)
  
})

plot_grid(plotlist = plot_list, ncol=2)
ggsave(paste0(out_dir, to_snake_case(cell_type), ".markers_expr_umap.png"),
       width=12, height=10)

# Neutro -------------------------------------------------------------------------

cell_type <- "Neutro"

markers <- markers_list[[cell_type]]

plot_list <- lapply(markers, function(gene) {
  
  p1 <- FeaturePlot(int_seu, features = gene, reduction="umap.harmony_filtered")
  p2 <- LabelClusters(plot = p1, id = "ident", repel=F)
  
  print(p2)
  ggsave(paste0(umap_dir, to_snake_case(cell_type), ".", gene, ".expr_umap.png"),
         width=6, height=5)
  
  return(p2)
  
})

plot_grid(plotlist = plot_list, ncol=3)
ggsave(paste0(out_dir, to_snake_case(cell_type), ".markers_expr_umap.png"),
       width=18, height=10, bg="white")

# Lymphoids -------------------------------------------------------------------------

cell_type <- "Lymphoids"

markers <- markers_list[[cell_type]]

plot_list <- lapply(markers, function(gene) {
  
  p1 <- FeaturePlot(int_seu, features = gene, reduction="umap.harmony_filtered")
  p2 <- LabelClusters(plot = p1, id = "ident", repel=F)
  
  print(p2)
  ggsave(paste0(umap_dir, to_snake_case(cell_type), ".", gene, ".expr_umap.png"),
         width=6, height=5)
  
  return(p2)
  
})

plot_grid(plotlist = plot_list, ncol=4)
ggsave(paste0(out_dir, to_snake_case(cell_type), ".markers_expr_umap.png"),
       width=24, height=15)

# Myeloid -------------------------------------------------------------------------

cell_type <- "Myeloid"

markers <- markers_list[[cell_type]]

plot_list <- lapply(markers, function(gene) {
  
  p1 <- FeaturePlot(int_seu, features = gene, reduction="umap.harmony_filtered")
  p2 <- LabelClusters(plot = p1, id = "ident", repel=F)
  
  print(p2)
  ggsave(paste0(umap_dir, to_snake_case(cell_type), ".", gene, ".expr_umap.png"),
         width=6, height=5)
  
  return(p2)
  
})

plot_grid(plotlist = plot_list, ncol=4)
ggsave(paste0(out_dir, to_snake_case(cell_type), ".markers_expr_umap.png"),
       width=24, height=20, bg="white")

# Astrocytes -------------------------------------------------------------------------

cell_type <- "Astrocytes"

markers <- markers_list[[cell_type]]

plot_list <- lapply(markers, function(gene) {
  
  p1 <- FeaturePlot(int_seu, features = gene, reduction="umap.harmony_filtered")
  p2 <- LabelClusters(plot = p1, id = "ident", repel=F)
  
  print(p2)
  ggsave(paste0(umap_dir, to_snake_case(cell_type), ".", gene, ".expr_umap.png"),
         width=6, height=5)
  
  return(p2)
  
})

plot_grid(plotlist = plot_list, ncol=4)
ggsave(paste0(out_dir, to_snake_case(cell_type), ".markers_expr_umap.png"),
       width=24, height=15, bg="white")

# OPC and oligodendrocytes -------------------------------------------------------------------------

cell_type <- "OPC and oligodendrocytes"

markers <- markers_list[[cell_type]]

plot_list <- lapply(markers, function(gene) {
  
  p1 <- FeaturePlot(int_seu, features = gene, reduction="umap.harmony_filtered")
  p2 <- LabelClusters(plot = p1, id = "ident", repel=F)
  
  print(p2)
  ggsave(paste0(umap_dir, to_snake_case(cell_type), ".", gene, ".expr_umap.png"),
         width=6, height=5)
  
  return(p2)
  
})

plot_grid(plotlist = plot_list, ncol=4)
ggsave(paste0(out_dir, to_snake_case(cell_type), ".markers_expr_umap.png"),
       width=24, height=15, bg="white")


# Neurons -------------------------------------------------------------------------

cell_type <- "Neurons"

markers <- markers_list[[cell_type]]

plot_list <- lapply(markers, function(gene) {
  
  p1 <- FeaturePlot(int_seu, features = gene, reduction="umap.harmony_filtered")
  p2 <- LabelClusters(plot = p1, id = "ident", repel=F)
  
  print(p2)
  ggsave(paste0(umap_dir, to_snake_case(cell_type), ".", gene, ".expr_umap.png"),
         width=6, height=5)
  
  return(p2)
  
})

plot_grid(plotlist = plot_list, ncol=3)
ggsave(paste0(out_dir, to_snake_case(cell_type), ".markers_expr_umap.png"),
       width=18, height=15, bg="white")

# Endothelial, Pericytes and ependymal -------------------------------------------------------------------------

cell_type <- "Endothelial, Pericytes and ependymal"

markers <- markers_list[[cell_type]]

plot_list <- lapply(markers, function(gene) {
  
  p1 <- FeaturePlot(int_seu, features = gene, reduction="umap.harmony_filtered")
  p2 <- LabelClusters(plot = p1, id = "ident", repel=F)
  
  print(p2)
  ggsave(paste0(umap_dir, to_snake_case(cell_type), ".", gene, ".expr_umap.png"),
         width=6, height=5)
  
  return(p2)
  
})

plot_grid(plotlist = plot_list, ncol=4)
ggsave(paste0(out_dir, to_snake_case(cell_type), ".markers_expr_umap.png"),
       width=24, height=15, bg="white")


# DC -------------------------------------------------------------------------

cell_type <- "DC"

markers <- markers_list[[cell_type]]

plot_list <- lapply(markers, function(gene) {
  
  p1 <- FeaturePlot(int_seu, features = gene, reduction="umap.harmony_filtered")
  p2 <- LabelClusters(plot = p1, id = "ident", repel=F)
  
  print(p2)
  ggsave(paste0(umap_dir, to_snake_case(cell_type), ".", gene, ".expr_umap.png"),
         width=6, height=5)
  
  return(p2)
  
})

plot_grid(plotlist = plot_list, ncol=2)
ggsave(paste0(out_dir, to_snake_case(cell_type), ".markers_expr_umap.png"),
       width=12, height=5, bg="white")


FeaturePlot(int_seu, features = "Matn4", reduction="umap.harmony_filtered")
FeaturePlot(int_seu, features = "Enpp6", reduction="umap.harmony_filtered")

