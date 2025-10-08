# Mostly cribbed from https://satijalab.org/seurat/articles/integration_introduction

library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(cowplot)

out_dir <- "results/singlecell_seurat_integrated_ifnb/"
dir.create(out_dir, recursive = T, showWarnings = F)

# install dataset
InstallData("ifnb")

# load dataset
ifnb <- LoadData("ifnb")
# split the RNA measurements into two layers one for control cells, one for stimulated cells

ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb

# first, to see the difference, run the analysis without integration
# run standard analysis workflow
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb) # This reduction is key for integration

ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))
ggsave(paste0(out_dir, "unintegrated_umap.png"), width=12, height=5)

# Now we run integration
# Goal is to identify common cell types between the samples
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, 
                        orig.reduction = "pca", 
                        new.reduction = "integrated.cca",
                        verbose = FALSE)

# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

# same Umap steps as before, this time integrated
ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")

# Visualization
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))
ggsave(paste0(out_dir, "integrated_umap.png"), width=12, height=5)

DimPlot(ifnb, reduction = "umap", split.by = "stim")
ggsave(paste0(out_dir, "integrated_umap.stim_split.png"), width=10, height=5)

# find conserved markers, irrespecitive of condition
# (useful for cell type determination)
Idents(ifnb) <- "seurat_annotations" # cell types are already included!
nk.markers <- FindConservedMarkers(ifnb, ident.1 = "NK", grouping.var = "stim", verbose = FALSE)
head(nk.markers)

# example dot plot of conserved markers
# NEEDS TO BE FIXED AND SET ORDER CORRECTLY
Idents(ifnb) <- factor(Idents(ifnb), levels = c("pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono",
                                                "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T"))

markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
DotPlot(ifnb, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
  RotatedAxis()
ggsave(paste0(out_dir, "cell_marker_dot_plot.png"), width=11, height=8, bg="white")

# differential expressed genes across conditions
theme_set(theme_cowplot())

# aggregate cells to create "pseudobulk" 
aggregate_ifnb <- AggregateExpression(ifnb, group.by = c("seurat_annotations", "stim"), 
                                      return.seurat = TRUE)
Idents(aggregate_ifnb)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")

p1 <- CellScatter(aggregate_ifnb, "CD14 Mono_CTRL", "CD14 Mono_STIM", highlight = genes.to.label)
p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

p3 <- CellScatter(aggregate_ifnb, "CD4 Naive T_CTRL", "CD4 Naive T_STIM", highlight = genes.to.label)
p4 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)

p2 + p4
ggsave(paste0(out_dir, "example_pseudobulk_expression_plots.png"),
       width=10, height=5, bg="white")

# since we only have one replicate, we can do typical find markers

ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"
b.interferon.response <- FindMarkers(ifnb, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 15)

# can do feature plot to compare cell markers (CD3D and GNLY) to affected gene IFI6
FeaturePlot(ifnb, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3, 
            cols = c("grey","red"), reduction = "umap")
ggsave(paste0(out_dir, "example_features_plot.stim_split.png"), width=8, height=10,
       bg="white")

# example 
plots <- VlnPlot(ifnb, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "seurat_annotations",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
ggsave(paste0(out_dir, "example_violin_gene_expression_plots.png"), width=8, height=10,
       bg="white")

# save the integrated seurat object
SeuratDisk::SaveH5Seurat(ifnb, paste0(out_dir, "ifnb_processed.h5Seurat"))

SaveSeuratRds(ifnb, paste0(out_dir, "ifnb_processed.seurat_rds.RDS"))
saveRDS(ifnb, paste0(out_dir, "ifnb_processed.RDS"))

# how would DE look with average expression?
# adapted code from here 
# https://satijalab.org/seurat/articles/parsebio_sketch_integration#compare-healthy-and-diabetic-samples

# recommend splitting up seurat object based on cell type

# Note: Due to lack of replicates, this does not work

# aggregate_ifnb_cd14 <- subset(aggregate_ifnb, subset = seurat_annotations=="CD14 Mono")
# Idents(aggregate_ifnb_cd14) <- "stim"
# 
# de_markers <- FindMarkers(aggregate_ifnb_cd14, 
#                           ident.1 = "STIM", ident.2 = "CTRL", 
#                           slot = "counts", test.use = "DESeq2",
#                           verbose = F)
# de_markers$gene <- rownames(de_markers)
# ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
#   ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,
#                                                                           "")), colour = "red", size = 3)
