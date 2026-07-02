library(openxlsx)
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(plyr)
library(dplyr)

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/GSE133283/"


out_dir <- paste0(root_dir, "results/endothelial_gene_exp_analysis/")
dir.create(out_dir, showWarnings = F, recursive = T)

# read in Seurat object

seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/integrate_seurat_samples/all_samples.integrated_seurat.RDS"))

# ensure correct assay
Idents(seu_obj) <- "harmony_clusters"

# set assay to RNA and join layers
DefaultAssay(seu_obj) <- "RNA"
seu_obj <- JoinLayers(seu_obj)

seu_obj <- NormalizeData(seu_obj)
seu_obj <- FindVariableFeatures(seu_obj)
seu_obj <- ScaleData(seu_obj)



# read in cassia results
cassia_results <- readRDS(paste0(root_dir, "results/cassia_results_analysis/cassia_results_comparison.RDS"))




# make some UMAPS
DimPlot(seu_obj, reduction="umap.harmony", group.by= "harmony_clusters", label=T)
ggsave(paste0(out_dir, "total_cell_umap.png"), width=7, height=5)

# subset the endothelial cells
endo_cells <- cassia_results[grepl("[Ee]ndothelial", cassia_results$gemini_prediction),]

endo_clusters <- gsub("cluster", "", endo_cells$harmony_clusters)

# subset the seurat
endo_seu <- subset(seu_obj, subset = harmony_clusters %in% endo_clusters)


# clean it up
umap_coords <- Embeddings(endo_seu, reduction = "umap.harmony")
umap_coords <- as.data.frame(umap_coords)

outliers <- umap_coords[umap_coords$umapharmony_1 < 0 |
                          umap_coords$umapharmony_2 < -5,]

endo_seu$cell_id <- colnames(endo_seu)

endo_seu <- subset(endo_seu, subset = !cell_id %in% rownames(outliers))

# plot the UMAP
DimPlot(endo_seu, reduction="umap.harmony", group.by= "harmony_clusters", label=T, raster = F, label.box=T)
ggsave(paste0(out_dir, "endothelial_cell_umap.png"), width=6, height=5)

# exam the genes of interest
group1_genes <- c("Il1r1",
                "Tnfrsf1a", 
                "Ifngr1",
                "Ifngr2")

group1_genes %in% rownames(endo_seu)

VlnPlot(endo_seu, features = group1_genes)
ggsave(paste0(out_dir, "group1_violin_plots.png"), width=9, height=5)

DotPlot(endo_seu, features=group1_genes) + labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "group1_dot_plots.png"), width=6, height=5, bg="white")


FeaturePlot(endo_seu, features = group1_genes, reduction="umap.harmony")
ggsave(paste0(out_dir, "group1_umap.png"), width=6, height=5)

# look at the next set of genes of interest
group2_genes <- c("B2m",
                  "Tap1",
                  "Cxcl10",
                  "H2-Ab1")

group2_genes %in% rownames(endo_seu)

VlnPlot(endo_seu, features = group2_genes, ncol = 2)
ggsave(paste0(out_dir, "group2_violin_plots.png"), width=7, height=7)

DotPlot(endo_seu, features=group2_genes) + labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "group2_dot_plots.png"), width=6, height=5, bg="white")


FeaturePlot(endo_seu, features = group2_genes, reduction="umap.harmony")
ggsave(paste0(out_dir, "group2_umap.png"), width=6, height=5)


# determine the unique markers for each endothelial cluster
endo_markers <- FindAllMarkers(endo_seu)
endo_markers$pct_diff <- endo_markers$pct.1 - endo_markers$pct.2
endo_markers <- endo_markers[order(endo_markers$pct_diff, decreasing=T),]

total_markers <- readRDS(paste0(root_dir, "results/find_all_markers/all_markers.RDS"))
total_markers$pct_diff <- total_markers$pct.1 - total_markers$pct.2
total_markers <- total_markers[order(total_markers$pct_diff, decreasing=T),]

# pull out top endomarkers
top_endo_markers <- lapply(endo_clusters, function(x) {
  
  head(endo_markers[endo_markers$cluster == x &
                      endo_markers$pct_diff > 0,])
  
})
top_endo_markers <- bind_rows(top_endo_markers)

DotPlot(endo_seu, features=unique(top_endo_markers$gene)) + RotatedAxis() +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "top_endothelial_cluster_markers.inter_endothelial.png"), width=12, height=5,
       bg="white")

# pull out top endomarkers
top_total_markers <- lapply(endo_clusters, function(x) {
  
  head(total_markers[total_markers$cluster == x &
                      total_markers$pct_diff > 0,])
  
})
top_total_markers <- bind_rows(top_total_markers)


DotPlot(endo_seu, features=unique(top_total_markers$gene)) + RotatedAxis() +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "top_endothelial_cluster_markers.all_clusters.png"), width=12, height=5,
       bg="white")

# split things by condition

endo_seu$condition <- gsub("\\-[0-9]", "", endo_seu$orig.ident)

# subset even more!

endo_seu_ae <- subset(endo_seu, subset = condition %in% c("Adult", "EAE"))

# further splits
endo_seu_adult <- subset(endo_seu, subset = condition == "Adult")
endo_seu_eae <- subset(endo_seu, subset = condition == "EAE")


# group 1 plots

FeaturePlot(endo_seu_ae, features = group1_genes, reduction="umap.harmony", split.by = "condition")
ggsave(paste0(out_dir, "group1_umap.condition_split.png"), width=6, height=8)


VlnPlot(endo_seu_ae, features = group1_genes, split.by = "condition", ncol=1)
ggsave(paste0(out_dir, "group1_violin_plots.condition_split.png"), width=6, height=10)

DotPlot(endo_seu_ae, features=group1_genes, split.by = "condition") + labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "group1_dot_plots.condition_split.png"), width=6, height=6, bg="white")

DotPlot(endo_seu_adult, features=group1_genes) + labs(x=NULL, y=NULL, title="Adult")
DotPlot(endo_seu_eae, features = group1_genes) + labs(x=NULL, y=NULL, title="EAE")

# group 2 plots

FeaturePlot(endo_seu_ae, features = group2_genes, reduction="umap.harmony", split.by = "condition")
ggsave(paste0(out_dir, "group2_umap.condition_split.png"), width=6, height=8)


VlnPlot(endo_seu_ae, features = group2_genes, split.by = "condition", ncol=1)
ggsave(paste0(out_dir, "group2_violin_plots.condition_split.png"), width=6, height=10)

DotPlot(endo_seu_ae, features=group2_genes, split.by = "condition") + labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "group2_dot_plots.condition_split.png"), width=6, height=6, bg="white")

# group 3!
group3_genes <- c("Stat1","Ciita", "Icam1", "Vcam1")

group3_genes %in% rownames(endo_seu)

VlnPlot(endo_seu, features = group3_genes, ncol = 2)
ggsave(paste0(out_dir, "group3_violin_plots.png"), width=7, height=7)

DotPlot(endo_seu, features=group3_genes) + labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "group3_dot_plots.png"), width=6, height=5, bg="white")


FeaturePlot(endo_seu, features = group3_genes, reduction="umap.harmony")
ggsave(paste0(out_dir, "group3_umap.png"), width=6, height=5)


# split

FeaturePlot(endo_seu_ae, features = group3_genes, reduction="umap.harmony", split.by = "condition")
ggsave(paste0(out_dir, "group3_umap.condition_split.png"), width=6, height=8)


VlnPlot(endo_seu_ae, features = group3_genes, split.by = "condition", ncol=1)
ggsave(paste0(out_dir, "group3_violin_plots.condition_split.png"), width=6, height=10)

DotPlot(endo_seu_ae, features=group3_genes, split.by = "condition") + labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "group3_dot_plots.condition_split.png"), width=6, height=6, bg="white")

