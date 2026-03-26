library(Seurat)
library(SeuratObject)
library(ggplot2)
library(cowplot)
library(openxlsx)
library(tibble)
library(dplyr)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/find_marker_clusters/")
dir.create(out_dir, showWarnings = F)


# read in integrated seurat
int_seu <- LoadSeuratRds(paste0(root_dir,
                        "results/integrate_seurat_samples/all_samples.integrated_seurat.RDS"))

# prepare data for using find markers

# set assay to RNA and join layers
DefaultAssay(int_seu) <- "RNA"
int_seu <- JoinLayers(int_seu)

int_seu <- NormalizeData(int_seu)
int_seu <- FindVariableFeatures(int_seu)
int_seu <- ScaleData(int_seu)

Idents(int_seu) <- "harmony_clusters"

# use the usual markers
lab_cluster_markers <- list("B Cells"=c("Cd19", "Cd79a", "Ms4a1"),
                            "Monocytes"=c("Ccr2", "Cd44"),
                            "Neutrophils"=c("Ly6g", "Ngp", "Mmp8"),
                            "Macrophages"=c("Pf4", "Mrc1", "Ms4a7"),
                            "T Cells"=c("Trbc2", "Cd3d", "Lck"),
                            "Microglia"=c("Sall1", "Hexb", "P2ry12"),
                            "Proliferating Cells"=c("Mki67", "Ccnb1", "Tpx2"),
                            "Pericytes"=c("Pdgfrb", "Rgs5"))

# cow dot plot for the markers
dotplot_list <- lapply(names(lab_cluster_markers), function(cluster_name) {
  
  markers <- lab_cluster_markers[[cluster_name]]
  
  DotPlot(int_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) 
  
})

plot_grid(plotlist = dotplot_list, nrow = 2)

ggsave(paste0(out_dir, "lab_cluster_marker_dotplots.png"), width=16, height=12, bg="white")

# subset for agnostic cluster markers

set.seed(42)

sampled_cells <- sample(colnames(int_seu), size=10000, replace = F)

subset_seu <- subset(int_seu, cells=sampled_cells)

clusters <- sort(unique(subset_seu$harmony_clusters))

cluster_markers <- lapply(clusters, function(cluster) {
  
  markers <- FindMarkers(subset_seu,
                         ident.1 = cluster,
                         test.use = "wilcox")
  markers <- rownames_to_column(markers, var = "gene")
  
  markers$pct_delta <- markers$pct.1 - markers$pct.2
  
  markers$cluster <- cluster
    
  return(markers[order(markers$pct_delta, decreasing = T),])
  
})

# save to file
write.xlsx(cluster_markers,
          paste0(out_dir, "cluster_markers.xlsx"),
          colWidths="auto")

# merge up for looking
cluster_markers_df <- bind_rows(cluster_markers)

# output top 10 for each cluster
cluster_dir <- paste0(out_dir, "cluster_marker_dot_plots/")

dir.create(cluster_dir, showWarnings = F)

for (data in cluster_markers) {
  
  top_markers <- data$gene[1:10]
  
  cluster_name <- unique(data$cluster)
  
  DotPlot(int_seu, features = top_markers) + 
    RotatedAxis() + 
    labs(x=NULL, y=NULL, title=cluster_name) 
  ggsave(paste0(cluster_dir, "cluster", cluster_name, ".top_markers_dots.png"),
         width=6, height=10, bg="white")
  
}


# distribuition of clusters
metadata <- int_seu@meta.data

ggplot(metadata,
       aes(x=orig.ident,
           fill=orig.ident)) +
  geom_bar(color="black") +
  theme_bw() +
  guides(fill="none") +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "cell_counts_per_sample.png"), width=5, height=4)

ggplot(metadata,
       aes(x=orig.ident,
           fill=orig.ident)) +
  geom_bar(color="black") +
  theme_bw() +
  guides(fill="none") +
  facet_wrap(~ harmony_clusters, scales="free_y", ncol=6) +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "cell_counts_per_sample.per_cluster.png"), width=12, height=11)

# umap

DimPlot(int_seu, reduction="umap.harmony", group.by= "harmony_clusters",
        raster = F) +
  labs(x="UMAP 1", y="UMAP 2")
ggsave(paste0(out_dir, "cluster_umap.png"), width=6, height=5)

DimPlot(int_seu, reduction="umap.harmony", group.by= "harmony_clusters",
        raster = F, label=T) +
  labs(x="UMAP 1", y="UMAP 2")
ggsave(paste0(out_dir, "cluster_umap.with_labels.png"), width=6, height=5)

DimPlot(int_seu, reduction="umap.harmony", group.by= "harmony_clusters",
        raster = F, split.by="orig.ident", ncol=3) +
  labs(x="UMAP 1", y="UMAP 2")
ggsave(paste0(out_dir, "cluster_umap.per_sample.png"), width=12, height=9)

DimPlot(int_seu, reduction="umap.harmony", group.by= "harmony_clusters",
        raster = F, split.by="orig.ident", ncol=3, label=T) +
  labs(x="UMAP 1", y="UMAP 2")
ggsave(paste0(out_dir, "cluster_umap.per_sample.with_labels.png"), width=12, height=9)




