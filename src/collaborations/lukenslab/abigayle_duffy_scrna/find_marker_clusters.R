library(Seurat)
library(SeuratObject)
library(ggplot2)
library(cowplot)
library(openxlsx)
library(tibble)
library(dplyr)

root_dir <- "~/Documents/projects/lukenslab/abigayle_duffy/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/find_marker_clusters/")
dir.create(out_dir, showWarnings = F)


# read in integrated seurat
int_seu <- LoadSeuratRds(paste0(root_dir,
                                "results/split_and_integrate/all_samples.integrated_seurat.RDS"))

# prepare data for using find markers

# set assay to RNA and join layers
DefaultAssay(int_seu) <- "RNA"
int_seu <- JoinLayers(int_seu)

int_seu <- NormalizeData(int_seu)
int_seu <- FindVariableFeatures(int_seu)
int_seu <- ScaleData(int_seu)

Idents(int_seu) <- "harmony_clusters"

known_markers <- list("Astrocyte"=c("Slc6a11", "Ntsr2","Htra1","Aqp4"),
                      "Microglia"=c("Cx3cr1","Ly86","P2ry12",
                                    "Gpr34","Hexb"))

dotplot_list <- lapply(names(known_markers), function(cluster_name) {
  
  markers <- known_markers[[cluster_name]]
  
  DotPlot(int_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) 
  
})

dotplot_list[[1]]
dotplot_list[[2]]

plot_grid(plotlist = dotplot_list, nrow = 1)
ggsave(paste0(out_dir, "astrocyte_microglia_markers.dotplot.png"), width=12, height=9, bg="white")

# subset and find agnostic markers for each cluster

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

# partially name the clusters

metadata <- int_seu@meta.data

metadata$cell_type_partial <- as.character(metadata$harmony_clusters)

metadata[metadata$harmony_clusters %in% c("6","7","10","16"),]$cell_type_partial <- "Astrocyte"
metadata[metadata$harmony_clusters %in% c("4","22","25","27"),]$cell_type_partial <- "Microglia"

int_seu$cell_type_partial <- metadata$cell_type_partial

DimPlot(int_seu, reduction="umap.harmony", 
        group.by= c("harmony_clusters","cell_type_partial"),
        label=T)
ggsave(paste0(out_dir, "harmony_clusters_v_celltype_partial.umap.png"),
       width=12, height=5)

SaveSeuratRds(int_seu, paste0(out_dir, "partially_celltype_named.seurat.RDS"))


