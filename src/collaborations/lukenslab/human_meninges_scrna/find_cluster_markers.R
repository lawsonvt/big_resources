library(Seurat)
library(SeuratObject)
library(ggplot2)
library(openxlsx)

root_dir <- "~/Documents/projects/lukenslab/ashley_bolte/human_meninges/"

out_dir <- paste0(root_dir, "results/find_cluster_markers/")
dir.create(out_dir, showWarnings = F)

seu_obj <- LoadSeuratRds(paste0(root_dir, "results/integrate_seurat_samples/all_samples.integrated_seurat.RDS"))

clusters <- sort(unique(seu_obj$harmony_clusters))

Idents(seu_obj) <- "harmony_clusters"

DefaultAssay(seu_obj) <- "RNA"

# scale and normalize the data
seu_obj <- NormalizeData(seu_obj)
seu_obj <- ScaleData(seu_obj)

seu_obj <- JoinLayers(seu_obj)


cluster_markers <- lapply(clusters, function(cluster) {
  
  print(paste0("Cluster ", cluster))
  
  markers <- FindMarkers(seu_obj,
                        ident.1 = cluster,
                        test.use="wilcox")
  markers$gene <- rownames(markers)
  markers$cluster <- cluster
  markers$pct_delta <- markers$pct.1 - markers$pct.2
  
  return(markers)
  
})
names(cluster_markers) <- paste0("Cluster ", clusters)

# filter it down to only positive ones

cluster_markers_f <- lapply(cluster_markers, function(data) {
  
  data <- data[data$avg_log2FC > 0,]
  
  data <- data[order(data$pct_delta, decreasing=T),]
  
  return(data)
})


# write to excel
write.xlsx(cluster_markers_f, file=paste0(out_dir, "cluster_markers.xlsx"),
           colWidths="auto")

# use previous markers

lab_cluster_markers <- list("B Cells"=c("Cd19", "Cd79a", "Ms4a1"),
                            "Monocytes"=c("Ccr2", "Cd44"),
                            "Neutrophils"=c("Ly6g", "Ngp", "Mmp8"),
                            "Macrophages"=c("Pf4", "Mrc1", "Ms4a7"),
                            "T Cells"=c("Trbc2", "Cd3d", "Lck"),
                            "Microglia"=c("Sall1", "Hexb", "P2ry12"),
                            "Proliferating Cells"=c("Mki67", "Ccnb1", "Tpx2"),
                            "Pericytes"=c("Pdgfrb", "Rgs5"))






