library(Seurat)
library(SeuratObject)
library(ggplot2)
library(cowplot)
library(openxlsx)
library(tibble)
library(dplyr)
library(stringr)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/marker_dot_plots/")
dir.create(out_dir, showWarnings = F)


# read in integrated seurat
int_seu <- LoadSeuratRds(paste0(root_dir,
                                "results/integrate_seurat_samples/all_samples.integrated_seurat.RDS"))

# set assay to RNA and join layers
DefaultAssay(int_seu) <- "RNA"
int_seu <- JoinLayers(int_seu)

int_seu <- NormalizeData(int_seu)
int_seu <- FindVariableFeatures(int_seu)
int_seu <- ScaleData(int_seu)

Idents(int_seu) <- "harmony_clusters"

# read in given marker file
marker_file <- paste0(root_dir, "marker gene.xlsx")

marker_raw <- read.xlsx(marker_file)

marker_list <- apply(marker_raw, 2, function(column) {
  
  str_to_title(column[!is.na(column)])
  
})

# correct some gene names
marker_list$Erythocyte <- c("Hbb-bt",
                            "Hbb-bs")

# make dotplots
dotplot_list <- lapply(names(marker_list), function(cluster_name) {
  
  markers <- marker_list[[cluster_name]]
  
  DotPlot(int_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) 
  
})
names(dotplot_list) <- names(marker_list)

plot_grid(plotlist = dotplot_list, nrow = 3)
ggsave(paste0(out_dir, "marker_dot_plots.total.png"),
       width=26, height=17, bg="white")

dot_dir <- paste0(out_dir, "dotplots/")

dir.create(dot_dir, showWarnings = F)



# output to individ files
for (marker in names(dotplot_list)) {
  print(dotplot_list[[marker]])
  ggsave(paste0(dot_dir, marker, ".dot_plot.png"), width=5, height=8, bg="white")
}



gene_names <- rownames(int_seu)



