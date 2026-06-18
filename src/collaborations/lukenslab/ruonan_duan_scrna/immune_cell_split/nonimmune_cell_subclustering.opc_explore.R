library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tibble)
library(openxlsx)
library(cowplot)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/immune_cell_split/nonimmune_cell_subclustering.opc_explore/")
dir.create(out_dir, showWarnings = F, recursive = T)

# read in nonimmune seurat
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/immune_cell_split/nonimmune_cell_subclustering/subset_nonimmune.seurat.RDS"))

# from here http://xteam.xbio.top/CellMarker/search.jsp?species=Mouse&tissue=Brain&cellname=Progenitor%20cell
opc_markers <- c("Pdgfra",
                 "Cacng4",
                 "Epn2",
                 "Megf11",
                 "Myt1",
                 "Nxph1",
                 "Pcdh15",
                 "Ppfibp1",
                 "Tnr",
                 "Vcan")

DotPlot(seu_obj, features = opc_markers) + 
  RotatedAxis() + 
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "opc_markers.dot_plot.png"), width=7, height=6, bg="white")

samples <- levels(seu_obj$sample)

dotplot_list <- lapply(samples, function(sample_id) {
  
  p <- DotPlot(subset(seu_obj, subset = sample == sample_id), 
               features = opc_markers) + 
    RotatedAxis() + 
    labs(x=NULL, y=NULL, title=sample_id)
  
  return(p)
  
})

plot_grid(plotlist = dotplot_list, nrow = 2)
ggsave(paste0(out_dir, "opc_markers.per_sample.all_clusters.dot_plot.png"), width=16, height=10, bg="white")

# focus just on the OPC clusters
opc_clusters <- c("2","10","18","19","20")

dotplot_list <- lapply(samples, function(sample_id) {
  
  p <- DotPlot(subset(seu_obj, subset = sample == sample_id,), 
               features = opc_markers,
               idents=opc_clusters,
               dot.min = 0) + 
    RotatedAxis() + 
    labs(x=NULL, y=NULL, title=sample_id)
  
  return(p)
  
})

plot_grid(plotlist = dotplot_list, nrow = 2)

