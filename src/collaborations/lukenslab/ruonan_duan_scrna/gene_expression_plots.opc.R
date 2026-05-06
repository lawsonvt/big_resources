library(Seurat)
library(ggplot2)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/gene_expression_plots.opc/")
dir.create(out_dir, showWarnings = F)


# read in integrated seurat
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/celltype_naming/all_samples.celltype_named.seurat.RDS"))

# group the conditions
metadata <- seu_obj@meta.data

metadata$sample <- factor(metadata$orig.ident,
                          levels=c("X1161_1NEG",
                                   "X1176-WT1",
                                   "X1177-WT2",
                                   "X1162_2Het",
                                   "X1178-KO1",
                                   "X1179-KO2"))

metadata$condition <- "WT"

metadata[metadata$sample %in% c("X1162_2Het",
                                "X1178-KO1",
                                "X1179-KO2"),]$condition <- "KO"
metadata$condition <- factor(metadata$condition,
                             levels=c("WT","KO"))

# add to seurat object
seu_obj$sample <- metadata$sample
seu_obj$condition <- metadata$condition

# subset down to tcell
subset_seu <- subset(seu_obj, subset = cell_type %in% c("OPC",
                                                        "Mature OL",
                                                        "Newly Formed OL"))

# set assay to RNA and join layers
DefaultAssay(subset_seu) <- "RNA"
subset_seu <- JoinLayers(subset_seu)

subset_seu <- NormalizeData(subset_seu)
subset_seu <- FindVariableFeatures(subset_seu)
subset_seu <- ScaleData(subset_seu)

genelist <- c("Gpx4",
              "Cd74",
              "H2-D1",
              "H2-K1")

for (gene in genelist) {

  VlnPlot(subset_seu,
          features = gene,
          group.by = "cell_type",
          split.by = "condition") +
    labs(x=NULL)
  ggsave(paste0(out_dir, gene, ".opc_ol.violin_plot.png"), width=6, height=4)
  
}







