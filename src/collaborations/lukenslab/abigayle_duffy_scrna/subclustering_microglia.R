library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tibble)
library(scDblFinder)
library(BiocParallel)
library(openxlsx)
library(harmony)
library(gtools)
library(stringr)

root_dir <- "~/Documents/projects/lukenslab/abigayle_duffy/"

out_dir <- paste0(root_dir, "results/subclustering_microglia/")
dir.create(out_dir, showWarnings = F)

# read in integrated seurat
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/find_marker_clusters/partially_celltype_named.seurat.RDS"))

# clean up metadata
metadata <- seu_obj@meta.data

metadata$sex <- factor(gsub("[Pp0-9]+-", "", metadata$MULTI_ID))
metadata$age <- toupper(gsub("-(female|male)", "", metadata$MULTI_ID))
metadata$age <- factor(metadata$age, levels=mixedsort(unique(metadata$age)))

metadata$sample_id <- str_to_title(metadata$MULTI_ID)
metadata$sample_id <- factor(metadata$sample_id,
                             levels=mixedsort(unique(metadata$sample_id)))

# add to seurat object
seu_obj$sample_id <- metadata$sample_id
seu_obj$sex <- metadata$sex
seu_obj$age <- metadata$age

# subset down to microglia
#subset_seu <- subset(seu_obj, subset = cell_type_partial == "Microglia")
subset_seu <- subset(seu_obj, subset = harmony_clusters == "4")

# clean up and free memory
rm(seu_obj)
gc()

# determine the optimal number of dimensions through an elbow plot

subset_seu <- RunPCA(subset_seu, npcs = 50)

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="Microglia Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "microglia_subset.pre_doublet.pca_elbow_plot.png"), width=6, height=5, bg="white")

max_pc_dim <- 10

# cluster the harmonized data
subset_seu <- FindNeighbors(subset_seu, dims = 1:max_pc_dim, reduction = "pca")
subset_seu <- FindClusters(subset_seu, cluster.name = "microglia_clusters")


# set assay to RNA and join layers
DefaultAssay(subset_seu) <- "RNA"
subset_seu <- JoinLayers(subset_seu)

subset_seu <- NormalizeData(subset_seu)
subset_seu <- FindVariableFeatures(subset_seu)
subset_seu <- ScaleData(subset_seu)

# do some doublet finding
# run doublet finder
subset_sce <- as.SingleCellExperiment(subset_seu)

bp <- MulticoreParam(2, RNGseed=1234) # equivalent to set seed, for reproducibility
subset_sce <- scDblFinder(subset_sce, clusters="microglia_clusters", BPPARAM=bp)

# add doublet calls to seurat object

subset_seu$scDblFinder.class <- subset_sce$scDblFinder.class

VlnPlot(subset_seu, group.by="sample_id", split.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
ggsave(paste0(out_dir, "microglia.doublet_qc_plots.png"), width=9, height=6)

doublet_counts <- as.data.frame(table(subset_seu$scDblFinder.class))

ggplot(doublet_counts, 
       aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="grey", color="black") +
  geom_text(aes(label=Freq), vjust=-0.4) +
  theme_bw() +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "microglia.doublet_counts.png"), width=4, height=5)

dbl_meta <- as.data.frame(colData(subset_sce))

saveRDS(dbl_meta, file=paste0(out_dir, "microglia.doublet_output.RDS"))

# remove doublets, redo clustering
subset_seu <- subset(subset_seu, scDblFinder.class == "singlet")

# set Assay back to SCT for clustering
#DefaultAssay(subset_seu) <- "SCT"

subset_seu <- SCTransform(subset_seu, vars.to.regress = c("percent.mt"), verbose = F)

subset_seu <- RunPCA(subset_seu, npcs = 50)

subset_seu <- RunHarmony(subset_seu, group.by.vars="sample_id")

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="Microglia Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "microglia_subset.post_doublet.pca_elbow_plot.png"), width=6, height=5, bg="white")

max_pc_dim <- 20

# cluster the harmonized data
subset_seu <- FindNeighbors(subset_seu, dims = 1:max_pc_dim, reduction = "pca")
subset_seu <- FindClusters(subset_seu, cluster.name = "microglia_clusters", 
                           resolution = 0.2)

# create umap
subset_seu <- RunUMAP(subset_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.microglia_pca")


DimPlot(subset_seu, reduction="umap.microglia_pca", group.by= "microglia_clusters",
        label=T) 
ggsave(paste0(out_dir, "microglia_subcluster.umap.png"), width=7, height=5)

DimPlot(subset_seu, reduction="umap.microglia_pca", group.by= "microglia_clusters",
        label=T, split.by = "sample_id", ncol=2) 



DimPlot(subset_seu, reduction="umap.microglia_pca", group.by= "harmony_clusters",
        label=T, split.by = "sample_id", ncol=2) 



DimPlot(subset_seu, reduction="umap.microglia_pca", group.by= "microglia_clusters",
        label=T, split.by = "age", ncol=3) 

DimPlot(subset_seu, reduction="umap.harmony", group.by= "harmony_clusters",
        label=T, split.by = "sample_id", ncol=2) 

DimPlot(subset_seu, reduction="umap.harmony", group.by= "harmony_clusters",
        label=T) 


subset_metadata <- subset_seu@meta.data

table(subset_metadata$sample_id,
      subset_metadata$microglia_clusters)


