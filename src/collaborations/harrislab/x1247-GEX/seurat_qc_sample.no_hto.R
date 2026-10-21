library(Seurat)
library(SeuratObject)
library(ggplot2)
library(scDblFinder)
library(BiocParallel)

root_dir <- "~/Documents/projects/harrislab/x1247-GEX/"

out_dir <- paste0(root_dir, "results/seurat_qc_samples.no_hto/")
dir.create(out_dir, showWarnings = F)


# read in data

h5_file <- paste0(root_dir, "cellranger_outputs/filtered_feature_bc_matrix.h5")

raw_data <- Read10X_h5(h5_file)

seu_obj <- CreateSeuratObject(counts=raw_data, 
                              project="x1247-GEX", 
                              min.cells = 3, 
                              min.features = 200)

# filtering
seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^mt-")

VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha=0.4) 
ggsave(paste0(out_dir,"pre_filter_count_plot.png"), width=8, height=6)

# filter and replot
feature_min <- 200
feature_max <- 6000 # neuronal cells express MANY genes, so we can make this high
percent_mt_max <- 5

# make a data frame of 
filter_counts <- data.frame(type=c("Raw Cell Count",
                                   paste0("nFeature_RNA > ", feature_min),
                                   paste0("nFeature_RNA < ", feature_max),
                                   paste0("percent.mt < ", percent_mt_max)),
                              count=c(nrow(seu_obj@meta.data),
                                    sum(seu_obj@meta.data$nFeature_RNA > feature_min),
                                    sum(seu_obj@meta.data$nFeature_RNA < feature_max),
                                    sum(seu_obj@meta.data$percent.mt < percent_mt_max)))

# filter
seu_obj <- subset(seu_obj, 
                 subset=nFeature_RNA > feature_min &
                   nFeature_RNA < feature_max &
                   percent.mt < 5)

VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha=0.4) 
ggsave(paste0(out_dir, "post_filter_count_plot.png"), width=8, height=6)

filter_counts <- rbind(filter_counts,
                       data.frame(type="Filtered Count",
                                  count=nrow(seu_obj@meta.data)))

filter_counts$type <- factor(filter_counts$type,
                             levels=filter_counts$type)

ggplot(filter_counts,
       aes(x=type,
           y=count)) +
  geom_bar(stat="identity",
           color="black",
           fill="maroon") +
  geom_text(aes(label=count),
            vjust=-0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "filter_parameter_counts.png"), width=4, height=5)
    



# inital PCA for doublet detection

seu_obj <- SCTransform(seu_obj, vars.to.regress = "percent.mt")

seu_obj <- RunPCA(seu_obj)

ElbowPlot(seu_obj, ndims=30) + 
  scale_x_continuous(breaks=seq(0,30,5)) +
  scale_y_continuous(breaks=seq(0,30,5), limits=c(0,NA))
ggsave(paste0(out_dir, "pre_doublet_filtering.pca_elbow_plot.png"), width=6, height=5, bg="white")

max_pc_dim <- 20

# doublet finding

# find clusters, purely for doublet finding
seu_obj <- FindNeighbors(seu_obj, dims = 1:max_pc_dim)
seu_obj <- FindClusters(seu_obj)

sce <- as.SingleCellExperiment(seu_obj)

bp <- MulticoreParam(3, RNGseed=1234) # equivalent to set seed, for reproducibility
sce <- scDblFinder(sce, clusters="seurat_clusters", BPPARAM=bp)

seu_obj$scDblFinder.class <- sce$scDblFinder.class

VlnPlot(seu_obj, group.by="orig.ident", split.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
ggsave(paste0(out_dir, "doublet_qc_plots.png"), width=9, height=6)

doublet_counts <- as.data.frame(table(seu_obj$scDblFinder.class))

ggplot(doublet_counts, 
       aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="grey", color="black") +
  geom_text(aes(label=Freq), vjust=-0.4) +
  theme_bw() +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "doublet_counts.png"), width=4, height=5)

dbl_meta <- as.data.frame(colData(sce))

saveRDS(dbl_meta, file=paste0(out_dir, "doublet_output.RDS"))

seu_obj <- subset(seu_obj, scDblFinder.class == "singlet")

# redo PCA and determine clusters

seu_obj <- SCTransform(seu_obj, vars.to.regress = "percent.mt")

seu_obj <- RunPCA(seu_obj)

ElbowPlot(seu_obj, ndims=30) + 
  scale_x_continuous(breaks=seq(0,30,5)) +
  scale_y_continuous(breaks=seq(0,30,5), limits=c(0,NA))
ggsave(paste0(out_dir, "post_doublet_filtering.pca_elbow_plot.png"), width=6, height=5, bg="white")

max_pc_dim <- 20

seu_obj <- FindNeighbors(seu_obj, dims = 1:max_pc_dim)
seu_obj <- FindClusters(seu_obj)

seu_obj <- RunUMAP(seu_obj, dims=1:max_pc_dim)

DimPlot(seu_obj, reduction="umap", label=T)
ggsave(paste0(out_dir, "clusters.umap.png"), width=7, height=5)

# set assay to RNA and join layers
DefaultAssay(seu_obj) <- "RNA"
seu_obj <- JoinLayers(seu_obj)

seu_obj <- NormalizeData(seu_obj)
seu_obj <- FindVariableFeatures(seu_obj)
seu_obj <- ScaleData(seu_obj)


all_markers <- FindAllMarkers(seu_obj)

saveRDS(all_markers, paste0(out_dir, "all_markers.RDS"))

SaveSeuratRds(seu_obj, paste0(out_dir, "seurat.post_qc.RDS"))

