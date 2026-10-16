library(Seurat)
library(SeuratObject)
library(openxlsx)

root_dir <- "~/Documents/projects/harrislab/x1247-GEX/"

# test comment

out_dir <- paste0(root_dir, "results/hash_out_samples/")
dir.create(out_dir, showWarnings = F)

h5_file <- paste0(root_dir, "cellranger_outputs/filtered_feature_bc_matrix.h5")

sample_metadata <- read.xlsx(paste0(root_dir, "hashtag_meta.xlsx"))

hto_dir <- paste0(root_dir, "fastq_files/x1248-FB-ds.79d9b0df3b254687839612ddc0479427/hto_counts/umi_count/")

# read in h5

total_data <- Read10X_h5(h5_file)

hto <- Read10X(hto_dir, gene.column = 1)

dim(hto)
rownames(hto)
rowSums(hto)
head(colnames(hto))

colnames(hto) <- paste0(colnames(hto), "-1")
hto <- hto[rownames(hto) != "unmapped", ]

common <- intersect(colnames(total_data), colnames(hto))
length(common); ncol(total_data)

seu <- CreateSeuratObject(counts = total_data[, common])
seu[["HTO"]] <- CreateAssayObject(counts = hto[, common])

seu <- NormalizeData(seu, assay = "HTO", normalization.method = "CLR", margin = 2)
seu <- HTODemux(seu, assay = "HTO", positive.quantile = 0.99)


hto_counts <- GetAssayData(seu, assay = "HTO", layer = "counts") 

sum(colSums(hto_counts) == 0)
ncol(seu)


hto_raw <- Read10X(hto_dir, gene.column = 1)
dim(hto_raw)
rownames(hto_raw)
rowSums(hto_raw)
head(colnames(hto_raw))

# grep("hash|HTO|TotalSeq|B025", rownames(total_data), ignore.case = TRUE, value = TRUE)
# 
# h5ls(h5_file)
# table(h5read(h5_file, "matrix/features/feature_type"))
