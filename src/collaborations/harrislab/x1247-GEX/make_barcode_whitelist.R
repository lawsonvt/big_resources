library(Seurat)

root_dir <- "~/Documents/projects/harrislab/x1247-GEX/"

# test comment

out_dir <- paste0(root_dir, "results/make_barcode_whitelist/")
dir.create(out_dir, showWarnings = F, recursive = T)

gex <- Read10X_h5(paste0(root_dir, "cellranger_outputs/filtered_feature_bc_matrix.h5"))
barcodes <- colnames(gex)

head(barcodes)   # e.g. "AAACCTGAGAAACCAT-1"

# Strip the "-1" suffix for CITE-seq-Count
barcodes_nosuffix <- sub("-1$", "", barcodes)

writeLines(barcodes_nosuffix, paste0(out_dir, "whitelist.txt"))
