library(Seurat)
library(SeuratObject)
library(DESeq2)
library(plyr)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(snakecase)
library(cowplot)
library(gtools)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/immune_cell_split/immune_cell.filtered.diff_exp_analysis/")
dir.create(out_dir, showWarnings = F, recursive = T)

# read in nonimmune seurat
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/immune_cell_split/immune_cell_subclustering/subset_immune.seurat.RDS"))

# filter down samples
samples <- levels(seu_obj$sample)
samples <- grep("(NEG)|(Het)", samples, value=T, invert=T)

seu_obj <- subset(seu_obj, subset = sample %in% samples)

metadata <- seu_obj@meta.data

# add in cluster word to clusters
levels(metadata$immune_clusters) <- paste0("cluster", levels(metadata$immune_clusters))

seu_obj$immune_clusters <- metadata$immune_clusters

# clusters to analyze
# clusters <- c("cluster2",
#               "cluster10",
#               "cluster18",
#               "cluster19",
#               "cluster20")

# just do all of them
clusters <- levels(seu_obj$immune_clusters)[1:20]

# set assay to RNA and join layers
DefaultAssay(seu_obj) <- "RNA"
seu_obj <- JoinLayers(seu_obj)

seu_obj <- NormalizeData(seu_obj)
seu_obj <- FindVariableFeatures(seu_obj)
seu_obj <- ScaleData(seu_obj)

seu_pseudo <- AggregateExpression(
  seu_obj,
  assays = "RNA",
  slot = "counts",
  group.by = c("immune_clusters","sample")
)

# Extract the count matrix
pseudobulk_matrix <- seu_pseudo$RNA

# fix sample names
colnames(pseudobulk_matrix) <- gsub("\\-1NEG", "_1NEG", colnames(pseudobulk_matrix))
colnames(pseudobulk_matrix) <- gsub("\\-2Het", "_2Het", colnames(pseudobulk_matrix))

# Extract unique sample-level metadata
sample_metadata <- metadata %>%
  select(sample, condition) %>%
  distinct() %>%
  as.data.frame()

rownames(sample_metadata) <- sample_metadata$sample

# WT v KO

condition1 <- "KO"
condition2 <- "WT"

ko_minus_wt.de_results <- lapply(clusters, function(cluster) {
  
  cat("Processing:", cluster, "\n")
  
  cluster_cols <- grep(paste0("^", cluster, "_"), 
                        colnames(pseudobulk_matrix), 
                        value = TRUE)
  
  
  # Extract sample names from column names
  sample_names <- sub(paste0("^", cluster, "_"), "", cluster_cols)
  
  # Subset counts
  counts <- pseudobulk_matrix[, cluster_cols]
  colnames(counts) <- sample_names
  
  meta <- sample_metadata[sample_names, , drop = FALSE]
  
  # Filter to only the two conditions being compared
  samples_to_keep <- meta$condition %in% c(condition1, condition2)
  counts <- counts[, samples_to_keep, drop = FALSE]
  meta <- meta[samples_to_keep, , drop = FALSE]
  
  # Keep genes with at least 10 counts in at least 2 samples
  keep <- rowSums(counts >= 10) >= 2
  counts <- counts[keep, ]
  
  meta$condition <- factor(as.character(meta$condition), levels = c(condition1, condition2))
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ condition
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds, contrast = c("condition", condition1, condition2))
  
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(pvalue) %>%
    mutate(cluster = cluster,
           comparison = paste0(condition1, "_m_", condition2))
  
  # filter out NAs
  res_df <- res_df[!is.na(res_df$padj),]
  
  return(res_df)
})
names(ko_minus_wt.de_results) <- clusters

write.xlsx(ko_minus_wt.de_results, 
           file=paste0(out_dir, "ko_minus_wt.clusters.de_results.xlsx"),
           colWidths="auto")

saveRDS(ko_minus_wt.de_results,
        file=paste0(out_dir, "ko_minus_wt.clusters.de_results.RDS"))

# Make Volcano plots -----------------------------------------------------------

volcano_dir <- paste0(out_dir, "volcano_plots.clusters/")

dir.create(volcano_dir, showWarnings = F)


top_genes <- 25

volcano_plot_list <-  lapply(names(ko_minus_wt.de_results), function(cluster) {
  
  print(cluster)
  
  subset <- ko_minus_wt.de_results[[cluster]]
  
  subset$log_p <- -log10(subset$pvalue)
  
  subset_sig <- subset[subset$padj < 0.05 &
                         abs(subset$log2FoldChange) > 0.5,]
  
  subset_top <- subset_sig[1:top_genes,]
  
  if (any(is.na(subset_top$gene))) {
    subset_top <- subset_top[!is.na(subset_top$gene),]
  }
  
  logp_thresh <- min(subset_sig$log_p)
  
  p1 <- ggplot(subset,
               aes(x=log2FoldChange,
                   y=log_p)) +
    geom_point(alpha=0.4, color="black") +
    geom_hline(yintercept = logp_thresh,
               color="red", linetype=2) +
    geom_vline(xintercept = 0.5,
               color="red", linetype=2) +
    geom_vline(xintercept = -0.5,
               color="red", linetype=2) +
    geom_point(data=subset_sig,
               color="red") +
    geom_text_repel(data=subset_top,
                    aes(label=gene),
                    color="red", size=2.5,
                    max.overlaps = 50) +
    theme_bw() +
    labs(x="Log2 Fold Change", y="-log10(P-Value)", title=cluster,
         subtitle="KO - WT")
  p1
  ggsave(paste0(volcano_dir, to_snake_case(cluster), ".volcano_plot.png"), width=5, height=6)
  
  return(p1)
})


plot_grid(plotlist = volcano_plot_list, nrow = 4)
ggsave(paste0(out_dir, "cluster_deg.volcanoes.png"), width=14, height=11, bg="white")




