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

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/x1247-GEX/"

out_dir <- paste0(root_dir, "results/microglia_subcluster.diff_exp_analysis/")
dir.create(out_dir, showWarnings = F)

# read in Seurat object

seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/microglia.cell_type_naming.with_hto/seurat.microglia_cell_named.RDS"))

metadata <- seu_obj@meta.data

# just do all cells
cells <- unique(metadata$cell_name)

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
  group.by = c("cell_name","sample_id")
)

# Extract the count matrix
pseudobulk_matrix <- seu_pseudo$RNA

# Extract unique sample-level metadata
sample_metadata <- metadata %>%
  select(sample_id, condition) %>%
  distinct() %>%
  as.data.frame()

rownames(sample_metadata) <- sample_metadata$sample

# ST-5X - MG-5X

condition1 <- "ST-5X"
condition2 <- "MG-5X"

cell_de_results <- lapply(cells, function(cell) {
  
  cat("Processing:", cell, "\n")
  
  cell_cols <- grep(paste0("^", cell, "_"), 
                    colnames(pseudobulk_matrix), 
                    value = TRUE)
  
  
  # Extract sample names from column names
  sample_names <- sub(paste0("^", cell, "_"), "", cell_cols)
  
  # Subset counts
  counts <- pseudobulk_matrix[, cell_cols]
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
    mutate(cell = cell,
           comparison = paste0(condition1, " - ", condition2))
  
  # filter out NAs
  res_df <- res_df[!is.na(res_df$padj),]
  
  return(res_df)
})
names(cell_de_results) <- cells

#lapply(cell_de_results, head)

write.xlsx(cell_de_results, 
           file=paste0(out_dir, "microglia_cells.de_results.xlsx"),
           colWidths="auto")

saveRDS(cell_de_results,
        file=paste0(out_dir, "microglia_cells.de_results.RDS"))

# Make Volcano plots -----------------------------------------------------------

volcano_dir <- paste0(out_dir, "volcano_plots.cell_types/")

dir.create(volcano_dir, showWarnings = F)


top_genes <- 25

volcano_plot_list <-  lapply(names(cell_de_results), function(cell) {
  
  print(cell)
  
  subset <- cell_de_results[[cell]]
  
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
    labs(x="Log2 Fold Change", y="-log10(P-Value)", title=cell,
         subtitle=unique(subset$comparison))
  p1
  ggsave(paste0(volcano_dir, to_snake_case(cell), ".volcano_plot.png"), width=5, height=6)
  
  return(p1)
})


plot_grid(plotlist = volcano_plot_list, ncol=3)
ggsave(paste0(out_dir, "cell_deg.volcanoes.png"), width=12, height=8, bg="white")



