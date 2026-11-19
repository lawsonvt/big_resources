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

root_dir <- "/Users/mjl3p/Documents/projects/gaultierlab/shreyas_gujar/snRNASeq/"
#root_dir <- "/home/mjl3p/projects/gaultierlab/shreyas_gujar/snRNASeq/"

out_dir <- paste0(root_dir, "results/deg_analysis.cells/")
dir.create(out_dir, showWarnings = F)

seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/cell_type_naming/seurat.cell_named.RDS"))

metadata <- seu_obj@meta.data

# ggplot(metadata,
#        aes(x=condition,
#            y=percent.mt)) +
#   geom_boxplot(outlier.shape = NA) +
#   #geom_jitter(width=0.1) +
#   theme_bw() +
#   facet_wrap(~ harmony_clusters, ncol=8, scales = "free_y")


# just do all of them
cells <- unique(seu_obj$cell_name)

# set assay to RNA and join layers
DefaultAssay(seu_obj) <- "RNA"
seu_obj <- JoinLayers(seu_obj)

seu_obj <- NormalizeData(seu_obj)
seu_obj <- FindVariableFeatures(seu_obj)
seu_obj <- ScaleData(seu_obj)

# aggregate data
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

# specify conditions / order
condition1 <- "S"
condition2 <- "O"

de_results <- lapply(cells, function(cell) {
  
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
names(de_results) <- cells

de_results_out <- de_results

names(de_results_out) <- substr(names(de_results_out), 1, 31)

# save the results
write.xlsx(de_results_out, 
           file=paste0(out_dir, "total_cell.de_results.xlsx"),
           colWidths="auto")

saveRDS(de_results,
        file=paste0(out_dir, "total_cell.de_results.RDS"))

# Make Volcano plots -----------------------------------------------------------

volcano_dir <- paste0(out_dir, "volcano_plots.cell_types/")

dir.create(volcano_dir, showWarnings = F)


top_genes <- 25

volcano_plot_list <-  lapply(names(de_results), function(cell) {
  
  print(cell)
  
  subset <- de_results[[cell]]
  
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

plot_grid(plotlist = volcano_plot_list, ncol = 5)
ggsave(paste0(out_dir, "total_volcanoes.png"), width=20, height=15)

# # inspect the results
# 
# de_results_df <- bind_rows(de_results)
# 
# degs <- de_results_df[de_results_df$padj < 0.05 &
#                         abs(de_results_df$log2FoldChange) > 0.5,]
# 
# 
# # the mitochondrial genes
# mt_degs <- degs[grepl("^mt", degs$gene),]
# 
# cat(unique(mt_degs$gene), sep="\n")
# 
# cluster_cell_xref <- unique(metadata[,c("harmony_clusters","cell_name")])
# cluster_cell_xref <- cluster_cell_xref[order(cluster_cell_xref$harmony_clusters),]
# 
# View(cluster_cell_xref[cluster_cell_xref$harmony_clusters %in% mt_degs$cluster,])
# 
# 
# VlnPlot(seu_obj, features="percent.mt", split.by="condition")
# 
# ggplot(metadata,
#        aes(x=condition,
#            y=percent.mt)) +
#   geom_boxplot(outlier.shape = NA) +
#   #geom_jitter(width=0.1) +
#   theme_bw() +
#   scale_y_log10() +
#   facet_wrap(~ harmony_clusters, ncol=8, scales = "free_y")
# 
# hist(metadata$percent.mt)
# 
# nrow(metadata[metadata$percent.mt < 1,])
# 
# # plot out the mitochondrial effect size vs baseline (Claude is helping with this)
# 
# mt_genes <- unique(mt_degs$gene)
# 
# # 1. Baseline mt per cluster, from controls only (average the 3 samples equally)
# baseline <- seu_obj@meta.data %>%
#   filter(condition == "O") %>%
#   group_by(harmony_clusters, sample_id) %>%
#   summarise(mt = mean(percent.mt), .groups = "drop") %>%
#   group_by(harmony_clusters) %>%
#   summarise(baseline_mt = mean(mt))
# 
# # 2. Treatment effect on mt-genes, per cluster (all 37, significant or not)
# effect <- bind_rows(lapply(names(de_results), function(cl) {
#   r <- as.data.frame(de_results[[cl]])
#   r <- r[r$gene %in% mt_genes, ]
#   data.frame(harmony_clusters = cl,
#              mt_lfc  = median(r$log2FoldChange, na.rm = TRUE))
# }))
# 
# # 3. Nuclei per cluster, for point size
# n_nuc <- seu_obj@meta.data %>% count(harmony_clusters, name = "n_nuclei")
# 
# # 4. Absolute shift in percent.mt (treatment - control), as a second view
# delta <- seu_obj@meta.data %>%
#   group_by(harmony_clusters, condition, sample_id) %>%
#   summarise(mt = mean(percent.mt), .groups = "drop") %>%
#   group_by(harmony_clusters, condition) %>% summarise(mt = mean(mt), .groups = "drop") %>%
#   tidyr::pivot_wider(names_from = condition, values_from = mt) %>%
#   mutate(delta_mt = S - O) %>%
#   select(harmony_clusters, delta_mt)
# 
# df <- effect %>%
#   left_join(baseline, by = "harmony_clusters") %>%
#   left_join(n_nuc,   by = "harmony_clusters") %>%
#   left_join(delta,   by = "harmony_clusters")
# 
# # add in mt sig cells
# df$mt_degs <- df$harmony_clusters %in% mt_degs$cluster
# 
# # Plot A: log2FC vs baseline
# ggplot(df, aes(baseline_mt, mt_lfc)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_point(aes(size = n_nuclei, color=mt_degs)) +
#   geom_smooth(method = "lm", se = TRUE) +
#   geom_text_repel(aes(label = harmony_clusters), size = 3) +
#   labs(x = "Baseline mean percent.mt (controls)",
#        y = "Median mt-gene log2FC (treatment vs control)") +
#   theme_bw()
# 
# # Plot B: absolute shift in percent.mt vs baseline
# ggplot(df, aes(baseline_mt, delta_mt)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_point(aes(size = n_nuclei, color=mt_degs)) +
#   geom_smooth(method = "lm") +
#   geom_text_repel(aes(label = harmony_clusters), size = 3) +
#   labs(x = "Baseline mean percent.mt (controls)",
#        y = "Δ percent.mt (treatment - control)") +
#   theme_bw()
# 
# cor.test(df$baseline_mt, df$mt_lfc,  method = "spearman")
# cor.test(df$baseline_mt, df$delta_mt, method = "spearman")
# 
# oxphos_pat <- "^(Ndufa|Ndufb|Ndufs|Ndufv|Cox4|Cox5|Cox6|Cox7|Atp5|Uqcr)"
# ox <- bind_rows(lapply(names(de_results), function(cl) {
#   r <- as.data.frame(de_results[[cl]])
#   r <- r[grepl(oxphos_pat, r$gene), ]
#   data.frame(harmony_clusters = cl, ox_lfc = median(r$log2FoldChange, na.rm = TRUE))
# }))
# df <- left_join(df, ox, by = "harmony_clusters")
# sum(df$ox_lfc < 0, na.rm = T); cor.test(df$mt_lfc, df$ox_lfc, method = "spearman")
# 
# summary(df$mt_lfc); summary(df$ox_lfc)
# wilcox.test(df$ox_lfc, mu = 0)
# 
# ribo_pat <- "^(Rpl|Rps)"
# cyto <- bind_rows(lapply(names(de_results), function(cl) {
#   r <- as.data.frame(de_results[[cl]])
#   data.frame(cluster = cl,
#              ribo_lfc   = median(r$log2FoldChange[grepl(ribo_pat, r$gene)], na.rm = TRUE),
#              malat1_lfc = r[r$gene == "Malat1",]$log2FoldChange)
# }))
# colMeans(cyto[, -1], na.rm = TRUE)
# 
# 
