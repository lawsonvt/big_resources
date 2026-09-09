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

root_dir <- "/Users/mjl3p/Documents/projects/lukenslab/ashley_bolte/human_meninges/"

out_dir <- paste0(root_dir, "results/cell_category.diff_exp_analysis/")
dir.create(out_dir, showWarnings = F)

# read in nonimmune seurat
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/annotate_cell_names/seurat.named.RDS"))

metadata <- seu_obj@meta.data

# find clusters that combine more than one cell type
cell_cats <- unique(metadata[,c("final_cell_name","cell_category")])

dupe_cats <- unique(cell_cats[duplicated(cell_cats$cell_category),]$cell_category)

# pull out cells
cells <- dupe_cats

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
  group.by = c("cell_category","Sample_name")
)

# Extract the count matrix
pseudobulk_matrix <- seu_pseudo$RNA

# Extract unique sample-level metadata
sample_metadata <- metadata %>%
  select(Sample_name, condition) %>%
  distinct() %>%
  as.data.frame()

rownames(sample_metadata) <- sample_metadata$Sample_name

contrasts <- list("trauma_aged-control_aged"=c("condition","chronic_trauma_aged","control_aged"),
                  "trauma_aged-control_young"=c("condition","chronic_trauma_aged","control_young"),
                  "control_aged-control_young"=c("condition","control_aged","control_young"))

cell_results_list <- lapply(cells, function(cell) {
  
  print(cell)
  
  # filter down to correct cells
  cell_cols <- grep(paste0("^", cell, "_"), 
                       colnames(pseudobulk_matrix), 
                       value = TRUE)
  
  
  # Extract sample names from column names
  sample_names <- sub(paste0("^", cell, "_"), "", cell_cols)
  
  # Subset counts
  counts <- pseudobulk_matrix[, cell_cols]
  colnames(counts) <- sample_names
  
  meta <- sample_metadata[sample_names, , drop = FALSE]
  
  # Keep genes with at least 10 counts in at least 2 samples
  keep <- rowSums(counts >= 10) >= 2
  counts <- counts[keep, ]
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ condition
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
 
  contrast_results_list <- lapply(names(contrasts), function(contrast_name) {
    
    print(contrast_name)
    
    contrast <- contrasts[[contrast_name]]
    
    # Get results
    res <- results(dds, contrast = contrast)
    
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene") %>%
      arrange(pvalue) %>%
      mutate(cell = cell,
             comparison = contrast_name)
    
    # filter out NAs
    res_df <- res_df[!is.na(res_df$padj),]
    
    return(res_df)
    
  })
  names(contrast_results_list) <- names(contrasts)
  
  return(contrast_results_list)
   
})
names(cell_results_list) <- cells

cell_results_df <- bind_rows(lapply(cell_results_list, bind_rows))

table(cell_results_df$cell, cell_results_df$comparison)

sig_cell_results_df <- cell_results_df[cell_results_df$padj < 0.05 &
                                         abs(cell_results_df$log2FoldChange) > 0.5,]
table(sig_cell_results_df$cell, sig_cell_results_df$comparison)

up_sig_cell_results_df <- sig_cell_results_df[sig_cell_results_df$log2FoldChange > 0,]
down_sig_cell_results_df <- sig_cell_results_df[sig_cell_results_df$log2FoldChange < 0,]

sig_cell_results_df$type <- "Up-regulated"
sig_cell_results_df[sig_cell_results_df$log2FoldChange < 0,]$type <- "Down-regulated"

ggplot(sig_cell_results_df,
       aes(x=comparison,
           fill=type)) +
  geom_bar(color="black") +
  scale_fill_manual(values=c("dodgerblue","firebrick1")) +
  theme_bw() +
  facet_wrap(~ cell, ncol=3, scales="free_y") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.position = "bottom") +
  labs(x=NULL, y="DEG Count at adjusted p-value < 0.05", fill=NULL)
ggsave(paste0(out_dir, "deg_counts.bar_plots.png"), width=10, height=6)

# output excel files
for (cell in cells) {
  
  out <- cell_results_list[[cell]]
  
  write.xlsx(out, paste0(out_dir, to_snake_case(cell), ".deg_results.xlsx"),
             colWidths="auto")
  
}

# volcano plots

# one plot per cell

top_genes <- 25

for (cell_results in cell_results_list) {
  
  plot_list <- lapply(cell_results, function(data) {
  
    cell <<- unique(data$cell)
    
    data$log_p <- -log10(data$pvalue)
    
    sig_data <- data[data$padj < 0.05 &
                       abs(data$log2FoldChange) > 0.5,]
    
    sig_top <- sig_data[1:top_genes,]
    
    if (any(is.na(sig_top$gene))) {
      sig_top <- sig_top[!is.na(sig_top$gene),]
    }
    
    logp_thresh <- min(sig_data$log_p)
    
    comparison <- unique(data$comparison)
    
    p1 <- ggplot(data,
                 aes(x=log2FoldChange,
                     y=log_p)) +
      geom_point(alpha=0.4, color="black") +
      geom_hline(yintercept = logp_thresh,
                 color="red", linetype=2) +
      geom_vline(xintercept = 0.5,
                 color="red", linetype=2) +
      geom_vline(xintercept = -0.5,
                 color="red", linetype=2) +
      geom_point(data=sig_data,
                 color="red") +
      geom_text_repel(data=sig_top,
                      aes(label=gene),
                      color="red", size=2.5,
                      max.overlaps = 50) +
      theme_bw() +
      labs(x="Log2 Fold Change", y="-log10(P-Value)", title=comparison)
    p1
  })
  plot_grid(plotlist = plot_list, nrow=1)
  ggsave(paste0(out_dir, to_snake_case(cell), ".volcano_plot.png"), width=11, height=5)
  
}

# output results
saveRDS(cell_results_list, file=paste0(out_dir, "deg_results.RDS"))






