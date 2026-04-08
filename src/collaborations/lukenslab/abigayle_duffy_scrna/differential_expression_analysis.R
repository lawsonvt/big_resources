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
library(stringr)

root_dir <- "~/Documents/projects/lukenslab/abigayle_duffy/"

out_dir <- paste0(root_dir, "results/differential_expression_analysis/")
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

# get the cell types

cells <- c("Astrocyte",
           "Microglia")

seu_pseudo <- AggregateExpression(
  seu_obj,
  assays = "RNA",
  slot = "counts",
  group.by = c("cell_type_partial","sample_id")
)

# Extract the count matrix
pseudobulk_matrix <- seu_pseudo$RNA

# Extract unique sample-level metadata
sample_metadata <- metadata %>%
  select(sample_id, sex, age) %>%
  distinct() %>%
  as.data.frame()

rownames(sample_metadata) <- sample_metadata$sample_id

# define contrasts
contrasts <- list("P14-P7"=c("age","P14","P7"),
                  "P40-P14"=c("age","P40","P14"),
                  "P40-P7"=c("age","P40","P7"),
                  "female-male"=c("sex","female","male"))

de_results <- lapply(cells, function(cell) {
  
  cat("Processing:", cell, "\n")
  
  celltype_cols <- grep(paste0("^", cell, "_"), 
                        colnames(pseudobulk_matrix), 
                        value = TRUE)
  
  # Extract sample names from column names
  sample_names <- sub(paste0("^", cell, "_"), "", celltype_cols)
  
  # Subset counts
  counts <- pseudobulk_matrix[, celltype_cols]
  colnames(counts) <- sample_names
  
  meta <- sample_metadata[sample_names, , drop = FALSE]
  
  # Keep genes with at least 10 counts in at least 2 samples
  keep <- rowSums(counts >= 10) >= 2
  counts <- counts[keep, ]
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ age + sex
  )
  
  dds <- DESeq(dds)
  
  results_list <- lapply(contrasts, function(contrast) {
    
    # calculate results
    res <- results(dds, contrast=contrast) 
    
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene") %>%
      arrange(padj) %>%
      mutate(cell_type = cell,
             comparison = paste0(contrast[2], "-", contrast[3]))
    
    # filter out NAs
    res_df <- res_df[!is.na(res_df$padj),]
    
    return(res_df)
  })
  names(results_list) <- names(contrasts)
  
  return(results_list)
})
names(de_results) <- cells

saveRDS(de_results,
        file=paste0(out_dir, "cell_types.age_sex_contrasts.de_results.RDS"))


# output to excel
for (cell in cells) {
  
  results_list <- de_results[[cell]]
  
  write.xlsx(results_list, file=paste0(out_dir, to_snake_case(cell), ".degs.xlsx"),
             colWidths="auto")
  
}

# Make Volcano plots -----------------------------------------------------------

volcano_dir <- paste0(out_dir, "volcano_plots.cell_types/")

dir.create(volcano_dir, showWarnings = F)


top_genes <- 25

volcano_plot_list <-  lapply(names(de_results), function(cell) {
  
  results_list <- de_results[[cell]]
  
  lapply(names(results_list), function(contrast) {
    
    subset <- results_list[[contrast]]
  
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
           subtitle=contrast)
    p1
    ggsave(paste0(volcano_dir, to_snake_case(cell), ".",
                  to_snake_case(contrast), ".volcano_plot.png"), width=5, height=6)
    
    return(p1)
  })
})

plot_grid(plotlist = c(volcano_plot_list[[1]],
                       volcano_plot_list[[2]]), nrow = 2)
ggsave(paste0(out_dir, "cell_types.volcanoes.png"), width=12, height=8, bg="white")

