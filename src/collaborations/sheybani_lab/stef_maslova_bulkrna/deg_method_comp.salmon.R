library(ggplot2)
library(ggrepel)
library(snakecase)

root_dir <- "~/Documents/projects/sheybanilab/stef_maslova/"

out_dir <- paste0(root_dir, "results/deg_method_comp.salmon/")
dir.create(out_dir, showWarnings = F)

# load in data sets
total_results <- readRDS(paste0(root_dir, "results/salmon_deg_analysis/deg_results.RDS"))
direct_results <- readRDS(paste0(root_dir, "results/salmon_deg_analysis.direct_comparison/deg_results.RDS"))

lapply(total_results, nrow)
lapply(direct_results, nrow)

# merge em up
merged_results <- lapply(names(total_results), function(contrast) {
  
  merged <- merge(total_results[[contrast]],
                  direct_results[[contrast]],
                  by=c("gene_id","gene_symbol","contrast"),
                  suffixes=c(".total",".direct"),
                  all=T)
  
  return(merged)
  
})

min_fdr <- 0.1

# make scatter plots for each comparison
for (result in merged_results) {
  
  result$sig <- factor(as.numeric(result$padj.total < 0.1) +
    as.numeric(result$padj.direct < 0.1) * 2)
  levels(result$sig) <- c("No DEG",
                          "DEG in Total",
                          "DEG in Direct",
                          "DEG in both")
  
  result <- result[!is.na(result$sig),]
  
  result <- result[order(result$sig),]
  
  result_top <- result[result$sig != "No DEG",]
  
  comparison <- unique(result$contrast)
  
  ggplot(result,
         aes(x=log2FoldChange.total,
             y=log2FoldChange.direct,
             color=sig)) +
    geom_point(alpha=0.5, size=2) +
    geom_hline(yintercept = 0, linetype=2, color="black") +
    geom_vline(xintercept = 0, linetype=2, color="black") +
    geom_text_repel(data=result_top,
                    aes(label=gene_symbol),
                    size=3,
                    max.overlaps = 50) +
    scale_color_manual(values=c("grey","red","blue","purple")) +
    theme_bw() +
    labs(x="Log2FC in Total Analysis",
         y="Log2FC in Direct Comparison Analysis",
         color=NULL,
         title=comparison) 
  ggsave(paste0(out_dir, to_snake_case(comparison), ".scatter_fc.png"), width=7, height=5)
  
  
} 



