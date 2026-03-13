library(openxlsx)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(stringr)
library(DESeq2)
library(snakecase)

# output folder
out_dir <- "~/Documents/projects/gaultierlab/sam_wachamo/bulkRNASeq/results/pathway_plots/"

dir.create(out_dir, showWarnings = F)

# read in data
deg_results <- readRDS("~/Documents/projects/gaultierlab/sam_wachamo/bulkRNASeq/results/deg_workflow/IAA_v_Vehicle.degs.RDS")
gsea_results <- readRDS("~/Documents/projects/gaultierlab/sam_wachamo/bulkRNASeq/results/gsea_workflow/IAA_v_Vehicle.gsea_pathway_results.RDS")
dds <- readRDS("~/Documents/projects/gaultierlab/sam_wachamo/bulkRNASeq/results/deg_workflow/dds.RDS")

gene_xref <- deg_results[,c("gene_id","gene_symbol")]
gene_xref$gene_id <- gsub("\\.[0-9]+", "", gene_xref$gene_id)
rownames(gene_xref) <- gene_xref$gene_id

# find the pathways of interest

gsea_results_df <- as.data.frame(bind_rows(gsea_results))

grep("MYELINATION", gsea_results_df$pathway, value=T)
grep("ARYL_HYDROCARBON", gsea_results_df$pathway, value=T)
grep("OXIDATIVE_STRESS", gsea_results_df$pathway, value=T)
grep("INFLAMMATORY", gsea_results_df$pathway, value=T)


gsea_results_df$pathway_pretty <- sapply(gsea_results_df$pathway, function(p) {
  
  # drop DB name and replace underscores
  p <- paste0(unlist(strsplit(p, "_"))[-1], collapse=" ")
  
  # make title
  p <- str_to_title(p)
  
  # drop WP IDs
  p <- gsub("Wp[0-9]+", "", p)
  p <- trimws(p)
  
  return(p)
  
})

# add the logp
gsea_results_df$logp <- -log10(gsea_results_df$pval)

pathway_list <- c("GOBP_REGULATION_OF_MYELINATION",
                  "WP_ARYL_HYDROCARBON_RECEPTOR_PATHWAY_WP2586",
                  "WP_OXIDATIVE_STRESS_RESPONSE",
                  "HALLMARK_INFLAMMATORY_RESPONSE",
                  "KEGG_TRYPTOPHAN_METABOLISM",
                  "REACTOME_PEROXISOMAL_LIPID_METABOLISM")

gsea_subset <- gsea_results_df[gsea_results_df$pathway %in% pathway_list,]
gsea_subset <- gsea_subset[order(gsea_subset$pval),]
gsea_subset$pathway_pretty <- factor(gsea_subset$pathway_pretty,
                                     levels=rev(gsea_subset$pathway_pretty))

ggplot(gsea_subset,
       aes(x=logp, y=pathway_pretty,fill=NES)) +
  geom_point(pch=21, size=5) +
  theme_bw() +
  xlim(0, 7) +
  scale_fill_gradient2(
    low  = "blue",
    mid  = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(x="-log10(PValue)", y=NULL)
ggsave(paste0(out_dir, "pathway_dot_plot.png"), width=7, height=5)

# make heatmaps!

# pull out vsd
vsd <- vst(dds, blind=F)

counts <- assay(vsd)
# fix rownames
rownames(counts) <- gsub("\\.[0-9]+", "", rownames(counts))

for (pathway in pathway_list) {
  
  pathway_gsea <- gsea_results_df[gsea_results_df$pathway == pathway,]
  
  pathway_gene_list <- unlist(pathway_gsea$leadingEdge)
  
  subset_mat <- counts[pathway_gene_list,]
  
  # scale it
  subset_mat <- t(scale(t(subset_mat)))
  
  # replace gene IDs with names
  rownames(subset_mat) <- gene_xref[pathway_gene_list,]$gene_symbol
  
  # height adjustment
  height_val <- 5
  
  if (nrow(subset_mat) > 40) {
    height_val <- 8
  }
  
  pdf(paste0(out_dir, to_snake_case(pathway_gsea$pathway_pretty),".expr_heatmap.pdf"),
      width=6, height=height_val)
  print(Heatmap(subset_mat,
          col=colorRamp2(c(-1.5,0,1.5), c("blue","white","red")),
          name = "Z-Score",
          cluster_columns = F))
  dev.off()
  
  
}







