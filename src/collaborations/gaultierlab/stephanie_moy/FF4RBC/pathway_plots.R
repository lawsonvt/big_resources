library(openxlsx)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(stringr)
library(DESeq2)
library(snakecase)

root_dir <- "~/Documents/projects/gaultierlab/stephanie_moy/"

out_dir <- paste0(root_dir, "FF4RBC/results/pathway_plots.no_sva/")
dir.create(out_dir, showWarnings = F)

# pull in data
deg_results <- readRDS(paste0(root_dir, "FF4RBC/results/deseq2_workflow.no_sva/deg_results.RDS"))
dds <- readRDS(paste0(root_dir, "FF4RBC/results/deseq2_workflow.no_sva/dds.RDS"))
gsea_results <- readRDS(paste0(root_dir, "FF4RBC/results/gsea_workflow.no_sva/total.gsea_results.RDS"))

pathways <- readRDS(paste0(root_dir, "FF4RBC/results/gsea_workflow.no_sva/msigdb_genesets.RDS"))

pathways <- unlist(pathways, recursive = F, use.names=T)
names(pathways) <- gsub("[a-z]+\\.", "", names(pathways))

gene_xref <- deg_results[,c("gene_id","gene_name")]
gene_xref$gene_id <- gsub("\\.[0-9]+", "", gene_xref$gene_id)
rownames(gene_xref) <- gene_xref$gene_id


# pull out top pathways

gsea_results_df <- as.data.frame(bind_rows(gsea_results$`treatment-control`$gsea_results))

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

gsea_results_top <- gsea_results_df[gsea_results_df$padj < 0.05,]

# get a list of desired pathways from Stephanie

pathway_list <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
                  "REACTOME_INTERFERON_SIGNALING",
                  "REACTOME_INTERFERON_GAMMA_SIGNALING",
                  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                  "KEGG_OXIDATIVE_PHOSPHORYLATION",
                  "HALLMARK_FATTY_ACID_METABOLISM",
                  "KEGG_FATTY_ACID_METABOLISM",
                  "REACTOME_FATTY_ACID_METABOLISM",
                  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                  "KEGG_TGF_BETA_SIGNALING_PATHWAY",
                  "REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS",
                  "REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX",
                  "HALLMARK_APOPTOSIS",
                  "REACTOME_REGULATION_OF_APOPTOSIS")

# gsea_results_top[grepl("interferon", tolower(gsea_results_top$pathway)),]$pathway
# gsea_results_top[grepl("oxidative", tolower(gsea_results_top$pathway)),]$pathway
# gsea_results_top[grepl("fatty_acid", tolower(gsea_results_top$pathway)),]$pathway
# gsea_results_top[grepl("cholesterol", tolower(gsea_results_top$pathway)),]$pathway
# gsea_results_top[grepl("tgf", tolower(gsea_results_top$pathway)),]$pathway
# gsea_results_top[grepl("apoptosis", tolower(gsea_results_top$pathway)),]$pathway

# make GSEA plots for each

enrich_dir <- paste0(out_dir, "enrichment_plots/")
dir.create(enrich_dir, showWarnings = F)

for (pathway in pathway_list) {
  
  pathway_pretty <- gsea_results_df[gsea_results_df$pathway == pathway,]$pathway_pretty
  
  plotEnrichment(
    pathways[[pathway]],
    gsea_results$`treatment-control`$ranked_genes
  ) + labs(title=pathway_pretty,
           x="Rank", y="Enrichment Score")
  ggsave(paste0(enrich_dir, to_snake_case(pathway), ".enrichment.png"), width=7, height=4)
  
}

# heatmaps!

vsd <- vst(dds, blind=F)

counts <- assay(vsd)
meta <- colData(dds)

column_names <- paste0(meta$condition,"|",
                       meta$sample_id)

colnames(counts) <- column_names

heatmap_dir <- paste0(out_dir, "heatmap_plots/")
dir.create(heatmap_dir, showWarnings = F)

for (pathway in pathway_list) {
  
  pathway_pretty <- gsea_results_df[gsea_results_df$pathway == pathway,]$pathway_pretty
  
  edge_genes <- gsea_results_df[gsea_results_df$pathway == pathway,]$leadingEdge
  edge_genes <- unlist(strsplit(edge_genes[[1]], ","))
  
  gene_ids <- gene_xref[gene_xref$gene_name %in% edge_genes,]$gene_id
  
  subset_mat <- counts[gene_ids,]
  
  # scale it
  subset_mat <- t(scale(t(subset_mat)))
  
  # replace gene IDs with names
  rownames(subset_mat) <- gene_xref[gene_ids,]$gene_name
  
  # height adjustment
  height_val <- 7
  
  if (nrow(subset_mat) > 40) {
    height_val <- 12
  }
  
  pdf(paste0(heatmap_dir, to_snake_case(pathway),".expr_heatmap.pdf"),
      width=6, height=height_val)
  print(Heatmap(subset_mat,
                col=colorRamp2(c(-2,0,2), c("blue","white","red")),
                name = "Z-Score",
                cluster_columns = F, column_title = pathway_pretty))
  dev.off()
  
  
}



