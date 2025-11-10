library(Seurat) 
library(msigdbr)  # For gene sets
library(dplyr)
library(fgsea)
library(openxlsx)
library(snakecase)
library(ggplot2)
library(cowplot)



out_dir <- "~/Documents/projects/lukenslab/ben_wendell/gsea_results/"

dir.create(out_dir, showWarnings = F)

# read in Seurat object

seu <- LoadSeuratRds("~/Documents/projects/lukenslab/ben_wendell/seurat_clusters_labeled_res.1_10_23_2025.RDS")



DimPlot(seu, reduction="umap", group.by=c("condition","cell_type"), alpha = 0.3)

meta <- seu@meta.data
table(meta$condition, meta$sample_id)
table(meta$cell_type)

# load in msigdbr
hallmark_gene_sets <- msigdbr(species = "Mus musculus", collection = "H")
reactome_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
wikipath_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:WIKIPATHWAYS")

gobp_gene_sets <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")
gomf_gene_sets <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:MF")

# function to convert data frames to lists
list_convert <- function(gene_sets) {
  gene_sets %>%
    split(x = .$gene_symbol, f = .$gs_name)
}

# put into one big list
total_gene_sets <- list(hallmark=list_convert(hallmark_gene_sets),
                        reactome=list_convert(reactome_gene_sets),
                        wikipathways=list_convert(wikipath_gene_sets),
                        gobp=list_convert(gobp_gene_sets),
                        gomf=list_convert(gomf_gene_sets))

# pull out all cell types
cell_types <- unique(seu$cell_type)
# pull out conditions
conditions <- unique(seu$condition)

deg_method <- "MAST"


results <- lapply(cell_types, function(cell_type) {
  
  print(cell_type)
  
  seu_subset <- subset(seu, idents = cell_type)
  
  Idents(seu_subset) <- "condition"
  
  print("Determine DEGs ...")
  
  # determine DEGs for cell type
  degs <- FindMarkers(seu_subset,
                      ident.1="FAD_KO",
                      ident.2="FAD_WT",
                      test.use=deg_method,
                      logfc.threshold=0.1,
                      min.pct=0.01,
                      verbose=F)
  degs$gene_name <- rownames(degs)
  
  # find minimum non-zero p-value
  min_nz_pval <- min(degs$p_val[degs$p_val > 0])
  
  # create "stat" value for ranking
  degs$stat <- -log10(pmax(degs$p_val, min_nz_pval)) * sign(degs$avg_log2FC)
  
  # Create ranked gene list (by log fold change or stat)
  # Important: remove NAs and sort
  ranked_genes <- degs %>%
    filter(!is.na(stat) &
             !is.infinite(stat)) %>%
    arrange(desc(stat)) %>%
    pull(stat, name = gene_name)
  
  print("Run GSEA on ranked DEGs ...")
  
  total_gsea_results <- lapply(names(total_gene_sets), function(gs_name) {
    
    print(gs_name)
    
    gene_sets <- total_gene_sets[[gs_name]]
    
    fgsea_results <- fgsea(
      pathways = gene_sets,
      stats = ranked_genes,
      minSize = 15, # Minimal size of a gene set to test
      maxSize = 500 # Maximal size of a gene set to test
    )
    
    return(fgsea_results[order(fgsea_results$pval),])
  })
  names(total_gsea_results) <- names(total_gene_sets)
  
  return(list(ranked_genes=ranked_genes,
              gsea_results=total_gsea_results))
  
})
names(results) <- cell_types


# save results
saveRDS(results, file=paste0(out_dir, "gsea_results.RDS"))

# output results to excel files

for (cell_type in names(results)) {
  
  cell_results <- results[[cell_type]]
  
  write.xlsx(cell_results$gsea_results,
             file = paste0(out_dir, 
                           to_snake_case(cell_type), 
                           ".", paste0(conditions, collapse="_minus_"),
                           ".gsea_results.xlsx"),
             colWidths="auto")
  
  
}

# output some plots

top_count <- 10

for (cell_type in names(results)) {
  
  cell_results <- results[[cell_type]]
  
  # make a plot dir
  plot_dir <- paste0(out_dir, to_snake_case(cell_type), "_plots/")
  dir.create(plot_dir, showWarnings = F)
  
  for (gene_set_name in names(cell_results$gsea_results)) {
    
    top_gene_sets <- cell_results$gsea_results[[gene_set_name]]$pathway[1:top_count]
    
    if (length(top_gene_sets) == 0) {
      next
    }
  
    plotGseaTable(
      total_gene_sets[[gene_set_name]][top_gene_sets],
      cell_results$ranked_genes,
      cell_results$gsea_results[[gene_set_name]],
      gseaParam = 0.5,
      pathwayLabelStyle = list(size=8)
    )
    ggsave(paste0(plot_dir, 
                  to_snake_case(cell_type), 
                  ".", paste0(conditions, collapse="_minus_"),
                  ".", gene_set_name,
                  ".top_pathway_gsea_table.png"), width=8, height=6, bg = "white")
      
  }
  
}

# GSEA plots
top_count <- 9

for (cell_type in names(results)) {
  
  cell_results <- results[[cell_type]]
  
  # make a plot dir
  plot_dir <- paste0(out_dir, to_snake_case(cell_type), "_plots/")
  dir.create(plot_dir, showWarnings = F)
  
  for (gene_set_name in names(cell_results$gsea_results)) {
    
    top_gene_sets <- cell_results$gsea_results[[gene_set_name]]$pathway[1:top_count]
    
    plot_list <- lapply(top_gene_sets, function(pathway) {
      
      plotEnrichment(
        total_gene_sets[[gene_set_name]][[pathway]],
        cell_results$ranked_genes
      ) + labs(title = pathway)
      
    })
    plot_grid(plotlist=plot_list, ncol=3)
    ggsave(paste0(plot_dir, 
                  to_snake_case(cell_type), 
                  ".", paste0(conditions, collapse="_minus_"),
                  ".", gene_set_name,
                  ".top_pathway_gsea_enrichment.png"), width=14, height=10, bg = "white")
  }
}







