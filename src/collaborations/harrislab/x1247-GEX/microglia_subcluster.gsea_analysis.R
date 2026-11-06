library(Seurat) # single cell analysis
library(msigdbr)  # for gene sets / pathways
library(dplyr) # data manipulation
library(fgsea) # GSEA test
library(openxlsx) # excel output
library(snakecase) # better file names
library(ggplot2) # plotting
library(cowplot) # combining plots

# Initial setup ----------------------------------------------------------------

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/x1247-GEX/"

# create output directory for results
out_dir <- paste0(root_dir, "results/microglia_subcluster.gsea_analysis/")
dir.create(out_dir, showWarnings = F)

# Load in msigdbr and get gene sets --------------------------------------------
hallmark_gene_sets <- msigdbr(species = "Mus musculus", collection = "H")
reactome_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
wikipath_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:WIKIPATHWAYS")
kegg_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY")

gobp_gene_sets <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")
gomf_gene_sets <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:MF")

# function to convert data frames to lists
list_convert <- function(gene_sets) {
  gene_sets %>%
    split(x = .$gene_symbol, f = .$gs_name)
}

# put into one big list
total_gene_sets <- list(hallmark=list_convert(hallmark_gene_sets),
                        kegg=list_convert(kegg_gene_sets),
                        reactome=list_convert(reactome_gene_sets),
                        wikipathways=list_convert(wikipath_gene_sets),
                        gobp=list_convert(gobp_gene_sets),
                        gomf=list_convert(gomf_gene_sets))

# Read in data -----------------------------------------------------------------

results_file <- paste0(root_dir, "results/microglia_subcluster.diff_exp_analysis/microglia_cells.de_results.RDS")

results <- readRDS(results_file)

for (cell in names(results)) {
  
  degs <- results[[cell]]
  
  # Create ranked gene list (by log fold change or stat)
  # Important: remove NAs and sort
  ranked_genes <- degs %>%
    filter(!is.na(stat) &
             !is.infinite(stat)) %>%
    arrange(desc(stat)) %>%
    pull(stat, name = gene)
  
  total_gsea_results <- lapply(names(total_gene_sets), function(gs_name) {
    
    print(gs_name)
    
    gene_sets <- total_gene_sets[[gs_name]]
    # run GSEA
    fgsea_results <- fgsea(
      pathways = gene_sets,
      stats = ranked_genes,
      minSize = 15, # Minimal size of a gene set to test
      maxSize = 500 # Maximal size of a gene set to test
    )
    
    return(fgsea_results[order(fgsea_results$pval),])
  })
  names(total_gsea_results) <- names(total_gene_sets)
  
  output_stump <- paste0(out_dir, cell, ".gsea_results")
  
  # save GSEA results as well as ranked list
  saveRDS(list(ranked_genes=ranked_genes,
               gsea_results=total_gsea_results,
               degs=degs), file=paste0(output_stump, ".RDS"))
  
  # save results to excel file
  write.xlsx(total_gsea_results,
             file=paste0(output_stump, ".xlsx"),
             colWidths="auto")
  
  
}


