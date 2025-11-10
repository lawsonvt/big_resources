library(Seurat) # single cell analysis
library(msigdbr)  # for gene sets / pathways
library(dplyr) # data manipulation
library(fgsea) # GSEA test
library(openxlsx) # excel output
library(snakecase) # better file names
library(ggplot2) # plotting
library(cowplot) # combining plots

# Initial setup ----------------------------------------------------------------

# create output directory for results
out_dir <- "~/Documents/projects/lukenslab/ben_wendell/gsea_results/"
dir.create(out_dir, showWarnings = F)

# read in Seurat object
seu <- LoadSeuratRds("~/Documents/projects/lukenslab/ben_wendell/seurat_clusters_labeled_res.1_10_23_2025.RDS")

# Load in msigdbr and get gene sets --------------------------------------------
hallmark_gene_sets <- msigdbr(species = "Mus musculus", collection = "H")
reactome_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
wikipath_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:WIKIPATHWAYS")

gobp_gene_sets <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")
gomf_gene_sets <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:MF")

# NOTE: More pathways can be added, check out msigdbr_collections() for more info

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


# Calculate GSEA results per cell type -----------------------------------------

# pull out all cell types
cell_types <- unique(seu$cell_type)
# pull out conditions
conditions <- unique(seu$condition)

# parameters for the analysis
deg_method <- "MAST" # can be changed to any method supported by Seurat

results <- lapply(cell_types, function(cell_type) {
  
  print(cell_type)
  # subset seurat to just the cell type of interest 
  # (speeding up the DEG analysis)
  seu_subset <- subset(seu, idents = cell_type)
  
  # reset idents to condition for differential expression analysis
  Idents(seu_subset) <- "condition"
  
  print("Determine DEGs ...")
  
  # determine DEGs for cell type
  degs <- FindMarkers(seu_subset,
                      ident.1="FAD_KO", # knock out conditon
                      ident.2="FAD_WT", # wildtype condition
                      test.use=deg_method, # DEG method
                      logfc.threshold=0.1, # log fold change threshold
                      min.pct=0.01, # minimum percentage of cells for gene expression
                      verbose=F)
  degs$gene_name <- rownames(degs)
  
  # create a stat value that rank genes on directionality and significance
  # find minimum non-zero p-value
  min_nz_pval <- min(degs$p_val[degs$p_val > 0])
  # create "stat" value for ranking
  degs$stat <- -log10(pmax(degs$p_val, min_nz_pval)) * sign(degs$avg_log2FC)
  
  # NOTE: The ranking can easily be replaced by using log2FC
  
  # Create ranked gene list (by log fold change or stat)
  # Important: remove NAs and sort
  ranked_genes <- degs %>%
    filter(!is.na(stat) &
             !is.infinite(stat)) %>%
    arrange(desc(stat)) %>%
    pull(stat, name = gene_name)
  
  print("Run GSEA on ranked DEGs ...")
  
  # for each set of gene sets, run GSEA
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
  
  # return GSEA results as well as ranked list
  return(list(ranked_genes=ranked_genes,
              gsea_results=total_gsea_results,
              degs=degs))
  
})
names(results) <- cell_types

# save results
saveRDS(results, file=paste0(out_dir, "gsea_results.RDS"))

# Output results to excel files ------------------------------------------------

# for each cell type, pull out results and write to excel
for (cell_type in names(results)) {
  
  cell_results <- results[[cell_type]]
  
  write.xlsx(cell_results$gsea_results,
             file = paste0(out_dir, 
                           to_snake_case(cell_type), 
                           ".", paste0(conditions, collapse="_minus_"),
                           ".gsea_results.xlsx"),
             colWidths="auto")
  
  
}


# GSEA table results plots -----------------------------------------------------
top_count <- 10 # top number of results to print out

# for each cell type, print out a table showing the top results
for (cell_type in names(results)) {
  
  cell_results <- results[[cell_type]]
  
  # make a plot dir for these results
  plot_dir <- paste0(out_dir, to_snake_case(cell_type), "_plots/")
  dir.create(plot_dir, showWarnings = F)
  
  # a separate plot for each set of gene sets
  for (gene_set_name in names(cell_results$gsea_results)) {
    
    # pull out the top ones
    top_gene_sets <- cell_results$gsea_results[[gene_set_name]]$pathway[1:top_count]
    
    # if no results, skip
    if (length(top_gene_sets) == 0) {
      next
    }
  
    # command to plot table
    plotGseaTable(
      total_gene_sets[[gene_set_name]][top_gene_sets], # list of top gene sets
      cell_results$ranked_genes, # the ranked genes
      cell_results$gsea_results[[gene_set_name]], # the GSEA results
      pathwayLabelStyle = list(size=8)
    )
    ggsave(paste0(plot_dir, 
                  to_snake_case(cell_type), 
                  ".", paste0(conditions, collapse="_minus_"),
                  ".", gene_set_name,
                  ".top_pathway_gsea_table.png"), width=8, height=6, bg = "white")
      
  }
  
}

# GSEA Enrichment plots that visualize the enrichment score calculation -------
top_count <- 9

# for each cell type, make enrichment plots for the top gene sets
for (cell_type in names(results)) {
  
  cell_results <- results[[cell_type]]
  
  # make a plot dir
  plot_dir <- paste0(out_dir, to_snake_case(cell_type), "_plots/")
  dir.create(plot_dir, showWarnings = F)
  
  for (gene_set_name in names(cell_results$gsea_results)) {
    
    top_gene_sets <- cell_results$gsea_results[[gene_set_name]]$pathway[1:top_count]
    
    # create plots for each gene set, put into list
    plot_list <- lapply(top_gene_sets, function(pathway) {
      
      plotEnrichment(
        total_gene_sets[[gene_set_name]][[pathway]],
        cell_results$ranked_genes
      ) + labs(title = pathway)
      
    })
    # plot the list as grid into a single plot
    plot_grid(plotlist=plot_list, ncol=3)
    ggsave(paste0(plot_dir, 
                  to_snake_case(cell_type), 
                  ".", paste0(conditions, collapse="_minus_"),
                  ".", gene_set_name,
                  ".top_pathway_gsea_enrichment.png"), width=14, height=10, bg = "white")
  }
}







