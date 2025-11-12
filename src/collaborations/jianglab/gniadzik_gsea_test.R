library(msigdbr)  # for gene sets / pathways
library(dplyr) # data manipulation
library(fgsea) # GSEA test
library(openxlsx) # excel output
library(ggplot2) # plotting
library(cowplot) # combining plots
library(data.table)

out_dir <- "~/Documents/projects/jianglab/weronika_gniadzik/gsea_results/"
dir.create(out_dir, showWarnings = F)


# Read in the files ------------------------------------------------------------
in_dir <- "~/Documents/projects/jianglab/weronika_gniadzik/"

csv_files <- list.files(in_dir, ".csv", full.names=T)

data_list <- lapply(csv_files, function(file) {
  
  data <- read.csv(file, sep = ";", dec = ",")
  colnames(data) <- c("gene", "log2FC", "neg_logFDR")
  
  data$signed_neg_logFDR <- data$neg_logFDR *
    sign(data$log2FC)
  
  # remove duplicates
  data <- data[order(data$neg_logFDR, decreasing=T),]
  data <- data[!duplicated(data$gene),]
  
  return(data)
})
names(data_list) <- gsub(".csv", "", basename(csv_files))

# Get pathways -----------------------------------------------------------------

reactome_gene_sets <- msigdbr(species = "human", collection = "C2", subcollection = "CP:REACTOME")
gobp_gene_sets <- msigdbr(species = "human", collection = "C5", subcollection = "GO:BP")

# function to convert data frames to lists
list_convert <- function(gene_sets) {
  gene_sets %>%
    split(x = .$gene_symbol, f = .$gs_name)
}

# put into one big list
total_gene_sets <- list(reactome=list_convert(reactome_gene_sets),
                        gobp=list_convert(gobp_gene_sets))

# Perform GSEA Test ------------------------------------------------------------
gsea_results <- lapply(data_list, function(data) {
  
  ranked_genes <- data %>%
    filter(!is.na(signed_neg_logFDR)) %>%
    arrange(desc(signed_neg_logFDR)) %>%
    pull(signed_neg_logFDR, name = gene)
  
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
    
    fgsea_results$avg_logfc <- lapply(fgsea_results$leadingEdge, function(gene_list) {
      
      mean(data[data$gene %in% gene_list,]$log2FC)
      
    })
    
    # reorder columns
    fgsea_results <- setcolorder(fgsea_results,
                                 c(colnames(fgsea_results)[colnames(fgsea_results) != "leadingEdge"],
                                     "leadingEdge"))
    
    return(fgsea_results[order(fgsea_results$pval),])
  })
  names(total_gsea_results) <- names(total_gene_sets)
  
  # return GSEA results as well as ranked list
  return(list(ranked_genes=ranked_genes,
              gsea_results=total_gsea_results))
})

# Write results to excel files
for (data_name in names(gsea_results)) {
  
  gsea_data <- gsea_results[[data_name]]
  
  write.xlsx(gsea_data$gsea_results,
             file=paste0(out_dir, data_name, ".gsea_results.xlsx"))
  
}







