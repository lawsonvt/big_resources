library(ggplot2)
library(plyr)
library(dplyr)

# load in data
logp_gsea <- readRDS("~/Documents/projects/lukenslab/ben_wendell/gsea_results.signed_logp/gsea_results.RDS")
logfc_gsea <- readRDS("~/Documents/projects/lukenslab/ben_wendell/gsea_results.log2fc/gsea_results.RDS")

cell_types <- names(logfc_gsea)
geneset_types <- names(logfc_gsea[[1]]$gsea_results)

merged_gsea <- lapply(cell_types, function(cell_type) {
  
  lapply(geneset_types, function(geneset_type) {
    
    logp_results <- logp_gsea[[cell_type]]$gsea_results[[geneset_type]]
    logfc_results <- logfc_gsea[[cell_type]]$gsea_results[[geneset_type]]
    
    merged_results <- merge(logp_results,
                            logfc_results,
                            by="pathway",
                            suffixes=c(".logp", ".logfc"))
    merged_results$cell_type <- cell_type
    merged_results$geneset_type <- geneset_type
    
    return(merged_results)
  })
  
})
merged_gsea_df <- bind_rows(merged_gsea)

table(merged_gsea_df$geneset_type,
      merged_gsea_df$cell_type)

cor(merged_gsea_df$NES.logp,
    merged_gsea_df$NES.logfc,
    use="complete.obs")

# 





