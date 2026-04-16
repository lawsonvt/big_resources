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

out_dir <- paste0(root_dir, "results/differential_expression_analysis.sex_split/")
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

# the cell types
cells <- c("Astrocyte",
           "Microglia")

# define contrasts
contrasts <- list("P14-P7"=c("age","P14","P7"),
                  "P40-P14"=c("age","P40","P14"),
                  "P40-P7"=c("age","P40","P7"))

Idents(seu_obj) <- "age"

# find the markers
total_results <- lapply(unique(metadata$sex), function(sex_type) {
  
  print(sex_type)
  
  sex_results <- lapply(cells, function(cell) {
    
    print(cell)
    
    # subset the object
    subset_seu <- subset(seu_obj, subset = sex == sex_type & cell_type_partial == cell)
    
    cell_results <- lapply(contrasts, function(contrast) {
      
      print(paste0(contrast[2], "-", contrast[3]))
      
      results <- FindMarkers(subset_seu,
                             ident.1 = contrast[2],
                             ident.2 = contrast[3],
                             test.use = "MAST")
      results <- tibble::rownames_to_column(results, "gene")
      
      results$sex <- sex_type
      results$cell_type <- cell
      results$contrast <- paste0(contrast[2], "-", contrast[3])
      
      return(results)
      
    })
    
    return(cell_results)
  })
  
  names(sex_results) <- cells
  
  return(sex_results)
})
names(total_results) <- unique(metadata$sex)

# output the results
saveRDS(total_results, file=paste0(out_dir, "total_results.RDS"))

# output to excel files
for (sex in names(total_results)) {
  
  sex_results <- total_results[[sex]]
  
  for (cell in names(sex_results)) {
    
    cell_results <- sex_results[[cell]]
    
    write.xlsx(cell_results, file=paste0(out_dir, cell, ".", sex, ".deg_results.xlsx"),
               colWidths="auto")
    
  }
  
}





