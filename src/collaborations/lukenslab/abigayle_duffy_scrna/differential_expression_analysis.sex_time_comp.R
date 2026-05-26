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

out_dir <- paste0(root_dir, "results/differential_expression_analysis.sex_time_comp/")
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

contrasts <- list("female-male"=c("sex","female","male"))

Idents(seu_obj) <- "sex"

# find the markers
total_results <- lapply(unique(metadata$age), function(age_number) {
  
  print(age_number)
  
  age_results <- lapply(cells, function(cell) {
    
    print(cell)
    
    # subset the object
    subset_seu <- subset(seu_obj, subset = age == age_number & cell_type_partial == cell)
    
    cell_results <- lapply(contrasts, function(contrast) {
      
      print(paste0(contrast[2], "-", contrast[3]))
      
      results <- FindMarkers(subset_seu,
                             ident.1 = contrast[2],
                             ident.2 = contrast[3],
                             test.use = "MAST")
      results <- tibble::rownames_to_column(results, "gene")
      
      results$age <- age_number
      results$cell_type <- cell
      results$contrast <- paste0(contrast[2], "-", contrast[3])
      
      return(results)
      
    })
    
    return(cell_results)
  })
  
  names(age_results) <- cells
  
  return(age_results)
})
names(total_results) <- unique(metadata$age)

# output the results
saveRDS(total_results, file=paste0(out_dir, "total_results.RDS"))

# output to excel files
for (age in names(total_results)) {
  
  age_results <- total_results[[age]]
  
  for (cell in names(age_results)) {
    
    cell_results <- age_results[[cell]]
    
    write.xlsx(cell_results, file=paste0(out_dir, cell, ".", age, ".deg_results.xlsx"),
               colWidths="auto")
    
  }
  
}

# Make Volcano plots -----------------------------------------------------------

volcano_dir <- paste0(out_dir, "volcano_plots.cell_types/")

dir.create(volcano_dir, showWarnings = F)


top_genes <- 25

volcano_plot_list <-  lapply(mixedsort(names(total_results)), function(age) {
  
  age_results <- total_results[[age]]
  
  lapply(names(age_results), function(cell) {
    
    subset <- age_results[[cell]][[1]]
    
    subset$log_p <- -log10(subset$p_val)
    
    # cap Inf
    if (any(is.infinite(subset$log_p))) {
      
      cap_val <- max(subset[!is.infinite(subset$log_p),]$log_p, na.rm=T)
      
      subset[is.infinite(subset$log_p),]$log_p <- cap_val
      
    }
    
    subset_sig <- subset[subset$p_val_adj < 0.05 &
                           abs(subset$avg_log2FC) > 0.5,]
    
    subset_top <- subset_sig[1:top_genes,]
    
    if (any(is.na(subset_top$gene))) {
      subset_top <- subset_top[!is.na(subset_top$gene),]
    }
    
    logp_thresh <- min(subset_sig$log_p)
    
    p1 <- ggplot(subset,
                 aes(x=avg_log2FC,
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
      labs(x="Log2 Fold Change", y="-log10(P-Value)", title=age,
           subtitle=cell)
    p1
    ggsave(paste0(volcano_dir, to_snake_case(age), ".",
                  to_snake_case(cell), ".volcano_plot.png"), width=5, height=6)
    
    return(p1)
  })
})

plot_grid(plotlist = c(volcano_plot_list[[1]],
                       volcano_plot_list[[2]],
                       volcano_plot_list[[3]]), nrow = 3)
ggsave(paste0(out_dir, "cell_types.volcanoes.png"), width=7, height=10, bg="white")




