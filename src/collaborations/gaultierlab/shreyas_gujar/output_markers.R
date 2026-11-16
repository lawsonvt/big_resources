library(openxlsx)

root_dir <- "/Users/mjl3p/Documents/projects/gaultierlab/shreyas_gujar/snRNASeq/"
#root_dir <- "/home/mjl3p/projects/gaultierlab/shreyas_gujar/snRNASeq/"

out_dir <- paste0(root_dir, "results/output_markers/")
dir.create(out_dir, showWarnings = F, recursive = T)

# load in markers
all_markers <- readRDS(paste0(root_dir, "results/find_all_markers/all_markers.RDS"))

column_order <- c("cluster","gene","avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")

marker_list <- lapply(unique(all_markers$cluster), function(cluster) {
  
  subset <- all_markers[all_markers$cluster == cluster,]
  
  subset <- subset[,column_order]
  
  subset$delta_pct <- subset$pct.1 - subset$pct.2
  
  return(subset)
  
})
names(marker_list) <- unique(all_markers$cluster)


top_marker_list <- lapply(marker_list, function(subset) {
  
  subset <- subset[subset$p_val_adj < 0.05 &
                  subset$delta_pct > 0.1,]

  subset <- subset[order(subset$avg_log2FC, 
                         decreasing = T),]
  
  return(subset[1:20,])
  
})

# output em
write.xlsx(marker_list, file=paste0(out_dir, "total_cluster_markers.xlsx"), colWidths="auto")
write.xlsx(top_marker_list, file=paste0(out_dir, "top_cluster_markers.xlsx"), colWidths="auto")



