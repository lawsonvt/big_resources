library(CASSIA)
library(openxlsx)
library(gtools)
library(Seurat)
library(SeuratObject)
library(stringr)

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/GSE133283/"


out_dir <- paste0(root_dir, "results/cassia_endothelial_subclustering_results_analysis.pct_diff/")
dir.create(out_dir, showWarnings = F, recursive = T)


# read in CASSIA results
gemini_results <- read.csv(paste0(root_dir, "results/cassia_endothelial_subclustering.pct_diff/results_gemini_lc_summary.csv"))
claude_results <- read.csv(paste0(root_dir, "results/cassia_endothelial_subclustering.pct_diff/results_claude_sonnet_summary.csv"))
#olgpt_results <- read.csv(paste0(root_dir, "results/cassia_annotate/results_ollama_gpt_summary.csv"))

# create smaller versions
gemini_simple <- gemini_results[,c("Cluster.ID","Predicted.General.Cell.Type",
                                   "Predicted.Detailed.Cell.Type",
                                   "Possible.Mixed.Cell.Types")]
colnames(gemini_simple) <- c("harmony_clusters",
                             "gemini_prediction",
                             "gemini_detailed_prediction",
                             "gemini_possible_mix")

claude_simple <- claude_results[,c("Cluster.ID","Predicted.General.Cell.Type",
                                   "Predicted.Detailed.Cell.Type",
                                   "Possible.Mixed.Cell.Types")]

colnames(claude_simple) <- c("harmony_clusters",
                             "claude_prediction",
                             "claude_detailed_prediction",
                             "claude_possible_mix")

# olgpt_simple <- olgpt_results[,c("Cluster.ID","Predicted.General.Cell.Type",
#                                  "Predicted.Detailed.Cell.Type",
#                                  "Possible.Mixed.Cell.Types")]
# 
# 
# colnames(olgpt_simple) <- c("harmony_clusters",
#                             "ollama_prediction",
#                             "ollama_detailed_prediction",
#                             "ollama_possible_mix")




# merge em up for a comparison


cell_comp <- merge(claude_simple,
                   gemini_simple,
                   by="harmony_clusters")

# cell_comp <- merge(cell_comp,
#                    olgpt_simple,
#                    by="harmony_clusters",
#                    all.x=T)

cell_comp <- cell_comp[mixedorder(cell_comp$harmony_clusters),]


# output the comp

wb <- createWorkbook()
addWorksheet(wb, "comparison")
writeData(wb, "comparison", cell_comp)

n_rows <- nrow(cell_comp) + 1  # +1 for header row

setColWidths(wb, "comparison", cols=1, width=20)
setColWidths(wb, "comparison", cols=2:ncol(cell_comp), width=40)

setRowHeights(wb, "comparison", rows = 2:n_rows, heights = 80)

wrap_style <- createStyle(wrapText = TRUE)
addStyle(wb, "comparison",
         style = wrap_style,
         rows = 2:n_rows,
         cols = 2:ncol(cell_comp),       # whichever columns should wrap
         gridExpand = TRUE)    # applies style to every row/col combination

saveWorkbook(wb, paste0(out_dir, "cassia_results_comparison.xlsx"), overwrite = T)

saveRDS(cell_comp, paste0(out_dir, "cassia_results_comparison.RDS"))

# move report htmls into this folder

html_reports <- list.files(paste0(root_dir, "results/cassia_endothelial_subclustering.pct_diff/"),
                           "_report.html", full.names = T)
html_reports <- html_reports[!grepl("scored", html_reports)]

file.copy(from=html_reports, to=out_dir)

# make some Seurat plots

endo_seu <- LoadSeuratRds(paste0(root_dir,
                                 "results/endothelial_subclustering/endothelial_subcluster.seurat.RDS"))


metadata <- endo_seu@meta.data

# first rename clustees
levels(metadata$endothelial_clusters) <- paste0("cluster", str_pad(levels(metadata$endothelial_clusters), 
                                                                   width=2, side="left", pad="0"))

colnames(cell_comp)[1] <- "endothelial_clusters"

# merge em in
metadata$cell_id <- rownames(metadata)

metadata <- merge(metadata,
                  cell_comp, 
                  by="endothelial_clusters")

rownames(metadata) <- metadata$cell_id

# adjust string width for long predictions

metadata$claude_detail_pretty <- str_wrap(metadata$claude_detailed_prediction, width=50)
metadata$gemini_detail_pretty <- str_wrap(metadata$gemini_detailed_prediction, width=50)

metadata$claude_cluster <- paste0(metadata$endothelial_clusters, "\n",
                                  metadata$claude_prediction)
metadata$gemini_cluster <- paste0(metadata$endothelial_clusters, "\n",
                                  metadata$gemini_prediction)

# add metadata back to seurat object
metadata <- metadata[colnames(endo_seu),]

endo_seu@meta.data <- metadata

# original UMAP
DimPlot(endo_seu, reduction="umap.endothelial_pca", group.by= "endothelial_clusters",
        label=T, label.box = T, label.size = 3, stroke.size = 1) + theme(legend.position = "none") +
  labs(x="UMAP1", y="UMAP2", title="Predicted Clusters")
ggsave(paste0(out_dir, "endothelial_subcluster.umap.png"), width=14, height=10)


DimPlot(endo_seu, reduction="umap.endothelial_pca", group.by= "gemini_cluster",
        label=T, label.box = T, label.size=3, stroke.size = 1) + theme(legend.position = "none") +
  labs(x="UMAP1", y="UMAP2", title="Gemini Simple Prediction")
ggsave(paste0(out_dir, "endothelial_subcluster.umap.gemini_basic.png"), width=14, height=10)


DimPlot(endo_seu, reduction="umap.endothelial_pca", group.by= "claude_cluster",
        label=T, label.box = T, label.size=3, stroke.size = 1) + theme(legend.position = "none") +
  labs(x="UMAP1", y="UMAP2", title="Claude Simple Prediction")
ggsave(paste0(out_dir, "endothelial_subcluster.umap.claude_basic.png"), width=14, height=10)


DimPlot(endo_seu, reduction="umap.endothelial_pca", group.by= "gemini_detail_pretty",
        label=T, label.box = T, label.size=3, stroke.size = 1) + theme(legend.position = "none") +
  labs(x="UMAP1", y="UMAP2", title="Gemini Detailed Prediction")
ggsave(paste0(out_dir, "endothelial_subcluster.umap.gemini_detail.png"), width=14, height=10)

DimPlot(endo_seu, reduction="umap.endothelial_pca", group.by= "claude_detail_pretty",
        label=T, label.box = T, label.size=3, stroke.size = 1) + theme(legend.position = "none") +
  labs(x="UMAP1", y="UMAP2", title="Claude Detailed Prediction")
ggsave(paste0(out_dir, "endothelial_subcluster.umap.claude_detail.png"), width=14, height=10)






