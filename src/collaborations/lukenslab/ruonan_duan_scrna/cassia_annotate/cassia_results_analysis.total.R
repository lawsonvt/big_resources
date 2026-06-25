library(CASSIA)
library(openxlsx)

# read in cell annotations
root_dir <- "/Users/mjl3p/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/cassia_annotate/cassia_results_analysis.total/")
dir.create(out_dir, showWarnings = F, recursive = T)

# originally determined cell annotations
cell_annot <- read.xlsx(paste0(root_dir, "cell type 20260329.xlsx"))
colnames(cell_annot) <- c("harmony_clusters","annotated_cell_type")

# read in CASSIA results
gemini_results <- read.csv(paste0(root_dir, "results/cassia_annotate/cassia_annotate_total/results_gemini_lc_summary.csv"))
claude_results <- read.csv(paste0(root_dir, "results/cassia_annotate/cassia_annotate_total/results_claude_sonnet_summary.csv"))
olgpt_results <- read.csv(paste0(root_dir, "results/cassia_annotate/cassia_annotate_total/results_ollama_gpt_summary.csv"))

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

olgpt_simple <- olgpt_results[,c("Cluster.ID","Predicted.General.Cell.Type",
                          "Predicted.Detailed.Cell.Type",
                          "Possible.Mixed.Cell.Types")]


colnames(olgpt_simple) <- c("harmony_clusters",
                             "ollama_prediction",
                             "ollama_detailed_prediction",
                             "ollama_possible_mix")




# merge em up for a comparison

cell_comp <- merge(cell_annot,
                   claude_simple,
                   by="harmony_clusters")

cell_comp <- merge(cell_comp,
                   gemini_simple,
                   by="harmony_clusters")

cell_comp <- merge(cell_comp,
                   olgpt_simple,
                   by="harmony_clusters",
                   all.x=T)


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

saveWorkbook(wb, paste0(out_dir, "cassia_results_comparison.total.xlsx"), overwrite = T)

# move report htmls into this folder

html_reports <- list.files(paste0(root_dir, "results/cassia_annotate/cassia_annotate_total/"),
                           "_report.html", full.names = T)
html_reports <- html_reports[!grepl("scored", html_reports)]

file.copy(from=html_reports, to=out_dir)






# output the results

high_out <- merge(cell_annot,
                  high_results,
                  by.x="harmony_clusters",
                  by.y = "Cluster.ID")

# output to excel file
wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", high_out)

n_rows <- nrow(high_out) + 1  # +1 for header row

setColWidths(wb, "Sheet1", cols=c(1,6,8,9,10,11,12), width = 20)
setColWidths(wb, "Sheet1", cols = c(2,3,4,5,7), widths = 40)

setRowHeights(wb, "Sheet1", rows = 2:n_rows, heights = 80)

wrap_style <- createStyle(wrapText = TRUE)
addStyle(wb, "Sheet1",
         style = wrap_style,
         rows = 1:n_rows,
         cols = 1:ncol(high_out),       # whichever columns should wrap
         gridExpand = TRUE)    # applies style to every row/col combination

saveWorkbook(wb, paste0(out_dir, "cassia_results_claude_sonnet.total.xlsx"), overwrite = T)







