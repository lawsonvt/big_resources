library(CASSIA)
library(openxlsx)
library(stringr)

# read in cell annotations
root_dir <- "/Users/mjl3p/Documents/projects/lukenslab/ashley_bolte/human_meninges/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/cassia_results_analysis.immune/")
dir.create(out_dir, showWarnings = F, recursive = T)


# read in CASSIA results
gemini_results <- read.csv(paste0(root_dir, "results/cassia_annotate_immune/results_gemini_lc_summary.csv"))
claude_results <- read.csv(paste0(root_dir, "results/cassia_annotate_immune/results_claude_sonnet_summary.csv"))

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

# merge em up

cell_comp <- merge(claude_simple,
                   gemini_simple,
                   by="harmony_clusters")

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

saveWorkbook(wb, paste0(out_dir, "cassia_results_comparison.immune.xlsx"), overwrite = T)

# move report htmls into this folder

html_reports <- list.files(paste0(root_dir, "results/cassia_annotate_immune/"),
                           "_report.html", full.names = T)
html_reports <- html_reports[!grepl("scored", html_reports)]

file.copy(from=html_reports, to=out_dir)



