library(ggplot2)
library(stringr)
library(snakecase)
library(openxlsx)

root_dir <- "~/Documents/projects/lukenslab/abigayle_duffy/"

out_dir <- paste0(root_dir, "results/gsea_plots/")
dir.create(out_dir, showWarnings = F)

# function to clean up pathway names
pathway_pretty <- function(p) {
  
  # drop DB name and replace underscores
  p <- paste0(unlist(strsplit(p, "_"))[-1], collapse=" ")
  
  # make title
  p <- str_to_title(p)
  
  # drop WP IDs
  p <- gsub("Wp[0-9]+", "", p)
  p <- trimws(p)
  
  # wrap text
  p <- str_wrap(p, width=40)
  
  return(p)
  
}

# read in excel files
gsea_analysis_dir <- paste0(root_dir, "results/gsea_analysis/")


gsea_result_files <- list.files(gsea_analysis_dir,
                                ".xlsx", full.names = T)

db <- "gobp"
top <- 10

for (file in gsea_result_files) {
  
  data <- read.xlsx(file, sheet = db)
  
  #print(summary(data$padj))
  
  data$pathway_pretty <- sapply(data$pathway, pathway_pretty)
  
  data$logp <- -log10(data$pval)
  # ensure proper order
  data <- data[order(data$pval),]
  data$pathway_pretty <- factor(as.character(data$pathway_pretty),
                                levels=rev(as.character(data$pathway_pretty)))
  
  ggplot(data[1:top,],
         aes(x=logp, y=pathway_pretty,fill=NES)) +
    geom_point(pch=21, size=5) +
    theme_bw() +
    xlim(0, max(data$logp)) +
    scale_fill_gradient2(
      low  = "blue",
      mid  = "white",
      high = "red",
      midpoint = 0
    ) +
    labs(x="-log10(PValue)", y=NULL)  
  
  outfile <- paste0(out_dir, gsub(".xlsx", paste0(".", db, ".png"), basename(file)))
  
  ggsave(outfile, width=7, height=6)
}

