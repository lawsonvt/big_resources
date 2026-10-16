library(stringr)
library(VennDiagram)

root_dir <- "~/Documents/projects/harrislab/hto_test/"

# test comment

out_dir <- paste0(root_dir, "results/hto_test/")
dir.create(out_dir, showWarnings = F, recursive = T)

whitelist_10x <- read.delim(paste0("~/Documents/projects/harrislab/x1247-GEX/fastq_files/3M-3pgex-may-2023.txt"),
                            header=F)[,1]

top_gex <- read.delim(paste0(root_dir, "top_GEX.txt"),
                      header=F)[,1]
top_fb <- read.delim(paste0(root_dir, "top_FB.txt"),
                     header=F)[,1]

sum(top_gex %in% whitelist_10x)
sum(top_fb %in% whitelist_10x)


sum(top_gex %in% top_fb) # yup, BIG overlap


barcodes <- list(
  
  FB_subset=top_fb,
  GEX_subset=top_gex
)

venn_object <- venn.diagram(
  x = barcodes,
  category.names = names(barcodes),
  filename = NULL,              # Displays directly in RStudio instead of saving a file
  
  # Proportional scaling arguments
  euler.d = TRUE,               # Enables Euler diagrams (drops empty intersections)
  scaled = TRUE,                # Scales circle sizes proportional to element counts
  
  # Styling and aesthetics
  fill = c("#E69F00", "#56B4E9"), # Soft grey, orange, and blue colors
  alpha = c(0.5, 0.5),                  # Semi-transparent so intersections show up
  lty = "blank",                             # Removes harsh black borders around circles
  cat.pos = c(-20, 180),                 # Position of labels relative to circles
  cat.dist = c(0.05, 0.05)             # Distance of labels from circles
)

pdf(paste0(out_dir, "gex_fb.venn.pdf"), width=7, height=5)
grid.newpage()
grid.draw(venn_object)
dev.off()

full_fb_counts <- read.delim(paste0(root_dir, "fb_barcode_counts.txt"),
                             header=F, sep=" ")
colnames(full_fb_counts) <- c("raw_count","cell_id")

# filter
full_fb_counts <- full_fb_counts[grepl("[0-9]+", full_fb_counts$raw_count),]
full_fb_counts$raw_count <- as.numeric(full_fb_counts$raw_count)

full_fb_counts$percent <- full_fb_counts$raw_count / sum(full_fb_counts$raw_count)
full_fb_counts$cell_id <- factor(full_fb_counts$cell_id,
                                 levels=full_fb_counts$cell_id)

ggplot(full_fb_counts[1:10,],
       aes(x=cell_id,
           y=percent * 100)) +
  geom_bar(stat="identity", fill="maroon", color="black") +
  facet_wrap(~ "Top FB Reads", ncol=1, scales="free_x") +
  geom_text(aes(label=round(percent * 100, digits=3)),
            vjust=-0.5, size=3) +
  theme_bw() +
  ylim(0, 21) +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  labs(x="Top Cell Barcode based on Total Read Percentage", 
       y="Percentage of Reads with Cell Barcode")
ggsave(paste0(out_dir, "top_fb_cell_barcode_percentages.png"), width=4, height=5)


