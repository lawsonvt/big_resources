library(Seurat)
library(SeuratObject)
library(gtools)
library(stringr)
library(ggplot2)
library(reshape2)
library(openxlsx)

root_dir <- "~/Documents/projects/lukenslab/abigayle_duffy/"

out_dir <- paste0(root_dir, "results/cell_breakdown_plots/")
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

# add cell category
metadata$cell_category <- metadata$cell_type_partial
metadata[!metadata$cell_type_partial %in% c("Microglia","Astrocyte"),]$cell_category <- "Other"

ggplot(metadata,
       aes(x=orig.ident,
           fill=cell_category)) +
  geom_bar() +
  theme_bw() +
  scale_fill_manual(values=c("blue","orange","grey")) +
  labs(x=NULL, y="Cell Count", fill=NULL) 
ggsave(paste0(out_dir, "cell_category_per_library.png"), width=6, height=4)
  
cell_counts <- as.data.frame(table(metadata$orig.ident,
                                    metadata$cell_category))
cell_counts <- dcast(cell_counts, Var1 ~ Var2, value.var = "Freq")
colnames(cell_counts)[1] <- "Library"

write.xlsx(cell_counts, file=paste0(out_dir, "cell_category_per_library.xlsx"),
           colWidths="auto")




