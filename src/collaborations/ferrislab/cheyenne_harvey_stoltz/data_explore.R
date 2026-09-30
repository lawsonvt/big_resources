library(openxlsx)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(snakecase)
library(mixOmics)
library(pheatmap)
library(RColorBrewer)
library(gtools)
library(pmp)

root_dir <- "~/Documents/projects/ferrislab/cheyenne_harvey_stolz/"

out_dir <- paste0(root_dir, "results/data_explore/")
dir.create(out_dir, showWarnings = F, recursive = T)

# read in data
metadata <- read.xlsx(paste0(root_dir, "BAF_CH_4681_Random number.xlsx"))
colnames(metadata)[1] <- "baf"
rownames(metadata) <- metadata$baf

# addded run data
run_data <- read.xlsx(paste0(root_dir, "BAF_CH_4681_Random number.xlsx"),
                      sheet = "Run Batch")
colnames(run_data)[1] <- "baf"


# add run data in
metadata <- merge(metadata,
                  run_data[,c("baf","LC-MS_pos","LC-MS_neg","run_batch")],
                  by="baf")
rownames(metadata) <- metadata$baf

metadata <- metadata[mixedorder(metadata$baf),]


n_matrix <- read.xlsx(paste0(root_dir, "BAF_CH_4681_Data.xlsx"), 
                      sheet = "p-value", startRow = 5, cols=1:51)

# add in sample type for ease

metadata$sample_type <- apply(metadata, 1, function(r) {
  
  if(grepl("Control 1", r[[4]])) {
    paste0(r[[3]], "c1")
  } else if (grepl("Control 2", r[[4]])) {
    paste0(r[[3]], "c2")
  } else if (grepl("Experimental 1", r[[4]])) {
    paste0(r[[3]], "e1")
  } else {
    paste0(r[[3]], "e2")
  }
  
})

# reformat n matrix

metabolite_xref <- n_matrix[,1:2]

metabolite_xref$metabolite_id <- 1:nrow(metabolite_xref)

dupes <- metabolite_xref[duplicated(metabolite_xref$Metabolite.name),]$Metabolite.name

# there are duplicates, so I can't simplify the matrix by using rownames
# WHY ARE THERE DUPLICATE METABOLITE NAMES?!?!


rownames(n_matrix) <- metabolite_xref$metabolite_id
rownames(metabolite_xref) <- metabolite_xref$metabolite_id

n_matrix$Metabolite.name <- NULL
n_matrix$Ontology <- NULL

# fix column names to match metadata baf
colnames(n_matrix) <- sapply(colnames(n_matrix), function(x) {
  as.numeric(unlist(strsplit(x, "_"))[1])
})

metadata[!metadata$baf %in% colnames(n_matrix),]

# remove from metadata missing sample
metadata <- metadata[metadata$baf %in% colnames(n_matrix),]

# ensure same order
n_matrix <- n_matrix[,rownames(metadata)]

# PCA analysis

pseudocount <- min(n_matrix[n_matrix > 0], na.rm = TRUE) / 2

n_matrix_log <- log2(n_matrix + pseudocount)

n_pca <- prcomp(t(n_matrix_log), scale. = T)

pca_vals <- data.frame(baf=rownames(n_pca$x),
                       PC1=n_pca$x[,1],
                       PC2=n_pca$x[,2])

pca_vals <- merge(pca_vals,
                  metadata,
                  by="baf") 

rownames(pca_vals) <- pca_vals$baf

ggplot(pca_vals,
       aes(x=PC1, y=PC2, fill=sample_type)) +
  geom_point(shape=21, size=3) +
  theme_bw() +
  geom_text_repel(aes(label=sample_type),
                  alpha=0.5) +
  scale_fill_brewer(palette="Set3")
ggsave(paste0(out_dir, "pca_plot.sample_type.png"), width=7, height=5)

ggplot(pca_vals,
       aes(x=PC1, y=PC2, fill=Group)) +
  geom_point(shape=21, size=3) +
  theme_bw() +
  geom_text_repel(aes(label=sample_type),
                  alpha=0.5) +
  scale_fill_brewer(palette="Set3")
ggsave(paste0(out_dir, "pca_plot.group.png"), width=7, height=5)

ggplot(pca_vals,
       aes(x=PC1, y=PC2, fill=Internal.Treatment)) +
  geom_point(shape=21, size=3) +
  theme_bw() +
  geom_text_repel(aes(label=sample_type),
                  alpha=0.5) +
  scale_fill_brewer(palette="Set3")
ggsave(paste0(out_dir, "pca_plot.internal_treatment.png"), width=7, height=5)

ggplot(pca_vals,
       aes(x=PC1, y=PC2, fill=run_batch)) +
  geom_point(shape=21, size=3) +
  theme_bw() +
  geom_text_repel(aes(label=sample_type),
                  alpha=0.5) +
  scale_fill_brewer(palette="Set1")
ggsave(paste0(out_dir, "pca_plot.run_batch.png"), width=7, height=5)

# output PCA split
pca_vals <- pca_vals[mixedorder(pca_vals$baf),]

pca_vals_left <- pca_vals[pca_vals$PC1 < 0,]
pca_vals_right <- pca_vals[pca_vals$PC1 > 0,]

write.xlsx(list(PC1_LEFT=pca_vals_left,
                PC1_RIGHT=pca_vals_right),
           file=paste0(out_dir, "metadata.pc1_split.xlsx"),
           colWidths="auto")

# correlation between PC1 and injection order
cor.test(pca_vals$PC1, pca_vals$`LC-MS_pos`, method="spearman")

ggplot(pca_vals,
       aes(y=PC1,
           x=`LC-MS_pos`)) +
  geom_point() +
  theme_bw()

# NO CORRELATION

# check to see correlation between PC1 and total intensity

total_intensity <- colSums(n_matrix, na.rm = T)
max_intensity <- apply(n_matrix, 2, max)

pca_vals$total_intensity <- total_intensity[rownames(pca_vals)]
pca_vals$max_intensity <- max_intensity[rownames(pca_vals)]

cor.test(pca_vals$PC1, pca_vals$total_intensity, method = "pearson")
cor.test(pca_vals$PC1, pca_vals$max_intensity, method="pearson")

cor.test(pca_vals$total_intensity, pca_vals$max_intensity)

ggplot(pca_vals,
       aes(x=total_intensity,
           y=PC1)) +
  geom_point(shape=21, size=3, color="black", fill="grey") +
  theme_bw() +
  labs(x="Total Intensity")
ggsave(paste0(out_dir, "pc1_total_intensity_correlation.png"), width=6, height=5)



# sample distance matrix

sample_dist <- dist(t(n_matrix_log), method = "euclidean")
sample_dist_mat <- as.matrix(sample_dist)

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

annotation_col <- data.frame(Group=metadata[colnames(sample_dist_mat),]$Group)

colnames(sample_dist_mat) <- paste0(metadata[colnames(sample_dist_mat),]$baf, "|",
                                    metadata[colnames(sample_dist_mat),]$sample_type)

rownames(sample_dist_mat) <- paste0(metadata[rownames(sample_dist_mat),]$baf, "|",
                                    metadata[rownames(sample_dist_mat),]$sample_type)

rownames(annotation_col) <- rownames(sample_dist_mat)


pheatmap(
  sample_dist_mat,
  clustering_distance_rows = sample_dist,
  clustering_distance_cols = sample_dist,
  clustering_method = "complete",   # or "average", "ward.D2"
  color = colors,
  main = "Sample Distance Heatmap (Euclidean, log2-transformed)",
  annotation_col = annotation_col,   # uncomment if using metadata
  annotation_row = annotation_col,
  fontsize_row = 8,
  fontsize_col = 8,
  filename = paste0(out_dir, "sample_distance_heatmap.group.pdf"),
  width = 8, height = 7
)
dev.off()

# trying out pqn normalization

sample_classes <- metadata[colnames(n_matrix),]$sample_type

# this seems to transpose the matrix
n_matrix_pqn <- pqn_normalisation(n_matrix,
                                    sample_classes,
                                    qc_label="all")

pqn_pseudocount <- min(n_matrix_pqn[n_matrix_pqn > 0], na.rm = TRUE) / 2

n_matrix_pqn_log <- log2(n_matrix_pqn + pseudocount)

n_pca_pqn <- prcomp(t(n_matrix_pqn_log), scale. = T)

pqn_pca_vals <- data.frame(baf=rownames(n_pca_pqn$x),
                       PC1=n_pca_pqn$x[,1],
                       PC2=n_pca_pqn$x[,2])

pqn_pca_vals <- merge(pqn_pca_vals,
                  metadata,
                  by="baf") 

rownames(pqn_pca_vals) <- pqn_pca_vals$baf

ggplot(pqn_pca_vals,
       aes(x=PC1, y=PC2, fill=sample_type)) +
  geom_point(shape=21, size=3) +
  theme_bw() +
  geom_text_repel(aes(label=sample_type),
                  alpha=0.5, max.overlaps = 100) +
  scale_fill_brewer(palette="Set3")
ggsave(paste0(out_dir, "pqn_pca_plot.sample_type.png"), width=7, height=5)

ggplot(pqn_pca_vals,
       aes(x=PC1, y=PC2, fill=run_batch)) +
  geom_point(shape=21, size=3) +
  theme_bw() +
  geom_text_repel(aes(label=sample_type),
                  alpha=0.5, max.overlaps = 100) +
  scale_fill_brewer(palette="Set1")

total_intensity_pqn <- colSums(n_matrix_pqn, na.rm = TRUE)

cor.test(pqn_pca_vals[names(total_intensity_pqn),]$PC1,
         total_intensity_pqn, method="pearson")

pqn_pca_vals$total_intensity <- total_intensity_pqn[rownames(pqn_pca_vals)]

ggplot(pqn_pca_vals,
       aes(x=total_intensity,
           y=PC1)) +
  geom_point(shape=21, size=3, color="black", fill="grey") +
  theme_bw() +
  labs(x="Total Intensity")
ggsave(paste0(out_dir, "pc1_total_intensity_correlation.post_pqn.png"), width=6, height=5)

ggplot(pqn_pca_vals,
       aes(x=run_batch,
           y=total_intensity)) +
  geom_boxplot() + theme_bw()

lipid_loadings <- n_pca_pqn$rotation

cor.test(pqn_pca_vals$PC1, pqn_pca_vals$`LC-MS_pos`, method="pearson")

loadings_pc1 <- n_pca_pqn$rotation[, "PC1"]
sorted <- sort(abs(loadings_pc1), decreasing = TRUE)
plot(sorted, type = "b", xlab = "Rank", ylab = "Abs(Loading)",
     main = "PC1 loading magnitude by rank")
head(sorted, 20)


summary(n_pca_pqn)$importance[2, "PC1"]
summary(n_pca)$importance[2, "PC1"]

cor.test(n_pca$x[,"PC1"], n_pca_pqn$x[,"PC1"])

hist(n_pca_pqn$x[, "PC1"], breaks = 30)
plot(density(n_pca_pqn$x[, "PC1"]))

pqn_pca_vals$pc1_clump <- ifelse(pqn_pca_vals$PC1 < 0, "ClumpA", "ClumpB")

table(pqn_pca_vals$pc1_clump,
      pqn_pca_vals$Group)

table(pqn_pca_vals$pc1_clump,
      pqn_pca_vals$sample_type)

table(pqn_pca_vals$pc1_clump,
      pqn_pca_vals$Internal.Treatment)

table(pqn_pca_vals$pc1_clump,
      pqn_pca_vals$run_batch)

ggplot(pqn_pca_vals,
       aes(x=PC1,
           y=`LC-MS_pos`)) +
  geom_point() + theme_bw()

# do fold change analysis

#View(unique(metadata[,c("Group-detail","Internal.Treatment","sample_type")]))

contrasts <- c("E4e1-E4c1",
               "E4e2-E4c2",
               "E3e1-E3c1",
               "E2e1-E2c1",
               "KOe1-KOc1",
               "E4e2-E4c1",
               "E4c2-E4c1",
               "KOe1-E4e1",
               "KOe1-E3e1",
               "KOe1-E2e1",
               "KOc1-E4c1",
               "KOc1-E3c1",
               "KOc1-E2c1",
               "E3e1-E4e1",
               "E2e1-E4e1",
               "E3e1-E2e1")

results_list <- lapply(contrasts, function(contrast) {
  
  print(contrast)
  
  group1 <- unlist(strsplit(contrast, "-"))[1]
  group2 <- unlist(strsplit(contrast, "-"))[2]
  
  group1_n <- n_matrix_pqn_log[,rownames(metadata[metadata$sample_type == group1,])]
  group2_n <- n_matrix_pqn_log[,rownames(metadata[metadata$sample_type == group2,])]
  
  group1_avg <- rowMeans(group1_n)
  group2_avg <- rowMeans(group2_n)
  
  log2fc <- group1_avg - group2_avg # because it's already log, subtract instead of divide
  
  pvals <- sapply(rownames(n_matrix_pqn_log), function(i) {
    
    if (length(unique(c(group1_n[i,],group2_n[i,]))) == 1) {
      return(NA)
    }
    
    t.test(log2(group1_n[i,]),
           log2(group2_n[i,]),
           alternative = "two.sided",
           var.equal = T)$p.value
    
  })
  
  adj_pvals <- p.adjust(pvals, method="fdr")
  
  # output the results
  out <- data.frame(metabolite_xref[rownames(n_matrix_pqn_log),],
             contrast = contrast,
             group1_avg=group1_avg,
             group2_avg=group2_avg,
             log2fc=log2fc,
             pvalue=pvals,
             adj_pvalue=adj_pvals)
  out[order(out$pvalue),]
})
names(results_list) <- contrasts

write.xlsx(results_list, file=paste0(out_dir, "comparison_results.xlsx"),
           colWidths="auto")

write.xlsx(metadata, file=paste0(out_dir, "metadata.xlsx"))

# volcano plots

# one for each contrast

for (results in results_list) {
  
  contrast <- unique(results$contrast)
  
  results$log_p <- -log10(results$pvalue)
  
  sig_results <- results[results$pvalue < 0.05 &
                         abs(results$log2fc) > 0.5,]
  
  ggplot(results,
         aes(x=log2fc,
             y=log_p)) +
    geom_point(alpha=0.4, color="black") +
    geom_hline(yintercept = -log10(0.05),
               color="red", linetype=2) +
    geom_vline(xintercept = 0.5,
               color="red", linetype=2) +
    geom_vline(xintercept = -0.5,
               color="red", linetype=2) +
    geom_point(data=sig_results,
               color="red",
               alpha=0.4) +
    geom_text_repel(data=sig_results,
                    aes(label=Metabolite.name),
                    color="red", size=2.5,
                    max.overlaps = 50) +
    theme_bw() +
    labs(x="Log2 Fold Change", y="-log10(P-Value)", 
         title=contrast)
  ggsave(paste0(out_dir, contrast, ".volcano.png"), width=8, height=6)
  
}

# volcano for each ontological group

for (results in results_list) {
  
  contrast <- unique(results$contrast)
  
  results$log_p <- -log10(results$pvalue)
  sig_results <- results[results$pvalue < 0.05 &
                           abs(results$log2fc) > 0.5,]
  
  # make a directory for the volcanoes
  volcano_dir <- paste0(out_dir, contrast, ".ontology_volcanoes/")
  dir.create(volcano_dir, showWarnings = F)
  
  ontologies <- unique(results$Ontology)
  
  for (ontology in ontologies) {
    
    subset <- results[results$Ontology == ontology,]
    subset_sig <- sig_results[sig_results$Ontology == ontology,]
    
    ggplot(subset,
           aes(x=log2fc,
               y=log_p)) +
      geom_point(alpha=0.4, color="black") +
      geom_hline(yintercept = -log10(0.05),
                 color="red", linetype=2) +
      geom_vline(xintercept = 0.5,
                 color="red", linetype=2) +
      geom_vline(xintercept = -0.5,
                 color="red", linetype=2) +
      geom_point(data=subset_sig,
                 color="red",
                 alpha=0.4) +
      geom_text_repel(data=subset_sig,
                      aes(label=Metabolite.name),
                      color="red", size=2.5,
                      max.overlaps = 50) +
      theme_bw() +
      labs(x="Log2 Fold Change", y="-log10(P-Value)", 
           title=contrast, subtitle=ontology)
    ggsave(paste0(volcano_dir, contrast, ".", ontology, ".volcano.png"), width=5, height=4)
    
  }
    
}



