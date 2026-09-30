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
library(limma)
library(reshape2)
library(dplyr)

root_dir <- "~/Documents/projects/ferrislab/cheyenne_harvey_stolz/"

out_dir <- paste0(root_dir, "results/comparisons_with_limma.PC1_batch_modeled/")
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



n_matrix <- read.xlsx(paste0(root_dir, "BAF_CH_4681_Data.xlsx"), 
                      sheet = "p-value", startRow = 5, cols=1:51)


# reformat n matrix

metabolite_xref <- n_matrix[,1:2]

metabolite_xref$metabolite_id <- 1:nrow(metabolite_xref)

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

# PQN normalization (since intensity was correlated with PC1)

sample_classes <- metadata[colnames(n_matrix),]$sample_type

mat_pqn <- pqn_normalisation(n_matrix,
                                  sample_classes,
                                  qc_label="all")

pseudocount <- min(mat_pqn[mat_pqn > 0], na.rm = TRUE) / 2

log_mat_pqn <- log2(mat_pqn + pseudocount)

# PCA

pca_pqn <- prcomp(t(log_mat_pqn), scale. = T)

pqn_pca_vals <- data.frame(baf=rownames(pca_pqn$x),
                           PC1=pca_pqn$x[,1],
                           PC2=pca_pqn$x[,2])

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
       aes(x=PC1, y=PC2, fill=Group)) +
  geom_point(shape=21, size=3) +
  theme_bw() +
  geom_text_repel(aes(label=sample_type),
                  alpha=0.5, max.overlaps = 100) +
  scale_fill_brewer(palette="Set1")
ggsave(paste0(out_dir, "pqn_pca_plot.group.png"), width=7, height=5)

ggplot(pqn_pca_vals,
       aes(x=PC1, y=PC2, fill=run_batch)) +
  geom_point(shape=21, size=3) +
  theme_bw() +
  geom_text_repel(aes(label=sample_type),
                  alpha=0.5, max.overlaps = 100) +
  scale_fill_brewer(palette="Set1")

# add in PC1 "clump"
pqn_pca_vals$pc1_clump <- factor(ifelse(pqn_pca_vals$PC1 < 0, "ClumpA", "ClumpB"))

# output to excel
pqn_pca_vals <- pqn_pca_vals[mixedorder(pqn_pca_vals$baf),]

write.xlsx(list("ClumpA"=pqn_pca_vals[pqn_pca_vals$pc1_clump == "ClumpA",],
                "ClumpB"=pqn_pca_vals[pqn_pca_vals$pc1_clump == "ClumpB",]),
           file=paste0(out_dir, "samples.pc1_clump_split.xlsx"),
           colWidths="auto")

# add clump to metadata
metadata$pc1_clump <- pqn_pca_vals[rownames(metadata),]$pc1_clump

# clump plot
ggplot(pqn_pca_vals,
       aes(x=PC1, y=PC2, fill=pc1_clump)) +
  geom_point(shape=21, size=3) +
  theme_bw() +
  geom_text_repel(aes(label=sample_type),
                  alpha=0.5, max.overlaps = 100) 
ggsave(paste0(out_dir, "pqn_pca_plot.pc1_clump.png"), width=7, height=5)

# try PCA after removing batch effect
design <- model.matrix(~ sample_type, data = metadata)

log_mat_pqn_corrected <- removeBatchEffect(
  log_mat_pqn,
  batch = metadata$pc1_clump,   # your pc1_clump variable
  design = design
)

pca_pqn_c <- prcomp(t(log_mat_pqn_corrected), scale. = T)

pqn_pca_c_vals <- data.frame(baf=rownames(pca_pqn_c$x),
                           PC1=pca_pqn_c$x[,1],
                           PC2=pca_pqn_c$x[,2])

pqn_pca_c_vals <- merge(pqn_pca_c_vals,
                      metadata,
                      by="baf") 

rownames(pqn_pca_c_vals) <- pqn_pca_c_vals$baf

ggplot(pqn_pca_c_vals,
       aes(x=PC1, y=PC2, fill=sample_type)) +
  geom_point(shape=21, size=3) +
  theme_bw() +
  geom_text_repel(aes(label=sample_type),
                  alpha=0.5, max.overlaps = 100) +
  scale_fill_brewer(palette="Set3")
ggsave(paste0(out_dir, "pqn_pca_plot.batch_corrected.sample_type.png"), width=7, height=5)

ggplot(pqn_pca_c_vals,
       aes(x=PC1, y=PC2, fill=Group)) +
  geom_point(shape=21, size=3) +
  theme_bw() +
  geom_text_repel(aes(label=sample_type),
                  alpha=0.5, max.overlaps = 100) +
  scale_fill_brewer(palette="Set1")
ggsave(paste0(out_dir, "pqn_pca_plot.batch_corrected.group.png"), width=7, height=5)

ggplot(pqn_pca_c_vals,
       aes(x=PC1, y=PC2, fill=run_batch)) +
  geom_point(shape=21, size=3) +
  theme_bw() +
  geom_text_repel(aes(label=sample_type),
                  alpha=0.5, max.overlaps = 100) +
  scale_fill_brewer(palette="Set1")

# clump plot
ggplot(pqn_pca_c_vals,
       aes(x=PC1, y=PC2, fill=pc1_clump)) +
  geom_point(shape=21, size=3) +
  theme_bw() +
  geom_text_repel(aes(label=sample_type),
                  alpha=0.5, max.overlaps = 100) 
ggsave(paste0(out_dir, "pqn_pca_plot.batch_corrected.pc1_clump.png"), width=7, height=5)


# begin differential expression analysis

# ensure same order
metadata <- metadata[colnames(log_mat_pqn),] 

all_results <- list()

# ---- 1. Groups WITH batch covariate: E2, KO ----
groups_with_batch <- c("E2", "KO")

for (grp in groups_with_batch) {
  
  meta_sub <- metadata[metadata$Group == grp, ]
  mat_sub  <- log_mat_pqn[, rownames(meta_sub)]
  
  # Ensure treatment and batch are factors with sensible reference levels
  meta_sub$treatment <- factor(meta_sub$sample_type,
                               levels = sort(unique(meta_sub$sample_type)))   
  meta_sub$batch      <- factor(meta_sub$pc1_clump)
  
  contrast <- paste0(levels(meta_sub$treatment)[2], "-", levels(meta_sub$treatment)[1])
  
  design <- model.matrix(~ treatment + batch, data = meta_sub)
  
  fit  <- lmFit(mat_sub, design)
  fit  <- eBayes(fit)
  
  # Coefficient of interest = treatment effect (2nd column, after intercept)
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  res$lipid_id <- rownames(res)
  res$group    <- grp
  res$contrast <- contrast
  res$model    <- "treatment + batch"
  
  res <- res[order(res$P.Value),]
  
  all_results[[contrast]] <- res
}

# E3 ------------------------------------------------

grp <- "E3"

meta_sub <- metadata[metadata$Group == grp, ]
mat_sub  <- log_mat_pqn[, rownames(meta_sub)]

# Ensure treatment and batch are factors with sensible reference levels
meta_sub$treatment <- factor(meta_sub$sample_type,
                             levels = sort(unique(meta_sub$sample_type)))   

contrast <- paste0(levels(meta_sub$treatment)[2], "-", levels(meta_sub$treatment)[1])

design <- model.matrix(~ treatment, data = meta_sub)

fit  <- lmFit(mat_sub, design)
fit  <- eBayes(fit)

# Coefficient of interest = treatment effect (2nd column, after intercept)
res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
res$lipid_id <- rownames(res)
res$group    <- grp
res$contrast <- contrast
res$model    <- "treatment"

res <- res[order(res$P.Value),]

all_results[[contrast]] <- res

# E4 -----------------------------------------
grp <- "E4"

meta_sub <- metadata[metadata$Group == grp, ]
mat_sub  <- log_mat_pqn[, rownames(meta_sub)]

# factorize
meta_sub$treatment <- factor(meta_sub$sample_type)

# No-intercept design: one coefficient per condition, cleanest for multiple contrasts
design <- model.matrix(~ 0 + treatment, data = meta_sub)
colnames(design) <- levels(meta_sub$treatment)

fit <- lmFit(mat_sub, design)

contrast_matrix <- makeContrasts(
  "E4e1-E4c1" = E4e1 - E4c1,
  "E4e2-E4c2" = E4e2 - E4c2,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

res_trt1 <- topTable(fit2, coef = "E4e1-E4c1", number = Inf, sort.by = "none")
res_trt2 <- topTable(fit2, coef = "E4e2-E4c2", number = Inf, sort.by = "none")

res_trt1$lipid_id <- rownames(res_trt1)
res_trt1$group    <- grp
res_trt1$contrast <- "E4e1-E4c1"
res_trt1$model    <- "treatment"

res_trt1 <- res_trt1[order(res_trt1$P.Value),]

res_trt2$lipid_id <- rownames(res_trt2)
res_trt2$group    <- grp
res_trt2$contrast <- "E4e2-E4c2"
res_trt2$model    <- "treatment"

res_trt2 <- res_trt2[order(res_trt2$P.Value),]

# add to total results
all_results[["E4e1-E4c1"]] <- res_trt1
all_results[["E4e2-E4c2"]] <- res_trt2

# merge in metabolite names
all_results <- lapply(all_results, function(x) {
  
  x <- merge(metabolite_xref,
             x,
             by.x="metabolite_id",
             by.y="lipid_id")
  
  x <- x[order(x$P.Value),]
  
  return(x)
  
})


# save results

saveRDS(all_results, file=paste0(out_dir, "comparison_results.RDS"))
write.xlsx(all_results, file=paste0(out_dir, "comparison_results.xlsx"))

# merge it all up
all_results_df <- bind_rows(all_results)

# spot checks

# those with batch effects
for (grp in groups_with_batch) {
  
  meta_sub <- metadata[metadata$Group == grp, ]
  mat_sub  <- log_mat_pqn[, rownames(meta_sub)]
  
  results <- all_results_df[all_results_df$group == grp,]
  results$data_print <- paste0(results$Metabolite.name,
                               "\nlogFC: ", round(results$logFC, digits=3),
                               "\nP-Value: ", formatC(results$P.Value, format = "e", digits=2),
                               "\nAdj P: ", formatC(results$adj.P.Val, format = "e", digits=2))
  results$data_print <- factor(results$data_print, levels=results$data_print)
  
  contrast <- unique(results$contrast)
  
  design <- model.matrix(~ sample_type, data = meta_sub)
  
  mat_corrected <- removeBatchEffect(
    mat_sub,
    batch = meta_sub$pc1_clump,   # your pc1_clump variable
    design = design
  )
  
  mat_sub_long <- melt(as.matrix(mat_sub))
  colnames(mat_sub_long) <- c("metabolite_id", "baf", "value")
  
  mat_corrected_long <- melt(as.matrix(mat_corrected))
  colnames(mat_corrected_long) <- c("metabolite_id", "baf", "value")
  
  # merge in metadata and DE results
  mat_sub_long <- merge(mat_sub_long,
                        metadata, by="baf")
  
  mat_sub_long <- merge(mat_sub_long,
                        results, by="metabolite_id")
  
  mat_corrected_long <- merge(mat_corrected_long,
                        metadata, by="baf")
  
  mat_corrected_long <- merge(mat_corrected_long,
                        results, by="metabolite_id")
  
  # start the plots
  top_metabolites <- results$metabolite_id[1:9]
  
  ggplot(mat_sub_long[mat_sub_long$metabolite_id %in% top_metabolites,],
         aes(x=sample_type,
             y=value)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    theme(legend.position = "bottom") +
    facet_wrap(~ data_print, ncol=3, scales="free_y") +
    geom_jitter(shape=21, color="black", aes(fill=pc1_clump), width = 0.1, size=3) +
    labs(x=NULL, y="Log Scaled PQN Expression", fill="PC1 Clump")
  ggsave(paste0(out_dir, contrast, ".top_metabolites.raw_expression.png"), width=9, height=8)
  
  ggplot(mat_corrected_long[mat_corrected_long$metabolite_id %in% top_metabolites,],
         aes(x=sample_type,
             y=value)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    theme(legend.position = "bottom") +
    facet_wrap(~ data_print, ncol=3, scales="free_y") +
    geom_jitter(shape=21, color="black", aes(fill=pc1_clump), width = 0.1, size=3) +
    labs(x=NULL, y="Log Scaled PQN Expression\nBatch Corrected", fill="PC1 Clump")
  ggsave(paste0(out_dir, contrast, ".top_metabolites.batch_corrected_expression.png"), width=9, height=8)
  
  
}


# do it for the remaining contrasts as well
non_batch_contrasts <- unique(all_results_df[!all_results_df$group %in% groups_with_batch,]$contrast)

for (contrast in non_batch_contrasts) {
  
  samples <- unlist(strsplit(contrast, "-"))
  
  meta_sub <- metadata[metadata$sample_type %in% samples, ]
  mat_sub  <- log_mat_pqn[, rownames(meta_sub)]
  
  results <- all_results_df[all_results_df$contrast == contrast,]
  results$data_print <- paste0(results$Metabolite.name,
                               "\nlogFC: ", round(results$logFC, digits=3),
                               "\nP-Value: ", formatC(results$P.Value, format = "e", digits=2),
                               "\nAdj P: ", formatC(results$adj.P.Val, format = "e", digits=2))
  results$data_print <- factor(results$data_print, levels=results$data_print)
  
  mat_sub_long <- melt(as.matrix(mat_sub))
  colnames(mat_sub_long) <- c("metabolite_id", "baf", "value")
  
  
  # merge in metadata and DE results
  mat_sub_long <- merge(mat_sub_long,
                        metadata, by="baf")
  
  mat_sub_long <- merge(mat_sub_long,
                        results, by="metabolite_id")

  
  # start the plots
  top_metabolites <- results$metabolite_id[1:9]
  
  ggplot(mat_sub_long[mat_sub_long$metabolite_id %in% top_metabolites,],
         aes(x=sample_type,
             y=value)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    facet_wrap(~ data_print, ncol=3, scales="free_y") +
    geom_jitter(shape=21, color="black", fill="grey", width = 0.1, size=3) +
    labs(x=NULL, y="Log Scaled PQN Expression",)
  ggsave(paste0(out_dir, contrast, ".top_metabolites.raw_expression.png"), width=9, height=8)
  
  
  
}


for (results in all_results) {
  
  contrast <- unique(results$contrast)
  
  results$log_p <- -log10(results$P.Value)
  
  sig_results <- results[results$adj.P.Val < 0.05 &
                           abs(results$logFC) > 0.5,]
  
  fdr_line <- min(sig_results$log_p)
  
  ggplot(results,
         aes(x=logFC,
             y=log_p)) +
    geom_point(alpha=0.4, color="black") +
    geom_hline(yintercept = fdr_line,
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







