library(DESeq2)
library(tximport)
library(ggplot2)
library(openxlsx)
library(data.table)
library(ggrepel)
library(readxl)
library(sva)
library(limma)
library(tibble)
library(reshape2)
library(snakecase)
library(cowplot)
library(pbayes)
library(dplyr)
library(purrr)

root_dir <- "~/Documents/projects/sheybanilab/stef_maslova/"

out_dir <- paste0(root_dir, "results/salmon_deg_analysis.direct_comparison/")
dir.create(out_dir, showWarnings = F, recursive = T)

# read in data
metadata <- as.data.frame(read_excel(paste0(root_dir, "260130_Bulk_Seq_List.xls")))
colnames(metadata)[1] <- "sample_id"

metadata$sample_id <- paste0("Mouse_", metadata$sample_id)
rownames(metadata) <- metadata$sample_id

metadata$condition <- factor(metadata$condition,
                             levels=unique(metadata$condition))

# get quant files
quant_files <- list.files(paste0(root_dir, "quants/"),
                          "quant.sf",
                          recursive=T, full.names = T)

names(quant_files) <- sapply(quant_files, function(x) {
  
  dirs <- unlist(strsplit(dirname(x), "\\/"))
  
  sample_number <- as.numeric(unlist(strsplit(last(dirs), "\\-"))[1])
  
  return(paste0("Mouse_", sample_number))
  
})

# ensure that files and metadata have same order
all(names(quant_files) %in% rownames(metadata))
all(rownames(metadata) %in% names(quant_files) )

metadata_total <- metadata[names(quant_files),]


# loop through treatments, comparing each to control

# treatments <- levels(metadata$condition)
# treatments <- treatments[!grepl("Control", treatments)]

# create the contrasts
contrasts <- list("FUS-Control"=c("condition","FUS","Control"),
                  "CAR-Control"=c("condition","CAR","Control"),
                  "CAR_FUS-Control"=c("condition","CAR_FUS","Control"),
                  "CAR_FUS-CAR"=c("condition","CAR_FUS","CAR"),
                  "CAR_FUS-FUS"=c("condition","CAR_FUS","FUS"),
                  "CAR-FUS"=c("condition","CAR","FUS"))

results_list <- lapply(contrasts, function(contrast) {
  
  # filter down
  metadata <- metadata_total[metadata_total$condition %in% contrast[c(3,2)],]
  
  # refactor conditions
  metadata$condition <- factor(as.character(metadata$condition),
                               levels=contrast[c(3,2)])
  
  subset_quant <- quant_files[rownames(metadata)]
  # get tx2gene
  tx2gene <- read.csv("~/Documents/projects/gaultierlab/sam_wachamo/bulkRNASeq/Mus_musculus.GRCm39.transcript2gene.csv")
  
  # read in quant files
  txi <- tximport(subset_quant, type="salmon", 
                  tx2gene = tx2gene[,c("transcript_id","gene_id")])
  
  # create DESeq object
  ddsTxi <- DESeqDataSetFromTximport(txi,
                                     colData = metadata,
                                     design = ~ condition)
  
  # pre filter counts data (for plotting purposes)
  smallestGroupSize <- 3
  keep <- rowSums(counts(ddsTxi) >= 10) >= smallestGroupSize
  ddsTxi <- ddsTxi[keep,]
  
  ddsTxi <- DESeq(ddsTxi)
  
  # QC plots
  vsd <- vst(ddsTxi)
  
  pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA") +
    geom_text_repel(aes(label=Name)) +
    theme_bw()
  ggsave(paste0(out_dir, "pca_plot.", paste0(contrast[2], "-",
                                             contrast[3]), ".png"), width=6, height=5)
  
  # try out sva
  
  # Get normalized counts for SVA
  dds_norm <- estimateSizeFactors(ddsTxi)
  norm_counts <- counts(dds_norm, normalized = TRUE)
  
  
  # dat <- vst(dds_norm, blind = FALSE)
  # dat_matrix <- assay(dat)
  
  mod <- model.matrix(~ 0 + condition, colData(ddsTxi))
  mod0 <- model.matrix(~1, colData(ddsTxi))
  
  n_sv <- num.sv(norm_counts, mod)
  
  # Calculate surrogate variables
  svobj <- svaseq(norm_counts, mod, mod0, n.sv = n_sv)
  
  
  vsd_corrected <- vst(ddsTxi, blind = FALSE)
  assay(vsd_corrected) <- limma::removeBatchEffect(
    assay(vsd_corrected),
    covariates = svobj$sv,
    design = mod
  )
  
  pcaData_after <- plotPCA(vsd_corrected, intgroup = "condition", returnData = TRUE)
  percentVar_after <- round(100 * attr(pcaData_after, "percentVar"))
  
  ggplot(pcaData_after, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar_after[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar_after[2], "% variance")) +
    ggtitle("PCA, SVA applied") +
    geom_text_repel(aes(label=Name)) +
    theme_bw()
  ggsave(paste0(out_dir, "pca_sva_plot.", paste0(contrast[2], "-",
                                                 contrast[3]), ".png"), width=6, height=5)
  
  # SVA seems to show a difference ...
  
  # try doing differential expression analysis
  
  # Add SVs to colData
  for (i in 1:n_sv) {
    colData(ddsTxi)[, paste0("SV", i)] <- svobj$sv[, i]
  }
  
  # Create new design formula including SVs
  # Build the SV terms dynamically
  sv_terms <- paste0("SV", 1:n_sv, collapse = " + ")
  
  design_formula <- as.formula(paste("~", sv_terms, "+ condition"))
  
  # Update the design
  design(ddsTxi) <- design_formula 
  
  ddsTxi <- DESeq(ddsTxi)
  
  # # create the contrasts
  # contrast <- c("condition",treatment,"Control")
  
  print(paste0(contrast[2], "-",
               contrast[3]))
  
  # calculate results
  res <- results(ddsTxi, contrast) 
  
  res <- as.data.frame(res)
  # remove NAs
  res <- res[!is.na(res$padj),]
  
  # posterior probability
  res$post_p <- pbayes(res$pvalue, n_cores=2, level_pvals = T)$posterior_prob
  
  # merge in gene names
  res <- rownames_to_column(res, var = "gene_id")
  res <- merge(unique(tx2gene[,c("gene_id","gene_symbol")]),
               res, by="gene_id")
  
  # sort the data
  res <- res[order(res$pvalue),]
  
  res$contrast <- paste0(contrast[2], "-",
                         contrast[3])
  
  # fix NAs in gene name
  res[is.na(res$gene_symbol),]$gene_symbol <- res[is.na(res$gene_symbol),]$gene_id
  
  return(res)
  
})
names(results_list) <- names(contrasts)

# drop results with low counts
results_list <- lapply(results_list, function(data) {
  
  return(data[data$baseMean > 30,])
  
})


write.xlsx(results_list, paste0(out_dir, "deg_results.xlsx"), colWidths="auto")

# save results to RDS
saveRDS(results_list, file=paste0(out_dir, "deg_results.RDS"))

# make some plots

# Make Volcano plots -----------------------------------------------------------

volcano_dir <- paste0(out_dir, "volcano_plots/")

dir.create(volcano_dir, showWarnings = F)


top_genes <- 25

cap <- 4

volcano_plot_list <-  lapply(names(results_list), function(contrast) {
  
  print(contrast)
  
  subset <- results_list[[contrast]]
  
  subset$log_p <- -log10(subset$pvalue)
  
  # cap foldchange
  if (any(subset$log2FoldChange > cap)) {
    subset[subset$log2FoldChange > cap,]$log2FoldChange <- cap 
  }
  
  if (any(subset$log2FoldChange < -cap)) {
    subset[subset$log2FoldChange < -cap,]$log2FoldChange <- -cap 
  }
    
  subset_sig <- subset[subset$padj < 0.1,]
  
  subset_top <- subset_sig[1:top_genes,]
  
  if (any(is.na(subset_top$gene_symbol))) {
    subset_top <- subset_top[!is.na(subset_top$gene_symbol),]
  }
  
  logp_thresh <- min(subset_sig$log_p)
  
  p1 <- ggplot(subset,
               aes(x=log2FoldChange,
                   y=log_p)) +
    geom_point(alpha=0.4, color="black") +
    geom_hline(yintercept = logp_thresh,
               color="red", linetype=2) +
    # geom_vline(xintercept = 0.5,
    #            color="red", linetype=2) +
    # geom_vline(xintercept = -0.5,
    #            color="red", linetype=2) +
    geom_point(data=subset_sig,
               color="red") +
    geom_text_repel(data=subset_top,
                    aes(label=gene_symbol),
                    color="red", size=2.5,
                    max.overlaps = 50) +
    theme_bw() +
    labs(x="Log2 Fold Change", y="-log10(P-Value)", title=contrast)
  p1
  ggsave(paste0(volcano_dir, to_snake_case(contrast), ".volcano_plot.png"), width=5, height=6)
  
  return(p1)
})


plot_grid(plotlist = volcano_plot_list, nrow = 2)
ggsave(paste0(out_dir, "contrast.volcanoes.png"), width=10, height=8, bg="white")

# identify trends

# simplify the results
results_list_simple <- lapply(results_list, function(x) {
  
  contrast <- unique(x$contrast)
  
  subset <- x[,c("gene_id","gene_symbol","log2FoldChange","pvalue","padj", "post_p")]
  colnames(subset)[3:6] <- paste0(colnames(subset)[3:6], ".", to_snake_case(contrast))
  
  return(subset)
})

# merge em up

total_df <- reduce(results_list_simple, inner_join, by=c("gene_id","gene_symbol"))

max_fdr <- 0.1

sum(total_df$padj.fus_control < max_fdr)
sum(total_df$padj.car_control < max_fdr)
sum(total_df$padj.car_fus_control < max_fdr)

sum(total_df$padj.car_fus < max_fdr)
sum(total_df$padj.car_fus_car < max_fdr)
sum(total_df$padj.car_fus_fus < max_fdr)

# identify combos
sig_both_ind <- total_df[total_df$padj.car_control < max_fdr &
                           total_df$padj.fus_control < max_fdr,]

sig_car_fus_control <- total_df[total_df$padj.car_fus_control < max_fdr,]
summary(sig_car_fus)



# combo effects
sig_fus <- total_df[total_df$padj.fus_control < max_fdr &
                      total_df$padj.car_fus_fus < max_fdr,]


sig_car <- total_df[total_df$padj.car_control < max_fdr &
                      total_df$padj.car_fus_car < max_fdr,]


sig_fus_car <- total_df[total_df$padj.fus_control < max_fdr &
                          total_df$padj.car_fus < max_fdr,]

sig_car_fus <- total_df[total_df$padj.car_control < max_fdr &
                          total_df$padj.car_fus < max_fdr,]


