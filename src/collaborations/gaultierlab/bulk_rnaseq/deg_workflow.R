library(DESeq2)
library(tximport)
library(ggplot2)
library(openxlsx)

out_dir <- "~/Documents/projects/gaultierlab/sam_wachamo/bulkRNASeq/results/deg_workflow/"

dir.create(out_dir, recursive = T, showWarnings = F)

# get metadata
metadata <- read.xlsx("~/Documents/projects/gaultierlab/sam_wachamo/bulkRNASeq/metadata.xlsx")
rownames(metadata) <- metadata$sample_id

# factorize things
metadata$condition <- factor(metadata$condition,
                             levels=c("Vehicle","IAA"))
metadata$sex <- factor(metadata$sex)

# get the quant files
quant_files <- list.files("~/Documents/projects/gaultierlab/sam_wachamo/bulkRNASeq/quants/",
                          "quant.sf",
                          recursive=T, full.names = T)

names(quant_files) <- sapply(quant_files, function(x) {
  
  dirs <- unlist(strsplit(dirname(x), "\\/"))
  
  return(gsub("quant", "", last(dirs)))
  
})

# ensure that files and metadata have same order
all(names(quant_files) %in% rownames(metadata))
all(rownames(metadata) %in% names(quant_files) )

metadata <- metadata[names(quant_files),]

# get tx2gene
tx2gene <- read.csv("~/Documents/projects/gaultierlab/sam_wachamo/bulkRNASeq/Mus_musculus.GRCm39.transcript2gene.csv")

# read in quant files
txi <- tximport(quant_files, type="salmon", 
                tx2gene = tx2gene[,c("transcript_id","gene_id")])

# create DESeq object
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = metadata,
                                   design = ~ condition + sex)

# pre filter counts data (for plotting purposes)
smallestGroupSize <- 5
keep <- rowSums(counts(ddsTxi) >= 10) >= smallestGroupSize
ddsTxi <- ddsTxi[keep,]

ddsTxi <- DESeq(ddsTxi)

# QC plots
vsd <- vst(ddsTxi)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape=sex)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA") +
  scale_color_manual(values=c("orange","blue")) +
  theme_bw()
ggsave(paste0(out_dir, "pca_plot.png"), width=5, height=4)

# everything looks good, no SVA needed

# get DEG results
res <- results(ddsTxi, name="condition_IAA_vs_Vehicle")

res <- as.data.frame(res)
# remove NAs
res <- res[!is.na(res$padj),]

res$gene_id <- rownames(res)

# annotate
res <- merge(unique(tx2gene[,c("gene_id","gene_symbol")]),
             res,
             by="gene_id", all.y=T)

# order it
res <- res[order(res$pvalue),]

# save to file
write.xlsx(res, paste0(out_dir, "IAA_v_Vehicle.degs.xlsx"), colWidths="auto")
saveRDS(res, paste0(out_dir, "IAA_v_Vehicle.degs.RDS"))

# make a volcano plot
res$logp <- -log10(res$pvalue)

# fix "Inf" values
max_logp <- max(res[!is.infinite(res$logp),]$logp)

res[is.infinite(res$logp),]$logp <- max_logp

# pull out significant ones
sig_res <- res[res$padj < 0.01 &
                 abs(res$log2FoldChange) > 0.5,]

sig_thresh <- min(sig_res$logp)

ggplot(res,
       aes(x=log2FoldChange,
           y=logp)) +
  geom_point(alpha=0.2) +
  geom_hline(yintercept = sig_thresh, color="red", linetype=2) +
  geom_vline(xintercept = 0.5, color="red", linetype=2) +
  geom_vline(xintercept = -0.5, color="red", linetype=2) +
  geom_point(data=sig_res, color="red", alpha=0.4) +
  theme_bw() +
  labs(x="Log2 Fold Change", y="-Log10 P-Value")
ggsave(paste0(out_dir, "volcano_plot.png"), width=7, height=6)







