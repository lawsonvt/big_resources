library(DESeq2)
library(tximport)
library(ggplot2)
library(openxlsx)
library(data.table)
library(ggrepel)
library(reshape2)
library(betareg)
library(openxlsx)

out_dir <- "~/Documents/projects/coutinholab/flyseq/results/isoform_expression/"

dir.create(out_dir, recursive = T, showWarnings = F)

# get the quant files
quant_files <- list.files("~/Documents/projects/coutinholab/flyseq/quants/",
                          "quant.sf",
                          recursive=T, full.names = T)

names(quant_files) <- sapply(quant_files, function(x) {
  
  dirs <- unlist(strsplit(dirname(x), "\\/"))
  
  return(gsub("quant", "", last(dirs)))
  
})

# get tx2gene
tx2gene <- read.csv("~/Documents/projects/coutinholab/flyseq/Drosophila_melanogaster.BDGP6.54.transcript2gene.csv")

# read in quant files
trans_txi <- tximport(quant_files, type="salmon", txOut = T,
                tx2gene = tx2gene[,c("transcript_id","gene_id")])


trans_tpm <- trans_txi$abundance
trans_counts <- trans_txi$counts

# plot out orion TPM scores

orion_transcripts <- tx2gene[tx2gene$gene_symbol == "orion",]$transcript_id

orion_tpm <- trans_tpm[orion_transcripts,]

orion_isoform_fraction <- orion_tpm[orion_transcripts[1],] / (orion_tpm[orion_transcripts[1],] +
                                                                orion_tpm[orion_transcripts[2],])

orion_tpm_long <- melt(orion_tpm)
colnames(orion_tpm_long) <- c("transcript_id","sample","tpm")

orion_tpm_long$condition <- gsub("[0-9]+", "", orion_tpm_long$sample)

ggplot(orion_tpm_long,
       aes(x=condition,
           y=tpm)) +
  facet_wrap(~ transcript_id, ncol=2) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  geom_text_repel(aes(label=sample)) +
  theme_bw() +
  labs(x=NULL, y="TPM")
ggsave(paste0(out_dir, "orion_tpm_comp.png"), width=7, height=5)

orion_if_long <- melt(orion_isoform_fraction)
orion_if_long$sample <- rownames(orion_if_long)
orion_if_long$condition <- gsub("[0-9]+", "", orion_if_long$sample)

ggplot(orion_if_long,
       aes(x=condition,
           y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  theme_bw() +
  geom_hline(yintercept = 0.5, linetype=2) +
  ylim(0, 1) +
  labs(x=NULL, y=paste0(orion_transcripts[1]," fraction"))
ggsave(paste0(out_dir, "orion_isoform_fraction.png"), width=4, height=5)

# do a differential expression analysis



fit <- betareg(value ~ condition, data=orion_if_long)
summary(fit)

# save to excel file
write.xlsx(list(tpm=orion_tpm_long,
     orionb_fraction=orion_if_long),
     paste0(out_dir, "orion_isoform_expression.xlsx"),
     colWidths="auto")








