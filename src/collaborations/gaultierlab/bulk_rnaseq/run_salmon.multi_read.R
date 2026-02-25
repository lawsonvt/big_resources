# parameters
root_dir <- "~/projects/gaultierlab/stephanie_moy/bulkRNASeq/raw_reads/"
salmon_index <- "~/projects/gaultierlab/stephanie_moy/bulkRNASeq/mm_GRCm39_index/"

fastq1_end <- "R1_001.fastq.gz"
fastq2_end <- "R2_001.fastq.gz"

# get the fastq files
fastq1_files <- list.files(root_dir,
                           fastq1_end, full.names=T, recursive = T)

fastq2_files <- list.files(root_dir,
                           fastq2_end, full.names=T, recursive = T)

# deal with multiple reads per sample
sample_dirs <- unique(dirname(fastq1_files))

# output directories
quant_dirs <- paste0(sample_dirs, "quant/")

# run the salmon command
for (sample_dir in sample_dirs) {
  
  print(basename(sample_dir))
  
  sample_fastq1 <- list.files(sample_dir,
                              fastq1_end, full.names=T)
  
  sample_fastq2 <- list.files(sample_dir,
                              fastq2_end, full.names = T)
  
  quant_dir <- paste0(sample_dir, "quant/")
  
  command <- paste0("salmon quant",
                    " -i ", salmon_index,
                    " -l A",
                    " -1 ", paste0(sample_fastq1, collapse=" "),
                    " -2 ", paste0(sample_fastq2, collapse=" "),
                    " -p 8 --validateMappings",
                    " -o ", quant_dir)
  
  print(command)
  system(command)
}


