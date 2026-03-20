# parameters
root_dir <- "~/projects/coutinholab/flyseq/raw_reads/"
salmon_index <- "~/projects/coutinholab/flyseq/dm_BDGP6_index/"

fastq1_end <- "R1_001.fp_out.fastq.gz"
fastq2_end <- "R2_001.fp_out.fastq.gz"

# get the fastq files
fastq1_files <- list.files(root_dir,
                           fastq1_end, full.names=T, recursive = T)

fastq2_files <- list.files(root_dir,
                           fastq2_end, full.names=T, recursive = T)

# output directories
quant_dir <- paste0(dirname(fastq1_files), "quant/")

# run the salmon command
for (i in 1:length(quant_dir)) {
  
  command <- paste0("salmon quant",
                    " -i ", salmon_index,
                    " -l A",
                    " -1 ", fastq1_files[i],
                    " -2 ", fastq2_files[i],
                    " -p 8 --validateMappings",
                    " -o ", quant_dir[i])
  
  print(command)
  system(command)
}


