# parameters
read_dir <- "/project/SheybaniLab/Stefanyda_Maslova/RNA\ Seq/Exp57_D7/SheybaniLab-483044923/BCLConvert_03_12_2026_05_55_30Z-906849967/"
fastp_location <- "~/fastp"

fastq1_end <- "R1_001.fastq.gz"
fastq2_end <- "R2_001.fastq.gz"

# get the fastq files
fastq1_files <- list.files(read_dir,
                           fastq1_end, full.names=T, recursive = T)

fastq2_files <- list.files(read_dir,
                           fastq2_end, full.names=T, recursive = T)


# make fastp out files
fp1_files <- gsub(".fastq.gz", ".fp_out.fastq.gz", fastq1_files)
fp2_files <- gsub(".fastq.gz", ".fp_out.fastq.gz", fastq2_files)


for (i in 1:length(fastq1_files)) {
  
  command <- paste0(fastp_location, 
                    " -i ", fastq1_files[i], 
                    " -I ", fastq2_files[i],
                    " -o ", fp1_files[i],
                    " -O ", fp2_files[i])
  
  print(command)
  
  system(command)
  
}
