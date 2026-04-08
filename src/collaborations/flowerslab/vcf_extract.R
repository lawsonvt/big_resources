library(vcfR)
library(openxlsx)

vcf_file <- read.vcfR("~/OneDrive - University of Virginia/projects/flowerslab/dana_may/KOLF2_SNPs/ABCA7_snps.vcf")


vcf_file <- read.delim("~/OneDrive - University of Virginia/projects/flowerslab/dana_may/KOLF2_SNPs/ABCA7_snps.vcf",
                       header=F)

vcf_header <- "CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HPSI0114i-kolf_2"
vcf_header <- unlist(strsplit(vcf_header, "\\s+"))

colnames(vcf_file) <- vcf_header

id_ones <- vcf_file[vcf_file$ID != ".",]

non_non_ref <- vcf_file[vcf_file$ALT != "<NON_REF>",]



write.xlsx(list(total=vcf_file,
                alternate=non_non_ref,
                rsids=id_ones), file="~/OneDrive - University of Virginia/projects/flowerslab/dana_may/KOLF2_SNPs/ABCA7_snps.vcf.xlsx",
           colWidths="auto")

