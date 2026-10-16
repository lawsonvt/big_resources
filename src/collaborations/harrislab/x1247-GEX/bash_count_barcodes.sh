gzip -dc x1247-GEX-ds.*/x1247-GEX_S1_L001_R1_001.fastq.gz | \
  awk 'NR%4==2{print substr($0,1,16)}' | \
  sort | uniq -c | sort -nr > gex_barcode_counts.txt
  
  gzip -dc fastq/SRR8281307_1.fastq.gz | \
  awk 'NR%4==2{print substr($0,1,16)}' | \
  sort | uniq -c | sort -nr > fb_barcode_counts.txt