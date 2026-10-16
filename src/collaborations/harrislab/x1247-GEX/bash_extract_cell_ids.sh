G=fastq_subset/SRR8281306_1.fastq.gz
F=fastq/SRR8281307_1.fastq.gz

for pair in "GEX:$G" "FB:$F"; do
name=${pair%%:*}; file=${pair#*:}
  gzip -dc "$file" | head -n 8000000 | awk 'NR%4==2{print substr($0,1,16)}' \
  | sort | uniq -c | sort -nr | head -n 20000 | awk '{print $2}' > top_${name}.txt
  done