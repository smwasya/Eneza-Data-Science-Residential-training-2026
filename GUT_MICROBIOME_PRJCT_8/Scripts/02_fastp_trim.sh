mkdir -p /home/sam/ENEZA/results/fastp_results

for R1 in /home/sam/gut_uganda/FASTQ/*_1.fastq.gz
do
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"
    SAMPLE=$(basename "$R1" _1.fastq.gz)

    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "/home/sam/ENEZA/results/fastp_results/${SAMPLE}_1.fastq.gz" \
        -O "/home/sam/ENEZA/results/fastp_results/${SAMPLE}_2.fastq.gz" \
        -h "/home/sam/ENEZA/results/fastp_results/${SAMPLE}.html" \
        -j "/home/sam/ENEZA/results/fastp_results/${SAMPLE}.json" \
        --thread 3

done
