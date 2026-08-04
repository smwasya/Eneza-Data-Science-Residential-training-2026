while read SRR
do
    echo "Downloading $SRR"
    prefetch "$SRR"
    fasterq-dump "$SRR" --split-files -e 8
    pigz -p 8 ${SRR}*.fastq
done < ids.txt
