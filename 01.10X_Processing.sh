#!/bin/bash

# Set environment variables here
X10_ROOT=""
OUT_DATA=""

# Create necessary directories
for dir in 00_clean 01_remove_bar 02_remove_bar_lambda 03_filter_reads 05_regen_clean; do
    mkdir -p "${OUT_DATA}/${dir}"
done

# Quality control with fastp
cd "${X10_ROOT}" || exit
for file in ./${CELL}*/*_L*_1.fq.gz; do
    p=$(dirname "$file")
    i=$(basename "$file" _1.fq.gz)
    new_file_name=$(echo "$i" | awk -F'_' '{print $1"_"$2"_"$NF}')
    
    echo "Processing ${new_file_name}"
    /home/users/ddu/software/Fastp/fastp -A -G -q 19 -u 50 -l 80 -n 7 -M 0 -w 10 -c \
        -i "${p}/${i}_1.fq.gz" -o "${OUT_DATA}/00_clean/${new_file_name}_1.fastq.gz" \
        -I "${p}/${i}_2.fq.gz" -O "${OUT_DATA}/00_clean/${new_file_name}_2.fastq.gz" \
        -h "${OUT_DATA}/00_clean/${i}.html" -j "${OUT_DATA}/00_clean/${i}.json"
    echo "fastp completed for ${new_file_name}"
done

echo "All fastp processing done."

# Removing barcode
cd "${OUT_DATA}" || exit
source activate python27

for i in $(ls 00_clean/*_L*_1.fastq.gz | xargs -n1 basename | sed 's/_1.fastq.gz//'); do
    /home/users/ddzhang/miniconda3/envs/python2/bin/python \
        /home/users/xyxu/tools/proc10xG/proc10xG/process_10xReads.py -a \
        -o "01_remove_bar/${i}" -1 "00_clean/${i}_1.fastq.gz" -2 "00_clean/${i}_2.fastq.gz"
    echo "Barcode removal done for ${i}"
done

echo "All barcode removals completed."

# Aligning reads to lambda genome
for i in $(ls 01_remove_bar/*_R1_001.fastq.gz | xargs -n1 basename | sed 's/_R1_001.fastq.gz//'); do
    /home/apps/bwa-0.7.15/bwa mem -t 20 /home/users/xyxu/assembly/rawdata/lambda/genome.fa \
        "01_remove_bar/${i}_R1_001.fastq.gz" "01_remove_bar/${i}_R2_001.fastq.gz" | \
        samtools sort --output-fmt BAM -T "02_remove_bar_lambda/${i}" -o "02_remove_bar_lambda/${i}.sort.bam"
    echo "Alignment done for ${i}"
done

echo "Lambda filtering completed."

# Extract unmapped reads
for i in $(ls 02_remove_bar_lambda/*.sort.bam | xargs -n1 basename | sed 's/.sort.bam//'); do
    samtools view -f4 "02_remove_bar_lambda/${i}.sort.bam" | awk '{print $1}' | sort -u > "02_remove_bar_lambda/${i}.unmapped.reads"
    echo "Unmapped reads extracted for ${i}"

done

echo "Unmapped read extraction done."

# Filter reads based on extracted unmapped reads
for i in $(ls 02_remove_bar_lambda/*.unmapped.reads | xargs -n1 basename | sed 's/.unmapped.reads//'); do
    /home/users/ddu/software/seqkit grep --threads 40 -f "02_remove_bar_lambda/${i}.unmapped.reads" \
        "01_remove_bar/${i}_R1_001.fastq.gz" | gzip - > "03_filter_reads/${i}_R1_001.fastq.gz"
    /home/users/ddu/software/seqkit grep --threads 40 -f "02_remove_bar_lambda/${i}.unmapped.reads" \
        "01_remove_bar/${i}_R2_001.fastq.gz" | gzip - > "03_filter_reads/${i}_R2_001.fastq.gz"
    echo "Filtering completed for ${i}"
done

echo "All filtering steps completed."

# Compute barcode statistics
cat 02_remove_bar_lambda/*unmapped.reads | awk -F':' '{print $1}' | sort | uniq -c | awk '$1>=3 {print $2}' > "03_filter_reads/barcode_filtered.txt"

# Compute enrichment percentage
a=$(cat 02_remove_bar_lambda/*unmapped.reads | wc -l)
b=$(cat 02_remove_bar_lambda/*.mapping.reads | wc -l)
c=$(echo "scale=4; $a / ($a + $b) * 100" | bc)
echo -e "${CELL}\t$a\t$b\t$c" >> "03_filter_reads/sample.mhc.reads.percent.txt"

echo "Enrichment calculation completed."

# Restore 10X format
for i in $(ls 03_filter_reads/*_R1_001.fastq.gz | xargs -n1 basename | sed 's/_R1_001.fastq.gz//'); do
    /home/users/ddzhang/miniconda3/envs/python2/bin/python /home/users/xyxu/tools/proc10xG/proc10xG/filter_10xReads.py \
        -L "03_filter_reads/barcode_filtered.txt" \
        -1 "03_filter_reads/${i}_R1_001.fastq.gz" \
        -2 "03_filter_reads/${i}_R2_001.fastq.gz" \
        -o "04_filter_low/${i}"
    echo "Filtering completed for ${i}"
    
    /home/users/ddzhang/miniconda3/envs/python2/bin/python /home/users/xyxu/tools/proc10xG/proc10xG/regen_10xReads.py \
        -o "05_regen_clean/${i}" \
        -1 "04_filter_low/${i}_R1_001.fastq.gz" \
        -2 "04_filter_low/${i}_R2_001.fastq.gz"
    echo "Reconstruction completed for ${i}"
done

echo "Filtering and reconstruction completed for all samples."

# Align 10X reads to GRCh38 reference for calling variants
longranger wgs \
    --id=longranger_wgs \
    --reference=/home/users/xyxu/assembly/rawdata/ref/refdata-GRCh38-2.1.0/ \
    --fastqs=${OUT_DATA}/05_regen_clean \
    --sample=${CELL} \
    --nopreflight \
    --localcores=40 --localmem=400 \
    --vcmode=gatk:/home/users/xyxu/tools/GATK/gatk3.8/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
    >Extract_reads_by_barcode_all_wgs.log