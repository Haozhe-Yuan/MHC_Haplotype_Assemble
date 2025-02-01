#!/bin/bash

set -e

# Define paths
HIFI_30_ROOT="/home/users/Liulab_data/xyxu/hifi/Q30/X101SC23065068-Z01-F001/Data-X101SC23065068-Z01-F001/${CELL}"
ADD_HIFI_30_ROOT="/home/users/Liulab_data/xyxu/addhifi_Clinical_Q30/X101SC23065068-Z01-F003/Data-X101SC23065068-Z01-F003/${CELL}"
OUT_DATA_ROOT="/home/users/hzyuan/workspace/cellsHLA_refine/processed_data/${CELL}"

# Create directories
for dir in hifi/q30A; do
    mkdir -p "${OUT_DATA_ROOT}/${dir}"
done

OUT_30_DATA="${OUT_DATA_ROOT}/hifi/q30A"

# Set paths based on QU value
RAW_PATH=${HIFI_30_ROOT}
RAW_ADD_PATH=${ADD_HIFI_30_ROOT}
OUT_PATH=${OUT_30_DATA}
SUFFIX=filt.hifi.bam


# Create necessary subdirectories
for subdir in 10X_hifi 10X_hifi/ref 10X_hifi/CCS_variants/vcf 10X_hifi/bam 00_Remove_Dup 01_Align imgs/evaluation; do
    mkdir -p "${OUT_PATH}/${subdir}"
done

source activate whatshap-env

# Merge original and additional HiFi data
pbmerge -o "${OUT_PATH}/00_Remove_Dup/${CELL}_Q${QU}.dup.bam" \
        "${RAW_PATH}/${CELL}.${SUFFIX}" \
        "${RAW_ADD_PATH}/${CELL}.${SUFFIX}" && \
pbindex "${OUT_PATH}/00_Remove_Dup/${CELL}_Q${QU}.dup.bam"

echo "Merge completed at $(date)"

# Remove duplicates
pbmarkdup -f -j 20 --cross-library \
    "${OUT_PATH}/00_Remove_Dup/${CELL}_Q${QU}.dup.bam" \
    "${OUT_PATH}/00_Remove_Dup/${CELL}_Q${QU}.hifi_reads.bam" \
    --dup-file "${OUT_PATH}/00_Remove_Dup/${CELL}_Q${QU}.hifi_reads.bam.fasta"

echo "Duplicates removed at $(date)"

# Alignment using minimap2
minimap2 -a -k 19 -O 5,56 -E 4,1 -B 5 -z 400,50 -r 2k -t 60 \
    -R "@RG\tSM:${CELL}\tID:${CELL}" --eqx --secondary=no \
    /home/users/xyxu/assembly/rawdata/10X_hifi/ref/genome.fa \
    <(samtools bam2fq -@ 60 "${OUT_PATH}/00_Remove_Dup/${CELL}_Q${QU}.hifi_reads.bam") | \
samtools sort -@ 60 --output-fmt BAM -o "${OUT_PATH}/10X_hifi/${CELL}.hg38.bam"

samtools index -@ 60 "${OUT_PATH}/10X_hifi/${CELL}.hg38.bam"
echo "Alignment completed at $(date)"

# Extract MHC region
samtools view -@ 60 -bh "${OUT_PATH}/10X_hifi/${CELL}.hg38.bam" chr6:28903952-33268517 > "${OUT_PATH}/10X_hifi/${CELL}.hg38.mhc.bam"
samtools index -@ 60 "${OUT_PATH}/10X_hifi/${CELL}.hg38.mhc.bam"
echo "MHC extraction completed at $(date)"

# Run DeepVariant
mkdir -p "${OUT_PATH}/10X_hifi/ref"
cp /home/users/ddu/software/Genome10X/longranger-2.2.2/refdata-GRCh38-2.1.0/fasta/genome.fa "${OUT_PATH}/10X_hifi/ref/hg38"
samtools faidx "${OUT_PATH}/10X_hifi/ref/hg38"

docker run -v "${OUT_PATH}/10X_hifi:/input" -v "${OUT_PATH}/10X_hifi/CCS_variants:/output" \
    google/deepvariant:0.10.0 /opt/deepvariant/bin/run_deepvariant \
    --model_type=PACBIO \
    --ref=/input/ref/hg38 \
    --reads=/input/${CELL}.hg38.bam \
    --output_vcf=/output/vcf/pacbioccs.chr6.vcf.gz \
    --regions "chr6" \
    --num_shards=11

echo "DeepVariant processing completed at $(date)"

# Run Whatshap genotype analysis
whatshap genotype --ignore-read-groups --chromosome 6 \
    --reference /home/users/xyxu/assembly/rawdata/10X_hifi/ref/genome.fa \
    -o "${OUT_PATH}/10X_hifi/CCS_variants/vcf/pacbioccs.chr6.WG.vcf.gz" \
    "${OUT_PATH}/10X_hifi/CCS_variants/vcf/pacbioccs.chr6.vcf.gz" \
    "${OUT_PATH}/10X_hifi/${CELL}.hg38.bam"
echo "Whatshap genotype analysis completed at $(date)"
