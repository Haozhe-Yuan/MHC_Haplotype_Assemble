#!/bin/bash

# Set up the environment for tools

# Cell and input paths
CELL=
FINAL1=~/WGS-${CELL}.Switched/00.ref/${CELL}.Switched.hap1.fa
FINAL2=~/WGS-${CELL}.Switched/00.ref/${CELL}.Switched.hap2.fa
REF=~/ref/mhc/mhc-hg38.fa
WKDIR=~/results/quast

# Check file size
du -h $FINAL1 $FINAL2

# Create working directory and run QUAST
mkdir -p $WKDIR/$CELL && cd $WKDIR && pwd
quast -t 40 -r $REF -o ./$CELL -l "Hap1,Hap2" $FINAL1 $FINAL2

# Busco analysis for different cells
function run_busco {
    CELL=$1
    INPUT1=$2
    INPUT2=$3
    DATASET=~/ref/busco/primates_odb10

    for i in 1 2; do
        INPUT_VAR="INPUT$i"
        OPTH="~/results/busco/${CELL}-${i}"
        mkdir -p $OPTH && cd $OPTH || exit
        nohup busco -i ${!INPUT_VAR} -c 60 -o ${CELL}-${i} -f -m geno -l $DATASET --offline &
    done
}


# Prepare and generate plots
python3 generate_plot.py -wd .

# Meryl and Merqury setup
cd ~/tools
tar -xJf meryl-1.4.*.tar.xz
cd meryl-1.4.1/bin || exit
export PATH=$PWD:$PATH

cd ~/tools && git clone https://github.com/marbl/merqury.git
cd merqury || exit
export MERQURY=$PWD

# Generate meryl files for fastq
function generate_meryl {
    CELL=$1
    BAM=~/workspace/cellsHLA_refine/processed_data/${CELL}/hifi/q30NA/10X_hifi/bam/pacbioccs.mhc.new_splitWP.bam
    samtools bam2fq $BAM > ~/merqury_fqgz/${CELL}.fastq
    gzip -c ~/merqury_fqgz/${CELL}.fastq > ~/merqury_fqgz/${CELL}.fastq.gz
}

for CELL in A549 Hela HepG2 K562 ${CELL}; do
    generate_meryl $CELL
done

# Merqury analysis
mkdir -p ~/merqury_results && cd ~/merqury_results && pwd
$MERQURY/best_k.sh 4548949

# Run Hapmer and Merqury
export PATH="~/miniconda3/envs/chip10/bin/bedtools:$PATH"
export PATH="~/miniconda3/envs/chip10/bin/samtools:$PATH"
export PATH="~/tools/meryl-1.4.1/bin:$PATH"
export MERQURY=~/tools/merqury

for CELL in A549 Hela HepG2 K562 U2OS; do
    samtools bam2fq ~/merqury_fqgz/${CELL}.fastq.gz > ./${CELL}_HP1.fastq
    samtools bam2fq ~/merqury_fqgz/${CELL}.fastq.gz > ./${CELL}_HP2.fastq
    ~/tools/meryl-1.4.1/bin/meryl k=16.0407 count output ${CELL}.meryl ./${CELL}.fastq.gz
done

$MERQURY/trio/hapmers.sh ${CELL}_HP1.meryl ${CELL}_HP2.meryl ${CELL}.meryl
$MERQURY/merqury.sh ${CELL}.meryl ${CELL}_HP1.meryl ${CELL}_HP2.meryl $READ1 $READ2 ${CELL}

# HLA Typing analysis
CELL=
QU=30
OUT_PATH=/home/users/hzyuan/workspace/cellsHLA_refine/processed_data/${CELL}/hifi/q30NA

mkdir -p ~/hla-typing-raw && cd ~/hla-typing-raw && pwd
samtools view -@ 60 -bh ${OUT_PATH}/10X_hifi/${CELL}.hg38.bam chr6:28903952-33268517 | samtools bam2fq - | seqkit fq2fa -w 0 > ./${CELL}.hg38.mhc.bam.reads

hifihla call-consensus --fasta ./${CELL}.hg38.mhc.bam.reads --threads 1 --outdir ./${CELL} -vv
samtools bam2fq -@ 60 ${OUT_PATH}/00_Remove_Dup/${CELL}_Q${QU}.hifi_reads.bam | less -S

# Immuannot analysis
cd ~/tools/Immuannot && pwd
for CELL in A549 Hela HepG2 K562 U2OS; do
    read1=~/WGS-Switched/00.ref/${CELL}.Switched.hap1.fa
    read2=~/WGS-Switched/00.ref/${CELL}.Switched.hap2.fa
    for CTG in $read1 $read2; do
        script=scripts/immuannot.sh
        refdir=refData-2023Jun05
        outpref=~/hla-typing-raw/${CELL}-HP1
        bash ${script} -t 40 -c $CTG -r $refdir -o $outpref
    done
done

# =======================Variants calling
minimap2=~/minimap2/minimap2

# Define cell types and reference paths
CELL=U20S
final1=~/WGS-U20S.Switched/00.ref/U20S.Switched.hap1.finally.fa
final2=~/WGS-U20S.Switched/00.ref/U20S.Switched.hap2.finally.fa

hap1fa=$final1
hap2fa=$final2

hg38=
# mhchg38=~/ref/mhc/mhc-hg38.fa

# Create working directory
mkdir -p results/variants/$CELL && wkdir=results/variants/$CELL
cd $wkdir && pwd

# Minimap2 alignment for hap1 and hap2
echo "Running minimap2 and sorting BAM files..."
$minimap2 --paf-no-hit -a -x asm20 -z200000,10000 --cs -r2k -t 60 $hg38 $hap1fa | samtools sort -@ 60 --output-fmt BAM -o hifiasm_WP_10X_H1_new3.50kb.aligned.bam
samtools index -@ 60 hifiasm_WP_10X_H1_new3.50kb.aligned.bam

$minimap2 --paf-no-hit -a -x asm20 -z200000,10000 --cs -r2k -t 60 $hg38 $hap2fa | samtools sort -@ 60 --output-fmt BAM -o hifiasm_WP_10X_H2_new3.50kb.aligned.bam
samtools index -@ 60 hifiasm_WP_10X_H2_new3.50kb.aligned.bam
echo "========Minimap2 and Samtools Done at $(date)"

# =============== SV Calling =======================
echo "Running SV calling using svim-asm..."
svim-asm diploid ./ hifiasm_WP_10X_H1_new3.50kb.aligned.bam hifiasm_WP_10X_H2_new3.50kb.aligned.bam $hg38

# Extract homozygous and heterozygous variants from VCF
echo "Extracting SVs from VCF..."
perl -lane 'print if /1\/1/' variants.vcf | perl -lane 'print "$F[0]\t$F[1]\t$F[-3]\t$F[-1]"' > home_sv.tsv
perl -lane 'print if /1\/0/' variants.vcf | perl -lane 'print "$F[0]\t$F[1]\t$F[-3]\t$F[-1]"' > hap1_sv.tsv
perl -lane 'print if /0\/1/' variants.vcf | perl -lane 'print "$F[0]\t$F[1]\t$F[-3]\t$F[-1]"' > hap2_sv.tsv

# Summarize SV types
echo "Summarizing SV types..."
{
  echo "home"
  awk '{print $3}' home_sv.tsv | cut -d ';' -f1 | sort | uniq -c
  echo "hap1"
  awk '{print $3}' hap1_sv.tsv | cut -d ';' -f1 | sort | uniq -c
  echo "hap2"
  awk '{print $3}' hap2_sv.tsv | cut -d ';' -f1 | sort | uniq -c
} > sv.txt

# Concatenate SV files for plotting
cat hap1_sv.tsv home_sv.tsv > hap1-home.tsv
cat hap2_sv.tsv home_sv.tsv > hap2-home.tsv
wc -l *.tsv

software='dipcall'

# Generate Makefile for diploid variant calling
echo "Running dipcall variant calling..."
run-dipcall -t 60 prefix.$software $hg38 $hap1fa $hap2fa > prefix.$software'.mak'

# Modify Makefile to set correct parameters
sed -i 's/\-xasm5/ \-z200000,10000/g' prefix.$software'.mak'
sed -i 's/samflt/samflt \-L 10000/g' prefix.$software'.mak'

# Run Makefile to generate diploid variants
make -j2 -f prefix.$software'.mak'

# Filter and replace syndip with cell name, then extract specific region from VCF
echo "Filtering and extracting VCF for specific region..."
zcat prefix.$software'.dip.vcf.gz' | perl -le 'while(<>){chomp;if(/^#/){print;}else{@F=split/\s+/,$_;next if $F[-1]=~/\./;print if ($F[0] eq "chr6") && ($F[1]>=28903952) && ($F[1]<=33268517);};}' \
  | sed 's/syndip/${CELL}/g' > prefix.$CELL.$software'.dip.vcf'

# Generate BCFTools stats
echo "Generating BCFTools stats..."
bcftools stats prefix.$CELL.$software'.dip.vcf' > dipcall.bcftools.stats

# Count the number of different genotypes in the VCF
echo "Counting genotype occurrences..."
awk '{
    genotype[$10]++
}
END {
    for (g in genotype) {
        print g, genotype[g]
    }
}' prefix.${CELL}.dipcall.dip.vcf

# Extract heterozygous variants for hap1, hap2, and common homozygous variants
echo "Extracting heterozygous and homozygous variants..."
awk 'BEGIN {
    hap1_het = 0; hap2_het = 0; common_hom = 0;
}
{
    if ($10 ~ /^1\|0/ || $10 ~ /^1\|2/ || $10 ~ /^2\|1/) {
        hap1_het++;
        print > "hap1.var.tsv";
    }
    if ($10 ~ /^0\|1/ || $10 ~ /^1\|2/ || $10 ~ /^2\|1/) {
        hap2_het++;
        print > "hap2.var.tsv";
    }
    if ($10 ~ /^1\|1/ || $10 ~ /^2\|2/) {
        common_hom++;
        print > "home.var.tsv";
    }
}
END {
    print "hap1 heterozygous: " hap1_het > "hap1.var.tsv";
    print "hap2 heterozygous: " hap2_het > "hap2.var.tsv";
    print "common homozygous: " common_hom > "home.var.tsv";
}' prefix.*.dipcall.dip.vcf

# Count variants by file
wc -l hap1.var.tsv hap2.var.tsv home.var.tsv

# Count homozygous variants
awk '{ if ($10 ~ /^1\|1/ || $10 ~ /^2\|2/ || $10 ~ /^0\|0/) count++; } END { print count; }' prefix.*.dipcall.dip.vcf > home.var.count

# Extract all variants for plotting
echo "Extracting all variants for plotting..."
awk '{
    if ($10 ~ /^1\|0/ || $10 ~ /^1\|2/ || $10 ~ /^2\|1/ || $10 ~ /^1\|1/)
        print > "plot.hap1.all.var.tsv";
    if ($10 ~ /^0\|1/ || $10 ~ /^1\|2/ || $10 ~ /^2\|1/ || $10 ~ /^1\|1/)
        print > "plot.hap2.all.var.tsv";
}' prefix.${CELL}.dipcall.dip.vcf

# Count lines in all TSV files
wc -l *.tsv

# ========== Output for Plotting =================
# hap1-home.tsv
# hap2-home.tsv
# plot.hap1.all.var.tsv
# plot.hap2.all.var.tsv

# ========== Other Outputs =======================
# variants.vcf
# sv.txt
# prefix.$CELL.$software'.dip.vcf'