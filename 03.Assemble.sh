#!/bin/bash

set -e

# Define paths
HIFI_30_ROOT="/home/users/Liulab_data/xyxu/hifi/Q30/X101SC23065068-Z01-F001/Data-X101SC23065068-Z01-F001/${CELL}"
ADD_HIFI_30_ROOT="/home/users/Liulab_data/xyxu/addhifi_Clinical_Q30/X101SC23065068-Z01-F003/Data-X101SC23065068-Z01-F003/${CELL}"
OUT_DATA_ROOT="/home/users/hzyuan/workspace/cellsHLA_refine/processed_data/${CELL}"

# Create directories
mkdir -p "${OUT_DATA_ROOT}/hifi/q30A"
OUT_30_DATA="${OUT_DATA_ROOT}/hifi/q30A"

# Set paths
RAW_PATH=${HIFI_30_ROOT}
RAW_ADD_PATH=${ADD_HIFI_30_ROOT}
OUT_PATH=${OUT_30_DATA}
SUFFIX=filt.hifi.bam

# Extract phased variants from VCF files
perl -le 'open(IN,$ARGV[0]);open(IN1,$ARGV[1]);open(IN2,$ARGV[2]);
while(<IN>){chomp;next if /^#/ or /1\/1/ or /0\/0/ or /\.\/\./;@F=split; $hash{$F[1]}=1;}
while(<IN1>){chomp;next if /^#/ or /1\/1/ or /0\/0/ or /\.\/\./;@F=split; $hash2{$F[1]}=1;}
while(<IN2>){chomp;next if /1\|1/ or /0\|0/;print if /^#/;@F=split;@M=split/:/,$F[-1];next if $M[0]=~/\./;if($hash{$F[1]} && $hash2{$F[1]}){print;}}' \
    <(zcat ./10X_hifi/CCS_variants/vcf/pacbioccs.chr6.vcf.gz) \
    <(zcat ./10X_hifi/CCS_variants/vcf/pacbioccs.chr6.WG.vcf.gz) \
    ${OUT_DATA_ROOT}/10X/L30ul_7/L30ul_7_all.vcf.gz.mhc.new.vcf \
    > ./10X_hifi/CCS_variants/vcf/phased_variants.vcf.mhc.ccs.10X.WG.heter.vcf

# Process VCF
sed -i "s/$(perl -lane 'print $F[9] if /^#CHROM/' ./10X_hifi/CCS_variants/vcf/phased_variants.vcf.mhc.ccs.10X.WG.heter.vcf)/${CELL}/g" \
    ./10X_hifi/CCS_variants/vcf/phased_variants.vcf.mhc.ccs.10X.WG.heter.vcf
bgzip -f ./10X_hifi/CCS_variants/vcf/phased_variants.vcf.mhc.ccs.10X.WG.heter.vcf
tabix -p vcf -f ./10X_hifi/CCS_variants/vcf/phased_variants.vcf.mhc.ccs.10X.WG.heter.vcf.gz

# Run Whatshap haplotag
whatshap haplotag --reference /home/users/ddu/software/Genome10X/longranger-2.2.2/refdata-GRCh38-2.1.0/fasta/genome.fa \
    -o ./10X_hifi/bam/pacbioccs.mhc.new_splitWP.bam \
    ./10X_hifi/CCS_variants/vcf/phased_variants.vcf.mhc.ccs.10X.WG.heter.vcf.gz \
    ./10X_hifi/${CELL}.hg38.mhc.bam
echo "===== MarkHap done at $(date)"

# Extract and filter reads
samtools view -h $OUT_PATH/bam/pacbioccs.mhc.new_splitWP.bam | grep -e "HP:i:1" -e "^@" | samtools bam2fq - | seqkit fq2fa -w 0 > $OUT_PATH/bam/pacbioccs.mhc.new_splitWP.bam.HP1.new.reads
samtools view -h $OUT_PATH/bam/pacbioccs.mhc.new_splitWP.bam | grep -e "HP:i:2" -e "^@" | samtools bam2fq - | seqkit fq2fa -w 0 > $OUT_PATH/bam/pacbioccs.mhc.new_splitWP.bam.HP2.new.reads
samtools view -h $OUT_PATH/bam/pacbioccs.mhc.new_splitWP.bam | grep -v -e "HP:i:1" -e "HP:i:2" | samtools bam2fq - | seqkit fq2fa -w 0 > $OUT_PATH/bam/pacbioccs.mhc.new_splitWP.bam.untagged.new.reads

# Run Hifiasm for phased reads
for hp in HP1 HP2; do
    hifiasm -o $OUT_PATH/asm_Q${QU}/hifiasm_${hp} -t 40 \
        $OUT_PATH/bam/pacbioccs.mhc.new_splitWP.bam.${hp}.new.reads \
        $OUT_PATH/bam/pacbioccs.mhc.new_splitWP.bam.untagged.new.reads
gfatools gfa2fa $OUT_PATH/asm_Q${QU}/hifiasm_${hp}.bp.p_ctg.gfa > $OUT_PATH/asm_Q${QU}/hifiasm_${hp}.bp.p_ctg.fa
seqkit seq -w 0 -m 50000 $OUT_PATH/asm_Q${QU}/hifiasm_${hp}.bp.p_ctg.fa > $OUT_PATH/asm_Q${QU}/hifiasm_${hp}.bp.p_ctg.50kb.fa
done
echo "======= Hifiasm and 50kb extraction done at $(date)"

# Run Hifiasm for phased reads
for hp in HP1 HP2; do
    hifiasm -o $OUT_PATH/asm_Q${QU}/hifiasm_${hp}_Cov30 -t 40 \
        $OUT_PATH/bam/pacbioccs.mhc.new_splitWP.${hp}.Cov30.bam.new.reads \
        $OUT_PATH/bam/pacbioccs.mhc.new_splitWP.bam.untagged.new.reads
    gfatools gfa2fa $OUT_PATH/asm_Q${QU}/hifiasm_${hp}_Cov30.bp.p_ctg.gfa > $OUT_PATH/asm_Q${QU}/hifiasm_${hp}_Cov30.bp.p_ctg.fa
    seqkit seq -w 0 -m 50000 $OUT_PATH/asm_Q${QU}/hifiasm_${hp}_Cov30.bp.p_ctg.fa > $OUT_PATH/asm_Q${QU}/hifiasm_${hp}_Cov30.bp.p_ctg.50kb.fa
done
echo "======= Hifiasm and extract 50kb from 30X done at $(date)"

# Merge haplotypes into a single FASTA
for hp in HP1 HP2; do
    perl -lane 'if(/^>/){print $_,"_${hp}"}else{print;}' $OUT_PATH/asm_Q${QU}/hifiasm_${hp}_Cov30.bp.p_ctg.50kb.fa >> $OUT_PATH/asm_Q${QU}/hifiasm_Cov30.bp.p_ctg.50kb.fa
done

# Align reads using minimap2
minimap2 -t 40 -R '@RG\tID:1\tSM:${CELL}' -ax asm20 \
    $OUT_PATH/asm_Q${QU}/hifiasm_Cov30.bp.p_ctg.50kb.fa \
    <(samtools bam2fq -@30 ./00_Remove_Dup/${CELL}_Q${QU}.hifi_reads.bam) | \
    samtools sort -@40 --output-fmt BAM -o $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam
samtools index -@40 $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam
echo "========= Minimap2 and samtools Done at $(date)"

# Filter supplementary alignments
samtools view -q 60 -f 0x800 $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam \
    | perl -lane 'if(/SA:Z:(.*),?-?,/){$b=$1}; print "$F[0]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$b"' \
    | perl -lane '@M=split/;/,$F[-1]; $i=0; @H=split/_/,$F[1]; for(@M){if(/$H[-1]/){$i+=1;}}; if($i==($#M+1)){print $_}' \
    > $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam.supplement60.txt

samtools view -q 60 -F 0x800 $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam > ./pacbioccs.mhc_asm.new.filtered.bam
echo "========= Filtering completed at $(date)"

# Filter supplementary alignments with specific match percent and score
perl -MList::Util -lane 'use List::Util qw/sum/;@MM=split/;/,$F[-1];$i=0;$j=0;
(@M)=($F[4]=~/\d+M/g);(@I)=($F[4]=~/\d+I/g);(@S)=($F[4]=~/\d+S/g);
(@eq)=($F[4]=~/\d+=/g);(@X)=($F[4]=~/\d+X/g);
$total=sum(@M)+sum(@I)+sum(@S)+sum(@eq)+sum(@X);
$total_m=sum(@M);$rate=$total_m/$total;
if(($rate>=0.5) && ($F[3]>=60)){$i+=1;};print if $i>=1;'
    $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam.supplement60.txt |
perl -MList::Util -lane '@MM=split/;/,$F[-1];$i=0;$j=0;
for my $m (@MM){chomp;use List::Util qw/sum/;@NN=split/,/,$m;
(@M)=($NN[3]=~/\d+M/g);(@I)=($NN[3]=~/\d+I/g);(@S)=($NN[3]=~/\d+S/g);
(@eq)=($NN[3]=~/\d+=/g);(@X)=($NN[3]=~/\d+X/g);
$total=sum(@M)+sum(@I)+sum(@S)+sum(@eq)+sum(@X);
$total_m=sum(@M);$rate=$total_m/$total;
if(($rate>=0.5) && ($NN[-1]>=60)){$j+=1;};};print if $j>=1;'
    > $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam.supplement60_0.5match_60score.txt

# Extract read IDs for haplotype 1 and 2
perl -lane 'print $F[0] if $F[1]=~/hp1/' \
    $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam.supplement60_0.5match_60score.txt \
    | sort -u > $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam.supplement60_0.5match_60score.txt.hp1.read.id

perl -lane 'print $F[0] if $F[1]=~/hp2/' \
    $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam.supplement60_0.5match_60score.txt \
    | sort -u > $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam.supplement60_0.5match_60score.txt.hp2.read.id

echo "===== Filtering and haplotype read extraction completed at $(date)"

# Extract filtered reads
seqkit grep --threads 40 -f \
    $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam.supplement60_0.5match_60score.txt.hp1.read.id \
    <(samtools bam2fq -@40 ./00_Remove_Dup/${CELL}_Q${QU}.hifi_reads.bam | seqkit fq2fa -w 0) \
    > $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam.supplement60_0.5match_60score.txt.hp1.read.id.reads

seqkit grep --threads 40 -f \
    $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam.supplement60_0.5match_60score.txt.hp2.read.id \
    <(samtools bam2fq -@40 ./00_Remove_Dup/${CELL}_Q${QU}.hifi_reads.bam | seqkit fq2fa -w 0) \
    > $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam.supplement60_0.5match_60score.txt.hp2.read.id.reads

echo "========Sort and seqkit Done at $(date)"

# Run Hifiasm with filtered reads
for hp in HP1 HP2; do
    hifiasm -o $OUT_PATH/asm_Q${QU}/hifiasm_${hp}_Cov30_S -t 40 \
        $OUT_PATH/bam/pacbioccs.mhc.new_splitWP.${hp}.Cov30.bam.new.reads \
        $OUT_PATH/bam/pacbioccs.mhc.new_splitWP.bam.untagged.new.reads \
        $OUT_PATH/bam/pacbioccs.mhc_asm_Cov30.new.bam.supplement60_0.5match_60score.txt.${hp}.read.id.reads

    gfatools gfa2fa $OUT_PATH/asm_Q${QU}/hifiasm_${hp}_Cov30_S.bp.p_ctg.gfa > $OUT_PATH/asm_Q${QU}/hifiasm_${hp}_Cov30_S.bp.p_ctg.fa
    seqkit seq -w 0 -m 50000 $OUT_PATH/asm_Q${QU}/hifiasm_${hp}_Cov30_S.bp.p_ctg.fa > $OUT_PATH/asm_Q${QU}/hifiasm_${hp}_Cov30_S.bp.p_ctg.50kb.fa
done
echo "======= Hifiasm and 50kb extraction completed at $(date)"

# Create final assembly directory
mkdir -p $OUT_PATH/finalfa

# Normalize Contig orientations based on HG38 reference
for hp in HP1 HP2; do
    awk '$NF == "True" { print }' $OUT_PATH/QUAST_HG38/contigs_reports/all_alignments_${hp}.tsv > $OUT_PATH/finalfa/${hp}_1.txt
    awk '{print $0, $4-$3}' $OUT_PATH/finalfa/${hp}_1.txt > $OUT_PATH/finalfa/${hp}_2.txt
    awk '{print $6, $NF}' $OUT_PATH/finalfa/${hp}_2.txt > $OUT_PATH/finalfa/${hp}_3.txt
    
    awk '
    {
        abs_value = ($2 >= 0) ? $2 : -$2;
        if (!($1 in max_values) || abs_value > max_values[$1]) {
            max_values[$1] = abs_value;
            lines[$1] = $0;
        }
    }
    END {
        for (key in lines) {
            print lines[key];
        }
    }' $OUT_PATH/finalfa/${hp}_3.txt > $OUT_PATH/finalfa/${hp}_4.txt
    
    awk '{if ($2 < 0) {a=$1; system("samtools faidx ./10X_hifi/asm_Q**/hifiasm_${hp}_Cov30_S.bp.p_ctg.50kb.fa " a " |seqkit seq -r -p --seq-type DNA >> ./10X_hifi/finalfa/hifiasm_${hp}_Cov30_S.bp.p_ctg.50kb.finally.fa")} else if ($2 > 0) {b=$1; system("samtools faidx ./10X_hifi/asm_Q**/hifiasm_${hp}_Cov30_S.bp.p_ctg.50kb.fa " b " >> ./10X_hifi/finalfa/hifiasm_${hp}_Cov30_S.bp.p_ctg.50kb.finally.fa")}}' $OUT_PATH/finalfa/${hp}_4.txt
    
    echo "===============Cov30_S.bp.p_ctg.50kb.finally.fa of ${hp} Generated"
done
