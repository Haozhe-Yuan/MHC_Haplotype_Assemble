#!/bin/bash

# Set working directory
wkdir=~/WGS-U20S.Switched
mkdir -p $wkdir && cd $wkdir && pwd

# Generate reference sequences
mkdir -p $wkdir/00.ref && cd $wkdir/00.ref || exit

final1=
final2=
# Check if the files exist and have content
grep '^>' $final1
grep '^>' $final2

# Produce and index hap-specific reference
hg38_nonMHC=~/ref/mixhg38/hg38.nonMHC.fa

# Generate bed files for fasta sequences
samtools faidx $final1
samtools faidx $final2

# Create bed files
awk '{print $1, 1, $2}' $final1.fai | sed 's/ /\t/g' > $final1.bed
awk '{print $1, 1, $2}' $final2.fai | sed 's/ /\t/g' > $final2.bed

# Sort bed files
bedtools sort -chrThenSizeA -i $final1.bed > $final1.pos
bedtools sort -chrThenSizeA -i $final2.bed > $final2.pos

# Combine with hg38 reference
hap1mixhg38=hg38.${CELL}.Switched.HP1.mix.fa
hap2mixhg38=hg38.${CELL}.Switched.HP2.mix.fa
cat $hg38_nonMHC $final1 > $hap1mixhg38
cat $hg38_nonMHC $final2 > $hap2mixhg38

# Generate bwa index for both haplotypes
nohup bwa index $hap1mixhg38 &
nohup bwa index $hap2mixhg38 &

# Setup alignment
cd $wkdir || exit
mkdir -p 02.align

# Set data path for fastq files
CELL=U20S
output_path=~/WGS-U20S.Switched/02.align

# Update paths to trimmed fastq files
fq1=~/WGS-cells/${CELL}/01.trim/fq1_1_val_1.fq.gz
fq2=~/WGS-cells/${CELL}/01.trim/fq1_2_val_2.fq.gz
fq3=~/WGS-cells/${CELL}/01.trim/fq2_1_val_1.fq.gz
fq4=~/WGS-cells/${CELL}/01.trim/fq2_2_val_2.fq.gz

# Reference paths
ref1=hg38.Switched.HP1.mix.fa
ref2=hg38.Switched.HP2.mix.fa

# Align reads to hap1 and hap2 references
sbatch -p fat align.sh $output_path $fq1 $fq2 $fq3 $fq4 $ref1 $ref2
sbatch -p fat align.sh $output_path $fq1 $fq2 $fq3 $fq4 $ref1 $ref2

# ===================== align.sh =====================

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --job-name=U20S-hap1-align

# Initialize Conda environment

# Set variables
CELL=U20S
output_path=$1
fq1=$2
fq2=$3
fq3=$4
fq4=$5
ref1=$6
ref2=$7

cd ${output_path} || exit

# Align reads for hap1
bwa mem -t 72 $ref1 $fq1 $fq2 | samtools view -Sb - > ${CELL}.hap1.bam
bwa mem -t 72 $ref1 $fq3 $fq4 | samtools view -Sb - > ${CELL}.hap1.l2.bam

# Sort and index bam files
samtools sort -@ 72 -O BAM ${CELL}.hap1.bam > ${CELL}.hap1.sort.bam && samtools index -@ 72 ${CELL}.hap1.sort.bam
samtools sort -@ 72 -O BAM ${CELL}.hap1.l2.bam > ${CELL}.hap1.l2.sort.bam && samtools index -@ 72 ${CELL}.hap1.l2.sort.bam

# Use GATK to mark duplicates
conda activate gatk
gatk MarkDuplicates \
  -I ${CELL}.hap1.sort.bam \
  -I ${CELL}.hap1.l2.sort.bam \
  -O ${CELL}.hap1.sort.markdup.bam \
  -M ${CELL}.hap1.sort.markdup_metrics.txt
conda deactivate 

# Generate perfect version
input=${CELL}.hap1.sort.markdup.bam
samtools view -h ${input} | perl -ne 'print if /^@/ || /NM:i:0/' | samtools sort -@ 60 -O BAM -o ${CELL}.hap1.perfect.sort.bam && samtools index -@ 60 ${CELL}.hap1.perfect.sort.bam

# Generate imperfect version
samtools sort -@ 60 -O BAM -o ${CELL}.hap1.sort.markdup.sort.bam ${CELL}.hap1.sort.markdup.bam && samtools index -@ 70 ${CELL}.hap1.sort.markdup.sort.bam

# Repeat for hap2 if necessary 

# ==================== Coordinate transform
# ---- Helper Functions ----

# Process dnadiff results (same steps for hap1 and hap2)
process_dnadiff() {
    local INFOR_FILE=$1
    local FINAL_FILE=$2
    local HAP=$3
    
    for i in $(grep '>' $FINAL_FILE | sed 's/>//g'); do
        show-aligns ${SAMPLE}.${HAP}.filter.delta chr6:28903952-33268517 $i | \
        perl -le 'while(<>){chomp;next if /^(\s+)|^$/;next if $.<=3;@F=split/\s+/,$_;if(/-/){print}else{print ;};}' | \
        perl -le '$n=$m=0;while(<>){chomp;chomp($G=<>);@A=@B=@AA=@BB=();if (/^--/ or /^=/){if(/Alignments between (\S+) and (\S+)/){$ref=$1;$denovo=$2;};if($G=~/BEGIN alignment \[ (.*)\]/){$M=$1;};@MM=split/\s+/,$M;print $_;print $G;}else{@A=split/\s+/,$_;@AA=split//,$A[1];@B=split/\s+/,$G;@BB=split//,$B[1];};$a=$A[0];$b=$B[0];for my $i (0..$#AA){if ($AA[$i] ne "."){if($MM[0] eq "+1"){$aa=$a+$i-$n}elsif($MM[0] eq "-1"){$aa=$a-$i+$n};}else{$aa ="N";$n+=1;};if ($BB[$i] ne "."){if($MM[5] eq "+1"){$bb=$b+$i-$m}elsif($MM[5] eq "-1"){$bb=$b-$i+$m};}else{$bb ="N";$m+=1;};if($aa eq "N"){$chr6 = "N"}else{$chr6=$aa+28903952-1};print "chr6\t$ref\t$denovo\t$chr6\t$aa\t$bb\t$AA[$i]\t$BB[$i]\t$MM[0]\t$MM[5]\t\t---$_\t$G";};$n=$m=0; ;}' \
        >> $INFOR_FILE
    done
}

# Process results after dnadiff (sorting, filtering, etc.)
process_results() {
    local INFOR_FILE=$1
    
    perl -le 'print "hg38\tmhc\tdenovo\thg38POS\tmhcPOS\tdenovoPOS";while(<>){chomp;@F=split/\s+/,$_;next if /^--|^=/;if($F[9] eq "+1"){print join("\t",@F[0..5]);}else{$m=$F[5]-1;print join("\t",@F[0..4]),"\t$m";};}' $INFOR_FILE | \
    perl -le 'while(<>){chomp;@F=split/\s+/,$_;push @{$hash{$F[3]}},$_;}END{foreach my $key (keys %hash){print "$key\t",join("\t",@{$hash{$key}});};}' | \
    sort -k1 | \
    perl -lane 'next if /hg38|N/;print join("\t",@F[1..6]) if $#F<=6' | \
    perl -le 'while(<>){chomp;@F=split/\s+/,$_;push @{$hash{"$F[2]\t$F[5]"}},$_;}END{foreach my $key (keys %hash){print "$key\t",join("\t",@{$hash{$key}});};}' | \
    perl -lane 'next if /hg38|N/;print join("\t",@F[2..7]) if $#F<=7' | \
    sort -k4 > ${INFOR_FILE}.bed.unique.pos.new.bed
}

# Set paths to commonly used files and directories
WORKDIR=
PYTHON_ENV=
MHC_REF="~/ref/mhc/mhc-hg38.fa"
SAMPLE="K562"  # You can change this as needed
FINAL1="path_to_final1.bam"  # Set this appropriately
FINAL2="path_to_final2.bam"  # Set this appropriately

# 1. Disk usage and Samtools depth
du -h ${SAMPLE}.common.perfect.common.bam ${SAMPLE}.unique.perfect.hap1.bam ${SAMPLE}.unique.perfect.hap2.bam
samtools depth -a ${SAMPLE}.unique.perfect.hap1.bam > ${SAMPLE}.unique.perfect.hap1.depth.txt
samtools depth -a ${SAMPLE}.unique.perfect.hap2.bam > ${SAMPLE}.unique.perfect.hap2.depth.txt
du -h ./*depth.txt

# 2. Create directories and dnadiff
mkdir -p HP1_dnadiff && cd HP1_dnadiff && pwd
dnadiff $MHC_REF $FINAL1
INFOR1="${SAMPLE}.hap1.filter.delta.coords.align.detail.infor"

process_dnadiff $INFOR1 $FINAL1 "hap1"  # Create a function for this part to avoid redundancy

# 3. Process and Sort Results
process_results $INFOR1

# 4. Repeat for HP2
mkdir -p HP2_dnadiff && cd HP2_dnadiff && pwd
dnadiff $MHC_REF $FINAL2
INFOR2="${SAMPLE}.hap2.filter.delta.coords.align.detail.infor"
process_dnadiff $INFOR2 $FINAL2 "hap2"  # Use the same function for hap2

# 5. Process and Sort Results for HP2
process_results $INFOR2
