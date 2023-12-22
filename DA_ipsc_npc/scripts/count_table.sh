#!/bin/bash

dir_bam=/data/share/htp/ATACseq/Novogene_ATAC/02_Preprocessing/bwa_mapping/
base=/data/share/htp/pleiotropy/paper_data/DA_ipsc_npc/
hum_bed=$base/bed/hg38_DHS.bed
mac_bed=$base/bed/macFas6_DHS.bed
dir_out=$base/count_tables



#declare -a hg=("01" "02" "05" "06" "09" "10")
#declare -a mc=("03" "04" "07" "08" "11" "12")
#for all samples
#sbatch --output=$base/slurm-%J.out --wrap="bedtools multicov -bams $dir_bam/ATAC_b1_sample*_hg38.bam -bed $hum_bed > $dir_out/hg38_counts.bed"
#sbatch --output=$base/slurm-%J.out --wrap="bedtools multicov -bams $dir_bam/ATAC_b1_sample*_macFas6.bam -bed $mac_bed > $dir_out/macFas6_counts.bed"


#for self-species samples
sbatch --output=$base/slurm-%J.out --wrap="bedtools multicov \
-bams $dir_bam/ATAC_b1_sample01_hg38.bam $dir_bam/ATAC_b1_sample02_hg38.bam \
$dir_bam/ATAC_b1_sample05_hg38.bam $dir_bam/ATAC_b1_sample06_hg38.bam \
-bed $hum_bed > $dir_out/hg38_counts.bed"

# $dir_bam/ATAC_b1_sample04_macFas6.bam  sample 4 somehow missing from macfas mapping?!?!?!
sbatch --output=$base/slurm-%J.out --wrap="bedtools multicov \
-bams $dir_bam/ATAC_b1_sample03_macFas6.bam $dir_bam/ATAC_b1_sample04_macFas6.bam \
$dir_bam/ATAC_b1_sample07_macFas6.bam $dir_bam/ATAC_b1_sample08_macFas6.bam \
-bed $mac_bed > $dir_out/macFas6_counts.bed"  



