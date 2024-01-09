#!/bin/bash
#SBATCH --error=cs.%J.err
#SBATCH --output=cs.%J.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G

target=/data/share/htp/EBgrant/genome_data/macFas6/fasta/genome.fa
reference=/data/home/EBgrant/genome_data_new/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa
gtf=/data/share/htp/hack_GRN/NPC_diff_network_analysis/remapping/genomes/hg38/genes.gtf

conda activate liftoff

liftoff -g $gtf -o mf6_liftoff_polished.gtf -u mf6_notfound.gtf  -polish -cds $target $reference 


