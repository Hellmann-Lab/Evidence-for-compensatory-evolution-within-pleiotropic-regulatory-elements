#!/bin/bash

workdir=/data/share/htp/pleiotropy/paper_data/INSIGHT
slurmdir=$workdir/slurm
mkdir -p $slurmdir
slurmpfx="$slurmdir/slurm.%J"
out=all_per_tissue


## declare an array variable
declare -a arr=("adrenal_gland" "brain" "heart" "kidney" "large_intestine" "lung" "muscle" "stomach" "thymus") 


## now loop through the above array
for i in "${arr[@]}"
do
mkdir -p $workdir/$out/$i/emOutput


for n in {1..9}
do

sbatch -J INS"$n"_"$i" --error=$slurmpfx$i$n.err --output=$slurmpfx$i$n.out --workdir=$workdir \
--wrap="$workdir/scripts/INSIGHT_filtering.sh \
$workdir/full_original/emInput/sp$n.ins \
$workdir/$out/$i/sp"$n".bed \
$workdir/$out/$i/sp"$n"_tmp.ins \
$workdir/$out/$i/sp"$n".ins \
$workdir/full_original/emInput/sp$n.flankPoly.forBetas.ins \
sp"$n"_"$i" \
$workdir/$out/$i/emOutput/ \
$workdir/$out/$i/sp$n.flanking.filtered.ins"

   
done
done
