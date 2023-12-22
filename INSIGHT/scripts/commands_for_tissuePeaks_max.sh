#!/bin/bash

workdir=/data/share/htp/pleiotropy/paper_data/INSIGHT
slurmdir=$workdir/slurm
mkdir -p $slurmdir
slurmpfx="$slurmdir/slurm.%J"
out=max_tissuePeak/all_tissues2
in=$workdir/max_tissuePeak_original


## declare an array variable
declare -a arr2=("maxPeak" "noMaxPeak")


for s in "${arr2[@]}"
do
mkdir -p $workdir/$out/$s/emOutput

for n in 1.downs 2.downs 3 4 5 6 7 8 9
do

sbatch -J INS"$n"_"$s" --error=$slurmpfx$s$n.err --output=$slurmpfx$s$n.out --workdir=$workdir --wrap="$workdir/scripts/INSIGHT_filtering.sh \
$in/emInput/sp$n.ins \
$workdir/$out/$s/sp"$n"_"$s".bed \
$workdir/$out/$s/sp"$n"_tmp.ins \
$workdir/$out/$s/sp"$n"_"$s".ins \
$in/emInput/sp$n.flankPoly.forBetas.ins \
sp"$n"_$s \
$workdir/$out/$s/emOutput/ \
$workdir/$out/$s/sp$n.flanking.filtered.ins"

   
done
done
