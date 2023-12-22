#!/bin/bash

workdir=/data/share/htp/pleiotropy/paper_data/INSIGHT
slurmdir=$workdir/slurm
mkdir -p $slurmdir
slurmpfx="$slurmdir/slurm.%J"
out=PeaksAntiPeaks/
in=$workdir/PeaksAntiPeaks_original


## declare an array variable
declare -a arr2=("brain_anti" "brain_peaks" "thymus_anti" "thymus_peaks")


for s in "${arr2[@]}"
do
mkdir -p $workdir/$out/$s/emOutput

for n in {2..9}
do

sbatch -J INS"$n"_"$s" --error=$slurmpfx$s$n.err --output=$slurmpfx$s$n.out --workdir=$workdir --wrap="$workdir/scripts/INSIGHT_filtering.sh \
$in/emInput/sp"$n"_"$s".ins \
$workdir/$out/$s/sp"$n"_"$s".bed \
$workdir/$out/$s/sp"$n"_"$s"_tmp.ins \
$workdir/$out/$s/sp"$n"_"$s".ins \
$in/emInput/sp"$n"_"$s".flankPoly.forBetas.ins \
sp"$n"_$s \
$workdir/$out/$s/emOutput/ \
$workdir/$out/$s/sp$n.flanking.filtered.ins"

   
done
done
