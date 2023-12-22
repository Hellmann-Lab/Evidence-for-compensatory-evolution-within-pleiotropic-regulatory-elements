#!/bin/bash

# combine full_original/split/emInput flanks and .ins in the same order into 2 files per sp only
# then run the EM algorithm to test if everything works.

workdir=/data/share/htp/pleiotropy/paper_data/INSIGHT/full_original/split
cd $workdir/emInput

cat sp1_part1.ins > ../../emInput/sp1.ins
cat sp1_part1.flankPoly.forBetas.ins > ../../emInput/sp1.flankPoly.forBetas.ins

for n in {2..4}
do
cat sp1_part$n.ins | tail -n +2 >> ../../emInput/sp1.ins
cat sp1_part$n.flankPoly.forBetas.ins | tail -n +2 >> ../../emInput/sp1.flankPoly.forBetas.ins
done


cat sp2_part1.ins > ../../emInput/sp2.ins
cat sp2_part1.flankPoly.forBetas.ins > ../../emInput/sp2.flankPoly.forBetas.ins

cat sp2_part2.ins | tail -n +2 >> ../../emInput/sp2.ins
cat sp2_part2.flankPoly.forBetas.ins | tail -n +2 >> ../../emInput/sp2.flankPoly.forBetas.ins



# run insight (actually just a check that it runs through.. not really necessary)
# slurmpfx="$workdir/combined/slurm.%J"
# workdir=/data/share/htp/pleiotropy/paper_data/INSIGHT/full_original/split
# inp = $workdir/combined/emInput/
# sbatch -o $slurmpfx --wrap="cd /home/zane/INSIGHT/; bash scripts/runINSIGHT-EM.sh sp1 $inp/sp1.ins $inp/sp1.flankPoly.forBetas.ins $workdir/combined/emOutput 15"

