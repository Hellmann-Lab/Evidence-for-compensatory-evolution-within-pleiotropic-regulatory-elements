#!/bin/bash

workdir=/data/share/htp/pleiotropy/paper_data
slurm=$workdir/DA_ipsc_npc/cbust/RDS/TFBS_position/slurms
slurmscripts=$workdir/DA_ipsc_npc/cbust/RDS/TFBS_position/slurmscripts
motif_filt_file=$workdir/DA_ipsc_npc/cbust/top10perc_motifs_exprTissues_dt.csv

mkdir -p $slurm
mkdir -p $slurmscripts


# need to select only IDs where both species got SOME binding, otherwise it will throw an error..
cd $workdir/DA_ipsc_npc/cbust/RDS/TFBS_position
echo "$(ls cbust_sep_exprTissues/hg38/)" > all_cbust_h.txt
echo "$(ls cbust_sep_exprTissues/macFas6/)" > all_cbust_m.txt
awk 'NR==FNR{a[$0];next}($0 in a)' all_cbust_h.txt all_cbust_m.txt > all_cbust.txt
sed -i -e "s/\.rds//" all_cbust.txt
grep -Fwf all_cbust.txt selected_ids.txt > selected_ids_cbust.txt
num="$(cat selected_ids_cbust.txt | wc -l)"

#this is set up for 90 nodes
for i in {1..389757..4330}
  do
    echo '#!/bin/bash' > $slurmscripts/$i.sh
    echo 'workdir=/data/share/htp/pleiotropy/paper_data' >> $slurmscripts/$i.sh
    echo 'slurm='$workdir'/DA_ipsc_npc/cbust/RDS/TFBS_position/slurms' >> $slurmscripts/$i.sh
     echo 'sed -n '$i','$((i+4330))'p '$workdir'/DA_ipsc_npc/cbust/RDS/TFBS_position/selected_ids_cbust.txt |  while read m j k;
    do Rscript '$workdir'/DA_ipsc_npc/scripts/calculatePositionCons.R $m $j $k $motif_filt_file; done' >> $slurmscripts/$i.sh
    sbatch -o $slurm/slurm-$i-%A.out $slurmscripts/$i.sh
  done