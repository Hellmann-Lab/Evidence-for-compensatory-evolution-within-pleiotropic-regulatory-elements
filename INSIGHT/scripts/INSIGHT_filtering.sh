#!/bin/bash


#change 1st & 2nd column, split up in two, hashtag blocks

sp=$1
filter_by=$2  # distance_to_TSS/distal_regions2.txt
sp_filtered=$3
sp_filtered_final=$4
flanking_regions=$5
dataID=$6
emOutput=$7 # distance_to_TSS/sp_promoter/emOutput/
flanking_regions_filtered=$8

perl_script=/data/share/htp/pleiotropy/paper_data/INSIGHT/scripts/filter_sites.pl
workdir=/data/share/htp/pleiotropy/paper_data/INSIGHT/

cd $workdir

/usr/bin/perl -F\\t -nlae 'print join("\t", @F[1,0,2..$#F])' $sp > tmp$dataID.ins
tr ':' $'\t' < tmp$dataID.ins > tmp2$dataID.ins
sed -i '/block/s/^/#/' tmp2$dataID.ins
 
#now filter 
/usr/bin/perl $perl_script $filter_by tmp2$dataID.ins $sp_filtered
 
sed -i 's/#//g' $sp_filtered

awk '{$1=$1":"$2;$2=""}1' $sp_filtered > tmp"$dataID"_filtered.ins

awk '{t=$1; $1=$2; $2=t; print; } '  tmp"$dataID"_filtered.ins > tmp2"$dataID"_filtered.ins

tr ' ' $'\t' < tmp2"$dataID"_filtered.ins > tmp3"$dataID"_filtered.ins

sed '/M/ s/$/\t\t\t/' tmp3"$dataID"_filtered.ins > $sp_filtered_final

sed -i '1 i\samples 108' $sp_filtered_final

#diff -q $sp_filtered_final sp1_distal_myComp/sp1.distal.ins

#selecting flanking stuff
awk -F'\t' 'NR==FNR{a[$0]; next} $0 in a {print $0;p=NR} NR==p+1 {while ($0 ~ "site"){print $0,p=NR, NR==p+1; next}}' $sp_filtered_final  $flanking_regions > temp_flanking$dataID.ins
sed 's/ .*//' temp_flanking$dataID.ins > $flanking_regions_filtered
sed -i '1s/.*/samples 108/' $flanking_regions_filtered
rm tmp*$dataID*
rm temp*$dataID*

#awk -F'\t' 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' $sp_filtered_final $flanking_regions > $flanking_regions_filtered

cd /home/zane/INSIGHT/
bash scripts/runINSIGHT-EM.sh $dataID $sp_filtered_final $flanking_regions_filtered $emOutput 15




#examples for each group

# zane@gorilla4:~/INSIGHT$ bash filtering3.sh emInput/sp1.ins distance_to_TSS/distal_regions2.txt 
# distance_to_TSS/sp_distal/sp1_tmp.ins distance_to_TSS/sp_distal/sp1_distal.ins emInput/sp1.flankPoly.forBetas.ins sp1_distal distance_to_TSS/sp_distal/emOutput/ distance_to_TSS/sp_distal/sp1.flanking.filtered.ins

# zane@gorilla4:~/INSIGHT$ bash filtering3.sh emInput/sp1.ins CpG_content/lowCpG_regions2.txt 
# CpG_content/low/sp1_tmp.ins CpG_content/low/sp1_lowCpG.ins emInput/sp1.flankPoly.forBetas.ins sp1_lowCpG CpG_content/low/emOutput/ CpG_content/low/sp1.flanking.filtered.ins





