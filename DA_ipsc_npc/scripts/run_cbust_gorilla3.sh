#!/bin/bash

#adjusted for running on the other server

genome=$1
TF_subset=$2

extension=500
base=/data/share/htp/pleiotropy/paper_data/DA_ipsc_npc/cbust
PFMs=$base/PWMs/$TF_subset
fasta=$base/fastas/$genome #modify the path for different subsets
out_folder=$base/cbust_out/$genome/$TF_subset #modify the path for different subsets


mkdir -p $out_folder
mkdir -p $out_folder/slurm
cd $fasta

declare -a arr=(*.fa)

for i in "${arr[@]}"


do
        # HEADER
        echo '#!/bin/bash' > $out_folder/slurm/"$i".sh  
        echo '#SBATCH -n 1' >> $out_folder/slurm/"$i".sh  
        echo '#SBATCH --error='$i'.%J.err' >> $out_folder/slurm/"$i".sh
        echo '#SBATCH --output='$i'.%J.out' >> $out_folder/slurm/"$i".sh
        echo '#SBATCH --job-name=cb_'$i'' >> $out_folder/slurm/"$i".sh

        # CBUST COMMAND
        echo "/home/zane/cluster-buster/cbust -c0 -m0 -r10000 -b$extension -f5 -G1 $PFMs $fasta/$i > $out_folder/"$i".txt" >>  $out_folder/slurm/"$i".sh

        # SUBMIT THE TASK 
        cd $out_folder/slurm/
        sbatch $out_folder/slurm/"$i".sh

done

# bash run_cbust_gorilla3.sh hg38 exprIpscNpc
# bash run_cbust_gorilla3.sh macFas6 exprIpscNpc 
# bash run_cbust_gorilla3.sh macFas6 exprTissues
# bash run_cbust_gorilla3.sh hg38 exprTissues 


