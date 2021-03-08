#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH -p clincloud
#SBATCH -A ACCOUNTNAME-CCLD-SL2-CPU
#SBATCH --time=0-02:59:59
#SBATCH -e slurm-%j.%N.err

# -------------------------------------------------------
#title           :citeseq_count_grinder.sh
#author          :hpb29
#date            :20200818
#version         :0.9
# --------------------------------------------------------


echo "Started processing library $1 on `date`"


singularity exec jupyter.mictlan.10j.sif CITE-seq-Count -R1 $1 \
                                                        -R2 $2 \
                                                        --tags $3 \
                                                        --expected_cells $4 \
                                                        --output $5 \
                                                        -cbf 1 \
                                                        -cbl 16 \
                                                        -umif 17 \
                                                        -umil 28

echo "Finishing processing library $1 on `date`"
