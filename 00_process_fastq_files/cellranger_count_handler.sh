#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH -p clincloud
#SBATCH -A ACCOUNTNAME-SL2-CPU
#SBATCH --time=0-9:59:59
#SBATCH -e slurm.%N.%j.err

# -------------------------------------------------------
#title           :cellranger_count_handler.sh
#author          :hpb29
#date            :20180404
#version         :0.9

#description     :This script submits 10x libraries via
#                 SLURM to 'cellranger count' and expects 
#                 to be fed four positional arguments:
#                 $1 - cellranger version
#                 $2 - library/sample id
#                 $3 - path to 10x reference
#                 $4 - path to fastq files dir
#                 $5 - expected number of cells
#                 $6 - chemistry

#                 and was designed to be used via a
#                 controller script. 
#                 (run_cellranger_count.sh)
# --------------------------------------------------------

module load cellranger/$1


echo "Started processing library $2 on `date`"

cellranger count --id=$2 \
                 --sample=$2 \
                 --transcriptome=$3\
                 --fastqs=$4 \
                 --expect-cells=$5\
                 --chemistry=$6

echo "Finishing processing library $2 on `date`"

