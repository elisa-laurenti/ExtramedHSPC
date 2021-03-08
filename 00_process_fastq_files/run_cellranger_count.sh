#!/bin/bash
# -------------------------------------------------------
#title           :run_cellranger_count.sh
#author          :hpb29
#date            :20180404
#version         :0.9    
#description     :Controls the submition of 10x libraries
#                 via SLURM to 'cellranger count'
#usage           :Edit variables below as needed and run
# --------------------------------------------------------

# Note: Some path may have changed with cluster migration.

cellranger=2.1.0
chemistry=SC3Pv2


refdir=/home/USSR/codex-pipeline/Programs/
ref=refdata-cellranger-hg19-1.2.0

# -----------------------------------------------------------------------------

# SLX-12978 (DOD1)

declare -A library
library[SIGAD9]=4319
library[SIGAE9]=4319
library[SIGAF9]=4319
library[SIGAG9]=4139
library[SIGAH9]=4319

fastqdir=/servers/bigscratch/laurenti-scratch/hb/data/nm/SLX-12978


for i in "${!library[@]}"
do
  sbatch cellranger_count_handler.sh ${cellranger} ${i} \
                                     ${refdir}${ref} ${fastqdir}${targetdir} \
                                     ${library["$i"]} ${chemistry}
  sleep 1
done

# -----------------------------------------------------------------------------

# SLX-14831 (DOD2)

declare -A library
library[SIGAB10]=9998
library[SIGAC10]=8576
library[SIGAD10]=9998

fastqdir=/servers/bigscratch/laurenti-scratch/hb/data/SLX-14831/fastq


for i in "${!library[@]}"
do
  sbatch cellranger_count_handler.sh ${cellranger} ${i} \
                                     ${refdir}${ref} ${fastqdir}${targetdir} \
                                     ${library["$i"]} ${chemistry}
  sleep 1
done


# =============================================================================


cellranger=3.1.0
chemistry=SC3Pv3

refdir=$RDS/$USER/rds-bg200-hphi-gottgens/references/10x/cellranger3x/
ref=refdata-cellranger-hg19-3.0.0


# SLX-18808 (TQ198, BP62j, BP37d)

declare -A library
library[SIGAF6]=7000
library[SIGAG6]=10000
library[SIGAH6]=10000

fastqdir=$RDS/$USER/EL/NM/2020/SLX18808/fastq/


for i in "${!library[@]}"
do
  sbatch cellranger_count_handler.sh ${cellranger} ${i} \
                                     ${refdir}${ref} ${fastqdir}${targetdir} \
                                     ${library["$i"]} ${chemistry}
  sleep 1
done

# -----------------------------------------------------------------------------

# SLX-19286 (BP74)

declare -A library
library[SIGAE2]=9000

fastqdir=/rds/user/$USER/rds-bg200-hphi-gottgens/users/$USER/EL/NM/2020/SLX19286/fastq/


for i in "${!library[@]}"
do
  sbatch cellranger_count_handler.sh ${cellranger} ${i} \
                                     ${refdir}${ref} ${fastqdir}${targetdir} \
                                     ${library["$i"]} ${chemistry}
  sleep 1
done


# -----------------------------------------------------------------------------

# SLX-19841 (BP1c, BP59h, DOD3, DOD4)

declare -A library
library[SIGAG4]=10000
library[SIGAH4]=10000
library[SIGAA12]=10000
library[SIGAB12]=10000
library[SIGAD12]=10000

fastqdir=/rds/user/$USER/rds-bg200-hphi-gottgens/users/$USER/EL/NM/2020/SLX19841/fastq/


for i in "${!library[@]}"
do
  sbatch cellranger_count_handler.sh ${cellranger} ${i} \
                                     ${refdir}${ref} ${fastqdir}${targetdir} \
                                     ${library["$i"]} ${chemistry}
  sleep 5
done
