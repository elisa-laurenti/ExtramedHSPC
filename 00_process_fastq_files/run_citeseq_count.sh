#!/bin/bash
# -------------------------------------------------------
#title           :run_citeseq_count.sh
#author          :hpb29
#date            :20200818
#version         :0.9
#description     :Controls the submition of 10x libraries
#                 via SLURM to 'citeseq_count_handler'
#usage           :Edit variables below accordingly and run
# --------------------------------------------------------


# Note: Verify input/output pathways each time before running


cd $RDS/users/$USER/EL/NM/2020/SLX19841/fastq/


barcodes=(SIGAB11 SIGAD11)


# R1

for b in "${barcodes[@]}"; do find -name "${b}_S*_L001_R1*" -type f -print0 | while IFS= read -r -d '' filename; do printf '%s\n' ${filename}; done |  sort | tr '\n' '\0' | xargs -0 cat > merged/${b}_L001_R1_001.fastq.gz ;  done;


for b in "${barcodes[@]}"; do find -name "${b}_S*_L002_R1*" -type f -print0 | while IFS= read -r -d '' filename; do printf '%s\n' ${filename}; done |  sort | tr '\n' '\0' | xargs -0 cat > merged/${b}_L002_R1_001.fastq.gz ;  done;

cd merged

for b in "${barcodes[@]}"; do find -name "${b}*_R1_*" -type f -print0 | while IFS= read -r -d '' filename; do printf '%s\n' ${filename}; done |  sort | tr '\n' '\0' | xargs -0 cat > ultimate/${b}_R1.fastq.gz ;  done;

cd ..

# R2

for b in "${barcodes[@]}"; do find -name "${b}_S*_L001_R2*" -type f -print0 | while IFS= read -r -d '' filename; do printf '%s\n' ${filename}; done |  sort | tr '\n' '\0' | xargs -0 cat > merged/${b}_L001_R2_001.fastq.gz ;  done;


for b in "${barcodes[@]}"; do find -name "${b}_S*_L002_R2*" -type f -print0 | while IFS= read -r -d '' filename; do printf '%s\n' ${filename}; done |  sort | tr '\n' '\0' | xargs -0 cat > merged/${b}_L002_R2_001.fastq.gz ;  done;


cd merged

for b in "${barcodes[@]}"; do find -name "${b}*_R2_*" -type f -print0 | while IFS= read -r -d '' filename; do printf '%s\n' ${filename}; done |  sort | tr '\n' '\0' | xargs -0 cat > ultimate/${b}_R2.fastq.gz ;  done;




declare -A library
library[SIGAB11]=15364
library[SIGAD11]=14264

fastqdir=$RDS/users/$USER/EL/NM/2020/SLX19841/fastq/merged/ultimate/


# tags file
tags=$RDS/users/hpb29/EL/NM/2020/COVID19/SLX19618/project_features_CITEseqCount.csv

outdir=$RDS/users/hpb29/EL/NM/2020/SLX19841/processed/

# =================================================================


echo "Current working directory is `pwd`"

for i in "${!library[@]}"
do
  sbatch ${outdir}citeseq_count_grinder.sh ${fastqdir}${targetdir}${i}_R1.fastq.gz ${fastqdir}${targetdir}${i}_R2.fastq.gz ${tags} ${library["$i"]} ${outdir}${i}
  sleep 3
done
