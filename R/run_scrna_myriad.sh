#!/bin/bash -l
#$ -S /bin/bash
#$ -N scrna_leptomeningeal
#$ -l mem=16G
#$ -l h_rt=48:00:00
#$ -pe smp 8
#$ -wd /home/sejkmor/Scratch/taylordataset/taylor25
#$ -o /home/sejkmor/Scratch/taylordataset/taylor25/logs/job.out
#$ -e /home/sejkmor/Scratch/taylordataset/taylor25/logs/job.err
#$ -m bea
#$ -M sejkmor@ucl.ac.uk

mkdir -p /home/sejkmor/Scratch/taylordataset/taylor25/logs

module --force purge
module load r/4.5.1-openblas/gnu-10.2.0

echo "Job started: $(date)"
echo "Running on: $(hostname)"
echo "Cores: $NSLOTS"

Rscript /home/sejkmor/Scratch/taylordataset/taylor25/scrna_leptomeningeal_analysis.R

echo "Job finished: $(date)"
