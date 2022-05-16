#!/bin/bash
#SBATCH --job-name=scclone_job4
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=/gpfs/research/fangroup/rk18g/longitudinal/clustering/scclone/output/exp2_7_D_12/scclone_output4.out
#SBATCH --error=/gpfs/research/fangroup/rk18g/longitudinal/clustering/scclone/output/exp2_7_D_12/scclone_error4.out
time ../../scclone-1.0/bin/scclone -i ./input/exp2_7_D_12.tsv -a 0.01 -b 0.1 -o ./output/exp2_7_D_12/exp2_7_D_12
