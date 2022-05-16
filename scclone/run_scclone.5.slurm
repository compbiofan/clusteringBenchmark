#!/bin/bash
#SBATCH --job-name=scclone_job5
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=/gpfs/research/fangroup/rk18g/longitudinal/clustering/scclone/output/exp2_8_D_4/scclone_output5.out
#SBATCH --error=/gpfs/research/fangroup/rk18g/longitudinal/clustering/scclone/output/exp2_8_D_4/scclone_error5.out
time ../../scclone-1.0/bin/scclone -i ./input/exp2_8_D_4.tsv -a 0.01 -b 0.1 -o ./output/exp2_8_D_4/exp2_8_D_4
