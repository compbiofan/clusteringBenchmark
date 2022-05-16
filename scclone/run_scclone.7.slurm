#!/bin/bash
#SBATCH --job-name=scclone_job7
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=/gpfs/research/fangroup/rk18g/longitudinal/clustering/scclone/output/exp2_2_D/scclone_output7.out
#SBATCH --error=/gpfs/research/fangroup/rk18g/longitudinal/clustering/scclone/output/exp2_2_D/scclone_error7.out
time ../../scclone-1.0/bin/scclone -i ./input/exp2_2_DMatrix.tsv -a 0.01 -b 0.1 -o ./output/exp2_2_D/exp2_2_D
