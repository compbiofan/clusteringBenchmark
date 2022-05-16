#!/bin/bash
#SBATCH --job-name=scclone_job0
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=/gpfs/research/fangroup/rk18g/longitudinal/clustering/scclone/output/exp2_11_D_20/scclone_output0.out
#SBATCH --error=/gpfs/research/fangroup/rk18g/longitudinal/clustering/scclone/output/exp2_11_D_20/scclone_error0.out
time ../../scclone-1.0/bin/scclone -i ./input/exp2_11_D_20.tsv -a 0.01 -b 0.1 -o ./output/exp2_11_D_20/exp2_11_D_20
