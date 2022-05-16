#!/bin/bash
#SBATCH --job-name=bnpc_job7
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=/gpfs/research/fangroup/rk18g/longitudinal/clustering/bnpc/output/exp2_2_D_bnpc/bnpc_output7.out
#SBATCH --error=/gpfs/research/fangroup/rk18g/longitudinal/clustering/bnpc/output/exp2_2_D_bnpc/bnpc_error7.out
time python ../../BnpC/run_BnpC.py ./input/exp2_2_DMatrix_val_bnpc.tsv -pp 0.75 0.75 -o ./output/exp2_2_D_bnpc
