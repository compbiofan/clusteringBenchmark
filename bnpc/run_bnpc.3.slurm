#!/bin/bash
#SBATCH --job-name=bnpc_job3
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=/gpfs/research/fangroup/rk18g/longitudinal/clustering/bnpc/output/exp2_2_D_10_bnpc/bnpc_output3.out
#SBATCH --error=/gpfs/research/fangroup/rk18g/longitudinal/clustering/bnpc/output/exp2_2_D_10_bnpc/bnpc_error3.out
time python ../../BnpC/run_BnpC.py ./input/exp2_2_D_10_bnpc.tsv -pp 0.75 0.75 -o ./output/exp2_2_D_10_bnpc
