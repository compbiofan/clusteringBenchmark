#!/bin/bash
#SBATCH --job-name=bnpc_job1
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=/gpfs/research/fangroup/rk18g/longitudinal/clustering/bnpc/output/exp2_13_D_9_bnpc/bnpc_output1.out
#SBATCH --error=/gpfs/research/fangroup/rk18g/longitudinal/clustering/bnpc/output/exp2_13_D_9_bnpc/bnpc_error1.out
time python ../../BnpC/run_BnpC.py ./input/exp2_13_D_9_bnpc.tsv -pp 0.75 0.75 -o ./output/exp2_13_D_9_bnpc
