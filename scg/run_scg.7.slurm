#!/bin/bash
#SBATCH --job-name=scg_job7
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=/gpfs/research/fangroup/rk18g/longitudinal/clustering/scg/output/exp2_2_D/scg_output7.out
#SBATCH --error=/gpfs/research/fangroup/rk18g/longitudinal/clustering/scg/output/exp2_2_D/scg_error7.out
python save_multipleSCGresults.py -opDir output/exp2_2_D -input ./input/exp2_2_DMatrix.tsv.gz -scg_config ../../SCG_Roth/scg/examples/config.yaml -niters 10000
