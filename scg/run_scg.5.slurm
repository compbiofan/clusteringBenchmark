#!/bin/bash
#SBATCH --job-name=scg_job5
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=/gpfs/research/fangroup/rk18g/longitudinal/clustering/scg/output/exp2_8_D_4/scg_output5.out
#SBATCH --error=/gpfs/research/fangroup/rk18g/longitudinal/clustering/scg/output/exp2_8_D_4/scg_error5.out
python save_multipleSCGresults.py -opDir output/exp2_8_D_4 -input ./input/exp2_8_D_4.tsv.gz -scg_config ../../SCG_Roth/scg/examples/config.yaml -niters 10000
