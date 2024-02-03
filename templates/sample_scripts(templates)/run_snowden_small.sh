#!/bin/bash
#SBATCH --job-name=sS-vs-S
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=oneside

export JULIA_DEPOT_PATH=/cluster/userdata/luigi/.julia


cd $SLURM_SUBMIT_DIR

echo running job k=$1 i=$2
srun julia run_psmall_snowden.jl 4 $1 $2
