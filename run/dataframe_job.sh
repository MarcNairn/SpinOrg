#!/bin/bash
#SBATCH --job-name=sims2df
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=250G
#SBATCH --time=0-01:00:00

srun run/phase_diagram_sims2df.sh 800