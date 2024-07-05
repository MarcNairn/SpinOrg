#!/bin/bash
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=4G
# Give job a reasonable name
#SBATCH --job-name=sS-$1
# File name for standard output (%j will be replaced by job id)
#SBATCH --output=serial_job-%j.out
# File name for error output
#SBATCH --error=serial_job-%j.err

julia run/JUSTUS_draft_run/run_parallel_justus_phi.jl 50 $1
