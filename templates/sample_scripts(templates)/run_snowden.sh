#!/bin/bash
#SBATCH --job-name=sS-$1
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=oneside
#SBATCH --mem=0

export JULIA_DEPOT_PATH=/cluster/userdata/luigi/.julia
julia run_parallel_snowden.jl 50 $1
