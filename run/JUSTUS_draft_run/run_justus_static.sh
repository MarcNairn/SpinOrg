#!/bin/bash

#TO BE USED ONLY IN CONJUNCTION WITH THE MAIN PARAMETER LOOP, DUE TO THE $1, $2 INITIALIZATION 

# Number of nodes to allocate, always 1
#SBATCH --nodes=1
# Number of MPI instances (ranks) to be executed per node, always 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
# Configure array parameters, split job in parts labeled 0-x. (only one job x=0)
# Give job a reasonable name
#SBATCH --job-name=static_pump_Nmc
# File name for standard output (%j will be replaced by job id)
#SBATCH --output=static_pump_Nmc-%a_%A.out
# File name for error output
#SBATCH --error=static_pump_Nmc-%a_%A.err

srun julia run/JUSTUS_draft_run/run_parallel_justus_static.jl $SLURM_ARRAY_TASK_ID $1 $2 $3 $4

#ARGS[1] is the file label, running from 0 to 499, with 2 trajectories per file
#ARGS[2] is the pumping strength
#ARGS[3] is the temperature
