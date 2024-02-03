#!/bin/bash
for i in {1..31}; do
    echo submitting job $i
		sbatch run_snowden.sh $i
		sleep 2s
done
