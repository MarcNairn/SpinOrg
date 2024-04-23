#!/bin/bash

# Wait until squeue has only one line (meaning only the header is present)
while [ $(squeue | wc -l) -gt 1 ]; do
    sleep 1  # Sleep for 1 second before checking again
done

# Run Julia file with argument passed from the command line
julia run/phase_diagram_sims2df.jl "$1"
