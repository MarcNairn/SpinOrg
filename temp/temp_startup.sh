#!/bin/bash

# Define the list of arguments
S_ARGS="28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47"
TEMP_ARGS="12 14 16 18" #5 10
RUN_INDEX="0 1 2 3 4 5 6 7 8 9"

# Set initial time
time_allocation="00:30:00"

# Flag to check if there are pending jobs
pending_jobs=false

for temp in $TEMP_ARGS; do
    # Iterate over the temperature arguments
    for run in $RUN_ARGS; do
        # Iterate over pump strengths
        for S in $S_ARGS; do
            # Count the number of lines in the output (excluding the header)
            num_jobs=$(echo "$(squeue)" | wc -l)
            ((num_jobs--))
            # Print the total number of pending jobs
            echo "Total number of pending jobs: $num_jobs"
            # Check if there are no pending jobs
            if [ "$num_jobs" -le 10 ]; then
                echo "Room for new jobs. Submitting a new job, run=$run, with arguments: temp=$temp..."
                # Check the value of S and adjust --time accordingly
                if [ "$S" -le 34 ]; then
                    time_allocation="00:30:00"
                elif [ "$S" -le 39 ]; then
                    # Steeper increase from 30 minutes to 2.5 hours
                    additional_minutes=$((($S - 34) * 30))
                    time_allocation="00:$(($additional_minutes + 20)):00"
                elif [ "$S" -le 44 ]; then
                    time_allocation="00:240:00"
                else
                    # Finish with a fixed 6 hours for S>44
                    time_allocation="00:360:00"
                fi
                sbatch -a 0-49 --time $time_allocation run/JUSTUS_draft_run/run_justus_static.sh $S $temp $run
                echo "New job submitted, run=$run, with arguments: S=$S, temp=$temp."
            else
                echo "Job queue full. No action taken for arguments: temp=$temp, S=$S, run=$run."
                pending_jobs=true
            fi
        done
        # Wait for all jobs to finish within the S loop
        echo "Waiting on queue to free up..."
        while [ $(echo "$(squeue)" | wc -l) -ge 10 ]; do
            sleep 1m  # Adjust the sleep interval as needed
        done
    done
done

echo "Reached end!"
