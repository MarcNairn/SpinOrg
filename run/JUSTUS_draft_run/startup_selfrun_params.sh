#!/bin/bash

## ALTERNATE SCRIPT
echo "Waiting on queue to empty..."
    while [ $(echo "$(squeue)" | wc -l) -gt 1 ]; do
        sleep 30s  # Adjust the sleep interval as needed
    done
echo "Starting..."

# Define the list of arguments
S_ARGS=$(seq 37 65) # 11 12 13 14 15 16 17 18 19 46 47 48 49 50 51 52 53 54 55 56 57 58 59" #29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47"
TEMP_ARGS="20" #20" # 30 40"
DELTAE_ARGS="2"
ATOM_NUMBER="100" #"10 50 100 1000"


# Set initial time
time_allocation="00:30:00"

for Nat in $ATOM_NUMBER; do
#ierate over atom numbers
    for temp in $TEMP_ARGS; do
    #iterate (if needed) over temperatures
        for Delta_e in $DELTAE_ARGS; do
        # Iterate over the atomic frequency arguments
            for S in $S_ARGS; do
                num_jobs=$(echo "$(squeue)" | wc -l)
                ((num_jobs--))
                # Print the total number of pending jobs
                echo "Total number of pending jobs: $num_jobs"

                # Check if there are no pending jobs
                if [ "$num_jobs" -lt 1 ]; then
                    echo "No pending jobs"
                    if [ "$Nat" -le 100 ]; then
                        # Check the value of S and adjust --time accordingly
                        if [ "$S" -le 20 ]; then
                            time_allocation="00:10:00"
                        elif [ "$S" -le 25 ]; then
                            additional_minutes=$((($S - 20) * 10))
                            time_allocation="00:$(($additional_minutes + 10)):00"
                        else
                            time_allocation="00:180:00"
                        fi
                    else
                        #increase time allowance when atom number is increased
                        if [ "$S" -le 20 ]; then
                            time_allocation="05:00:00"
                        elif [ "$S" -le 25 ]; then
                            # additional_minutes=$((($S - 20) * 20))
                            # time_allocation="00:$(($additional_minutes + 30)):00"
                            time_allocation="10:00:00"
                        else
                            time_allocation="15:00:00"
                        fi
                    fi
                    # Record start time before submitting the job
                    start_time=$(date +"%s")

                    sbatch -a 0-999 --time $time_allocation run/JUSTUS_draft_run/run_justus_static.sh $S $temp $Delta_e $Nat
                    echo "New job submitted, with arguments: N_at=$Nat, Delta_e=$Delta_e, S=$S, temp=$temp."

                else
                    echo "ERROR: Job queue overran! Current job (S=$S, Delta_e=$Delta_e, temp=$temp) will be skipped..."
                fi

            echo "Waiting on jobs to finish..."
                while [ $(echo "$(squeue)" | wc -l) -gt 1 ]; do
                    sleep 1s  # Adjust the sleep interval as needed
                done

                # Calculate and print elapsed time for the job
                end_time=$(date +"%s")
                elapsed_time=$((end_time - start_time))
                echo "Elapsed time for S=$S, temp=$temp: $elapsed_time seconds"
                
            done
        done
    done
done
echo "Reached end!"



