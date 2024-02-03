#!/bin/bash
for k in {1..23}; do
    for i in {1..26}; do
        # if [ $k == 1 ] && [ $i -lt 6 ]; then
        #     echo skipping job  k=$k and i=$i
        # else
              echo submitting job  k=$k and i=$i
              sbatch --time 00:$(($k*10+1200)):00 run_snowden_small.sh $k $i
        # fi
        sleep 0.5s
    done
done

# for k in {1..66}; do
#     echo submitting job k=$k
#     sbatch --time 1-06:00:00 run_snowden_small.sh $k
#     sleep 0.5s
# done
