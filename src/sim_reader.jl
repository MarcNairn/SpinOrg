#THIS SCRIPT EXTRACTS ALL SIMULATIONS STORED AS .JLD FILES INTO A SIMULATION ARRAY FROM WHICH WE MAY EXTRACT THE OBSERVABLE sim_data
#ALSO USED TO CLEAN UP ALL THE USED JLD2 FILES AND STORE THEM IN N_TRAJECTORY SIZED BATCHES
using JLD2

include("selforg_core.jl")

#TO EXTRACT ALL SIMULATION DATA INTO A COMMON ARRAY
function load_datalll(directory::AbstractString)
    sim = []

    for file in readdir(directory)
        if endswith(file, ".jld2")
            filetemp = joinpath(directory, file)
            println("loading $filetemp...")
            push!(sim, load_datal(filetemp))
        end
    end

    merge_sim(sim...)
end

#RUN THROUGH ALL SIM DATA, CATEGORISE IT IN BATCHES OF THE SAME PARAMETERS AND SAVE IT IN MANY TRAJECTORY FILES
function sort_save_datal(directory::AbstractString, set_traj_number::Int) #take array of array of solutions, i.e. array of many trajectory simulations per entry
    #set_traj_number to some integer we are aware should match with the desired number of trajectories 

    sims_to_sort =  load_datalll(directory)
    sorted_sims = split_sim_from_par(sims_to_sort)


    for batch in sorted_sims 
        g = Int(batch[1].p.S₁) #to sort by pumping strength
        trajs = length(batch) #make sure all have same length otherwise discard current iteration
        temp = Int(batch[1].p.temp)

        #Check for correct number of trajectories:
        if trajs == set_traj_number
            # Construct the directory path
            target_directory = "sorted_sims"
        
            # Create the directory if it doesn't exist
            isdir(target_directory) || mkdir(target_directory)
            
             # Construct the file path
             file_path = joinpath(target_directory, "S=$(g), temp=$(temp), Ntraj=$(trajs)-sim_data.jld2")
            
             # Check if the file already exists
            if !isfile(file_path)
                # save data
                save_datal(file_path, batch)
                println("Creating file for g=$(g), temp=$(temp) and $(trajs) trajs in $target_directory")
            else
                println("File for g=$(g), temp=$(temp) and $(trajs) trajs already exists. Skipping...")
            end
        else
            println("Skipping batch for g=$(g), temp=$(temp) and $(trajs) trajs as it does not have the required length.")
        end
    end
end