#module CustomFunctions

#REMOVED MODULE CALLS

using DifferentialEquations: EnsembleSolution, RODESolution
using JLD, JLD2

export extract_solution, save_datal, load_datal, Sol,merge_sol, merge_sim
export intervallize_array, load_datall
export yourpart

struct Sol
    u::Array{Float64,2}
    p
    t::Array{Float64,1}
    alg::String
end

Sol(u,p,t) = Sol(u,p,t,"")

function extract_solution(sol::RODESolution)
    [Sol(hcat(sol.u...),sol.prob.p,sol.t,repr(sol.alg))]
end

function extract_solution(sim::EnsembleSolution)
    trajectories = Array{Sol,1}()
    for sol in sim
        append!(trajectories,extract_solution(sol))
    end
    trajectories
end

function extract_solution(sim::Array{Sol,1})
    sim
end

# use vcat to merge Array{Sol,1}

function merge_sim(sim1::Array{Sol,1}...)
    vcat(sim1...)
end

function merge_sol(sol1::Sol,sol2::Sol)
    # perform checks
    if sol1.p != sol2.p
        error("solutions do not have same parameters")
    end
    if sol1.u[:,end] != sol2.u[:,1]
        error("endpoint in sol1 differ from initial point int sol2")
    end
    t = vcat(sol1.t,sol2.t[2:end])
    Sol(hcat(sol1.u,sol2.u[:,2:end]),sol1.p,t,sol1.alg*sol2.alg)
end

function merge_sol(sol1::Sol,sol2::Sol,sols::Sol...)
    sol = merge_sol(sol1,sol2)
    for x in sols
        sol = merge_sol(sol,x)
    end
    sol
end

function merge_sol(sim1::Array{Sol,1},sim2::Array{Sol,1})
    if size(sim1) != size(sim2)
        error("sim1 and sim2 must have the same size")
    end
    sim = Sol[]
    for i in 1:size(sim1)[1]
        push!(sim, merge_sol(sim1[i],sim2[i]))
    end
    sim
end

function merge_sol(sim1::Array{Sol,1},sim2::Array{Sol,1},sims::Array{Sol,1}...)
    sim = merge_sol(sim1,sim2)
    for x in sims
        sim = merge_sol(sim,x)
    end
    sim
end


function keep_input_time(sim::Array{Sol,1}, required_t::Int)
    sorted_sim = Sol[]
    for i in 1:length(sim)
        if sim[i].t[end]==required_t
            push!(sorted_sim, sim[i])
        end
    end
    return sorted_sim
end


function save_datal(filename::String,solorsim)
    print("saving data in $(filename) ..")
    elt = @elapsed JLD.jldopen(filename, "w") do file
        write(file, "solution", extract_solution(solorsim))
    end
    println("done in $elt seconds.")
end

function load_datal(filename::String)
    # jldopen(filename, true, true, true, IOStream) do file
    try
        JLD.jldopen(filename, "r") do file
            read(file, "solution")
        end
    catch
        @warn("using JLD2")
        JLD2.jldopen(filename, "r") do file
            file["solution"]
        end
    end
end

function load_datal(filename::String,n::Int)
    sim = []
    for i in 1:n
        filetemp = filename * string(i) * ".jld2"
        println("loading "* filetemp* "...")
        push!(sim,load_datal(filetemp))
    end
    merge_sim(sim...)
end

function load_datall(path::String)
    sim = []

    dir = dirname(path)
    filename = basename(path)

    for i in readdir(dir)
        if length(i) >= length(filename) && i[1:length(filename)] == filename
            filetemp = joinpath(dir,i)
            println("loading "* filetemp* "...")
            push!(sim,load_datal(filetemp))
        end
    end
    merge_sim(sim...)
end

#new addition below
#=
function load_datalll(path::String)
    sim = []
    dir = dirname(path)
    filename = basename(path)

    files = filter(f -> startswith(f, filename), readdir(dir))

    for file in files
        filetemp = joinpath(dir, file)
        println("loading $filetemp...")
        push!(sim, load_datal(filetemp))
    end

    merge_sim(sim...)
end
=#
function intervallize_array(arr::Array{Float64,1})
    y = zeros(length(arr)+1)
    y[2:end] .= arr
    y[1:end-1] .+= arr
    y = y/2
    y[1] = 2arr[1]-y[2]
    y[end] = 2arr[end]-y[end-1]
    return y
end

function yourpart(x::AbstractArray, n::Int)
    [x[i:min(i+n-1,length(x))] for i in 1:n:length(x)]
end


#THE FOLLOWING EXTRACTS ALL SIMULATIONS STORED AS .JLD FILES INTO A SIMULATION ARRAY FROM WHICH WE MAY EXTRACT THE OBSERVABLE sim_data
#ALSO USED TO CLEAN UP ALL THE USED JLD2 FILES AND STORE THEM IN N_TRAJECTORY SIZED BATCHES

#TO EXTRACT ALL SIMULATION DATA IN DIRECTORY INTO A COMMON ARRAY
function load_datalll(directory::AbstractString)
    sim = []
    println("loading data...")
    for file in readdir(directory)
        if endswith(file, ".jld2")
            filetemp = joinpath(directory, file)
            #println("loading $filetemp...")
            push!(sim, load_datal(filetemp))
        end
    end
    println("all files loaded!")
    merge_sim(sim...)
end

#load all simulations of a given temperature value
function load_data_temp(directory::AbstractString, temp_value::Int)
    sim = []
    println("loading data...")
    for file in readdir(directory)
        if endswith(file, ".jld2") && contains(file, "temp=$temp_value,") #need the comma so that only unique integer values are picked and not stringed, i.e. 1 and not 1,10,12...
            filetemp = joinpath(directory, file)
            #println("loading $filetemp...")
            push!(sim, load_datal(filetemp))
        end
    end
    println("all files loaded!")
    merge_sim(sim...)
end

#load all simulations of a given atomic frequency value
function load_data_delta(directory::AbstractString, Delta_e::Int)
    sim = []
    println("loading data...")
    for file in readdir(directory)
        if endswith(file, ".jld2") && contains(file, "Delta_e=$Delta_e,") #need the comma so that only unique integer values are picked and not stringed, i.e. 1 and not 1,10,12...
            filetemp = joinpath(directory, file)
            #println("loading $filetemp...")
            push!(sim, load_datal(filetemp))
        end
    end
    println("all files loaded!")
    merge_sim(sim...)
end

# function load_data_deltatemp(directory::AbstractString, Delta_e::Int, temp::Int)
#     sim = []
#     println("loading data...")
#     elt=@elapsed for file in readdir(directory)
#         if endswith(file, ".jld2") && contains(file, "Delta_e=$Delta_e,") && contains(file, "temp=$temp,") #need the comma so that only unique integer values are picked and not stringed, i.e. 1 and not 1,10,12...
#             filetemp = joinpath(directory, file)
#             println("loading $filetemp...")
#             push!(sim, load_datal(filetemp))
#         end
#     end
#     println("all files loaded in $elt seconds!")
#     merge_sim(sim...)
# end


function load_data_deltatemp_alt(directory::AbstractString, temp::Int) #Delta_e=10 files do not contain the string Delta_e, should fix later
    sim = []
    println("loading data...")
    elt=@elapsed for file in readdir(directory)
        if endswith(file, ".jld2") && !contains(file, "Delta_e=") && contains(file, "temp=$temp,") #need the comma so that only unique integer values are picked and not stringed, i.e. 1 and not 1,10,12...
            filetemp = joinpath(directory, file)
            println("loading $filetemp...")
            push!(sim, load_datal(filetemp))
        end
    end
    println("all files loaded in $elt seconds!")
    merge_sim(sim...)
end


function load_data_deltatemp_alt(directory::AbstractString, S::Int, temp::Int) #Delta_e=10 files do not contain the string Delta_e, should fix later
    sim = []
    println("loading data...")
    elt=@elapsed for file in readdir(directory)
        if endswith(file, ".jld2") && !contains(file, "Delta_e=") && contains(file, "temp=$temp,") && contains(file, "S=$S,") #need the comma so that only unique integer values are picked and not stringed, i.e. 1 and not 1,10,12...
            filetemp = joinpath(directory, file)
            println("loading $filetemp...")
            push!(sim, load_datal(filetemp))
        end
    end
    println("all files loaded in $elt seconds!")
    merge_sim(sim...)
end


#RUN THROUGH ALL SIM DATA, CATEGORISE IT IN BATCHES OF THE SAME PARAMETERS AND SAVE IT IN MANY TRAJECTORY FILES
function sort_save_datal(directory::AbstractString, Delta_e::Int, temps::Vector{Int}) #take array of array of solutions, i.e. array of many trajectory simulations per entry
    #set_traj_number to some integer we are aware should match with the desired number of trajectories 

    for temp in temps
        data = load_data_deltatemp(directory::AbstractString, Delta_e::Int, temp::Int)
            for batch in data 
                g = Int(batch[1].p.S₁) #to sort by pumping strength

                # Construct the directory path
                target_directory = "sorted_sims"
            
                # Create the directory if it doesn't exist
                isdir(target_directory) || mkdir(target_directory)
                
                # Construct the file path
                file_path = joinpath(target_directory, "complete_S=$(g), Delta_e=$Delta_e, temp=$(temp).jld2")
                
                # Check if the file already exists
                if !isfile(file_path)
                    # save data
                    save_datal(file_path, batch)
                    println("Creating file for g=$(g), temp=$(temp) in $target_directory")
                else
                    println("File for g=$(g), temp=$(temp) already exists. Skipping...")
                end
            end
    end
end
#end

################################
#NEW INITIAL CONDITIONS; BINORND SAMPLED spins

function load_data_binordn(directory::AbstractString, temp::Int) # to import all temp data
    sim = []
    println("loading data...")
    for file in readdir(directory)
        if endswith(file, ".jld2") && contains(file, "binomial_init_spin") && contains(file, "temp=$temp,") #need the comma so that only unique integer values are picked and not stringed, i.e. 1 and not 1,10,12...
            filetemp = joinpath(directory, file)
            #println("loading $filetemp...")
            push!(sim, load_datal(filetemp))
        end
    end
    println("all files loaded!")
    merge_sim(sim...)
end

function load_data_binordn(directory::AbstractString, temp::Int, S1::Int, S2::Int, Delta_e::Int) #for two point Normal/SR phase plotting
    sim = []
    println("loading data...")
    for file in readdir(directory)
        if endswith(file, ".jld2") && contains(file, "binomial_init_spin") && contains(file, "temp=$temp,") && contains(file, "Delta_e=$Delta_e,") && (contains(file, "S=$S1") || contains(file, "S=$S2"))#need the comma so that only unique integer values are picked and not stringed, i.e. 1 and not 1,10,12...
            filetemp = joinpath(directory, file)
            #println("loading $filetemp...")
            push!(sim, load_datal(filetemp))
                    end
    end
    println("all files loaded!")
    merge_sim(sim...)
end


function load_data_binordn(directory::AbstractString, temp::Int, Slist::Tuple, Delta_e::Int) #for S range plotting
    sim = []
    Srange = Slist[1]:1:Slist[end]
    println("loading data...")
    for file in readdir(directory)
        for S in Srange
            if endswith(file, ".jld2") && contains(file, "binomial_init_spin") && contains(file, "temp=$temp,") && contains(file, "S=$S,") && contains(file, "Delta_e=$Delta_e,")#need the comma so that only unique integer values are picked and not stringed, i.e. 1 and not 1,10,12...
                filetemp = joinpath(directory, file)
                #println("loading $filetemp...")
                push!(sim, load_datal(filetemp))
            end
        end
    end
    println("all files loaded!")
    merge_sim(sim...)
end

function load_data_binordn(directory::AbstractString, Nat::Int, temp::Int, S1::Int, S2::Int, Delta_e::Int) #for normal/SR data for a given atom number
    sim = []
    println("loading data...")
    for file in readdir(directory)
        if endswith(file, ".jld2") && contains(file, "binomial_init_spin") && contains(file, "temp=$temp,") && contains(file, "Nat=$Nat,") && contains(file, "Delta_e=$Delta_e,") && (contains(file, "S=$S1") || contains(file, "S=$S2"))#need the comma so that only unique integer values are picked and not stringed, i.e. 1 and not 1,10,12...
            filetemp = joinpath(directory, file)
            #println("loading $filetemp...")
            push!(sim, load_datal(filetemp))
        end
    end
    println("all files loaded!")
    merge_sim(sim...)
end

function load_data_binordn(directory::AbstractString, Nat::Int, temp::Int, Slist::Tuple, Delta_e::Int) #for normal/SR data for a given atom number
    sim = []
    println("loading data...")
    for file in readdir(directory)
        for S in Slist[1]:1:Slist[end]
            if endswith(file, ".jld2") && contains(file, "binomial_init_spin") && contains(file, "temp=$temp,") && contains(file, "Nat=$Nat,") && contains(file, "Delta_e=$Delta_e,") && contains(file, "S=$S,") #need the comma so that only unique integer values are picked and not stringed, i.e. 1 and not 1,10,12...
                filetemp = joinpath(directory, file)
                #println("loading $filetemp...")
                push!(sim, load_datal(filetemp))
            end
        end
    end
    println("all files loaded!")
    merge_sim(sim...)
end


function load_data_pm(directory::AbstractString, Nat::Int, temp::Int, Slist::Tuple, Delta_e::Int) #for normal/SR data for a given atom number
    sim = []
    println("loading data...")
    for file in readdir(directory)
        for S in Slist[1]:1:Slist[end]
            if endswith(file, ".jld2") && contains(file, "pm_init_spin_badcav") && contains(file, "temp=$temp,") && contains(file, "Nat=$Nat,") && contains(file, "Delta_e=$Delta_e,") && contains(file, "S=$S,") #need the comma so that only unique integer values are picked and not stringed, i.e. 1 and not 1,10,12...
                filetemp = joinpath(directory, file)
                push!(sim, load_datal(filetemp))
            else
                nothing
            end
        end
    end
    println("all files loaded!")
    merge_sim(sim...)
end

function load_data_pm(directory::AbstractString, Nat::Int, temp::Int, S::Int, Delta_e::Int)
    files_to_load = filter(file ->
            endswith(file, ".jld2") &&
            contains(file, "pm_init_spin_badcav") &&
            contains(file, "temp=$temp,") &&
            contains(file, "Nat=$Nat,") &&
            contains(file, "Delta_e=$Delta_e,") &&
            contains(file, "S=$S,"),
        readdir(directory)
    )    
    sim = [load_datal(joinpath(directory, file)) for file in files_to_load]
    
    println("All files loaded!")
    
    merge_sim(sim...)
end

function load_data_pm_generic(directory::AbstractString, custom_string::AbstractString)
    sim = []
    println("loading data...")
    for file in readdir(directory)
        if endswith(file, ".jld2") && contains(file, "$(custom_string)")
             filetemp = joinpath(directory, file)
             push!(sim, load_datal(filetemp))
        else
            nothing
        end
    end
    println("all files loaded!")
    merge_sim(sim...)
end

function load_merge_trajs(directory::AbstractString, Nat::Int, temp::Int, Slist::Tuple, Delta_e::Int) #take array of array of solutions, i.e. array of many trajectory simulations per entry

    # Construct the directory path
    target_directory = "sorted_sims_pm"

    # Create the directory if it doesn't exist
    isdir(target_directory) || mkdir(target_directory)

    for S in Slist[1]:1:Slist[end]
        batch = load_data_pm(directory::AbstractString, Nat::Int, temp::Int, S::Int, Delta_e::Int)

        # Construct the file path
        file_path = joinpath(target_directory, "pm_complete_badcav(S=$(g), Delta_e=$Delta_e, temp=$(temp), Nat=$(Nat)).jld2")
        
        # Check if the file already exists
        if !isfile(file_path)
            # save data
            save_datal(file_path, batch)
            println("Creating file for g=$(g), temp=$(temp), Delta_e=$Delta_e in $target_directory")
        else
            println("File for g=$(g), temp=$(temp), Delta_e=$Delta_e already exists. Skipping...")
        end
    end
end