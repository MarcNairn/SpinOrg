### THIS SCRIPT USES THE FUNCTIONALITIES DEFINED IN selforg_core.jl TO RUN THE NUMERICS OVER A RANGE OF PARAMETERS ###

using JLD2
using Distributed


@everywhere include("../src/selforg_core.jl")

g_init_value = 30
g_final_value = 35
t_range = 800
N_MC = 100
N_spin = 100

p_array = vcat([fill(System_p(0.0, 0.0, i , i , 10.0, -100.0, 100.0, 10.0, N_spin, (0.0, t_range), N_MC)) for i in g_init:g_final]...)


sim_array = vcat([fill(many_trajectory_solver(p_array[i], saveat=10.0, seed=abs(rand(Int)), maxiters=Int(1e9))) for i in 1:length(p_array)]...)

# this object may be used for PLOTTING
sim_extract_array = vcat([extract_solution(sim_array[i]) for i in 1:length(sim_array)]...)


# SAVE THE sim_array
save_datal("sim_data/pump_range/g($(g_init)-$(g_final))_range_data(t_range=$(t_range),N_MC=($N_MC).jld", sim_extract_array)
