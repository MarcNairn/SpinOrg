using JLD2
include("../../src/sim_reader.jl")


directory = ARGS[1]
set_traj_number = parse(Int, ARGS[2])



sort_save_datal(directory, set_traj_number)

