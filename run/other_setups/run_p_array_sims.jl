include("lib_parallel.jl")
include("../src/custom_functions.jl")

@everywhere include("./parameters_SDe.jl")

sim = many_trajectory_solver(ps)

save_datal("../sim_data/data_vs_SDe_2.jld",sim)
