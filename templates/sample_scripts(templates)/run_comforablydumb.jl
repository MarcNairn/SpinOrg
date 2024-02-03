include("lib_parallel.jl")

@everywhere include("./parameters_SDe.jl")

sim = run_selfSynch(ps)

save_datal("../data/data_vs_SDe_2.jld2",sim)
