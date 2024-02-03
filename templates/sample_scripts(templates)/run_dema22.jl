include("lib_parallel.jl")
using Random

@everywhere include("./parameters_phi.jl")

Random.seed!(232344)

for (i, psi) in enumerate(ps)
    sim = run_selfSynch(psi,abs(rand(Int)))
    save_datal("../data/data_vs_phi_S4_part" * string(i) * ".jld2",sim)
end
