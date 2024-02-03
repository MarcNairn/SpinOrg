include("lib_parallel.jl")
using Random

@everywhere include("./parameters_phi.jl")

Random.seed!(234)

for (i, psi) in enumerate(ps)
    sim = run_selfSynch(psi,abs(rand(Int)))
    save_datal("/localdisk/luigi/data/selfSynch/system1/data_vs_phi_orc_part" * string(i) * ".jld",sim)
end

