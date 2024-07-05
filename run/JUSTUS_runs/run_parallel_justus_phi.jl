using Pkg
Pkg.activate("..")

using Distributed
using ClusterManagers

np = parse(Int,ARGS[1])
i = parse(Int,ARGS[2])

addprocs(SlurmManager(np))

@everywhere push!(LOAD_PATH,"../src")
@everywhere include("../load.jl")

using Random

# Random.seed!(232344)

@everywhere include("../parameters_phi.jl")

sim = many_trajectory_solver(ps[i],saveat=10.0, seed=abs(rand(Int), maxiters=Int(1.0e9)))
save_datal("data_vs_phi_S4_part" * string(i) * ".jld2",sim)


# for (i, psi) in enumerate(ps)
#     sim = run_selfSynch(psi,abs(rand(Int)))
#     save_datal("../data/data_vs_phi_S4_part" * string(i) * ".jld2",sim)
# end

for i in workers()
	rmprocs(i)
end
