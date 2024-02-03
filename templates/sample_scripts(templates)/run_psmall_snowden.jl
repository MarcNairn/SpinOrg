using Pkg
Pkg.activate("..")

using Distributed

np = parse(Int,ARGS[1])
k = parse(Int,ARGS[2])
i = parse(Int,ARGS[3])

addprocs(np)

@everywhere push!(LOAD_PATH,"../src")
@everywhere include("./load.jl")

using Random

# Random.seed!(232344)

@everywhere include("./parameters.jl")

# sim = run_selfSynch(p,maxiters=Int(1e11))
# save_datal("../data/ZZZ_data_De10_S5_Tk_phipi3_part" * string(k) * ".jld",sim)


sim = run_selfSynch(ps[k])
save_datal("../data/newinitialconditions/vs_S/data_De1_Tk" * string((k-1)*26+i) * ".jld",sim)



# for (i, psi) in enumerate(ps)
#     sim = run_selfSynch(psi,abs(rand(Int)))
#     save_datal("../data/data_vs_phi_S4_part" * string(i) * ".jld2",sim)
# end

for i in workers()
	rmprocs(i)
end
