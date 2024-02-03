addprocs([("dema70",8),("dema48",8),("dema62",8)])
# addprocs([("ichthys2",44),("ichthys4",44)])

include("lib_parallel.jl")
using Random

@everywhere include("./parameters_S_2.jl")

Random.seed!(241)

for (i, psi) in enumerate(ps)
    sim = run_selfSynch(psi,abs(rand(Int)))
    save_datal("/localdisk/luigi/data/selfSynch/system1/data_vs_S_2_Tk_part" * string(i) * ".jld2",sim)
end
