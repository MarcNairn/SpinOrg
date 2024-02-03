include("lib_parallel.jl")
using Random

@everywhere include("./parameters_S.jl")

# sim = run_selfSynch(ps,1847)
# save_datal("/localdisk/luigi/data/selfSynch/system1/data_vs_S_ordering.jld2",sim)

# maybe not useful
# ps_parts = yourpart(ps,num_monte)

Random.seed!(2341)

for (i, psi) in enumerate(ps)
    if i<=42
        println("random number "*string(i)*" is "*string(abs(rand(Int))))
    else
        sim = run_selfSynch(psi,abs(rand(Int)))
        save_datal("/localdisk/luigi/data/selfSynch/system1/data_vs_S_Tk_part" * string(i) * ".jld2",sim)
    end
end

# using Distributed

# Base.@ccallable function julia_main(ARGS::Vector{String})::Cint

#     Distributed.@everywhere push!(LOAD_PATH,"../src")
#     Distributed.@everywhere include("./load.jl")

#     # include("lib_parallel.jl")
#     @everywhere include("./parameters.jl")
#     sim = run_selfSynch(p)
#     save_datal("/localdisk/luigi/data/selfSynch/system1/data1.jld2",sim)
#     return 0
# end
