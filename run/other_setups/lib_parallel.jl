using Distributed

@everywhere push!(LOAD_PATH,"../src")
@everywhere include("./load.jl")

# @time sim = solve(monte_prob,SOSRA(),num_monte=100,saveat=1.,abstol=1e-3,reltol=1e-3,parallel_type=:pmap,pmap_batch_size=1)
