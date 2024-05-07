using JLD2
using DataFrames

include("../src/selforg_core.jl")
include("../src/custom_functions.jl")

N = parse(Int, ARGS[1])

# X = ∑ⱼ(σₓʲ cos(xⱼ))^2/N^2
function traj_X2(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [mean((sol.u[2N+1:3N,j].*cos.(sol.u[1:N,j]))).^2 for j=1:length(sol.t)]
end
X2 = Observable(traj_X2,L"N^{-1}\sum_j (\sigma_x^j \cos(x_j))^2",L"X",L"order parameter $X^2$");

### SIMS TO DATAFRAME
function sim2df(sim::Array{Sol,1}, obs_dict::Dict)
    value_names = keys(obs_dict) 
    dfs = DataFrame[]

    dict = Dict{Symbol,Any}()
    for (o_, o) in obs_dict
        dict[o_] = expect(o,sim)
    end
    #dict[:timestamp] = sim[1].t

    push!(dfs,DataFrame(dict))

    vcat(dfs...)
end;

### DEFINE OBSERVABLE DICTIONARY OF INTEREST ###
dict= Dict(
    :X2 => X2,
    :adaga => adaga,
    :Sz => Sz,
    :kurt => kurt,
    :Cos2 => Cos2,
);

### READ SIMULATIONS FROM DIRECTORY
sim =  load_data_pm_generic("ode_phasediagram_N=100_",99)
### SORT THROUGH PARAMETERS
sorted_sim =  split_sim_from_par(sim);
### CONVERT TO DATAFRAME ARRAY
simsdf = [sim2df(sorted_sim[i],dict) for i in 1:length(sorted_sim)]
### SAVE DATAFRAME ARRAY
JLD2.jldopen("short_time_ode_dfsN=$(N).jld2", "w") do file
    write(file, "solution", simsdf)
end