include("../src/selforg_core.jl")
include("../src/custom_functions.jl")

sim_index = parse(Int, ARGS[1])

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
    :Ekin => Ekin,
    :adaga => adaga,
    :Sz => Sz,
    :kurt => kurt,
    :Cos2 => Cos2,
);


### PARAMETERS ###
N = parse(Int, ARGS[2]); ### ATOM NUMBER (read from console)

kappa = 100
omegauv = 5:2:60; ### LIST OF ATOMIC FREQUENCIES
Deltac = -kappa;
Ntraj = 1 ### TRAJECTORY NUMBER
tlist = (0.0, 10.0) ### SIM TIME
deltaD = 10; ### DOPPLER WIDTH ~ temperature 
gvrange = 5:2:70; ### UNNORMALISED COUPLING

gv = sqrt.(-1 * (kappa^2 + Deltac^2) / (2 * Deltac) / N) * sqrt.(gvrange); ### RESCALE COUPLING STRENGTH

### PARAMETER ARRAY
p_array = [
    System_p(0.0, 0.0, g, g, omegau, Deltac, kappa, deltaD/2, N, tlist, Ntraj)
    for omegau in omegauv for g in gv
]


sim = extract_solution(many_trajectory_solver(p_array, saveat=0.05, seed=abs(rand(Int)), dt=1e-5));

JLD2.jldopen("short_time_sim_N=$(N)_i=$(sim_index).jld2", "w") do file
        write(file, "solution", sim)
    end

# sorted_sim = split_sim_from_par(sim);


# ### DATA MANIPULATION & STORAGE

# simsdf = [sim2df(sorted_sim[i],dict) for i in 1:length(sorted_sim)];
# JLD2.jldopen("short_time_sim_dfs_N=$(N).jld2", "w") do file
#     write(file, "solution", simsdf)
# end
