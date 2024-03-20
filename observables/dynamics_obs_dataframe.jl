using CSV
include("../src/selforg_core.jl")

function sim2df(sim::Array{Sol,1}, obs_dict::Dict)
    value_names = keys(obs_dict) 
    dfs = DataFrame[]

    dict = Dict{Symbol,Any}()
    for (o_, o) in obs_dict
        dict[o_] = expect(o,sim)
    end

    push!(dfs,DataFrame(dict))

    vcat(dfs...)
end 

# X = ∑ⱼ(σₓʲ cos(xⱼ))^2/N
function traj_X2(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [mean((sol.u[2N+1:3N,j].*cos.(sol.u[1:N,j])).^2) for j=1:length(sol.t)]
end
X2 = Observable(traj_X2,L"N^{-1}\sum_j (\sigma_x^j \cos(x_j))^2",L"X",L"order parameter $X^2$")

dict= Dict(
    :X2 => X2,
    :Sz => Sz,
    :Ekin => Ekin,
    :adaga => adaga,
)

U₁ = 0.0
U₂ = 0.0
g = 10
S₁ = g
S₂ = g

κ = 200
omegau = κ/10
Δc = -κ
deltaD = 10
N = 100
N_MC = 100

tspan = (0.0, 10.0)

p = System_p(U₁, U₂,S₁, S₂, omegau, Δc, κ, deltaD/2, N, tspan,N_MC)

sim = many_trajectory_solver(p,saveat=0.02, seed=abs(rand(Int)))

df = sim2df(sim, dict)

# Save DataFrame to CSV file with parameters in the filename
filename = "observables_g=$(g)_omegau=$(omegau)_kappa=$(κ)_deltaD=$(deltaD).csv"
CSV.write(filename, df)