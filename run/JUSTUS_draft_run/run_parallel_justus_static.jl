using ClusterManagers
using Random

include("../load.jl")
include("../parameters_static.jl")

np = parse(Int, ARGS[1])
Sre = parse(Int, ARGS[2])
Temp = parse(Int, ARGS[3])
Delta_e = parse(Int, ARGS[4])
Nat = parse(Int, ARGS[5])




sim = many_trajectory_solver(p,saveat=10.0, seed=abs(rand(Int)), maxiters=Int(1.0e9))
#save_datal("full_sims_traj1000/sim_data_static_pump(S=$(Sre), Delta_e=$(Delta_e), temp=$(Temp), t=800, $np).jld2", sim)
save_datal("full_sims_traj1000/pm_init_spin_badcav(S=$(Sre), Nat=$Nat, Delta_e=$(Delta_e), temp=$(Temp), t=800, $np).jld2", sim)
