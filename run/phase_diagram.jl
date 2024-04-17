include("../src/selforg_core.jl")
include("../src/custom_functions.jl")

# X = ∑ⱼ(σₓʲ cos(xⱼ))^2/N^2
function traj_X2(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [mean((sol.u[2N+1:3N,j].*cos.(sol.u[1:N,j]))).^2 for j=1:length(sol.t)]
end
X2 = Observable(traj_X2,L"N^{-1}\sum_j (\sigma_x^j \cos(x_j))^2",L"X",L"order parameter $X^2$");


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
N = 200; ### ATOM NUMBER (READ FROM TERMINAL)

kappa = 100
omegauv = 5:2:60; ### LIST OF ATOMIC FREQUENCIES
Deltac = -kappa;
Ntraj = 100 ### TRAJECTORY NUMBER
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
sorted_sim = split_sim_from_par(sim);


