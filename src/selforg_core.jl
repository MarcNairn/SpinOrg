    # THIS SCRIPT SOLVES THE SDE OF THE MAIN TEXT THROUGH A MONTE CARLO SAMPLING TECHNIQUE
    #
    # for more info look up documentation on 
    # https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble


    #REMOVED MODULE CALLS FOR NOW TO SIMPLIFY ANALYSIS
    
    #module Selforg 

    #using .CustomFunctions # .imports module without running 

	include("custom_functions.jl")

    using Random
    using Distributed
    using Statistics
    using LaTeXStrings
    using Printf
    using DataFrames
    using DifferentialEquations 

    export f_det, f_noise
    export Ekin,X,Y,Sx,Sy,Sz,Sx2,Sy2,Sz2,Cos2,Z,adaga,ar,absar,ai,absai,Cos
    export absX,absY,absSx,absSy,absSz,absZ,adiabaticar,adiabaticadaga
    export System_p
    export Observable, observable_list, observable_dict
    export expect, expect_full
    export generate_spin,initial_conditions
    export define_prob_from_parameters,regenerate_prob,get_last_steps
    export join_trajectories,sim2df
    export many_trajectory_solver
    export categorize_traj, split_sim, split_simPhi
    export split_sim_from_par, sort_sim_from_S


    #define parameter struct
    struct System_p
        U₁::Float64 
        U₂::Float64
        S₁::Complex{Float64}
        S₂::Complex{Float64}
        Δₑ::Float64
        Δc::Float64
        κ::Float64
        temp::Float64
        N::Int
        tspan::Tuple{Float64,Float64}
        N_MC::Int
    end

    ##############################################################################
    ################################SYSTEM FUNCTIONS##############################
    ##############################################################################
    function f_det(du,u,p,t)
    # t is in unit of w_R

    # u = [ xⱼ,     pⱼ,      σˣⱼ,       σʸⱼ,      σᶻⱼ,    aᵣ,    aᵢ]
    #     [1:N, N+1:2N, 2N+1..3N, 3N+1..4N, 4N+1..5N, 5N+1, 5N+2]
    # x_j in units of 1/kc (x'_j = kc x_j)
    # p_j in units of hbar*kc (p'_j = p_j/(hbar*kc))

    # ancilla variables
    N::Int = p.N
    REP::Float64 = real(p.S₁) + real(p.S₂) # ℜ(S₁) + ℜ(S₂)
    REM::Float64 = real(p.S₁) - real(p.S₂) # ℜ(S₁) - ℜ(S₂)
    IMP::Float64 = imag(p.S₁) + imag(p.S₂) # ℑ(S₁) + ℑ(S₂)
    IMM::Float64 = imag(p.S₁) - imag(p.S₂) # ℑ(S₁) - ℑ(S₂)
    aa::Float64 = (u[5N+1]^2 + u[5N+2]^2 - 0.5) # aᵣ² + aᵢ² - 1/2
    bb::Float64 = 2p.Δc
    cc::Float64 = 0.0
    dd::Float64 = 0.0

    for j = 1:N
        # ancilla variables
        sinuj::Float64, cosuj::Float64 = sincos(u[j])
        bb -= ((1-u[4N+j])p.U₁ + (1+u[4N+j])p.U₂)cosuj^2
        cc -= ( IMM*u[2N+j] + REM*u[3N+j] )cosuj
        dd += ( REP*u[2N+j] - IMP*u[3N+j] )cosuj

        # positions x_j
        du[j] = 2u[N+j]
        # momenta p_j
        du[N+j] = sinuj * ( cosuj*( (1-u[4N+j])p.U₁ + (1+u[4N+j])p.U₂ ) * aa + ( REP*u[2N+j]u[5N+1] + REM*u[3N+j]u[5N+2] + IMM*u[2N+j]u[5N+2] - IMP*u[3N+j]u[5N+1] ) )
        # sigmax_j
        du[2N+j] = -u[3N+j]*(p.Δₑ-(p.U₁-p.U₂)*cosuj^2*aa) - 2cosuj*u[4N+j] * (IMP*u[5N+1]-REM*u[5N+2])
        # sigmay_j
        du[3N+j] = u[2N+j]*(p.Δₑ-(p.U₁-p.U₂)*cosuj^2*aa) - 2cosuj*u[4N+j] * (REP*u[5N+1]+IMM*u[5N+2])
        # sigmaz_j
        du[4N+j] = 2cosuj*(IMP*u[2N+j]u[5N+1] + IMM*u[3N+j]u[5N+2] - REM*u[2N+j]u[5N+2] + REP*u[3N+j]u[5N+1])
    end

    # a_r
    du[5N+1] = bb/2 * u[5N+2] + cc/2 - p.κ*u[5N+1]
    # a_i
    du[5N+2] = -bb/2 * u[5N+1] + dd/2 - p.κ*u[5N+2]

end

    ###Stochastic part of the SDE (cavity noise alone)####

    #@everywhere 
    function f_noise(du,u,p,t)
        N::Int = p.N
        du[1:5N] .= 0.0
        du[5N+1] = sqrt(2*p.κ)
        du[5N+2] = sqrt(2*p.κ)
    end

###INITIAL CONDITIONS###

const σ₁ = Complex[0. 1.; 1. 0.]
const σ₂ = Complex[0. -im; im 0.]
const σ₃ = Complex[1. 0.; 0. -1.]
    
    function generate_spin(θ, ϕ)
        v = [cos(θ/2) sin(θ/2)*cis(ϕ)]
        real(v*σ₁*v')[], real(v*σ₂*v')[], real(v*σ₃*v')[]
    end

    function initial_conditions(p::System_p, seed=abs(rand(Int)))
        N::Int = p.N
        Random.seed!(seed) # random number generator
        u0 = zeros(5N + 2)

        u0[1:N] = 2pi.*rand(N) # generate random positions
        u0[N+1:2N] = p.temp .* randn(N) # generate random momenta, temp refers to the square root of a temperature unit

        u0[2N+1:4N] = 2bitrand(2N) .- 1 # σˣⱼ and σʸⱼ are 1 or -1 
        u0[4N+1:5N] .= -1. # σᶻⱼ = -1, atoms in the ground state

        # #evaluate spin expectation values as initial conditions for σx,σy,σz. 
        # #i.e. we define σᵢ≡<ψ(0)|̂σᵢ|ψ(0)> where ̂σᵢ is the i-th Pauli matrix 
        # for j = 1:N
        #     u0[2N+j], u0[3N+j], u0[4N+j] = generate_spin(pi+rand()*0.1pi,2pi*rand())
        # end

        u0[5N+1:end] .= 0.0 # cavity empty
        return u0
    end

    #@everywhere 
    function define_prob_from_parameters(p::System_p,seed=abs(rand(Int)))
        # initial conditions
        Random.seed!(seed) # random number generator
        u0 = initial_conditions(p,abs(rand(Int)))
        u0_arr = [initial_conditions(p,abs(rand(Int))) for j=1:p.N_MC]

    # NEED TO DEFINE THE FUNCTION TO REDRAW INITIAL CONDITIONS AT EVERY TRAJECTORY
        function prob_func(prob,i,repeat)
            # @. prob.u0 = initial_conditions(N,κ,rng)
            @. prob.u0 = u0_arr[i]
        prob
        end

        prob = SDEProblem(f_det,f_noise,u0,p.tspan,p)
        monte_prob = EnsembleProblem(prob, prob_func = prob_func)

        return prob, monte_prob
    end

     function define_prob_from_parameters(ps::Array{System_p,1}, seed=abs(rand(Int)))
        # initial conditions
        Random.seed!(seed) # random number generator
    
        p_arr = [p for p in ps for i in 1:p.N_MC]
        u0 = initial_conditions(ps[1],abs(rand(Int)))
        u0_arr = [initial_conditions(p,abs(rand(Int))) for p in p_arr]
    
        function prob_func(prob,i,repeat)
            SDEProblem(f_det,f_noise,u0_arr[i],p_arr[i].tspan,p_arr[i])
        end
    
        prob = SDEProblem(f_det,f_noise,u0,ps[1].tspan,ps[1])
        monte_prob = EnsembleProblem(prob, prob_func = prob_func)
    
        return prob, monte_prob
    end
    



    function many_trajectory_solver(p::System_p;saveat::Float64=10.0,seed::Int=abs(rand(Int)),maxiters::Int=Int(1e9))#, reltol::Float64=1e-4)
        prob, monte_prob = define_prob_from_parameters(p,seed)
        #print("calculating $(p.N_MC) trajectories on $(gethostname()) with $(nworkers()) workers..")
        elt = @elapsed sim = solve(monte_prob::EnsembleProblem, EM(), EnsembleThreads(), trajectories=p.N_MC, saveat=saveat, maxiters=maxiters, progress=true)

        # EnsembleDistributed() recommended here when each trajectory is not very quick (like here)
        println("done in $elt seconds.")
        return sim
    end
    


    function many_trajectory_solver(ps::Array{System_p,1};saveat::Float64=10.0,seed::Int=abs(rand(Int)),maxiters::Int=Int(1e12), dt::Float64)
        prob, monte_prob = define_prob_from_parameters(ps,seed)
        N_MC = sum([p.N_MC for p in ps])
       #print("calculating $N_MC trajectories on $(gethostname()) with $(nworkers()) workers..")
        elt = @elapsed sim = solve(monte_prob::EnsembleProblem, EM(), EnsembleThreads(), trajectories=N_MC, saveat=saveat, maxiters=maxiters; dt=dt)
        println("done in $elt seconds.")
        return sim
    end



    

function get_last_steps(sol::Sol)
    Sol(sol.u[:,end-1:end],sol.p,sol.t[end-1:end],sol.alg)
end

function get_last_steps(sim::Array{Sol,1})
    [get_last_steps(sol) for sol in sim]
end

function regenerate_prob(sol::Sol,deltat::Real)
    u0 = sol.u[:,end]
    p = sol.p
    tspan = (sol.t[end], sol.t[end]+deltat)
    prob = SDEProblem(f,f_noise,u0,tspan,p)
end

function regenerate_prob(sim::Array{Sol,1},deltat::Real)
    prob =  regenerate_prob(sim[1],deltat)
    function prob_func(prob,i,repeat)
        regenerate_prob(sim[i],deltat)
    end

    monte_prob = EnsembleProblem(prob,prob_func=prob_func)
end

function stringfromp(p)
    U₁,U₂,S₁,S₂,Δₑ,Δc,κ,N = parsfromp(p)
    str1 = @sprintf("U1%.4f_U2%.4f_S1%.4f%+.4fim_S2%.4f%+.4fim_De%.4f_Dc%.4f_k%.4f_N%d", U₁,U₂,real(S₁),imag(S₁),real(S₂),imag(S₂),Δₑ,Δc,κ,N)
end

function pfromstring(str::String)
    x = split(str,'_')
    U₁ = parse(Float64,x[1][3:end])
    U₂ = parse(Float64,x[2][3:end])
    S₁ = parse(Complex{Float64},x[3][3:end])
    S₂ = parse(Complex{Float64},x[4][3:end])
    Δₑ = parse(Float64,x[5][3:end])
    Δc = parse(Float64,x[6][3:end])
    κ = parse(Float64,x[7][2:end])
    N = parse(Int,x[8][2:end])
    return U₁,U₂,S₁,S₂,Δₑ,Δc,κ,N
end

function join_trajectories(sim::Array{Sol,1})
    idx = size(sim[1].t)[1]
    join_trajectories(sim,idx)
end

function join_trajectories(sim::Array{Sol,1},idx)
    N::Int = try sim[1].p.N catch; sim[1].p[10] end # for backward compatibility

    ntraj = size(sim)[1]

    u0 = zeros(ntraj*size(sim[1].u)[1])

    for i in 1:ntraj
        u0[range((i-1)*N+1,length=N)] = sim[i].u[1:N,idx]
        u0[range((ntraj+i-1)*N+1,length=N)] = sim[i].u[N+1:2N,idx]
        u0[range((2ntraj+i-1)*N+1,length=N)] = sim[i].u[2N+1:3N,idx]
        u0[range((3ntraj+i-1)*N+1,length=N)] = sim[i].u[3N+1:4N,idx]
        u0[range((4ntraj+i-1)*N+1,length=N)] = sim[i].u[4N+1:5N,idx]
        u0[range(5ntraj*N+(i-1)*2+1,length=2)] = sim[i].u[5N+1:5N+2,idx]
    end
    u0
end

# function sim2df(sim::Array{Sol,1})
#     nvars = size(sim[1].u)[1]
#     N::Int = (nvars-2)/5

#     value_names = Symbol[]
#     for var in ["x_","p_","sx_","sy_","sz_"]
#         for i = 1:N
#             push!(value_names,Symbol(var*"$i"))
#         end
#     end
#     push!(value_names, :a_r)
#     push!(value_names, :a_i)
#     # push!(value_names, :timestamp)
#     # push!(value_names, :traj)

#     dfs = DataFrame[]
#     for (k,sol) in enumerate(sim)
#         dict = Dict{Symbol,Any}(value_names[i] => sol.u[i,:] for i = 1:nvars)
#         dict[:a_r] = sol.u[5N+1,:]
#         dict[:a_i] = sol.u[5N+2,:]
#         for (o_, o) in observable_dict
#             dict[o_] = expect(o,sol)
#         end
#         dict[:timestamp] = sol.t
#         dict[:traj] = k
#         push!(dfs,DataFrame(dict))
#     end
#     vcat(dfs...)
# end
######################################################################

############################## OBSERVABLES ###########################

######################################################################

Base.@kwdef struct Observable
    s_traj
    formula::String
    short_name::String = formula
    name::String = short_name
end

Observable(s_traj,formula,short_name) = Observable(s_traj=s_traj,formula=formula,short_name=short_name)
Observable(s_traj,formula) = Observable(s_traj=s_traj,formula=formula)

function expect(o::Observable,sol::RODESolution)
    sol_ = extract_solution(sol)[1]
    expect(o,sol_)
end

function expect(o::Observable,sim::EnsembleSolution)
    sim_ = extract_solution(sim)
    expect(o,sim_)
end

function expect(o::Observable,sol::Sol)
    o.s_traj(sol)
end

function expect_full(o::Observable,sim::Array{Sol,1})
    [expect(o,sim[j]) for j=1:length(sim)]
end

function expect(o::Observable,sim::Array{Sol,1})
    Os = expect_full(o,sim)
    Omean = mean(Os)
    Ostd = stdm(Os,Omean)
    bb = hcat(Os...)
    Oq90 = [quantile(bb[i,:],[0.05,0.95]) for i in 1:size(bb)[1]]

    return Omean, Ostd, Oq90
end

# kinetic energy/N
function traj_kinetic_energy(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end # for backward compatibility
    [mean(sol.u[N+1:2N,j].^2) for j=1:length(sol.t)]
end
Ekin = Observable(traj_kinetic_energy,L"\sum_j p_j^2",L"E_\mathrm{kin}","kinetic energy")

# kurtosis of momentum distribution
function traj_kurt(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end # for backward compatibility
    [mean(sol.u[N+1:2N,j].^4)./mean(sol.u[N+1:2N,j].^2).^2 for j=1:length(sol.t)]
end
kurt = Observable(traj_kurt,L"\sum_j p_j^4 / (p^2)^2",L"\mathrm{kurtosis}","kurtosis")

# X = ∑ⱼσₓʲ cos(xⱼ)/N
function traj_X(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [mean(sol.u[2N+1:3N,j].*cos.(sol.u[1:N,j])) for j=1:length(sol.t)]
end
X = Observable(traj_X,L"N^{-1}\sum_j \sigma_x^j \cos(x_j)",L"X",L"order parameter $X$")

# absX = |∑ⱼσₓʲ cos(xⱼ)|/N
function traj_absX(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [abs(mean(sol.u[2N+1:3N,j].*cos.(sol.u[1:N,j]))) for j=1:length(sol.t)]
end
absX = Observable(traj_absX,L"N^{-1}\abs{\sum_j \sigma_x^j \cos(x_j)}",
                  L"\abs{X}",L"order parameter $\abs{X}$")

# Y = ∑ⱼσ\_yʲ cos(xⱼ)/N
function traj_Y(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [mean(sol.u[3N+1:4N,j].*cos.(sol.u[1:N,j])) for j=1:length(sol.t)]
end
Y = Observable(traj_Y,L"N^{-1}\sum_j \sigma_y^j \cos(x_j)",L"Y",L"order parameter $Y$")

# absY = |∑ⱼσ\_yʲ cos(xⱼ)|/N
function traj_absY(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [abs(mean(sol.u[3N+1:4N,j].*cos.(sol.u[1:N,j]))) for j=1:length(sol.t)]
end
absY = Observable(traj_absY,L"N^{-1}\abs{\sum_j \sigma_y^j \cos(x_j)}",
                  L"\abs{Y}",L"order parameter $\abs{Y}$")

# Z =  ∑ⱼσ\_zʲ⋅cos(xⱼ)/N
function traj_Z(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [mean(sol.u[4N+1:5N,j].*cos.(sol.u[1:N,j])) for j=1:length(sol.t)]
end
Z = Observable(traj_Z,L"N^{-1}\sum_j\sigma_z\cos(x_j)",L"Z")

# |Z| =  |∑ⱼσ\_zʲ⋅cos(xⱼ)|/N
function traj_absZ(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [abs(mean(sol.u[4N+1:5N,j].*cos.(sol.u[1:N,j]))) for j=1:length(sol.t)]
end
absZ = Observable(traj_absZ,L"N^{-1}\abs{\sum_j\sigma_z\cos(x_j)}",
                       L"\abs{Z}")

# sigmax = ∑ⱼσₓʲ/N
function traj_sigmax(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [mean(sol.u[2N+1:3N,j]) for j=1:length(sol.t)]
end
Sx = Observable(traj_sigmax,L"N^{-1}\sum_j \sigma_x^j",L"$\langle \hat J_x \rangle$",L"\sigma_x")

# absSx = |∑ⱼσₓʲ|/N
function traj_absSx(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [abs(mean(sol.u[2N+1:3N,j])) for j=1:length(sol.t)]
end
absSx = Observable(traj_absSx,L"N^{-1}\abs{\sum_j \sigma_x^j}",L"\abs{\sigma_x}")

# sigmay = ∑ⱼσ\_yʲ/N
function traj_sigmay(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [mean(sol.u[3N+1:4N,j]) for j=1:length(sol.t)]
end
Sy = Observable(traj_sigmay,L"N^{-1}\sum_j \sigma_y^j",L"$\langle \hat J_y \rangle$",L"\sigma_y")

# absSy = |∑ⱼσ\_yʲ|/N
function traj_absSy(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [abs(mean(sol.u[3N+1:4N,j])) for j=1:length(sol.t)]
end
absSy = Observable(traj_absSy,L"N^{-1}\abs{\sum_j \sigma_y^j}",L"\abs{\sigma_y}")


# sigmaz = ∑ⱼσ\_zʲ/N
function traj_sigmaz(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [mean(sol.u[4N+1:5N,j]) for j=1:length(sol.t)]
end
Sz = Observable(traj_sigmaz,L"N^{-1}\sum_j \sigma_z^j",L"$\langle \hat J_z \rangle$",L"\sigma_z")

# absSz = |∑ⱼσ\_zʲ|/N
function traj_absSz(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [abs(mean(sol.u[4N+1:5N,j])) for j=1:length(sol.t)]
end
absSz = Observable(traj_absSz,L"N^{-1}\abs{\sum_j \sigma_z^j}",L"\abs{\sigma_z}")

# sigmax² = ∑ⱼ(σₓʲ)²/N
function traj_sigmax2(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [mean(sol.u[2N+1:3N,j].^2) for j=1:length(sol.t)]
end
Sx2 = Observable(traj_sigmax2,L"N^{-1}\sum_j (\sigma_x^j)^2",L"\sigma_x^2")

# sigmay² = ∑ⱼ(σ\_yʲ)²/N
function traj_sigmay2(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [mean(sol.u[3N+1:4N,j].^2) for j=1:length(sol.t)]
end
Sy2 = Observable(traj_sigmay2,L"N^{-1}\sum_j (\sigma_y^j)^2",L"\sigma_y^2")

# sigmaz² = ∑ⱼ(σ\_zʲ)²/N
function traj_sigmaz2(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [mean(sol.u[4N+1:5N,j].^2) for j=1:length(sol.t)]
end
Sz2 = Observable(traj_sigmaz2,L"N^{-1}\sum_j (\sigma_z^j)^2",L"\sigma_z^2")

# cos(x)^2 = ∑ⱼcos(xⱼ)²/N
function traj_cos2(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [mean(cos.(sol.u[1:N,j]).^2) for j=1:length(sol.t)]
end
Cos2 = Observable(traj_cos2,L"N^{-1}\sum_j\cos(x_j)^2",L"\mathcal{B}",L"\cos^2")

# cos(x) = ∑ⱼcos(xⱼ)/N
function traj_cos(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [mean(cos.(sol.u[1:N,j])) for j=1:length(sol.t)]
end
Cos = Observable(traj_cos,L"N^{-1}\sum_j\cos(x_j)",L"\cos")

# cavity population = aᵣ² + aᵢ² - 0.5
function traj_adaga(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [sum(sol.u[5N+1:5N+2,j].^2) for j=1:length(sol.t)]
end
adaga = Observable(traj_adaga,L"a_r^2+a_i^2-0.5",L"\langle a^\dag a\rangle","cavity population")

# real part cavity field aᵣ
function traj_ar(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [sol.u[5N+1,j] for j=1:length(sol.t)]
end
ar = Observable(traj_ar, L"\langle \frac{a+a^\dag}{2}\rangle", L"$\langle \hat{a}_{\mathrm{r}} \rangle$", "real part of cavity field")

# absolute value of real part cavity field aᵣ
function traj_absar(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [abs(sol.u[5N+1,j]) for j=1:length(sol.t)]
end
absar = Observable(traj_absar, L"\abs{a_r}", L"\langle \frac{\abs{a+a^\dag}}{2}\rangle", "absolute value of real part of cavity field")

# imaginary part cavity field aᵢ
function traj_ai(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [sol.u[5N+2,j] for j=1:length(sol.t)]
end
ai = Observable(traj_ai, L"\langle i\frac{a-a^\dag}{2}\rangle", L"$\langle \hat{a}_{\mathrm{r}} \rangle$", "imaginary part of cavity field")

# absolute value of imaginary part cavity field aᵢ
function traj_absai(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [abs(sol.u[5N+2,j]) for j=1:length(sol.t)]
end
absai = Observable(traj_absai, L"\abs{a_i}", L"\langle \frac{\abs{a-a^\dag}}{2}\rangle", "absolute value of imaginary part of cavity field")

# adiabatic cavity field aᵣ
function traj_adiabaticar(sol::Sol)
    N::Int = sol.p.N
    S₁::Complex{Float64} = sol.p.S₁
    S₂::Complex{Float64} = sol.p.S₂
    U₁::Float64 = sol.p.U₁
    U₂::Float64 = sol.p.U₂
    Δc::Float64 = sol.p.Δc
    κ::Float64 = sol.p.κ
    [real(sum((S₁.*(sol.u[2N+1:3N,j]-sol.u[3N+1:4N,j]*im)./2 + S₂.*(sol.u[2N+1:3N,j]+sol.u[3N+1:4N,j]*im)./2).*cos.(sol.u[1:N,j]))/(Δc*im-κ)) for j=1:length(sol.t)]
end
adiabaticar = Observable(traj_adiabaticar, L"\tilde{a}", L"\langle \tilde{a_r} \rangle", "adiabatic value of the real part of the cavity field")

# adiabatic cavity population
function traj_adiabaticadaga(sol::Sol)
    N::Int = sol.p.N
    S₁::Complex{Float64} = sol.p.S₁
    S₂::Complex{Float64} = sol.p.S₂
    U₁::Float64 = sol.p.U₁
    U₂::Float64 = sol.p.U₂
    Δc::Float64 = sol.p.Δc
    κ::Float64 = sol.p.κ
    [abs2(sum((S₁.*(sol.u[2N+1:3N,j]-sol.u[3N+1:4N,j]*im)./2 + S₂.*(sol.u[2N+1:3N,j]+sol.u[3N+1:4N,j]*im)./2).*cos.(sol.u[1:N,j]))/(Δc*im-κ)) for j=1:length(sol.t)]
end
adiabaticadaga = Observable(traj_adiabaticadaga, L"\tilde{a}", L"\langle \tilde{a^\dag a} \rangle", "adiabatic value of the cavity population")

# single_sigmax² = (σₓʲ)²
function traj_single_sx2(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [(sol.u[2N+50,j].^2) for j=1:length(sol.t)]
end
single_sx2 = Observable(traj_single_sx2,L"(\sigma_x)^2",L"\sigma_x^2")

function traj_single_sy2(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [(sol.u[3N+50,j].^2) for j=1:length(sol.t)]
end
single_sy2 = Observable(traj_single_sy2,L"(\sigma_y)^2",L"\sigma_y^2")

function traj_single_sz2(sol::Sol)
    N::Int = try sol.p.N catch; sol.p[10] end
    [(sol.u[4N+50,j].^2) for j=1:length(sol.t)]
end
single_sz2 = Observable(traj_single_sz2,L"(\sigma_z)^2",L"\sigma_z^2")


observable_list = [Ekin,kurt,X,Y,Sx,Sy,Sz,Cos2,Z,ar,adaga,absX,absY,absSx,absSy,absSz,
                   absZ, absar, adiabaticar, adiabaticadaga, single_sx2, single_sy2, single_sz2]

observable_dict = Dict( :Ekin => Ekin,
                        :kurt => kurt,
                        :X => X,
                        :Y => Y,
                        :Z => Z,
                        :Sx => Sx,
                        :Sy => Sy,
                        :Sz => Sz,
                        :Sx2 => Sx2,
                        :Sy2 => Sy2,
                        :Sz2 => Sz2,
                        :Cos2 => Cos2,
                        :adaga => adaga,
                        :ar => ar,
                        :ai => ai,
                        :absX => absX,
                        :absY => absY,
                        :absZ => absZ,
                        :absSx => absSx,
                        :absSy => absSy,
                        :absSz => absSz,
                        :absar => absar,
                        :absai => absai,
                        :adiabaticar => adiabaticar,
                        :adiabaticadaga => adiabaticadaga,
                        :single_sx2 => single_sx2,
                        :single_sy2 => single_sy2,
                        :single_sz2 => single_sz2,
                        )


                        

                     
####################################################################################
############################ EXTRACT SOLUTIONS FROM TRAJECTORIES ###################
####################################################################################


function categorize_traj(sol::Sol) #categorise the broken Z₂ symmetry X,Y states
    expect(X,sol)[end] >=0 ? x=1 : x=-1
    expect(Y,sol)[end] >=0 ? y=1 : y=-1
    return x,y
end

function split_sim(sim::Array{Sol,1})  
    a,b,c,d = Sol[],Sol[],Sol[],Sol[]
    # a = (-1,-1), b = (-1,1), c=(1,-1), d=(1,1)
    for sol in sim
        catsol = categorize_traj(sol)
        catsol[1]==-1 ? catsol[2]==-1 ? push!(a,sol) : push!(b,sol) :
            catsol[2]==-1 ? push!(c,sol) : push!(d,sol)
    end
    return a,b,c,d
end

function split_simPhi(sim::Array{Sol,1})
    a,b = Sol[],Sol[]
    # a = Φ>0, b = Φ<0
    for sol in sim
        ϕ = -1//2*angle(sol.p.S₂/sol.p.S₁)
        Φ = cos(ϕ)*expect(X,sol)[end] + sin(ϕ)*expect(Y,sol)[end]
        Φ>0 ? push!(a,sol) : push!(b,sol)
    end
    return a,b
end

function split_sim_from_par(sim::Array{Sol,1}, ret_pars::Bool=false)
    par_list = System_p[]
    categories = Array{Sol,1}[]

    for sol in sim
        if !(sol.p in par_list)
            push!(par_list,sol.p)
            push!(categories,sim[findall(x -> x.p == sol.p,sim)])
        end
    end
    if ret_pars
        return categories, par_list
    else
        return categories
    end
end

function sort_sim_from_S(sim::Array{Sol,1}, ret_pars::Bool=false)
    categories, par_list = split_sim_from_par(sim,true)
    sp = sortperm(par_list,lt=(x,y)->isless(abs(x.S₁)^2+abs(x.S₂)^2,abs(y.S₁)^2+abs(y.S₂)^2))
    if ret_pars
        return categories[sp], par_list[sp]
    else
        return categories[sp]
    end
end

#####################

#end
