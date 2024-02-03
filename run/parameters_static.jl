using DifferentialEquations
include("../src/selforg_core.jl")

#Read pump strength from parameters

temp = parse(Int, ARGS[3]) #same goes for here
Δₑ = parse(Int, ARGS[4])
N = parse(Int, ARGS[5]) # number of particles
g = parse(Int, ARGS[2])*sqrt(100/N)*(1/100*10) #in the case kappa=!100 or N=!100 need to renormalise the coupling, rescale with the standard N=100.


# parameters
#N = 100 
U₁ = U₂ = 0.0
κ = 10. # decay rate of the cavity field (units of ω\_R)
S₁ = g 
S₂ = g
Δc = -κ
num_monte = 1
#recall per trajectory need around 20 minutes

tspan = (0.0, 5000.0)

p = System_p(U₁,U₂,S₁,S₂,Δₑ,Δc,κ,temp,N,tspan,num_monte)
