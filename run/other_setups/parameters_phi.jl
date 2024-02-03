# parameters
N = 100 # number of particles
U₁ = U₂ = 0.0
κ = 100. # decay rate of the cavity field (units of ω\_R)
S₁ = 4κ/sqrt(N)
S₂ = 4κ/sqrt(N)
Δₑ = κ/100
Δc = -κ
num_monte = 3

tspan = (0.0, 800.0)

temp = κ
p = System_p(U₁,U₂,S₁,S₂,Δₑ,Δc,κ,temp,N,tspan,num_monte)

step = 2π/31
ps = System1_p[]
for phi in range(-pi,stop=pi-step,step=step)
    S₂ = S₁ * exp(im*phi)
    p = System1_p(U₁,U₂,S₁,S₂,Δₑ,Δc,κ,temp,N,tspan,num_monte)
    push!(ps,p)
end
