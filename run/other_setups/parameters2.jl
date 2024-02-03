# parameters
N = 100 # number of particles
U₁ = U₂ = 0.0
κ = 100. # decay rate of the cavity field (units of ω\_R)
S₁ = 5κ/sqrt(N)
S₂ = 5κ/sqrt(N)*exp(im*π/3)
Δₑ = 0.0
Δc = -κ
num_monte = 3

tspan = (0.0, 400.0)

temp = 1κ
p = System_p(U₁,U₂,S₁,S₂,Δₑ,Δc,κ,temp,N,tspan,num_monte)

ps = System_p[]
for S₁ in range(κ/sqrt(N)/2,stop=5κ/sqrt(N),length=30)
    S₂ = S₁
    p = System1_p(U₁,U₂,S₁,S₂,Δₑ,Δc,κ,temp,N,tspan,num_monte)
    push!(ps,p)
end
