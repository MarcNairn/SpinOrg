# parameters
N = 100 # number of particles
U₁ = U₂ = 0.0
κ = 100. # decay rate of the cavity field (units of ω\_R)
S₁ = 4κ/sqrt(N)
S₂ = 4κ/sqrt(N)
Δₑ = κ/100
Δc = -κ
num_monte = 1

tspan = (0.0, 800.0)

temp = 10

p = System_p(U₁,U₂,S₁,S₂,Δₑ,Δc,κ,temp,N,tspan,num_monte)

ps = System_p[]
for S₁ in range(2κ/sqrt(N),stop=8κ/sqrt(N),length=33)
    S₂ = S₁
    for Δₑ in range(κ/1000,stop=κ/2,length=33)
        p = System_p(U₁,U₂,S₁,S₂,Δₑ,Δc,κ,temp,N,tspan,num_monte)
        push!(ps,p)
    end
end
