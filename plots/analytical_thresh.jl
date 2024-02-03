using QuadGK
using PyPlot


# Define parameters
#selected_deltaDv = [5, 10, 20, 100]
omegauv = range(0, stop=100, length=201)
delta_Dvv = range(1, stop=100, length=201)
selected_omegau = [1,10,100]

gc_values = Vector{Vector{Float64}}(undef, length(selected_omegau))

opacities = range(1.0, stop=0.2, length=length(selected_omegau)) # Vary opacities

colors = ["#edf8b1", "#7fcdbb", "#2c7fb8"]

# Store points where omegau * delta_Dv = 1
# intersection_points = Dict{Int, Float64}()

for (idx, omegau) in enumerate(selected_omegau)
    gc = zeros(Float64, length(omegauv))

    for (bb, deltaDv) in enumerate(delta_Dvv)
        N = 100
        kappa = 100
        Deltac = 100
        omegar = 1 / 100 * deltaDv / (2 * π)

        fun(t) = exp(-deltaDv^2 * t^2 / 2) * (sin(omegau * t) + cos(omegau * t) * omegar * t)
        
        # Numerical integration using QuadGK
        Int, _ = quadgk(fun, 0, Inf;order=9)

        gc[bb] = sqrt((kappa^2 + Deltac^2 / (2 * Deltac) / Int) / N)
    end

    gc_values[idx] = gc

end

# Create the plot with varying opacities
figure()

for (idx, omegauv) in enumerate(selected_omegau)
    plot(delta_Dvv, gc_values[idx], label="omegau=$omegauv", linewidth=3, color=colors[idx], alpha=0.7)
end


xlabel(L"Initial ensemble temperature $\delta_D$", fontsize=16)
#xlabel(L"Atomic frequency $\omega_\uparrow$", fontsize=20)
ylabel(L"Coupling strength $g_c$", fontsize=16)
grid(false)
legend(loc="best", fontsize=14)
ylim(9, 40)
fontsize = 12
xticks(fontsize=fontsize)
yticks(fontsize=fontsize)
gca().set_position([0.125, 0.125, 0.775, 0.775])  # Adjust the position of the axis as needed

# Set the DPI for the figure
PyPlot.gcf().set_size_inches(8, 6)  # Adjust the width and height as needed
# savefig("threshold_Doppler_inset.png", dpi=600)  # Save the plot with x DPI
show()
