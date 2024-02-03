using QuadGK
using PyPlot


# Define parameters
omegauv = range(0.1, stop=5, length=201)
delta_Dvv = range(0.1, stop=125, length=201)

gc_values = zeros(Float64, length(omegauv), length(delta_Dvv))

# Calculate the integral for all combinations of omegauv and delta_Dvv
for (i, omegau) in enumerate(omegauv)
    for (j, deltaDv) in enumerate(delta_Dvv)
        N = 100
        kappa = 100
        Deltac = 100
        omegar =  1/100 * deltaDv/(2*pi)

        fun(t) = exp(-deltaDv^2 * t^2 / 2) * (sin(omegau * t) + cos(omegau * t) * omegar * t)
        
        # Numerical integration using QuadGK
        Int, _ = quadgk(fun, 0, Inf; order=9)

        gc_values[i, j] = sqrt((kappa^2 + Deltac^2 / (2 * Deltac) / Int) / N)
    end
end

# Calculate values for the line \omega_u * 1 / \delta_D = 1
line_values = ones(length(omegauv)) ./ omegauv

# Create the heatmap
figure()
pcolormesh(delta_Dvv, omegauv, gc_values, cmap="inferno")
colorbar(label=L"Critical coupling strength $g_c$")
xlabel(L"Initial ensemble temperature $\delta_D$")
ylabel(L"Atomic frequency $\Delta_e$")
clim(10,70)
grid(false)

# Plot the line \omega_u * 1 / \delta_D = 1
plot(delta_Dvv, line_values, label=L"$\omega_\uparrow \cdot \frac{1}{\delta_D} = 1$", color="white", linestyle="-", linewidth=3)


# Restrict the y-axis range
ylim(0.1, 20)

