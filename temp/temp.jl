using Interpolations
using JLD2
using Plots
using LaTeXStrings

include("../src/custom_functions.jl")
include("../src/selforg_core.jl")
# script to generate the numerical value of the critical coupling in the normal-SR phase transition of the main script   

#FIND NEW FUNCTION AS "plot_interpolate_threshold"


# include("../src/selforg_core.jl")
# include("../src/sim_reader.jl")


# sim_array = load_datalll("full_sims_traj1000/") 
# #would need tweaking when additional temperatures are in the folder, this works for only one temperatur now.
# sorted_sims = split_sim_from_par(sim_array)

# S = Float64[]
# y = (Float64[], Float64[])

# for sim in sorted_sims
#     push!(S, abs(sim[1].p.S₁)) #extract coupling strengths from sorted_sims

#     m, s, q = expect(adaga, sim) #extract mean, standard dev and 90 quantile
#     push!(y[1], m[end])
#     push!(y[2], s[end])
# end 


# #Need to make S into an equally spaced Tuple
# x = S[1]:1.0:S[end]

# itp_linear = linear_interpolation(x, y[1])
# itp_cubic = cubic_spline_interpolation(x, y[1])
# # Interpolation functions
# f_linear(x) = itp_linear(x)
# f_cubic(x) = itp_cubic(x)
# # Plots
# width, height = 1500, 800 

# x_new = S[1]:0.1:S[end] # smoother interval, necessary for cubic spline

# scatter(S, y[1], markersize=10,label="Data points")
# plot!(f_linear, x_new, w=3,label="Linear interpolation")
# plot!(f_cubic, x_new, linestyle=:dash, w=3, label="Cubic Spline interpolation")

# xlabel!(L"S", fontsize=14)
# ylabel!(L"⟨a^\dagger a\rangle", fontsize=14)

# plot!(size = (width, height))
# plot!(legend = :topleft, legendfontsize=15)


# function find_interp_crossover(sim::Array{Sol, 1})

#     sorted_sims = split_sim_from_par(sim)
    
#     S = Float64[]
#     y = (Float64[], Float64[])
    
#     for sorted_sim in sorted_sims
#         push!(S, abs(sorted_sim[1].p.S₁)) #extract coupling strengths from sorted_sims
    
#         m, s, q = expect(adaga, sorted_sim) #extract mean, standard dev and 90 quantile
#         push!(y[1], m[end])
#         push!(y[2], s[end])
#     end 
#     #Need to make S into an equally spaced Tuple
#     x = S[1]:1.0:S[end]

#     itp_cubic = cubic_spline_interpolation(x, y[1])
#     # Interpolation function
#     f_cubic(x) = itp_cubic(x)

#     x_new = S[1]:0.1:S[end] # smoother interval, necessary for cubic spline
    
#     # Calculate second derivative using Central Finite Difference
#     delta_x = x_new[2] - x_new[1]
#     f_cubic_second_derivative(x) = (f_cubic(x - delta_x) - 2 * f_cubic(x) + f_cubic(x + delta_x)) / delta_x^2
    
#     inset_x_range = 35:0.01:37.5 # !!!!!need to fix and automise this range
#     inset_y_values = f_cubic_second_derivative.(inset_x_range)


#     closest_to_zero_index = argmin(abs.(inset_y_values))
#     closest_to_zero_x = inset_x_range[closest_to_zero_index]
#     closest_to_zero_y = inset_y_values[closest_to_zero_index]
#     #returns interpolated numerical approximation of the second derivative crossover point
#     #closest_to_zero_y should be 0 in an ideal case
#     return closest_to_zero_x, closest_to_zero_y
# end

sort_save_datal("full_sims_traj1000", 1, [5,10,20])


##Delta_e=1

sim_array_5d1 = load_data_deltatemp("full_sims_traj1000", 1, 5)

sim_array_10d1 = load_data_deltatemp("full_sims_traj1000", 1, 10)

sim_array_20d1 = load_data_deltatemp("full_sims_traj1000", 1, 20)

#Delta_e=10

sim_array_5_d10 = load_data_deltatemp_alt("full_sims_traj1000", 5)

sim_array_10_d10 = load_data_deltatemp_alt("full_sims_traj1000", 10)

sim_array_20_d10 = load_data_deltatemp_alt("full_sims_traj1000", 20)

