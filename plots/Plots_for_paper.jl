#HERE WE LIST ADDITIONAL PLOTS FOR THE SPIN-MOTION ORGANISATION PAPER IN PREPARATION 
#THAT ARE NOT IMMEDIATELY LISTED IN THE LUIGIPLOTTING MODULE 

using Interpolations
using JLD2
using DiffEqBase
using LaTeXStrings
using LinearAlgebra
using LsqFit
using Measurements
using Random: MersenneTwister
using PyPlot
using StatsBase
#using EasyFit #needed for exponential fitting
using Distributions

function plot_adaga(sim::Array{Sol,1})
    y,y_std,y_q90 = expect(adaga,sim)
    y_q90 = hcat(y_q90...)
    tlist = sim[1].t

    matplotlib[:rc]("axes", labelpad=0.5)

    fig, ax = subplots(1,1,figsize=[3.25, 2.])

    color="C2"

    ax[:set_ylabel](L"cavity population $|\alpha|^2/N_\mathrm{at}$")
    ax[:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
    # ax[:set_xscale]("log")
    ax[:plot](tlist.+1,y.*1/sim[1].p.N,color=color)
    ax[:fill_between](tlist.+1,y_q90[1,:]*1/sim[1].p.N,y_q90[2,:]*1/sim[1].p.N,color=color,alpha=0.2)
    #ax[:fill_between](tlist.+1,y.+y_std,y.-y_std,color=color,alpha=0.5)

    fig[:tight_layout]()

    return fig, ax
end

function plot_adaga(sim::Array{Sol,1},filename::String)
    fig, ax = plot_adaga(sim)
    fig[:savefig](filename)
end

function plot_adaga_derivative(sim::Array{Sol,1})
    y,y_std,y_q90 = expect(adaga,sim)
    y_q90 = hcat(y_q90...)
    tlist = sim[1].t
    fig, ax =subplots()
    left, bottom, width, height = [0.25, 0.6, 0.5, 0.2]
    axinset = fig.add_axes([left, bottom, width, height])
    
    t_range = tlist[1]:10.0:tlist[end]
    # Cubic Spline Interpolation
    itp_cubic_adaga = cubic_spline_interpolation(t_range, y./sim[1].p.N)

    # Interpolation functions
    f_cubic_adaga(x) = itp_cubic_adaga(x)

    # Create Smooth Time Vector
    delta_t = 0.01
    t_smooth = tlist[2]:delta_t:tlist[end-1]
    t_inset = Int(tlist[end]/2):delta_t:tlist[end-1]

    # Evaluate First Derivative with Interpolated Data (Central Finite Difference)

    f_cubic_first_derivative_adaga(x) = (f_cubic_adaga(x + delta_t) - f_cubic_adaga(x - delta_t)) / (2 * delta_t)
    mean_sq_diff(x) = (f_cubic_adaga(x) - f_cubic_adaga(x + delta_t))^2

    # Plotting First Derivatives
    ax[:plot](t_smooth, f_cubic_first_derivative_adaga.(t_smooth), linestyle="--", color="C3", label="∂ₜ|α|²/N")
    axinset[:plot](t_inset, mean_sq_diff.(t_inset)./delta_t, linestyle="-", color="grey", label="MSD")
    axinset[:set_xscale]("log")
    axinset[:legend]()
    ax[:legend](handlelength=2.5,loc="upper left",bbox_to_anchor=(0.2, 1.08),framealpha=1)  
    fig[:tight_layout]()
    ax[:set_xscale]("log")
    return fig, ax
end

function plot_adaga_vs_S(sim::Array{Sol, 1})
    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")

    categories, par_list = split_sim_from_par(sim, true)

    println(par_list)

    S = Float64[]
    y = (Float64[], Float64[])
    for x in categories
        push!(S, abs(x[1].p.S₁))

        m, s, q = expect(adaga, x)
        push!(y[1], m[end])
        push!(y[2], q[end])
    end

    matplotlib[:rc]("axes", labelpad=1)

    fig, ax = subplots(figsize=[6.2, 4.])

    ax[:set_ylabel](L"\langle a^\dag a\rangle")
    
    # Plot the line
    ax[:plot](S, y[1], color="C0")

    ax[:fill_between](S, y[1] + y[2], y[1] - y[2], color="C0", alpha=0.2)

    ax[:set_xlabel](L"{|S_1|}={|S_2|}")

    fig[:tight_layout](h_pad=0., w_pad=-0.)
    return fig, ax
end

function plot_adaga_vs_S(solorsim,filename::String)
    fig, ax = plot_adaga_vs_S(solorsim)
    fig[:savefig](filename)
end


##### FIX THIS; IMPLEMENT IN PYPLOT
function plot_ai_ar_scatter(sim::Array{Sol, 1}, bins=50) 
    sorted_sim = split_sim_from_par(sim)
    fig, ax = plt.subplots()
    N::Int = sorted_sim[1][1].p.N  # Extract atom number of simulations

    # Extracting values for scatter plot
    scatter_values_sr = hcat(
        [sorted_sim[end][traj].u[5N + 1, end] for traj in 1:length(sorted_sim[end])],
        [sorted_sim[end][traj].u[5N + 2, end] for traj in 1:length(sorted_sim[end])]
    )

    scatter_values_normal = hcat(
        [sorted_sim[1][traj].u[5N + 1, end] for traj in 1:length(sorted_sim[end])],
        [sorted_sim[1][traj].u[5N + 2, end] for traj in 1:length(sorted_sim[end])]
    )

    mean_x_normal = mean(scatter_values_normal[:, 1])
    mean_y_normal = mean(scatter_values_normal[:, 2])

    mean_x_sr = mean(scatter_values_sr[:, 1])
    mean_y_sr = mean(scatter_values_sr[:, 2])

    ax.scatter(scatter_values_sr[:, 1], scatter_values_sr[:, 2], color="blue", alpha=0.5, label="SR", s=25)
    ax.scatter(scatter_values_normal[:, 1], scatter_values_normal[:, 2], color="red", alpha=0.5, label="Normal", s=25)

    ax.set_label(L"a_r")#, fontsize=15)
    ax.set_ylabel(L"a_i")#, fontsize=15)

    ax.legend()

    #ax.xlim(mean_x_sr - 2*std(scatter_values_sr[:, 1]), mean_x_sr + 2*std(scatter_values_sr[:, 1]))
    #ax.ylim(mean_y_sr - 2*std(scatter_values_sr[:, 2]), mean_y_sr + 2*std(scatter_values_sr[:, 2]))

end


function plot_ordering_vs_S(sims::Array{Sol,1}...) #can include multiple sims

    colorlist = ["C1","C2","C3","C4","C5","C6"]
    linelist = ["-","--",":","-.","-."]

    # fig, ax = subplots(3,1,figsize=[3.4, 5.3],sharex=true)
    fig, ax = subplots(4,1,figsize=[3.4, 7.3],sharex=true)

    ax[1].set_ylabel(L"cavity population $|\alpha|^2/N_\mathrm{at}$")
    ax[2].set_ylabel(L"order parameter $\vert\Phi\vert$")
    ax[3].set_ylabel(L"bunching parameter $\mathcal{B}$")
    ax[4].set_ylabel(L"final kinetic energy $E_\mathrm{kin}/\hbar\kappa$")
    ax[end].set_xlabel(L"pump strength $\sqrt{N}S/\kappa$")

    for (i,sim) in enumerate(sims)
        println("found $(size(sim)[end]) trajectories")
        println("-----------------------------------------")

        categories, par_list = split_sim_from_par(sim,true)

        S = Float64[]
        y1 = Float64[]
        y2 = Array{Float64,1}[]
        y3 = Float64[]
        y4 = Array{Float64,1}[]
        y5 = Float64[]
        y6 = Array{Float64,1}[]
        y7 = Float64[]
        y8 = Array{Float64,1}[]
        for x in categories
            push!(S,abs(x[1].p.S₁)*sqrt(x[1].p.N)/x[1].p.κ)

            m,s,q = expect(adaga,x)
            push!(y1,m[end])
            push!(y2,q[end])

            m,s,q = expect(absX,x)
            push!(y3,m[end])
            push!(y4,q[end])

            m,s,q = expect(Cos2,x)
            push!(y5,m[end])
            push!(y6,q[end])

            m,s,q = expect(Ekin,x)./x[1].p.κ
            push!(y7,sqrt.(m[end]))
            push!(y8,sqrt.(q[end]))
        end

        A = sortslices(hcat(S,y1,vcat(y2'...),y3,vcat(y4'...),y5,vcat(y6'...),y7,vcat(y8'...)),dims=1)
        S = A[:,1]
        y1 = A[:,2]
        y2 = A[:,3:4]
        y3 = A[:,5]
        y4 = A[:,6:7]
        y5 = A[:,8]
        y6 = A[:,9:10]
        y7 = A[:,11]
        y8 = A[:,12:13]

        # matplotlib[:rc]("axes", labelpad=1)

        # if par_list[1].temp == 0.0
        #     label = "\$temp=0\$"
        # elseif par_list[1].temp < par_list[1].κ
        #     label = "\$temp=\\kappa/"*string(trunc(Int,par_list[1].κ/par_list[1].temp))*"\$"
        # else
        #     label = "\$temp="*string(trunc(Int,par_list[1].temp/par_list[1].κ))*"\\kappa\$"
        # end

        
        # if par_list[1].Δₑ == 0.0
        #     label = "\$Δₑ=0\$"
        # elseif par_list[1].Δₑ < par_list[1].κ
        #     label = "\$Δₑ=\\kappa/"*string(trunc(Int,par_list[1].κ/par_list[1].Δₑ))*"\$"
        # else
        #     label = "\$Δₑ="*string(trunc(Int,par_list[1].Δₑ/par_list[1].κ))*"\\kappa\$"
        # end

                
        if par_list[1].temp == 0.0
            label = "\$δD=0\$"
        elseif par_list[1].Δₑ < par_list[1].κ
            label = "\$δD=\\kappa/"*string(trunc(Int,par_list[1].κ/par_list[1].temp))*", Δₑ=\\kappa/"*string(trunc(Int,par_list[1].κ/par_list[1].Δₑ))*", N_{at}="*string(trunc(Int,par_list[1].N))*"\$"
        else
            label = "\$δD="*string(trunc(Int,par_list[1].temp/par_list[1].κ))*"\\kappa\$"
        end

        ax[1].plot(S,y1.*1/sim[1].p.N,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        ax[1].fill_between(S,y2[:,1]*1/sim[1].p.N,y2[:,2]*1/sim[1].p.N,color=colorlist[i],alpha=0.2,linewidth=0.1)

        ax[2].plot(S,y3,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        ax[2].fill_between(S,y4[:,1],y4[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)

        ax[3].plot(S,y5,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        ax[3].fill_between(S,y6[:,1],y6[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)

        ax[4].plot(S,y7,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        ax[4].fill_between(S,y8[:,1],y8[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)
    end

    ax[1].legend(handlelength=2.5,loc="upper left",bbox_to_anchor=(0.2, 1.08),framealpha=1)
    fig.tight_layout(h_pad=0.)

    for i in 1:4
        letter = Char(Int('a')+i-1)
        ax[i].text(0.05,0.87,"("*letter*")",transform=ax[i].transAxes)
    end
    return fig, ax
end

function plot_ordering_vs_S(filename::String,sims::Array{Sol,1}...)
    fig, ax = plot_ordering_vs_S(sims...)
    fig.savefig(filename,dpi=1200)
end


# CHANGE TO AVOID USING PLOTS; CHANGE TO USE PYPLOT INSTEAD

function plot_interp_threshold_adaga(sims::Array{Sol, 1}...)

    colorlist = ["C1","C2","C3","C4","C5","C6"]
    linelist = ["-","--",":","-.","-."]
    fig, ax1 = plt.subplots()
    left, bottom, width, height = [0.25, 0.6, 0.2, 0.2]
    ax2 = fig.add_axes([left, bottom, width, height])

    for (i,sim) in enumerate(sims)
        println("found $(size(sim)[end]) trajectories")
        println("-----------------------------------------")

        categories, par_list = split_sim_from_par(sim, true)
        
        S = Float64[]
        y1 = Float64[]
        y2 = Array{Float64,1}[]
        
        for x in categories
            push!(S,abs(x[1].p.S₁))#*sqrt(x[1].p.N)/x[1].p.κ) #extract coupling strengths from sorted_sims
        
            m, s, q = expect(adaga, x) #extract mean, standard dev and 90 quantile
            push!(y1, m[end])
            push!(y2, q[end])
        end 

        A = sortslices(hcat(S,y1,vcat(y2'...)),dims=1)
        S = A[:,1]
        y1 = A[:,2]
        y2 = A[:,3:4]

        if par_list[1].Δₑ != 0.0
            label = "\$\\Delta_e=\\kappa/"*string(trunc(Int,par_list[1].κ/par_list[1].Δₑ))*"\$"
        else
            label = "\$\\Delta_e=0\$"
        end

        #Need to make S into an equally spaced Tuple
        x = S[1]:1.0:S[end]
            
        itp_linear = linear_interpolation(x, y1)
        itp_cubic = cubic_spline_interpolation(x, y1)
        # Interpolation functions
        f_linear(x) = itp_linear(x)
        f_cubic(x) = itp_cubic(x)
        
        x_smooth = S[1]:0.1:S[end] # smoother interval, necessary for cubic spline



        ax1.plot(x, y1, ls=linelist[1],color=colorlist[i], label=latexstring(label), alpha=0.8) #data
        ax1.fill_between(x,y2[:,1],y2[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)
        #ax1.plot(x_smooth, f_linear.(x_smooth), ls=linelist[i+1],color=colorlist[i+1], label="linear", alpha=0.8, markersize=0.5) #linear
        ax1.plot(x_smooth, f_cubic.(x_smooth), ls=linelist[2],color=colorlist[i], label="cubic", alpha=0.5, markersize=0.5)

        # Customize plot
        ax1.set_xlabel(L"|S₁|=|S₂|")
        ax1.set_ylabel(L"\langle a^\dagger a\rangle")
        ax1.legend(loc="best")
        #Inset 
        # Calculate second derivative using Central Finite Difference
        delta_x = x_smooth[2] - x_smooth[1]
        f_cubic_second_derivative(x) = (f_cubic(x - delta_x) - 2 * f_cubic(x) + f_cubic(x + delta_x)) / delta_x^2

        inset_x_range = 34.5:0.01:41
        inset_y_values = f_cubic_second_derivative.(inset_x_range)
        ax2.plot(inset_x_range, inset_y_values, color="black", ls=linelist[i], alpha=0.5, label=latexstring(label))#, inset=(1, bbox(0.1, 0.5, 0.3, 0.3)), subplot=2, linestyle=:dash, linewidth=5,, label=nothing, grid=false, xtickfontsize=10, ytickfontsize=10, xlabelfontsize=15, ylabelfontsize=15)
        ax2.set_xlabel(L"S")
        ax2.set_ylabel(L"\partial_S^2⟨a^\dagger a\rangle")
        # Draw a horizontal line at y=0 in the inset plot
        ax2.axhline(y=0, color=:black, alpha=0.5, label=nothing)
        ax2.legend(loc="best")

        closest_to_zero_index = argmin(abs.(inset_y_values))
        closest_to_zero_x = inset_x_range[closest_to_zero_index]
        closest_to_zero_y = inset_y_values[closest_to_zero_index]
        ax2.scatter([closest_to_zero_x], [closest_to_zero_y], color=:black, alpha=0.8)
    end
end

function plot_interp_threshold_adaga(solorsim,filename::String)
    fig, ax = plot_interp_threshold_adaga(solorsim)
    fig[:savefig](filename)
end

function plot_interp_threshold_order_param_delta(sims::Array{Sol, 1}...)

    colorlist = ["C1","C2","C3","C4","C5","C6"]
    linelist = ["-","--",":","-.","-."]
    fig, ax1 = plt.subplots()
    left, bottom, width, height = [0.25, 0.6, 0.2, 0.2]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax1.set_xlabel(L"|S₁|=|S₂|")
    ax1.set_ylabel(L"order parameter $\vert\Phi\vert$")
    ax2.set_xlabel(L"|S₁|=|S₂|")
    ax2.set_ylabel(L"\partial_S^2\vert\Phi\vert")



    for (i,sim) in enumerate(sims)
        println("found $(size(sim)[end]) trajectories")
        println("-----------------------------------------")

        categories, par_list = split_sim_from_par(sim, true)
        
        S = Float64[]
        y1 = Float64[]
        y2 = Array{Float64,1}[]
        
        for x in categories
            push!(S,abs(x[1].p.S₁))#*sqrt(x[1].p.N)/x[1].p.κ) #extract coupling strengths from sorted_sims
        
            m, s, q = expect(absX, x) #extract mean, standard dev and 90 quantile
            push!(y1, m[end])
            push!(y2, q[end])
        end 

        A = sortslices(hcat(S,y1,vcat(y2'...)),dims=1)
        S = A[:,1]
        y1 = A[:,2]
        y2 = A[:,3:4]

        if par_list[1].Δₑ != 0.0
            label = "\$\\Delta_e=\\kappa/"*string(trunc(Int,par_list[1].κ/par_list[1].Δₑ))*"\$"
        else
            label = "\$\\Delta_e=0\$"
        end

        #Need to make S into an equally spaced Tuple
        x = S[1]:1.0:S[end]
            
        itp_linear = linear_interpolation(x, y1)
        itp_cubic = cubic_spline_interpolation(x, y1)
        # Interpolation functions
        f_linear(x) = itp_linear(x)
        f_cubic(x) = itp_cubic(x)
        
        x_smooth = S[1]:0.1:S[end] # smoother interval, necessary for cubic spline



        ax1.plot(x, y1, ls=linelist[1],color=colorlist[i], label=latexstring(label), alpha=0.8) #data
        ax1.fill_between(x,y2[:,1],y2[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)
        #ax1.plot(x_smooth, f_linear.(x_smooth), ls=linelist[i+1],color=colorlist[i+1], label="linear", alpha=0.8, markersize=0.5) #linear
        ax1.plot(x_smooth, f_cubic.(x_smooth), ls=linelist[2],color=colorlist[i], label="cubic", alpha=0.5, markersize=0.5)
        ax1.legend(loc="best")

        #Inset 
        # Calculate second derivative using Central Finite Difference
        delta_x = x_smooth[2] - x_smooth[1]
        f_cubic_second_derivative(x) = (f_cubic(x - delta_x) - 2 * f_cubic(x) + f_cubic(x + delta_x)) / delta_x^2

        inset_x_range = 33:0.01:37
        inset_y_values = f_cubic_second_derivative.(inset_x_range)
        ax2.plot(inset_x_range, inset_y_values, color=colorlist[i], ls=linelist[i], alpha=0.5, label=latexstring(label))#, inset=(1, bbox(0.1, 0.5, 0.3, 0.3)), subplot=2, linestyle=:dash, linewidth=5,, label=nothing, grid=false, xtickfontsize=10, ytickfontsize=10, xlabelfontsize=15, ylabelfontsize=15)
        # Draw a horizontal line at y=0 in the inset plot
        ax2.axhline(y=0, color=:black, alpha=0.5, label=nothing)
        ax2.legend(loc="best")

        closest_to_zero_index = argmin(abs.(inset_y_values))
        closest_to_zero_x = inset_x_range[closest_to_zero_index]
        closest_to_zero_y = inset_y_values[closest_to_zero_index]
        ax2.scatter([closest_to_zero_x], [closest_to_zero_y], color=:black, alpha=0.8)
    end
end

function plot_interp_threshold_order_param_temp(sims::Array{Sol, 1}...)

    colorlist = ["C1","C2","C3","C4","C5","C6"]
    linelist = ["-","--",":","-.","-."]
    fig, ax1 = plt.subplots()
    left, bottom, width, height = [0.25, 0.6, 0.2, 0.2]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax1.set_xlabel(L"|S₁|=|S₂|")
    ax1.set_ylabel(L"order parameter $\vert\Phi\vert$")
    ax2.set_xlabel(L"|S₁|=|S₂|")
    ax2.set_ylabel(L"\partial_S^2\vert\Phi\vert")



    for (i,sim) in enumerate(sims)
        println("found $(size(sim)[end]) trajectories")
        println("-----------------------------------------")

        categories, par_list = split_sim_from_par(sim, true)
        
        S = Float64[]
        y1 = Float64[]
        y2 = Array{Float64,1}[]
        
        for x in categories
            push!(S,abs(x[1].p.S₁))#*sqrt(x[1].p.N)/x[1].p.κ) #extract coupling strengths from sorted_sims
        
            m, s, q = expect(absX, x) #extract mean, standard dev and 90 quantile
            push!(y1, m[end])
            push!(y2, q[end])
        end 

        A = sortslices(hcat(S,y1,vcat(y2'...)),dims=1)
        S = A[:,1]
        y1 = A[:,2]
        y2 = A[:,3:4]

        if par_list[1].temp != 0.0
            label = "\$temp=\\kappa/"*string(trunc(Int,par_list[1].κ/par_list[1].temp))*"\$"
        else
            label = "\$temp=0\$"
        end

        #Need to make S into an equally spaced Tuple
        x = S[1]:1.0:S[end]
            
        itp_linear = linear_interpolation(x, y1)
        itp_cubic = cubic_spline_interpolation(x, y1)
        # Interpolation functions
        f_linear(x) = itp_linear(x)
        f_cubic(x) = itp_cubic(x)
        
        x_smooth = S[1]:0.1:S[end] # smoother interval, necessary for cubic spline



        ax1.plot(x, y1, ls=linelist[1],color=colorlist[i], label=latexstring(label), alpha=0.8) #data
        ax1.fill_between(x,y2[:,1],y2[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)
        #ax1.plot(x_smooth, f_linear.(x_smooth), ls=linelist[i+1],color=colorlist[i+1], label="linear", alpha=0.8, markersize=0.5) #linear
        ax1.plot(x_smooth, f_cubic.(x_smooth), ls=linelist[2],color=colorlist[i], label="cubic", alpha=0.5, markersize=0.5)
        ax1.legend(loc="best")

        #Inset 
        # Calculate second derivative using Central Finite Difference
        delta_x = x_smooth[2] - x_smooth[1]
        f_cubic_second_derivative(x) = (f_cubic(x - delta_x) - 2 * f_cubic(x) + f_cubic(x + delta_x)) / delta_x^2

        inset_x_range = 31:0.01:36
        inset_y_values = f_cubic_second_derivative.(inset_x_range)
        ax2.plot(inset_x_range, inset_y_values, color="black", ls=linelist[i], alpha=0.5, label=latexstring(label))#, inset=(1, bbox(0.1, 0.5, 0.3, 0.3)), subplot=2, linestyle=:dash, linewidth=5,, label=nothing, grid=false, xtickfontsize=10, ytickfontsize=10, xlabelfontsize=15, ylabelfontsize=15)
        # Draw a horizontal line at y=0 in the inset plot
        ax2.axhline(y=0, color=:black, alpha=0.5, label=nothing)
        ax2.legend(loc="best")

        closest_to_zero_index = argmin(abs.(inset_y_values))
        closest_to_zero_x = inset_x_range[closest_to_zero_index]
        closest_to_zero_y = inset_y_values[closest_to_zero_index]
        ax2.scatter([closest_to_zero_x], [closest_to_zero_y], color=:black, alpha=0.8)
    end
end



####################################################################################





####################################################################################

function plot_initial_conditions(sim::Array{Sol,1})

    u0 = join_trajectories(sim,1)

    N::Int = sum(try sol.p.N catch; sol.p[10] end for sol in sim) # for backwards compatibility
    nbins = trunc(Int, sqrt(N))

    rangex = range(0,length=trunc(Int,sqrt(N)),stop=2pi)
    pdfx0 = fit(Histogram,mod2pi.(u0[1:N]),rangex)
    pdfx0 = normalize(pdfx0)

    minp = minimum(u0[N+1:2N])
    maxp = maximum(u0[N+1:2N])
    rangep = range(minp,length=trunc(Int,sqrt(N)),stop=maxp)
    pdfp0 = fit(Histogram,u0[N+1:2N],rangep)
    pdfp0 = normalize(pdfp0)

    ranges = range(-1.1,length=trunc(Int,sqrt(N)),stop=1.1)
    pdfsx = fit(Histogram,u0[2N+1:3N],ranges)
    pdfsy = fit(Histogram,u0[3N+1:4N],ranges)
    pdfsz = fit(Histogram,u0[4N+1:5N],ranges)
    pdfsx = normalize(pdfsx)
    pdfsy = normalize(pdfsy)
    pdfsz = normalize(pdfsz)


    matplotlib[:rc]("axes", labelpad=2.)

    fig, ax = subplots(2,2,figsize=[6.2, 4.6])

    ax[1, 1][:set_ylabel]("distribution")
    ax[1, 1][:set_xlabel](L"atom position $x$")
    ax[1, 1][:set_ylim]([0.,maximum(pdfx0.weights)*1.1])
    ax[1, 1][:set_xticks](pi/2*collect(0:4))
    ax[1, 1][:set_xticklabels]([L"0",L"\pi/2",L"\pi",L"3\pi/2",L"2\pi"])
    ax[1, 1][:step](pdfx0.edges[1][1:end-1],pdfx0.weights,label="initial",where="post")

    ax[1, 2][:set_ylabel]("distribution")
    ax[1, 2][:set_xlabel](L"atom momentum $p$")
    ax[1, 2][:step](pdfp0.edges[1][1:end-1],pdfp0.weights,label="initial",where="post")

    ax[2, 1][:set_ylabel]("distribution")
    ax[2, 1][:set_xlabel](L"spins $\sigma^x$ and $\sigma^y$")
    ax[2, 1][:set_xlabel](L"spins $\sigma^x$ and $\sigma^y$")
    ax[2, 1][:step](pdfsx.edges[1][1:end-1],pdfsx.weights,label=L"\sigma^x",where="post")
    ax[2, 1][:step](pdfsy.edges[1][1:end-1],pdfsy.weights,label=L"\sigma^y",where="post")
    ax[2, 1][:legend]()

    ax[2, 2][:set_ylabel]("distribution")
    ax[2, 2][:set_xlabel](L"spins $\sigma^z$")
    ax[2, 2][:step](pdfsz.edges[1][1:end-1],pdfsz.weights,label=L"\sigma^z",where="post")


    fig[:tight_layout](h_pad=0., w_pad=-0.)

    return fig, ax
end

function plot_initial_conditions(sim::Array{Sol,1},filename::String)
    fig, ax = plot_initial_conditions(sim)
    fig[:savefig](filename)
end


function plot_initial_spin(sim::Array{Sol,1})
    u0 = join_trajectories(sim,1)

    N::Int = sum(try sol.p.N catch; sol.p[10] end for sol in sim) # for backwards compatibility
    nbins = trunc(Int, sqrt(N))

    ranges = range(-1.1,length=trunc(Int,sqrt(N)),stop=1.1)
    pdfsx = fit(Histogram,u0[2N+1:3N],ranges)
    pdfsy = fit(Histogram,u0[3N+1:4N],ranges)
    pdfsz = fit(Histogram,u0[4N+1:5N],ranges)
    pdfsx = normalize(pdfsx)
    pdfsy = normalize(pdfsy)
    pdfsz = normalize(pdfsz)

    matplotlib[:rc]("axes", labelpad=2.)

    fig, ax = subplots()

    ax[:set_ylabel]("distribution")
    ax[:set_xlabel](L"spins $\sigma^x$ and $\sigma^y$")
    ax[:set_xlabel](L"spins $\sigma^x$ and $\sigma^y$")
    ax[:step](pdfsx.edges[1][1:end-1],pdfsx.weights,label=L"\sigma^x",where="post")
    ax[:step](pdfsy.edges[1][1:end-1],pdfsy.weights,label=L"\sigma^y",where="post")
    ax[:legend]()

    # ax[1, 2][:set_ylabel]("distribution")
    # ax[1, 2][:set_xlabel](L"spins $\sigma^z$")
    # ax[1, 2][:step](pdfsz.edges[1][1:end-1],pdfsz.weights,label=L"\sigma^z",where="post")

    fig[:tight_layout](h_pad=0., w_pad=-0.)

    return fig, ax
end

function plot_initial_spin(sim::Array{Sol,1},filename::String)
    fig, ax = plot_initial_spin(sim)
    fig[:savefig](filename)
end


function plot_spinspositionhisto(sim::Array{Sol,1})
    # matplotlib[:rc]("axes", labelpad=1)
    matplotlib.rc("image", cmap="inferno")

    sim_list = Array{Sol,1}[]

    for i in split_sim(sim)
        if length(i)>0
            push!(sim_list,i)
        end
    end

    nrows = length(sim_list)

    # fig, ax = subplots(nrows,2,figsize=[4.65, 2.3*nrows],sharex=true,sharey=false)
    fig, ax = subplots(nrows,2,figsize=[3.7, 1.5*nrows],sharex=true,sharey=false)

    for (i,a) in enumerate(sim_list)
        u1 = join_trajectories(a)
        N::Int = sum(try sol.p.N catch; sol.p[10] end for sol in a) # for backwards compatibility
        x = mod2pi.(u1[1:N])/(2pi)

        nbins = trunc(Int, ^(N,1//3))
        println("number of bins: "*string(nbins))

        ax1 = ax[i,1]
        ax1.set_ylabel(L"spin $\langle\sigma_x\rangle$")
        ax1.set_yticks(collect(-1:0.5:1))
        ax1.set_yticklabels([L"-1",L"-\frac{1}{2}",L"0",L"\frac{1}{2}",L"1"])
        ax1.hist2d(x,u1[2N+1:3N],bins=nbins,density=true)# norm=matplotlib.colors.LogNorm())

        ax2 = ax[i,2]
        ax2.set_ylabel(L"spin $\langle\sigma_y\rangle$")
        ax2.set_yticks(collect(-1:0.5:1))
        ax2.set_yticklabels([L"-1",L"-\frac{1}{2}",L"0",L"\frac{1}{2}",L"1"])
        cs = ax2.hist2d(x,u1[3N+1:4N], bins=nbins,density=true)

        ax2.set_xticks([cs[2][1],0.5,cs[2][end]])
        ax2.set_xticklabels(["0","0.5","1"])

    end

    fig.tight_layout()

    fig.subplots_adjust(top=0.85,bottom=0.13)
    cbar_ax = fig.add_axes([0.2, 0.92, 0.7, 0.02])
    cbar = fig.colorbar(matplotlib.cm.ScalarMappable(),cax=cbar_ax,orientation="horizontal")
    # cbar_ax.xaxis.set_ticks_position("top")
    cbar_ax.xaxis.set_label_position("top")
    cbar_ax.xaxis.set_tick_params(pad=1)
    cbar.set_label("density (a.u.)",labelpad=3.)

    fig.text(0.5, 0.03, L"atom position mod $\lambda_\mathrm{c}$ (units of $\lambda_\mathrm{c}$)", ha="center")

    return fig, ax

end

function plot_spinspositionhisto(filename::String,sim::Array{Sol,1})
    fig, ax = plot_spinspositionhisto(sim)
    fig.savefig(filename)
end

function plot_XY(sim::Array{Sol,1})

    a,b,c,d = split_sim(sim)

    tlist = sim[1].t

    matplotlib[:rc]("axes", labelpad=1)

    sf = 0.83
    fig, ax = subplots(1,2,figsize=[4.25*sf, 2*sf],sharey="row")

    ax[1][:set_ylabel]("X",labelpad=-2)
    ax[1][:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
    ax[2][:set_ylabel]("Y",labelpad=7)
    ax[2][:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
    # ax[1][:set_xscale]("log")
    # ax[2][:set_xscale]("log")

    for (i,part) in enumerate([a,d])
        if isempty(part)==true
            println("Array $i is empty")
        else
            Xs,X_std,X_q90 = expect(X,part)
            X_q90 = hcat(X_q90...)
            Ys,Y_std,Y_q90 = expect(Y,part)
            Y_q90 = hcat(Y_q90...)

            ax[1][:plot](tlist.+1,Xs,color="C0")
            # ax[1][:fill_between](tlist.+1,X_q90[1,:],X_q90[2,:],color="C0",alpha=0.2)
            ax[1][:fill_between](tlist.+1,Xs.+X_std,Xs.-X_std,color="C0",alpha=0.5)

            ax[2][:plot](tlist.+1,Ys,color="C1")
            # ax[2][:fill_between](tlist.+1,Y_q90[1,:],Y_q90[2,:],color="C1",alpha=0.2)
            ax[2][:fill_between](tlist.+1,Ys.+Y_std,Ys.-Y_std,color="C1",alpha=0.5)
        end
    end

    fig[:tight_layout](h_pad=0.4, w_pad=0.)

    return fig, ax
end

function plot_XY(sim::Array{Sol,1},filename::String)
    fig, ax = plot_XY(sim)
    fig[:savefig](filename)
end

function plot_XY_derivative(sim::Array{Sol, 1})
    a, b, c, d = split_sim(sim)
    tlist = sim[1].t

    matplotlib[:rc]("axes", labelpad=1)

    sf = 0.83
    fig, ax = subplots(1, 2, figsize=[4.25 * sf, 2 * sf], sharey="row")

    ax[1][:set_ylabel]("X", labelpad=-2)
    ax[1][:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
    ax[2][:set_ylabel]("Y", labelpad=7)
    ax[2][:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")

    for (i, part) in enumerate([a, d])
        Xs, X_std, X_q90 = expect(X, part)
        X_q90 = hcat(X_q90...)
        Ys, Y_std, Y_q90 = expect(Y, part)
        Y_q90 = hcat(Y_q90...)

        t_range = tlist[1]:10.0:tlist[end]
        # Cubic Spline Interpolation
        itp_cubic_X = cubic_spline_interpolation(t_range, Xs)
        itp_cubic_Y = cubic_spline_interpolation(t_range, Ys)

        # Interpolation functions
        f_cubic_X(x) = itp_cubic_X(x)
        f_cubic_Y(x) = itp_cubic_Y(x)

        # Create Smooth Time Vector
        delta_t = 0.01
        t_smooth = tlist[2]:delta_t:tlist[end-1]

        # Evaluate First Derivative with Interpolated Data (Central Finite Difference)

        f_cubic_first_derivative_X(x) = (f_cubic_X(x + delta_t) - f_cubic_X(x - delta_t)) / (2 * delta_t)
        f_cubic_first_derivative_Y(x) = (f_cubic_Y(x + delta_t) - f_cubic_Y(x - delta_t)) / (2 * delta_t)

        # Plotting First Derivatives
        ax[1][:plot](t_smooth, f_cubic_first_derivative_X.(t_smooth), linestyle="--", color="C2", label="dX/dt (spline)")
        ax[2][:plot](t_smooth, f_cubic_first_derivative_Y.(t_smooth), linestyle="--", color="C3", label="dY/dt (spline)")
    end

    fig[:tight_layout](h_pad=0.4, w_pad=0.)
    ax[1][:legend]()
    ax[2][:legend]()
    # ax[1][:set_xscale]("log")
    # ax[2][:set_xscale]("log")
    # ax[1][:set_yscale]("log")
    # ax[2][:set_yscale]("log")
    return fig, ax
end


function plot_Ekin(phase::Int, sims::Array{Sol,1}...)
    #sorted_sims = split_sim_from_par(sims) 
    fig, ax = subplots()
    ax[:set_xscale]("log")

    for u1 in sims
        sim = split_sim_from_par(u1)[phase] #1 to pick normal phase, end to pick SR phase
        temp = sim[1].p.temp
        Delta_e = sim[1].p.Δₑ
        S = real(sim[1].p.S₁)
        y,y_std,y_q90 = expect(Ekin,sim)
        y_q90 = hcat(y_q90...)
        tlist = sim[1].t
        
        matplotlib[:rc]("axes", labelpad=1.5)

        ax[:set_ylabel](L"$K_B T$ (units of ${\hbar\omega_\mathrm{r}}$ )")
        ax[:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
        #ax[:set_xscale]("log")
        ax[:plot](tlist,y, label = "\$N_{at}= "*string(trunc(Int, sim[1].p.N))*"\$")# temp=$temp, |S|=$S, Delta_e=$Delta_e")
        ax[:fill_between](tlist,y_q90[1,:],y_q90[2,:],alpha=0.2)

        ax[:legend](handlelength=2.5,loc="upper left",bbox_to_anchor=(0.2, 1.08),framealpha=1)  
        fig[:tight_layout]()
    end
    return fig, ax
end


function plot_Ekin(sim::Array{Sol,1},filename::String)
    fig, ax = plot_Ekin(sim)
    fig[:savefig](filename)
end

# function plot_Ekin(sim::Array{Sol,1})
#     #sorted_sims = split_sim_from_par(sims) 
#     fig, ax = subplots(1,1,figsize=[3.4, 2.3])

#         temp = sim[1].p.temp
#         S = real(sim[1].p.S₁)
#         D_e = sim[1].p.Δₑ
#         y,y_std,y_q90 = expect(Ekin,sim)
#         y_q90 = hcat(y_q90...)
#         tlist = sim[1].t
        
#         expression(t) = sim[1].p.Δc^2+sim[1].p.κ^2*1/(-4*sim[1].p.Δc) #Simons approx

#         analytic = [expression(t) for t in tlist]

#         efit =  fitexp(tlist,y, n=1) #fit with order n exponential

#         a = round(efit.a, digits=2)
#         b = round(efit.b, digits=2)
#         c = round(efit.c, digits=2)
#         matplotlib[:rc]("axes", labelpad=1.5)

#         ax[:set_ylabel](L"temperature")
#         ax[:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
#         #ax[:set_xscale]("log")
#         ax[:plot](tlist,y, label = "temp=$temp, |S|=$S, D_e = $D_e")
#         #ax[:plot](efit.x,efit.y, label = "y=$a*exp(-$Γ t)+$c)")
#         ax[:plot](efit.x,efit.y, label =  "fit, Γ=$b, δ^f=$c")
#         ax[:plot](tlist, analytic, label= "analytic approx")
#         ax[:fill_between](tlist,y_q90[1,:],y_q90[2,:],alpha=0.2)
#         #ax[:fill_between](tlist.+1,y.+y_std,y.-y_std,color=red,alpha=0.5)
#         println("A = $(efit.a) B = $(efit.b)) C = $(efit.c))")
#         ax[:legend](handlelength=2.5,loc="upper left",bbox_to_anchor=(0.2, 1.08),framealpha=1)  
#         fig[:tight_layout]()
#     return fig, ax
# end

function plot_Ekin_vs_S(sims::Array{Sol,1}...)

    colorlist = ["C1","C2","C3","C4","C5","C6"]
    linelist = ["-","--",":","-.","-."]

    # fig, ax = subplots(3,1,figsize=[3.4, 5.3],sharex=true)
    fig, ax = subplots()

    ax.set_ylabel(L"final ensemble temperature $K_B T$(units of $\hbar ω_R)$")
    ax.set_xlabel(L"pump strength $\sqrt{N}S/\kappa$")

    for (i,sim) in enumerate(sims)
        println("found $(size(sim)[end]) trajectories")
        println("-----------------------------------------")

        categories, par_list = split_sim_from_par(sim,true)

        S = Float64[]
        y1 = Float64[]
        y2 = Array{Float64,1}[]
        for x in categories
            push!(S,abs(x[1].p.S₁)*sqrt(x[1].p.N)/x[1].p.κ)
            m,s,q = expect(Ekin,x)
            push!(y1,m[end]./x[1].p.κ)
            push!(y2,q[end]./x[1].p.κ)
        end

        A = sortslices(hcat(S,y1,vcat(y2'...)),dims=1)
        S = A[:,1]
        y1 = A[:,2]
        y2 = A[:,3:4]
                
        if par_list[1].temp == 0.0
            label = "Initial temp\$=0\$"
        elseif par_list[1].temp < par_list[1].κ
            label = "Initial temp\$=\\kappa/"*string(trunc(Int,par_list[1].κ/par_list[1].temp))*", Δₑ=\\kappa/"*string(trunc(Int,par_list[1].κ/par_list[1].Δₑ))*"\$"
        else
            label = "Initial temp\$="*string(trunc(Int,par_list[1].temp/par_list[1].κ))*"\\kappa\$"
        end

        ax.plot(S,y1,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        ax.fill_between(S,y2[:,1],y2[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)
    end

    ax.legend(handlelength=2.5,loc="upper left",bbox_to_anchor=(0.2, 1.08),framealpha=1)
    fig.tight_layout(h_pad=0.)

    return fig, ax
end

function plot_init_final_momentum(sim::Array{Sol,1})
    colorlist = ["C1","C2","C3","C4","C5","C6"]
    linelist = ["-","--",":","-.","-."]
    linewidth=2

    final = Int(sim[end].t[end]./10)
    u0 = join_trajectories(sim,1)
    uf = join_trajectories(sim,final)
    N::Int = sum(try sol.p.N catch; sol.p[10] end for sol in sim) # for backwards compatibility
    nbins = trunc(Int, sqrt(N))

    fig, ax = subplots()

    minp = minimum(u0[N+1:2N])
    maxp = maximum(u0[N+1:2N])
    rangep = range(minp,length=trunc(Int,sqrt(N)),stop=maxp)
    pdfp0 = fit(Histogram,u0[N+1:2N],rangep)
    pdfpf = fit(Histogram,uf[N+1:2N],rangep)

    # #Gaussian fitting

    dist_init = fit(Normal, u0[N+1:2N])
    dist_final = fit(Normal, uf[N+1:2N])
    
    fit_init = pdf.(dist_init, rangep)
    fit_final = pdf.(dist_final, rangep)

    # # q-Gaussian fitting

    # @. model(x, p) = p[1]*(1-(1-p[2])*p[3]*x^2)^(1/(1-p[2]))#cos(x*p[2])+p[3]
    # xdata = collect(pdfpf.edges[1][1:end-1])
    # ydata = pdfpf.weights
    # p0 = [10, 1, 0.]
    # fitted = curve_fit(model, xdata, ydata, p0)
    # println("Fitted parameters", fitted.param)

    matplotlib[:rc]("axes", labelpad=2.)

    ax[:set_ylabel]("distribution")
    ax[:set_yscale]("log")
    ax[:set_ylim](1e-5, 0.1)
    ax[:set_xlabel](L"atom momentum $p$")
    ax[:step](normalize(pdfp0).edges[1][1:end-1],normalize(pdfp0).weights,label=L"t_0",where="post", color=colorlist[3], alpha=0.7, linestyle=linelist[1], linewidth=linewidth)
    ax[:step](normalize(pdfpf).edges[1][1:end-1],normalize(pdfpf).weights,label=L"t_{ss}",where="post", color=colorlist[2], alpha=0.7, linestyle=linelist[1], linewidth=linewidth)
    #ax[:plot](xdata,normalize(model(xdata,fitted.param)), label="q-Gaussian fit")
    ax[:plot](rangep, fit_init, label=L"Gaussian $t_0$ fit",color=colorlist[3] , alpha=0.7, linestyle=linelist[2], linewidth=linewidth)
    ax[:plot](rangep, fit_final, label=L"Gaussian $t_{ss}$ fit", color=colorlist[2], alpha=0.7, linestyle=linelist[2], linewidth=linewidth)
    fig[:tight_layout](h_pad=0., w_pad=-0.)
    ax[:legend](loc="best")
    return fig, ax
end

function plot_Ekin_vs_Deltae(sims::Array{Sol,1}...)
    color = ["C2","C3"]#,"C2","C3","C4","C5","C6"]
    line = "-"#,"--",":","-.","-."]

    fig, ax = subplots()

    ax.set_ylabel(L"final ensemble temperature $K_B T$(units of $\hbar ω_R)$")
    ax.set_xlabel(L"atomic detuning $Δₑ/\kappa$")
    

    for (i,sim) in enumerate(sims)
        println("found $(size(sim)[end]) trajectories")
        println("-----------------------------------------")
        Deltae = Float64[]
        y1 = Float64[]
        y2 = Float64[]#Array{Float64,1}[]
        y1_end = Float64[]
        y2_end = Float64[]#Array{Float64,1}[]

        categories, par_list = split_sim_from_par(sim,true)
        x = categories[1]
        x_end = categories[end]
            push!(Deltae,abs(x[1].p.Δₑ)/x[1].p.κ)
            m,s,q = expect(Ekin,x)
            push!(y1,m[end])
            push!(y2,s[end])
            
            m_end,s_end,q_end = expect(Ekin,x_end)
            push!(y1_end,m_end[end])
            push!(y2_end,s_end[end])
                
        # if par_list[1].temp == 0.0
        #     label = "Initial temp\$=0\$"
        # elseif par_list[1].temp < par_list[1].κ
        #     label = "Initial temp\$=κ/" * string(trunc(Int, par_list[1].κ / par_list[1].temp)) * "\\, S=" * string(trunc(Int, real(par_list[1].S₁))) * "\$"
        # else
        #     label = "Initial temp\$=" * string(trunc(Int, par_list[1].temp / par_list[1].κ)) * "κ\$"
        # end      

        ax[:errorbar](Deltae,y1,yerr=y2 ,ls=line,color=color[1], label="Normal", fmt="o", alpha=0.7)
        ax[:errorbar](Deltae,y1_end, yerr=y2_end, ls=line, color=color[2], label = "SR", fmt="o", alpha=0.7)
        #ax.fill_between(S,y2[:,1],y2[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)
    end

    #ax.legend(handlelength=2.5,loc="upper left",bbox_to_anchor=(0.2, 1.08),framealpha=1)
    fig.tight_layout(h_pad=0.)

    return fig, ax
end

function plot_total_spin(sim::Array{Sol, 1})
    sf = 0.83
    fig, ax = subplots(1,4,figsize=[4.45*sf, 2*sf],sharex=true, sharey="row")

    ax[1][:set_ylabel](L"$⟨S_x^2⟩$")
    #ax[1][:set_xscale]("log")
    ax[2][:set_ylabel](L"$⟨S_y^2⟩$")
    #ax[2][:set_xscale]("log")
    ax[3][:set_ylabel](L"$⟨S_z^2⟩$")
    #ax[3][:set_xscale]("log")
    ax[4][:set_ylabel](L"$\sum_i⟨S_i^2⟩$")
    #ax[4][:set_xscale]("log")

    # ax[2,1][:set_ylabel](L"$⟨S_x⟩$")
    # ax[2,1][:set_xscale]("log")
    # ax[2,2][:set_ylabel](L"$⟨S_y⟩$")
    # ax[2,2][:set_xscale]("log")
    # ax[2,3][:set_ylabel](L"$⟨S_z⟩$")
    # ax[2,3][:set_xscale]("log")

    tlist = sim[1].t

    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")

    S2x,Sx_err,q = expect(Sx2,sim) 
    S2y,Sy_err,q = expect(Sy2,sim)
    S2z,Sz_err,q = expect(Sz2,sim)

    S2tot = S2x .+ S2y.+S2z
    S2_err = Sx_err.+Sz_err.+Sz_err
    
    Ssx,Ssx_err,q = expect(Sx,sim) 
    Ssy,Ssy_err,q = expect(Sy,sim)
    Ssz,Ssz_err,q = expect(Sz,sim)

    ax[1][:plot](tlist.+1, S2x)
    ax[1][:fill_between](tlist.+1, S2x.+Sx_err, S2x.-Sx_err, color="C1", alpha=0.5)
    ax[1][:hlines](sqrt(3), minimum(tlist .+ 1), maximum(tlist .+ 1), linestyle="--", color="C3", alpha=0.8,  label=L"$\sqrt{3}$")#, fontsize=14)
    # ax[1][:text](maximum(tlist .- 25), sqrt(3) .- 0.1, L"$\sqrt{3}$", va="center", ha="right", color="C3",)
    ax[1][:legend]()
    ax[2][:plot](tlist.+1, S2y)
    ax[2][:fill_between](tlist.+1, S2y.+Sy_err, S2y.-Sy_err, color="C2", alpha=0.5)

    ax[3][:plot](tlist.+1, S2z)
    ax[3][:fill_between](tlist.+1, S2z.+Sz_err, S2z.-Sz_err, color="C3", alpha=0.5)

    ax[4][:plot](tlist.+1, S2tot)
    ax[4][:fill_between](tlist.+1, S2tot.+S2_err, S2tot.-S2_err, color="C4", alpha=0.5)
    ax[4][:hlines](3, minimum(tlist .+ 1), maximum(tlist .+ 1), linestyle="--", color="C3", alpha=0.8, label="3")
    ax[4][:legend]()
    # ax[2,1][:plot](tlist.+1, Ssx)
    # ax[2,1][:fill_between](tlist.+1, Ssx.+Ssx_err, Ssx.-Ssx_err, color="C1", alpha=0.5)

    # ax[2,2][:plot](tlist.+1, S2y)
    # ax[2,2][:fill_between](tlist.+1, Ssy.+Ssy_err, Ssy.-Ssy_err, color="C2", alpha=0.5)

    # ax[2,3][:plot](tlist.+1, S2z)
    # ax[2,3][:fill_between](tlist.+1, Ssz.+Ssz_err, Ssz.-Ssz_err, color="C3", alpha=0.5)
    return fig, ax
end

function plot_spin_length(sim::Array{Sol, 1})
    #sf = 0.83
    fig, ax = subplots()

    ax[:set_ylabel](L"$⟨S_x^2⟩$")
    ax[1][:set_xscale]("log")


    tlist = sim[1].t

    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")
    
    Sx,Sxe,q = expect(Sx,sim) 
    Sy,Sye,q = expect(Sy,sim)
    Sz,Sze,q = expect(Sz,sim)

    S2tot = S2x.^2 .+ S2y.^2 .+ S2z.^2
    S2_err = Sx_err.+Sz_err.+Sz_err

    ax[:plot](tlist.+1, S2tot)
    ax[:fill_between](tlist.+1, S2tot.+S2_err, S2tot.-S2_err, color="C1", alpha=0.5)
    ax[:hlines](3, minimum(tlist .+ 1), maximum(tlist .+ 1), linestyle="--", color="C3", alpha=0.8,  label=L"$3$")#, fontsize=14)
    # ax[1][:text](maximum(tlist .- 25), sqrt(3) .- 0.1, L"$\sqrt{3}$", va="center", ha="right", color="C3",)
    ax[:legend]()

    return fig, ax
end

function plot_spin_squeezing(sims::Array{Sol, 1}...)
    fig, ax = subplots()

    ax[:set_ylabel](L"$\xi^2=N\frac{(\Delta J_x)^2}{|⟨J_z⟩|^2}$")
    ax[:set_xscale]("log")

    for (i, sim) in enumerate(sims)
        println("found $(size(sim)[end]) trajectories")
        println("-----------------------------------------")

        categories, par_list = split_sim_from_par(sim,false)

        S = Float64[]
        ξ2 = Float64[]
        ξerr = Array{Float64,1}[]
        for x in categories
            push!(S,abs(x[1].p.S₁)*sqrt(x[1].p.N)/x[1].p.κ)
            Sxt,Sxe,q = expect(Sx,x) 
            Syt,Sye,q = expect(Sy,x)
            Szt,Sze,q = expect(Sz,x)

            Sx2t,Sxe2,q2 = expect(Sx2,x) 
            Sxvar = sqrt.(Sx2t[end]-Sxt[end].^2)

            push!(ξ2, N .* Sxvar.^2)
            push!(ξerr, N .* Sxe[end].^2)
        end
            # ξ2 = N .* Sxvar.^2# ./ (abs.(Sxt.^2+Syt.^2+Szt.^2))
            # ξe = N .* Sxe.^2
            #ax[:set_ylim](0,N)
            ax[:plot](S, ξ2)
            ax[:fill_between](S, ξ2.+ξerr, ξ2.-ξerr, color="C1", alpha=0.5)
    end

    fig.tight_layout(h_pad=0.)
    return fig, ax
end

function plot_single_spin_stats(sim::Array{Sol, 1})
    #sf = 0.83
    fig, ax = subplots(3,1,sharex=true)

    ax[1][:set_ylabel](L"$⟨s_x^2⟩$")
    ax[2][:set_ylabel](L"$⟨s_y^2⟩$")
    ax[3][:set_ylabel](L"$⟨s_z^2⟩$")
    #ax[:set_xscale]("log")


    tlist = sim[1].t

    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")
    
    sx2,sx2e,qx = expect(single_sx2,sim)
    qx = hcat(qx...)
    sy2,sy2e,qy = expect(single_sy2,sim) 
    qy = hcat(qy...)
    sz2,sz2e,qz = expect(single_sz2,sim)  
    qz = hcat(qz...)

    ax[1][:plot](tlist.+1, sx2)
    ax[1][:fill_between](tlist.+1, sx2.+sx2e, sx2.-sx2e, color="C1", alpha=0.5)
    ax[2][:plot](tlist.+1, sy2)
    ax[2][:fill_between](tlist.+1, sy2.+sy2e, sy2.-sy2e, color="C1", alpha=0.5)
    ax[3][:plot](tlist.+1, sz2)
    ax[3][:fill_between](tlist.+1, sz2.+sz2e, sz2.-sz2e, color="C1", alpha=0.5)
    #ax[:hlines](3, minimum(tlist .+ 1), maximum(tlist .+ 1), linestyle="--", color="C3", alpha=0.8,  label=L"$3$")#, fontsize=14)
    # ax[1][:text](maximum(tlist .- 25), sqrt(3) .- 0.1, L"$\sqrt{3}$", va="center", ha="right", color="C3",)
    ax[1][:legend]()

    return fig, ax
end

function plot_bunching(sim::Array{Sol,1})
    y,y_std,y_q90 = expect(Cos2,sim)
    y_q90 = hcat(y_q90...)
    tlist = sim[1].t

    matplotlib[:rc]("axes", labelpad=0.5)

    fig, ax = subplots(1,1,figsize=[3.25, 2.])

    color="C2"

    ax[:set_ylabel](L"bunching parameter $\mathcal{B}$")
    ax[:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
    # ax[:set_xscale]("log")
    ax[:plot](tlist.+1,y,color=color)
    ax[:fill_between](tlist.+1,y_q90[1,:],y_q90[2,:],color=color,alpha=0.2)
    #ax[:fill_between](tlist.+1,y.+y_std,y.-y_std,color=color,alpha=0.5)

    fig[:tight_layout]()

    return fig, ax
end

function plot_kurtosis(sim::Array{Sol,1})
    y,y_std,y_q90 = expect(kurt,sim)
    y_q90 = hcat(y_q90...)
    tlist = sim[1].t

    matplotlib[:rc]("axes", labelpad=0.5)

    fig, ax = subplots(1,1,figsize=[3.25, 2.])

    color="C2"

    ax[:set_ylabel](L"kurtosis $\mathcal{K}$")
    ax[:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
    # ax[:set_xscale]("log")
    ax[:plot](tlist.+1,y,color=color)
    ax[:fill_between](tlist.+1,y_q90[1,:],y_q90[2,:],color=color,alpha=0.2)
    #ax[:fill_between](tlist.+1,y.+y_std,y.-y_std,color=color,alpha=0.5)

    fig[:tight_layout]()

    return fig, ax
end