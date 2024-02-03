#module LPlots

#REMOVED MODULE CALLS


using DiffEqBase
using LaTeXStrings
using LinearAlgebra
using LsqFit
using Measurements
using Random: MersenneTwister
using PyPlot
using StatsBase
using PyCall


include("../src/selforg_core.jl")
#include("custom_functions.jl")
#using Selforg
#using CustomFunctions: extract_solution, Sol, intervallize_array

export plot_initial_conditions, plot_random_spins
export plot_observable, plot_observables, plot_resume
export plot_status, video_status, plot_status2
export plot_vs_phi, plot_vs_S, plot_vs_SDe
export plot_position, plot_XY, plot_Ekin, plot_adaga, plot_spinspositionhisto, plot_threshold, plot_vs_phi_talk, plot_Ekin_vs_S # plots for talks


function plot_Bloch()
    qt = pyimport("qutip")
    b = qt.Bloch()


    b.show()
end


function plot_random_spins(sol::Sol, seed=0)
    # choose 2 random spins
    N::Int = sol.p.N
    rng = MersenneTwister(seed) # random number generator
    s1, s2 = rand(rng, 1:N, 2)

    fig, ax = subplots(1, 3, figsize=[11, 3])
    ax[1][:set_ylabel](L"$\sigma_x$")
    ax[1][:plot](sol.t, sol.u[2N+s1,:], label="spin $s1")
    ax[1][:plot](sol.t, sol.u[2N+s2,:], label="spin $s2")
    ax[2][:set_ylabel](L"$\sigma_y$")
    ax[2][:plot](sol.t, sol.u[3N+s1,:])
    ax[2][:plot](sol.t, sol.u[3N+s2,:])
    ax[3][:set_ylabel](L"$\sigma_x$")
    ax[3][:plot](sol.t, sol.u[4N+s1,:])
    ax[3][:plot](sol.t, sol.u[4N+s2,:])

    for a in ax
        a[:set_xlabel]("time")
        a[:plot]([sol.t[1], sol.t[end]], [1, 1], ls=":", color="k")
        a[:plot]([sol.t[1], sol.t[end]], [-1, -1], ls=":", color="k")
    end

    ax[1][:legend](loc=1)
    fig[:tight_layout]()
end

function plot_random_spins(sol::RODESolution, seed=0)
    plot_random_spins(extract_solution(sol)[1], seed)
end

function plot_status2(sim::Array{Sol,1})
    N::Int = sum(try sol.p.N catch; sol.p[10] end for sol in sim) # for backwards compatibility
    nbins = trunc(Int, ^(N,1//3))
    plot_status2(sim,nbins)
end

function plot_status2(sim::Array{Sol,1},nbins::Int,idx::Int=size(sim[1].t)[1])
    u1 = join_trajectories(sim,idx)
    u0 = join_trajectories(sim,1)

    N::Int = sum(try sol.p.N catch; sol.p[10] end for sol in sim) # for backwards compatibility

    println("number of bins: "*string(nbins))

    matplotlib[:rc]("axes", labelpad=1)
    matplotlib[:rc]("image", cmap="cividis")

    x = mod2pi.(u1[1:N])

    fig, ax = subplots(4,2,figsize=[6.2, 9.])

    ax[1, 1][:set_ylabel](Sx.short_name)
    ax[1, 1][:hist2d](x,u1[2N+1:3N], bins=nbins)    # norm=matplotlib[:colors][:LogNorm]()

    ax[1, 2][:set_ylabel](Sy.short_name)
    ax[1, 2][:hist2d](x,u1[3N+1:4N], bins=nbins)

    ax[2, 1][:set_ylabel](Sz.short_name)
    ax[2, 1][:hist2d](x,u1[4N+1:5N], bins=nbins)

    ax[2, 2][:set_xlabel](Sx.short_name)
    ax[2, 2][:set_ylabel](Sy.short_name)
    ax[2, 2][:hist2d](u1[2N+1:3N],u1[3N+1:4N], bins=nbins)

    ax[3, 1][:set_xlabel](Sx.short_name)
    ax[3, 1][:set_ylabel](Sz.short_name)
    ax[3, 1][:hist2d](u1[2N+1:3N],u1[4N+1:5N], bins=nbins)

    ax[3, 2][:set_xlabel](Sz.short_name)
    ax[3, 2][:set_ylabel](Sy.short_name)
    ax[3, 2][:hist2d](u1[4N+1:5N],u1[3N+1:4N], bins=nbins)

    ax[4, 1][:set_ylabel]("total spin")
    ax[4, 1][:set_xlabel]("atom")
    ax[4, 1][:plot](u1[2N+1:3N].^2 .+ u1[3N+1:4N].^2 .+ u1[4N+1:5N].^2)
    # ax[4, 1][:plot](u0[4N+1:5N])

    ax[4, 2][:set_ylabel]("total spin")
    ax[4, 2][:set_xlabel]("atom")
    ax[4, 2][:hist](u1[2N+1:3N].^2 .+ u1[3N+1:4N].^2 .+ u1[4N+1:5N].^2)

    for i in [ax[1,1],ax[1,2],ax[2,1]]
        i[:set_xlabel](L"atom position mod $2\pi$")
    end

    fig[:tight_layout](h_pad=0., w_pad=-0.)

    return fig, ax

end

function plot_status2(sim::Array{Sol,1},filename::String)
    fig, ax = plot_status2(sim)
    fig[:savefig](filename)
end

function plot_status2(sim::Array{Sol,1},nbins::Int,filename::String)
    fig, ax = plot_status2(sim,nbins)
    fig[:savefig](filename)
end

function plot_status(sim::Array{Sol,1},idx::Int=size(sim[1].t)[1])

    u0 = join_trajectories(sim,1)
    u1 = join_trajectories(sim,idx)

    N::Int = sum(try sol.p.N catch; sol.p[10] end for sol in sim) # for backwards compatibility
    nbins = trunc(Int, sqrt(N))

    rangex = range(0,length=trunc(Int,sqrt(N)),stop=2pi)
    pdfx0 = fit(Histogram,mod2pi.(u0[1:N]),rangex)
    pdfx1 = fit(Histogram,mod2pi.(u1[1:N]),rangex)

    minp = minimum(vcat(u0[N+1:2N],u1[N+1:2N]))
    maxp = maximum(vcat(u0[N+1:2N],u1[N+1:2N]))
    rangep = range(minp,length=trunc(Int,sqrt(N)),stop=maxp)
    pdfp0 = fit(Histogram,u0[N+1:2N],rangep)
    pdfp1 = fit(Histogram,u1[N+1:2N],rangep)

    matplotlib[:rc]("axes", labelpad=1)

    fig, ax = subplots(4,2,figsize=[6.2, 8.6])

    ax[1, 1][:set_ylabel]("histogram")
    ax[1, 1][:set_xlabel](L"atom position mod $2\pi$")
    ax[1, 1][:set_xticks](pi/2*collect(0:4))
    ax[1, 1][:set_xticklabels]([L"0",L"\pi/2",L"\pi",L"3\pi/2",L"2\pi"])
    ax[1, 1][:step](pdfx0.edges[1][1:end-1],pdfx0.weights,label="initial",where="post")
    ax[1, 1][:step](pdfx1.edges[1][1:end-1],pdfx1.weights,label="final",where="post")
    ax[1,1][:legend](loc=8,ncol=2)

    ax[1, 2][:set_xlabel]("time")
    ax[1, 2][:set_ylabel]("positions of 5 random atoms")
    sol = rand(sim)
    for i in rand(1:trunc(Int,try sol.p.N catch; sol.p[10] end),3) # for backwards compatibility
        ax[1, 2][:plot](sol.t,sol.u[i,:])
    end

    ax[2, 1][:set_ylabel]("histogram")
    ax[2, 1][:set_xlabel](L"atom momentum $p$")
    ax[2, 1][:step](pdfp0.edges[1][1:end-1],pdfp0.weights,label="initial",where="post")
    ax[2, 1][:step](pdfp1.edges[1][1:end-1],pdfp1.weights,label="final",where="post")
    ax[2,1][:legend]()

    ax[3, 1][:set_ylabel]("histogram")
    ax[3, 1][:set_xlabel](L"spins $\sigma_x$")
    ax[3, 1][:hist]((u0[2N+1:3N],u1[2N+1:3N]),histtype="step",
                    bins=nbins,label=["initial", "final"])
    ax[3, 1][:legend]()

    ax[3,2][:set_ylabel]("histogram")
    ax[3,2][:set_xlabel](L"spins $\sigma_y$")
    ax[3,2][:hist]((u0[3N+1:4N],u1[3N+1:4N]),histtype="step",
                   bins=nbins,label=["initial", "final"])
    ax[3,2][:legend]()

    ax[2,2][:set_ylabel]("histogram")
    ax[2,2][:set_xlabel](L"spins $\sigma_z$")
    ax[2,2][:hist]((u0[4N+1:5N],u1[4N+1:5N]),histtype="step",
                   bins=nbins,label=["initial", "final"])
    ax[2,2][:legend]()


    fig[:tight_layout](h_pad=0., w_pad=-0.)

    return fig, ax
end

function plot_status(sim::Array{Sol,1},filename::String)
    fig, ax = plot_status(sim,size(sim[1].t)[1])
    fig[:savefig](filename)
end

function plot_status(sim::Array{Sol,1},idx::Int,filename::String)
    fig, ax = plot_status(sim,idx)
    fig[:savefig](filename)
end

function video_status(sim::Array{Sol,1})
    for i in 1:size(sim[1].t)[1]
        filename = "./video/fig_"*lpad(repr(i),4,"0")
        fig, ax = plot_status(sim,i)
        fig[:savefig](filename)
    end
end

function plot_initial_conditions_old(N::Int, u0)
    fig, ax = subplots(3, 2, figsize=[10, 10])
    nbins = trunc(Int, sqrt(N))

    ax[1, 1][:set_ylabel]("histogram")
    ax[1, 1][:set_xlabel](L"atom position $x$")
    ax[1, 1][:hist](u0[1:N], bins=nbins)
    ax[1, 2][:set_xlabel](L"atom number $i$")
    ax[1, 2][:set_ylabel](L"atoms positions $x_i$")
    ax[1, 2][:plot](u0[1:N], ".")

    ax[2, 1][:set_ylabel]("histogram")
    ax[2, 1][:set_xlabel](L"atom momentum $p$")
    ax[2, 1][:hist](u0[N+1:2N], bins=nbins)
    ax[2, 2][:set_xlabel](L"atom number $i$")
    ax[2, 2][:set_ylabel](L"atom momentum $p_i$")
    ax[2, 2][:plot](u0[N+1:2N], ".")

    ax[3, 1][:set_ylabel]("histogram")
    ax[3, 1][:set_xlabel](L"spins $\sigma_x$, $\sigma_y$ and $\sigma_z$")
    ax[3, 1][:hist]((u0[2N+1:3N],u0[3N+1:4N],u0[4N+1:5N]),bins=3nbins)
    ax[3, 2][:set_xlabel](L"atom number $i$")
    ax[3, 2][:set_ylabel]("spins")
    ax[3, 2][:plot](u0[2N+1:3N], label=L"$\sigma_x$")
    ax[3, 2][:plot](u0[3N+1:4N], label=L"$\sigma_y$")
    ax[3, 2][:plot](u0[4N+1:5N], label=L"$\sigma_z$")
    ax[3, 2][:plot](@. u0[2N+1:3N]^2 + u0[3N+1:4N]^2 + u0[4N+1:5N]^2)

    ax[3, 2][:legend]()
    fig[:tight_layout]()

    return fig, ax
end

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



function plot_observable(o::Observable, sol::DiffEqBase.RODESolution)
    N::Int = sol.prob.p.N
    y = o.s_traj(sol)

    fig, ax = subplots(1,1)
    ax[:set_xlabel]("time")
    ax[:set_ylabel](o.name)
    ax[:set_title](o.name)
    ax[:plot](sol.t, y)
    fig[:tight_layout]()

    return fig, ax
end

function plot_observable(o::Observable, sim::EnsembleSolution)
    N::Int = sim[1].prob.p.N
    y, ystd = o.montec(sim)

    fig, ax = subplots(1,1)
    ax[:set_xlabel]("time")
    ax[:set_ylabel](o.name)
    ax[:set_title](o.name)
    ax[:plot](sim[1].t, y)
    ax[:fill_between](sim[1].t,y.+ystd,y.-ystd,alpha=0.5)
    fig[:tight_layout]()

    return fig, ax
end

function plot_observable(obs::Array{Observable}, solorsim)
    for o in obs
        plot_observable(o,solorsim)
    end
end

function plot_observable(o,solorsim,filename::String)
    fig, ax = plot_observable(o,solorsim)
    fig[:savefig](filename)
end

function plot_observables(obs::Array{Observable}, sol::DiffEqBase.RODESolution)
    fig, ax = subplots(1,1)
    ax[:set_xlabel]("time")
    ax[:set_ylabel]("observables")
    title = ""
    for o in obs
        ax[:plot](sol.t, o.s_traj(sol), label=o.name)
        title *= o.name * ", "
    end
    ax[:set_title](title[1:end-2])
    ax[:legend]()
    fig[:tight_layout]()

    return fig, ax
end

function plot_observables(obs::Array{Observable}, sim::EnsembleSolution)
    fig, ax = subplots(1,1)
    ax[:set_xlabel]("time")
    ax[:set_ylabel]("observables")
    title = ""
    for o in obs
        y, ystd = o.montec(sim)
        ax[:plot](sim[1].t, y, label=o.name)
        ax[:fill_between](sim[1].t,y.+ystd,y.-ystd,alpha=0.5)
        title *= o.name * ", "
    end
    ax[:set_title](title[1:end-2])
    ax[:legend]()
    fig[:tight_layout]()

    return fig, ax
end

function plot_observables(obs::Array{Observable},solorsim,filename::String)
    fig, ax = plot_observables(obs,solorsim)
    fig[:savefig](filename)
end

function plot_resume(sol::RODESolution)
    sol_ = extract_solution(sol)
    plot_resume(sol_)
end

function plot_resume(sim)
    plot_grid = [[Ekin,Cos2,Cos],
                 # [X,Y,Z],
                 [absX,absY,absZ],
                 [Sx,Sy,Sz],
                 [absSx,absSy,absSz],
                 [adaga,ar,ai],
                 [Sx2,Sy2,Sz2]]

    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")

    matplotlib[:rc]("axes", labelpad=1)

    fig, ax = subplots(size(plot_grid)[1],size(plot_grid[2])[1],figsize=[6.2, 8.6],sharex="col")

    for (i,row) in enumerate(plot_grid)
        color="C"*repr(i-1)
        for (j, o) in enumerate(row)
            println("elaborating "*o.name)
            y,y_std,y_q90 = expect(o,sim)
            y_q90 = hcat(y_q90...)
            tlist = sim[1].t
            ax[i,j][:set_ylabel](o.short_name)
            ax[i,j][:plot](tlist,y,color=color,label=o.short_name)
            ax[i,j][:fill_between](tlist,y_q90[1,:],y_q90[2,:],color=color,alpha=0.2)
            ax[i,j][:fill_between](tlist,y.+y_std,y.-y_std,color=color,alpha=0.5)
        end
    end
    println("-----------------------------------------")

    ax[end,1][:set_xlabel]("time")
    ax[end,2][:set_xlabel]("time")
    ax[end,3][:set_xlabel]("time")

    # ax[end,1][:set_xscale]("log")
    # ax[end,2][:set_xscale]("log")
    # ax[end,3][:set_xscale]("log")

    fig[:tight_layout](h_pad=0., w_pad=-0.)
    return fig, ax
end

function plot_resume(solorsim,filename::String)
    fig, ax = plot_resume(solorsim)
    fig[:savefig](filename)
end


########################## MISCELLANEOUS PLOTS ################################
function plot_vs_phi(sim::Array{Sol,1})
    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")

    categories, par_list = split_sim_from_par(sim,true)

    phi = Float64[]
    # y1 = (Float64[], Float64[], Array{Float64,1}[])
    y1 = (Float64[], Float64[])
    y2 = (Float64[], Float64[])
    y3 = (Float64[], Float64[])
    for x in categories
        push!(phi,angle(x[1].p.S₂))

        m,s,q = expect(adaga,x)
        # push!(y1[1],m[end])
        # push!(y1[2],s[end])
        # push!(y1[3],q[end])
        push!(y1[1],m[end])
        push!(y1[2],s[end])

        m1,s1,q1 = expect(absX,x)
        m2,s2,q2 = expect(absY,x)
        a = m1[end] ± s1[end]
        b = m2[end] ± s2[end]
        y = angle(a+im*b)
        push!(y2[1],y.val)
        push!(y2[2],y.err)

        m1,s1,q1 = expect(absar,x)
        m2,s2,q2 = expect(absai,x)
        a = m1[end] ± s1[end]
        b = m2[end] ± s2[end]
        y = angle(a+im*b)
        push!(y3[1],y.val)
        push!(y3[2],y.err)


    end

    matplotlib[:rc]("axes", labelpad=1)

    fig, ax = subplots(2,2,figsize=[6.2, 4.],sharex="col")

    ax[1,1][:set_ylabel](L"\langle a^\dag a\rangle")
    # ax[1,1][:plot](phi,y1[1],color="C0",".")
    # ax[1,1][:fill_between](phi,hcat(y1[3]...)[1,:],hcat(y1[3]...)[2,:],color="C0",alpha=0.2)
    # ax[1,1][:fill_between](phi,y1[1].+y1[2],y1[1].-y1[2],color="C0",alpha=0.5)
    ax[1,1][:errorbar](phi,y1[1],yerr=y1[2],color="C0",fmt="o")

    ax[2,1][:set_ylabel](L"\atan(\abs{Y}/\abs{X})")
    ax[2,1][:errorbar](phi,y2[1],yerr=y2[2],color="C1",fmt="o")
    ax[2,1][:set_yticks](pi/8*collect(0:4))
    ax[2,1][:set_yticklabels]([L"0",L"\pi/8",L"\pi/4",L"3\pi/8",L"\pi/2"])

    ax[2,2][:set_ylabel](L"phase of cavity field $a$")
    ax[2,2][:errorbar](phi,y3[1],yerr=y3[2],color="C2",fmt="o")
    ax[2,2][:set_yticks](pi/8*collect(0:4))
    ax[2,2][:set_yticklabels]([L"0",L"\pi/8",L"\pi/4",L"3\pi/8",L"\pi/2"])

    ax[end, 1][:set_xlabel](L"\phi")
    ax[end, 1][:set_xticks](pi/2*collect(-2:2))
    ax[end, 1][:set_xticklabels]([L"-\pi",L"-\pi/2",L"0",L"\pi/2",L"\pi"])

    ax[end, 2][:set_xlabel](L"\phi")
    ax[end, 2][:set_xticks](pi/2*collect(-2:2))
    ax[end, 2][:set_xticklabels]([L"-\pi",L"-\pi/2",L"0",L"\pi/2",L"\pi"])


    fig[:tight_layout](h_pad=0., w_pad=-0.)
    return fig, ax
end

function plot_vs_phi(solorsim,filename::String)
    fig, ax = plot_vs_phi(solorsim)
    fig[:savefig](filename)
end

function plot_vs_S(sim::Array{Sol,1})
    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")

    categories, par_list = split_sim_from_par(sim,true)

    println(par_list)

    S = Float64[]
    # y1 = (Float64[], Float64[], Array{Float64,1}[])
    y1 = (Float64[], Float64[])
    y2 = (Float64[], Float64[])
    y3 = (Float64[], Float64[])
    y4 = (Float64[], Float64[])
    for x in categories
        push!(S,abs(x[1].p.S₁))

        m,s,q = expect(adaga,x)
        push!(y1[1],m[end])
        push!(y1[2],s[end])

        m1,s1,q1 = expect(absX,x)
        m2,s2,q2 = expect(absY,x)
        a = m1[end] ± s1[end]
        b = m2[end] ± s2[end]
        y = angle(a+im*b)
        push!(y2[1],y.val)
        push!(y2[2],y.err)

        m1,s1,q1 = expect(absar,x)
        m2,s2,q2 = expect(absai,x)
        a = m1[end] ± s1[end]
        b = m2[end] ± s2[end]
        y = angle(a+im*b)
        push!(y3[1],y.val)
        push!(y3[2],y.err)

        m,s,q = expect(Cos2,x)
        push!(y4[1],m[end])
        push!(y4[2],s[end])
    end

    matplotlib[:rc]("axes", labelpad=1)

    fig, ax = subplots(2,2,figsize=[6.2, 4.],sharex="col")

    ax[1,1][:set_ylabel](L"\langle a^\dag a\rangle")
    ax[1,1][:errorbar](S,y1[1],yerr=y1[2],color="C0",fmt="o")

    ax[1,2][:set_ylabel](Cos2.short_name)
    ax[1,2][:errorbar](S,y4[1],yerr=y4[2],color="C3",fmt="o")

    ax[2,1][:set_ylabel](L"\atan(\abs{Y}/\abs{X})")
    ax[2,1][:errorbar](S,y2[1],yerr=y2[2],color="C1",fmt="o")
    ax[2,1][:set_yticks](pi/8*collect(0:4))
    ax[2,1][:set_yticklabels]([L"0",L"\pi/8",L"\pi/4",L"3\pi/8",L"\pi/2"])

    ax[2,2][:set_ylabel](L"phase of cavity field $a$")
    ax[2,2][:errorbar](S,y3[1],yerr=y3[2],color="C2",fmt="o")
    ax[2,2][:set_yticks](pi/8*collect(0:4))
    ax[2,2][:set_yticklabels]([L"0",L"\pi/8",L"\pi/4",L"3\pi/8",L"\pi/2"])

    ax[end, 1][:set_xlabel](L"|S_1|=|S_2|")
    ax[end, 2][:set_xlabel](L"|S_1|=|S_2|")

    fig[:tight_layout](h_pad=0., w_pad=-0.)
    return fig, ax
end

function plot_vs_S(solorsim,filename::String)
    fig, ax = plot_vs_S(solorsim)
    fig[:savefig](filename)
end

function plot_Ekin_vs_S(sim::Array{Sol,1})
    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")

    categories, par_list = split_sim_from_par(sim,true)

    println(par_list)

    S = Float64[]
    # y1 = (Float64[], Float64[], Array{Float64,1}[])
    y1 = (Float64[], Float64[])
    y2 = (Float64[], Float64[])
    y3 = (Float64[], Float64[])
    y4 = (Float64[], Float64[])
    for x in categories
        push!(S,abs(x[1].p.S₁))

        m,s,q = expect(Ekin,x)
        push!(y1[1],m[end])
        push!(y1[2],s[end])

    end

    # matplotlib[:rc]("axes", labelpad=1)

    fig, ax = subplots(1,1,figsize=[3.4, 2.3])

    ax[:set_ylabel](L"kinetic energy (units of $\hbar\omega_\mathrm{r}$)")
    ax[:errorbar](S,y1[1],yerr=y1[2],color="C0",fmt="o")

    ax[:set_xlabel](L"pump strength $|S_1|=|S_2|$")

    fig[:tight_layout]()
    return fig, ax
end

function plot_Ekin_vs_S(solorsim,filename::String)
    fig, ax = plot_Ekin_vs_S(solorsim)
    fig[:savefig](filename)
end

function plot_vs_SDe(sim::Array{Sol,1})
    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")

    categories, par_list = split_sim_from_par(sim,true)

    S_list = unique([abs(p.S₁) for p in par_list])
    De_list = unique([p.Δₑ for p in par_list])
    if length(De_list)*length(S_list) > length(categories)
        @info("Not all data have been plotted.")
        S_list = S_list[1:end-1]
        categories = reshape(categories[1:length(De_list)*length(S_list)],length(De_list),length(S_list))
    else
        categories = reshape(categories,length(De_list),length(S_list))
    end

    S_int = intervallize_array(S_list)
    De_int = intervallize_array(De_list)

    y1 = map(x->expect(adaga,x)[1][end],categories)

    y2 = (Float64[],Float64[])
    S = Float64[]
    De = Float64[]
    for x in categories
        push!(S,abs(x[1].p.S₁))
        push!(De,abs(x[1].p.Δₑ))

        m,s,q = expect(adaga,x)
        push!(y2[1],m[end])
        push!(y2[2],s[end])
    end

    matplotlib[:rc]("axes", labelpad=1)
    matplotlib[:rc]("image", cmap="cividis")

    fig, ax = subplots(2,2,figsize=[6.2, 4.])

    ax[1,1][:set_ylabel](L"\Delta_e")
    ax[1,1][:set_xlabel](L"|S_1|=|S_2|")
    cs11 = ax[1,1][:pcolormesh](S_int,De_int,y1)
    cbar11 = fig[:colorbar](cs11,ax=ax[1,1])
    cbar11[:set_label](L"\langle a^\dag a\rangle")

    ax[1,2][:set_ylabel](L"\Delta_e")
    ax[1,2][:set_xlabel](L"|S_1|=|S_2|")
    cs12 = ax[1,2][:contourf](S_list,De_list,y1,101)
    cbar12 = fig[:colorbar](cs12,ax=ax[1,2])
    cbar12[:set_label](L"\langle a^\dag a\rangle")

    ax[2,1][:set_ylabel](L"\Delta_e")
    ax[2,1][:set_xlabel](L"|S_1|=|S_2|")
    cs21 = ax[2,1][:pcolormesh](S_int,De_int,y1,norm=matplotlib[:colors][:LogNorm]())
    cbar21 = fig[:colorbar](cs21,ax=ax[2,1])
    cbar21[:set_label](L"\langle a^\dag a\rangle")

    ax[2,2][:set_ylabel](L"\Delta_e")
    ax[2,2][:set_xlabel](L"|S_1|=|S_2|")
    cs22 = ax[2,2][:contourf](S_list,De_list,y1,norm=matplotlib[:colors][:LogNorm]())
    cbar22 = fig[:colorbar](cs22,ax=ax[2,2])
    cbar22[:set_label](L"\langle a^\dag a\rangle")

    # ax[3,1][:set_ylabel](L"\Delta_e")
    # ax[3,1][:set_xlabel](L"|S_1|=|S_2|")
    # cs31 = ax[3,1][:tricontourf](S,De,y2[1],101,norm=matplotlib[:colors][:LogNorm]())
    # cbar31 = fig[:colorbar](cs31,ax=ax[3,1])
    # cbar31[:set_label](L"\langle a^\dag a\rangle")

    # ax[3,2][:set_ylabel](L"\Delta_e")
    # ax[3,2][:set_xlabel](L"|S_1|=|S_2|")
    # cs32 = ax[3,2][:tricontourf](S,De,y2[1],101)
    # cbar32 = fig[:colorbar](cs32,ax=ax[3,2])
    # cbar32[:set_label](L"\langle a^\dag a\rangle")


    # ax[3,1][:set_ylabel](L"phase of cavity field $a$")
    # ax[3,1][:errorbar](S,y3[1],yerr=y3[2],color="C2",fmt="o")
    # ax[3,1][:set_yticks](pi/8*collect(0:4))
    # ax[3,1][:set_yticklabels]([L"0",L"\pi/8",L"\pi/4",L"3\pi/8",L"\pi/2"])

    ax[end, 1][:set_xlabel](L"|S_1|=|S_2|")

    fig[:tight_layout](h_pad=0., w_pad=-0.)
    return fig, ax
end

function plot_vs_SDe(solorsim,filename::String)
    fig, ax = plot_vs_SDe(solorsim)
    fig[:savefig](filename)
end

########################## PLOTS FOR TALKS ################################

function plot_position(sim::Array{Sol,1},idx::Int=size(sim[1].t)[1])

    u0 = join_trajectories(sim,1)
    u1 = join_trajectories(sim,idx)

    N::Int = sum(try sol.p.N catch; sol.p[10] end for sol in sim) # for backwards compatibility
    nbins = trunc(Int, sqrt(N))

    rangex = range(0,length=trunc(Int,sqrt(N)),stop=1.)
    pdfx0 = fit(Histogram,mod2pi.(u0[1:N])./(2pi),rangex)
    pdfx1 = fit(Histogram,mod2pi.(u1[1:N])./(2pi),rangex)

    pdfx0 = normalize(pdfx0)
    pdfx1 = normalize(pdfx1)

    minp = minimum(vcat(u0[N+1:2N],u1[N+1:2N]))
    maxp = maximum(vcat(u0[N+1:2N],u1[N+1:2N]))
    rangep = range(minp,length=trunc(Int,sqrt(N)),stop=maxp)
    pdfp0 = fit(Histogram,u0[N+1:2N],rangep)
    pdfp1 = fit(Histogram,u1[N+1:2N],rangep)


    @. model(x, p) = p[1]*cos(x*p[2])+p[3]
    xdata = collect(pdfx1.edges[1][1:end-1])
    ydata = pdfx1.weights
    p0 = [600., 1., 0.]
    fitted = curve_fit(model, xdata, ydata, p0)
    println("Fitted parameters", fitted.param)
    matplotlib[:rc]("axes", labelpad=2)

    fig, ax = subplots(1,1,figsize=[4.2, 2.6])

    ax[:set_ylabel]("distribution")
    ax[:set_xlabel](L"atom position mod $\lambda_\mathrm{c}$ (units of $\lambda_\mathrm{c}$)")
    # ax[:set_xticks](pi/2*collect(0:4))
    # ax[:set_xticklabels]([L"0",L"\pi/2",L"\pi",L"3\pi/2",L"2\pi"])
    ax[:step](pdfx0.edges[1][1:end-1],pdfx0.weights,label="initial",where="post")
    ax[:step](pdfx1.edges[1][1:end-1],pdfx1.weights,label="final",where="post")
    # ax[:plot](xdata,model(xdata,fitted.param))
    ax[:legend](loc=8,ncol=1)

    fig[:tight_layout]()

    return fig, ax
end

function plot_position(sim::Array{Sol,1},filename::String)
    fig, ax = plot_position(sim)
    fig[:savefig](filename)
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

    fig[:tight_layout](h_pad=0.4, w_pad=0.)

    # sf = 0.83
    # fig, ax = subplots(2,2,figsize=[4.25*sf, 3*sf],sharex="col",sharey="row")

    # ax[1][1][:set_ylabel]("X",labelpad=-0)
    # ax[2][1][:set_ylabel]("X",labelpad=6)
    # ax[2][1][:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
    # ax[1][2][:set_ylabel]("Y",labelpad=3)
    # ax[2][2][:set_ylabel]("Y",labelpad=3)
    # ax[2][2][:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
    # # ax[1][:set_xscale]("log")
    # # ax[2][:set_xscale]("log")

    # for (i,part) in enumerate([a,d])
    #     Xs,X_std,X_q90 = expect(X,part)
    #     X_q90 = hcat(X_q90...)
    #     Ys,Y_std,Y_q90 = expect(Y,part)
    #     Y_q90 = hcat(Y_q90...)

    #     ax[i][1][:plot](tlist.+1,Xs,color="C0")
    #     # ax[1][:fill_between](tlist.+1,X_q90[1,:],X_q90[2,:],color="C0",alpha=0.2)
    #     ax[i][1][:fill_between](tlist.+1,Xs.+X_std,Xs.-X_std,color="C0",alpha=0.5)

    #     ax[i][2][:plot](tlist.+1,Ys,color="C1")
    #     # ax[2][:fill_between](tlist.+1,Y_q90[1,:],Y_q90[2,:],color="C1",alpha=0.2)
    #     ax[i][2][:fill_between](tlist.+1,Ys.+Y_std,Ys.-Y_std,color="C1",alpha=0.5)
    # end

    # fig[:tight_layout](h_pad=0.4, w_pad=0.4)

    return fig, ax
end

function plot_XY(sim::Array{Sol,1},filename::String)
    fig, ax = plot_XY(sim)
    fig[:savefig](filename)
end

function plot_Ekin(sim::Array{Sol,1})

    y,y_std,y_q90 = expect(Ekin,sim)
    y_q90 = hcat(y_q90...)
    tlist = sim[1].t

    matplotlib[:rc]("axes", labelpad=1.5)

    fig, ax = subplots(1,1,figsize=[3.4, 2.3])

    color="C3"

    ax[:set_ylabel](L"kinetic energy (units of $\hbar\omega_\mathrm{r}$)")
    ax[:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
    # ax[:set_xscale]("log")
    ax[:plot](tlist,y,color=color)
    ax[:fill_between](tlist,y_q90[1,:],y_q90[2,:],color=color,alpha=0.2)
    # ax[:fill_between](tlist.+1,y.+y_std,y.-y_std,color=color,alpha=0.5)

    fig[:tight_layout]()

    return fig, ax
end

function plot_Ekin(sim::Array{Sol,1},filename::String)
    fig, ax = plot_Ekin(sim)
    fig[:savefig](filename)
end

function plot_adaga(sim::Array{Sol,1})
    y,y_std,y_q90 = expect(adaga,sim)
    y_q90 = hcat(y_q90...)
    tlist = sim[1].t

    matplotlib[:rc]("axes", labelpad=0.5)

    fig, ax = subplots(1,1,figsize=[3.25, 2.])

    color="C2"

    ax[:set_ylabel]("cavity population")
    ax[:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
    # ax[:set_xscale]("log")
    ax[:plot](tlist.+1,y,color=color)
    # ax[:fill_between](tlist.+1,y_q90[1,:],y_q90[2,:],color=color,alpha=0.2)
    ax[:fill_between](tlist.+1,y.+y_std,y.-y_std,color=color,alpha=0.5)

    fig[:tight_layout]()

    return fig, ax
end

function plot_adaga(sim::Array{Sol,1},filename::String)
    fig, ax = plot_adaga(sim)
    fig[:savefig](filename)
end

function plot_spinspositionhisto(sim::Array{Sol,1})
    u1 = join_trajectories(sim,size(sim[1].t)[1])

    N::Int = sum(try sol.p.N catch; sol.p[10] end for sol in sim) # for backwards compatibility

    nbins = trunc(Int, ^(N,1//3))
    println("number of bins: "*string(nbins))

    matplotlib[:rc]("axes", labelpad=1)
    matplotlib[:rc]("image", cmap="viridis")

    x = mod2pi.(u1[1:N])/(2pi)

    sf = 0.78
    fig, ax = subplots(1,2,figsize=[4.25*sf, 2.0*sf],sharex=true,sharey=true)

    ax[1][:set_ylabel](Sx.short_name,labelpad=-5.)
    ax[1][:hist2d](x,u1[2N+1:3N],bins=nbins,density=true)   #,norm=matplotlib[:colors][:LogNorm]())

    ax[2][:set_ylabel](Sy.short_name,labelpad=4.)
    cs = ax[2][:hist2d](x,u1[3N+1:4N], bins=nbins,density=true)

    for i in ax
        # i[:set_xlabel](L"atom position mod $2\pi$ (units of $\lambda_\mathrm{c}$)")
        i[:set_xticks](0:0.5:1)
        # i[:set_xticklabels]([L"0",L"\pi/2",L"\pi",L"3\pi/2",L"2\pi"])
        i[:set_yticks](collect(-1:0.5:1))
        # i[:set_yticklabels]([L"-1",L"-\frac{1}{2}",L"0",L"\frac{1}{2}",L"1"])
    end

    fig[:tight_layout](h_pad=0., w_pad=-0.3)
    fig[:text](0.5, 0.02, L"atom position mod $\lambda_\mathrm{c}$ (units of $\lambda_\mathrm{c}$)", ha="center")

    fig[:subplots_adjust](right=0.86)
    cbar_ax = fig[:add_axes]([0.88, 0.2, 0.02, 0.7])
    cbar = fig[:colorbar](cs[4],cax=cbar_ax)
    # cbar[:set_label](L"\langle a^\dag a\rangle")


    return fig, ax

end

function plot_spinspositionhisto(sim::Array{Sol,1},filename::String)
    fig, ax = plot_spinspositionhisto(sim)
    fig[:savefig](filename)
end

function plot_threshold(sim::Array{Sol,1})
    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")

    categories, par_list = split_sim_from_par(sim,true)

    println(par_list)

    S = Float64[]
    # y1 = (Float64[], Float64[], Array{Float64,1}[])
    y1 = (Float64[], Float64[])
    y2 = (Float64[], Float64[])
    y3 = (Float64[], Float64[])
    y4 = (Float64[], Float64[])
    for x in categories
        push!(S,abs(x[1].p.S₁))

        m,s,q = expect(adaga,x)
        push!(y1[1],m[end])
        push!(y1[2],s[end])

        m1,s1,q1 = expect(absX,x)
        m2,s2,q2 = expect(absY,x)
        a = m1[end] ± s1[end]
        b = m2[end] ± s2[end]
        y = angle(a+im*b)
        push!(y2[1],m1[end])
        push!(y2[2],s1[end])

        a,b,c,d = split_sim(x)
        m3,s3,q3 = expect(X,a)
        m4,s4,q4 = expect(X,d)
        push!(y3[1],m3[end])
        push!(y3[2],s3[end])
        push!(y4[1],m4[end])
        push!(y4[2],s4[end])

    end

    unitS = 100.0

    # matplotlib[:rc]("axes", labelpad=1)

    y5 = (Float64[], Float64[])
    S1 = Float64[]

    for i in 1:length(S)
        if S[i]<30
            push!(S1,S[i])
            push!(y5[1],y2[1][i])
            push!(y5[2],y2[2][i])
        else
            push!(S1,S[i])
            push!(y5[1],y3[1][i])
            push!(y5[2],y3[2][i])
            push!(S1,S[i])
            push!(y5[1],y4[1][i])
            push!(y5[2],y4[2][i])
        end
    end




    fig, ax = subplots(2,1,figsize=[2.8, 3.],sharex="col")

    ax[1][:set_ylabel](L"$\langle a^\dag a\rangle$",labelpad=1)
    ax[2][:set_xlabel](L"pump strength $\Omega$ (units of $\kappa$)")
    ln1 = ax[1][:errorbar](S/unitS,y1[1],yerr=y1[2],color="C2",fmt="o",markersize=3, label="cavity population")

    ax[2][:set_ylabel](L"$\Phi$",labelpad=0)
    # ln2 = ax[2][:errorbar](S/unitS,y2[1],yerr=y2[2],color="C0",fmt="o", label="order parameter")
    ln2 = ax[2][:errorbar](S1/unitS,y5[1],yerr=y5[2],color="C0",fmt="o", markersize=3,label="order parameter")
    # ln2 = ax[1][:errorbar](S/unitS,y4[1],yerr=y4[2],color="C1",fmt="o", label="order parameter")
    # ax2.tick_params("y", colors="C0")

    # lns = [ln1,ln2]
    # labs = ["cavity population","order parameter"]

    # ax[:legend](lns,labs,loc=2)

    fig[:tight_layout](h_pad=0.)
    return fig, ax
end

function plot_threshold(solorsim,filename::String)
    fig, ax = plot_threshold(solorsim)
    fig[:savefig](filename)
end

function plot_vs_phi_talk(sim::Array{Sol,1})
    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")

    categories, par_list = split_sim_from_par(sim,true)

    phi = Float64[]
    # y1 = (Float64[], Float64[], Array{Float64,1}[])
    y1 = (Float64[], Float64[])
    y2 = (Float64[], Float64[])
    y3 = (Float64[], Float64[])
    y4 = ((Float64[], Float64[]),(Float64[], Float64[]),(Float64[], Float64[]),(Float64[], Float64[]))
    for x in categories
        push!(phi,angle(x[1].p.S₂))

        m,s,q = expect(adaga,x)
        # push!(y1[1],m[end])
        # push!(y1[2],s[end])
        # push!(y1[3],q[end])
        push!(y1[1],m[end])
        push!(y1[2],s[end])

        m1,s1,q1 = expect(absX,x)
        m2,s2,q2 = expect(absY,x)
        a = m1[end] ± s1[end]
        b = m2[end] ± s2[end]
        y = angle(a+im*b)
        push!(y2[1],y.val)
        push!(y2[2],y.err)

        m1,s1,q1 = expect(absar,x)
        m2,s2,q2 = expect(absai,x)
        a = m1[end] ± s1[end]
        b = m2[end] ± s2[end]
        y = angle(a+im*b)
        push!(y3[1],y.val)
        push!(y3[2],y.err)


        for (i,c) in enumerate(split_sim(x))
            if isempty(c)
                push!(y4[i][1],NaN)
                push!(y4[i][2],NaN)
            else
                m1,s1,q1 = expect(X,c)
                m2,s2,q2 = expect(Y,c)
                a = m1[end] ± s1[end]
                b = m2[end] ± s2[end]
                y = angle(a+im*b)
                push!(y4[i][1],y.val)
                push!(y4[i][2],y.err)
            end
        end
    end
    phi = -phi/2

    fig, ax = subplots(1,2,figsize=[4., 1.8],sharey=true)

    ax[1][:set_ylabel](L"spin phase $\phi_s$",labelpad=1)
    ax[1][:errorbar](phi,y2[1],yerr=y2[2],color="C1",fmt="o")
    ax[1][:set_yticks](pi/8*collect(0:4))
    ax[1][:set_yticklabels]([L"0",L"\pi/8",L"\pi/4",L"3\pi/8",L"\pi/2"])

    # for c in y4
    #     ax[1][:errorbar](phi,c[1],yerr=c[2],fmt="o")
    # end


    ax[2][:set_ylabel]("phase of cavity field",labelpad=6)
    ax[2][:errorbar](phi,y3[1],yerr=y3[2],color="C2",fmt="o")
    # ax[2][:set_yticks](pi/8*collect(0:4))
    # ax[2][:set_yticklabels]([L"0",L"\pi/8",L"\pi/4",L"3\pi/8",L"\pi/2"])

    ax[1][:set_xlabel](L"relative lasers phase $\phi$")
    ax[1][:set_xticks](pi/4*collect(-2:2))
    ax[1][:set_xticklabels]([L"-\pi/2",L"-\pi/4",L"0",L"\pi/4",L"\pi/2"])

    ax[2][:set_xlabel](L"relative lasers phase $\phi$")
    ax[2][:set_xticks](pi/4*collect(-2:2))
    ax[2][:set_xticklabels]([L"-\pi/2",L"-\pi/4",L"0",L"\pi/4",L"\pi/2"])

    fig[:tight_layout](h_pad=0,w_pad=0.5)
    return fig, ax
end

function plot_vs_phi_talk(solorsim,filename::String)
    fig, ax = plot_vs_phi_talk(solorsim)
    fig[:savefig](filename)
end



#end
