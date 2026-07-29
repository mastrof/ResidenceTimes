using DrWatson
using ResidenceTimes
using CSV
using DataFrames
using KernelDensity
using StatsBase
using LaTeXStrings
using Printf
using Colors
using GLMakie
set_theme!(theme_latexfonts();
    fontsize=22,
    Axis=(
        xgridvisible=false,
        ygridvisible=false,
        xticksize=-5,
        yticksize=-5,
        xticksmirrored=true,
        yticksmirrored=true,
        titlesize=18,
    ),
    Legend=(
        framecolor=:transparent,
        backgroundcolor=RGBAf(1,1,1,0.6),
        patchsize=(40, 16),
        titlesize=18,
        titlefont=:italic,
        labelsize=18,
    ),
    Lines=(
        linewidth=6,
    )
)

it(s) = rich(s; font=:italic)

function plot_exposurepdf!(ax, g, x, reducer, color)
    U = first(g.U)
    PER = first(g.PER)
    Aphy = first(g.Aphy)
    α = first(g.α)
    c = [reducer(a.c) for a in groupby(g, :id)] # nm
    K = kde(c)
    y = pdf(K, x)
    for i in eachindex(y)
        if y[i] <= 1e-5#2
            y[i] = NaN
        end
    end
    # label = join((U, PER, Aphy, α), ", ")
    label = string(U) * " μm/s"
    lines!(ax, x, y; label, color)
end

function plot_cumulexposure!(ax, g, reducer, color)
    ps = 1:2:99
    Cs = percentile([reducer(a.c) for a in groupby(g,:id)], ps)
    Cs_print = percentile([reducer(a.c) for a in groupby(g,:id)], [15, 95])
    @printf "/** %sμm/s %s **/\n" string(first(g.U)) string(nameof(reducer))
    @printf "bottom 15%%: %f\ntop 5%%: %f\n" Cs_print[1] Cs_print[2]
    lines!(ax, Cs, 1 .- ps ./ 100; color)
end

df = collect_datasets(datadir("comm_7classes");
    # U=40,
    Aphy=30000,
    PER=0.3,
    α=1.0,
)
gdf = groupby(df, [:U, :PER, :Aphy, :α])
# x = logrange(1e-4, maximum(df.c); length=80)
df.c .*= 1e3 # μM → nM
x = logrange(1e-1, maximum(df.c); length=80)

fig = Figure(size=(1400,400))
colors = [
    "#B133DE",
    "#D57D19",
    "#205C80",
]
# panel A - mean
let
    ax = Axis(fig[1:2,1];
        xlabel="Mean cell exposure (nM)",
        ylabel="Probability density",
        # xticks=logrange(1e-4, 1; length=5)
        xticks=logrange(1e-1, 1e3; length=5)
    )
    # xlims!(ax, 5e-5, 6e-1)
    xlims!(ax, 5e-2, 6e2)
    text!(ax, 0.75, 0.75+0.075;
        text="Mean:",
        space=:relative,
        fontsize=18,
        color=:black
    )
    for (i,g) in enumerate(gdf)
        plot_exposurepdf!(ax, g, x, mean, colors[i])
        avg = mean(mean(a.c) for a in groupby(g, :id))
        text!(ax, 0.75, 0.75-0.075*(i-1);
            text=string(round(avg; sigdigits=2)) * " nM",
            space=:relative,
            fontsize=18,
            color=colors[i]
        )
    end
    # axislegend(ax, "U(μm/s), PER, Nphy, α"; position=:lb)
    axislegend(ax, "U"; position=:cb)
    ax.yscale = log10
    ax.xscale = log10
    ylims!(ax, 8e-6, 0.9)
    text!(ax, 0.9, 0.88; text="A",
        space=:relative,
        font=:bold,
        fontsize=28
    )
end
# panel B - median
let
    ax = Axis(fig[1:2,2];
        xlabel="Median cell exposure (nM)",
        # ylabel="probability density",
        # xticks=logrange(1e-4, 1; length=5),
        xticks=logrange(1e-1, 1e3; length=5),
    )
    # xlims!(ax, 5e-5, 6e-1)
    xlims!(ax, 5e-2, 6e2)
    text!(ax, 0.75, 0.75+0.075;
        text="Mean:",
        space=:relative,
        fontsize=18,
        color=:black
    )
    for (i,g) in enumerate(gdf)
        plot_exposurepdf!(ax, g, x, median, colors[i])
        avg = mean(median(a.c) for a in groupby(g, :id))
        text!(ax, 0.75, 0.75-0.075*(i-1);
            text=string(round(avg; sigdigits=2)) * " nM",
            space=:relative,
            fontsize=18,
            color=colors[i]
        )
    end
    ax.yscale = log10
    ax.xscale = log10
    ylims!(ax, 8e-6, 0.9)
    text!(ax, 0.9, 0.88; text="B",
        space=:relative,
        font=:bold,
        fontsize=28
    )
end
# # panel C and D - mean and median percentiles
let
    ax1 = Axis(fig[1,3];
        # xlabel="C (nM)",
        ylabel=rich("P(mean > ", it("c"), ")"),
        # xticks=logrange(1e-4, 1; length=5),
        xticks=logrange(1e-1, 1e3; length=5),
        xticklabelsvisible=false,
        yticks=logrange(1e-2, 1; length=3),
    )
    # xlims!(ax1, 5e-5, 6e-1)
    # xlims!(ax1, 5e-2, 6e2)
    xlims!(ax1, 9e-2, 6e1)
    ax2 = Axis(fig[2,3];
        xlabel=rich(it("c"), " (nM)"),
        ylabel=rich("P(median > ", it("c"), ")"),
        # xticks=logrange(1e-4, 1; length=5),
        xticks=logrange(1e-1, 1e3; length=5),
        yticks=logrange(1e-2, 1; length=3),
    )
    # xlims!(ax2, 5e-5, 6e-1)
    xlims!(ax2, 9e-2, 6e1)
    for (i,g) in enumerate(gdf)
        plot_cumulexposure!(ax1, g, mean, colors[i])
        plot_cumulexposure!(ax2, g, median, colors[i])
    end
    ax1.xscale = log10
    ax2.xscale = log10
    ax1.yscale = log10
    ax2.yscale = log10
    text!(ax1, 0.9, 0.75; text="C",
        space=:relative,
        font=:bold,
        fontsize=28
    )
    text!(ax2, 0.9, 0.75; text="D",
        space=:relative,
        font=:bold,
        fontsize=28
    )
end
fig
