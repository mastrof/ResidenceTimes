using DrWatson
push!(LOAD_PATH, projectdir("analysis")) # Makie lives in the analysis workspace project
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
        linewidth=4,
    )
)

it(s) = rich(s; font=:italic)

function plot_exposurepdf!(ax, g, x, color)
    U = first(g.U)
    e = g.exposure # μM s
    K = kde(log.(e))
    y = pdf(K, log.(x)) ./ x
    label = string(U) * " μm/s"
    lines!(ax, x, y; label, color)
end
function plot_exposurepdf_mot!(ax, g, x, color)
    mot = first(g.mot)
    e = g.exposure # μM s
    K = kde(log.(e))
    y = pdf(K, log.(x)) ./ x
    label = string(mot)
    lines!(ax, x, y; label, color)
end

df_oligotrophic_rrf = collect_datasets(datadir("comm_7classes_exposure");
    Aphy=30000,
    PER=0.5,
    α=1.0,
    mot="RRF",
)
df_eutrophic_rrf = collect_datasets(datadir("comm_7classes_exposure");
    Aphy=300000,
    PER=0.5,
    α=1.0,
    mot="RRF",
)
df_slow = collect_datasets(datadir("comm_7classes_exposure");
    Aphy=30000,
    PER=0.5,
    α=1.0,
    U=10,
)
sort!(df_slow, :mot)
df_fast = collect_datasets(datadir("comm_7classes_exposure");
    Aphy=30000,
    PER=0.5,
    α=1.0,
    U=40,
)
sort!(df_fast, :mot)

fig = Figure(size=(800,800))
colors = [
    "#B133DE",
    "#D57D19",
    "#205C80",
]
# panel A - oligotrophic
let df=df_oligotrophic_rrf
    gdf = groupby(df, [:U, :PER, :Aphy, :α])
    x = logrange(extrema(df.exposure)...; length=100)
    ax = Axis(fig[1,1];
        xlabel="Nutrient exposure (μM⋅s)",
        ylabel="Probability density",
        xticks=LogTicks(-1:3),
        yticks=LogTicks(-5:1),
    )
    text!(ax, 0.52, 0.75+0.075;
        text="Mean:",
        space=:relative,
        fontsize=18,
        color=:black
    )
    text!(ax, 0.78, 0.75+0.075;
        text="99.9%:",
        space=:relative,
        fontsize=18,
        color=:black
    )
    for (i,g) in enumerate(gdf)
        plot_exposurepdf!(ax, g, x, colors[i])
        avg = mean(g.exposure)
        text!(ax, 0.55, 0.75-0.075*(i-1);
            text=string(round(avg; sigdigits=2)),# * " μM⋅s",
            space=:relative,
            fontsize=18,
            color=colors[i]
        )
        p99 = percentile(g.exposure, 99.9)
        text!(ax, 0.80, 0.75-0.075*(i-1);
            text=string(round(p99; sigdigits=2)),# * " μM⋅s",
            space=:relative,
            fontsize=18,
            color=colors[i]
        )
    end
    # axislegend(ax, "U(μm/s), PER, Nphy, α"; position=:lb)
    axislegend(ax, "U"; position=:cb)
    ax.yscale = log10
    ax.xscale = log10
    xlims!(ax, 8e-2, 205)
    ylims!(ax, 1e-6, 3)
    text!(ax, 0.9, 0.88; text="A",
        space=:relative,
        font=:bold,
        fontsize=28
    )
end
# panel B - eutrophic
let df=df_eutrophic_rrf
    gdf = groupby(df, [:U, :PER, :Aphy, :α])
    x = logrange(extrema(df.exposure)...; length=100)
    ax = Axis(fig[1,2];
        xlabel="Nutrient exposure (μM⋅s)",
        # ylabel="Probability density",
        xticks=LogTicks(-1:3),
        yticks=LogTicks(-5:1),
    )
    text!(ax, 0.52, 0.75+0.075;
        text="Mean:",
        space=:relative,
        fontsize=18,
        color=:black
    )
    text!(ax, 0.78, 0.75+0.075;
        text="99.9%:",
        space=:relative,
        fontsize=18,
        color=:black
    )
    for (i,g) in enumerate(gdf)
        plot_exposurepdf!(ax, g, x, colors[i])
        avg = mean(g.exposure)
        text!(ax, 0.55, 0.75-0.075*(i-1);
            text=string(round(avg; sigdigits=2)),# * " μM⋅s",
            space=:relative,
            fontsize=18,
            color=colors[i]
        )
        p99 = percentile(g.exposure, 99.9)
        text!(ax, 0.80, 0.75-0.075*(i-1);
            text=string(round(p99; sigdigits=2)),# * " μM⋅s",
            space=:relative,
            fontsize=18,
            color=colors[i]
        )
    end
    ax.yscale = log10
    ax.xscale = log10
    xlims!(ax, 3, 205)
    ylims!(ax, 1e-6, 3)
    text!(ax, 0.9, 0.88; text="B",
        space=:relative,
        font=:bold,
        fontsize=28
    )
end
# panel C - all motilities, U=10
let df=df_slow
    gdf = groupby(df, [:mot, :PER, :Aphy, :α])
    x = logrange(extrema(df.exposure)...; length=100)
    ax = Axis(fig[2,1];
        xlabel="Nutrient exposure (μM⋅s)",
        ylabel="Probability density",
        xticks=LogTicks(-1:3),
        yticks=LogTicks(-5:1),
    )
    text!(ax, 0.52, 0.75+0.075;
        text="Mean:",
        space=:relative,
        fontsize=18,
        color=:black
    )
    text!(ax, 0.78, 0.75+0.075;
        text="99.9%:",
        space=:relative,
        fontsize=18,
        color=:black
    )
    for (i,g) in enumerate(gdf)
        plot_exposurepdf_mot!(ax, g, x, colors[i])
        avg = mean(g.exposure)
        text!(ax, 0.55, 0.75-0.075*(i-1);
            text=string(round(avg; sigdigits=2)),# * " μM⋅s",
            space=:relative,
            fontsize=18,
            color=colors[i]
        )
        p99 = percentile(g.exposure, 99.9)
        text!(ax, 0.80, 0.75-0.075*(i-1);
            text=string(round(p99; sigdigits=2)),# * " μM⋅s",
            space=:relative,
            fontsize=18,
            color=colors[i]
        )
    end
    axislegend(ax, "Motility"; position=:cb, titlefont=:regular)
    ax.yscale = log10
    ax.xscale = log10
    xlims!(ax, 8e-2, 205)
    ylims!(ax, 1e-6, 3)
    text!(ax, 0.9, 0.88; text="C",
        space=:relative,
        font=:bold,
        fontsize=28
    )
end
# panel D - all motilities, U=40
let df=df_fast
    gdf = groupby(df, [:mot, :PER, :Aphy, :α])
    x = logrange(extrema(df.exposure)...; length=100)
    ax = Axis(fig[2,2];
        xlabel="Nutrient exposure (μM⋅s)",
        # ylabel="Probability density",
        xticks=LogTicks(-1:3),
        yticks=LogTicks(-5:1),
    )
    text!(ax, 0.52, 0.75+0.075;
        text="Mean:",
        space=:relative,
        fontsize=18,
        color=:black
    )
    text!(ax, 0.78, 0.75+0.075;
        text="99.9%:",
        space=:relative,
        fontsize=18,
        color=:black
    )
    for (i,g) in enumerate(gdf)
        plot_exposurepdf_mot!(ax, g, x, colors[i])
        avg = mean(g.exposure)
        text!(ax, 0.55, 0.75-0.075*(i-1);
            text=string(round(avg; sigdigits=2)),# * " μM⋅s",
            space=:relative,
            fontsize=18,
            color=colors[i]
        )
        p99 = percentile(g.exposure, 99.9)
        text!(ax, 0.80, 0.75-0.075*(i-1);
            text=string(round(p99; sigdigits=2)),# * " μM⋅s",
            space=:relative,
            fontsize=18,
            color=colors[i]
        )
    end
    ax.yscale = log10
    ax.xscale = log10
    xlims!(ax, 8e-2, 205)
    ylims!(ax, 1e-6, 3)
    text!(ax, 0.9, 0.88; text="D",
        space=:relative,
        font=:bold,
        fontsize=28
    )
end
fig
