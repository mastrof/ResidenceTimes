using DrWatson
using ResidenceTimes
using CSV
using DataFrames
using KernelDensity
using LaTeXStrings
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
        labelsize=18
    ),
    Lines=(
        linewidth=4,
    )
)

df = collect_datasets(datadir("comm_7classes");
    # U=10,
    Aphy=80000,
    PER=0.3,
    α=1.0,
)
gdf = groupby(df, [:U, :PER, :Aphy, :α])

fig = Figure(size=(1400,400))
# panel A - mean
let
    ax = Axis(fig[1,1];
        xlabel="Mean cell exposure (μM)",
        ylabel="probability density",
        xticks=logrange(1e-4, 1; length=5)
    )
    x = logrange((extrema(df.c) .+ 1e-4)...; length=80)
    for g in gdf
        U = first(g.U)
        PER = first(g.PER)
        Aphy = first(g.Aphy)
        α = first(g.α)
        c = [mean(a.c) for a in groupby(g, :id)] # μM
        # @show mean(c) median(c) std(c)
        K = kde(c)
        y = pdf(K, x)
        for i in eachindex(y)
            if y[i] <= 1e-2
                y[i] = NaN
            end
        end
        label = join((U, PER, Aphy, α), ", ")
        lines!(ax, x, y; label)
    end
    axislegend(ax, "U(μm/s), PER, Nphy, α"; position=:rt)
    ax.yscale = log10
    ax.xscale = log10
end
# panel B - median
let
    ax = Axis(fig[1,2];
        xlabel="Median cell exposure (μM)",
        # ylabel="probability density",
        xticks=logrange(1e-4, 1; length=5)
    )
    x = logrange((extrema(df.c) .+ 1e-4)...; length=80)
    for g in gdf
        U = first(g.U)
        PER = first(g.PER)
        Aphy = first(g.Aphy)
        α = first(g.α)
        c = [median(a.c) for a in groupby(g, :id)] # μM
        # @show mean(c) median(c) std(c)
        K = kde(c)
        y = pdf(K, x)
        for i in eachindex(y)
            if y[i] <= 1e-2
                y[i] = NaN
            end
        end
        label = join((U, PER, Aphy, α), ", ")
        lines!(ax, x, y; label)
    end
    ax.yscale = log10
    ax.xscale = log10
end
# panel C - minimum
let
    ax = Axis(fig[1,3];
        xlabel="Minimum cell exposure (μM)",
        # ylabel="probability density",
        xticks=logrange(1e-4, 1; length=5)
    )
    x = logrange((extrema(df.c) .+ 1e-4)...; length=80)
    for g in gdf
        U = first(g.U)
        PER = first(g.PER)
        Aphy = first(g.Aphy)
        α = first(g.α)
        c = [minimum(a.c) for a in groupby(g, :id)] # μM
        # @show mean(c) median(c) std(c)
        K = kde(c)
        y = pdf(K, x)
        for i in eachindex(y)
            if y[i] <= 1e-2
                y[i] = NaN
            end
        end
        label = join((U, PER, Aphy, α), ", ")
        lines!(ax, x, y; label)
    end
    ax.yscale = log10
    ax.xscale = log10
end
fig
