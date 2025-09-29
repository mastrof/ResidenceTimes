using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using DataFrames
using StatsBase
using KernelDensity
using JLD2
using GLMakie
using LaTeXStrings
set_theme!(theme_latexfonts();
    fonts=(
        title="Arial Bold",
        annotation="Arial",
        legendtitle="Arial Bold",
        legend="Arial",
    ),
    fontsize=22,
    Text=(
        font=:annotation,
    ),
    Axis=(
        xgridvisible=false,
        ygridvisible=false,
        xticksize=-5,
        yticksize=-5,
        xticksmirrored=true,
        yticksmirrored=true,
        titlefont=:title,
    ),
    Legend=(
        framevisible=false,
        backgroundcolor=nothing,
        patchsize=(40, 16),
        titlefont=:legendtitle,
        # labelfont=:legend,
    ),
    Lines=(
        linewidth=4,
    )
)

function makeplot(radial_distributions, U, R; L=1e3, npoints=100)
    w = 400
    r = range(1, L/2; length=npoints)
    gdf = groupby(
        subset(
            radial_distributions,
            :U => u -> u .== U,
            :R => r -> r .== R,
            :Cs => Cs -> Cs .> 0;
            view = true
        ),
        [:mot, :Cs]
    )
    cs = sort(unique(radial_distributions.Cs))[2:end] # skip 0
    mots = unique(radial_distributions.mot)
    nc = length(cs)
    fig = Figure(size=(nc*w, w))
    xlabel = L"r\ \mathrm{(μm)}"
    ylabel = L"g(r)"
    axs = [Axis(fig[1,i]; xlabel, ylabel) for i in 1:nc]
    linkxaxes!(axs...)
    xlims!.(axs, 0, 190)
    for i in 1:nc
        Cs = cs[i]
        ax = axs[i]
        for mot in mots
            g = gdf[(mot, Cs)]
            lines!(ax, r, g.g[1]; label=mot)
        end
        hlines!(ax, 1.0; color=:black, linestyle=:dash)
        if i == 1
            axislegend(ax, "Motility"; position=:rt)
        end
        ax.title = LaTeXString("U=$U\\,\\mathrm{μm/s}\\ R=$R\\,\\mathrm{μm}\\ C_\\mathrm{S}=$(Cs)\\,\\mathrm{μM}")
    end
    fig
end

radial_distributions = jldopen(datadir("brumley", "rdf.jld2"), "r")["rdf"]
L = 1e3
r = range(1, L/2; length=100)

# R=1, U=25, 3 motility patterns at all Cs
Us = unique(radial_distributions.U)
Rs = unique(radial_distributions.R)
for U in Us, R in Rs
    fig = makeplot(radial_distributions, U, R)
    save(plotsdir("rdf_R=$(R)_U=$(U).png"), fig)
end
