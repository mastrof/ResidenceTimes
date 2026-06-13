using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using DataFrames
using StatsBase
using KernelDensity
using LinearAlgebra
using HypothesisTests
using LaTeXStrings
using GLMakie
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

function fithist(data::AbstractVector; nbins=50)
    h = normalize(fit(Histogram, data; nbins))
    x = midpoints(h.edges[1])
    y = h.weights
    Point2f.(x, y)
end

function makeplot(exposure, R, Cs, mot, λ, Us)
    fig = Figure()
    ax = Axis(fig[1,1];
        xlabel="log relative exposure",
        ylabel="PDF"
    )
    for i in eachindex(exposure)
        e = log10.(exposure[i])
        scatter!(ax, fithist(e; nbins=40))
        vlines!(ax, [mean(e)])
        x = range(-5.1,0.5;length=400)
        label = string(Us[i])
        # ks = ApproximateTwoSampleKSTest(exposure[1], exposure[i])
        # sd = pvalue(ks; tail=:right) < 0.05
        # if sd
        #     label *= "  ⃰"
        # end
        lines!(ax, x, pdf(kde(e; bandwidth=0.09), x);
            label
        )
    end
    ylims!(ax, 1e-4, 2)
    ax.title = LaTeXString(
        "C_\\mathrm{S}=$(Cs)\\,μM\\ R=$(R)\\,μm\\ λ=$(λ)\\,s^{-1}\\ \\mathrm{$(mot)}"
    )
    ax.yscale = log10
    axislegend(ax, L"U\,\mathrm{(μm/s)}";
        position=:rt,
        framevisible=true,
        framecolor=RGBAf(1,1,1,0.6),
        backgroundcolor=RGBAf(1,1,1,0.6)
    )
    fig
end

Rs = [1, 3, 10]
Css = [0, 0.1, 1.0]
mots = ["RR", "RRF", "RT"]
λs = [2.2, 0.5]
for R in Rs, Cs in Css, mot in mots, λ in λs
    datasets = collect_datasets(datadir("brumley");
        R=R,
        Cs=Cs,
        mot=mot,
        λ=λ,
    )
    sort!(datasets, [:mot, :U])
    gdata = groupby(datasets, [:Cs, :R, :U, :λ, :mot])
    exposure = map(collect(gdata)) do dataset
        gdf = groupby(dataset, :id)
        [sum(g.c) / nrow(g) for g in gdf]
    end
    Us = sort(unique(datasets.U))
    fig = makeplot(exposure, R, Cs, mot, λ, Us)
    save(plotsdir(savename("exposure", (@strdict R Cs λ mot), "png")), fig)
end

# let
#     fig = Figure()
#     ax = Axis(fig[1,1];
#         xlabel="time",
#         ylabel="relative exposure"
#     )
#     for n in 1:gdata.ngroups
#         g = gdata[n]
#         ts = unique(g.time)[2:end]
#         e0 = [mean(gdata[1][g.time .== t, :c]) for t in ts]
#         e = [mean(g[g.time .== t, :c]) for t in ts]
#         lines!(ax, ts, e ./ e0)
#     end
#     # ax.yscale = log10
#     fig
# end
