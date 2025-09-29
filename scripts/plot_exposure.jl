using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using DataFrames
using StatsBase
using LinearAlgebra
using HypothesisTests
using LaTeXStrings
using GLMakie

function fithist(data::AbstractVector; nbins=50)
    h = normalize(fit(Histogram, data; nbins))
    x = midpoints(h.edges[1])
    y = h.weights
    Point2f.(x, y)
end

function makeplot(exposure, R, U, mot, Cs)
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
        label = string(Cs[i])
        ks = ApproximateTwoSampleKSTest(exposure[1], exposure[i])
        sd = pvalue(ks; tail=:right) < 0.05
        if sd
            label *= "  ⃰"
        end
        lines!(ax, x, pdf(kde(e; bandwidth=0.09), x);
            label
        )
    end
    ylims!(ax, 1e-4, 2)
    ax.title = LaTeXString("U=$(U)\\,μm/s\\ R=$(R)\\,μm\\ \\mathrm{$(mot)}")
    ax.yscale = log10
    axislegend(ax, L"C_\mathrm{S}\,\mathrm{(μM)}";
        position=:lt,
        framevisible=true,
        framecolor=RGBAf(1,1,1,0.6),
        backgroundcolor=RGBAf(1,1,1,0.6)
    )
    fig
end

Rs = [1, 3, 10]
Us = [25, 45]
mots = ["RR", "RRF", "RT"]
for R in Rs, U in Us, mot in mots
    datasets = collect_datasets(datadir("brumley");
        R=R,
        U=U,
        mot=mot
    )
    sort!(datasets, [:mot, :Cs])
    gdata = groupby(datasets, [:Cs, :R, :U, :mot])
    exposure = map(collect(gdata)) do dataset
        gdf = groupby(dataset, :id)
        [sum(g.c) / nrow(g) for g in gdf]
    end
    Cs = sort(unique(datasets.Cs))
    fig = makeplot(exposure, R, U, mot, Cs)
    save(plotsdir(savename("exposure", (@strdict R U mot), "png")), fig)
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
