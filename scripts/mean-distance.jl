using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using DataFrames
using StatsBase
using LinearAlgebra
using KernelDensity
using GLMakie

function mindist(data::AbstractDataFrame)
    [minimum(g.r) for g in groupby(data, :id)]
end

function entrances(data::AbstractDataFrame, R::Real)
    gdf = groupby(data, :id)
    counts = zeros(Int, gdf.ngroups)
    for (i, g) in enumerate(gdf)
        for t in 2:nrow(g)
            if g.r[t-1] > R && g.r[t] <= R
                counts[i] += 1
            end
        end
    end
    counts
end

function residence(data::AbstractDataFrame, R::Real)
    gdf = groupby(data, :id)
    n = gdf.ngroups
    τ = [Int[] for _ in 1:n]
    for (i, g) in enumerate(gdf)
        for t in 2:nrow(g)
            if g.r[t-1] > R && g.r[t] <= R
                push!(τ[i], 1)
            elseif g.r[t-1] <= R && g.r[t] <= R
                if isempty(τ[i]) # inside at time 0
                    push!(τ[i], 1)
                end
                τ[i][end] += 1
            end
        end
    end
    idxs = findall(!isempty, τ)
    # idxs .=> τ[idxs]
    vcat(τ[idxs]...)
end

datasets = collect_datasets(
    datadir("brumley");
    R=1,
    U=10,
    λ=2.2,
    mot="RRF",
)
# filter!(:Cs => C -> C .!= 1, datasets)
sort!(datasets, [:mot, :Cs])
gdata = groupby(datasets, [:Cs, :R, :U, :mot])

dist = map(collect(gdata)) do dataset
    gdf = groupby(dataset, :time)
    [mean(g.r) for g in gdf[2:end]] # skip t=0
end
let
    fig = Figure()
    ax = Axis(fig[1,1];
        xlabel="time (min)",
        ylabel="mean distance from source (μm)",
    )
    n = gdata.ngroups
    for i in 1:n
        dt = gdata[i].dt[1]
        t = unique(gdata[i].time)[2:end] .* (dt/60)
        label = string(gdata[i].Cs[1])
        lines!(ax, t, dist[i]; label)
    end
    axislegend(ax, "Cs"; position=:rc)
    fig
end

let
    fig = Figure()
    ax = Axis(fig[1,1];
        xlabel="initial distance from source (μm)",
        ylabel="relative gross exposure"
    )
    for i in 1:gdata.ngroups-1
        gdf = groupby(gdata[i], :id)
        p = [Point2f(first(g.r), sum(g.c)) for g in gdf]
        scatter!(ax, p; alpha=0.2)
    end
    ax.xscale = log10
    ax.yscale = log10
    fig
end

let
    fig = Figure()
    ax = Axis(fig[1,1];
        xlabel="gross displacement (μm)",
        ylabel="relative gross exposure"
    )
    for i in 1:gdata.ngroups-1
        gdf = groupby(gdata[i], :id)
        p = [Point2f(sum(abs.(diff(g.r))), sum(g.c)) for g in gdf]
        scatter!(ax, p; alpha=0.2)
    end
    ax.xscale = log10
    ax.yscale = log10
    fig
end

let R = 35
    entr = map(dataset -> entrances(dataset, R), collect(gdata))
    fig = Figure()
    ax = Axis(fig[1,1];
        xlabel="# entrances in $(R)μm phycosphere",
        ylabel="PDF",
    )
    for i in eachindex(entr)
        h = normalize(fit(Histogram, entr[i], -0.5:10.5); mode=:pdf)
        x = midpoints(h.edges[1])
        y = h.weights
        j = findall(y .> 0)
        scatterlines!(ax, x[j], y[j])
    end
    xlims!(ax, -0.25, 10.25)
    ax.xticks = 0:10
    ax.yscale = log10
    fig
end

let R = 35
    res = map(dataset -> residence(dataset, R), collect(gdata))
    fig = Figure()
    ax = Axis(fig[1,1];
        xlabel="residence time within $(R)μm",
        ylabel="PDF"
    )
    for i in eachindex(res)[1:end-1]
        h = normalize(fit(Histogram, res[i], 0.5:50.5); mode=:pdf)
        x = midpoints(h.edges[1])
        y = h.weights
        j = findall(y .> 0)
        scatterlines!(ax, x[j], y[j])
    end
    ax.yscale = log10
    fig
end

# let
#     fig = Figure()
#     ax = Axis(fig[1,1];
#         xlabel="minimum distance from source (μm)",
#         ylabel="PDF"
#     )
#     n = gdata.ngroups
#     for i in 1:n
#         d = mindist(gdata[i])
#         k = kde(d; boundary=(0,500), bandwidth=25)
#         lines!(ax, k.x, k.density)
#     end
#     fig
# end
