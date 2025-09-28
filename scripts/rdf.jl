#================================================================
Evaluate radial distribution function from ABM simulations.

For each simulation, all timepoints are pooled together.
This is valid assuming the simulations are equilibrated.
The distribution of bacteria as a function of distance
from the source is evaluated using kernel density estimation.
The KDE of a chemotactic system is then divided by the
KDE of the equivalent non-chemotactic system, providing the RDF.
================================================================#

using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using DataFrames
using StatsBase
using KernelDensity

function p_r(dataset::AbstractDataFrame; bandwidth=25)
    InterpKDE(kde(dataset.r; bandwidth))
end

function rdf(
    kde_chemo::InterpKDE,
    kde_random::InterpKDE,
    r::AbstractVector,
)
    pdf(kde_chemo, r) ./ pdf(kde_random, r)
end

datasets = collect_datasets(datadir("brumley"); rinclude=[r"abm"])
gdata = groupby(datasets, [:R, :U, :mot, :L]) # don't group Cs here

npoints = 100 # points for the rdf sampling
L = 1e3 # L is always 1 mm
r = range(1, L/2; length=npoints)
radial_distributions = vcat(map(collect(gdata)) do dataset
    gdf = groupby(dataset, :Cs)
    data_random = gdf[(0,)] # take random walk as reference
    k0 = p_r(data_random)
    # for each Cs estimate the RDF by comparing to random walk
    vcat(map(keys(gdf)) do key
        data_chemotax = gdf[key]
        k = p_r(data_chemotax)
        DataFrame(;
            R=dataset.R[1],
            U=dataset.U[1],
            mot=dataset.mot[1],
            Cs=key.Cs,
            k=k,
            g=[rdf(k, k0, r)]
        )
    end...)
end...)
wsave(
    datadir("brumley", "rdf.jld2"),
    Dict("rdf" => radial_distributions)
)
