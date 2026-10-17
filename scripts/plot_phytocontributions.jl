using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using CSV
using DataFrames
using LinearAlgebra
using StatsBase
using KernelDensity
using LaTeXStrings
using GLMakie
set_theme!(theme_latexfonts();
    fontsize=18,
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
        titlefont=:regular,
        labelsize=18,
    ),
    Lines=(
        linewidth=6,
    )
)
it(s) = rich(s; font=:italic)

function individual_encounters(d, trk; threshold=0.063)
    idx = findall(>=(threshold), trk.c)
    DataFrame(; c=trk[idx,:c], idp=d[idx,:idp])
end

function idptor(g::AbstractDataFrame, comm)
    [comm[id,:R] for id in g.idp]
end
function idptor!(g::AbstractDataFrame, comm)
    R = idptor(g, comm)
    g.R = R
    g
end

function stackedbarplot!(ax, dfs, i)
    df = vcat(dfs...)
    sums = combine(groupby(df, :R), :c => sum)
    sort!(sums, :R)
    n = nrow(sums)
    stackidx = 1:n
    tot = sum(sums.c_sum)
    sums.c_sum ./= tot
    color = categorical_colors(Makie.wong_colors(), n)
    barplot!(ax, fill(i, n), sums.c_sum;
        stack=stackidx,
        color=color[stackidx],
    )
end

filenames = filter(
    s -> startswith(s, "dist"),
    readdir(datadir("comm_7classes_proc"))
)[3:4]
comm = let
    ks = ["Aphy", "L", "α"]
    parameters = parse_savename(filenames[1])[2]
    comm_parameters = Dict(ks .=> get.(Ref(parameters), ks, nothing))
    CSV.read(datadir(savename("phytoplankton", comm_parameters, "csv")), DataFrame)
end

data = map(filenames) do fname
    ids, es = eachcol(CSV.read(
        datadir(
            "comm_7classes_ids",
            replace(fname, "dist"=>"ids1-99")
        ),
        DataFrame
    ))
    q = sortperm(es)
    lo = q[1:1000]
    hi = q[1001:end]
    dist = CSV.read(
        datadir("comm_7classes_proc", fname),
        DataFrame;
    )
    n = length(unique(dist.id))
    gdist = groupby(dist, :id)
    dt = 10 # distances are downsampled every dt steps
    tracks = CSV.read(
        datadir("comm_7classes_tracks", replace(fname, "dist"=>"comm")),
        DataFrame;
    )
    tmin, tmax = extrema(tracks.time)
    subset!(tracks, :time => ByRow(t -> t ∈ tmin:dt:tmax))
    gtracks = groupby(tracks, :id)
    E_lo = [
        individual_encounters(gdist[(i,)], gtracks[(i,)])
        for i in ids[lo] if i in unique(dist.id)
    ]
    E_hi = [
        individual_encounters(gdist[(i,)], gtracks[(i,)])
        for i in ids[hi] if i in unique(dist.id)
    ]
    Dict(
        :lo => E_lo,
        :hi => E_hi,
    )
end
df_lowU = data[1]
df_hiU = data[2]

idptor!.(df_lowU[:hi], Ref(comm))
idptor!.(df_lowU[:lo], Ref(comm))
idptor!.(df_hiU[:hi], Ref(comm))
idptor!.(df_hiU[:lo], Ref(comm))

let
    fig = Figure(size=(600,500))
    Rs = sort(unique(comm.R))
    xlabels = [
        "worst\n10 μm/s",
        "best\n10 μm/s",
        "worst\n40 μm/s",
        "best\n40 μm/s",
    ]
    ax = Axis(fig[1,1];
        xticks=1:4,
        yticks=0:0.25:1,
        xtickformat = xs -> map(i -> xlabels[i], eachindex(xs)),
        ylabel="Relative contribution to exposure",
    )
    xlims!(ax, 0.5, 4.5)
    ylims!(ax, 0, 1)
    stackedbarplot!(ax, df_lowU[:lo], 1)
    stackedbarplot!(ax, df_lowU[:hi], 2)
    stackedbarplot!(ax, df_hiU[:lo], 3)
    stackedbarplot!(ax, df_hiU[:hi], 4)
    colors = categorical_colors(Makie.wong_colors(), length(Rs))
    Legend(fig[1,0],
        [PolyElement(color=c) for c in colors],
        string.(Rs) .* " μm";
    )
    fig
end
