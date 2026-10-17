using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using CSV
using DataFrames
using StatsBase
using Distributions
using KernelDensity
using Random
using StaticArrays
using Agents: ContinuousSpace
import BubbleBath
using GLMakie
set_theme!(theme_latexfonts();
    # fontsize=22,
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
        titlefont=:italic,
        labelsize=18,
    ),
    Lines=(
        linewidth=1,
    )
)
it(s) = rich(s; font=:italic)

function threshold_crossings(x::AbstractVector{<:Real}, threshold::Real)
    above = x .> threshold
    n = length(x)
    durations = Int64[]
    i = 1
    while i <= n
        if above[i]
            j = i
            while j <= n && above[j]
                j += 1
            end
            push!(durations, j - i)
            i = j
        else
            i += 1
        end
    end
    return length(durations), durations
end
function hotspot_encounters(gdf, threshold)
    encounters = zeros(Int, gdf.ngroups)
    residencetimes = [Int[] for _ in gdf]
    for i in eachindex(encounters)
        g = gdf[i]
        counts, durations = threshold_crossings(g.c, threshold)
        encounters[i] = counts
        residencetimes[i] = durations
    end
    return encounters, residencetimes
end

function commstructure(L, α, Aphy; Rmin=0.65, Rmax=10.0, rng=Xoshiro(1))
    space = ContinuousSpace(fill(L, SVector{3}); periodic=true)
    comm = generate_community(space, α, Aphy, Rmin, Rmax, 0; rng)[1]
    # periodic images
    offsets = filter(
        x -> sum(abs.(x))==1,
        collect(Iterators.product(-1:1, -1:1, -1:1))
    )
    comm_images = [deepcopy(comm) for _ in offsets]
    for i in eachindex(offsets)
        d = offsets[i]
        for j in eachindex(comm_images[i])
            p = comm_images[i][j].pos
            r = comm_images[i][j].radius
            comm_images[i][j] = BubbleBath.Sphere(p .+ (d.*L), r)
        end
    end
    [comm; comm_images...]
end
function spheresinrange(comm, xs, ys, zs)
    idxs = findall(
        s -> (xs[1]<s.pos[1]<xs[2]) && (ys[1]<s.pos[2]<ys[2]) &&(zs[1]<s.pos[3]<zs[2]),
        comm
    )
    comm[idxs]
end
dist(a, b) = hypot((a.-b)...)


filenames = readdir(datadir("comm_7classes_tracks"))
fname = filenames[3]
α = parse_savename(fname)[2]["α"]
L = parse_savename(fname)[2]["L"]
Aphy = parse_savename(fname)[2]["Aphy"]
comm = commstructure(L, α, Aphy)
df = CSV.read(
    datadir("comm_7classes_tracks", fname),
    DataFrame;
)
gdf = groupby(df, :id)
ids, es = eachcol(CSV.read(datadir("comm_7classes_ids", replace(fname, "comm"=>"ids1-99")), DataFrame))
q = sortperm(es)
lo = q[1:1000]
hi = q[1001:end]
hotspot_threshold = 0.063 # μM
encounters, residencetimes = hotspot_encounters(gdf, hotspot_threshold)
restime_per_encounter = sum.(residencetimes) ./ encounters

i = 1711#1823
let
    hotspot_threshold = 0.063 # μM
    g = gdf[q[i]]
    s = 1:10:27000
    ts = s .* 0.1 ./ 60
    fig = Figure(size=(467,680))
    a = Axis3(fig[1,1];
        xlabelvisible=false,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        zticklabelsvisible=false,
        ylabelvisible=false,
        zlabelvisible=false,
        aspect=:data
    )
    b = Axis(fig[2,1];
        ylabel=rich(it("E"), " (μM)"),
        xticklabelsvisible=false,
        yscale=log10,
        yticks=logrange(1e-3, 1; length=4),
    )
    c = Axis(fig[3,1];
        ylabel=rich(it("d"), " (μm)"),
        xlabel="time (min)",
        yticks=0:60:300,
    )
    color = [c > hotspot_threshold ? :purple2 : :black for c in g.c[s]]
    # trajectory extrema
    xs = extrema(g.x[s])
    ys = extrema(g.y[s])
    zs = extrema(g.z[s])
    xlims!(a, xs .+ (-10, 10))
    ylims!(a, ys .+ (-10, 10))
    zlims!(a, zs .+ (-10, 10))
    a.xticks = xs[1]:50:xs[2]
    a.yticks = ys[1]:50:ys[2]
    a.zticks = zs[1]:50:zs[2]
    w = (-25, 25)
    phytoplankton = spheresinrange(comm, xs.+w, ys.+w, zs.+w)
    for p in phytoplankton
        mesh!(a, Sphere(Point3f(p.pos), p.radius*6);
            color=:green,
        )
    end
    # trajectory
    lines!(a, g.x[s], g.y[s], g.z[s]; color, alpha=0.5)
    # exposure timeseries
    lines!(b, ts, g.c[s]; color)
    # distance from closest phytoplankton
    hlines!(b, hotspot_threshold; linestyle=:dash, color=:purple2, linewidth=0.5)
    r = [
        minimum(
            dist([g.x[t], g.y[t], g.z[t]], p.pos)
            for p in phytoplankton
        )
        for t in s
    ]
    lines!(c, ts, r; color)
    rowsize!(fig.layout, 2, Relative(1/6))
    rowsize!(fig.layout, 3, Relative(1/6))
    fig
end

save(
    plotsdir("trajectories",
        savename("track_id=$(i)", parse_savename(fname)[2], "png")
    ),
    current_figure();
    px_per_unit=3,
    backend=GLMakie
)
