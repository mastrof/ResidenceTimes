using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using CSV
using DataFrames
using StatsBase

dist(a, b) = hypot((a.-b)...)
function inrange(comm, xs, ys, zs)
    idxs = findall(
        p -> (xs[1]<p.x<xs[2]) && (ys[1]<p.y<ys[2]) &&(zs[1]<p.z<zs[2]),
        eachrow(comm)
    )
    comm[idxs, :]
end
function mindist(gdf::GroupedDataFrame, comm, dt=10)
    vcat((mindist(g, comm, dt) for g in gdf)...)
end
function mindist(g::SubDataFrame, comm, dt=10)
    ts = 1:dt:nrow(g)
    nt = length(ts) # number of time points
    d = fill(Inf, nt) # distance from closest phytoplankton
    idp = zeros(Int, nt) # id of closest phytoplankton
    # retain only phytoplankton within explored volume
    w = (-1, +1) .* 100
    xs = extrema(g.x) .+ w
    ys = extrema(g.y) .+ w
    zs = extrema(g.z) .+ w
    subcomm = inrange(comm, xs, ys, zs)
    if isempty(subcomm)
        return DataFrame(; time=empty(g.time), id=empty(g.id), d=empty(d), idp=empty(idp))
    end
    np = nrow(subcomm)
    for j in eachindex(ts)
        t = ts[j]
        for i in 1:np
            b = Vector(g[t,[:x, :y, :z]])
            p = Vector(subcomm[i,[:x, :y, :z]])
            new_d = dist(b, p)
            if new_d < d[j]
                d[j] = new_d
                idp[j] = subcomm[i, :id]
            end
        end
    end
    DataFrame(; time=g.time[ts], id=g.id[ts], d, idp)
end

filenames = readdir(datadir("comm_7classes_tracks"))
for fname in filenames
    println(fname)
    ks = ["Aphy", "L", "α"]
    parameters = parse_savename(fname)[2]
    comm_parameters = Dict(ks .=> get.(Ref(parameters), ks, nothing))
    comm = CSV.read(datadir(savename("phytoplankton", comm_parameters, "csv")), DataFrame)
    df = CSV.read(
        datadir("comm_7classes_tracks", fname),
        DataFrame;
    )
    gdf = groupby(df, :id)
    dt = 10 # sample every 10 frames (=1s)
    outdf = mindist(gdf, comm, dt)
    CSV.write(
        datadir("comm_7classes_proc", savename("dist", parameters, "csv")),
        outdf
    )
end
