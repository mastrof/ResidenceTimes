using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using CSV
using DataFrames

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

filenames = readdir(datadir("comm_7classes_tracks"))

for fname in filenames
    println(fname)
    df = CSV.read(
        datadir("comm_7classes_tracks", fname),
        DataFrame;
        select=[:time, :id, :c]
    )
    gdf = groupby(df, :id)
    ids, es = eachcol(CSV.read(datadir("comm_7classes_ids", replace(fname, "comm"=>"ids1-99")), DataFrame))
    q = sortperm(es)
    lo = q[1:1000]
    hi = q[1001:end]
    hotspot_threshold = 0.063 # μM
    encounters, residencetimes = hotspot_encounters(gdf, hotspot_threshold)
    restime_per_encounter = sum.(residencetimes) ./ encounters
    #
    CSV.write(
        datadir("comm_7classes_proc", replace(fname, "comm"=>"encounters")),
        DataFrame(; id=ids, encounters=encounters, residencetime=restime_per_encounter)
    )
    #
    fout = replace(fname, "comm"=>"residencetimes", "csv"=>"txt")
    open(datadir("comm_7classes_proc", fout), "w") do io
        for i in eachindex(ids)
            id = ids[i]
            ts = residencetimes[i]
            write(io, join([id; ts], " ")*"\n")
        end
    end
end
