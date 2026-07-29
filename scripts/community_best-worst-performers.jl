using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using CSV
using DataFrames
using StatsBase

function find_ids_extremal_performers(x, lo=1, hi=99)
    thrs_low, thrs_hi = percentile(x, [lo, hi])
    findall(id -> id<=thrs_low || id >= thrs_hi, x)
end

allfiles = readdir(datadir("comm_7classes"))
for fname in allfiles
    fout = replace(fname, "comm" => "ids_$(lo)-$(hi)")
    isfile(datadir("comm_7classes_ids", fout)) && continue
    df = CSV.read(datadir("comm_7classes", fname), DataFrame)
    gdf = groupby(df, :id)
    mean_exposures = [mean(g.c) for g in gdf]
    lo, hi = 1, 99
    ids = find_ids_extremal_performers(mean_exposures, lo, hi)
    CSV.write(
        datadir("comm_7classes_ids", fout),
        (id=ids, e=mean_exposures[ids])
    )
end
