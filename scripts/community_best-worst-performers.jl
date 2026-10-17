using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using CSV
using DataFrames
using Random

# k lowest and k highest exposures, plus k random ids as a reference group;
# the random group is drawn from all ids, so it may overlap with the tails
function select_ids(df, k; rng=Xoshiro(1))
    order = sortperm(df.exposure)
    bottom = df.id[order[1:k]]
    top = df.id[order[end-k+1:end]]
    random = df.id[randperm(rng, nrow(df))[1:k]]
    ids = sort(union(bottom, top, random))
    selection = DataFrame(id=ids)
    leftjoin!(selection, select(df, :id, :exposure); on=:id)
    sort!(selection, :id)
    transform!(selection, :id => ByRow(in(Set(bottom))) => :bottom)
    transform!(selection, :id => ByRow(in(Set(top))) => :top)
    transform!(selection, :id => ByRow(in(Set(random))) => :random)
    selection
end

on_cluster = haskey(ENV, "SCRATCH")
k = 1000
for fname in readdir(datadir("comm_7classes_exposure"; on_cluster))
    fout = datadir("comm_7classes_selection", replace(fname, r"^comm" => "sel"); on_cluster)
    isfile(fout) && continue
    df = CSV.read(datadir("comm_7classes_exposure", fname; on_cluster), DataFrame)
    wsave(fout, select_ids(df, k))
end
