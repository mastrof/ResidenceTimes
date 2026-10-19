using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using ResidenceTimes: leaked_concentration
using CSV
using DataFrames
using Random
using StaticArrays
using Agents: ContinuousSpace
import BubbleBath

# Reproduce the phytoplankton community of the community simulations (same generator,
# same seed, same defaults as `setup_abm_community`) and write, for every configuration:
#   1. phytoplankton_<params>.csv   position, radius and leakage of every cell
#      (the 6 face-adjacent periodic images are included, flagged by `image=true`)
#   2. cfield_<params>.jld2         the interpolated concentration field; load it with
#          field = load_community_field(file);  field(x);  field_gradient(field, x)

# must match the defaults of setup_abm_community and the `abm_community.jl` script
# (PER is read from the simulation filenames)
const Rmin, Rmax = 0.65, 10.0
const Cb = 0.0
const γ = 150.0
const σ = 40.0
const SEED = 1

function commstructure(L, α, Aphy, PER; rng=Xoshiro(SEED))
    space = ContinuousSpace(fill(L, SVector{3}); periodic=true)
    generate_community(space, α, Aphy, Rmin, Rmax, PER; rng)
end

function commtable(spatialcfg, leakage, L)
    offsets = filter(x -> sum(abs.(x)) == 1, collect(Iterators.product(-1:1, -1:1, -1:1)))
    shifts = [(0, 0, 0); offsets]
    ids = eachindex(spatialcfg)
    reduce(vcat, map(shifts) do d
        DataFrame(;
            id=collect(ids),
            x=[s.pos[1] + d[1]*L for s in spatialcfg],
            y=[s.pos[2] + d[2]*L for s in spatialcfg],
            z=[s.pos[3] + d[3]*L for s in spatialcfg],
            R=[s.radius for s in spatialcfg],
            Cs=leakage,
            image=fill(d != (0, 0, 0), length(ids)),
        )
    end)
end

filenames = readdir(datadir("comm_7classes_exposure"))
dicts = getindex.(parse_savename.(filenames), 2)
ks = ["Aphy", "L", "α", "PER"]
pars = unique(map(d -> get.(Ref(d), ks, nothing), dicts))

for p in pars
    d = Dict(ks .=> p)
    L = d["L"]
    spatialcfg, leakage = commstructure(L, d["α"], d["Aphy"], d["PER"])
    CSV.write(datadir(savename("phytoplankton", d, "csv")), commtable(spatialcfg, leakage, L))
    field = CommunityField(
        getproperty.(spatialcfg, :pos), getproperty.(spatialcfg, :radius), leakage,
        L, γ, σ, Cb
    )
    save_community_field(datadir(savename("cfield", d, "jld2")), field)
    @info "wrote community" d length(spatialcfg)
end
