using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using CSV
using DataFrames
using Random
using StaticArrays
using Agents: ContinuousSpace
import BubbleBath

# reproduce community used in experiments using same seed
# include also periodic images at each boundary of the domain
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
    c = [comm; comm_images...]
    ids = repeat(eachindex(comm), length(offsets)+1)
    getx(s) = s.pos[1]
    gety(s) = s.pos[2]
    getz(s) = s.pos[3]
    getr(s) = s.radius
    DataFrame(;
        id=ids,
        x=getx.(c),
        y=gety.(c),
        z=getz.(c),
        R=getr.(c),
    )
end

filenames = readdir(datadir("comm_7classes_tracks"))
dicts = getindex.(parse_savename.(filenames), 2)
ks = ["Aphy", "L", "α"]
pars = unique(map(d -> get.(Ref(d), ks, nothing), dicts))

for i in eachindex(pars)
    p = pars[i]
    d = Dict(ks .=> p)
    comm = commstructure(d["L"], d["α"], d["Aphy"])
    CSV.write(
        datadir(savename("phytoplankton", d, "csv")),
        comm
    )
end
