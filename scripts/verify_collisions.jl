#================================================================
Verify that no bacterium ever penetrates a phytoplankton, in both
the single-source and community models, and plot a few trajectories
grazing a sphere. Run via the julia-mcp `julia_eval` tool.
================================================================#

using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using MicrobeAgents
using Agents
using StaticArrays
using LinearAlgebra
using Random
using DataFrames
using CairoMakie

const TOL = 1e-4   # μm; clamped cells sit at the surface within rounding
const WARMUP = 50  # steps to let any random initial overlaps eject before asserting

"No swimmer centre is inside any sphere (centre distance ≥ R + microbe radius)."
function assert_no_penetration(model, centres, radii)
    for a in allagents(model)
        for j in eachindex(radii)
            d = distance(a, centres[j], model)
            @assert d ≥ radii[j] + radius(a) - TOL "penetration: id=$(a.id) d=$d R=$(radii[j])"
        end
    end
end

# --- single-source model ---
model = setup_abm(; n=2_000, L=1000.0, R=10.0, U=40.0, mot="RRF", dt=0.1, Cs=1.0)
centres = [chemoattractant(model).origin]
radii = [chemoattractant(model).radius]
for _ in 1:WARMUP            # eject any cell randomly initialized inside the sphere
    step!(model, 1)
end
for _ in 1:500
    step!(model, 1)
    assert_no_penetration(model, centres, radii)
end
println("single-source: no penetration over 500 steps ✓")

# --- community model ---
cmodel = setup_abm_community(; n=5_000, L=1000.0, Aphy=1e6, mot="RRF",
                             U=40.0, dt=0.1, rng=Xoshiro(7))
ccentres = cmodel.neighborlist.ypositions
cradii = cmodel.phytoplankton_radii
for _ in 1:WARMUP            # eject any cell randomly initialized inside a phytoplankton
    step!(cmodel, 1)
end
for _ in 1:500
    step!(cmodel, 1)
    assert_no_penetration(cmodel, ccentres, cradii)
end
println("community: no penetration over 500 steps ✓")

# --- trajectory plot: a few swimmers grazing the central sphere ---
rng = Xoshiro(9)
tmodel = setup_abm(;
    n=120, L=400.0, R=25.0, U=40.0,
    mot="RRF", dt=0.1, Cs=1.0, λ=0.0,
    Drot=0.04,
    rng
)
origin = chemoattractant(tmodel).origin
# seed swimmers heading at the sphere from one side
for a in allagents(tmodel)
    newpos = origin .- SVector(60.0, 0.0, 0.0) .+
        SVector(0.0, 5.0*(rand(abmrng(tmodel))-0.5), 0.0) .+
        SVector(0.0, 0.0, 5.0*(rand(abmrng(tmodel))-0.5))
    move_agent!(a, newpos, tmodel)
    a.vel = SVector(1.0, 0.0, 0.0)
    a.speed = 40.0
end
tracks = [Point3f[] for _ in 1:nagents(tmodel)]
for _ in 1:40
    step!(tmodel, 1)
    assert_no_penetration(tmodel, [origin], [25.0])
    for (k, a) in enumerate(allagents(tmodel))
        push!(tracks[k], Point3f(position(a)))
    end
end
fig = Figure()
ax = Axis3(fig[1,1]; aspect=:data,
    xlabel="x (μm)", ylabel="y (μm)", zlabel="z (μm)",
    title="Swimmers grazing a sphere (R=25 μm)"
)
mesh!(ax, Sphere(Point3f(origin), 25.0); color=(:black, 0.2))
# arc!(ax, Point2f(origin[1], origin[2]), 20.0, 0, 2π; color=:black)
for t in tracks
    scatterlines!(ax, t; color=(:steelblue, 0.5))
end
save(plotsdir("verify_collisions.png"), fig)
println("trajectory plot written to ", plotsdir("verify_collisions.png"))
