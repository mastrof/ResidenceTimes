using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using DataFrames
using MicrobeAgents
using Random
using StaticArrays
using GLMakie

U = 10
PER = 0.3
α = 0.75
mot = "RRF"
λ = 2.2
Aphy = 300_000
L = 1500
Cb = 0.0
dt = 0.1
γ = 90.0
n = 0
rng = Xoshiro(123)
model = setup_abm_community(; mot, dt, U, λ, L, α, Aphy, PER, Cb, γ, n, rng)

positions = model.neighborlist.ypositions
radii = model.phytoplankton_radii
spheres = [
    Sphere(Point(pos), r) for (pos, r) in zip(positions, radii)
]

l = 300
xs = ys = zs = range(0, L; length=l)
c = zeros(l, l, l)
p0 = SVector{3,Float64}(fill(L/2, 3))
r0 = 100
idxs = eachindex(xs)
Threads.@threads for ijk in collect(Iterators.product(idxs, idxs, idxs))
    i, j, k = ijk
    x = xs[i]
    y = ys[j]
    z = zs[k]
    pos = SVector{3,Float64}(x,y,z)
    c[i,j,k] = concentration(model)(pos, model)
    r = distance(pos, p0, model)
    if r >= L/2-r0
        h = (r-(L/2-r0))/(1.5*r0)
        c[i,j,k] *= exp(-h)
    end
end

colormap = cgrad(
    [:black, :purple, :magenta, :pink],
    [0.03, 0.1, 0.25, 0.95]
)
let
    fig = Figure(backgroundcolor=:black, size=(600,600))
    ax = Axis3(fig[1,1], viewmode=:fit, aspect=:data)
    hidespines!(ax)
    hidedecorations!(ax)
    volume!(ax, 0..L, 0..L, 0..L, c;
        colormap,
    )
    fig
end

let
    fig = Figure(backgroundcolor=:black, size=(900,900))
    ax = Axis3(fig[1,1]; viewmode=:fit)
    hidespines!(ax)
    hidedecorations!(ax)
    ax.aspect = :data
    volume!(ax, 0..L, 0..L, 0..L, c;
        colormap,
    )
    display(fig)
    # @async while events(fig).window_open[]
    #     ax.azimuth[] += 0.005
    #     sleep(1/60)
    # end
    n_frames = 200
    θ_step = 2π / n_frames
    record(fig, "community.mp4", 1:n_frames; framerate=15, compression=10) do frame
        ax.azimuth[] += θ_step
    end
end
