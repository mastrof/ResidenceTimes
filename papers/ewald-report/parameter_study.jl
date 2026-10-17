# Accuracy and cost of the Ewald community field as a function of its numerical parameters.
# Produces the numbers reported in ewald-report.tex.
#
#   include("papers/ewald-report/parameter_study.jl")   # in a project-activated session
using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes, MicrobeAgents, StaticArrays, LinearAlgebra, Random, Statistics
using CellListMap: map_pairwise!
using Agents: allids, step!
using NonuniformFFTs: HalfSupport
using ResidenceTimes: EwaldFar, far_gradient, ewald_near, point_strength,
    _concentration, _gradient, generate_community

const γ = 150.0

## test community and probe points
function test_community(; L=800.0, Aphy=3e5, α=1.0, rng=Xoshiro(3))
    space = ContinuousSpace(fill(L, SVector{3}); periodic=true)
    cfg, Cs = generate_community(space, α, Aphy, 0.65, 10.0, 0.5; rng)
    P = SVector{3}.(getproperty.(cfg, :pos))
    R = getproperty.(cfg, :radius)
    (; P, R, Cs, Q=point_strength.(R, Cs, γ), L)
end

wrap(d, L) = d - L * round(d / L)
minimage(v, L) = SVector(wrap(v[1], L), wrap(v[2], L), wrap(v[3], L))

# half of the probes 0.5-5 μm outside a random cell surface, half uniform outside all cells
function probes(src; n=400, rng=Xoshiro(11))
    (; P, R, L) = src
    outside(x) = all(j -> norm(minimage(x - P[j], L)) > R[j], eachindex(P))
    surface = map(1:n÷2) do _
        j = rand(rng, eachindex(P))
        mod.(P[j] + normalize(randn(rng, SVector{3,Float64})) * (R[j] + 0.5 + 4.5rand(rng)), L)
    end
    bulk = map(1:n÷2) do _
        x = SVector{3}(rand(rng, 3) .* L)
        while !outside(x)
            x = SVector{3}(rand(rng, 3) .* L)
        end
        x
    end
    (; surface, bulk)
end

## reference: direct sum over all sources and periodic images closer than rmax
function brute_force(x, src; rmax=18γ)
    (; P, R, Cs, L) = src
    n = ceil(Int, rmax / L)
    c, g = 0.0, zero(SVector{3,Float64})
    for j in eachindex(P)
        d0 = minimage(x - P[j], L)
        for i1 in -n:n, i2 in -n:n, i3 in -n:n
            dvec = d0 + L * SVector(i1, i2, i3)
            d = norm(dvec)
            d < rmax || continue
            c += _concentration(d, R[j], Cs[j], γ)
            g += _gradient(d, dvec, R[j], Cs[j], γ)
        end
    end
    c, g
end

## Ewald evaluation with explicit parameters (same kernels as the model)
function ewald(x, src, far, σ, rc)
    (; P, R, Q, L) = src
    c, g = far(x), far_gradient(far, x)
    for j in eachindex(P)
        dvec = minimage(x - P[j], L)
        d = norm(dvec)
        d < rc || continue
        G, dG = ewald_near(max(d, R[j]), γ, σ)
        c += Q[j] * G
        g += Q[j] * dG / d * dvec
    end
    c, g
end

function errors(src, xs, ref; σ=40.0, h=σ/2.4, rc=3.5σ, m=8)
    far = EwaldFar(src.P, src.Q, src.L, γ, σ; h, m=HalfSupport(m))
    eC = Float64[]; eG = Float64[]; eGs = Float64[]
    for (x, (c0, g0)) in zip(xs, ref)
        c, g = ewald(x, src, far, σ, rc)
        push!(eC, abs(c - c0) / c0)
        push!(eG, norm(g - g0) / norm(g0))
        push!(eGs, norm(g - g0) / (c0 / γ)) # relative to the natural gradient scale C/γ
    end
    (; eC_med=median(eC), eC_max=maximum(eC), eG_med=median(eG), eG_max=maximum(eG),
        eGs_max=maximum(eGs), h=far.h)
end

function report(label, rows)
    println("\n# ", label)
    println(rpad("param", 12), rpad("h (μm)", 10), "C med      C max      ∇C med     ∇C max     ∇C max (C/γ)")
    for (p, e) in rows
        println(rpad(string(p), 12), rpad(string(round(e.h; digits=2)), 10),
            join((rpad(string(round(v; sigdigits=2)), 11) for v in (e.eC_med, e.eC_max, e.eG_med, e.eG_max, e.eGs_max))))
    end
end

function accuracy_study()
    src = test_community()
    pr = probes(src)
    xs = vcat(pr.surface, pr.bulk)
    println("sources: ", length(src.P), ", probes: ", length(xs))
    ref = [brute_force(x, src) for x in xs]
    # convergence of the reference itself
    ref24 = [brute_force(x, src; rmax=24γ) for x in xs[1:40]]
    println("reference 18γ vs 24γ: max rel. diff C = ",
        maximum(abs(a[1] - b[1]) / b[1] for (a, b) in zip(ref, ref24)),
        ", ∇C = ", maximum(norm(a[2] - b[2]) / norm(b[2]) for (a, b) in zip(ref, ref24)))
    # near-surface probes only, for comparison
    nsurf = length(pr.surface)
    e_surf = errors(src, xs[1:nsurf], ref[1:nsurf])
    e_bulk = errors(src, xs[nsurf+1:end], ref[nsurf+1:end])
    report("default parameters, surface vs bulk probes", ["surface" => e_surf, "bulk" => e_bulk])

    report("mesh spacing σ/h (σ = 40, rc = 3.5σ, m = 8)",
        [r => errors(src, xs, ref; h=40.0 / r) for r in (1.2, 1.6, 2.0, 2.4, 3.2, 4.0, 4.8)])
    report("near cutoff rc/σ (σ = 40, h = σ/4.8, m = 8)",
        [a => errors(src, xs, ref; rc=a * 40.0, h=40.0 / 4.8) for a in (2.0, 2.5, 3.0, 3.5, 4.0)])
    report("NUFFT half-support m (σ = 40, h = σ/4.8, rc = 3.5σ)",
        [m => errors(src, xs, ref; m, h=40.0 / 4.8) for m in (2, 3, 4, 6, 8)])
    report("splitting length σ (h = σ/2.4, rc = 3.5σ, m = 8)",
        [σ => errors(src, xs, ref; σ) for σ in (20.0, 30.0, 40.0, 60.0, 80.0, 100.0)])
    return nothing
end

## cost of the field evaluation in a production-size system
function timing(; L, Aphy, σ, n=100_000, α=0.85, reps=5)
    t_setup = @elapsed model = setup_abm_community(; n, L, Aphy, α, σ, rng=Xoshiro(1))
    far = model.ewald_far
    t_far_setup = @elapsed EwaldFar(model.neighborlist.ypositions,
        point_strength.(model.phytoplankton_radii, model.phytoplankton_leakage, γ), L, γ, σ)
    nlist = model.neighborlist
    ids = collect(allids(model))
    near!() = map_pairwise!((x, y, i, j, d2, out) -> comm_out!(x, y, i, j, d2, out, model), nlist)
    function far!()
        out = nlist.measurements
        for i in ids
            x = position(model[i])
            out.c[i] += far(x)
            out.dc[i] += far_gradient(far, x)
        end
    end
    near!(); far!()
    t_near = minimum(@elapsed(near!()) for _ in 1:reps)
    t_far = minimum(@elapsed(far!()) for _ in 1:reps)
    step!(model, 1)
    t_step = minimum(@elapsed(step!(model, 1)) for _ in 1:reps)
    M = round(Int, L / far.h)
    (; L, σ, nphyto=length(model.phytoplankton_radii), M, mesh_MB=8M^3 / 1e6,
        t_setup, t_far_setup, t_near, t_far, t_step)
end

function timing_study()
    rows = [timing(; L, Aphy, σ) for (L, Aphy) in ((1000.0, 3e5), (2150.0, 3e4))
            for σ in (20.0, 30.0, 40.0, 60.0, 80.0)]
    println("\n# cost at 100 000 bacteria (seconds; near part serial)")
    println("L       σ     Nphy  M     mesh MB  far setup  near/step  far/step  full step")
    for r in rows
        println(rpad(Int(r.L), 8), rpad(Int(r.σ), 6), rpad(r.nphyto, 6), rpad(r.M, 6),
            rpad(round(r.mesh_MB; digits=1), 9), rpad(round(r.t_far_setup; digits=2), 11),
            rpad(round(r.t_near; digits=3), 11), rpad(round(r.t_far; digits=3), 10),
            round(r.t_step; digits=3))
    end
    return rows
end
