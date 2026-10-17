using ResidenceTimes
using MicrobeAgents
using Agents: step!, allagents
using StaticArrays
using LinearAlgebra
using Random
using Test

@testset "ResidenceTimes collisions" begin
    @testset "single-source collision step" begin
        # one swimmer, central sphere R=10 at the domain centre
        model = setup_abm(; n=1, L=1000.0, R=10.0, U=40.0, mot="RT", dt=0.1, Cs=1.0)
        origin = chemoattractant(model).origin
        R_eff = chemoattractant(model).radius + radius(model[1]) # sphere + microbe radius

        # head-on: start 12 μm from centre, step 4 μm inward → pinned at the surface
        a = model[1]
        a.pos = origin .- SVector(12.0, 0.0, 0.0)
        a.vel = SVector(1.0, 0.0, 0.0)
        a.speed = 40.0
        move_step!(a, model)
        @test isapprox(distance(a, origin, model), R_eff; atol=1e-6)
        @test a.speed == 0.0                                          # halted

        # a pinned cell does not leave the surface by rotational diffusion alone
        pos = position(a)
        for _ in 1:20
            move_step!(a, model)
        end
        @test position(a) == pos

        # on the surface, pointing outward → moves away, keeps speed
        b = model[1]
        b.pos = origin .- SVector(R_eff, 0.0, 0.0)
        b.vel = SVector(-1.0, 0.0, 0.0)
        b.speed = 40.0
        move_step!(b, model)
        @test b.speed == 40.0
        @test distance(b, origin, model) > R_eff
    end

    @testset "setup_abm uses collision stepping" begin
        # with ghost stepping the swimmer would move the full 4 μm and end up inside
        model = setup_abm(; n=1, L=1000.0, R=10.0, U=40.0, mot="RR", dt=0.1, Cs=1.0)
        origin = chemoattractant(model).origin
        a = model[1]
        a.pos = origin .- SVector(12.0, 0.0, 0.0)
        a.vel = SVector(1.0, 0.0, 0.0)
        a.speed = 40.0
        step!(model, 1)
        R_eff = chemoattractant(model).radius + radius(a)
        @test distance(a, origin, model) ≥ R_eff - 1e-6   # never penetrated the sphere
    end

    @testset "collision across periodic boundary" begin
        space = ContinuousSpace((100.0, 100.0, 100.0); periodic=true)
        C = SVector(2.0, 50.0, 50.0)       # sphere centre near the x=0 edge
        sphere = Obstacles(space, [ObstacleSphere(C, 1.0)]; contact=HaltUntilTumble())
        model = StandardABM(ResidenceTimes.microbe_type("RT"), space, 0.1;
            properties=Dict(:obstacles => (sphere,)), container=Vector)
        a = add_agent!(SVector(99.0, 50.0, 50.0), model; # 3 μm from C across the boundary
            motility=ResidenceTimes.make_motility("RT", 40.0, 2.2),
            radius=0.5, behaviors=(Brumley(),),
            vel=SVector(1.0, 0.0, 0.0), speed=40.0,       # step 4 μm, wrapping toward C
        )
        move_step!(a, model)
        @test isapprox(distance(a, C, model), 1.5; atol=1e-6)  # sphere + microbe radius
        @test a.speed == 0.0
    end

    @testset "community surface collision" begin
        # L must exceed 2*field_cutoff = 2*3.5σ = 280 μm (CellListMap unit-cell rule)
        model = setup_abm_community(; n=1, L=800.0, Aphy=1e6, mot="RRF",
                                    U=40.0, dt=0.1, rng=Xoshiro(1))
        @test length(model.phytoplankton_radii) ≥ 1

        # aim the single swimmer straight at the first phytoplankton, 5 μm outside it,
        # fast enough (step 8 μm) to overshoot the surface without clamping
        P = model.neighborlist.ypositions[1]
        R = model.phytoplankton_radii[1]
        a = model[1]
        dir = normalize(distancevector(position(a), P, model))
        a.pos = P .- dir .* (R + 5.0)
        a.vel = dir
        a.speed = 80.0
        step!(model, 1)
        R_eff = R + radius(a)
        @test isapprox(distance(a, P, model), R_eff; atol=1e-4)
    end
end

@testset "community field (Ewald)" begin
    model = setup_abm_community(; n=100, L=800.0, Aphy=3e5, α=1.0, mot="RT",
                                U=40.0, dt=0.1, rng=Xoshiro(3))
    P = model.neighborlist.ypositions
    Rs = model.phytoplankton_radii
    Qs = model.phytoplankton_leakage
    γ = chemoattractant(model).γ
    L = spacesize(model)[1]
    rng = Xoshiro(11)
    # half the probes 0.5–5 μm outside a phytoplankton surface, half anywhere outside
    outside(x) = all(j -> distance(x, P[j], model) > Rs[j], eachindex(P))
    for (k, a) in enumerate(allagents(model))
        if isodd(k)
            j = rand(rng, eachindex(P))
            a.pos = normalize(randn(rng, SVector{3,Float64})) * (Rs[j] + 0.5 + 4.5rand(rng)) + P[j]
            a.pos = mod.(a.pos, L)
        else
            x = SVector{3}(rand(rng, 3) .* L)
            while !outside(x)
                x = SVector{3}(rand(rng, 3) .* L)
            end
            a.pos = x
        end
    end
    ResidenceTimes.community_field!(model)
    # brute-force sum over all sources and periodic images within 18γ
    rmax = 18γ
    n = ceil(Int, rmax / L)
    for a in allagents(model)
        c, g = 0.0, zero(SVector{3,Float64})
        for j in eachindex(P)
            d0 = distancevector(P[j], position(a), model)
            for i1 in -n:n, i2 in -n:n, i3 in -n:n
                dvec = d0 + L * SVector(i1, i2, i3)
                d = norm(dvec)
                d < rmax || continue
                c += ResidenceTimes._concentration(d, Rs[j], Qs[j], γ)
                g += ResidenceTimes._gradient(d, dvec, Rs[j], Qs[j], γ)
            end
        end
        @test concentration(a, model) ≈ chemoattractant(model).Cb + c rtol=1e-3
        @test norm(gradient(a, model) - g) ≤ 1e-2 * norm(g)
    end
    # the convenience evaluator agrees with the stepping path
    a = model[2]
    @test ResidenceTimes.field_comm(position(a), model) ≈ concentration(a, model)
end

@testset "gradients match finite differences" begin
    h = 1e-5
    e(k) = SVector(ntuple(==(k), 3)) .* 1.0
    # community kernel
    C, G = ResidenceTimes._concentration, ResidenceTimes._gradient
    R, Cs, γ = 5.0, 1.0, 150.0
    for x in (SVector(12.0, 3.0, -4.0), SVector(-6.0, 0.5, 1.0), SVector(40.0, -80.0, 25.0))
        fd = SVector(ntuple(k -> (C(norm(x + h*e(k)), R, Cs, γ) - C(norm(x - h*e(k)), R, Cs, γ)) / 2h, 3))
        @test G(norm(x), x, R, Cs, γ) ≈ fd rtol=1e-6
    end
    # short-range Ewald kernel: derivative, and G_near + G_far → G
    γ, σ = 150.0, 40.0
    for r in (1.0, 12.0, 60.0, 120.0)
        fd = (ResidenceTimes.ewald_near(r + h, γ, σ)[1] - ResidenceTimes.ewald_near(r - h, γ, σ)[1]) / 2h
        @test ResidenceTimes.ewald_near(r, γ, σ)[2] ≈ fd rtol=1e-6
    end
    @test ResidenceTimes.ewald_near(3.5σ, γ, σ)[1] < 1e-5 * exp(-3.5σ/γ) / 3.5σ
    # single-source field
    model = setup_abm(; n=1, L=1000.0, R=5.0, Cs=1.0)
    a = model[1]
    origin = chemoattractant(model).origin
    cfield(p) = (a.pos = p; ResidenceTimes.field_exp(a, model))
    for x in (SVector(12.0, 3.0, -4.0), SVector(40.0, -80.0, 25.0))
        p = origin + x
        fd = SVector(ntuple(k -> (cfield(p + h*e(k)) - cfield(p - h*e(k))) / 2h, 3))
        a.pos = p
        @test ResidenceTimes.gradient_exp(a, model) ≈ fd rtol=1e-6
    end
end

@testset "community exposure and checkpoints" begin
    kw = (; n=200, L=800.0, Aphy=1e6, mot="RRF", U=40.0, dt=0.1)

    @testset "exposure accumulates only after equilibration" begin
        model = setup_abm_community(; kw..., teq=0.3, rng=Xoshiro(1))
        @test model.exposure_from == 3
        step!(model, 3)
        @test all(iszero, model.exposure)
        expected = zeros(kw.n)
        for _ in 1:4
            step!(model, 1)
            expected .+= kw.dt .* model.neighborlist.measurements.c
        end
        @test model.exposure == expected
        @test all(>(0), model.exposure)
    end

    @testset "restart from checkpoint is bitwise identical" begin
        model = setup_abm_community(; kw..., teq=0.2, rng=Xoshiro(1))
        step!(model, 3)
        file = save_checkpoint(joinpath(mktempdir(), "ckpt", "state.jls"), model)
        step!(model, 3)
        restarted = load_checkpoint(file)
        @test !isidentical(model, restarted)
        step!(restarted, 3)
        @test isidentical(model, restarted)
        # a 1-ulp change in a single coordinate is detected
        a = restarted[1]
        a.pos = SVector(nextfloat(a.pos[1]), a.pos[2], a.pos[3])
        @test !isidentical(model, restarted)
    end
end
