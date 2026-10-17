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
        # L must exceed 2*field_cutoff = 2*2.5γ = 750 μm (CellListMap unit-cell rule)
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

@testset "community field" begin
    model = setup_abm_community(; n=50, L=800.0, Aphy=3e5, α=1.0, mot="RT",
                                U=40.0, dt=0.1, rng=Xoshiro(3))
    step!(model, 3)
    P = model.neighborlist.ypositions
    Rs = model.phytoplankton_radii
    Qs = model.phytoplankton_leakage
    γ = chemoattractant(model).γ
    rc = model.neighborlist.cutoff
    # brute-force sum over the sources within the neighbor-list cutoff (minimum image)
    for a in allagents(model)
        c, g = 0.0, zero(SVector{3,Float64})
        for j in eachindex(P)
            dvec = distancevector(P[j], position(a), model)
            d = norm(dvec)
            d < rc || continue
            c += ResidenceTimes._concentration(d, Rs[j], Qs[j], γ)
            g += ResidenceTimes._gradient(d, dvec, Rs[j], Qs[j], γ)
        end
        @test concentration(a, model) ≈ chemoattractant(model).Cb + c
        @test gradient(a, model) ≈ g
    end
end

