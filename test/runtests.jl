using ResidenceTimes
using MicrobeAgents
using Agents: step!
using StaticArrays
using LinearAlgebra
using Random
using CellListMap
using Test

@testset "ResidenceTimes collisions" begin
    @testset "ray_sphere_fraction" begin
        rsf = ResidenceTimes.ray_sphere_fraction
        # head-on: start 3 to the -x of centre, step +4x, R=1 → cross surface at t=0.5
        @test rsf(SVector(-3.0, 0.0, 0.0), SVector(4.0, 0.0, 0.0), 1.0) ≈ 0.5
        # moving away from the sphere → no clamp
        @test rsf(SVector(-3.0, 0.0, 0.0), SVector(-4.0, 0.0, 0.0), 1.0) == 1.0
        # zero displacement (pinned cell) → no clamp
        @test rsf(SVector(-3.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 1.0) == 1.0
        # lateral miss (offset 2 > R) → no clamp
        @test rsf(SVector(-3.0, 2.0, 0.0), SVector(4.0, 0.0, 0.0), 1.0) == 1.0
        # exactly on the surface, pointing OUTWARD → free to leave (entry root only)
        @test rsf(SVector(-1.0, 0.0, 0.0), SVector(-4.0, 0.0, 0.0), 1.0) == 1.0
        # exactly on the surface, pointing INWARD → cannot advance
        @test rsf(SVector(-1.0, 0.0, 0.0), SVector(4.0, 0.0, 0.0), 1.0) == 0.0
    end

    @testset "single-source collision step" begin
        # one swimmer, central sphere R=10 at the domain centre
        model = setup_abm(; n=1, L=1000.0, R=10.0, U=40.0, mot="RT", dt=0.1, Cs=1.0)
        origin = chemoattractant(model).origin

        # head-on: start 12 μm from centre, step 4 μm inward → clamp at the surface (r=10)
        a = model[1]
        a.pos = origin .- SVector(12.0, 0.0, 0.0)
        a.vel = SVector(1.0, 0.0, 0.0)
        a.speed = 40.0
        ResidenceTimes.collision_move_step!(a, model)
        R_eff = chemoattractant(model).radius + radius(a)
        @test isapprox(distance(a, origin, model), R_eff; atol=1e-6)  # pinned at surface
        @test a.speed == 0.0                                          # halted

        # at the nominal sphere radius, pointing outward → not clamped, moves away, keeps speed
        b = model[1]
        b.pos = origin .- SVector(10.0, 0.0, 0.0)   # at the bare sphere radius (inside R_eff)
        b.vel = SVector(-1.0, 0.0, 0.0)             # away from centre
        b.speed = 40.0
        ResidenceTimes.collision_move_step!(b, model)
        @test b.speed == 40.0
        @test distance(b, origin, model) > R_eff
    end

    @testset "collision across periodic boundary" begin
        model = setup_abm(; n=1, L=100.0, R=1.0, U=40.0, dt=0.1)
        a = model[1]
        C = SVector(2.0, 50.0, 50.0)       # sphere centre near the x=0 edge
        a.pos = SVector(99.0, 50.0, 50.0)  # 3 μm from C across the boundary
        a.vel = SVector(1.0, 0.0, 0.0)     # heading +x, wrapping toward C
        a.speed = 40.0                     # step 4 μm
        f = distancevector(C, position(a), model)   # minimal image → (-3, 0, 0)
        d = velocity(a) .* abmtimestep(model)
        R_eff = 1.0 + radius(a)            # sphere radius + microbe radius
        α = ResidenceTimes.ray_sphere_fraction(f, d, R_eff)
        @test α ≈ (3.0 - R_eff) / 4.0     # entry fraction: (|f| - R_eff) / |d|
    end

    @testset "setup_abm uses collision stepping" begin
        # head-on at the central sphere: with the default (ghost) stepping the swimmer
        # would move the full 4 μm and end up inside (distance 8 < R_eff); once the
        # collision agent_step! is wired it is clamped to the surface instead.
        # NOTE: the Brumley microbe radius is 0.5, so R_eff = sphere radius + 0.5 = 10.5.
        model = setup_abm(; n=1, L=1000.0, R=10.0, U=40.0, mot="RR", dt=0.1, Cs=1.0)
        origin = chemoattractant(model).origin
        a = model[1]
        a.pos = origin .- SVector(12.0, 0.0, 0.0)
        a.vel = SVector(1.0, 0.0, 0.0)
        a.speed = 40.0
        step!(model, 1)                       # runs the wired agent_step!
        R_eff = chemoattractant(model).radius + radius(a)
        @test distance(a, origin, model) ≥ R_eff - 1e-6   # never penetrated the sphere
    end

    @testset "OutCollision reducer/reset" begin
        out = ResidenceTimes.OutCollision([0.4, 1.0, 0.8])
        # reducer keeps the element-wise minimum
        other = ResidenceTimes.OutCollision([0.7, 0.5, 0.9])
        r = CellListMap.reducer(out, other)
        @test r.αmin == [0.4, 0.5, 0.8]
        # reset returns every entry to 1.0 (no collision)
        CellListMap.reset_output!(r)
        @test all(r.αmin .== 1.0)
        # copy is independent of its source
        src = ResidenceTimes.OutCollision([0.2, 0.3])
        c = CellListMap.copy_output(src)
        c.αmin[1] = 9.0
        @test c.αmin == [9.0, 0.3]
        @test src.αmin == [0.2, 0.3]   # source untouched
    end

    @testset "community surface collision" begin
        # L must exceed 2*field_cutoff = 2*2.5γ = 750 μm (CellListMap unit-cell rule)
        model = setup_abm_community(; n=1, L=800.0, Aphy=1e6, mot="RRF",
                                    U=40.0, dt=0.1, rng=Xoshiro(1))
        @test length(model.phytoplankton_radii) ≥ 1

        # aim the single swimmer straight at the first phytoplankton, 5 μm outside it,
        # fast enough (step 8 μm) to overshoot the surface without clamping
        P = model.neighborlist.ypositions[1]
        R = model.phytoplankton_radii[1]                         # bare sphere radius
        a = model[1]
        dir = normalize(distancevector(position(a), P, model))   # unit vector toward P
        a.pos = P .- dir .* (R + 5.0)
        a.vel = dir
        a.speed = 80.0
        step!(model, 1)                                          # runs community_step!
        # clamped at the surface instead of penetrating (R_eff = sphere + microbe radius)
        R_eff = R + radius(a)
        @test isapprox(distance(a, P, model), R_eff; atol=1e-4)
    end
end
