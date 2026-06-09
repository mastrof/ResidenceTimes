using ResidenceTimes
using MicrobeAgents
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
end
