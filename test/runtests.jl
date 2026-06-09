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
end
