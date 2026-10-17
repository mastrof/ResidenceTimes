using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using CSV
using DataFrames
using LinearAlgebra
using StatsBase
using StaticArrays

gyration_tensor(g::SubDataFrame) = gyration_tensor(g.x, g.y, g.z)
"""
Gyration tensor S_ab = <(r_a - <r_a>)(r_b - <r_b>)> over the sampled points.
Returns a 3x3 SMatrix. Assumes uniform sampling in time.
"""
function gyration_tensor(x::AbstractVector, y::AbstractVector, z::AbstractVector)
    N = length(x)
    @assert length(y) == N == length(z)
    T = float(eltype(x))
    x̄ = mean(x)
    ȳ = mean(y)
    z̄ = mean(z)
    Sxx = Sxy = Sxz = Syy = Syz = Szz = zero(T)
    @inbounds @simd for i in 1:N
        dx = x[i] - x̄
        dy = y[i] - ȳ
        dz = z[i] - z̄
        Sxx += dx*dx
        Sxy += dx*dy
        Sxz += dx*dz
        Syy += dy*dy
        Syz += dy*dz
        Szz += dz*dz
    end
    SMatrix{3,3,T}(Sxx, Sxy, Sxz,
                   Sxy, Syy, Syz,
                   Sxz, Syz, Szz) ./ N
end
"""
λ sorted descending; Rg2 = tr(S); b = asphericity; c = acylindricity;
κ² = relative shape anisotropy ∈ [0,1] (0 = spherically symmetric, 1 = collinear).
V = rotation-invariant volume in μL
"""
function shape_descriptors(S::SMatrix{3,3})
    λ = reverse(eigvals(Symmetric(Array(S))))   # eigvals(Symmetric) is ascending
    λ1, λ2, λ3 = λ
    Rg2 = λ1 + λ2 + λ3
    b = λ1 - (λ2 + λ3)/2
    c = λ2 - λ3
    κ2 = (b^2 + 0.75c^2) / Rg2^2
    V = det(S) * 1e-9 # μL
    (; λ1, λ2, λ3, Rg2, Rg=sqrt(Rg2), b, c, κ2, V)
end

filenames = readdir(datadir("comm_7classes_tracks"))

for fname in filenames
    df = CSV.read(
        datadir("comm_7classes_tracks", fname),
        DataFrame;
        select=[:id, :x, :y, :z]
    )
    gdf = groupby(df, :id)
    ids, es = eachcol(CSV.read(datadir("comm_7classes_ids", replace(fname, "comm"=>"ids1-99")), DataFrame))
    S = [gyration_tensor(g) for g in gdf]
    sds = vcat([
        DataFrame(; id=ids[i], shape_descriptors(S[i])...)
        for i in eachindex(ids)
    ]...)
    CSV.write(
        datadir("comm_7classes_proc", replace(fname, "comm"=>"gyration")),
        sds
    )
end
