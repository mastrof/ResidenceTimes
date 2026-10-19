export CommunityField, field_gradient, save_community_field, load_community_field

"""
    CommunityField(positions, radii, leakage, L, γ, σ, Cb, far)

Standalone evaluator of the community concentration field used in the simulations
(same Ewald split as `setup_abm_community`): a direct near-field sum over the
phytoplankton (minimum image) plus the interpolated far field `far`.
`leakage` is the surface concentration `Cs` of each cell.
Call `field(x)` for the concentration (μM) and `field_gradient(field, x)` for its gradient.
"""
struct CommunityField{F<:EwaldFar}
    positions::Vector{SVector{3,Float64}}
    radii::Vector{Float64}
    leakage::Vector{Float64}
    L::Float64
    γ::Float64
    σ::Float64
    Cb::Float64
    far::F
end

function CommunityField(positions, radii, leakage, L, γ, σ, Cb; h=σ/2.4)
    pos = SVector{3,Float64}.(positions)
    Q = point_strength.(radii, leakage, γ)
    far = EwaldFar(pos, Q, L, γ, σ; h)
    CommunityField(pos, Float64.(radii), Float64.(leakage), Float64(L), Float64(γ), Float64(σ), Float64(Cb), far)
end

@inline _minimage(d, L) = d - L * round(d / L)

function (f::CommunityField)(x)
    rc = 3.5f.σ
    c = f.Cb + f.far(x)
    for j in eachindex(f.positions)
        d = norm(_minimage.(x .- f.positions[j], f.L))
        d < rc || continue
        R = f.radii[j]
        c += point_strength(R, f.leakage[j], f.γ) * first(ewald_near(max(d, R), f.γ, f.σ))
    end
    return c
end

function field_gradient(f::CommunityField, x)
    rc = 3.5f.σ
    g = far_gradient(f.far, x)
    for j in eachindex(f.positions)
        dvec = SVector{3,Float64}(_minimage.(x .- f.positions[j], f.L))
        d = norm(dvec)
        (d < rc && d >= f.radii[j]) || continue
        Q = point_strength(f.radii[j], f.leakage[j], f.γ)
        g += Q * ewald_near(d, f.γ, f.σ)[2] / d * dvec
    end
    return g
end

"""
    save_community_field(file, field::CommunityField)

Write the field (far-field mesh, sources, and parameters) to a JLD2 file.
"""
function save_community_field(file, f::CommunityField)
    jldsave(file;
        mesh=f.far.mesh,
        positions=[Vector(p) for p in f.positions],
        radii=f.radii, leakage=f.leakage,
        L=f.L, γ=f.γ, σ=f.σ, Cb=f.Cb,
    )
    return file
end

"""
    load_community_field(file) -> CommunityField

Inverse of [`save_community_field`](@ref); the result is callable, `field(x)`.
"""
function load_community_field(file)
    d = load(file)
    CommunityField(
        SVector{3,Float64}.(d["positions"]), d["radii"], d["leakage"],
        d["L"], d["γ"], d["σ"], d["Cb"], EwaldFar(d["mesh"], d["L"]),
    )
end
