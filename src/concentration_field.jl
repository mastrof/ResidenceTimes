export ExpField, field_exp, gradient_exp
export CommField, field_comm, gradient_comm
export comm_out!, OutCommField

@kwdef struct ExpField <: AbstractChemoattractant{3}
    origin::SVector{3,Float64}
    radius::Float64
    Cs::Float64
    Cb::Float64
    γ::Float64
    concentration_field::Function = field_exp
    concentration_gradient::Function = gradient_exp
    concentration_ramp::Function = (_, _) -> zero(Float64)
    diffusivity::Function = (_, _) -> 500.0 # μm²/s
end
function ExpField(
    origin::AbstractVector{<:Real},
    radius::Real,
    Cs::Real,
    Cb::Real,
    γ::Real,
)
    @assert length(origin) == 3
    @assert radius >= 0
    @assert Cs >= 0
    @assert Cb >= 0
    @assert γ >= 0
    T = Float64
    ExpField(;
        origin=SVector{3,T}(origin),
        radius=T(radius),
        Cs=T(Cs),
        Cb=T(Cb),
        γ=T(γ)
    )
end

@kwdef struct CommField <: AbstractChemoattractant{3}
    Cb::Float64
    γ::Float64
    concentration_field::Function = field_comm
    concentration_gradient::Function = gradient_comm
    concentration_ramp::Function = (_, _) -> zero(Float64)
    diffusivity::Function = (_, _) -> 500.0 # μm²/s
end
function CommField(Cb, γ)
    CommField(; Cb, γ)
end

function field_exp(microbe, model)
    chemo = chemoattractant(model)
    Cb = chemo.Cb
    Cs = chemo.Cs
    R = chemo.radius
    P = chemo.origin
    γ = chemo.γ
    pos = position(microbe)
    r = max(distance(pos, P, model), R)
    return Cb + Cs*R*exp(-(r-R)/γ) / r
end

function gradient_exp(microbe, model)
    chemo = chemoattractant(model)
    Cs = chemo.Cs
    R = chemo.radius
    P = chemo.origin
    γ = chemo.γ
    pos = position(microbe)
    rvec = distancevector(P, pos, model)
    r2 = dot(rvec, rvec)
    r = sqrt(r2)
    r3 = r * r2
    return SVector{3}(
        r >= R ? -(γ+r)*Cs*exp(-(r-R)/γ) / (γ*r3) * x : 0.0
        for x in rvec
    )
end

# keep this for convenience
function field_comm(pos::SVector{3,T}, model::ABM) where {T}
    chemo = chemoattractant(model)
    γ = chemo.γ
    Cb = chemo.Cb
    radii = model.phytoplankton_radii
    positions = model.neighborlist.ypositions
    intensities = model.phytoplankton_leakage
    sum_field = sum(zip(radii, positions, intensities)) do (R, P, Cs)
        r = max(distance(pos, P, model), R)
        Cs*R*exp(-(r-R)/γ) / r
    end
    return Cb + sum_field
end

mutable struct OutCommField
    c::Vector{Float64}
    dc::Vector{SVector{3,Float64}}
end
function CellListMap.copy_output(x::OutCommField)
    return OutCommField(copy(x.c), copy(x.dc))
end
function CellListMap.reset_output!(x::OutCommField)
    for i in eachindex(x.c)
        x.c[i] = 0.0
        x.dc[i] = zero(SVector{3,Float64})
    end
    return x
end
function CellListMap.reducer(x::OutCommField, y::OutCommField)
    for i in eachindex(x.c)
        x.c[i] += y.c[i]
        x.dc[i] += y.dc[i]
    end
    return OutCommField(x.c, x.dc)
end

function comm_out!(x,y,i,j,d2,out::OutCommField,model::ABM)
    P = abmproperties(model)[:neighborlist].ypositions[j]
    R = abmproperties(model)[:phytoplankton_radii][j]
    Cs = abmproperties(model)[:phytoplankton_leakage][j]
    chemo = chemoattractant(model)
    γ = chemo.γ
    d = sqrt(d2)
    dvec = distancevector(P, position(model[i]), model)
    out.c[i] += _concentration(d, R, Cs, γ)
    out.dc[i] += _gradient(d, dvec, R, Cs, γ)
    return out
end
function _concentration(d, R, Cs, γ)::Float64
    r = max(d, R)
    # no background here!
    Cs * R/r * exp(-(r-R)/γ)
end
function _gradient(r, rvec, R, Cs, γ)::SVector{3,Float64}
    r3 = r * r * r
    if r >= R
        SVector{3,Float64}(
            -(γ+r)*Cs*exp(-(r-R)/γ) / (γ*r3) * x
            for x in rvec
        )
    else
        zero(SVector{3,Float64})
    end
end

function field_comm(microbe::AbstractMicrobe, model::ABM)::Float64
    Csum = abmproperties(model)[:neighborlist].measurements.c[microbe.id]
    Cb = chemoattractant(model).Cb
    Cb + Csum
end
function gradient_comm(microbe::AbstractMicrobe, model::ABM)::SVector{3,Float64}
    abmproperties(model)[:neighborlist].measurements.dc[microbe.id]
end
