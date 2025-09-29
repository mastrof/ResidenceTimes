export ExpField, field_exp, gradient_exp

@kwdef struct ExpField <: AbstractChemoattractant{3,Float64}
    origin::SVector{3,Float64}
    radius::Float64
    Cs::Float64
    Cb::Float64
    γ::Float64
    concentration_field::Function = field_exp
    concentration_gradient::Function = gradient_exp
    concentration_ramp::Function = (pos, model) -> zero(Float64)
    diffusivity::Float64 = 500.0 # μm²/s
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

function field_exp(pos, model)
    chemo = chemoattractant(model)
    Cb = chemo.Cb
    Cs = chemo.Cs
    R = chemo.radius
    P = chemo.origin
    γ = chemo.γ
    r = max(distance(pos, P, model), R)
    return Cb + Cs*R*exp(-(r-R)/γ) / r
end

function gradient_exp(pos, model)
    chemo = chemoattractant(model)
    Cs = chemo.Cs
    R = chemo.radius
    P = chemo.origin
    γ = chemo.γ
    rvec = distancevector(P, pos, model)
    r2 = dot(rvec, rvec)
    r = max(sqrt(r2), R)
    r3 = r * r * r
    return SVector{3}(
        -(γ+r)*Cs*exp(-(r-R)/γ) / (γ*r3) * x
        for x in rvec
    )
end
