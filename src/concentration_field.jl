export ExpField, field_exp, gradient_exp

@kwdef struct ExpField <: AbstractChemoattractant{3,Float64}
    origin::SVector{3,Float64}
    radius::Float64
    Cs::Float64
    Cb::Float64
    λ::Float64
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
    λ::Real,
)
    @assert length(origin) == 3
    @assert radius >= 0
    @assert Cs >= 0
    @assert Cb >= 0
    @assert λ >= 0
    T = Float64
    ExpField(;
        origin=SVector{3,T}(origin),
        radius=T(radius),
        Cs=T(Cs),
        Cb=T(Cb),
        λ=T(λ)
    )
end

function field_exp(pos, model)
    chemo = chemoattractant(model)
    Cb = chemo.Cb
    Cs = chemo.Cs
    R = chemo.radius
    P = chemo.origin
    λ = chemo.λ
    r = distance(pos, P, model)
    return Cb + Cs*R*exp(-(r-R)/λ) / r
end

function gradient_exp(pos, model)
    chemo = chemoattractant(model)
    Cs = chemo.Cs
    R = chemo.radius
    P = chemo.origin
    λ = chemo.gl
    rvec = distancevector(P, pos, model)
    r2 = dot(rvec, rvec)
    r = sqrt(r2)
    r3 = r * r2
    return SVector{3}(
        -(λ+x)*Cs*exp(-(x-R)/λ) / (λ*r3) * x
        for x in rvec
    )
end
