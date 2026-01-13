export ExpField, field_exp, gradient_exp
export CommField, field_comm, gradient_comm

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

@kwdef struct CommField <: AbstractChemoattractant{3,Float64}
    radii::Vector{Float64}
    positions::Vector{SVector{3,Float64}}
    intensities::Vector{Float64}
    Cb::Float64
    γ::Float64
    concentration_field::Function = field_comm
    concentration_gradient::Function = gradient_comm
    concentration_ramp::Function = (pos, model) -> zero(Float64)
    diffusivity::Float64 = 500.0 # μm²/s
end
function CommField(community, Cb, γ)
    spatialcfg, intensities = community
    radii = getproperty.(spatialcfg, :radius)
    positions = SVector{3}.(getproperty.(spatialcfg, :pos))
    CommField(; radii, positions, intensities, Cb, γ)
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
    r = sqrt(r2)
    r3 = r * r2
    return SVector{3}(
        r >= R ? -(γ+r)*Cs*exp(-(r-R)/γ) / (γ*r3) * x : 0.0
        for x in rvec
    )
end

function field_comm(pos, model)
    chemo = chemoattractant(model)
    γ = chemo.γ
    Cb = chemo.Cb
    radii = chemo.radii
    positions = chemo.positions
    intensities = chemo.intensities
    sum_field = sum(zip(radii, positions, intensities)) do (R, P, Cs)
        r = max(distance(pos, P, model), R)
        Cs*R*exp(-(r-R)/γ) / r
    end
    return Cb + sum_field
end

function gradient_comm(pos, model)
    chemo = chemoattractant(model)
    γ = chemo.γ
    radii = chemo.radii
    positions = chemo.positions
    intensities = chemo.intensities
    return sum(zip(radii, positions, intensities)) do (R, P, Cs)
        rvec = distancevector(P, pos, model)
        r2 = dot(rvec, rvec)
        r = sqrt(r2)
        r3 = r*r2
        SVector{3}(
            r >= R ? -(γ+r)*Cs*exp(-(r-R)/γ) / (γ*r3) * x : 0.0
            for x in rvec
        )
    end
end
