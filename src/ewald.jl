export EwaldFar

#==
Ewald summation of the community field in a periodic box.

Each phytoplankton j is, outside its surface, the point source
    Cs_j R_j e^{-(r-R_j)/γ} / r = Q_j G(r),   G(r) = e^{-r/γ}/r,   Q_j = Cs_j R_j e^{R_j/γ}.
G is split with a Gaussian of width σ into G = G_near + G_far (Yukawa-Ewald split,
Salin & Caillol, J. Chem. Phys. 113, 10459 (2000)):
- G_near decays like erfc(r/σ) and is summed directly over sources within r_c (neighbor list);
- G_far is smooth; its sum over all sources and periodic images is computed once
  in Fourier space, where
      Ĉ_far(k) = (4π/V) e^{-(k² + γ⁻²)σ²/4} / (k² + γ⁻²) Σ_j Q_j e^{-ik·c_j},
  using a type-1 NUFFT (NonuniformFFTs.jl) for the structure factor Σ_j Q_j e^{-ik·c_j}.
  The far field is sampled on a uniform mesh by inverse FFT and evaluated at arbitrary
  points by periodic cubic B-spline interpolation (Interpolations.jl).
==#

"""
    EwaldFar(positions, strengths, L, γ, σ; h=σ/2.4)

Far (smooth) part of the Ewald-split screened field `Σ_j Q_j e^{-r/γ}/r` summed over
all sources and periodic images in a cubic box of side `L`, sampled on a mesh of
spacing ≲ `h` and interpolated with periodic cubic B-splines.
Call `far(x)` to get the far concentration at `x`, `far_gradient(far, x)` for its gradient.
"""
struct EwaldFar{I}
    itp::I
    h::Float64
    L::Float64
end
function EwaldFar(positions, strengths, L, γ, σ; h=σ/2.4)
    M = 2ceil(Int, L / (2h)) # even number of mesh nodes (= Fourier modes) per side
    plan = PlanNUFFT(Float64, (M, M, M); m=HalfSupport(8))
    set_points!(plan, ntuple(k -> [mod(p[k], L) * 2π / L for p in positions], 3))
    S = Array{ComplexF64}(undef, size(plan))
    exec_type1!(S, plan, collect(Float64, strengths))
    kx = rfftfreq(M, M) .* (2π / L)
    ky = fftfreq(M, M) .* (2π / L)
    κ2 = γ^-2
    for I in CartesianIndices(S)
        k2 = kx[I[1]]^2 + ky[I[2]]^2 + ky[I[3]]^2
        S[I] *= 4π / L^3 * exp(-(k2 + κ2) * σ^2 / 4) / (k2 + κ2)
    end
    # drop the Nyquist modes so that the truncated series is real
    S[M÷2+1, :, :] .= 0
    S[:, M÷2+1, :] .= 0
    S[:, :, M÷2+1] .= 0
    mesh = brfft(S, M, 1:3) # far field at nodes x = (i-1)h
    itp = extrapolate(interpolate(mesh, BSpline(Cubic(Periodic(OnCell())))), Periodic())
    EwaldFar(itp, L / M, Float64(L))
end
@inline _mesh_coords(far::EwaldFar, x) = ntuple(k -> mod(x[k], far.L) / far.h + 1, 3)
(far::EwaldFar)(x)::Float64 = far.itp(_mesh_coords(far, x)...)
far_gradient(far::EwaldFar, x)::SVector{3,Float64} =
    SVector{3,Float64}(Interpolations.gradient(far.itp, _mesh_coords(far, x)...)) ./ far.h

"""
    ewald_near(r, γ, σ)

Short-range part of the Ewald-split kernel, `G_near(r) = G(r) - G_far(r)` with
`G(r) = e^{-r/γ}/r`, and its radial derivative:

    G_near(r) = [e^{-r/γ} erfc(r/σ - σ/2γ) + e^{r/γ} erfc(r/σ + σ/2γ)] / 2r
"""
function ewald_near(r, γ, σ)
    a = σ / (2γ)
    u = r / σ
    em, ep = exp(-r / γ), exp(r / γ)
    E1, E2 = erfc(u - a), erfc(u + a)
    F = em * E1 + ep * E2
    dF = (ep * E2 - em * E1) / γ - 4 / (σ * sqrt(π)) * exp(-u^2 - a^2)
    G = F / (2r)
    dG = dF / (2r) - G / r
    return G, dG
end

point_strength(R, Cs, γ) = Cs * R * exp(R / γ)
