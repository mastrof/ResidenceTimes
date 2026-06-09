#==
Self-contained surface-collision handling for bacteria–phytoplankton encounters.
Deliberately independent of MicrobeAgents' own collision utilities
(HyperSphere / line_sphere_intersection / is_encounter), which will be phased out.
==#

"""
    ray_sphere_fraction(f, d, R) -> Float64

Fraction `t ∈ [0,1]` of the displacement `d` at which a point starting at `f`
(its position relative to the sphere centre) first crosses *into* the sphere of
radius `R`. Returns `1.0` if the step does not enter the sphere.

Only the entry intersection (smaller root) is considered: a point on or inside the
surface is free to move outward (returns `1.0`) but blocked from moving inward
(returns `0.0`). An initial overlap therefore resolves once the point's heading
turns outward, rather than by swimming deeper.
"""
function ray_sphere_fraction(f::SVector{D}, d::SVector{D}, R::Real)::Float64 where {D}
    a = dot(d, d)
    iszero(a) && return 1.0            # no displacement (e.g. a pinned cell)
    b = 2 * dot(f, d)
    c = dot(f, f) - R * R
    # On or inside the surface (c ≤ 0): block inward motion (0), allow outward (1).
    # Handled explicitly because the entry-root formula below is numerically
    # unstable for c ≈ 0 — a cell pinned exactly on the surface would otherwise
    # leak through (rounding makes c slightly negative and the entry root slightly
    # negative, which the [0,1] test would read as "no collision").
    c ≤ 0 && return b < 0 ? 0.0 : 1.0
    disc = b * b - 4 * a * c
    disc < 0 && return 1.0             # the ray misses the sphere
    t = (-b - sqrt(disc)) / (2a)       # entry point (smaller root)
    return (0.0 ≤ t ≤ 1.0) ? t : 1.0
end

"""
    resolve_collision!(microbe, model, α_hit)

Shared core for both models. Advance `microbe` along its velocity by the clamped
fraction `α_hit` of one timestep. If a collision occurred (`α_hit < 1`), pin the
cell by zeroing its speed; its heading is preserved, so the cell resumes motion
once it next reorients (`reorient_step!` restores speed on a motile-state switch).
"""
function resolve_collision!(microbe::AbstractMicrobe, model::ABM, α_hit::Real)
    move_agent!(microbe, model, α_hit * abmtimestep(model))
    if α_hit < 1
        microbe.speed = zero(microbe.speed)
    end
    return microbe
end

"""
    collision_move_step!(microbe, model)

Single-source variant of MicrobeAgents' `move_step!`: clamp the step against the
one chemoattractant sphere (`chemoattractant(model).origin`/`.radius`), then apply
rotational diffusion.
"""
function collision_move_step!(microbe::AbstractMicrobe, model::ABM)
    chemo = chemoattractant(model)
    f = distancevector(chemo.origin, position(microbe), model)
    d = velocity(microbe) .* abmtimestep(model)
    R = chemo.radius + radius(microbe)
    α = ray_sphere_fraction(f, d, R)
    resolve_collision!(microbe, model, α)
    rotational_diffusion!(microbe, model)
    return microbe
end

"""
    microbe_step_collision!(microbe, model)

Drop-in replacement for `microbe_step!` that resolves single-source surface
collisions during the translation substep.
"""
function microbe_step_collision!(microbe::AbstractMicrobe, model::ABM)
    collision_move_step!(microbe, model)
    affect_step!(microbe, model)
    reorient_step!(microbe, model)
    return microbe
end
