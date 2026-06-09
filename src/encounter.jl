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

Only the entry intersection (smaller root) is considered: a point already on the
surface moving outward is free to leave (returns `1.0`), while one moving inward
cannot advance (returns `0.0`).
"""
function ray_sphere_fraction(f::SVector{D}, d::SVector{D}, R::Real)::Float64 where {D}
    a = dot(d, d)
    iszero(a) && return 1.0            # no displacement (e.g. a pinned cell)
    b = 2 * dot(f, d)
    c = dot(f, f) - R * R
    disc = b * b - 4 * a * c
    disc < 0 && return 1.0             # the ray misses the sphere
    t = (-b - sqrt(disc)) / (2a)       # entry point (smaller root)
    return (0.0 ≤ t ≤ 1.0) ? t : 1.0
end
