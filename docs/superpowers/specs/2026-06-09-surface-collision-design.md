# Surface-collision handling for bacteria–phytoplankton encounters

**Date:** 2026-06-09
**Branch:** `stop-bacteria-at-surface`
**Status:** Approved design — ready for implementation planning

## Problem

On `main`, phytoplankton are "ghost particles": chemotactic bacteria feel their
chemoattractant field but swim straight through their bodies. We want bacteria to
**stop at the phytoplankton surface** on contact, accumulate exposure there, and
resume swimming when they reorient away. This must work for **both** model
constructors: the single central source (`setup_abm`) and the polydisperse
community (`setup_abm_community`).

The branch already contains a partial, non-working attempt in `src/encounter.jl`:

- `intersection_ray_sphere` is broken — its body references undefined `p`/`d`, and
  it is declared with three arguments but called with four.
- `collision_move_step!` moves the cell to the surface but never zeroes its speed,
  so the cell is not actually stopped (it penetrates on the next step). It also
  reads `phytoplankton_radii` / `neighborlist`, which exist **only** in the
  community model, and loops over *all* phytoplankton (`# PERF` note).
- `resolve_collision!` (which *does* zero speed) is dead code.
- Neither model is wired to any collision step: `setup_abm` uses the default
  `microbe_step!`; `setup_abm_community` uses `community_step!`, which never calls a
  collision routine.

Separately, the current `community_step!` (`src/model.jl:160`) has latent bugs
independent of collisions: it references an undefined `nlist` and calls the
non-mutating `map_pairwise` instead of `map_pairwise!`. These are fixed as part of
restructuring that function.

## Decisions

| Question | Decision |
| --- | --- |
| Stopping mechanics | **Pin + halt until tumble.** Clamp the cell to the surface and set `speed = 0`. It stays put (accumulating exposure) until its next reorientation event restores speed in a new direction; if that still points inward, it re-pins. |
| API scope | **Replace entirely.** Collisions are always on; ghost mode is removed. No toggle keyword. |
| Dependency on MicrobeAgents collision code | **None.** Reimplement the collision math in this repo (copying the quadratic from MicrobeAgents is fine). MicrobeAgents' `HyperSphere`/`line_sphere_intersection`/`is_encounter` will be phased out in future versions, so collision detection here must be self-contained and stable against those updates. |
| Verification | **Both** deterministic unit tests and an end-to-end verification script. |
| Architecture | **Approach A:** one shared `resolve_collision!` core; model-specific candidate finding (one sphere for single-source; a dedicated small-cutoff neighbor list for the community). |

### Why the pin mechanism is sound

In MicrobeAgents, `velocity(m) = direction(m) .* speed(m)`: the microbe's `vel`
field stores the unit *direction*; `speed` is separate. Therefore `speed = 0`
genuinely halts motion (`move_agent!` walks `velocity·dt = 0`) while preserving the
heading. Rotational diffusion keeps rotating the heading each step; the cell only
regains speed at a reorientation event (`reorient_step!` → `update_speed!`),
resuming in its reoriented direction. For run-reverse / run-reverse-flick motility a
reversal points the cell straight back out of the sphere, so it escapes cleanly.

## Architecture

### Core primitives (`src/encounter.jl`, rewritten)

- **`ray_sphere_fraction(f, d, R) → Float64`** (replaces `intersection_ray_sphere`).
  - `f`: displacement-start position relative to the sphere centre. The caller
    computes it with `distancevector(centre, pos, model)` so periodic boundaries
    are respected.
  - `d`: the step displacement `velocity(microbe) * dt`.
  - `R`: summed radius `R_sphere + radius(microbe)`.
  - Solves `a = d·d`, `b = 2 f·d`, `c = f·f − R²`; returns the smallest root in
    `[0, 1]`, otherwise `1.0` (step does not cross the sphere within this step).
  - Guards `a ≈ 0` (a pinned, zero-speed cell → zero displacement) by returning
    `1.0`. Self-contained; only depends on `LinearAlgebra.dot`.

- **`resolve_collision!(microbe, model, α_hit)`** — the shared core both models call:
  - `move_agent!(microbe, model, α_hit * dt)` — advances along the velocity by the
    clamped fraction (respects the periodic space via `walk!`).
  - if `α_hit < 1`: `microbe.speed = 0` (pin at the surface).

### Single-source model (`setup_abm`)

The only sphere is `chemoattractant(model).origin` with radius
`chemoattractant(model).radius`. Wire `agent_step! = microbe_step_collision!`:

```
microbe_step_collision!(microbe, model):
    collision_move_step!(microbe, model)
    affect_step!(microbe, model)
    reorient_step!(microbe, model)

collision_move_step!(microbe, model):
    dt = abmtimestep(model)
    f  = distancevector(origin, position(microbe), model)
    d  = velocity(microbe) * dt
    α  = ray_sphere_fraction(f, d, R + radius(microbe))
    resolve_collision!(microbe, model, α)
    rotational_diffusion!(microbe, model)
```

No neighbor list — a single sphere. `affect_step!` / `reorient_step!` are the
existing MicrobeAgents subroutines, reused unchanged.

### Community model (`setup_abm_community`)

There are ~hundreds of phytoplankton against ~5×10⁵ swimmers, so per-agent brute
force over all cells is infeasible. Collision detection needs the **pre-move**
positions, whereas the existing field `map_pairwise!` runs on **post-move**
positions — so the two cannot share a single pass without changing field timing.

Add a **dedicated collision neighbor list**: a second `CellListMap.ParticleSystem`
with cutoff `Rmax + microbe_radius + max_speed · dt`. With that small cutoff each
swimmer has ~0 candidates and the pass is nearly free. Its output records, per
microbe, the **earliest hit fraction `αmin`** across nearby spheres.

Restructured `community_step!` (also fixing the existing `nlist` /
`map_pairwise` bugs):

1. Refresh the collision list's `xpositions` to current (pre-move) positions;
   `map_pairwise!` → `αmin[i]` per microbe (min hit fraction over near spheres; the
   pairwise function uses each microbe's current velocity to form `d`).
2. Per microbe: `resolve_collision!(microbe, model, αmin[i])`;
   `rotational_diffusion!(microbe, model)`.
3. Refresh the field neighbor list to the new positions; run the existing field
   `map_pairwise!` (concentration + gradient) — timing unchanged from `main`.
4. Per microbe: `affect_step!(microbe, model)`; `reorient_step!(microbe, model)`.

The collision output is a small struct (in `src/concentration_field.jl`, beside
`OutCommField`) holding `αmin::Vector{Float64}`, with `CellListMap.copy_output`,
`reset_output!` (reset to `1.0`), and `reducer` (element-wise **min**) methods.

## Verification

### `test/runtests.jl` (new; add `Test` to `Project.toml` `[extras]` + `[targets]`)

Seeded, deterministic:

- `ray_sphere_fraction` unit cases: head-on hit returns the expected fraction;
  no-hit returns `1.0`; zero displacement returns `1.0`; grazing/tangent case.
- A cell aimed head-on at a sphere stops at the surface (centre distance ≈ `R + r`
  within float tolerance) and has `speed == 0`.
- A grazing / near-miss cell completes its full step unclamped (`speed` unchanged).
- A pinned cell escapes after a reversal (RunReverse): once reoriented outward it
  moves away and `speed > 0`.
- A hit across the periodic boundary clamps correctly (sphere near a domain edge).

### `scripts/verify_collisions.jl` (new)

Short run of each model with the **no-penetration invariant** asserted at every
sampled step (no swimmer centre inside any sphere beyond tolerance), plus a
trajectory plot of a few cells grazing a sphere for visual confirmation. Matches
the repo's script-driven workflow.

## Files touched

- `src/encounter.jl` — rewrite: `ray_sphere_fraction`, `resolve_collision!`,
  `collision_move_step!`, `microbe_step_collision!`.
- `src/model.jl` — wire `setup_abm` (`agent_step!`) and restructure
  `community_step!`; remove ghost-mode behavior.
- `src/concentration_field.jl` — collision output struct + CellListMap methods.
- `Project.toml` — add `Test` to `[extras]` and a `[targets]` test target.
- `test/runtests.jl` — new.
- `scripts/verify_collisions.jl` — new.

## Out of scope

- Bacterium–bacterium collisions (swimmers remain non-interacting).
- Hydrodynamic interactions or surface accumulation beyond the geometric pin.
- Changing the chemoattractant field model or its timing.
- Re-running / regenerating the published datasets and figures.
