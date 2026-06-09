# Surface-Collision Handling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make chemotactic bacteria stop at phytoplankton surfaces on contact (clamp to surface + zero speed, resume at next reorientation) in both the single-source (`setup_abm`) and community (`setup_abm_community`) models.

**Architecture:** A single self-contained ray–sphere intersection (`ray_sphere_fraction`) and a shared `resolve_collision!` core. The single-source model wires a collision-aware `agent_step!` that tests its one sphere; the community model computes the earliest hit fraction per swimmer through a dedicated small-cutoff CellListMap neighbor list, then calls the same `resolve_collision!`. No dependency on MicrobeAgents' own collision code (it will be phased out).

**Tech Stack:** Julia 1.11, MicrobeAgents 0.6.2, Agents.jl, CellListMap.jl, StaticArrays, Test (stdlib).

---

## Conventions for every task

- **Run Julia only through the `julia-mcp` `julia_eval` tool** with `env_path = "/home/riccardo/Science/ResidenceTimes"`. Never call `julia` via Bash. The first `julia_eval` in a fresh session must run `using Revise` before anything else, so later source edits are hot-reloaded.
- **Run the test suite** with this `julia_eval` code (absolute path; Revise picks up `src/` edits without a restart):
  ```julia
  include("/home/riccardo/Science/ResidenceTimes/test/runtests.jl")
  ```
- All new collision math lives in this repo. Do **not** call `MicrobeAgents.line_sphere_intersection`, `HyperSphere`, or `is_encounter`.

---

## File Structure

- `src/encounter.jl` — **rewritten.** Geometry + stepping: `ray_sphere_fraction`, `resolve_collision!`, `collision_move_step!` (single-source), `microbe_step_collision!`.
- `src/concentration_field.jl` — **appended.** Community collision output: `OutCollision` struct + CellListMap `copy_output`/`reset_output!`/`reducer`, and the `collision_out!` pairwise kernel.
- `src/model.jl` — **modified.** Wire `setup_abm` `agent_step!`; add the collision neighbor list to `setup_abm_community`; restructure `community_step!` (also fixing its pre-existing `nlist`/`map_pairwise` bugs).
- `Project.toml` — **modified.** Add `Test` to `[extras]` + a `[targets]` test target.
- `test/runtests.jl` — **new.** Grows one `@testset` per task.
- `scripts/verify_collisions.jl` — **new.** End-to-end no-penetration assertion + trajectory plot.

---

## Task 0: Test scaffolding

**Files:**
- Modify: `Project.toml`
- Create: `test/runtests.jl`

- [ ] **Step 1: Add the test target to `Project.toml`**

Append these two sections to the end of `Project.toml` (after the existing `[compat]` block):

```toml
[extras]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[targets]
test = ["Test"]
```

- [ ] **Step 2: Create the test entry point**

Create `test/runtests.jl`:

```julia
using ResidenceTimes
using MicrobeAgents
using StaticArrays
using LinearAlgebra
using Random
using CellListMap
using Test

@testset "ResidenceTimes collisions" begin
    # testsets are added by later tasks
end
```

- [ ] **Step 3: Run it to confirm the harness loads**

Run via `julia_eval`:
```julia
using Revise
include("/home/riccardo/Science/ResidenceTimes/test/runtests.jl")
```
Expected: prints a passing `Test Summary:` for `ResidenceTimes collisions` with 0 tests (no errors).

- [ ] **Step 4: Commit**

```bash
git add Project.toml test/runtests.jl
git commit -m "test: add Test target and runtests scaffolding"
```

---

## Task 1: `ray_sphere_fraction`

The fraction `t ∈ [0,1]` of a step at which a swimmer first crosses *into* a sphere; `1.0` if it does not cross this step. Returns only the **entry** root so that a cell on the surface pointing outward is free to leave.

**Files:**
- Modify: `src/encounter.jl` (rewrite — see below; in this task we replace the whole file with just the geometry function, then later tasks add to it)
- Test: `test/runtests.jl`

- [ ] **Step 1: Write the failing test**

Add to `test/runtests.jl`, inside the outer `@testset`, before its closing `end`:

```julia
@testset "ray_sphere_fraction" begin
    rsf = ResidenceTimes.ray_sphere_fraction
    # head-on: start 3 to the -x of centre, step +4x, R=1 → cross surface at t=0.5
    @test rsf(SVector(-3.0, 0.0, 0.0), SVector(4.0, 0.0, 0.0), 1.0) ≈ 0.5
    # moving away from the sphere → no clamp
    @test rsf(SVector(-3.0, 0.0, 0.0), SVector(-4.0, 0.0, 0.0), 1.0) == 1.0
    # zero displacement (pinned cell) → no clamp
    @test rsf(SVector(-3.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 1.0) == 1.0
    # lateral miss (offset 2 > R) → no clamp
    @test rsf(SVector(-3.0, 2.0, 0.0), SVector(4.0, 0.0, 0.0), 1.0) == 1.0
    # exactly on the surface, pointing OUTWARD → free to leave (entry root only)
    @test rsf(SVector(-1.0, 0.0, 0.0), SVector(-4.0, 0.0, 0.0), 1.0) == 1.0
    # exactly on the surface, pointing INWARD → cannot advance
    @test rsf(SVector(-1.0, 0.0, 0.0), SVector(4.0, 0.0, 0.0), 1.0) == 0.0
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```julia
include("/home/riccardo/Science/ResidenceTimes/test/runtests.jl")
```
Expected: FAIL — `UndefVarError: ray_sphere_fraction` (function not defined yet).

- [ ] **Step 3: Replace `src/encounter.jl` with the geometry function**

Overwrite the entire file `src/encounter.jl` with:

```julia
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
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```julia
include("/home/riccardo/Science/ResidenceTimes/test/runtests.jl")
```
Expected: PASS — `ray_sphere_fraction` testset green.

- [ ] **Step 5: Commit**

```bash
git add src/encounter.jl test/runtests.jl
git commit -m "feat: self-contained ray_sphere_fraction for surface collisions"
```

---

## Task 2: `resolve_collision!` and the single-source collision step

**Files:**
- Modify: `src/encounter.jl` (append)
- Test: `test/runtests.jl`

- [ ] **Step 1: Write the failing test**

Add to `test/runtests.jl` inside the outer `@testset`:

```julia
@testset "single-source collision step" begin
    # one swimmer, central sphere R=10 at the domain centre
    model = setup_abm(; n=1, L=1000.0, R=10.0, U=40.0, mot="RT", dt=0.1, Cs=1.0)
    origin = chemoattractant(model).origin

    # head-on: start 12 μm from centre, step 4 μm inward → clamp at the surface (r=10)
    a = model[1]
    a.pos = origin .- SVector(12.0, 0.0, 0.0)
    a.vel = SVector(1.0, 0.0, 0.0)
    a.speed = 40.0
    ResidenceTimes.collision_move_step!(a, model)
    @test isapprox(distance(a, origin, model), 10.0; atol=1e-6)  # pinned at surface
    @test a.speed == 0.0                                          # halted

    # on the surface pointing outward → not clamped, moves away, keeps speed
    b = model[1]
    b.pos = origin .- SVector(10.0, 0.0, 0.0)   # exactly on the surface
    b.vel = SVector(-1.0, 0.0, 0.0)             # away from centre
    b.speed = 40.0
    ResidenceTimes.collision_move_step!(b, model)
    @test b.speed == 40.0
    @test distance(b, origin, model) > 10.0
end

@testset "collision across periodic boundary" begin
    model = setup_abm(; n=1, L=100.0, R=1.0, U=40.0, dt=0.1)
    a = model[1]
    C = SVector(2.0, 50.0, 50.0)       # sphere centre near the x=0 edge
    a.pos = SVector(99.0, 50.0, 50.0)  # 3 μm from C across the boundary
    a.vel = SVector(1.0, 0.0, 0.0)     # heading +x, wrapping toward C
    a.speed = 40.0                     # step 4 μm
    f = distancevector(C, position(a), model)   # minimal image → (-3, 0, 0)
    d = velocity(a) .* abmtimestep(model)
    α = ResidenceTimes.ray_sphere_fraction(f, d, 1.0 + radius(a))
    @test α ≈ 0.5                      # travels 2 of 4 μm to reach the surface
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```julia
include("/home/riccardo/Science/ResidenceTimes/test/runtests.jl")
```
Expected: FAIL — `UndefVarError: collision_move_step!`.

- [ ] **Step 3: Append the stepping functions to `src/encounter.jl`**

Add to the end of `src/encounter.jl`:

```julia
"""
    resolve_collision!(microbe, model, α_hit)

Shared core for both models. Advance `microbe` along its velocity by the clamped
fraction `α_hit` of one timestep. If a collision occurred (`α_hit < 1`), pin the
cell by zeroing its speed; its heading is preserved, so a later reorientation
(`update_speed!`) resumes motion in a new direction.
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
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```julia
include("/home/riccardo/Science/ResidenceTimes/test/runtests.jl")
```
Expected: PASS — both new asserts green. (`collision_move_step!` also calls `rotational_diffusion!`, which changes only the heading, so the position/speed asserts are deterministic.)

- [ ] **Step 5: Commit**

```bash
git add src/encounter.jl test/runtests.jl
git commit -m "feat: resolve_collision! and single-source collision step"
```

---

## Task 3: Wire collisions into `setup_abm`

**Files:**
- Modify: `src/model.jl`
- Test: `test/runtests.jl`

- [ ] **Step 1: Write the failing test**

Add to `test/runtests.jl` inside the outer `@testset`:

```julia
@testset "setup_abm uses collision stepping" begin
    # head-on at the central sphere: with the default (ghost) stepping the swimmer
    # would move the full 4 μm and end up inside (distance 8 < R=10); once the
    # collision agent_step! is wired it is clamped to the surface instead.
    model = setup_abm(; n=1, L=1000.0, R=10.0, U=40.0, mot="RR", dt=0.1, Cs=1.0)
    origin = chemoattractant(model).origin
    a = model[1]
    a.pos = origin .- SVector(12.0, 0.0, 0.0)
    a.vel = SVector(1.0, 0.0, 0.0)
    a.speed = 40.0
    step!(model, 1)                       # runs the wired agent_step!
    @test distance(a, origin, model) ≥ 10.0 - 1e-6   # never penetrated the sphere
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```julia
include("/home/riccardo/Science/ResidenceTimes/test/runtests.jl")
```
Expected: FAIL — with the default `microbe_step!` the swimmer moves the full step and ends up at distance 8, so `distance ≥ 10 - 1e-6` fails.

- [ ] **Step 3: Add `agent_step!` to the `setup_abm` model constructor**

In `src/model.jl`, find this block inside `setup_abm` (it is the one **without** an `agent_step!`):

```julia
    model = StandardABM(Brumley{3,N}, space, dt;
        properties,
        container=Vector,
    )
```

Replace it with:

```julia
    model = StandardABM(Brumley{3,N}, space, dt;
        properties,
        container=Vector,
        agent_step! = microbe_step_collision!,
    )
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```julia
include("/home/riccardo/Science/ResidenceTimes/test/runtests.jl")
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/model.jl test/runtests.jl
git commit -m "feat: enable surface collisions in setup_abm"
```

---

## Task 4: Community collision output (`OutCollision` + `collision_out!`)

**Files:**
- Modify: `src/concentration_field.jl` (append)
- Test: `test/runtests.jl`

- [ ] **Step 1: Write the failing test**

Add to `test/runtests.jl` inside the outer `@testset`:

```julia
@testset "OutCollision reducer/reset" begin
    out = ResidenceTimes.OutCollision([0.4, 1.0, 0.8])
    # reducer keeps the element-wise minimum
    other = ResidenceTimes.OutCollision([0.7, 0.5, 0.9])
    r = CellListMap.reducer(out, other)
    @test r.αmin == [0.4, 0.5, 0.8]
    # reset returns every entry to 1.0 (no collision)
    CellListMap.reset_output!(r)
    @test all(r.αmin .== 1.0)
    # copy is independent
    c = CellListMap.copy_output(ResidenceTimes.OutCollision([0.2, 0.3]))
    c.αmin[1] = 9.0
    @test c.αmin == [9.0, 0.3]
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```julia
include("/home/riccardo/Science/ResidenceTimes/test/runtests.jl")
```
Expected: FAIL — `UndefVarError: OutCollision`.

- [ ] **Step 3: Append the collision output to `src/concentration_field.jl`**

Add to the end of `src/concentration_field.jl`:

```julia
export OutCollision, collision_out!

# Per-microbe earliest hit fraction across nearby phytoplankton, accumulated over
# a CellListMap pairwise pass. αmin[i] == 1.0 means microbe i hits nothing.
mutable struct OutCollision
    αmin::Vector{Float64}
end
function CellListMap.copy_output(x::OutCollision)
    return OutCollision(copy(x.αmin))
end
function CellListMap.reset_output!(x::OutCollision)
    fill!(x.αmin, 1.0)
    return x
end
function CellListMap.reducer(x::OutCollision, y::OutCollision)
    @inbounds for i in eachindex(x.αmin)
        x.αmin[i] = min(x.αmin[i], y.αmin[i])
    end
    return x
end

# Pairwise kernel: clamp microbe i's step against phytoplankton j and keep the
# earliest (minimum) hit fraction. Uses pre-move positions and the microbe's
# current velocity; periodicity is handled by `distancevector`.
function collision_out!(x, y, i, j, d2, out::OutCollision, model::ABM)
    microbe = model[i]
    P = abmproperties(model)[:collisionlist].ypositions[j]
    R = abmproperties(model)[:phytoplankton_radii][j] + radius(microbe)
    f = distancevector(P, position(microbe), model)
    d = velocity(microbe) .* abmtimestep(model)
    α = ray_sphere_fraction(f, d, R)
    @inbounds out.αmin[i] = min(out.αmin[i], α)
    return out
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```julia
include("/home/riccardo/Science/ResidenceTimes/test/runtests.jl")
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/concentration_field.jl test/runtests.jl
git commit -m "feat: OutCollision accumulator and collision_out! kernel"
```

---

## Task 5: Wire collisions into the community model

Add the dedicated small-cutoff collision neighbor list and restructure `community_step!` to (1) detect collisions on pre-move positions, (2) move with clamping, (3) recompute the field, (4) run chemotaxis + reorientation. This also fixes the existing `community_step!` bugs (undefined `nlist`; non-mutating `map_pairwise`; the `speed(microbe)*dt` double-speed move).

**Files:**
- Modify: `src/model.jl`
- Test: `test/runtests.jl`

- [ ] **Step 1: Write the failing test**

Add to `test/runtests.jl` inside the outer `@testset`:

```julia
@testset "community surface collision" begin
    model = setup_abm_community(; n=1, L=300.0, Aphy=1e6, mot="RRF",
                                U=40.0, dt=0.1, rng=Xoshiro(1))
    @test length(model.phytoplankton_radii) ≥ 1

    # aim the single swimmer straight at the first phytoplankton, 5 μm outside it,
    # fast enough (step 8 μm) to overshoot the surface without clamping
    P = model.neighborlist.ypositions[1]
    R = model.phytoplankton_radii[1]
    a = model[1]
    dir = normalize(distancevector(position(a), P, model))   # unit vector toward P
    a.pos = P .- dir .* (R + 5.0)
    a.vel = dir
    a.speed = 80.0
    step!(model, 1)                                          # runs community_step!
    # clamped at the surface instead of penetrating
    @test isapprox(distance(a, P, model), R; atol=1e-4)
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```julia
include("/home/riccardo/Science/ResidenceTimes/test/runtests.jl")
```
Expected: FAIL — the current `community_step!` does not clamp (and references an undefined `nlist`), so the swimmer penetrates or the step errors.

- [ ] **Step 3: Add the collision neighbor list to `setup_abm_community`**

In `src/model.jl`, find this block in `setup_abm_community`:

```julia
    system = ParticleSystem(
        xpositions=zeros(SVector{3,Float64}, n),
        ypositions=phytoplankton_positions,
        unitcell=spacesize(space),
        cutoff=2.5γ,
        output=OutCommField(zeros(n), zeros(SVector{3,Float64}, n)),
        output_name=:measurements,
        parallel=false
    )
    chemo = CommField(Cb, γ)
    properties = Dict(
        :neighborlist => system,
        :phytoplankton_radii => phytoplankton_radii,
        :phytoplankton_leakage => phytoplankton_leakage,
        :chemoattractant => chemo,
    )
```

Replace it with (adds `collisionsystem` and the `:collisionlist` property; cutoff = max sphere radius + max step, with a small margin, since microbe radius is 0):

```julia
    system = ParticleSystem(
        xpositions=zeros(SVector{3,Float64}, n),
        ypositions=phytoplankton_positions,
        unitcell=spacesize(space),
        cutoff=2.5γ,
        output=OutCommField(zeros(n), zeros(SVector{3,Float64}, n)),
        output_name=:measurements,
        parallel=false
    )
    collision_cutoff = Rmax + U*dt + 2.0  # μm: largest sphere + max step + margin
    collisionsystem = ParticleSystem(
        xpositions=zeros(SVector{3,Float64}, n),
        ypositions=phytoplankton_positions,
        unitcell=spacesize(space),
        cutoff=collision_cutoff,
        output=OutCollision(ones(n)),
        output_name=:measurements,
        parallel=false
    )
    chemo = CommField(Cb, γ)
    properties = Dict(
        :neighborlist => system,
        :collisionlist => collisionsystem,
        :phytoplankton_radii => phytoplankton_radii,
        :phytoplankton_leakage => phytoplankton_leakage,
        :chemoattractant => chemo,
    )
```

- [ ] **Step 4: Replace `community_step!`**

In `src/model.jl`, replace the entire existing `community_step!` function:

```julia
function community_step!(model)
    dt = abmtimestep(model)
    # move all cells to new positions
    for microbe in allagents(model)
        move_agent!(microbe, model, speed(microbe)*dt)
        MicrobeAgents.rotational_diffusion!(microbe, model)
    end
    # update neighbor list
    for i in allids(model)
        nlist.xpositions[i] = position(model[i])
    end
    # compute new values of concentration and gradient
    map_pairwise(
        (x,y,i,j,d2,out) -> comm_out!(x,y,i,j,d2,out,model),
        nlist
    )
    # continue with individual microbe steps as usual
    for microbe in allagents(model)
        model.affect!(microbe, model)
        if MicrobeAgents.can_turn(microbe)
            MicrobeAgents.turn!(microbe, model)
        end
        p = switching_probability(microbe, model)
        if rand(abmrng(model)) < p
            update_motilestate!(microbe, model)
            new_motilestate = motilestate(microbe)
            if variantof(new_motilestate) === TurnState && iszero(duration(new_motilestate))
                MicrobeAgents.turn!(microbe, model)
                update_motilestate!(microbe, model)
            end
            MicrobeAgents.update_speed!(microbe, model)
        end
    end
end
```

with:

```julia
function community_step!(model)
    collisionlist = abmproperties(model)[:collisionlist]
    nlist = abmproperties(model)[:neighborlist]

    # 1. detect collisions on the current (pre-move) positions
    for i in allids(model)
        collisionlist.xpositions[i] = position(model[i])
    end
    map_pairwise!(
        (x, y, i, j, d2, out) -> collision_out!(x, y, i, j, d2, out, model),
        collisionlist
    )
    αmin = collisionlist.measurements.αmin

    # 2. move every swimmer with surface clamping, then rotational diffusion
    for microbe in allagents(model)
        resolve_collision!(microbe, model, αmin[microbe.id])
        rotational_diffusion!(microbe, model)
    end

    # 3. recompute concentration + gradient at the new positions
    for i in allids(model)
        nlist.xpositions[i] = position(model[i])
    end
    map_pairwise!(
        (x, y, i, j, d2, out) -> comm_out!(x, y, i, j, d2, out, model),
        nlist
    )

    # 4. chemotactic update and reorientation
    for microbe in allagents(model)
        affect_step!(microbe, model)
        reorient_step!(microbe, model)
    end
end
```

- [ ] **Step 5: Run the test to verify it passes**

Run:
```julia
include("/home/riccardo/Science/ResidenceTimes/test/runtests.jl")
```
Expected: PASS — the swimmer is clamped to within `1e-4` of the surface. (Revise reloads `src/model.jl`; if the `ParticleSystem` struct change is not picked up, restart the session with `mcp__julia__julia_restart` for this `env_path`, then re-run.)

- [ ] **Step 6: Commit**

```bash
git add src/model.jl test/runtests.jl
git commit -m "feat: surface collisions in community model; fix community_step!"
```

---

## Task 6: End-to-end verification script + authoritative test gate

**Files:**
- Create: `scripts/verify_collisions.jl`
- Test: the packaged `]test` target

- [ ] **Step 1: Write the verification script**

Create `scripts/verify_collisions.jl`:

```julia
#================================================================
Verify that no bacterium ever penetrates a phytoplankton, in both
the single-source and community models, and plot a few trajectories
grazing a sphere. Run via the julia-mcp `julia_eval` tool.
================================================================#

using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using MicrobeAgents
using StaticArrays
using LinearAlgebra
using Random
using DataFrames
using CairoMakie

const TOL = 1e-4  # μm; clamped cells sit at the surface within rounding

"No swimmer centre is inside any sphere (centre distance ≥ R + microbe radius)."
function assert_no_penetration(model, centres, radii)
    for a in allagents(model)
        for j in eachindex(radii)
            d = distance(a, centres[j], model)
            @assert d ≥ radii[j] + radius(a) - TOL "penetration: id=$(a.id) d=$d R=$(radii[j])"
        end
    end
end

# --- single-source model ---
model = setup_abm(; n=2_000, L=1000.0, R=10.0, U=40.0, mot="RRF", dt=0.1, Cs=1.0)
centres = [chemoattractant(model).origin]
radii = [chemoattractant(model).radius]
for _ in 1:500
    step!(model, 1)
    assert_no_penetration(model, centres, radii)
end
println("single-source: no penetration over 500 steps ✓")

# --- community model ---
cmodel = setup_abm_community(; n=5_000, L=400.0, Aphy=1e6, mot="RRF",
                             U=40.0, dt=0.1, rng=Xoshiro(7))
ccentres = cmodel.neighborlist.ypositions
cradii = cmodel.phytoplankton_radii
for _ in 1:500
    step!(cmodel, 1)
    assert_no_penetration(cmodel, ccentres, cradii)
end
println("community: no penetration over 500 steps ✓")

# --- trajectory plot: a few swimmers grazing the central sphere ---
tmodel = setup_abm(; n=40, L=400.0, R=20.0, U=40.0, mot="RRF", dt=0.1, Cs=1.0)
origin = chemoattractant(tmodel).origin
# seed swimmers heading at the sphere from one side
for a in allagents(tmodel)
    a.pos = origin .- SVector(60.0, 0.0, 0.0) .+ SVector(0.0, 20.0*(rand()-0.5), 0.0)
    a.vel = SVector(1.0, 0.0, 0.0)
    a.speed = 40.0
end
tracks = [Point2f[] for _ in 1:nagents(tmodel)]
for _ in 1:300
    step!(tmodel, 1)
    for (k, a) in enumerate(allagents(tmodel))
        push!(tracks[k], Point2f(position(a)[1], position(a)[2]))
    end
end
fig = Figure()
ax = Axis(fig[1,1]; aspect=DataAspect(), xlabel="x (μm)", ylabel="y (μm)",
          title="Swimmers grazing a sphere (R=20 μm)")
arc!(ax, Point2f(origin[1], origin[2]), 20.0, 0, 2π; color=:black)
for t in tracks
    lines!(ax, t; color=(:steelblue, 0.5))
end
save(plotsdir("verify_collisions.png"), fig)
println("trajectory plot written to ", plotsdir("verify_collisions.png"))
```

- [ ] **Step 2: Run the verification script**

Run via `julia_eval`:
```julia
include("/home/riccardo/Science/ResidenceTimes/scripts/verify_collisions.jl")
```
Expected: prints both `no penetration … ✓` lines and the plot path, with no `AssertionError`.

- [ ] **Step 3: Run the authoritative packaged test target**

Run via `julia_eval` (uses the `[targets]` test env in a clean subprocess):
```julia
import Pkg
Pkg.test("ResidenceTimes")
```
Expected: PASS — `Testing ResidenceTimes tests passed`.

- [ ] **Step 4: Commit**

```bash
git add scripts/verify_collisions.jl
git commit -m "test: end-to-end no-penetration verification script"
```

---

## Notes / known edge cases

- **Initial overlaps:** swimmers are placed at random and a few may start *inside* a phytoplankton. `ray_sphere_fraction` returns `1.0` for an inside-start (the entry root is negative), so such cells move freely until they exit, then collide normally. The models discard an equilibration window before sampling, so this does not affect measured exposure. No rejection sampling is added (YAGNI).
- **Cutoff:** the community collision cutoff is `Rmax + U·dt + 2.0` μm. If a future sweep raises `U` or `dt` substantially, this is the value to revisit — it must stay ≥ `Rmax + microbe_radius + max_speed·dt`.
- **Field timing unchanged:** the concentration/gradient pass still runs on post-move positions, exactly as on `main`; only the move itself is now clamped.
```
