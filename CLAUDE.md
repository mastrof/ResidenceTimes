# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Julia research code that measures **residence times / chemoattractant exposure of chemotactic
bacteria near phytoplankton cells** via 3D agent-based simulation. Built on
[MicrobeAgents.jl](https://github.com/mastrof/MicrobeAgents.jl) (which wraps Agents.jl) and organized
as a [DrWatson](https://juliadynamics.github.io/DrWatson.jl/) project.

All physical quantities use μm (length), μM (concentration), s (time).

## Running Julia code

ALWAYS run Julia through the `julia-mcp` `julia_eval` tool (a persistent REPL session with state
preserved between calls), passing this repo as `env_path`. NEVER invoke `julia` via the command line /
Bash. The first time you operate in a session, run `using Revise` explicitly (do not assume it is
already loaded) — afterwards edits to source files are picked up without restarting the session.

## Running things

Always activate the project first; scripts do this with `@quickactivate "ResidenceTimes"`.

```julia
# REPL
julia --project=.
] instantiate          # first-time dependency install

# Run a parameter sweep (uses Distributed/pmap — launch with workers)
julia -p auto --project=. scripts/abm.jl            # single phytoplankton source
julia -p auto --project=. scripts/abm_community.jl  # full phytoplankton community

# Post-process raw simulation CSVs into exposure distributions + RDFs
julia --project=. scripts/process_simulations.jl
julia --project=. scripts/process_simulations_comm.jl

# Figures
julia --project=. scripts/plot_exposure.jl   # etc.
```

### Tests

`test/runtests.jl` covers surface collisions (both models), the Ewald community field (vs. a brute-force periodic sum) and kernel gradients. `Pkg.test()` does **not** work
for this DrWatson project (it fails on the package name/UUID check), so run the suite by including it
in a project-activated session instead:

```julia
include("test/runtests.jl")
```

### Environments (Julia 1.13, workspace, pinned)

The project requires Julia 1.13. `Project.toml` declares a workspace with one sub-project,
`analysis/`, holding the plotting-only deps (GLMakie, CairoMakie, Colors, LaTeXStrings); both share the
root `Manifest.toml`. Production (cluster) uses the root project only, so Makie is never loaded there.
Plot scripts run from the root project and do `push!(LOAD_PATH, projectdir("analysis"))` after
`@quickactivate`. Add plotting deps with `Pkg.activate("analysis"); Pkg.add(...)`.
Every registry package in the manifest is pinned (`Pkg.pin`); `Pkg.pin(all_pkgs=true)` fails on 1.13
because of the URL-tracked MicrobeAgents, so pin by name/uuid instead (MicrobeAgents is fixed by `[sources]`).
Checkpoints from Julia 1.11 runs are not readable on 1.13.

### Cluster vs. local

`src/utils.jl` overrides `DrWatson.datadir`: when the `SCRATCH` env var is set (cluster), output goes
to `$SCRATCH/ResidenceTimes/...`; otherwise to `./data/`. Scripts detect the cluster via
`haskey(ENV, "SCRATCH")`. `process_simulations.jl` additionally moves raw `data/abm` files to an
external drive (`/media/Elements/...`) after processing — that path is hardcoded.

## Data pipeline

1. **Simulate** (`scripts/abm*.jl`): define a parameter grid, expand with DrWatson's `dict_list`,
   run each config in parallel with `pmap`, save one CSV per config named by `savename`. Raw output
   is per-agent time series (position/radial distance `r` and normalized concentration `c`), collected
   only after an equilibration period and then sampled periodically (the `when` callback in `run!`).
2. **Process** (`scripts/process_simulations*.jl`): re-load raw CSVs, compute per-individual
   **exposure** (`combine(groupby(df,:id), :c => sum)`) and **radial distribution functions** (KDE of
   `r` for a given `Cs` divided by the matching `Cs=0` control). Writes to `data/exposure`, `data/rdf`.
   The **community** pipeline is different, since full per-step output is too large:
   - `abm_community.jl` (coarse pass) stores only each bacterium's exposure, accumulated at every
     step after equilibration (`teq` keyword, model property `:exposure`, in μM·s), to
     `comm_7classes_exposure/`, plus full-model checkpoints (`save_checkpoint`, plain
     `Serialization`) at start / end of equilibration / end of run in `comm_7classes_states/`.
   - `community_best-worst-performers.jl` picks the 1000 lowest, 1000 highest and 1000 random ids
     → `comm_7classes_selection/`.
   - `abm_community_track-tails.jl` (fine pass) reloads the equilibration checkpoint, replays the
     run recording the selected ids at every step (Arrow, `comm_7classes_replay/`), and checks
     bitwise (`isidentical`) that it ends on the stored final state (Arrow metadata
     `replay_identical`). Checkpoints are only readable with the same Julia/package versions:
     don't update the Manifest between coarse and fine pass.
   - Older analysis scripts still read the previous layout (`comm_7classes`, `_ids`, `_tracks`).
3. **Plot** (`scripts/plot_*.jl`): load processed CSVs with `collect_datasets` and render with Makie.

`src/reader.jl`'s `collect_datasets(folder; rinclude, rexclude, key=val...)` is the main loader: it
reads all CSVs in a folder, parses parameter values back out of the filenames (via DrWatson
`parse_savename`) into dataframe columns, and supports regex/`key=val` filtering.

`src/wsave.jl` makes `wsave` write `DataFrame`s as CSV.

## Model architecture

Built on **MicrobeAgents v1** (pinned in `Project.toml`/`Manifest.toml` to a commit of
`main` via the GitHub URL, since v1 is unreleased). Two model constructors in `src/model.jl`, both
build a `StandardABM` of the concrete type `microbe_type(mot) = Microbe{3,N,Tuple{Brumley}}`
(chemotaxis is the `Brumley` *behavior*, passed per agent via `behaviors=`) and accept the same
physical/motility parameters (see the keyword defaults there). Microbes get `radius=0.5` explicitly
(v1 defaults to 0, which would also disable Brumley's sensing noise check).

- **`setup_abm`** — single spherical source at the domain center. Concentration comes from `ExpField`
  (exponentially-decaying analytic field, `src/concentration_field.jl`). Stepping is the plain
  MicrobeAgents `microbe_step!`; the source is an `ObstacleSphere` in the `:obstacles` property.
- **`setup_abm_community`** — a polydisperse community of phytoplankton generated by
  `generate_community` (size classes via allometric scaling, packed with `BubbleBath`, per-cell
  leakage from `leaked_concentration`). The field is the **full periodic superposition** of all
  cells' contributions (all sources, all periodic images), computed by Ewald summation
  (`src/ewald.jl`, periodic domains only): the screened kernel `e^{-r/γ}/r` is split with a Gaussian
  of width `σ` (keyword, default 40 μm) into a short-range part (`ewald_near`, summed over a
  `CellListMap.ParticleSystem` neighbor list `:neighborlist` with cutoff `3.5σ`) and a smooth
  long-range part (`EwaldFar`, property `:ewald_far`), built once at setup in Fourier space (type-1
  NUFFT of the sources via NonuniformFFTs.jl, inverse FFT onto a mesh of spacing `σ/2.4`) and
  evaluated by periodic cubic B-spline interpolation (Interpolations.jl). Accuracy vs. the brute-force
  periodic sum: ~1e-4 relative on C, ~3e-3 on ∇C. This requires a custom
  `model_step! = community_step!` that, in order: moves every agent (`move_step!`, obstacle-aware,
  + rotational diffusion), recomputes concentration+gradient on the new positions
  (`community_field!`: `comm_out!` / `OutCommField` near part, then the far part), then runs
  `affect_step!` / `reorient_step!` per microbe.

Chemoattractant fields are `AbstractChemicalField{3}` structs (`ExpField`, `CommField`) carrying
`concentration_field` / `concentration_gradient` / `concentration_ramp` / `diffusivity` function
fields, which MicrobeAgents' generic accessors read. `CommField`'s field/gradient just read
precomputed values out of the neighbor-list `measurements`. Inside a phytoplankton (`r < R`, only
reachable by microbes initialized there) the community field is not physical: the near part is
clamped at `r = R`, the far part is not.

Ewald references: Yukawa-Ewald split, Salin & Caillol, J. Chem. Phys. 113, 10459 (2000); mesh-based
(particle-mesh) Ewald, Darden et al., J. Chem. Phys. 98, 10089 (1993) and Essmann et al.,
J. Chem. Phys. 103, 8577 (1995).

### Collision handling (`src/encounter.jl`)

Surface collisions are **always on** in both models and are detected by MicrobeAgents' `Obstacles`
interface: phytoplankton are an `Obstacles` set of `ObstacleSphere`s (effective contact radius is
`sphere_radius + radius(microbe)`). `src/encounter.jl` only defines the contact rule,
`HaltUntilTumble`: on contact the rest of the step is discarded and the speed is zeroed (heading
preserved); the cell resumes at its next motile-state switch, since `reorient_step!` resamples the
speed. This differs from MicrobeAgents' `Halt`, which keeps the speed and lets a cell slide off the
surface as soon as rotational diffusion turns it away.

## Conventions

- Filenames and parameter round-tripping go entirely through DrWatson `savename` / `parse_savename` /
  `dict_list` — keep parameter dict keys consistent so processing/plotting can recover them.
- Data collection in `run!` uses small named functions (`r`, `c`, `x`, `y`, `z`) as `adata`. Local
  closures don't carry stable column names, so scripts rename the trailing columns by position after
  `run!` (note the `# HACK`/`# rename because...` comments).
- Plotting scripts use `GLMakie` interactively; the `pub/` and figure scripts use `CairoMakie`.
