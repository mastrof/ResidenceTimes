# ABM data fetch/process/archive pipeline

## Problem

Raw ABM simulation output (`data/abm/*.csv`) is produced on the Euler cluster
(`/cluster/scratch/rfoffi/ResidenceTimes/abm/`), which is short-term scratch storage
(periodically cleared). It needs to be:

1. fetched to the local machine (`data/abm/`, limited to ~8GB at a time — disk is
   nearly full),
2. processed by `scripts/process_simulations.jl` into `data/exposure/` and `data/rdf/`,
3. archived long-term to `/media/Elements/ResidenceTimes/data/abm/`.

RDF computation pairs each `Cs>0` file with its matching `Cs=0` control (same
parameters otherwise), so files must be fetched and processed in **groups** —
all `Cs` variants of one parameter configuration together.

Currently 428 files (64GB) sit on the cluster; 255 (39GB) are already staged
locally (over budget); the archive is empty.

## State model

No separate status field. State is derived from three listings:

- **remote**: `ls` on the cluster scratch dir (ground truth for "exists upstream",
  but transient — may be cleared after archiving)
- **local**: `data/abm/` (staged, not yet processed)
- **archived**: `data/abm_archived.txt`, a plain-text append-only log of filenames,
  one per line, written by `process_simulations.jl` after each file is moved to
  `/media/Elements/...`

**Pending** (still need fetching) = `remote − archived − local`.

Using `archived` rather than re-listing `/media/Elements` avoids slow listings on
the external drive, and is robust to the cluster scratch being cleared after the
fact (point 1 above) — once a file is in `abm_archived.txt`, it's considered done
regardless of whether it still exists on remote.

## Components

### 1. `data/abm_archived.txt` (new, git-tracked)

Plain text, one filename per line, append-only. `/data` is gitignored, so add a
`!/data/abm_archived.txt` exception to `.gitignore`.

### 2. `scripts/process_simulations.jl` (small edit)

After the existing `mv` of a raw file to `/media/Elements/...`, append the
filename to `data/abm_archived.txt` (create if missing). No other behavioral
change.

### 3. `scripts/fetch_batch.sh` (new)

A bash script using `ssh`/`sshpass`/`rsync`/standard text tools (no Julia needed —
this is pure file/text manipulation):

1. `sshpass -f ~/Documents/clusterpw ssh ... ls /cluster/scratch/rfoffi/ResidenceTimes/abm/`
   → remote filenames.
2. Read `data/abm/` listing and `data/abm_archived.txt`.
3. `pending = remote − archived − local` (via `comm`/`grep -vFf`).
4. Group `pending` by stripping the `_Cs=...` token from each filename
   (`sed -E 's/_Cs=[^_]*//'`), so each group = one physical configuration with
   all its `Cs` variants.
5. Check current size of `data/abm/` (`du -sh`); if already ≥ 8GB, print a
   message to run `process_simulations.jl` first and exit.
6. Iterate pending groups in listing order; for each group, get remote file
   sizes (`ssh ... du -ch <files>` or `ls -la`), and if adding the group keeps
   `data/abm/` under 8GB, `rsync`/`scp` (via sshpass) all files of that group to
   `data/abm/`. Stop once the budget would be exceeded.
7. Print a summary: groups/files fetched, GB fetched, groups/GB still pending.

## Workflow

```sh
./scripts/fetch_batch.sh        # pulls up to 8GB of new groups
julia --project=. scripts/process_simulations.jl  # processes data/abm/, archives, logs
# repeat until fetch_batch.sh reports nothing pending
```

The existing 255-file/39GB local backlog is over budget, so the first run of
`fetch_batch.sh` will simply report "over budget, run process_simulations.jl
first" — the first real step is processing that backlog, which archives it and
frees space for new fetches.

## Out of scope

- Cs naming inconsistencies in already-staged files (resolved by user separately).
- Automatic cleanup of cluster scratch (user's responsibility / cluster policy).
- Resumability of partial `rsync` transfers beyond rsync's own behavior.
