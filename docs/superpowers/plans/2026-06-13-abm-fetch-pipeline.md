# ABM Fetch/Process/Archive Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a fetch/process/archive loop for ABM simulation data, so raw CSVs can be pulled from the Euler cluster in disk-budget-sized batches, processed into exposure/RDF outputs, and archived to external storage, with an append-only log tracking what's already been archived.

**Architecture:** A new `scripts/fetch_batch.sh` computes the "pending" file set as `remote − archived − local` (grouped by parameter configuration so `Cs` variants travel together for RDF pairing), and `rsync`s groups to `data/abm/` until an 8GB local budget would be exceeded. `scripts/process_simulations.jl` gets a one-function change so that every file it archives to `/media/Elements/...` is also appended to `data/abm_archived.txt`, which is git-tracked despite `/data` being gitignored.

**Tech Stack:** Bash (ssh/sshpass/rsync/comm/sed for fetch_batch.sh), Julia (process_simulations.jl edit), bash test script for the pure-text grouping/diff logic.

---

### Task 1: Track `data/abm_archived.txt` in git

**Files:**
- Create: `data/abm_archived.txt`

- [ ] **Step 1: Create the empty archive log**

```bash
touch data/abm_archived.txt
```

- [ ] **Step 2: Verify it's not ignored**

Run: `git status --porcelain data/abm_archived.txt`
Expected: `?? data/abm_archived.txt` (i.e. it shows up as untracked, not silently ignored — the `.gitignore` exception from the previous commit already covers this)

- [ ] **Step 3: Commit**

```bash
git add data/abm_archived.txt
git commit -m "Add empty abm_archived.txt log"
```

---

### Task 2: Log archived filenames from `process_simulations.jl`

**Files:**
- Modify: `scripts/process_simulations.jl`

The script currently has two separate `mv` blocks (lines 49-53 and 56-62) that
move files to `/media/Elements/ResidenceTimes/data/abm/`. Replace both with a
shared `archive!` function that moves the file *and* appends its name to
`data/abm_archived.txt`.

- [ ] **Step 1: Add the `archive!` helper and use it in both mv sites**

Current content (lines 18-63):

```julia
filenames = readdir(datadir("abm"))
L = 1e3 # L is always 1mm
npoints = 100 # points for the rdf sampling
r = range(1, L/2; length=npoints)
for filename in filenames
    prefix, config, suffix = parse_savename(filename)
    df = CSV.read(datadir("abm", filename), DataFrame)
    # distribution of individual exposures
    df_exposure = combine(groupby(df, :id), :c => sum)
    CSV.write(
        datadir("exposure", savename("exposure", config, "csv")),
        df_exposure
    )
    # radial distribution functions
    # each Cs value is compared to corresponding Cs=0 simulation
    iszero(config["Cs"]) && continue
    config_random = copy(config)
    config_random["Cs"] = 0.0
    filename_random = savename(prefix, config_random, suffix)
    !isfile(datadir("abm", filename_random)) && continue
    df_random = CSV.read(datadir("abm", filename_random), DataFrame)
    Pr = kde(df.r; bandwidth=25)
    Pr0 = kde(df_random.r; bandwidth=25)
    k = pdf(Pr, r)
    k0 = pdf(Pr0, r)
    g = k ./ k0
    df_rdf = DataFrame(; r, g)
    CSV.write(
        datadir("rdf", savename("rdf", config, "csv")),
        df_rdf
    )
    # move abm data to hard drive
    mv(
        datadir("abm", filename),
        joinpath("/media/Elements/ResidenceTimes/data/abm/", filename)
    )
end
# the Cs=0 files have been kept, remove them now
for filename in filenames
    !isfile(datadir("abm", filename)) && continue
    mv(
        datadir("abm", filename),
        joinpath("/media/Elements/ResidenceTimes/data/abm/", filename)
    )
end
```

Replace with:

```julia
function archive!(filename)
    mv(
        datadir("abm", filename),
        joinpath("/media/Elements/ResidenceTimes/data/abm/", filename)
    )
    open(datadir("abm_archived.txt"), "a") do io
        println(io, filename)
    end
end

filenames = readdir(datadir("abm"))
L = 1e3 # L is always 1mm
npoints = 100 # points for the rdf sampling
r = range(1, L/2; length=npoints)
for filename in filenames
    prefix, config, suffix = parse_savename(filename)
    df = CSV.read(datadir("abm", filename), DataFrame)
    # distribution of individual exposures
    df_exposure = combine(groupby(df, :id), :c => sum)
    CSV.write(
        datadir("exposure", savename("exposure", config, "csv")),
        df_exposure
    )
    # radial distribution functions
    # each Cs value is compared to corresponding Cs=0 simulation
    iszero(config["Cs"]) && continue
    config_random = copy(config)
    config_random["Cs"] = 0.0
    filename_random = savename(prefix, config_random, suffix)
    !isfile(datadir("abm", filename_random)) && continue
    df_random = CSV.read(datadir("abm", filename_random), DataFrame)
    Pr = kde(df.r; bandwidth=25)
    Pr0 = kde(df_random.r; bandwidth=25)
    k = pdf(Pr, r)
    k0 = pdf(Pr0, r)
    g = k ./ k0
    df_rdf = DataFrame(; r, g)
    CSV.write(
        datadir("rdf", savename("rdf", config, "csv")),
        df_rdf
    )
    # move abm data to hard drive
    archive!(filename)
end
# the Cs=0 files have been kept, remove them now
for filename in filenames
    !isfile(datadir("abm", filename)) && continue
    archive!(filename)
end
```

- [ ] **Step 2: Sanity-check the file parses**

In the persistent Julia session (per CLAUDE.md, `julia_eval` with `env_path` set
to this repo; run `using Revise` first if not already done this session):

```julia
using DrWatson
@quickactivate "ResidenceTimes"
include("scripts/process_simulations.jl")
```

Don't run this yet against the real backlog — this step is just to confirm the
file has no syntax errors. If you want to avoid accidentally processing the
255-file backlog while just checking syntax, instead run:

```julia
using JuliaSyntax
JuliaSyntax.parseall(Expr, read("scripts/process_simulations.jl", String))
```

Expected: returns an `Expr` with no error.

- [ ] **Step 3: Commit**

```bash
git add scripts/process_simulations.jl
git commit -m "Log archived filenames to data/abm_archived.txt"
```

---

### Task 3: `fetch_batch.sh` — pure grouping/diff logic + tests

**Files:**
- Create: `scripts/fetch_batch.sh`
- Create: `scripts/test_fetch_batch.sh`

This task implements the two pure-text-processing functions that don't need
network access, and tests them against fixture data. The ssh/rsync fetch logic
that uses these functions is added in Task 4.

- [ ] **Step 1: Write `scripts/fetch_batch.sh` with the two helper functions, guarded so `main` only runs when executed directly**

```bash
#!/usr/bin/env bash
set -euo pipefail

# group_key: strip the _Cs=<value> token from a filename, so all Cs variants
# of one configuration map to the same key.
group_key() {
    echo "$1" | sed -E 's/_Cs=[^_]*//'
}

# compute_pending: files that exist on remote but are neither already staged
# locally nor already archived.
#   $1 = path to file listing remote filenames (one per line)
#   $2 = local directory (e.g. data/abm)
#   $3 = path to archive log (e.g. data/abm_archived.txt), may not exist
compute_pending() {
    local remote_list="$1" local_dir="$2" archive_log="$3"
    local done_list
    done_list=$(mktemp)
    {
        ls "$local_dir" 2>/dev/null || true
        [[ -f "$archive_log" ]] && cat "$archive_log"
    } | sort -u > "$done_list"
    comm -23 <(sort -u "$remote_list") "$done_list"
    rm -f "$done_list"
}

main() {
    echo "main: not yet implemented (see Task 4)"
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
```

- [ ] **Step 2: Make it executable**

```bash
chmod +x scripts/fetch_batch.sh
```

- [ ] **Step 3: Write the test script**

```bash
#!/usr/bin/env bash
set -euo pipefail

# shellcheck source=fetch_batch.sh
source "$(dirname "$0")/fetch_batch.sh"

fail=0

# --- group_key ---
result=$(group_key "abm_Cb=0.03_Cs=0.0_L=1000_R=1_U=10_dt=0.1_mot=RR_λ=0.1.csv")
expected="abm_Cb=0.03_L=1000_R=1_U=10_dt=0.1_mot=RR_λ=0.1.csv"
if [[ "$result" != "$expected" ]]; then
    echo "FAIL group_key: got '$result', expected '$expected'"
    fail=1
else
    echo "PASS group_key"
fi

# --- compute_pending ---
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

mkdir -p "$tmpdir/local"
# fileA (Cs=0) and fileB (Cs=1) are one group, not yet fetched
# fileC is already staged locally
# fileD is already archived
cat > "$tmpdir/remote.txt" <<'EOF'
abm_Cb=0.03_Cs=0.0_L=1000.csv
abm_Cb=0.03_Cs=1.0_L=1000.csv
abm_Cb=0.05_Cs=0.0_L=1000.csv
abm_Cb=0.07_Cs=0.0_L=1000.csv
EOF
touch "$tmpdir/local/abm_Cb=0.05_Cs=0.0_L=1000.csv"
echo "abm_Cb=0.07_Cs=0.0_L=1000.csv" > "$tmpdir/archived.txt"

result=$(compute_pending "$tmpdir/remote.txt" "$tmpdir/local" "$tmpdir/archived.txt" | sort)
expected=$(printf '%s\n%s' \
    "abm_Cb=0.03_Cs=0.0_L=1000.csv" \
    "abm_Cb=0.03_Cs=1.0_L=1000.csv")

if [[ "$result" != "$expected" ]]; then
    echo "FAIL compute_pending: got:"
    echo "$result"
    echo "expected:"
    echo "$expected"
    fail=1
else
    echo "PASS compute_pending"
fi

exit $fail
```

- [ ] **Step 4: Make it executable and run it**

```bash
chmod +x scripts/test_fetch_batch.sh
./scripts/test_fetch_batch.sh
```

Expected output:
```
PASS group_key
PASS compute_pending
```

- [ ] **Step 5: Commit**

```bash
git add scripts/fetch_batch.sh scripts/test_fetch_batch.sh
git commit -m "Add grouping/diff logic for fetch_batch.sh with tests"
```

---

### Task 4: `fetch_batch.sh` — remote listing, grouping, budgeted rsync

**Files:**
- Modify: `scripts/fetch_batch.sh`

This task fills in `main()`. It is not covered by `test_fetch_batch.sh` because
it requires cluster/VPN access; verification is manual (Task 5).

- [ ] **Step 1: Replace the placeholder `main` with the full implementation**

```bash
REMOTE_HOST="rfoffi@euler.ethz.ch"
REMOTE_DIR="/cluster/scratch/rfoffi/ResidenceTimes/abm"
PW_FILE="$HOME/Documents/clusterpw"
LOCAL_DIR="data/abm"
ARCHIVE_LOG="data/abm_archived.txt"
BUDGET_BYTES=$((8 * 1024 * 1024 * 1024))

# group_files_by_key: read filenames from stdin, print "<group_key>\t<filename>"
# lines, one per input filename.
group_files_by_key() {
    while IFS= read -r f; do
        [[ -z "$f" ]] && continue
        printf '%s\t%s\n' "$(group_key "$f")" "$f"
    done
}

main() {
    local_size=$(du -sb "$LOCAL_DIR" 2>/dev/null | cut -f1)
    local_size=${local_size:-0}

    if (( local_size >= BUDGET_BYTES )); then
        echo "data/abm/ is already at or above the ${BUDGET_BYTES}-byte budget (current: ${local_size} bytes)."
        echo "Run scripts/process_simulations.jl first to free up space."
        exit 0
    fi

    remote_list=$(mktemp)
    sshpass -f "$PW_FILE" ssh "$REMOTE_HOST" "ls $REMOTE_DIR" > "$remote_list"

    pending=$(mktemp)
    compute_pending "$remote_list" "$LOCAL_DIR" "$ARCHIVE_LOG" > "$pending"

    if [[ ! -s "$pending" ]]; then
        echo "Nothing pending."
        rm -f "$remote_list" "$pending"
        exit 0
    fi

    # remote file sizes, "<size>\t<filename>"
    sizes=$(mktemp)
    sshpass -f "$PW_FILE" ssh "$REMOTE_HOST" \
        "cd $REMOTE_DIR && du -b $(paste -sd' ' "$pending")" \
        | sed -E 's#^([0-9]+)\t.*/([^/]+)$#\1\t\2#' > "$sizes"

    grouped=$(mktemp)
    group_files_by_key < "$pending" | sort > "$grouped"

    fetched_files=()
    fetched_groups=0
    fetched_bytes=0
    skipped_groups=0

    for key in $(cut -f1 "$grouped" | uniq); do
        files=$(awk -F'\t' -v k="$key" '$1==k {print $2}' "$grouped")
        group_bytes=0
        while IFS= read -r f; do
            sz=$(awk -F'\t' -v fn="$f" '$2==fn {print $1}' "$sizes")
            group_bytes=$((group_bytes + ${sz:-0}))
        done <<< "$files"

        if (( local_size + fetched_bytes + group_bytes > BUDGET_BYTES )); then
            skipped_groups=$((skipped_groups + 1))
            continue
        fi

        while IFS= read -r f; do
            fetched_files+=("${REMOTE_HOST}:${REMOTE_DIR}/${f}")
        done <<< "$files"
        fetched_bytes=$((fetched_bytes + group_bytes))
        fetched_groups=$((fetched_groups + 1))
    done

    if (( ${#fetched_files[@]} > 0 )); then
        sshpass -f "$PW_FILE" rsync -avz -e ssh "${fetched_files[@]}" "$LOCAL_DIR/"
    fi

    total_pending_groups=$(cut -f1 "$grouped" | uniq | wc -l)
    echo "Fetched ${fetched_groups} group(s), ${#fetched_files[@]} file(s), $((fetched_bytes / 1024 / 1024)) MB."
    echo "Skipped ${skipped_groups} group(s) due to budget (of ${total_pending_groups} pending groups)."

    rm -f "$remote_list" "$pending" "$sizes" "$grouped"
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
```

Note this replaces the earlier placeholder `main` AND the
`if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then main "$@"; fi` block at the end of
the file from Task 3 — the new content above includes its own copy of that
guard, so the file should end with this version (don't duplicate it).

- [ ] **Step 2: Run the existing tests again to make sure the helper additions didn't break `group_key`/`compute_pending`**

```bash
./scripts/test_fetch_batch.sh
```

Expected output (unchanged):
```
PASS group_key
PASS compute_pending
```

- [ ] **Step 3: Commit**

```bash
git add scripts/fetch_batch.sh
git commit -m "Implement remote listing, grouping, and budgeted rsync in fetch_batch.sh"
```

---

### Task 5: Manual end-to-end verification

**Files:** none (verification only)

This task requires VPN access to the ETH network and the cluster password at
`~/Documents/clusterpw`. Run it interactively, not as part of automated CI.

- [ ] **Step 1: Confirm the script refuses to fetch while over budget**

The local `data/abm/` currently holds ~39GB (255 backlog files), well over the
8GB budget.

```bash
./scripts/fetch_batch.sh
```

Expected: prints the "already at or above the budget... run
process_simulations.jl first" message and exits 0 without contacting the
cluster.

- [ ] **Step 2: Process the existing backlog**

In the persistent Julia session for this repo:

```julia
using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
include("scripts/process_simulations.jl")
```

Expected: `data/exposure/` and `data/rdf/` get populated, `data/abm/` empties
out (files moved to `/media/Elements/ResidenceTimes/data/abm/`), and
`data/abm_archived.txt` now has 255 lines:

```bash
wc -l data/abm_archived.txt
ls data/abm | wc -l
ls /media/Elements/ResidenceTimes/data/abm | wc -l
```

- [ ] **Step 3: Run a real fetch**

```bash
./scripts/fetch_batch.sh
```

Expected: prints something like `Fetched N group(s), M file(s), ~XXXX MB.` and
`data/abm/` now contains those files, each under the 8GB total.

```bash
du -sh data/abm
```

- [ ] **Step 4: Confirm re-running fetch doesn't refetch the same files**

```bash
./scripts/fetch_batch.sh
```

Expected: since `data/abm/` is now non-empty but likely still under 8GB, this
either fetches the *next* batch of pending groups, or (if the fetched files
already cover all remaining pending groups) prints `Nothing pending.` Either
way, none of the files fetched in Step 3 should be fetched again — confirm by
checking `data/abm/` doesn't contain duplicates:

```bash
ls data/abm | sort | uniq -d
```

Expected: empty output (no duplicates).

- [ ] **Step 5: Run process_simulations.jl on the new batch and confirm the archive log grows**

```julia
include("scripts/process_simulations.jl")
```

```bash
wc -l data/abm_archived.txt   # should be > 255 now
ls data/abm                    # should be empty again
```

---

## Self-review notes

- All three spec components (archive log, process_simulations.jl edit,
  fetch_batch.sh) are covered by Tasks 1, 2, and 3+4 respectively.
- The 8GB budget, group-by-non-Cs-params logic, and `remote − archived − local`
  pending computation are all implemented and unit-tested where feasible
  (Task 3) or manually verified against the real cluster (Task 5).
- The `.gitignore` exception and design doc were already committed during
  brainstorming; Task 1 just adds the actual (empty) file.
