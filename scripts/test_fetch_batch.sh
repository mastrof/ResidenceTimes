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
