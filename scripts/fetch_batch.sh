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
