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
