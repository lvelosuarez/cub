#!/usr/bin/env bash
# cub_cron.sh — scan for new sequencing runs and dispatch them
#
# Called by cron, e.g. every 30 minutes:
#   */30 * * * * /path/to/cub/cub_cron.sh >> /var/log/cub_cron.log 2>&1
#
# A run folder is processed when:
#   - It contains lib_type.txt           (marks it as a cub run)
#   - It does NOT contain .cub.done      (not yet completed)
#   - It does NOT contain .cub.lock      (not currently running)
#
# On success a .cub.done marker is written so the folder is not re-processed.

set -euo pipefail

SCRIPT_DIR="$(dirname "$(realpath "$0")")"
DISPATCHER="${SCRIPT_DIR}/cub_dispatcher.sh"

# ── Configure base directories to scan ─────────────────────────────────────
# Add all root directories where sequencing runs land.
SCAN_DIRS=(
    "/mnt/san/microbio/m_clinique"
    "/mnt/san/microbio/m_clinique/NOSCENDO"
    "/mnt/nas/R60-011/BACTERIOLOGIE/metagenomique_NextSeq/NOSCENDO/output"
)

echo "=== cub_cron start: $(date) ==="

for base_dir in "${SCAN_DIRS[@]}"; do
    [[ -d "$base_dir" ]] || continue

    for run_folder in "${base_dir}"/*/; do
        [[ -d "$run_folder" ]]               || continue
        [[ -f "${run_folder}/lib_type.txt" ]] || continue   # not a cub run
        [[ -f "${run_folder}/.cub.done" ]]   && continue   # already done
        [[ -f "${run_folder}/.cub.lock" ]]   && continue   # currently running

        echo "INFO: Dispatching ${run_folder}"
        if "$DISPATCHER" "$run_folder"; then
            touch "${run_folder}/.cub.done"
            echo "INFO: Completed ${run_folder}"
        else
            echo "ERROR: Failed ${run_folder} (exit $?)" >&2
        fi
    done
done

echo "=== cub_cron end: $(date) ==="
