#!/usr/bin/env bash
# cub_dispatcher.sh — detect library type and launch the correct Snakefile
#
# Usage: cub_dispatcher.sh <run_folder>
#
# Each run folder must contain a plain-text file 'lib_type.txt' with one of:
#   lib1  — paired-end 150bp (standard)
#   lib2  — single-end 75bp  (NOSCENDO)
#   lib3  — paired-end 75bp + poly-T tag
#
# The script:
#   1. Reads lib_type.txt
#   2. Resolves the correct working directory for Snakemake
#   3. Acquires a lock to prevent concurrent runs on the same folder
#   4. Runs Snakemake with --rerun-incomplete (safe to call repeatedly)
#   5. Writes a timestamped log inside the run folder

set -euo pipefail

RUN_FOLDER="${1:?Usage: cub_dispatcher.sh <run_folder>}"
RUN_FOLDER="$(realpath "$RUN_FOLDER")"
PIPELINE_DIR="$(dirname "$(dirname "$(realpath "$0")")")"
CONDA_ENV="cub"
CONFIGFILE="${PIPELINE_DIR}/config/config.yaml"

# ── 1. Read library type ────────────────────────────────────────────────────
LIB_TYPE_FILE="${RUN_FOLDER}/lib_type.txt"
if [[ ! -f "$LIB_TYPE_FILE" ]]; then
    echo "ERROR: ${LIB_TYPE_FILE} not found." \
         "Create it with content: lib1, lib2, or lib3" >&2
    exit 1
fi
LIB_TYPE=$(tr -d '[:space:]' < "$LIB_TYPE_FILE")

# ── 2. Resolve Snakefile and working directory ──────────────────────────────
case "$LIB_TYPE" in
    lib1)
        SNAKEFILE="${PIPELINE_DIR}/workflow/Snakefile"
        # FASTQs live under Alignment_1/{timestamp}/Fastq/ — resolve dynamically
        FASTQ_DIR=$(ls -d "${RUN_FOLDER}/Alignment_1"/*/Fastq 2>/dev/null | sort | tail -1)
        if [[ -z "$FASTQ_DIR" ]]; then
            echo "ERROR: Could not find Alignment_1/*/Fastq/ under ${RUN_FOLDER}" >&2
            exit 1
        fi
        WORK_DIR="$FASTQ_DIR"
        ;;
    lib2)
        SNAKEFILE="${PIPELINE_DIR}/workflow/Snakefile_SE"
        WORK_DIR="${RUN_FOLDER}/alignment_lou"
        if [[ ! -d "$WORK_DIR" ]]; then
            echo "ERROR: alignment_lou/ not found under ${RUN_FOLDER}" >&2
            exit 1
        fi
        ;;
    lib3)
        SNAKEFILE="${PIPELINE_DIR}/workflow/Snakefile_lib3"
        WORK_DIR="$RUN_FOLDER"
        ;;
    *)
        echo "ERROR: Unknown lib_type '${LIB_TYPE}'. Must be lib1, lib2, or lib3." >&2
        exit 1
        ;;
esac

# ── 3. Lock — prevent concurrent runs on the same folder ───────────────────
LOCK="${RUN_FOLDER}/.cub.lock"
if [[ -f "$LOCK" ]]; then
    echo "INFO: Run folder is locked (another instance may be running): ${LOCK}"
    exit 0
fi
touch "$LOCK"
trap 'rm -f "$LOCK"' EXIT

# ── 4. Run Snakemake ────────────────────────────────────────────────────────
LOG="${RUN_FOLDER}/cub_$(date +%Y%m%d_%H%M%S).log"
echo "INFO: lib_type=${LIB_TYPE}  work_dir=${WORK_DIR}  log=${LOG}"

conda run -n "$CONDA_ENV" \
    snakemake \
        -s "$SNAKEFILE" \
        --configfile "$CONFIGFILE" \
        --directory "$WORK_DIR" \
        --cores all \
        --rerun-incomplete \
        --latency-wait 60 \
        --printshellcmds \
    2>&1 | tee "$LOG"

echo "INFO: Done — $(date)"
