#!/usr/bin/env bash
# Cron wrapper for the NOSCENDO/lib2 pipeline (Snakefile_SE).
# Mirrors ~/enter (real $HOME over quota); runs the cub env via PATH.
# Emails are sent by the Snakefile's notify() on real-run events only.
set -uo pipefail

mkdir -p /mnt/san/microbio/conda/cache
export XDG_CACHE_HOME=/mnt/san/microbio/conda/cache
export HOME=/mnt/san/microbio/conda/0143813a
export PATH=/mnt/san/microbio/conda/envs/cub/bin:$PATH

# Startup CWD must be writable by lourdes: cron starts in the over-quota real home,
# and the shared cfdna_metagenomics/.snakemake is owned by another user (pe00936a).
# The Snakefile itself chdir's into each run's own folder for real work.
WORKDIR="$HOME/.cron_se_work"
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1

exec snakemake -s /mnt/san/microbio/apps/cub/workflow/Snakefile_SE \
  --configfile /mnt/san/microbio/apps/cub/config/config.yaml \
  --cores 50 --rerun-incomplete
