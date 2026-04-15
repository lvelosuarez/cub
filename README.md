# CUB — Clinical metagenomics pipeline
**CBAM (Centre Brestois Analyse Microbiota), CHRU Brest**  
Author: L. Velo Suarez — lourdes.velosuarez@chu-brest.fr

Snakemake pipeline for shotgun Illumina metagenomics from clinical samples.  
Handles three library types, two Kraken2 databases, and poly-T contamination tracking.

---

## Library types

| Type | Description | Snakefile |
|---|---|---|
| `lib1` | Paired-end 150bp × 2 (standard DNA/RNA) | `workflow/Snakefile` |
| `lib2` | Single-end 75bp (NOSCENDO) | `workflow/Snakefile_SE` |
| `lib3` | Paired-end 75bp + 10bp poly-A tag | `workflow/Snakefile_lib3` |

**lib3 design:** R1 carries the DNA fragment; R2 carries a poly-A tail added before library preparation, which appears as poly-T in sequencing. Reads where R2 is poly-T originated from the sample before library prep (genuine DNA). Reads where R2 is not poly-T were introduced during library preparation or from reagents (contamination). Both streams are kept and processed independently.

---

## Pipeline overview

```
raw FASTQs
    │
    ▼
QC              bbduk — adapter/artifact trimming, quality filtering
    │
    ▼ (lib3 only)
split_polyt     bbduk — split by R2 poly-T content → tagged / untagged streams
    │
    ▼
dedup           clumpify — optical/PCR duplicate removal
    │
    ▼
dehost          hostile + bowtie2 — human read removal
    │                               (index: human-t2t-hla-argos985-mycob140)
    │ (lib1 only: merge 4 lanes)
    ▼
kraken2         PlusPF database   → output/kraken2/{sample}.report + .pavian
kraken2         EuPathDB48        → output/kraken2_eupathdb/{sample}.report + .pavian
barlett         bowtie2 → GTDB    → output/bam/{sample}.bam
    │
    ▼
Nreads.csv      read counts at each stage
```

For lib3, kraken2 and barlett run on both `tagged` and `untagged` streams independently.

---

## Installation

```bash
# Clone the repository
git clone https://github.com/lvelosuarez/cub
cd cub

# Create conda environment
mamba env create -f cub.yaml
```

---

## Running an analysis

### 1. Mark the run folder with its library type

Create a plain-text `lib_type.txt` file inside the run folder:

```bash
echo "lib1" > /path/to/run_folder/lib_type.txt
# or
echo "lib2" > /path/to/run_folder/lib_type.txt
# or
echo "lib3" > /path/to/run_folder/lib_type.txt
```

### 2. Run the dispatcher

```bash
bash cub_dispatcher.sh /path/to/run_folder
```

The dispatcher reads `lib_type.txt`, resolves the correct FASTQ directory, and launches Snakemake. A timestamped log is written to the run folder.

**lib1 note:** FASTQs are expected under `Alignment_1/{timestamp}/Fastq/` — the dispatcher resolves this automatically.  
**lib2 note:** FASTQs are expected under `alignment_lou/` inside the run folder.  
**lib3 note:** FASTQs are expected in the run folder root; the `microbial/` subfolder is ignored.

### 3. Check outputs

All results are written to `output/` inside the working FASTQ directory.

| File | Description |
|---|---|
| `output/kraken2/{sample}.report` | Kraken2 report — PlusPF database |
| `output/kraken2/{sample}.pavian` | Pavian-compatible report — PlusPF |
| `output/kraken2_eupathdb/{sample}.report` | Kraken2 report — EuPathDB48 |
| `output/kraken2_eupathdb/{sample}.pavian` | Pavian-compatible report — EuPathDB48 |
| `output/bam/{sample}.bam` | Barlett (GTDB) alignment |
| `output/Nreads.csv` | Read counts at each pipeline stage |

For **lib3**, all files include `_tagged` or `_untagged` in the sample name:  
`output/kraken2/{sample}_tagged.report`, `output/bam/{sample}_untagged.bam`, etc.

`output/Nreads.csv` for lib3 includes: `raw, qc, tagged, untagged, tag_ratio_pct, dedup_tagged, dedup_untagged, dehost_tagged, dehost_untagged`

---

## Dry-run (recommended before first execution)

```bash
conda run -n cub snakemake \
    -s /path/to/cub/workflow/Snakefile_lib3 \
    --configfile /path/to/cub/config/config.yaml \
    --directory /path/to/run_folder \
    -n
```

Replace `Snakefile_lib3` with `Snakefile` or `Snakefile_SE` for lib1/lib2.

---

## Automated processing with cron

To process new runs automatically, set up a cron job that calls `cub_cron.sh` every 30 minutes:

```bash
crontab -e
```

Add:
```
*/30 * * * * /path/to/cub/cub_cron.sh >> /var/log/cub_cron.log 2>&1
```

Edit `cub_cron.sh` to add your run base directories to `SCAN_DIRS`. The cron job will process any folder containing `lib_type.txt` that has not already been completed (no `.cub.done` marker).

---

## Configuration

`config/config.yaml`:

```yaml
PlusPF:   "/mnt/san/microbio/fermion/tl/kraken2/PlusPF20251015/"
barlett:  "/mnt/san/microbio/fermion/tl/gtdb/bartlett_index/bartlett"
eupathdb: "/mnt/san/microbio/fermion/tl/kraken2/eupathdb48_20230407/"
threads:  15

lib3:
  polyt_literal: "TTTTTTTTT"   # poly-T sequence to detect in R2
  polyt_kmer:    9              # k-mer length (fits in 10bp poly-T region of 19bp R2)
  polyt_hdist:   1              # mismatches allowed (tolerates Illumina Q30 errors)
```

---

## Repository structure

```
cub/
├── workflow/
│   ├── Snakefile           # lib1 — paired-end 150bp
│   ├── Snakefile_SE        # lib2 — single-end 75bp
│   ├── Snakefile_lib3      # lib3 — paired-end 75bp + poly-T tag
│   └── rules/
│       ├── count.smk       # read counting for lib1/lib2
│       └── count_lib3.smk  # read counting for lib3 (per stream)
├── config/
│   └── config.yaml         # database paths and parameters
├── scr/
│   └── bonobo/             # statistical analysis suite (cub2.py, bam.py, post_analysis.py)
├── cub_dispatcher.sh        # launches the correct pipeline for a run folder
├── cub_cron.sh              # scans for new runs and dispatches them
└── cub.yaml                 # conda environment
```
