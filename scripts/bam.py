#!/usr/bin/env python3
"""
bam.py

Compare per-reference read counts across BAM files.

Statistics logic (harmonized with cub.py):
- If --neg-name is provided:
    * use negative-control–based statistics (per reference)
- Else:
    * use cohort-based statistics only if a reference is seen in >= 3 samples
    * otherwise z_score and p_value are set to null (NA)

Other features:
- Counts paired reads per reference with a MAPQ threshold.
- Can ignore specified samples.
- Annotates references with GTDB taxonomy table.
"""

import argparse
import math
import sys
from pathlib import Path

import polars as pl
import pysam
from scipy.stats import norm

# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Compare BAM files per reference, compute stats versus a negative "
            "control if provided, or across samples otherwise, and annotate "
            "references with GTDB taxonomy."
        )
    )
    p.add_argument(
        "--bam-dir",
        required=True,
        help="Directory containing BAM files.",
    )
    p.add_argument(
        "--output",
        required=True,
        help="Output TSV file.",
    )
    p.add_argument(
        "--taxonomy",
        required=True,
        help=("GTDB taxonomy table (.rds, .tsv, .txt or .csv) with columns: accessions, gtdb_taxonomy, ncbi_taxid."),
    )
    p.add_argument(
        "--neg-name",
        default=None,
        help=("Negative control BAM filename (e.g. 'MG-Run35-Ech-11.bam'). If provided, stats are computed versus this sample."),
    )
    p.add_argument(
        "--ignore-samples",
        default=None,
        help="Text file with one BAM filename per line to ignore.",
    )
    p.add_argument(
        "--min-mapq",
        type=int,
        default=25,
        help="Minimum mapping quality to count reads (default: 25).",
    )
    p.add_argument(
        "--pseudocount",
        type=float,
        default=1.0,
        help="Pseudocount added to negative counts (default: 1.0).",
    )
    p.add_argument(
        "--coverage-dir",
        default=None,
        help=(
            "Directory containing samtools coverage output files "
            "({sample}_coverage.txt). If provided, coverage_breadth and "
            "mean_depth columns are added to the output."
        ),
    )
    return p.parse_args()


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------


def load_taxonomy(taxonomy_path: Path) -> pl.DataFrame:
    """
    Load taxonomy table from .rds, .tsv, .txt or .csv into a Polars DataFrame.

    Expected columns:
      - accessions
      - gtdb_taxonomy
      - optionally: ncbi_taxid, gtdb_genome_representative, gtdb_representative
    """
    suffix = taxonomy_path.suffix.lower()

    if suffix == ".rds":
        try:
            import pyreadr  # type: ignore
        except ImportError as err:
            raise SystemExit("ERROR: pyreadr is required to read .rds files.\nInstall it with: pip install pyreadr") from err

        result = pyreadr.read_r(str(taxonomy_path))
        if not result:
            raise SystemExit(f"ERROR: No object found in RDS file: {taxonomy_path}")
        pandas_df = next(iter(result.values()))
        tax_df = pl.from_pandas(pandas_df)

    elif suffix in {".tsv", ".txt"}:
        tax_df = pl.read_csv(taxonomy_path, separator="\t")

    elif suffix == ".csv":
        tax_df = pl.read_csv(taxonomy_path)

    else:
        raise SystemExit(f"ERROR: Unsupported taxonomy file extension: {suffix}\nUse .rds, .tsv, .txt or .csv")

    if "accessions" not in tax_df.columns:
        raise SystemExit("ERROR: taxonomy table must contain a column named 'accessions'.")

    keep_cols = [
        c
        for c in [
            "accessions",
            "gtdb_taxonomy",
            "ncbi_taxid",
            "gtdb_genome_representative",
            "gtdb_representative",
        ]
        if c in tax_df.columns
    ]

    return tax_df.select(keep_cols).unique(subset=["accessions"])


def count_paired_reads_per_reference(bam_path: Path, min_mapq: int) -> dict[str, int]:
    """
    Count reads per reference.

    - Paired-end BAMs: count read1 of proper pairs passing MAPQ threshold.
    - Single-end BAMs: count all mapped reads passing MAPQ threshold.
    """
    counts: dict[str, int] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            if read.is_paired:
                if not read.is_read1 or read.mate_is_unmapped:
                    continue
            if read.mapping_quality < min_mapq:
                continue
            ref_name = bam.get_reference_name(read.reference_id)
            counts[ref_name] = counts.get(ref_name, 0) + 1
    return counts


def load_coverage(coverage_dir: Path) -> pl.DataFrame:
    """
    Load all samtools coverage TSVs from *coverage_dir*.

    samtools coverage columns (tab-separated, with header):
      #rname  startpos  endpos  numreads  covbases  coverage  meandepth  meanbaseq  meanmapq

    Returns a DataFrame with columns: sample, rname, coverage_breadth, mean_depth
    """
    frames = []
    for txt in sorted(coverage_dir.glob("*_coverage.txt")):
        # Derive sample name: strip trailing _coverage
        stem = txt.stem  # e.g. "BG_tagged_coverage" or "sample1_coverage"
        sample = stem[: -len("_coverage")] if stem.endswith("_coverage") else stem
        try:
            cov = pl.read_csv(txt, separator="\t", comment_prefix="#")
        except Exception as exc:
            print(f"[WARN] Could not read {txt}: {exc}", file=sys.stderr)
            continue
        # samtools coverage writes the header with leading '#'; polars strips it
        # Rename columns defensively in case the '#' prefix is kept
        col_map = {}
        for c in cov.columns:
            clean = c.lstrip("#")
            if clean != c:
                col_map[c] = clean
        if col_map:
            cov = cov.rename(col_map)

        needed = {"rname", "coverage", "meandepth"}
        if not needed.issubset(set(cov.columns)):
            print(f"[WARN] Unexpected columns in {txt}: {cov.columns}", file=sys.stderr)
            continue

        frames.append(
            cov.select([
                pl.lit(sample).alias("sample"),
                pl.col("rname"),
                pl.col("coverage").cast(pl.Float64).alias("coverage_breadth"),
                pl.col("meandepth").cast(pl.Float64).alias("mean_depth"),
            ])
        )

    if not frames:
        return pl.DataFrame({"sample": [], "rname": [], "coverage_breadth": [], "mean_depth": []})

    return pl.concat(frames)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------


def main():
    args = parse_args()

    bam_dir = Path(args.bam_dir)
    if not bam_dir.is_dir():
        sys.exit(f"--bam-dir does not exist or is not a directory: {bam_dir}")

    bam_files = sorted(bam_dir.glob("*.bam"))
    if not bam_files:
        sys.exit(f"No BAM files (*.bam) found in {bam_dir}")

    # Handle ignore list
    ignore_set: set[str] = set()
    if args.ignore_samples:
        ignore_path = Path(args.ignore_samples)
        if not ignore_path.exists():
            sys.exit(f"--ignore-samples file not found: {ignore_path}")
        for line in ignore_path.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            ignore_set.add(Path(line).stem)

    # Normalize negative control name (compare by stem)
    use_negative = args.neg_name is not None
    neg_sample_stem: str | None = None
    if use_negative:
        neg_sample_stem = Path(args.neg_name).stem

    # ------------------------------------------------------------------
    # Count per-reference reads for all samples
    # ------------------------------------------------------------------

    rows: list[dict[str, object]] = []
    for bam_path in bam_files:
        sample_id = bam_path.stem  # e.g. CasClinique3_ADN
        if sample_id in ignore_set:
            print(f"[INFO] Ignoring sample: {sample_id}", file=sys.stderr)
            continue

        print(f"[INFO] Counting paired reads for sample: {sample_id}", file=sys.stderr)
        ref_counts = count_paired_reads_per_reference(bam_path, args.min_mapq)
        for ref, n in ref_counts.items():
            rows.append(
                {
                    "sample": sample_id,
                    "reference": ref,
                    "read_count": n,
                }
            )

    if not rows:
        print("[WARN] No reads counted in any BAM; writing empty output.", file=sys.stderr)
        empty = pl.DataFrame(
            {
                "sample": [],
                "reference": [],
                "mapq_min": [],
                "read_count": [],
                "neg_read_count": [],
                "z_score": [],
                "p_value": [],
                "stat_mode": [],
                "coverage_breadth": [],
                "mean_depth": [],
                "gtdb_taxonomy": [],
                "ncbi_taxid": [],
            }
        )
        empty.write_csv(args.output, separator="\t")
        return

    df = pl.DataFrame(rows)
    df = df.with_columns(pl.lit(args.min_mapq).alias("mapq_min"))

    # ------------------------------------------------------------------
    # Load taxonomy and annotate references
    # ------------------------------------------------------------------

    taxonomy_path = Path(args.taxonomy)
    if not taxonomy_path.exists():
        sys.exit(f"Taxonomy file not found: {taxonomy_path}")

    print(f"[INFO] Loading taxonomy from: {taxonomy_path}", file=sys.stderr)
    tax_df = load_taxonomy(taxonomy_path)
    print(
        f"[INFO] Taxonomy table loaded with {tax_df.height} accessions; columns: {', '.join(tax_df.columns)}",
        file=sys.stderr,
    )

    # Join on accession / reference
    df = df.join(
        tax_df.select([c for c in ["accessions", "gtdb_taxonomy", "ncbi_taxid"] if c in tax_df.columns]),
        left_on="reference",
        right_on="accessions",
        how="left",
    ).drop("accessions", strict=False)

    # ------------------------------------------------------------------
    # Statistics
    # ------------------------------------------------------------------

    # Initialize neg_read_count for compatibility
    df = df.with_columns(pl.lit(0.0).alias("neg_read_count"))

    if use_negative:
        # -----------------------------
        # Negative-control–based stats
        # -----------------------------
        assert neg_sample_stem is not None
        if neg_sample_stem not in df.select("sample").to_series().unique():
            print(
                f"[WARN] Negative control sample '{neg_sample_stem}' not found among samples.",
                file=sys.stderr,
            )

        neg_df = df.filter(pl.col("sample") == neg_sample_stem)

        neg_bg = neg_df.group_by("reference").agg(pl.col("read_count").sum().alias("neg_read_count"))

        df = df.join(neg_bg, on="reference", how="left")
        df = df.with_columns(pl.col("neg_read_count").fill_null(0.0).alias("neg_read_count"))

        df = df.with_columns((pl.col("neg_read_count") + args.pseudocount).alias("mu"))
        df = df.with_columns(((pl.col("read_count") - pl.col("mu")) / pl.col("mu").sqrt()).alias("z_score"))

        df = df.with_columns(
            pl.col("z_score")
            .map_elements(
                lambda z: None if z is None or (isinstance(z, float) and math.isnan(z)) else 1.0 - float(norm.cdf(z)),
                return_dtype=pl.Float64,
            )
            .alias("p_value"),
            pl.lit("negative").alias("stat_mode"),
        )

    else:
        # -----------------------------
        # Cohort-based stats (n >= 3)
        # -----------------------------
        bg = df.group_by("reference").agg(
            pl.col("read_count").mean().alias("mean_count"),
            pl.col("read_count").std().alias("sd_count"),
            (pl.col("read_count") > 0).sum().alias("n_samples_with_ref"),
        )

        df = df.join(bg, on="reference", how="left")

        df = df.with_columns(
            ((pl.col("n_samples_with_ref") >= 3) & (pl.col("sd_count").is_not_null()) & (pl.col("sd_count") > 0)).alias("has_valid_stats")
        )

        df = df.with_columns(
            pl.when(pl.col("has_valid_stats"))
            .then((pl.col("read_count") - pl.col("mean_count")) / pl.col("sd_count"))
            .otherwise(None)
            .alias("z_score")
        )

        df = df.with_columns(
            pl.col("z_score")
            .map_elements(
                lambda z: None if z is None or (isinstance(z, float) and math.isnan(z)) else 1.0 - float(norm.cdf(z)),
                return_dtype=pl.Float64,
            )
            .alias("p_value")
        )

        df = df.with_columns(pl.when(pl.col("has_valid_stats")).then(pl.lit("cohort")).otherwise(pl.lit("none")).alias("stat_mode"))

    # ------------------------------------------------------------------
    # Join coverage (optional)
    # ------------------------------------------------------------------

    if args.coverage_dir:
        cov_dir = Path(args.coverage_dir)
        if not cov_dir.is_dir():
            print(f"[WARN] --coverage-dir does not exist: {cov_dir}", file=sys.stderr)
        else:
            print(f"[INFO] Loading coverage from: {cov_dir}", file=sys.stderr)
            cov_df = load_coverage(cov_dir)
            if cov_df.height > 0:
                df = df.join(
                    cov_df,
                    left_on=["sample", "reference"],
                    right_on=["sample", "rname"],
                    how="left",
                )
            else:
                print("[WARN] No coverage data loaded; skipping join.", file=sys.stderr)

    # ------------------------------------------------------------------
    # Final formatting & output
    # ------------------------------------------------------------------

    df = df.sort(["sample", "read_count"], descending=[False, True])

    out_cols = [
        "sample",
        "reference",
        "mapq_min",
        "read_count",
        "neg_read_count",
        "z_score",
        "p_value",
        "stat_mode",
        "coverage_breadth",
        "mean_depth",
        "gtdb_taxonomy",
        "ncbi_taxid",
    ]
    out_cols = [c for c in out_cols if c in df.columns]

    df.select(out_cols).write_csv(args.output, separator="\t")
    print(f"[INFO] Wrote: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
