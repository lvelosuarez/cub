#!/usr/bin/env python3
"""
cub2.py

Statistical analysis of Kraken2 reports using minimizer data.

Features
--------
- Reads Kraken2 standard reports produced with --report-minimizer-data.
- Reads Kraken2 DB inspect.txt (kraken2-inspect output).
- Host organism exclusion via NCBI taxonomy tree (names.dmp + nodes.dmp
  from the DB taxonomy/ directory) with fallback to simple name/taxID match.
- Species-level filtering: non-viruses kept at rank S; viruses at S/S1/S2
  where direct reads dominate (fragments_direct / fragments_clade >= threshold).
- Minimizer proportion per sample computed on non-host taxa only.
- genome_coverage_estimate = distinct_minimizers / db_taxon_minimizers.
- Statistics: cohort mode (n >= 3 samples) or negative-control mode.
- Benjamini–Hochberg FDR correction.
- Statistics logic harmonized with bam.py.

Statistics logic
----------------
- If --negative-control is provided:
    * z-score based on negative-control background minimizer proportion.
- Else:
    * cohort z-score if n_samples_with_taxon >= 3 and sd > 0; else null.
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import norm

# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------


def parse_args():
    p = argparse.ArgumentParser(
        description="Statistical analysis of Kraken2 reports using minimizer data."
    )
    p.add_argument(
        "--std-reports",
        required=True,
        help="Directory containing Kraken2 standard reports (--report-minimizer-data).",
    )
    p.add_argument(
        "--organism",
        default=None,
        help=(
            "Host organism: scientific name (e.g. 'Homo sapiens') or taxID (e.g. '9606'). "
            "Host taxa are excluded from minimizer proportion denominators and statistics."
        ),
    )
    p.add_argument(
        "--reference",
        required=True,
        help="Kraken2 DB inspect.txt (output of kraken2-inspect).",
    )
    p.add_argument(
        "--domain",
        default=None,
        help="Comma-separated list of domains of interest (currently unused).",
    )
    p.add_argument(
        "--metadata",
        default=None,
        help="Optional metadata CSV with sample-level information (currently unused).",
    )
    p.add_argument(
        "--sample-col",
        default="sample",
        help="Column in metadata containing sample IDs (currently unused).",
    )
    p.add_argument(
        "--columns",
        default=None,
        help="Comma-separated list of metadata columns to keep (currently unused).",
    )
    p.add_argument(
        "--samples-to-remove",
        default=None,
        help="Text file with one sample ID per line to exclude.",
    )
    p.add_argument(
        "--prefix",
        default="cub",
        help="Prefix for output files (default: cub).",
    )
    p.add_argument(
        "--outdir",
        required=True,
        help="Output directory.",
    )
    p.add_argument(
        "--min-prop",
        type=float,
        default=0.0,
        help="Minimum non-host minimizer proportion to keep rows (default: 0.0).",
    )
    p.add_argument(
        "--negative-control",
        default=None,
        help=(
            "Sample name corresponding to a negative control. "
            "If provided, background is estimated from this sample."
        ),
    )
    p.add_argument(
        "--pseudocount",
        type=float,
        default=1e-6,
        help="Pseudocount for negative-control statistics (default: 1e-6).",
    )
    p.add_argument(
        "--virus-direct-ratio",
        type=float,
        default=0.8,
        help=(
            "For virus taxa, minimum fragments_direct / fragments_clade ratio "
            "to retain the node (default: 0.8). Keeps only the deepest node "
            "with real read support."
        ),
    )
    return p.parse_args()


# ---------------------------------------------------------------------
# Inspect.txt reader
# ---------------------------------------------------------------------


def read_inspect(path: Path) -> pl.DataFrame:
    """
    Read kraken2-inspect output into a Polars DataFrame.

    Format (tab-separated, no header):
        db_pct  db_clade_minimizers  db_taxon_minimizers  db_rank_code  taxid  db_name
    """
    df = pl.read_csv(path, separator="\t", has_header=False)

    if df.width < 6:
        raise SystemExit(
            f"ERROR: Expected at least 6 columns in inspect file: {path}, got {df.width}"
        )

    cols = df.columns
    df = df.rename(
        {
            cols[0]: "db_pct",
            cols[1]: "db_clade_minimizers",
            cols[2]: "db_taxon_minimizers",
            cols[3]: "db_rank_code",
            cols[4]: "taxid",
            cols[5]: "db_name",
        }
    )
    return df.with_columns(pl.col("taxid").cast(pl.Utf8))


# ---------------------------------------------------------------------
# Report reader
# ---------------------------------------------------------------------


def read_std_report(path: Path, sample: str) -> pl.DataFrame:
    """
    Read a Kraken2 standard report (--report-minimizer-data) into Polars.

    Columns: pct_fragments, fragments_clade, fragments_direct,
             minimizers, distinct_minimizers, rank_code, taxid, name
    """
    df = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=[
            "pct_fragments",
            "fragments_clade",
            "fragments_direct",
            "minimizers",
            "distinct_minimizers",
            "rank_code",
            "taxid",
            "name",
        ],
    )
    return df.with_columns(
        [
            pl.col("taxid").cast(pl.Utf8),
            pl.col("name").str.strip_chars(),
            pl.lit(sample).alias("sample"),
        ]
    )


# ---------------------------------------------------------------------
# NCBI taxonomy tree (names.dmp + nodes.dmp)
# ---------------------------------------------------------------------


def _get_taxonomy_paths(reference_path: Path):
    """Return (names_path, nodes_path) from the DB directory."""
    taxonomy_dir = reference_path.parent / "taxonomy"
    return taxonomy_dir / "names.dmp", taxonomy_dir / "nodes.dmp"


def _find_host_taxid(organism: str, names_path: Path) -> str | None:
    """Resolve organism scientific name to a taxid via names.dmp."""
    target = organism.strip().lower()
    with open(names_path) as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue
            tax_id, name_txt, _, name_class = parts[:4]
            if name_class == "scientific name" and name_txt.lower() == target:
                return str(tax_id)
    return None


def _load_parent_children(nodes_path: Path) -> dict[str, list[str]]:
    """Parse nodes.dmp into a parent → children mapping."""
    tree: dict[str, list[str]] = {}
    with open(nodes_path) as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 2:
                continue
            tax_id, parent_id = str(parts[0]), str(parts[1])
            tree.setdefault(parent_id, []).append(tax_id)
    return tree


def _get_descendants(root: str, tree: dict[str, list[str]]) -> set[str]:
    """BFS over taxonomy tree, returns all descendants including root."""
    seen = {root}
    stack = [root]
    while stack:
        node = stack.pop()
        for child in tree.get(node, []):
            if child not in seen:
                seen.add(child)
                stack.append(child)
    return seen


def resolve_host_taxids(organism: str, reference_path: Path) -> set[str] | None:
    """
    Return the set of all taxids in the host clade (including root).

    Strategy:
    1. If organism is all-digits → treat as taxID directly.
    2. Otherwise look up names.dmp in the DB taxonomy/ directory.
    3. Walk nodes.dmp to collect all descendants.
    4. Falls back to None (simple name match) if taxonomy files are absent.
    """
    if organism is None:
        return None

    organism = organism.strip()
    names_path, nodes_path = _get_taxonomy_paths(reference_path)

    if not names_path.exists() or not nodes_path.exists():
        print(
            "[WARN] taxonomy/names.dmp or nodes.dmp not found — "
            "falling back to simple name/taxID host matching.",
            file=sys.stderr,
        )
        return None

    if organism.isdigit():
        root_taxid = organism
    else:
        root_taxid = _find_host_taxid(organism, names_path)
        if root_taxid is None:
            print(
                f"[WARN] Could not find '{organism}' in names.dmp — "
                "falling back to simple name host matching.",
                file=sys.stderr,
            )
            return None

    tree = _load_parent_children(nodes_path)
    taxids = _get_descendants(root_taxid, tree)
    print(
        f"[INFO] Host clade '{organism}' (taxid {root_taxid}): "
        f"{len(taxids)} taxids resolved from taxonomy tree.",
        file=sys.stderr,
    )
    return taxids


# ---------------------------------------------------------------------
# Host marking
# ---------------------------------------------------------------------


def mark_host_taxa(
    df: pl.DataFrame,
    organism: str | None,
    reference_path: Path,
) -> pl.DataFrame:
    """
    Add boolean column 'is_host'.

    Priority:
    1. Taxonomy-tree set (names.dmp + nodes.dmp) if available.
    2. Simple taxID match (if organism is all-digits).
    3. Case-insensitive name substring match.
    4. No organism → all False.
    """
    if organism is None:
        return df.with_columns(pl.lit(False).alias("is_host"))

    host_taxids = resolve_host_taxids(organism, reference_path)

    if host_taxids is not None:
        return df.with_columns(
            pl.col("taxid").is_in(list(host_taxids)).alias("is_host")
        )

    # Fallback
    org = organism.strip()
    if org.isdigit():
        return df.with_columns((pl.col("taxid") == org).alias("is_host"))
    return df.with_columns(
        pl.col("name").str.contains(f"(?i){org}", literal=False).alias("is_host")
    )


# ---------------------------------------------------------------------
# Species-level filtering
# ---------------------------------------------------------------------


def filter_to_best_kraken_taxa(
    df: pl.DataFrame,
    virus_direct_ratio: float = 0.8,
) -> pl.DataFrame:
    """
    Retain the most clinically relevant taxonomic level per organism.

    Non-viruses: keep only species level (rank_code == 'S').
    Viruses    : keep S / S1 / S2 where
                 fragments_direct / fragments_clade >= virus_direct_ratio.
                 This keeps the deepest node with real read support
                 (e.g. 'Human adenovirus 5' at S2) and drops ancestors.
    """
    is_virus = pl.col("name").str.contains("(?i)virus", literal=False)

    non_virus = df.filter(~is_virus).filter(pl.col("rank_code") == "S")

    virus = (
        df.filter(is_virus)
        .filter(pl.col("rank_code").is_in(["S", "S1", "S2"]))
        .with_columns(
            pl.when(pl.col("fragments_clade").cast(pl.Float64) > 0)
            .then(
                pl.col("fragments_direct").cast(pl.Float64)
                / pl.col("fragments_clade").cast(pl.Float64)
            )
            .otherwise(0.0)
            .alias("_direct_ratio")
        )
        .filter(pl.col("_direct_ratio") >= virus_direct_ratio)
        .drop("_direct_ratio")
    )

    return pl.concat([non_virus, virus])


# ---------------------------------------------------------------------
# Benjamini–Hochberg FDR
# ---------------------------------------------------------------------


def benjamini_hochberg(p_values: pl.Series) -> pl.Series:
    """
    Benjamini–Hochberg FDR correction.
    Null values are preserved as null in the output.
    """
    mask = p_values.is_not_null()
    valid_idx = [i for i, m in enumerate(mask.to_list()) if m]
    p_arr = np.array(p_values.drop_nulls().to_list(), dtype=float)

    n = len(p_arr)
    if n == 0:
        return p_values.cast(pl.Float64)

    order = np.argsort(p_arr)
    ranks = np.arange(1, n + 1, dtype=float)
    p_sorted = p_arr[order]
    p_adj = p_sorted * n / ranks
    # Enforce monotonicity (from right)
    p_adj = np.minimum.accumulate(p_adj[::-1])[::-1]
    p_adj = np.clip(p_adj, 0.0, 1.0)

    # Scatter back to original positions
    result_arr = np.full(len(p_values), None, dtype=object)
    for rank_pos, orig_pos in enumerate(order):
        result_arr[valid_idx[orig_pos]] = float(p_adj[rank_pos])

    return pl.Series(values=list(result_arr), dtype=pl.Float64)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    reference_path = Path(args.reference)
    use_negative = args.negative_control is not None

    # ------------------------------------------------------------------
    # Load inspect.txt
    # ------------------------------------------------------------------
    if not reference_path.exists():
        sys.exit(f"ERROR: --reference not found: {reference_path}")

    print(f"[INFO] Reading inspect.txt: {reference_path}", file=sys.stderr)
    inspect_df = read_inspect(reference_path)

    # ------------------------------------------------------------------
    # Load Kraken2 reports
    # ------------------------------------------------------------------
    std_dir = Path(args.std_reports)
    if not std_dir.is_dir():
        sys.exit(f"ERROR: --std-reports directory not found: {std_dir}")

    reports = sorted(std_dir.glob("*.report"))
    if not reports:
        sys.exit(f"ERROR: No *.report files found in {std_dir}")

    dfs: list[pl.DataFrame] = []
    for rpt in reports:
        dfs.append(read_std_report(rpt, rpt.stem))

    df = pl.concat(dfs, how="vertical")
    print(
        f"[INFO] Loaded {len(reports)} reports, {df.height} rows total.",
        file=sys.stderr,
    )

    # ------------------------------------------------------------------
    # Remove unwanted samples
    # ------------------------------------------------------------------
    if args.samples_to_remove:
        bad = {
            line.strip()
            for line in Path(args.samples_to_remove).read_text().splitlines()
            if line.strip()
        }
        if bad:
            df = df.filter(~pl.col("sample").is_in(list(bad)))
            print(f"[INFO] Removed samples: {bad}", file=sys.stderr)

    # ------------------------------------------------------------------
    # Join DB info
    # ------------------------------------------------------------------
    df = df.join(inspect_df, on="taxid", how="left")

    # ------------------------------------------------------------------
    # Species-level filtering
    # ------------------------------------------------------------------
    print("[INFO] Filtering to species-level taxa.", file=sys.stderr)
    df = filter_to_best_kraken_taxa(df, args.virus_direct_ratio)

    # ------------------------------------------------------------------
    # Host marking
    # ------------------------------------------------------------------
    print(f"[INFO] Marking host taxa (organism={args.organism!r}).", file=sys.stderr)
    df = mark_host_taxa(df, args.organism, reference_path)

    # ------------------------------------------------------------------
    # Minimizer proportions (non-host denominator)
    # ------------------------------------------------------------------
    non_host = df.filter(~pl.col("is_host"))

    tot = non_host.group_by("sample").agg(
        pl.col("distinct_minimizers").sum().alias("total_distinct_minimizers_non_host")
    )
    non_host = non_host.join(tot, on="sample", how="left").with_columns(
        (
            pl.col("distinct_minimizers") / pl.col("total_distinct_minimizers_non_host")
        ).alias("minimizer_proportion_sample")
    )

    # ------------------------------------------------------------------
    # genome_coverage_estimate
    # ------------------------------------------------------------------
    if "db_taxon_minimizers" in non_host.columns:
        non_host = non_host.with_columns(
            pl.when(pl.col("db_taxon_minimizers").cast(pl.Float64).fill_null(0.0) > 0)
            .then(
                pl.col("distinct_minimizers").cast(pl.Float64)
                / pl.col("db_taxon_minimizers").cast(pl.Float64)
            )
            .otherwise(None)
            .alias("genome_coverage_estimate")
        )
    else:
        non_host = non_host.with_columns(pl.lit(None).cast(pl.Float64).alias("genome_coverage_estimate"))

    # ------------------------------------------------------------------
    # min-prop filter
    # ------------------------------------------------------------------
    non_host = non_host.filter(pl.col("minimizer_proportion_sample") >= args.min_prop)

    # ------------------------------------------------------------------
    # Taxon descriptive stats (always computed)
    # ------------------------------------------------------------------
    stats = non_host.group_by("taxid").agg(
        pl.col("minimizer_proportion_sample").mean().alias("taxon_mean_prop"),
        pl.col("minimizer_proportion_sample").std().alias("taxon_sd_prop"),
        pl.col("sample").n_unique().alias("n_samples_with_taxon"),
    )
    non_host = non_host.join(stats, on="taxid", how="left")

    # ------------------------------------------------------------------
    # Statistics
    # ------------------------------------------------------------------

    if use_negative:
        # Negative-control–based stats
        neg = non_host.filter(pl.col("sample") == args.negative_control)
        if neg.is_empty():
            print(
                f"[WARN] Negative control '{args.negative_control}' not found.",
                file=sys.stderr,
            )

        neg_bg = neg.group_by("taxid").agg(
            pl.col("minimizer_proportion_sample").mean().alias("neg_prop")
        )
        non_host = (
            non_host.join(neg_bg, on="taxid", how="left")
            .with_columns(pl.col("neg_prop").fill_null(0.0))
            .with_columns(
                (pl.col("neg_prop") + args.pseudocount).alias("mu")
            )
            .with_columns(
                ((pl.col("minimizer_proportion_sample") - pl.col("mu")) / pl.col("mu").sqrt()).alias("z_score")
            )
            .with_columns(
                pl.col("z_score")
                .map_elements(
                    lambda z: None
                    if z is None or (isinstance(z, float) and math.isnan(z))
                    else 1.0 - float(norm.cdf(z)),
                    return_dtype=pl.Float64,
                )
                .alias("p_value"),
                pl.lit("negative").alias("stat_mode"),
            )
        )

    else:
        # Cohort-based stats (n >= 3)
        non_host = (
            non_host.with_columns(
                (
                    (pl.col("n_samples_with_taxon") >= 3)
                    & pl.col("taxon_sd_prop").is_not_null()
                    & (pl.col("taxon_sd_prop") > 0)
                ).alias("has_valid_stats")
            )
            .with_columns(
                pl.when(pl.col("has_valid_stats"))
                .then(
                    (pl.col("minimizer_proportion_sample") - pl.col("taxon_mean_prop"))
                    / pl.col("taxon_sd_prop")
                )
                .otherwise(None)
                .alias("z_score")
            )
            .with_columns(
                pl.col("z_score")
                .map_elements(
                    lambda z: None
                    if z is None or (isinstance(z, float) and math.isnan(z))
                    else 1.0 - float(norm.cdf(z)),
                    return_dtype=pl.Float64,
                )
                .alias("p_value")
            )
            .with_columns(
                pl.when(pl.col("has_valid_stats"))
                .then(pl.lit("cohort"))
                .otherwise(pl.lit("none"))
                .alias("stat_mode")
            )
        )

    # ------------------------------------------------------------------
    # Benjamini–Hochberg FDR
    # ------------------------------------------------------------------
    non_host = non_host.with_columns(
        benjamini_hochberg(non_host["p_value"]).alias("fdr")
    )

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------
    out = outdir / f"{args.prefix}.tsv"
    non_host.sort(
        by=["sample", "minimizer_proportion_sample"],
        descending=[False, True],
    ).write_csv(out, separator="\t")

    print(f"[INFO] Wrote {non_host.height} rows to: {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
