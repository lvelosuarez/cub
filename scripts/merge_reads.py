#!/usr/bin/env python3
"""
merge_reads.py — Clinical Metagenomics Multi-source Read Merger
===============================================================

PURPOSE
-------
This script merges per-read evidence from three independent sources:

  1. Kraken2 PlusPF classification  (.kraken2 per-read output)
  2. Kraken2 EuPathDB classification (.kraken2 per-read output)
  3. Bowtie2 BAM mapping against a curated pathogen index (~1100–1500 species)

It is designed for low-biomass clinical samples where contamination is a
major concern. Rather than relying on a single tool, it integrates:

  - Per-read kmer evidence quality (hit fraction, distinct ratio, hit groups)
  - Mapping quality and alignment metrics from the BAM
  - Coverage breadth across the reference genome (not just depth)
  - Human kmer conflict detection (reads that map to bacteria but have
    strong human kmer signal in PlusPF — common in conserved gene regions)

OUTPUT
------
A TSV file with one row per read, containing all merged evidence columns
plus a 'read_call' column classifying each read into one of:

  KRAKEN_CONFIRMED       BAM maps to pathogen, Kraken2 agrees (non-human taxid),
                         all quality metrics pass
  BAM_ONLY               BAM maps to pathogen, Kraken2 did not classify this read
                         (could be a region not well covered by PlusPF)
  HUMAN_CONFLICT         Kraken2 classifies the read as human with strong signal
                         (kmer_hit_frac >= threshold) AND breadth is low —
                         likely a human read hitting a conserved bacterial locus
  HUMAN_KMER_BROAD_COVERAGE
                         Kraken2 says human BUT the accession has broad genome
                         coverage — organism is likely truly present, the human
                         kmer signal is incidental (conserved domain)
  LOW_CONF_HUMAN_KMER    Kraken2 gives a weak human signal — borderline case,
                         kept but flagged for cohort-level review
  KRAKEN_LOW_QUALITY     Kraken2 classified the read but quality metrics are poor
                         (low hit fraction, low distinct ratio, few hit groups)

FILTERING LOGIC
---------------
BAM-side filters (applied before merge):
  --mapq-min        Minimum mapping quality (default: 30)
                    Reads below this threshold are discarded UNLESS the accession
                    has broad coverage (breadth rescue — see below)
  --mapq-floor      Minimum MAPQ for breadth-rescued reads (default: 10)
                    Reads with MAPQ between mapq-floor and mapq-min are kept
                    if the accession has enough breadth evidence
  --nm-rate-max     Maximum mismatch rate NM/read_length (default: 0.05 = 5%)
  --alen-min        Minimum aligned length in bp, excluding soft/hard clips
                    (default: 45 bp, appropriate for 75 bp SE reads)

Breadth parameters:
  --window-size     Size of genomic windows for breadth computation (default: 500 bp)
                    The genome is divided into non-overlapping windows; a window
                    is 'covered' if at least one read maps to it
  --breadth-min-windows   Minimum number of covered windows for an accession to
                          be considered 'broadly covered' (default: 3)
  --breadth-min-fraction  Minimum fraction of the genome covered (default: 0.01 = 1%)
                          C. acnes in a typical skin sample reaches ~2.6%

Kraken2 quality filters (applied post-merge):
  --kmer-hit-frac-min    Minimum fraction of read kmers that hit an informative
                         taxid (default: 0.15). Equivalent to --confidence in
                         Kraken2 runtime, but computed post-hoc so you do not
                         need to re-run Kraken2
  --distinct-ratio-min   Minimum fraction of classified kmers owned by the top
                         taxid (default: 0.50). Low value = kmers spread across
                         many taxa = likely conserved region noise
  --hit-groups-min       Minimum number of non-contiguous kmer runs on the read
                         (default: 3). A single conserved domain hit produces
                         1 hit group; a real organism read produces many
  --n-taxa-max           Maximum number of distinct informative taxids hit by
                         the read (default: 5). More than this = ambiguous

Human conflict thresholds:
  --human-khf-cutoff     kmer_hit_frac above which a human Kraken2 call is
                         considered a hard conflict (default: 0.50)

TYPICAL USAGE
-------------
Single sample:
  python merge_reads.py \\
      --k1  sample.pluspf.krk \\
      --k2  sample.eupath.krk \\
      --bam sample.bam \\
      --taxmap bartlett_kraken2_index.tsv \\
      --sample-id SAMPLE_001 \\
      --out sample_001_merged.tsv

Override thresholds:
  python merge_reads.py \\
      --k1 sample.pluspf.krk --k2 sample.eupath.krk \\
      --bam sample.bam --taxmap bartlett_kraken2_index.tsv \\
      --sample-id SAMPLE_001 --out sample_001_merged.tsv \\
      --mapq-min 20 --breadth-min-fraction 0.005 --hit-groups-min 2

NOTE ON DEDUPLICATION
---------------------
This script assumes reads have already been deduplicated BEFORE dehosting
(the recommended order for low-biomass clinical samples). The BAM dedup
step here only removes reads with identical read_id that mapped to multiple
references — it keeps the highest MAPQ hit per read_id.

NOTE ON BREADTH COMPUTATION
----------------------------
Breadth is computed on ALL primary mapped reads (before MAPQ/NM filtering).
This is intentional: breadth answers 'is there any coverage spread across
the genome?' which is a genome-level question independent of individual
read quality. MAPQ filtering controls which reads enter the evidence merge.
"""

import argparse
import logging
import re
import sys
from pathlib import Path

import polars as pl
import pysam


# ── logging setup ─────────────────────────────────────────────────────────────

def setup_logging(sample_id: str, log_file: str | None = None) -> logging.Logger:
    logger = logging.getLogger("merge_reads")
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    # screen handler
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    # file handler
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


# ── constants ──────────────────────────────────────────────────────────────────

UNINFORMATIVE_TAXIDS = frozenset({0, 1, 131567})
MATE_SUFFIX_RE       = re.compile(r"/[12]$")


# ── kmer string parser ────────────────────────────────────────────────────────

def parse_kmer_string(kmer_str: str) -> tuple[float, float, int, int]:
    """
    Parse the raw Kraken2 kmer hit string (column 5 of .kraken2 output).

    Returns
    -------
    kmer_hit_frac  : classified_kmers / total_kmers
                     proxy for Kraken2 --confidence threshold
    distinct_ratio : top_taxid_kmers / classified_kmers
                     1.0 = one taxid owns all classified kmers (clean)
                     low = kmers spread across many taxa (ambiguous)
    n_taxa_hit     : number of distinct informative taxids hit
    hit_groups     : number of non-contiguous runs of informative kmers
                     1 = single conserved domain hit
                     ≥3 = reads distributed across multiple loci (genuine)
    """
    if not kmer_str or kmer_str.strip() == "":
        return 0.0, 0.0, 0, 0

    total        = 0
    classified   = 0
    taxa_counts: dict[int, int] = {}
    hit_groups   = 0
    in_hit_run   = False

    for token in kmer_str.split():
        if ":" not in token:
            continue
        taxid_s, count_s = token.split(":", 1)
        try:
            taxid = int(taxid_s)
            count = int(count_s)
        except ValueError:
            continue

        total += count
        is_informative = taxid not in UNINFORMATIVE_TAXIDS

        if is_informative:
            classified += count
            taxa_counts[taxid] = taxa_counts.get(taxid, 0) + count
            if not in_hit_run:
                hit_groups += 1
                in_hit_run = True
        else:
            in_hit_run = False

    kmer_hit_frac  = classified / total if total > 0 else 0.0
    distinct_ratio = (
        max(taxa_counts.values()) / classified if classified > 0 else 0.0
    )
    return kmer_hit_frac, distinct_ratio, len(taxa_counts), hit_groups


# ── kraken2 loader ────────────────────────────────────────────────────────────

def load_kraken2_per_read(path: str, db_tag: str) -> pl.DataFrame:
    """Load a Kraken2 .kraken2 per-read file into a Polars DataFrame."""
    rows = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            status   = parts[0]
            rid      = MATE_SUFFIX_RE.sub("", parts[1])
            taxid_s  = parts[2]
            length_s = parts[3]
            kmer_s   = parts[4]

            try:
                taxid = int(taxid_s)
            except ValueError:
                taxid = 0
            try:
                read_len = int(length_s)
            except ValueError:
                read_len = None

            khf, dr, n_taxa, hg = parse_kmer_string(kmer_s)
            rows.append({
                "read_id":        rid,
                "status":         status,
                "taxid":          taxid,
                "read_len":       read_len,
                "kmer_hit_frac":  khf,
                "distinct_ratio": dr,
                "n_taxa_hit":     n_taxa,
                "hit_groups":     hg,
                "db":             db_tag,
            })

    return pl.DataFrame(rows, schema={
        "read_id":        pl.Utf8,
        "status":         pl.Utf8,
        "taxid":          pl.Int64,
        "read_len":       pl.Int32,
        "kmer_hit_frac":  pl.Float64,
        "distinct_ratio": pl.Float64,
        "n_taxa_hit":     pl.Int32,
        "hit_groups":     pl.Int32,
        "db":             pl.Utf8,
    })


def prefix_kr_cols(df: pl.DataFrame, tag: str) -> pl.DataFrame:
    """Rename all columns except read_id with a db-specific prefix."""
    return df.rename({
        c: f"{tag}_{c}" for c in df.columns if c != "read_id"
    })


# ── BAM loader + breadth ──────────────────────────────────────────────────────

def load_bam(
    path: str,
    window_size: int,
    min_depth: int,
    logger: logging.Logger,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    Single-pass BAM scan.
    Returns (per_read_df, breadth_df).
    Breadth is computed on ALL primary reads before quality filtering.
    """
    rows:        list[dict]           = []
    window_hits: dict[str, list[int]] = {}
    ref_lengths: dict[str, int]       = {}
    n_secondary = n_supplementary = n_unmapped = 0

    with pysam.AlignmentFile(path, "rb") as bam:
        for ref, length in zip(bam.references, bam.lengths):
            ref_lengths[ref] = length
            window_hits[ref] = [0] * max(1, length // window_size)

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                n_unmapped += 1
                continue
            if read.is_secondary:
                n_secondary += 1
                continue
            if read.is_supplementary:
                n_supplementary += 1
                continue

            rid  = MATE_SUFFIX_RE.sub("", read.query_name)
            rnam = read.reference_name

            rows.append({
                "read_id":   rid,
                "bam_rname": rnam,
                "bam_mapq":  read.mapping_quality,
                "bam_nm":    read.get_tag("NM") if read.has_tag("NM") else None,
                "bam_rlen":  read.query_length,
                "bam_alen":  read.query_alignment_length,
                "bam_flag":  read.flag,
            })

            n_win = len(window_hits[rnam])
            w     = min(read.reference_start // window_size, n_win - 1)
            window_hits[rnam][w] += 1

    logger.info(
        "BAM raw counts | primary_mapped=%d | unmapped=%d | "
        "secondary=%d | supplementary=%d",
        len(rows), n_unmapped, n_secondary, n_supplementary,
    )

    # breadth DataFrame
    breadth_rows = []
    for ref, hits in window_hits.items():
        ref_len   = ref_lengths[ref]
        n_total   = len(hits)
        n_covered = sum(1 for h in hits if h >= min_depth)
        breadth_rows.append({
            "bam_rname":         ref,
            "ref_length":        ref_len,
            "n_windows_total":   n_total,
            "n_windows_covered": n_covered,
            "breadth_fraction":  (n_covered * window_size) / ref_len
                                 if ref_len > 0 else 0.0,
        })

    breadth_df = pl.DataFrame(breadth_rows, schema={
        "bam_rname":         pl.Utf8,
        "ref_length":        pl.Int64,
        "n_windows_total":   pl.Int32,
        "n_windows_covered": pl.Int32,
        "breadth_fraction":  pl.Float64,
    })

    top10 = (
        breadth_df
        .sort("n_windows_covered", descending=True)
        .head(10)
    )
    logger.info("Top 10 accessions by breadth:\n%s", top10)

    # per-read DataFrame
    bam_raw = pl.DataFrame(rows, schema={
        "read_id":   pl.Utf8,
        "bam_rname": pl.Utf8,
        "bam_mapq":  pl.Int32,
        "bam_nm":    pl.Int32,
        "bam_rlen":  pl.Int32,
        "bam_alen":  pl.Int32,
        "bam_flag":  pl.Int32,
    }).with_columns(
        pl.when(pl.col("bam_rlen") > 0)
          .then(pl.col("bam_nm") / pl.col("bam_rlen"))
          .otherwise(None)
          .alias("bam_nm_rate")
    )

    return bam_raw, breadth_df


# ── BAM filter ────────────────────────────────────────────────────────────────

def filter_bam(
    bam_raw: pl.DataFrame,
    breadth_df: pl.DataFrame,
    mapq_min: int,
    mapq_floor: int,
    nm_rate_max: float,
    alen_min: int,
    breadth_min_windows: int,
    breadth_min_fraction: float,
    logger: logging.Logger,
) -> pl.DataFrame:
    """Apply breadth-conditional MAPQ/NM/alen filter and dedup."""
    # join breadth before filtering so we can use it in the condition
    bam = bam_raw.join(breadth_df, on="bam_rname", how="left")
    n_before = len(bam)

    high_confidence = (
        (pl.col("bam_mapq")    >= mapq_min)    &
        (pl.col("bam_nm_rate") <= nm_rate_max) &
        (pl.col("bam_alen")    >= alen_min)
    )
    low_mapq_broad = (
        (pl.col("bam_mapq")          >= mapq_floor)           &
        (pl.col("bam_mapq")          <  mapq_min)             &
        (pl.col("bam_nm_rate")       <= nm_rate_max)          &
        (pl.col("bam_alen")          >= alen_min)             &
        (pl.col("n_windows_covered") >= breadth_min_windows)  &
        (pl.col("breadth_fraction")  >= breadth_min_fraction)
    )

    bam = bam.with_columns(
        pl.when(high_confidence).then(pl.lit("HIGH_CONF"))
          .when(low_mapq_broad).then(pl.lit("BREADTH_RESCUED"))
          .otherwise(pl.lit("FAILED"))
          .alias("mapq_filter_reason")
    )

    n_hc      = (bam["mapq_filter_reason"] == "HIGH_CONF").sum()
    n_rescued = (bam["mapq_filter_reason"] == "BREADTH_RESCUED").sum()
    n_failed  = (bam["mapq_filter_reason"] == "FAILED").sum()

    logger.info(
        "BAM pre-filter | mapq>=%d OR [mapq>=%d + breadth>=%.3f] | "
        "high_conf=%d | breadth_rescued=%d | failed=%d | kept=%d",
        mapq_min, mapq_floor, breadth_min_fraction,
        n_hc, n_rescued, n_failed, n_hc + n_rescued,
    )

    bam = bam.filter(pl.col("mapq_filter_reason") != "FAILED")

    # dedup: keep highest MAPQ per read_id
    n_before_dedup = len(bam)
    bam = (
        bam
        .sort("bam_mapq", descending=True)
        .group_by("read_id")
        .first()
    )
    n_removed = n_before_dedup - len(bam)
    if n_removed > 0:
        logger.info("BAM dedup | removed=%d duplicate read_ids", n_removed)

    return bam


# ── taxmap join ───────────────────────────────────────────────────────────────

def join_taxmap(bam: pl.DataFrame, taxmap_path: str, logger: logging.Logger) -> pl.DataFrame:
    taxmap = pl.read_csv(taxmap_path, separator="\t")

    dupes = taxmap.group_by("accession").len().filter(pl.col("len") > 1)
    if len(dupes) > 0:
        logger.warning(
            "taxmap has %d duplicated accessions — keeping first occurrence", len(dupes)
        )
        taxmap = taxmap.unique(subset=["accession"], keep="first")

    taxmap_clean = (
        taxmap
        .rename({"accession": "bam_rname"})
        .with_columns([
            pl.col("ncbi_taxid_strain").cast(pl.Int64).alias("bam_taxid_accession"),
            pl.col("ncbi_taxid_species").cast(pl.Int64).alias("bam_taxid_species"),
            pl.coalesce([
                pl.when(
                    pl.col("species_name_kraken2").is_not_null() &
                    (pl.col("species_name_kraken2") != "")
                ).then(pl.col("species_name_kraken2")),
                pl.when(
                    pl.col("bartlett_species").is_not_null() &
                    (pl.col("bartlett_species") != "")
                ).then(pl.col("bartlett_species")),
                pl.col("gtdb_species"),
            ]).alias("taxon_name"),
        ])
        .select(["bam_rname", "bam_taxid_accession", "bam_taxid_species", "taxon_name"])
    )

    result = bam.join(taxmap_clean, on="bam_rname", how="left")

    no_taxmap = result.filter(pl.col("bam_taxid_species").is_null())
    if len(no_taxmap) > 0:
        top10 = (
            no_taxmap["bam_rname"]
            .value_counts()
            .sort("count", descending=True)
            .head(10)
        )
        logger.warning(
            "%d reads have no taxmap entry. Top unknown accessions:\n%s",
            len(no_taxmap), top10,
        )
    return result


# ── merge ─────────────────────────────────────────────────────────────────────

def merge_all(
    bam_df: pl.DataFrame,
    k1: pl.DataFrame,
    k2: pl.DataFrame,
    kmer_hit_frac_min: float,
    distinct_ratio_min: float,
    hit_groups_min: int,
    n_taxa_max: int,
    logger: logging.Logger,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Merge BAM + Kraken2 frames and apply Kraken quality filter."""
    k1p = prefix_kr_cols(k1, "k1")
    k2p = prefix_kr_cols(k2, "k2")

    merged = (
        bam_df
        .join(k1p, on="read_id", how="left")
        .join(k2p, on="read_id", how="left")
    )
    logger.info("Merged shape (pre-Kraken-filter): %s", merged.shape)

    def kr_pass(tag: str) -> pl.Expr:
        s   = pl.col(f"{tag}_status")
        khf = pl.col(f"{tag}_kmer_hit_frac")
        dr  = pl.col(f"{tag}_distinct_ratio")
        hg  = pl.col(f"{tag}_hit_groups")
        nt  = pl.col(f"{tag}_n_taxa_hit")
        return (
            s.is_null() | (s == "U") |
            (
                (s == "C") &
                (khf >= kmer_hit_frac_min) &
                (dr  >= distinct_ratio_min) &
                (hg  >= hit_groups_min) &
                (nt  <= n_taxa_max)
            )
        )

    passed  = merged.filter(kr_pass("k1") & kr_pass("k2"))
    rejected = merged.filter(
        ~(kr_pass("k1") & kr_pass("k2"))
    ).with_columns(pl.lit("kraken_prefilter").alias("filter_reason"))

    logger.info(
        "Kraken post-merge filter | kmer_hit_frac>=%.2f | distinct_ratio>=%.2f | "
        "hit_groups>=%d | n_taxa_max<=%d | before=%d | passed=%d | rejected=%d",
        kmer_hit_frac_min, distinct_ratio_min, hit_groups_min, n_taxa_max,
        len(merged), len(passed), len(rejected),
    )
    return passed, rejected


# ── read classifier ───────────────────────────────────────────────────────────

def classify_reads(
    df: pl.DataFrame,
    human_taxids: frozenset,
    human_khf_cutoff: float,
    kmer_hit_frac_min: float,
    distinct_ratio_min: float,
    hit_groups_min: int,
    n_taxa_max: int,
    breadth_min_windows: int,
    breadth_min_fraction: float,
) -> pl.DataFrame:
    """Assign read_call label to every read."""
    has_breadth = (
        (pl.col("n_windows_covered") >= breadth_min_windows) &
        (pl.col("breadth_fraction")  >= breadth_min_fraction)
    )
    is_human      = pl.col("k1_taxid").is_in(list(human_taxids))
    k1_classified = pl.col("k1_status") == "C"
    k1_quality = (
        (pl.col("k1_kmer_hit_frac")  >= kmer_hit_frac_min)  &
        (pl.col("k1_distinct_ratio") >= distinct_ratio_min) &
        (pl.col("k1_hit_groups")     >= hit_groups_min)     &
        (pl.col("k1_n_taxa_hit")     <= n_taxa_max)
    )

    return df.with_columns(
        pl.when(
            k1_classified & is_human &
            (pl.col("k1_kmer_hit_frac") >= human_khf_cutoff) &
            ~has_breadth
        ).then(pl.lit("HUMAN_CONFLICT"))

        .when(
            k1_classified & is_human &
            (pl.col("k1_kmer_hit_frac") >= human_khf_cutoff) &
            has_breadth
        ).then(pl.lit("HUMAN_KMER_BROAD_COVERAGE"))

        .when(
            k1_classified & is_human &
            (pl.col("k1_kmer_hit_frac") < human_khf_cutoff)
        ).then(pl.lit("LOW_CONF_HUMAN_KMER"))

        .when(
            k1_classified & ~is_human & k1_quality
        ).then(pl.lit("KRAKEN_CONFIRMED"))

        .when(
            pl.col("k1_status").is_null() | (pl.col("k1_status") == "U")
        ).then(pl.lit("BAM_ONLY"))

        .otherwise(pl.lit("KRAKEN_LOW_QUALITY"))
        .alias("read_call")
    )


# ── argument parser ───────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="merge_reads.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # required inputs
    req = p.add_argument_group("required inputs")
    req.add_argument("--k1",      required=True, metavar="FILE",
                     help="Kraken2 PlusPF per-read output (.kraken2 format)")
    req.add_argument("--k2",      required=True, metavar="FILE",
                     help="Kraken2 EuPathDB per-read output (.kraken2 format)")
    req.add_argument("--bam",     required=True, metavar="FILE",
                     help="Bowtie2 BAM aligned against pathogen index")
    req.add_argument("--taxmap",  required=True, metavar="FILE",
                     help="Taxmap TSV (accession → taxid + species name). "
                          "Columns: accession, ncbi_taxid_strain, ncbi_taxid_species, "
                          "species_name_kraken2, bartlett_species, gtdb_species")
    req.add_argument("--sample-id", required=True, metavar="STR",
                     help="Sample identifier added to every output row")

    # output
    out = p.add_argument_group("output")
    out.add_argument("--out", required=True, metavar="FILE",
                     help="Output TSV path (per-read merged table)")
    out.add_argument("--log", default=None, metavar="FILE",
                     help="Log file path (default: stdout only)")
    out.add_argument("--rejected-out", default=None, metavar="FILE",
                     help="Optional TSV for reads rejected by Kraken quality filter")

    # BAM filter parameters
    bam_grp = p.add_argument_group("BAM filter parameters")
    bam_grp.add_argument("--mapq-min", type=int, default=30,
                         help="Minimum MAPQ for high-confidence reads (default: 30)")
    bam_grp.add_argument("--mapq-floor", type=int, default=10,
                         help="Minimum MAPQ for breadth-rescued reads (default: 10)")
    bam_grp.add_argument("--nm-rate-max", type=float, default=0.05,
                         help="Maximum mismatch rate NM/read_length (default: 0.05)")
    bam_grp.add_argument("--alen-min", type=int, default=45,
                         help="Minimum aligned length in bp (default: 45)")

    # breadth parameters
    br_grp = p.add_argument_group("breadth parameters")
    br_grp.add_argument("--window-size", type=int, default=500,
                        help="Genomic window size in bp for breadth computation (default: 500)")
    br_grp.add_argument("--min-depth", type=int, default=1,
                        help="Minimum read depth for a window to be counted as covered "
                             "(default: 1)")
    br_grp.add_argument("--breadth-min-windows", type=int, default=3,
                        help="Minimum covered windows for breadth rescue (default: 3)")
    br_grp.add_argument("--breadth-min-fraction", type=float, default=0.01,
                        help="Minimum genome fraction covered for breadth rescue "
                             "(default: 0.01)")

    # Kraken2 quality parameters
    kr_grp = p.add_argument_group("Kraken2 quality filter parameters")
    kr_grp.add_argument("--kmer-hit-frac-min", type=float, default=0.15,
                        help="Minimum kmer_hit_frac (default: 0.15)")
    kr_grp.add_argument("--distinct-ratio-min", type=float, default=0.50,
                        help="Minimum distinct_ratio (default: 0.50)")
    kr_grp.add_argument("--hit-groups-min", type=int, default=3,
                        help="Minimum hit_groups (default: 3)")
    kr_grp.add_argument("--n-taxa-max", type=int, default=5,
                        help="Maximum n_taxa_hit (default: 5)")

    # human conflict parameters
    hum_grp = p.add_argument_group("human conflict parameters")
    hum_grp.add_argument("--human-khf-cutoff", type=float, default=0.50,
                         help="kmer_hit_frac threshold above which a human Kraken2 call "
                              "is considered a hard conflict (default: 0.50)")

    return p


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = build_parser()
    args   = parser.parse_args()

    logger = setup_logging(args.sample_id, args.log)
    logger.info("=" * 70)
    logger.info("merge_reads.py — sample: %s", args.sample_id)
    logger.info("=" * 70)
    logger.info("Inputs:")
    logger.info("  k1 (PlusPF)  : %s", args.k1)
    logger.info("  k2 (EuPath)  : %s", args.k2)
    logger.info("  bam          : %s", args.bam)
    logger.info("  taxmap       : %s", args.taxmap)
    logger.info("Parameters:")
    logger.info("  mapq_min=%d  mapq_floor=%d  nm_rate_max=%.3f  alen_min=%d",
                args.mapq_min, args.mapq_floor, args.nm_rate_max, args.alen_min)
    logger.info("  window_size=%d  breadth_min_windows=%d  breadth_min_fraction=%.4f",
                args.window_size, args.breadth_min_windows, args.breadth_min_fraction)
    logger.info("  kmer_hit_frac_min=%.2f  distinct_ratio_min=%.2f  "
                "hit_groups_min=%d  n_taxa_max=%d",
                args.kmer_hit_frac_min, args.distinct_ratio_min,
                args.hit_groups_min, args.n_taxa_max)
    logger.info("  human_khf_cutoff=%.2f", args.human_khf_cutoff)

    # ── load Kraken2 ──────────────────────────────────────────────────────────
    logger.info("Loading Kraken2 PlusPF: %s", args.k1)
    k1 = load_kraken2_per_read(args.k1, db_tag="pluspf")
    logger.info("  k1 rows: %d", len(k1))

    logger.info("Loading Kraken2 EuPathDB: %s", args.k2)
    k2 = load_kraken2_per_read(args.k2, db_tag="eupath")
    logger.info("  k2 rows: %d", len(k2))

    # ── load BAM + breadth ────────────────────────────────────────────────────
    logger.info("Loading BAM: %s", args.bam)
    bam_raw, breadth_df = load_bam(
        args.bam,
        window_size=args.window_size,
        min_depth=args.min_depth,
        logger=logger,
    )

    # ── filter BAM ────────────────────────────────────────────────────────────
    bam_filtered = filter_bam(
        bam_raw, breadth_df,
        mapq_min=args.mapq_min,
        mapq_floor=args.mapq_floor,
        nm_rate_max=args.nm_rate_max,
        alen_min=args.alen_min,
        breadth_min_windows=args.breadth_min_windows,
        breadth_min_fraction=args.breadth_min_fraction,
        logger=logger,
    )

    # ── join taxmap ───────────────────────────────────────────────────────────
    logger.info("Joining taxmap: %s", args.taxmap)
    bam_df = join_taxmap(bam_filtered, args.taxmap, logger)

    # ── merge ─────────────────────────────────────────────────────────────────
    logger.info("Merging BAM + Kraken2 frames")
    passed, rejected = merge_all(
        bam_df, k1, k2,
        kmer_hit_frac_min=args.kmer_hit_frac_min,
        distinct_ratio_min=args.distinct_ratio_min,
        hit_groups_min=args.hit_groups_min,
        n_taxa_max=args.n_taxa_max,
        logger=logger,
    )

    # ── classify reads ────────────────────────────────────────────────────────
    logger.info("Classifying reads")
    human_taxids = frozenset({9606, 9605, 207598, 9604})
    passed = classify_reads(
        passed,
        human_taxids=human_taxids,
        human_khf_cutoff=args.human_khf_cutoff,
        kmer_hit_frac_min=args.kmer_hit_frac_min,
        distinct_ratio_min=args.distinct_ratio_min,
        hit_groups_min=args.hit_groups_min,
        n_taxa_max=args.n_taxa_max,
        breadth_min_windows=args.breadth_min_windows,
        breadth_min_fraction=args.breadth_min_fraction,
    )

    # add sample_id to output
    passed = passed.with_columns(pl.lit(args.sample_id).alias("sample_id"))

    # ── summary log ───────────────────────────────────────────────────────────
    summary = (
        passed
        .group_by("read_call")
        .agg([
            pl.len().alias("n_reads"),
            pl.col("taxon_name").n_unique().alias("n_taxa"),
            pl.col("k1_kmer_hit_frac").mean().alias("mean_khf"),
            pl.col("k1_hit_groups").mean().alias("mean_hg"),
            pl.col("n_windows_covered").mean().alias("mean_windows"),
            pl.col("breadth_fraction").mean().alias("mean_breadth"),
        ])
        .sort("n_reads", descending=True)
    )
    logger.info("Read call summary:\n%s", summary)

    taxon_summary = (
        passed
        .group_by(["taxon_name", "read_call"])
        .agg([
            pl.len().alias("n_reads"),
            pl.col("n_windows_covered").first(),
            pl.col("breadth_fraction").first(),
            pl.col("k1_kmer_hit_frac").mean().alias("mean_khf"),
            pl.col("bam_mapq").mean().alias("mean_mapq"),
        ])
        .sort(["n_reads"], descending=True)
    )
    logger.info("Per-taxon summary:\n%s", taxon_summary)

    # ── write outputs ─────────────────────────────────────────────────────────
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    passed.write_csv(out_path, separator="\t")
    logger.info("Output written: %s  (%d rows, %d columns)",
                out_path, len(passed), len(passed.columns))

    if args.rejected_out:
        rej_path = Path(args.rejected_out)
        rej_path.parent.mkdir(parents=True, exist_ok=True)
        rejected.write_csv(rej_path, separator="\t")
        logger.info("Rejected reads written: %s  (%d rows)", rej_path, len(rejected))

    logger.info("Done.")


if __name__ == "__main__":
    main()
