#!/usr/bin/env python3
"""
analyse.py — CUB post-analysis pipeline (extracted from folder.py marimo notebook).

Runs evidence integration, RPM normalisation, and the full decontamination decision
chain for one sequencing run, then writes:
  output/analysis/cohort_norm.parquet   — patient rows with decontam_decision
  output/analysis/tpos_qc.parquet       — one row per TPOS sample (QC metrics)
  output/analysis/tpos_report.html      — self-contained HTML QC report

Usage (called by Snakemake rule 'analyse'):
    python /mnt/san/microbio/apps/cub/scripts/analyse.py [--run-id N]

The working directory must be the run folder (set by snakemake --directory).
"""
import argparse
import math
import os
import re
from pathlib import Path

import polars as pl
import pysam

from analyse_report import write_tpos_html

# ── Constants ─────────────────────────────────────────────────────────────────
UNINFORMATIVE_TAXIDS  = frozenset({0, 1, 131567, 28384, 2759})
HUMAN_TAXIDS          = frozenset({9606, 9605, 207598, 9604})
HUMAN_KHF_HARD_CUTOFF = 0.50
KMER_HIT_FRAC_MIN     = 0.15
DISTINCT_RATIO_MIN    = 0.50
HIT_GROUPS_MIN        = 3
N_TAXA_MAX            = 5
BREADTH_MIN_WINDOWS   = 3
BREADTH_MIN_FRACTION  = 0.01
MAPQ_MIN              = 30
MAPQ_FLOOR            = 10
NM_RATE_MAX           = 0.05
ALEN_MIN              = 45
WINDOW_SIZE           = 500
MIN_DEPTH             = 1
MATE_SUFFIX_RE        = re.compile(r"/[12]$")

REPORTABLE_RANKS = {"species", "subspecies", "strain", "genus", "no rank"}
NONREPORTABLE_RANKS = {
    "phylum", "class", "order", "family",
    "superphylum", "subphylum", "subclass", "suborder", "superfamily",
    "kingdom", "superkingdom",
}
GENUS_LEVEL_COLLAPSE = {
    3048414: 687331, 3048415: 687331, 3048417: 687331, 3048418: 687331,
    3048419: 687331, 3048420: 687331, 3048424: 687331, 3048427: 687331,
    3048432: 687331, 3048433: 687331, 3050298: 687331,
    3052014: 687332,
    1982268: 1982251, 1982262: 1982251, 1982269: 1982251,
    1982283: 1982251, 1982257: 1982251, 1654781: 1982251,
}

P_KRAKEN2_NODES  = "/mnt/san/microbio/database/kraken2/PlusPF20260226/nodes.dmp"
P_KRAKEN2_NAMES  = "/mnt/san/microbio/database/kraken2/PlusPF20260226/names.dmp"
P_EUPATH_INSPECT = "/mnt/san/microbio/database/kraken2/eupathdb48_20230407/inspect.txt"
P_TAXMAP         = "/mnt/san/microbio/apps/cub/resources/bartlett_kraken2_index.tsv"
BASE             = "/mnt/san/microbio/metagenomique_clinique"

P_CONTAMINANT_DB    = "/mnt/san/microbio/apps/cub/resources/decontam/contaminant_db.parquet"
EBV_TAXID           = 10376
MIN_BLANKS_N        = 2
FC_TIER1_BREADTH    = 5.0
FC_TIER1_NO_BREADTH = 10.0
FC_TIER2            = 10.0
FC_HELMINTH         = 3.0
FC_PRIORITY         = 2.0
RUN_PREVALENCE_MAX  = 0.5
RUN_PREVALENCE_RPM  = 300.0

ENVIRONMENTAL_FAMILIES = frozenset([
    "Weeksellaceae", "Sphaerotilaceae", "Deinococcaceae", "Azospirillaceae",
    "Pyriculariaceae", "Sclerotiniaceae", "Listeriaceae", "Comamonadaceae",
])
ENVIRONMENTAL_GENERA = frozenset([
    "Xanthomonas", "Caldimonas", "Diaphorobacter",
    "Exiguobacterium", "Skermanella", "Paracoccus",
])
PHAGE_PHYLA = frozenset(["Uroviricota", "Phixviricota"])
PRIORITY_TAXIDS = frozenset([
    1773, 1765, 78331, 746128, 5476, 5478, 5482,
    3071, 4909, 727, 1392, 637, 632, 777,
])

HUMAN_CALLS = frozenset([
    "HUMAN_CLASSIFIED", "HUMAN_CONFLICT",
    "HUMAN_KMER_BROAD_COVERAGE", "LOW_CONF_HUMAN_KMER",
])
UNINFORMATIVE_DISPLAY_TAXIDS = frozenset({
    2759, 2, 10239, 2157, 131567, 1, 33154, 4751,
})

TPOS_EXPECTED = {
    1639:  ("Listeria monocytogenes",           0.891,      "MUST_DETECT"),
    287:   ("Pseudomonas aeruginosa",            0.089,      "MUST_DETECT"),
    1423:  ("Bacillus subtilis",                 0.0089,     "MUST_DETECT"),
    4932:  ("Saccharomyces cerevisiae",          0.0089,     "MUST_DETECT"),
    562:   ("Escherichia coli",                  0.00089,    "SHOULD_DETECT"),
    28901: ("Salmonella enterica",               0.00089,    "SHOULD_DETECT"),
    1613:  ("Limosilactobacillus fermentum",     0.000089,   "BORDERLINE"),
    1351:  ("Enterococcus faecalis",             0.0000089,  "BORDERLINE"),
    5207:  ("Cryptococcus neoformans",           0.0000089,  "BORDERLINE"),
    1280:  ("Staphylococcus aureus",             0.00000089, "LOD_CHALLENGE"),
}

# ── Taxonomy databases (loaded once at module level) ──────────────────────────
parent_g: dict = {}
rank_g:   dict = {}
name_lookup_g: dict = {}
k2_name_df_g:  pl.DataFrame
taxmap_g:       pl.DataFrame

def _load_taxonomy_dbs() -> None:
    global parent_g, rank_g, name_lookup_g, k2_name_df_g, taxmap_g

    with open(P_KRAKEN2_NODES) as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if len(parts) < 3:
                continue
            try:
                t = int(parts[0].strip())
                parent_g[t] = int(parts[1].strip())
                rank_g[t]   = parts[2].strip()
            except ValueError:
                continue
    print(f"  nodes.dmp  : {len(parent_g):,} nodes")

    with open(P_KRAKEN2_NAMES) as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if len(parts) < 4:
                continue
            try:
                t = int(parts[0].strip())
            except ValueError:
                continue
            if parts[3].strip().rstrip("\t|").strip() == "scientific name":
                name_lookup_g[t] = parts[1].strip()
    print(f"  names.dmp  : {len(name_lookup_g):,} scientific names")

    eupath_rows = []
    with open(P_EUPATH_INSPECT) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                t = int(parts[4].strip())
            except ValueError:
                continue
            eupath_rows.append({"k2_taxid": t, "k2_sci_name": parts[5].strip()})
    k2_name_df_g = (
        pl.DataFrame(eupath_rows, schema={"k2_taxid": pl.Int64, "k2_sci_name": pl.Utf8})
        .unique(subset=["k2_taxid"], keep="first")
    )
    print(f"  EuPathDB   : {len(k2_name_df_g):,} taxa")

    taxmap_g = (
        pl.read_csv(P_TAXMAP, separator="\t")
        .unique(subset=["accession"], keep="first")
        .rename({"accession": "bam_rname"})
        .with_columns([
            pl.col("ncbi_taxid_strain").cast(pl.Int64).alias("bam_taxid_accession"),
            pl.col("ncbi_taxid_species").cast(pl.Int64).alias("bam_taxid_species"),
            pl.coalesce([
                pl.when(pl.col("species_name_kraken2").is_not_null() & (pl.col("species_name_kraken2") != "")).then(pl.col("species_name_kraken2")),
                pl.when(pl.col("bartlett_species").is_not_null()     & (pl.col("bartlett_species")     != "")).then(pl.col("bartlett_species")),
                pl.col("gtdb_species"),
            ]).alias("taxon_name"),
        ])
        .select(["bam_rname", "bam_taxid_accession", "bam_taxid_species", "taxon_name"])
    )
    print(f"  taxmap     : {len(taxmap_g):,} accessions")


# ── Helper functions ───────────────────────────────────────────────────────────
def _strip_mate_suffix(s):
    return MATE_SUFFIX_RE.sub("", s)


def _parse_kmer_string(kmer_str):
    if not kmer_str or kmer_str.strip() == "":
        return 0.0, 0.0, 0, 0
    total = 0; classified = 0; taxa_counts = {}
    hit_groups = 0; in_hit_run = False
    for token in kmer_str.split():
        if ":" not in token:
            continue
        taxid_s, count_s = token.split(":", 1)
        try:
            taxid = int(taxid_s); count = int(count_s)
        except ValueError:
            continue
        total += count
        is_informative = taxid not in UNINFORMATIVE_TAXIDS
        if is_informative:
            classified += count
            taxa_counts[taxid] = taxa_counts.get(taxid, 0) + count
            if not in_hit_run:
                hit_groups += 1; in_hit_run = True
        else:
            in_hit_run = False
    kmer_hit_frac  = classified / total if total > 0 else 0.0
    distinct_ratio = max(taxa_counts.values()) / classified if classified > 0 else 0.0
    return kmer_hit_frac, distinct_ratio, len(taxa_counts), hit_groups


def load_kraken2_per_read(path, db_tag):
    rows = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            status = parts[0]; rid = _strip_mate_suffix(parts[1])
            try:    taxid    = int(parts[2])
            except: taxid    = 0
            try:    read_len = int(parts[3])
            except: read_len = None
            khf, dr, n_taxa, hg = _parse_kmer_string(parts[4])
            rows.append({
                "read_id": rid, "status": status, "taxid": taxid,
                "read_len": read_len, "kmer_hit_frac": khf,
                "distinct_ratio": dr, "n_taxa_hit": n_taxa,
                "hit_groups": hg, "db": db_tag,
            })
    return pl.DataFrame(rows, schema={
        "read_id": pl.Utf8, "status": pl.Utf8, "taxid": pl.Int64,
        "read_len": pl.Int32, "kmer_hit_frac": pl.Float64,
        "distinct_ratio": pl.Float64, "n_taxa_hit": pl.Int32,
        "hit_groups": pl.Int32, "db": pl.Utf8,
    })


def classify_read(df, taxid_kingdom_df):
    df_k = df.join(
        taxid_kingdom_df.rename({"taxid": "k1_taxid", "kingdom": "k1_kingdom"}),
        on="k1_taxid", how="left",
    ).with_columns(pl.col("k1_kingdom").fill_null("Unknown"))

    UNRESOLVABLE = frozenset([
        0, 1, 2, 131567, 2759, 2157, 10239, 33154, 4751, 33090, 33208,
    ])
    _k3_classified = pl.col("k3_status") == "C" if "k3_status" in df_k.columns \
                     else pl.lit(False)
    df_k = df_k.filter(
        (pl.col("bam_rname").is_not_null() & (pl.col("k1_taxid") == 0)) |
        (~pl.col("k1_taxid").is_in(list(UNRESOLVABLE)))                  |
        _k3_classified
    )

    has_bam          = pl.col("bam_rname").is_not_null()
    has_k2           = pl.col("k2_status") == "C"
    has_k3           = pl.col("k3_status") == "C" if "k3_status" in df_k.columns \
                       else pl.lit(False)
    is_human         = pl.col("k1_taxid").is_in(list(HUMAN_TAXIDS))
    is_bacteria      = pl.col("k1_kingdom") == "Bacteria"
    is_eukaryote     = pl.col("k1_kingdom") == "Eukaryota"
    is_virus         = pl.col("k1_kingdom") == "Viruses"
    is_human_kingdom = pl.col("k1_kingdom") == "Human"
    has_breadth      = (
        (pl.col("n_windows_covered") >= BREADTH_MIN_WINDOWS) &
        (pl.col("breadth_fraction")  >= BREADTH_MIN_FRACTION)
    )
    k1_classified = pl.col("k1_status") == "C"

    k1_bacteria_quality = (
        (pl.col("k1_kmer_hit_frac")  >= KMER_HIT_FRAC_MIN)  &
        (pl.col("k1_distinct_ratio") >= DISTINCT_RATIO_MIN)  &
        (pl.col("k1_hit_groups")     >= HIT_GROUPS_MIN)      &
        (pl.col("k1_n_taxa_hit")     <= N_TAXA_MAX)
    )
    k1_eukaryote_quality = (
        (pl.col("k1_kmer_hit_frac")  >= 0.10) &
        (pl.col("k1_distinct_ratio") >= 0.30) &
        (pl.col("k1_hit_groups")     >= 1)    &
        (pl.col("k1_n_taxa_hit")     <= 8)
    )
    k1_virus_quality = (
        (pl.col("k1_kmer_hit_frac")  >= 0.10) &
        (pl.col("k1_distinct_ratio") >= 0.30) &
        (pl.col("k1_hit_groups")     >= 1)    &
        (pl.col("k1_n_taxa_hit")     <= 10)
    )
    k2_eukaryote_quality = (
        (pl.col("k2_kmer_hit_frac")  >= 0.10) &
        (pl.col("k2_distinct_ratio") >= 0.30) &
        (pl.col("k2_hit_groups")     >= 1)
    )
    k3_helminth_quality = (
        (pl.col("k3_kmer_hit_frac")  >= 0.10) &
        (pl.col("k3_distinct_ratio") >= 0.30) &
        (pl.col("k3_hit_groups")     >= 1)
    ) if "k3_kmer_hit_frac" in df_k.columns else pl.lit(False)

    return (
        df_k.with_columns(
            pl.when(k1_classified & (is_human | is_human_kingdom) & (pl.col("k1_kmer_hit_frac") >= HUMAN_KHF_HARD_CUTOFF) & ~has_breadth).then(pl.lit("HUMAN_CONFLICT"))
              .when(k1_classified & (is_human | is_human_kingdom) & (pl.col("k1_kmer_hit_frac") >= HUMAN_KHF_HARD_CUTOFF) & has_breadth).then(pl.lit("HUMAN_KMER_BROAD_COVERAGE"))
              .when(k1_classified & (is_human | is_human_kingdom) & (pl.col("k1_kmer_hit_frac") < HUMAN_KHF_HARD_CUTOFF) & (pl.col("k1_kmer_hit_frac") >= KMER_HIT_FRAC_MIN)).then(pl.lit("LOW_CONF_HUMAN_KMER"))
              .when(k1_classified & (is_human | is_human_kingdom)).then(pl.lit("HUMAN_CLASSIFIED"))
              .when(k1_classified & is_eukaryote & k1_eukaryote_quality & has_k3 & k3_helminth_quality).then(pl.lit("EUKARYOTE_helminth"))
              .when(k1_classified & is_eukaryote & k1_eukaryote_quality & has_k2 & k2_eukaryote_quality).then(pl.lit("EUKARYOTE_eupath"))
              .when(k1_classified & is_eukaryote & k1_eukaryote_quality & ~has_k2 & ~has_k3).then(pl.lit("EUKARYOTE_PlusPFOnly"))
              .when((pl.col("k1_status").is_null() | (pl.col("k1_status") == "U")) & ~has_bam & has_k3 & k3_helminth_quality).then(pl.lit("EUKARYOTE_K3_HELMINTH"))
              .when(k1_classified & is_virus & k1_virus_quality).then(pl.lit("VIRUS_K1_ONLY"))
              .when(k1_classified & is_bacteria & k1_bacteria_quality & has_bam).then(pl.lit("KRAKEN_CONFIRMED"))
              .when(k1_classified & is_bacteria & k1_bacteria_quality & ~has_bam).then(pl.lit("BACTERIA_K1_ONLY"))
              .when(k1_classified & is_bacteria & (pl.col("k1_kmer_hit_frac") >= KMER_HIT_FRAC_MIN) & (pl.col("k1_distinct_ratio") >= DISTINCT_RATIO_MIN) & has_bam).then(pl.lit("BAM_K1_LOW_CONF"))
              .when(k1_classified & is_bacteria & (pl.col("k1_kmer_hit_frac") >= KMER_HIT_FRAC_MIN) & (pl.col("k1_distinct_ratio") >= DISTINCT_RATIO_MIN) & ~has_bam).then(pl.lit("BACTERIA_LOW_CONF"))
              .when(has_bam & (pl.col("k1_status").is_null() | (pl.col("k1_status") == "U"))).then(pl.lit("BAM_ONLY"))
              .otherwise(pl.lit("KRAKEN_LOW_QUALITY"))
              .alias("read_call")
        )
        .with_columns(
            pl.when(pl.col("read_call") == "BAM_ONLY")
              .then(pl.lit("Bacteria"))
              .when(pl.col("read_call") == "EUKARYOTE_K3_HELMINTH")
              .then(pl.lit("Eukaryota"))
              .otherwise(pl.col("k1_kingdom"))
              .alias("kingdom")
        )
    )


_SUPERKINGDOM_ROOTS = {
    2: "Bacteria", 2157: "Archaea", 2759: "Eukaryota",
    10239: "Viruses", 9606: "Human",
}


def _resolve_kingdom(taxid):
    current = taxid; visited = set()
    while current not in visited:
        if current in _SUPERKINGDOM_ROOTS:
            return _SUPERKINGDOM_ROOTS[current]
        if current in (0, 1):
            break
        visited.add(current)
        current = parent_g.get(current, 1)
    return "Unknown"


_STOP_RANKS    = {"species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"}
_NO_COLLAPSE   = frozenset([32603, 32604, 10376, 1891762, 687331])
_LINEAGE_RANKS = {"superkingdom", "phylum", "class", "order", "family", "genus", "species"}


def _find_species(taxid):
    if taxid in _NO_COLLAPSE or rank_g.get(taxid, "no rank") in _STOP_RANKS:
        return taxid
    current = taxid; visited = set()
    while current not in visited and current not in (0, 1):
        if rank_g.get(current, "no rank") in _STOP_RANKS:
            return current
        visited.add(current)
        current = parent_g.get(current, 1)
    return taxid


def _get_lineage(taxid):
    lineage = {}; current = taxid; visited = set()
    while current not in visited and current not in (0, 1):
        r = rank_g.get(current, "no rank")
        if r in _LINEAGE_RANKS:
            lineage[r] = current
        visited.add(current)
        current = parent_g.get(current, 1)
    return lineage


def run_sample(s):
    sid = s["sample_id"]

    k1 = load_kraken2_per_read(s["p_k1"], db_tag="pluspf")
    k2 = load_kraken2_per_read(s["p_k2"], db_tag="eupath")
    k3 = load_kraken2_per_read(s["p_k3"], db_tag="helminths") \
         if Path(s["p_k3"]).exists() else None

    _frames = [k1.select("taxid"), k2.select("taxid")]
    if k3 is not None:
        _frames.append(k3.select("taxid"))
    _all_taxids = (
        pl.concat(_frames).unique().filter(pl.col("taxid") > 1)["taxid"].to_list()
    )
    taxid_kingdom_df = pl.DataFrame(
        [{"taxid": t, "kingdom": _resolve_kingdom(t)} for t in _all_taxids],
        schema={"taxid": pl.Int64, "kingdom": pl.Utf8},
    )

    rows = []; window_hits = {}; ref_lengths = {}
    with pysam.AlignmentFile(s["p_bam"], "rb") as bam:
        for ref, length in zip(bam.references, bam.lengths):
            ref_lengths[ref] = length
            window_hits[ref] = [0] * max(1, length // WINDOW_SIZE)
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            rid  = MATE_SUFFIX_RE.sub("", read.query_name)
            rnam = read.reference_name
            rows.append({
                "read_id": rid, "bam_rname": rnam,
                "bam_mapq": read.mapping_quality,
                "bam_nm":   read.get_tag("NM") if read.has_tag("NM") else None,
                "bam_rlen": read.query_length,
                "bam_alen": read.query_alignment_length,
                "bam_flag": read.flag,
            })
            n_win = len(window_hits[rnam])
            w = min(read.reference_start // WINDOW_SIZE, n_win - 1)
            window_hits[rnam][w] += 1

    breadth_df = pl.DataFrame(
        [{"bam_rname": ref, "ref_length": ref_lengths[ref],
          "n_windows_total":   len(hits),
          "n_windows_covered": sum(1 for h in hits if h >= MIN_DEPTH),
          "breadth_fraction":  (sum(1 for h in hits if h >= MIN_DEPTH) * WINDOW_SIZE)
                               / ref_lengths[ref] if ref_lengths[ref] > 0 else 0.0}
         for ref, hits in window_hits.items()],
        schema={"bam_rname": pl.Utf8, "ref_length": pl.Int64,
                "n_windows_total": pl.Int32, "n_windows_covered": pl.Int32,
                "breadth_fraction": pl.Float64},
    )

    bam_raw = (
        pl.DataFrame(rows, schema={
            "read_id": pl.Utf8, "bam_rname": pl.Utf8, "bam_mapq": pl.Int32,
            "bam_nm": pl.Int32, "bam_rlen": pl.Int32,
            "bam_alen": pl.Int32, "bam_flag": pl.Int32,
        })
        .with_columns(
            pl.when(pl.col("bam_rlen") > 0)
              .then(pl.col("bam_nm") / pl.col("bam_rlen"))
              .otherwise(None).alias("bam_nm_rate")
        )
        .join(breadth_df, on="bam_rname", how="left")
    )

    high_confidence = (
        (pl.col("bam_mapq")    >= MAPQ_MIN)    &
        (pl.col("bam_nm_rate") <= NM_RATE_MAX) &
        (pl.col("bam_alen")    >= ALEN_MIN)
    )
    low_mapq_broad_coverage = (
        (pl.col("bam_mapq")          >= MAPQ_FLOOR)          &
        (pl.col("bam_mapq")          <  MAPQ_MIN)            &
        (pl.col("bam_nm_rate")       <= NM_RATE_MAX)         &
        (pl.col("bam_alen")          >= ALEN_MIN)            &
        (pl.col("n_windows_covered") >= BREADTH_MIN_WINDOWS) &
        (pl.col("breadth_fraction")  >= BREADTH_MIN_FRACTION)
    )
    bam_raw = (
        bam_raw
        .with_columns(
            pl.when(high_confidence).then(pl.lit("HIGH_CONF"))
              .when(low_mapq_broad_coverage).then(pl.lit("BREADTH_RESCUED"))
              .otherwise(pl.lit("FAILED")).alias("mapq_filter_reason")
        )
        .filter(pl.col("mapq_filter_reason") != "FAILED")
        .sort("bam_mapq", descending=True)
        .group_by("read_id").first()
    )
    bam_df = bam_raw.join(taxmap_g, on="bam_rname", how="left")

    def _prefix(df, tag):
        return df.rename({c: f"{tag}_{c}" for c in df.columns if c != "read_id"})

    k1p = _prefix(k1, "k1")
    k2p = _prefix(k2, "k2")
    k3p = _prefix(k3, "k3") if k3 is not None else None

    merged_k1 = k1p.filter(pl.col("k1_status") == "C").join(k2p, on="read_id", how="left")
    if k3p is not None:
        merged_k1 = merged_k1.join(k3p, on="read_id", how="left")
    merged_k1 = merged_k1.join(bam_df, on="read_id", how="left")

    merged_bam_only = (
        bam_df.join(merged_k1.select("read_id"), on="read_id", how="anti")
        .join(k1p, on="read_id", how="left")
        .join(k2p, on="read_id", how="left")
    )
    if k3p is not None:
        merged_bam_only = merged_bam_only.join(k3p, on="read_id", how="left")

    _to_concat = [merged_k1, merged_bam_only]
    if k3p is not None:
        _seen = pl.concat([merged_k1.select("read_id"),
                           merged_bam_only.select("read_id")]).unique()
        merged_k3_only = (
            k3p.filter(pl.col("k3_status") == "C")
            .join(_seen, on="read_id", how="anti")
            .join(k1p,  on="read_id", how="left")
            .join(k2p,  on="read_id", how="left")
            .join(bam_df, on="read_id", how="left")
        )
        if len(merged_k3_only) > 0:
            _to_concat.append(merged_k3_only)

    merged = pl.concat(_to_concat, how="diagonal_relaxed")
    passed  = classify_read(merged, taxid_kingdom_df)
    nonhuman = passed.filter(~pl.col("read_call").is_in(list(HUMAN_CALLS)))

    _k1_taxids = set(nonhuman["k1_taxid"].drop_nulls().unique().to_list())
    taxid_name_df = pl.DataFrame(
        [{"k1_taxid": t, "k1_sci_name": name_lookup_g.get(t, str(t))} for t in _k1_taxids],
        schema={"k1_taxid": pl.Int64, "k1_sci_name": pl.Utf8},
    )
    _k3_name_df = None
    if "k3_taxid" in nonhuman.columns:
        _k3_taxids = set(nonhuman["k3_taxid"].drop_nulls().unique().to_list())
        _k3_name_df = pl.DataFrame(
            [{"k3_taxid": t, "k3_sci_name": name_lookup_g.get(t, str(t))} for t in _k3_taxids],
            schema={"k3_taxid": pl.Int64, "k3_sci_name": pl.Utf8},
        )

    _joined = (
        nonhuman
        .join(taxid_name_df, on="k1_taxid", how="left")
        .join(k2_name_df_g,  on="k2_taxid", how="left")
    )
    if _k3_name_df is not None:
        _joined = _joined.join(_k3_name_df, on="k3_taxid", how="left")
    else:
        _joined = _joined.with_columns(pl.lit(None).cast(pl.Utf8).alias("k3_sci_name"))

    summary = (
        _joined
        .with_columns([
            pl.when(pl.col("read_call") == "EUKARYOTE_eupath")
              .then(pl.col("k2_taxid"))
              .when(pl.col("read_call").is_in(["EUKARYOTE_helminth", "EUKARYOTE_K3_HELMINTH"]))
              .then(pl.col("k3_taxid").cast(pl.Int64))
              .when(pl.col("read_call") == "BAM_ONLY")
              .then(pl.col("bam_taxid_species"))
              .otherwise(pl.col("k1_taxid"))
              .alias("display_taxid"),
            pl.when(pl.col("read_call") == "EUKARYOTE_eupath")
              .then(pl.col("k2_sci_name"))
              .when(pl.col("read_call").is_in(["EUKARYOTE_helminth", "EUKARYOTE_K3_HELMINTH"]))
              .then(pl.col("k3_sci_name"))
              .when(pl.col("read_call") == "BAM_ONLY")
              .then(pl.col("taxon_name"))
              .otherwise(pl.col("k1_sci_name"))
              .alias("display_name"),
        ])
        .group_by(["read_call", "kingdom", "display_taxid", "display_name"])
        .agg([
            pl.len().alias("n_reads"),
            pl.col("k1_kmer_hit_frac").mean().alias("mean_k1_khf"),
            pl.col("k2_kmer_hit_frac").mean().alias("mean_k2_khf"),
            pl.col("k1_hit_groups").mean().alias("mean_k1_hg"),
            pl.col("k2_hit_groups").mean().alias("mean_k2_hg"),
            pl.col("n_windows_covered").first().alias("n_windows_covered"),
            pl.col("breadth_fraction").first().alias("breadth_fraction"),
        ])
    )

    summary = summary.filter(pl.col("kingdom") != "Unknown")
    summary = summary.filter(
        ~pl.col("display_taxid").is_in(list(UNINFORMATIVE_DISPLAY_TAXIDS))
    )

    _taxids   = summary["display_taxid"].drop_nulls().unique().to_list()
    _sp_map   = pl.DataFrame(
        [{"display_taxid": t, "species_taxid": _find_species(t)} for t in _taxids],
        schema={"display_taxid": pl.Int64, "species_taxid": pl.Int64},
    )
    _sp_names = pl.DataFrame(
        [{"species_taxid": t, "species_name": name_lookup_g.get(t, str(t))}
         for t in _sp_map["species_taxid"].unique().to_list()],
        schema={"species_taxid": pl.Int64, "species_name": pl.Utf8},
    )

    sc = (
        summary
        .join(_sp_map,   on="display_taxid", how="left")
        .join(_sp_names, on="species_taxid",  how="left")
        .with_columns([
            pl.coalesce([pl.col("species_taxid"), pl.col("display_taxid")]).alias("group_taxid"),
            pl.coalesce([pl.col("species_name"),  pl.col("display_name")]).alias("group_name"),
        ])
        .group_by(["kingdom", "group_taxid", "group_name"])
        .agg([
            pl.col("n_reads").sum().alias("n_reads_total"),
            pl.when(pl.col("read_call").is_in([
                "KRAKEN_CONFIRMED", "EUKARYOTE_eupath", "EUKARYOTE_helminth",
            ])).then(pl.col("n_reads")).otherwise(0).sum().alias("n_reads_tier1"),
            pl.when(pl.col("read_call").is_in([
                "BAM_K1_LOW_CONF", "BACTERIA_K1_ONLY",
                "EUKARYOTE_K1_ONLY", "VIRUS_K1_ONLY",
                "BAM_ONLY", "EUKARYOTE_K3_HELMINTH",
            ])).then(pl.col("n_reads")).otherwise(0).sum().alias("n_reads_tier2"),
            pl.when(pl.col("read_call").is_in([
                "BACTERIA_LOW_CONF", "HUMAN_KMER_BROAD_COVERAGE",
            ])).then(pl.col("n_reads")).otherwise(0).sum().alias("n_reads_tier3"),
            pl.col("read_call").str.split("|").list.explode().unique().sort().str.join("|").alias("evidence_calls"),
            pl.col("mean_k1_khf").mean(), pl.col("mean_k2_khf").mean(),
            pl.col("n_windows_covered").max(), pl.col("breadth_fraction").max(),
        ])
        .with_columns([
            pl.lit(sid).alias("sample_id"),
            pl.lit(s["is_blank"]).alias("is_blank"),
            pl.when(pl.col("n_reads_tier1") > 0).then(1)
              .when(pl.col("n_reads_tier2") > 0).then(2)
              .when(pl.col("n_reads_tier3") > 0).then(3)
              .otherwise(4).alias("best_tier"),
            pl.when(pl.col("group_taxid").map_elements(
                lambda t: rank_g.get(t, "no rank") == "genus", return_dtype=pl.Boolean)
            ).then(pl.lit("genus"))
            .otherwise(pl.col("group_taxid").map_elements(
                lambda t: rank_g.get(t, "no rank"), return_dtype=pl.Utf8)
            ).alias("tax_rank"),
        ])
        .with_columns([
            pl.when(pl.col("best_tier") == 4).then(0).otherwise(pl.col("n_windows_covered")).alias("n_windows_covered"),
            pl.when(pl.col("best_tier") == 4).then(0.0).otherwise(pl.col("breadth_fraction")).alias("breadth_fraction"),
        ])
        .filter(pl.col("group_taxid").map_elements(
            lambda t: rank_g.get(t, "no rank") not in NONREPORTABLE_RANKS,
            return_dtype=pl.Boolean,
        ))
    )

    if GENUS_LEVEL_COLLAPSE:
        _genus_needed = set(GENUS_LEVEL_COLLAPSE.values())
        _genus_names  = {t: name_lookup_g.get(t, str(t)) for t in _genus_needed}
        sc = (
            sc.with_columns(pl.col("group_taxid").replace(GENUS_LEVEL_COLLAPSE).alias("group_taxid"))
              .with_columns(
                  pl.when(pl.col("group_taxid").is_in(list(_genus_needed)))
                    .then(pl.col("group_taxid").map_elements(lambda t: _genus_names.get(t, ""), return_dtype=pl.Utf8))
                    .otherwise(pl.col("group_name")).alias("group_name")
              )
              .group_by(["kingdom", "group_taxid", "group_name", "sample_id", "is_blank"])
              .agg([
                  pl.col("n_reads_total").sum(), pl.col("n_reads_tier1").sum(),
                  pl.col("n_reads_tier2").sum(), pl.col("n_reads_tier3").sum(),
                  pl.col("evidence_calls").str.split("|").list.explode().unique().sort().str.join("|"),
                  pl.col("mean_k1_khf").mean(), pl.col("mean_k2_khf").mean(),
                  pl.col("n_windows_covered").max(), pl.col("breadth_fraction").max(),
                  pl.col("best_tier").min(), pl.col("tax_rank").first(),
              ])
        )

    _lin_rows = []
    for _t in sc["group_taxid"].drop_nulls().unique().to_list():
        _lin = _get_lineage(_t)
        _lin_rows.append({
            "group_taxid": _t,
            "tax_phylum":  name_lookup_g.get(_lin.get("phylum",  0), ""),
            "tax_class":   name_lookup_g.get(_lin.get("class",   0), ""),
            "tax_order":   name_lookup_g.get(_lin.get("order",   0), ""),
            "tax_family":  name_lookup_g.get(_lin.get("family",  0), ""),
            "tax_genus":   name_lookup_g.get(_lin.get("genus",   0), ""),
            "tax_species": name_lookup_g.get(_lin.get("species", 0), ""),
            "taxonomy": "|".join(filter(None, [
                name_lookup_g.get(_lin.get(r, 0), "")
                for r in ["phylum", "class", "order", "family", "genus", "species"]
            ])),
        })
    lineage_df = pl.DataFrame(_lin_rows, schema={
        "group_taxid": pl.Int64, "tax_phylum": pl.Utf8, "tax_class": pl.Utf8,
        "tax_order": pl.Utf8, "tax_family": pl.Utf8, "tax_genus": pl.Utf8,
        "tax_species": pl.Utf8, "taxonomy": pl.Utf8,
    })

    return (
        sc.join(lineage_df, on="group_taxid", how="left")
          .filter(
              ((pl.col("best_tier") == 1) & (pl.col("n_reads_total") >= 1)) |
              ((pl.col("best_tier") == 2) & (pl.col("n_reads_total") >= 3)) |
              ((pl.col("best_tier") == 3) & (pl.col("n_reads_total") >= 5))
          )
    )


def _paths(folder: str, sample: str) -> dict:
    return {
        "p_k1":  f"{BASE}/{folder}/output/kraken2/{sample}.krk",
        "p_k2":  f"{BASE}/{folder}/output/kraken2_eupathdb/{sample}.krk",
        "p_k3":  f"{BASE}/{folder}/output/kraken2_helminth/{sample}.krk",
        "p_bam": f"{BASE}/{folder}/output/bam/{sample}.bam",
    }


def cohort_from_samplesheet(folder: str, run: int, p_sheet: str | None = None) -> list[dict]:
    if p_sheet is None:
        p_sheet = f"{BASE}/{folder}/SampleSheet.csv"
    BLANK_PREFIXES = ("TENV", "TPOS")

    def _lib_type(name: str) -> str:
        return "RNA" if name.upper().endswith(("_ARN", "-ARN")) else "DNA"

    entries = []
    in_data = False
    with open(p_sheet, newline="", encoding="utf-8-sig") as fh:
        for line in fh:
            line = line.strip()
            if line == "[Data]":
                in_data = True; continue
            if not in_data or not line or line.startswith("Sample_ID"):
                continue
            parts = line.split(",")
            if not parts or not parts[0].strip():
                continue
            sample_id = parts[0].strip()
            patient   = parts[1] if len(parts) > 1 else ""
            qbit      = parts[2] if len(parts) > 2 else ""
            entries.append({
                "sample_id":   sample_id,
                "folder":      folder,
                "sample":      sample_id,
                "dna_ng":      qbit,
                "patient":     patient,
                "sample_type": "blood",
                "is_blank":    any(sample_id.startswith(p) for p in BLANK_PREFIXES),
                "run":         run,
                "lib_type":    _lib_type(sample_id),
            })
    for s in entries:
        s.update(_paths(s["folder"], s["sample"]))
    return entries


def run_track(samples: list, label: str) -> tuple[pl.DataFrame, list]:
    results = []; errors = []
    for i, s in enumerate(samples):
        print(f"  [{i+1:>2}/{len(samples)}] {s['sample_id']}", end="  →  ", flush=True)
        try:
            row = run_sample(s)
            results.append(row)
            print(f"OK  ({len(row)} rows)")
        except Exception as e:
            errors.append({"sample_id": s["sample_id"], "error": str(e)})
            print(f"ERROR  {e}")
    df = pl.concat(results, how="diagonal_relaxed") if results else pl.DataFrame()
    return df, errors


def _count_microbial_reads(path: str) -> int:
    n = 0
    with open(path) as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                taxid = int(parts[2])
            except ValueError:
                continue
            if taxid not in HUMAN_TAXIDS:
                n += 1
    return n


def _read_nreads(folder: str) -> pl.DataFrame:
    path = os.path.join(BASE, folder, "output", "Nreads.csv")
    return (
        pl.read_csv(path)
        .rename({"sample": "sample_id", "dedup_reads": "total_reads"})
        .with_columns(pl.col("sample_id").str.replace(r"_S[0-9]+.*$", ""))
    )


def run_tpos_qc(
    tpos_result:      pl.DataFrame,
    run:              int,
    lib_type:         str,
    tpos_total_reads: int,
    history:          pl.DataFrame | None = None,
) -> pl.DataFrame:
    expected_df = pl.DataFrame([
        {"group_taxid": taxid, "expected_name": name,
         "expected_fraction": frac, "detection_tier": tier}
        for taxid, (name, frac, tier) in TPOS_EXPECTED.items()
    ], schema={
        "group_taxid": pl.Int64, "expected_name": pl.Utf8,
        "expected_fraction": pl.Float64, "detection_tier": pl.Utf8,
    })

    observed = (
        tpos_result
        .with_columns(
            (pl.col("n_reads_total") / tpos_total_reads).alias("observed_fraction")
        )
        .select(["group_taxid", "group_name", "n_reads_total",
                 "observed_fraction", "best_tier"])
    )

    comparison = (
        expected_df
        .join(observed, on="group_taxid", how="left")
        .with_columns([
            pl.col("observed_fraction").fill_null(0.0),
            pl.col("n_reads_total").fill_null(0),
            pl.col("observed_fraction").is_not_null().alias("detected"),
            pl.when(pl.col("observed_fraction") > 0)
              .then((pl.col("observed_fraction") / pl.col("expected_fraction")).log(base=2.0))
              .otherwise(None)
              .alias("log2_fold_change"),
        ])
        .sort("expected_fraction", descending=True)
    )

    lod_fraction = comparison.filter(pl.col("detected"))["expected_fraction"].min()

    _detected = comparison.filter(pl.col("observed_fraction") > 0)
    spearman_rho = None
    if len(_detected) >= 3:
        _exp = [math.log10(x) for x in _detected["expected_fraction"].to_list()]
        _obs = [math.log10(x) for x in _detected["observed_fraction"].to_list()]
        n = len(_exp)
        def _rank(lst):
            sorted_idx = sorted(range(n), key=lambda i: lst[i])
            ranks = [0.0] * n
            for rank, idx in enumerate(sorted_idx):
                ranks[idx] = rank + 1
            return ranks
        _r_exp = _rank(_exp); _r_obs = _rank(_obs)
        _mean_r = (n + 1) / 2
        _num = sum((_r_exp[i]-_mean_r)*(_r_obs[i]-_mean_r) for i in range(n))
        _den = (sum((r-_mean_r)**2 for r in _r_exp) *
                sum((r-_mean_r)**2 for r in _r_obs)) ** 0.5
        spearman_rho = _num / _den if _den > 0 else 0.0

    must_detect_rate = (
        comparison.filter(pl.col("detection_tier") == "MUST_DETECT")["detected"].mean()
    )
    qc_flag = (
        "FAIL" if must_detect_rate < 1.0
        else "FAIL" if (spearman_rho or 0) < 0.6
        else "WARN" if (spearman_rho or 1) < 0.8
        else "WARN" if lod_fraction and lod_fraction > 0.001
        else "PASS"
    )

    print(f"\n{'═'*60}")
    print(f"  TPOS QC — ZymoBIOMICS Std II Log  |  Run {run}  |  {lib_type}")
    print(f"{'═'*60}")
    print(f"  Overall QC flag      : {qc_flag}")
    print(f"  MUST_DETECT rate     : {must_detect_rate:.0%}")
    if lod_fraction:
        print(f"  LOD (lowest detected): {lod_fraction:.2e}")
    if spearman_rho is not None:
        print(f"  Spearman rho (log)   : {spearman_rho:.3f}")

    return pl.DataFrame([{
        "run": run, "lib_type": lib_type, "qc_flag": qc_flag,
        "must_detect_rate": must_detect_rate, "lod_fraction": lod_fraction,
        "spearman_rho": spearman_rho,
        "n_detected": int(comparison["detected"].sum()),
        "n_expected": len(comparison),
    }]), comparison


def _apply_decontam(run_cohort_norm: pl.DataFrame) -> pl.DataFrame:
    blank_db_ref = (
        pl.read_parquet(P_CONTAMINANT_DB)
        .filter(pl.col("n_blanks_present") >= MIN_BLANKS_N)
        .select(["group_taxid", "max_rpm_blank", "n_blanks_present"])
        .rename({"max_rpm_blank": "blank_db_rpm"})
    )
    run_blank_rpm_ref = (
        run_cohort_norm
        .filter(pl.col("is_blank"))
        .group_by("group_taxid")
        .agg(pl.col("rpm").max().alias("run_blank_rpm"))
    )
    _patients = run_cohort_norm.filter(~pl.col("is_blank"))
    _n_pat    = _patients["sample_id"].n_unique()

    run_prevalence = (
        _patients
        .group_by("group_taxid")
        .agg([
            pl.len().alias("n_samples_with_taxid"),
            pl.col("rpm").median().alias("median_rpm_run"),
        ])
        .with_columns(
            (pl.col("n_samples_with_taxid") / _n_pat).alias("run_prevalence")
        )
    )

    decontam = (
        _patients
        .join(blank_db_ref,      on="group_taxid", how="left")
        .join(run_blank_rpm_ref, on="group_taxid", how="left")
        .join(run_prevalence,    on="group_taxid", how="left")
        .with_columns([
            pl.col("blank_db_rpm"   ).fill_null(0.0),
            pl.col("run_blank_rpm"  ).fill_null(0.0),
            pl.col("run_prevalence" ).fill_null(0.0),
            pl.col("median_rpm_run" ).fill_null(0.0),
        ])
        .with_columns(
            pl.max_horizontal("blank_db_rpm", "run_blank_rpm").alias("reference_rpm")
        )
        .with_columns(
            pl.when(pl.col("reference_rpm") == 0)
              .then(pl.lit(float("inf")))
              .otherwise(pl.col("rpm") / pl.col("reference_rpm"))
              .alias("fold_change")
        )
    )

    _is_ebv           = pl.col("group_taxid") == EBV_TAXID
    _is_env_family    = pl.col("tax_family").is_in(list(ENVIRONMENTAL_FAMILIES))
    _is_env_genus     = pl.col("tax_genus").is_in(list(ENVIRONMENTAL_GENERA))
    _is_phage         = pl.col("tax_phylum").is_in(list(PHAGE_PHYLA))
    _is_priority      = pl.col("group_taxid").is_in(list(PRIORITY_TAXIDS))
    _is_helminth      = pl.col("evidence_calls").str.contains("EUKARYOTE_helminth|EUKARYOTE_K3_HELMINTH")
    _no_blank         = pl.col("reference_rpm") == 0
    _has_breadth      = (
        (pl.col("breadth_fraction")  >= 0.01) &
        (pl.col("n_windows_covered") >= 3)
    )
    _is_run_prevalent = (
        (pl.col("run_prevalence")  >= RUN_PREVALENCE_MAX) &
        (pl.col("median_rpm_run")  <  RUN_PREVALENCE_RPM)
    )

    decontam = decontam.with_columns(
        pl.when(_is_ebv).then(pl.lit("FILTERED_EBV"))
        .when(_is_phage).then(pl.lit("FILTERED_ENVIRONMENTAL"))
        .when(_is_env_family | _is_env_genus).then(pl.lit("FILTERED_ENVIRONMENTAL"))
        .when(pl.col("fold_change") < 1.0).then(pl.lit("FILTERED_BACKGROUND"))
        .when(pl.col("best_tier") == 3).then(pl.lit("FILTERED_LOW_CONF"))
        .when(_is_run_prevalent).then(pl.lit("FILTERED_RUN_PREVALENT"))
        .when(_is_priority & (pl.col("fold_change") >= FC_PRIORITY)).then(pl.lit("PASS"))
        .when(_is_priority & _no_blank).then(pl.lit("PASS"))
        .when(_is_helminth & (pl.col("fold_change") >= FC_HELMINTH)).then(pl.lit("PASS"))
        .when(_is_helminth & _no_blank).then(pl.lit("PASS"))
        .when(_is_helminth).then(pl.lit("FILTERED_BACKGROUND"))
        .when(
            (pl.col("best_tier") == 1) & _has_breadth &
            ((pl.col("fold_change") >= FC_TIER1_BREADTH) | _no_blank)
        ).then(pl.lit("PASS"))
        .when(
            (pl.col("best_tier") == 1) & ~_has_breadth &
            ((pl.col("fold_change") >= FC_TIER1_NO_BREADTH) | _no_blank)
        ).then(pl.lit("PASS"))
        .when(
            (pl.col("best_tier") == 2) &
            ((pl.col("fold_change") >= FC_TIER2) | _no_blank)
        ).then(pl.lit("PASS"))
        .otherwise(pl.lit("FILTERED_BACKGROUND"))
        .alias("decontam_decision")
    )

    decontam = decontam.with_columns(
        pl.when(
            (pl.col("decontam_decision") == "FILTERED_BACKGROUND") &
            (pl.col("run_blank_rpm") > pl.col("blank_db_rpm"))
        )
        .then(pl.lit("FILTERED_RUN_SPECIFIC"))
        .otherwise(pl.col("decontam_decision"))
        .alias("decontam_decision")
    )
    return decontam


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description="CUB post-analysis pipeline")
    ap.add_argument("--run-id", type=int, default=0,
                    help="Integer run identifier (used in TPOS QC history tracking)")
    args = ap.parse_args()

    run_dir = Path(os.getcwd())
    folder  = run_dir.name
    out_dir = run_dir / "output" / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading taxonomy databases...")
    _load_taxonomy_dbs()

    print(f"\nParsing SampleSheet for run folder: {folder}")
    cohort = cohort_from_samplesheet(folder, args.run_id)

    patients_dna = [s for s in cohort if not s["is_blank"] and "TPOS" not in s["sample_id"] and s["lib_type"] == "DNA"]
    tenv_dna     = [s for s in cohort if s["is_blank"]     and "TPOS" not in s["sample_id"] and s["lib_type"] == "DNA"]
    tpos_dna     = [s for s in cohort if "TPOS" in s["sample_id"]                           and s["lib_type"] == "DNA"]

    print(f"  Patients DNA : {len(patients_dna)}")
    print(f"  TENV     DNA : {len(tenv_dna)}")
    print(f"  TPOS     DNA : {len(tpos_dna)}")

    # ── TPOS QC (isolated, never enters decontam) ─────────────────────────────
    tpos_qc_rows = []
    tpos_comparisons = []
    tpos_raw_results = []
    for s in tpos_dna:
        print(f"\n── TPOS QC: {s['sample_id']} ────────────────────────────────────")
        try:
            tpos_raw = run_sample(s)
            tpos_raw_results.append(tpos_raw)
            total = sum(1 for _ in open(s["p_k1"]) if _.strip())
            qc_row, comparison = run_tpos_qc(
                tpos_result      = tpos_raw,
                run              = args.run_id,
                lib_type         = "DNA",
                tpos_total_reads = total,
            )
            tpos_qc_rows.append(qc_row)
            tpos_comparisons.append(comparison)
        except Exception as e:
            print(f"  ERROR: {e}")

    if tpos_qc_rows:
        tpos_qc_df = pl.concat(tpos_qc_rows)
        tpos_qc_df.write_parquet(out_dir / "tpos_qc.parquet")
        print(f"\nWrote TPOS QC → {out_dir}/tpos_qc.parquet")
        if tpos_raw_results and tpos_comparisons:
            write_tpos_html(
                tpos_run_sample_df = tpos_raw_results[0],
                comparison_df      = tpos_comparisons[0],
                qc_row             = tpos_qc_rows[0],
                run_id             = args.run_id,
                out_path           = out_dir / "tpos_report.html",
            )
            print(f"Wrote TPOS HTML → {out_dir}/tpos_report.html")

    # ── Evidence integration: patients + blanks ───────────────────────────────
    print("\n── Evidence integration: patients + TENV ────────────────────────")
    run_cohort, errors = run_track(patients_dna + tenv_dna, label="cohort")
    if errors:
        print(f"  {len(errors)} samples failed:")
        for e in errors:
            print(f"    {e['sample_id']}: {e['error']}")

    # ── RPM normalisation ─────────────────────────────────────────────────────
    print("\nComputing microbial reads + RPM normalisation...")
    nreads_df = _read_nreads(folder)
    lib_sizes = []
    for s in patients_dna + tenv_dna:
        microbial = _count_microbial_reads(s["p_k1"])
        nrow = nreads_df.filter(pl.col("sample_id") == s["sample_id"])
        lib_sizes.append({
            "sample_id":       s["sample_id"],
            "microbial_reads": microbial,
            "Library_size":    nrow["total_reads"][0] if len(nrow) else None,
        })
    lib_df = pl.DataFrame(lib_sizes)

    run_cohort_norm = (
        run_cohort
        .join(lib_df.select(["sample_id", "microbial_reads", "Library_size"]),
              on="sample_id", how="left")
        .with_columns(
            (pl.col("n_reads_total") / pl.col("microbial_reads") * 1_000_000)
            .alias("rpm")
        )
    )

    # ── Decontam ──────────────────────────────────────────────────────────────
    print("\nApplying decontamination decision chain...")
    decontam = _apply_decontam(run_cohort_norm)

    _counts = decontam.group_by("decontam_decision").agg(pl.len().alias("n")).sort("n", descending=True)
    print("  decontam_decision summary:")
    for row in _counts.iter_rows(named=True):
        print(f"    {row['decontam_decision']:<30} {row['n']:>6}")

    decontam.write_parquet(out_dir / "cohort_norm.parquet")
    print(f"\nWrote {len(decontam)} rows → {out_dir}/cohort_norm.parquet")

    # Touch the HTML output file if no TPOS was present (Snakemake requires all outputs)
    html_out = out_dir / "tpos_report.html"
    if not html_out.exists():
        html_out.write_text("<html><body><p>No TPOS sample in this run.</p></body></html>")
    tpos_out = out_dir / "tpos_qc.parquet"
    if not tpos_out.exists():
        pl.DataFrame(schema={
            "run": pl.Int64, "lib_type": pl.Utf8, "qc_flag": pl.Utf8,
            "must_detect_rate": pl.Float64, "lod_fraction": pl.Float64,
            "spearman_rho": pl.Float64, "n_detected": pl.Int64, "n_expected": pl.Int64,
        }).write_parquet(tpos_out)


if __name__ == "__main__":
    main()
