import marimo

__generated_with = "0.23.5"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    import re
    import pysam
    from pathlib import Path

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
    P_OUTPUT_DIR     = "/mnt/san/microbio/apps/cub/resources/decontam"
    Path(P_OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    return (
        ALEN_MIN,
        BREADTH_MIN_FRACTION,
        BREADTH_MIN_WINDOWS,
        DISTINCT_RATIO_MIN,
        GENUS_LEVEL_COLLAPSE,
        HIT_GROUPS_MIN,
        HUMAN_KHF_HARD_CUTOFF,
        HUMAN_TAXIDS,
        KMER_HIT_FRAC_MIN,
        MAPQ_FLOOR,
        MAPQ_MIN,
        MATE_SUFFIX_RE,
        MIN_DEPTH,
        NM_RATE_MAX,
        NONREPORTABLE_RANKS,
        N_TAXA_MAX,
        P_EUPATH_INSPECT,
        P_KRAKEN2_NAMES,
        P_KRAKEN2_NODES,
        P_OUTPUT_DIR,
        P_TAXMAP,
        Path,
        UNINFORMATIVE_TAXIDS,
        WINDOW_SIZE,
        pl,
        pysam,
    )


@app.cell
def _(MATE_SUFFIX_RE, UNINFORMATIVE_TAXIDS, pl):
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

    return (load_kraken2_per_read,)


@app.cell
def _(
    BREADTH_MIN_FRACTION,
    BREADTH_MIN_WINDOWS,
    DISTINCT_RATIO_MIN,
    HIT_GROUPS_MIN,
    HUMAN_KHF_HARD_CUTOFF,
    HUMAN_TAXIDS,
    KMER_HIT_FRAC_MIN,
    N_TAXA_MAX,
    pl,
):
    def classify_read(df, taxid_kingdom_df):
        df_k = df.join(
            taxid_kingdom_df.rename({"taxid": "k1_taxid", "kingdom": "k1_kingdom"}),
            on="k1_taxid", how="left",
        ).with_columns(pl.col("k1_kingdom").fill_null("Unknown"))

        UNRESOLVABLE = frozenset([
            0, 1, 2, 131567, 2759, 2157, 10239, 33154, 4751, 33090, 33208,
        ])
        # K3-only reads have k1_taxid=0 and no BAM — must be explicitly passed through
        _k3_classified = pl.col("k3_status") == "C" if "k3_status" in df_k.columns \
                         else pl.lit(False)
        df_k = df_k.filter(
            (pl.col("bam_rname").is_not_null() & (pl.col("k1_taxid") == 0)) |
            (~pl.col("k1_taxid").is_in(list(UNRESOLVABLE)))                  |
            _k3_classified
        )

        has_bam          = pl.col("bam_rname").is_not_null()
        has_k2           = pl.col("k2_status") == "C"
        # K3 columns exist only when K3 was loaded for this sample
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
        # K3 quality gate — same thresholds as K2 eukaryote
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
                  # ── Eukaryota: K1+K3 double confirmation (tier1) ─────────────
                  .when(k1_classified & is_eukaryote & k1_eukaryote_quality & has_k3 & k3_helminth_quality).then(pl.lit("EUKARYOTE_helminth"))
                  # ── Eukaryota: K1+K2 confirmation (tier1) ────────────────────
                  .when(k1_classified & is_eukaryote & k1_eukaryote_quality & has_k2 & k2_eukaryote_quality).then(pl.lit("EUKARYOTE_eupath"))
                  .when(k1_classified & is_eukaryote & k1_eukaryote_quality & ~has_k2 & ~has_k3).then(pl.lit("EUKARYOTE_PlusPFOnly"))
                  # ── K3-only helminth (K1 unclassified, tier2) ────────────────
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

    return (classify_read,)


@app.cell
def _():
    BASE = "/mnt/san/microbio/metagenomique_clinique/noscendo_benchmarking"

    def _paths(folder, sample):
        return {
            "p_k1":  f"{BASE}/{folder}/output/kraken2/{sample}.krk",
            "p_k2":  f"{BASE}/{folder}/output/kraken2_eupathdb/{sample}.krk",
            "p_k3":  f"{BASE}/{folder}/output/kraken2_helminth/{sample}.krk",  # ← K3
            "p_bam": f"{BASE}/{folder}/output/bam/{sample}.bam",
        }

    COHORT = [
        {"sample_id": "LRG_r13",   "folder": "240911_NB501647_0270_AHNJNGBGXW",   "sample": "A1001", "dna_ng": 1.16, "is_blank": False, "run": 13},
        {"sample_id": "AM_r13",    "folder": "240911_NB501647_0270_AHNJNGBGXW",   "sample": "A1002", "dna_ng": 0.64, "is_blank": False, "run": 13},
        {"sample_id": "MRV_r13",   "folder": "240911_NB501647_0270_AHNJNGBGXW",   "sample": "A1003", "dna_ng": 0.20, "is_blank": False, "run": 13},
        {"sample_id": "KJM_r13",   "folder": "240911_NB501647_0270_AHNJNGBGXW",   "sample": "A1004", "dna_ng": 0.20, "is_blank": False, "run": 13},
        {"sample_id": "GR_r13",    "folder": "240911_NB501647_0270_AHNJNGBGXW",   "sample": "A1005", "dna_ng": 0.17, "is_blank": False, "run": 13},
        {"sample_id": "CP_r13",    "folder": "240911_NB501647_0270_AHNJNGBGXW",   "sample": "A1006", "dna_ng": 0.09, "is_blank": False, "run": 13},
        {"sample_id": "LG_r13",    "folder": "240911_NB501647_0270_AHNJNGBGXW",   "sample": "A1007", "dna_ng": 8.85, "is_blank": False, "run": 13},
        {"sample_id": "BLANK_r13", "folder": "240911_NB501647_0270_AHNJNGBGXW",   "sample": "A1008", "dna_ng": 0.94, "is_blank": True,  "run": 13},
        {"sample_id": "CL_r15",    "folder": "241024_NB501647_0277_AHNHGMBGXW",   "sample": "A1001", "dna_ng": 0.13, "is_blank": False, "run": 15},
        {"sample_id": "EHM_r15",   "folder": "241024_NB501647_0277_AHNHGMBGXW",   "sample": "A1002", "dna_ng": 0.22, "is_blank": False, "run": 15},
        {"sample_id": "OS_r15",    "folder": "241024_NB501647_0277_AHNHGMBGXW",   "sample": "A1003", "dna_ng": 0.12, "is_blank": False, "run": 15},
        {"sample_id": "GM_r15",    "folder": "241024_NB501647_0277_AHNHGMBGXW",   "sample": "A1004", "dna_ng": 0.80, "is_blank": False, "run": 15},
        {"sample_id": "LP_r15",    "folder": "241024_NB501647_0277_AHNHGMBGXW",   "sample": "A1005", "dna_ng": 0.26, "is_blank": False, "run": 15},
        {"sample_id": "BLANK_r15", "folder": "241024_NB501647_0277_AHNHGMBGXW",   "sample": "A1006", "dna_ng": 0.94, "is_blank": True,  "run": 15},
        {"sample_id": "CL_r16",    "folder": "241128_NB501647_0281_AHG7G7BGYW",   "sample": "A1001", "dna_ng": 0.10, "is_blank": False, "run": 16},
        {"sample_id": "OS_r16",    "folder": "241128_NB501647_0281_AHG7G7BGYW",   "sample": "A1002", "dna_ng": 0.11, "is_blank": False, "run": 16},
        {"sample_id": "LMP_r16",   "folder": "241128_NB501647_0281_AHG7G7BGYW",   "sample": "A1003", "dna_ng": 1.34, "is_blank": False, "run": 16},
        {"sample_id": "CF_r16",    "folder": "241128_NB501647_0281_AHG7G7BGYW",   "sample": "A1004", "dna_ng": 1.02, "is_blank": False, "run": 16},
        {"sample_id": "BLANK_r16", "folder": "241128_NB501647_0281_AHG7G7BGYW",   "sample": "A1005", "dna_ng": 0.88, "is_blank": True,  "run": 16},
        {"sample_id": "GP_r17",    "folder": "241205_NB501647_0283_AHYKTVBGXW",   "sample": "A1001", "dna_ng": 0.59, "is_blank": False, "run": 17},
        {"sample_id": "CE_r17",    "folder": "241205_NB501647_0283_AHYKTVBGXW",   "sample": "A1002", "dna_ng": 0.11, "is_blank": False, "run": 17},
        {"sample_id": "TM_r17",    "folder": "241205_NB501647_0283_AHYKTVBGXW",   "sample": "A1003", "dna_ng": 0.16, "is_blank": False, "run": 17},
        {"sample_id": "SA_r17",    "folder": "241205_NB501647_0283_AHYKTVBGXW",   "sample": "A1004", "dna_ng": 1.13, "is_blank": False, "run": 17},
        {"sample_id": "LF_r17",    "folder": "241205_NB501647_0283_AHYKTVBGXW",   "sample": "A1005", "dna_ng": 0.12, "is_blank": False, "run": 17},
        {"sample_id": "BLANK_r17", "folder": "241205_NB501647_0283_AHYKTVBGXW",   "sample": "A1006", "dna_ng": 1.00, "is_blank": True,  "run": 17},
        {"sample_id": "CE_r18",    "folder": "241212_NB501647_0284_AHYKJKBGXW",   "sample": "A1001", "dna_ng": 0.12, "is_blank": False, "run": 18},
        {"sample_id": "LF_r18",    "folder": "241212_NB501647_0284_AHYKJKBGXW",   "sample": "A1002", "dna_ng": 0.09, "is_blank": False, "run": 18},
        {"sample_id": "PA_r18",    "folder": "241212_NB501647_0284_AHYKJKBGXW",   "sample": "A1003", "dna_ng": 0.13, "is_blank": False, "run": 18},
        {"sample_id": "EK_r18",    "folder": "241212_NB501647_0284_AHYKJKBGXW",   "sample": "A1004", "dna_ng": 2.64, "is_blank": False, "run": 18},
        {"sample_id": "BLANK_r18", "folder": "241212_NB501647_0284_AHYKJKBGXW",   "sample": "A1005", "dna_ng": 1.00, "is_blank": True,  "run": 18},
        {"sample_id": "KG_r20",    "folder": "250108_NB501647_0288_AHJGJVBGYW",   "sample": "A1001", "dna_ng": 0.99, "is_blank": False, "run": 20},
        {"sample_id": "HR_r20",    "folder": "250108_NB501647_0288_AHJGJVBGYW",   "sample": "A1002", "dna_ng": 1.70, "is_blank": False, "run": 20},
        {"sample_id": "LGJ_r20",   "folder": "250108_NB501647_0288_AHJGJVBGYW",   "sample": "A1003", "dna_ng": 0.50, "is_blank": False, "run": 20},
        {"sample_id": "DI_r20",    "folder": "250108_NB501647_0288_AHJGJVBGYW",   "sample": "A1004", "dna_ng": 0.20, "is_blank": False, "run": 20},
        {"sample_id": "BLANK_r20", "folder": "250108_NB501647_0288_AHJGJVBGYW",   "sample": "A1005", "dna_ng": 0.98, "is_blank": True,  "run": 20},
        {"sample_id": "KG_r21",    "folder": "250213_NB501647_0293_AH32GMBGYW",   "sample": "A1001", "dna_ng": 0.95, "is_blank": False, "run": 21},
        {"sample_id": "EMB_r21",   "folder": "250213_NB501647_0293_AH32GMBGYW",   "sample": "A1002", "dna_ng": 1.16, "is_blank": False, "run": 21},
        {"sample_id": "LA_r21",    "folder": "250213_NB501647_0293_AH32GMBGYW",   "sample": "A1003", "dna_ng": 0.26, "is_blank": False, "run": 21},
        {"sample_id": "SJC_r21",   "folder": "250213_NB501647_0293_AH32GMBGYW",   "sample": "A1004", "dna_ng": 0.11, "is_blank": False, "run": 21},
        {"sample_id": "BLANK_r21", "folder": "250213_NB501647_0293_AH32GMBGYW",   "sample": "A1005", "dna_ng": 0.92, "is_blank": True,  "run": 21},
        {"sample_id": "SJC_r22",   "folder": "250220_NB501647_0294_AHH7VFBGYW",   "sample": "A1001", "dna_ng": 0.08, "is_blank": False, "run": 22},
        {"sample_id": "SS_r22",    "folder": "250220_NB501647_0294_AHH7VFBGYW",   "sample": "A1002", "dna_ng": 5.65, "is_blank": False, "run": 22},
        {"sample_id": "PG_r22",    "folder": "250220_NB501647_0294_AHH7VFBGYW",   "sample": "A1003", "dna_ng": 0.13, "is_blank": False, "run": 22},
        {"sample_id": "DV_r22",    "folder": "250220_NB501647_0294_AHH7VFBGYW",   "sample": "A1004", "dna_ng": 0.88, "is_blank": False, "run": 22},
        {"sample_id": "RM_r22",    "folder": "250220_NB501647_0294_AHH7VFBGYW",   "sample": "A1005", "dna_ng": 6.13, "is_blank": False, "run": 22},
        {"sample_id": "BLANK_r22", "folder": "250220_NB501647_0294_AHH7VFBGYW",   "sample": "A1006", "dna_ng": 0.92, "is_blank": True,  "run": 22},
        {"sample_id": "SS_r23",    "folder": "250304_NB501647_0296_AHJJFCBGYW",   "sample": "A1001", "dna_ng": 6.48, "is_blank": False, "run": 23},
        {"sample_id": "PG_r23",    "folder": "250304_NB501647_0296_AHJJFCBGYW",   "sample": "A1002", "dna_ng": 0.14, "is_blank": False, "run": 23},
        {"sample_id": "MK_r23",    "folder": "250304_NB501647_0296_AHJJFCBGYW",   "sample": "A1003", "dna_ng": 0.08, "is_blank": False, "run": 23},
        {"sample_id": "LA_r23",    "folder": "250304_NB501647_0296_AHJJFCBGYW",   "sample": "A1004", "dna_ng": 6.08, "is_blank": False, "run": 23},
        {"sample_id": "RF_r23",    "folder": "250304_NB501647_0296_AHJJFCBGYW",   "sample": "A1005", "dna_ng": 3.41, "is_blank": False, "run": 23},
        {"sample_id": "PF_r23",    "folder": "250304_NB501647_0296_AHJJFCBGYW",   "sample": "A1006", "dna_ng": 1.06, "is_blank": False, "run": 23},
        {"sample_id": "BLANK_r23", "folder": "250304_NB501647_0296_AHJJFCBGYW",   "sample": "A1007", "dna_ng": 0.88, "is_blank": True,  "run": 23},
        {"sample_id": "HD_r24",    "folder": "250326_NB501647_0299_AHTY5MBGYW",   "sample": "A1001", "dna_ng": 0.16, "is_blank": False, "run": 24},
        {"sample_id": "PM_r24",    "folder": "250326_NB501647_0299_AHTY5MBGYW",   "sample": "A1002", "dna_ng": 0.25, "is_blank": False, "run": 24},
        {"sample_id": "GT_r24",    "folder": "250326_NB501647_0299_AHTY5MBGYW",   "sample": "A1003", "dna_ng": 0.26, "is_blank": False, "run": 24},
        {"sample_id": "GL_r24",    "folder": "250326_NB501647_0299_AHTY5MBGYW",   "sample": "A1004", "dna_ng": 0.08, "is_blank": False, "run": 24},
        {"sample_id": "DJ_r24",    "folder": "250326_NB501647_0299_AHTY5MBGYW",   "sample": "A1005", "dna_ng": 0.52, "is_blank": False, "run": 24},
        {"sample_id": "LG_r24",    "folder": "250326_NB501647_0299_AHTY5MBGYW",   "sample": "A1006", "dna_ng": 0.05, "is_blank": False, "run": 24},
        {"sample_id": "SS_r24",    "folder": "250326_NB501647_0299_AHTY5MBGYW",   "sample": "A1007", "dna_ng": 6.22, "is_blank": False, "run": 24},
        {"sample_id": "BLANK_r24", "folder": "250326_NB501647_0299_AHTY5MBGYW",   "sample": "A1008", "dna_ng": 0.90, "is_blank": True,  "run": 24},
        {"sample_id": "HD_r25",    "folder": "250327_NB501647_0300_AHJGV3BGYW_1", "sample": "A1001", "dna_ng": 0.12, "is_blank": False, "run": 25},
        {"sample_id": "PM_r25",    "folder": "250327_NB501647_0300_AHJGV3BGYW_1", "sample": "A1002", "dna_ng": 0.23, "is_blank": False, "run": 25},
        {"sample_id": "GT_r25",    "folder": "250327_NB501647_0300_AHJGV3BGYW_1", "sample": "A1003", "dna_ng": 0.21, "is_blank": False, "run": 25},
        {"sample_id": "DJ_r25",    "folder": "250327_NB501647_0300_AHJGV3BGYW_1", "sample": "A1004", "dna_ng": 0.53, "is_blank": False, "run": 25},
        {"sample_id": "LG_r25",    "folder": "250327_NB501647_0300_AHJGV3BGYW_1", "sample": "A1005", "dna_ng": 0.05, "is_blank": False, "run": 25},
        {"sample_id": "SS_r25",    "folder": "250327_NB501647_0300_AHJGV3BGYW_1", "sample": "A1006", "dna_ng": 5.77, "is_blank": False, "run": 25},
        {"sample_id": "LGG_r25",   "folder": "250327_NB501647_0300_AHJGV3BGYW_1", "sample": "A1007", "dna_ng": 0.30, "is_blank": False, "run": 25},
        {"sample_id": "BLANK_r25", "folder": "250327_NB501647_0300_AHJGV3BGYW_1", "sample": "A1008", "dna_ng": 0.88, "is_blank": True,  "run": 25},
        {"sample_id": "LBR_r26",   "folder": "250411_NB501647_0303_AHWGT3BGYW",   "sample": "A1001", "dna_ng": 0.53, "is_blank": False, "run": 26},
        {"sample_id": "RA_r26",    "folder": "250411_NB501647_0303_AHWGT3BGYW",   "sample": "A1002", "dna_ng": 0.21, "is_blank": False, "run": 26},
        {"sample_id": "SL_r26",    "folder": "250411_NB501647_0303_AHWGT3BGYW",   "sample": "A1003", "dna_ng": 0.09, "is_blank": False, "run": 26},
        {"sample_id": "RN_r26",    "folder": "250411_NB501647_0303_AHWGT3BGYW",   "sample": "A1004", "dna_ng": 0.33, "is_blank": False, "run": 26},
        {"sample_id": "BLANK_r26", "folder": "250411_NB501647_0303_AHWGT3BGYW",   "sample": "A1005", "dna_ng": 0.96, "is_blank": True,  "run": 26},
        {"sample_id": "DB_r27",    "folder": "250423_NB501647_0304_AHWGTLBGYW",   "sample": "A1001", "dna_ng": 0.30, "is_blank": False, "run": 27},
        {"sample_id": "PR_r27",    "folder": "250423_NB501647_0304_AHWGTLBGYW",   "sample": "A1002", "dna_ng": 0.39, "is_blank": False, "run": 27},
        {"sample_id": "BR_r27",    "folder": "250423_NB501647_0304_AHWGTLBGYW",   "sample": "A1003", "dna_ng": 2.33, "is_blank": False, "run": 27},
        {"sample_id": "LGM_r27",   "folder": "250423_NB501647_0304_AHWGTLBGYW",   "sample": "A1004", "dna_ng": 0.07, "is_blank": False, "run": 27},
        {"sample_id": "CT_r27",    "folder": "250423_NB501647_0304_AHWGTLBGYW",   "sample": "A1005", "dna_ng": 0.47, "is_blank": False, "run": 27},
        {"sample_id": "PF_r27",    "folder": "250423_NB501647_0304_AHWGTLBGYW",   "sample": "A1006", "dna_ng": 1.55, "is_blank": False, "run": 27},
        {"sample_id": "BLANK_r27", "folder": "250423_NB501647_0304_AHWGTLBGYW",   "sample": "A1007", "dna_ng": 1.02, "is_blank": True,  "run": 27},
        {"sample_id": "PJP_r28",   "folder": "250522_NB501647_0309_AH2VL2BGYX",   "sample": "A1001", "dna_ng": 0.67, "is_blank": False, "run": 28},
        {"sample_id": "RJ_r28",    "folder": "250522_NB501647_0309_AH2VL2BGYX",   "sample": "A1002", "dna_ng": 0.55, "is_blank": False, "run": 28},
        {"sample_id": "CJM_r28",   "folder": "250522_NB501647_0309_AH2VL2BGYX",   "sample": "A1003", "dna_ng": 0.34, "is_blank": False, "run": 28},
        {"sample_id": "BLANK_r28", "folder": "250522_NB501647_0309_AH2VL2BGYX",   "sample": "A1004", "dna_ng": 0.91, "is_blank": True,  "run": 28},
        {"sample_id": "CJ_r29",    "folder": "250528_NB501647_0311_AH2VVFBGYX",   "sample": "A1001", "dna_ng": 0.06, "is_blank": False, "run": 29},
        {"sample_id": "KD_r29",    "folder": "250528_NB501647_0311_AH2VVFBGYX",   "sample": "A1002", "dna_ng": 0.09, "is_blank": False, "run": 29},
        {"sample_id": "LHG_r29",   "folder": "250528_NB501647_0311_AH2VVFBGYX",   "sample": "A1003", "dna_ng": 1.54, "is_blank": False, "run": 29},
        {"sample_id": "PJP2_r29",  "folder": "250528_NB501647_0311_AH2VVFBGYX",   "sample": "A1004", "dna_ng": 0.22, "is_blank": False, "run": 29},
        {"sample_id": "BLANK_r29", "folder": "250528_NB501647_0311_AH2VVFBGYX",   "sample": "A1005", "dna_ng": 0.28, "is_blank": True,  "run": 29},
        {"sample_id": "LBP_r30",   "folder": "250619_NB501647_0314_AH2TJGBGYX",   "sample": "A1001", "dna_ng": 0.29, "is_blank": False, "run": 30},
        {"sample_id": "GJ_r30",    "folder": "250619_NB501647_0314_AH2TJGBGYX",   "sample": "A1002", "dna_ng": 0.16, "is_blank": False, "run": 30},
        {"sample_id": "KD_r30",    "folder": "250619_NB501647_0314_AH2TJGBGYX",   "sample": "A1003", "dna_ng": 0.08, "is_blank": False, "run": 30},
        {"sample_id": "BL_r30",    "folder": "250619_NB501647_0314_AH2TJGBGYX",   "sample": "A1004", "dna_ng": 0.24, "is_blank": False, "run": 30},
        {"sample_id": "SAB_r30",   "folder": "250619_NB501647_0314_AH2TJGBGYX",   "sample": "A1005", "dna_ng": 0.26, "is_blank": False, "run": 30},
        {"sample_id": "BLANK_r30", "folder": "250619_NB501647_0314_AH2TJGBGYX",   "sample": "A1006", "dna_ng": 0.29, "is_blank": True,  "run": 30},
        {"sample_id": "LD_r31",    "folder": "250717_NB501647_0319_AH37VHBGYX",   "sample": "A1001", "dna_ng": 0.33, "is_blank": False, "run": 31},
        {"sample_id": "HA_r31",    "folder": "250717_NB501647_0319_AH37VHBGYX",   "sample": "A1002", "dna_ng": 1.52, "is_blank": False, "run": 31},
        {"sample_id": "AM_r31",    "folder": "250717_NB501647_0319_AH37VHBGYX",   "sample": "A1003", "dna_ng": 0.57, "is_blank": False, "run": 31},
        {"sample_id": "DJ_r31",    "folder": "250717_NB501647_0319_AH37VHBGYX",   "sample": "A1004", "dna_ng": 0.16, "is_blank": False, "run": 31},
        {"sample_id": "BLANK_r31", "folder": "250717_NB501647_0319_AH37VHBGYX",   "sample": "A1005", "dna_ng": 0.28, "is_blank": True,  "run": 31},
        {"sample_id": "EB_r32",    "folder": "250811_NB501647_0324_AH7C7HBGYX",   "sample": "A1001", "dna_ng": 0.17, "is_blank": False, "run": 32},
        {"sample_id": "LGY_r32",   "folder": "250811_NB501647_0324_AH7C7HBGYX",   "sample": "A1002", "dna_ng": 1.25, "is_blank": False, "run": 32},
        {"sample_id": "GF_r32",    "folder": "250811_NB501647_0324_AH7C7HBGYX",   "sample": "A1003", "dna_ng": 2.93, "is_blank": False, "run": 32},
        {"sample_id": "BL_r32",    "folder": "250811_NB501647_0324_AH7C7HBGYX",   "sample": "A1004", "dna_ng": 0.73, "is_blank": False, "run": 32},
        {"sample_id": "MW_r32",    "folder": "250811_NB501647_0324_AH7C7HBGYX",   "sample": "A1005", "dna_ng": 0.22, "is_blank": False, "run": 32},
        {"sample_id": "RA_r32",    "folder": "250811_NB501647_0324_AH7C7HBGYX",   "sample": "A1006", "dna_ng": 3.09, "is_blank": False, "run": 32},
        {"sample_id": "BLANK_r32", "folder": "250811_NB501647_0324_AH7C7HBGYX",   "sample": "A1007", "dna_ng": 0.29, "is_blank": True,  "run": 32},
        {"sample_id": "P1_r33",    "folder": "250905_NB501647_0328_AH7C52BGYX",   "sample": "A1001", "dna_ng": 1.45, "is_blank": False, "run": 33},
        {"sample_id": "P2_r33",    "folder": "250905_NB501647_0328_AH7C52BGYX",   "sample": "A1002", "dna_ng": 0.82, "is_blank": False, "run": 33},
        {"sample_id": "P3_r33",    "folder": "250905_NB501647_0328_AH7C52BGYX",   "sample": "A1003", "dna_ng": 0.14, "is_blank": False, "run": 33},
        {"sample_id": "BLANK_r33", "folder": "250905_NB501647_0328_AH7C52BGYX",   "sample": "A1004", "dna_ng": 0.24, "is_blank": True,  "run": 33},
    ]
    for _s in COHORT:
        _s.update(_paths(_s["folder"], _s["sample"]))

    BLANKS = [s for s in COHORT if s["is_blank"]]
    print(f"Cohort: {len(COHORT)} samples total | {len(BLANKS)} blanks")
    return (BLANKS,)


@app.cell
def _(P_EUPATH_INSPECT, P_KRAKEN2_NAMES, P_KRAKEN2_NODES, P_TAXMAP, pl):
    # ── nodes.dmp → parent + rank ─────────────────────────────────────────────
    parent_g = {}
    rank_g   = {}
    with open(P_KRAKEN2_NODES) as _fh:
        for _line in _fh:
            _parts = _line.split("\t|\t")
            if len(_parts) < 3:
                continue
            try:
                _t = int(_parts[0].strip())
                parent_g[_t] = int(_parts[1].strip())
                rank_g[_t]   = _parts[2].strip()
            except ValueError:
                continue
    print(f"  nodes.dmp  : {len(parent_g):,} nodes")

    # ── names.dmp → taxid → scientific name ──────────────────────────────────
    name_lookup_g = {}
    with open(P_KRAKEN2_NAMES) as _fh:
        for _line in _fh:
            _parts = _line.split("\t|\t")
            if len(_parts) < 4:
                continue
            try:
                _t = int(_parts[0].strip())
            except ValueError:
                continue
            if _parts[3].strip().rstrip("\t|").strip() == "scientific name":
                name_lookup_g[_t] = _parts[1].strip()
    print(f"  names.dmp  : {len(name_lookup_g):,} scientific names")

    # ── EuPathDB inspect.txt → k2_taxid → name ───────────────────────────────
    _eupath_rows = []
    with open(P_EUPATH_INSPECT) as _fh:
        for _line in _fh:
            _parts = _line.rstrip("\n").split("\t")
            if len(_parts) < 6:
                continue
            try:
                _t = int(_parts[4].strip())
            except ValueError:
                continue
            _eupath_rows.append({"k2_taxid": _t, "k2_sci_name": _parts[5].strip()})
    k2_name_df_g = (
        pl.DataFrame(_eupath_rows, schema={"k2_taxid": pl.Int64, "k2_sci_name": pl.Utf8})
        .unique(subset=["k2_taxid"], keep="first")
    )
    print(f"  EuPathDB   : {len(k2_name_df_g):,} taxa")

    # ── taxmap ────────────────────────────────────────────────────────────────
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
    return k2_name_df_g, name_lookup_g, parent_g, rank_g, taxmap_g


@app.cell
def _(
    ALEN_MIN,
    BREADTH_MIN_FRACTION,
    BREADTH_MIN_WINDOWS,
    GENUS_LEVEL_COLLAPSE,
    MAPQ_FLOOR,
    MAPQ_MIN,
    MATE_SUFFIX_RE,
    MIN_DEPTH,
    NM_RATE_MAX,
    NONREPORTABLE_RANKS,
    Path,
    WINDOW_SIZE,
    classify_read,
    k2_name_df_g,
    load_kraken2_per_read,
    name_lookup_g,
    parent_g,
    pl,
    pysam,
    rank_g,
    taxmap_g,
):
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

    HUMAN_CALLS = frozenset([
        "HUMAN_CLASSIFIED", "HUMAN_CONFLICT",
        "HUMAN_KMER_BROAD_COVERAGE", "LOW_CONF_HUMAN_KMER",
    ])
    UNINFORMATIVE_DISPLAY_TAXIDS = frozenset({
        2759, 2, 10239, 2157, 131567, 1, 33154, 4751,
    })

    def run_blank(s):
        sid = s["sample_id"]

        # ── Step 1: load K1 + K2 + K3 (K3 optional) ──────────────────────────
        k1 = load_kraken2_per_read(s["p_k1"], db_tag="pluspf")
        k2 = load_kraken2_per_read(s["p_k2"], db_tag="eupath")
        k3 = load_kraken2_per_read(s["p_k3"], db_tag="helminths") \
             if Path(s["p_k3"]).exists() else None

        # ── Step 2: kingdom resolution (K1 + K2 + K3 taxids) ─────────────────
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

        # ── Step 3: BAM scan + breadth + taxmap ───────────────────────────────
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

        # ── Step 4: merge K1 + K2 + K3 + BAM ─────────────────────────────────
        def _prefix(df, tag):
            return df.rename({c: f"{tag}_{c}" for c in df.columns if c != "read_id"})

        k1p = _prefix(k1, "k1")
        k2p = _prefix(k2, "k2")
        k3p = _prefix(k3, "k3") if k3 is not None else None

        # Path A: K1-classified
        merged_k1 = k1p.filter(pl.col("k1_status") == "C").join(k2p, on="read_id", how="left")
        if k3p is not None:
            merged_k1 = merged_k1.join(k3p, on="read_id", how="left")
        merged_k1 = merged_k1.join(bam_df, on="read_id", how="left")

        # Path B: BAM-only
        merged_bam_only = (
            bam_df.join(merged_k1.select("read_id"), on="read_id", how="anti")
            .join(k1p, on="read_id", how="left")
            .join(k2p, on="read_id", how="left")
        )
        if k3p is not None:
            merged_bam_only = merged_bam_only.join(k3p, on="read_id", how="left")

        # Path C: K3-only (helminth, K1 missed, not in BAM)
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

        # ── Step 5: classify_read ─────────────────────────────────────────────
        passed = classify_read(merged, taxid_kingdom_df)

        # ── Step 6: human filter ──────────────────────────────────────────────
        nonhuman = passed.filter(~pl.col("read_call").is_in(list(HUMAN_CALLS)))

        # ── Step 7: name resolution (K1 + K3 from pre-loaded dict) ───────────
        _k1_taxids = set(nonhuman["k1_taxid"].drop_nulls().unique().to_list())
        taxid_name_df = pl.DataFrame(
            [{"k1_taxid": t, "k1_sci_name": name_lookup_g.get(t, str(t))} for t in _k1_taxids],
            schema={"k1_taxid": pl.Int64, "k1_sci_name": pl.Utf8},
        )
        # K3 names from same pre-loaded dict
        _k3_name_df = None
        if "k3_taxid" in nonhuman.columns:
            _k3_taxids = set(nonhuman["k3_taxid"].drop_nulls().unique().to_list())
            _k3_name_df = pl.DataFrame(
                [{"k3_taxid": t, "k3_sci_name": name_lookup_g.get(t, str(t))} for t in _k3_taxids],
                schema={"k3_taxid": pl.Int64, "k3_sci_name": pl.Utf8},
            )

        # ── Step 8: summary aggregation ───────────────────────────────────────
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
                  .when(pl.col("read_call").is_in(["EUKARYOTE_helminth",
                                                   "EUKARYOTE_K3_HELMINTH"]))
                  .then(pl.col("k3_taxid").cast(pl.Int64))
                  .when(pl.col("read_call") == "BAM_ONLY")
                  .then(pl.col("bam_taxid_species"))
                  .otherwise(pl.col("k1_taxid"))
                  .alias("display_taxid"),
                pl.when(pl.col("read_call") == "EUKARYOTE_eupath")
                  .then(pl.col("k2_sci_name"))
                  .when(pl.col("read_call").is_in(["EUKARYOTE_helminth",
                                                   "EUKARYOTE_K3_HELMINTH"]))
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

        # ── Step 9: summarise_sample ───────────────────────────────────────────
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
                # tier1: K1+K2, K1+K3 double confirmation
                pl.when(pl.col("read_call").is_in([
                    "KRAKEN_CONFIRMED", "EUKARYOTE_eupath",
                    "EUKARYOTE_helminth",          # ← K3 tier1
                ])).then(pl.col("n_reads")).otherwise(0).sum().alias("n_reads_tier1"),
                # tier2: single-source + K3 helminth
                pl.when(pl.col("read_call").is_in([
                    "BAM_K1_LOW_CONF", "BACTERIA_K1_ONLY",
                    "EUKARYOTE_K1_ONLY", "VIRUS_K1_ONLY",
                    "BAM_ONLY", "EUKARYOTE_K3_HELMINTH",  # ← K3 tier2
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
                pl.lit(s["dna_ng"]).alias("dna_input_ng"),
                pl.lit(s["is_blank"]).alias("is_blank"),   # ← fixed (was hardcoded True)
                pl.lit(s["run"]).cast(pl.Int32).alias("run"),
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
                  .group_by(["kingdom", "group_taxid", "group_name",
                             "sample_id", "dna_input_ng", "is_blank", "run"])
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

    return (run_blank,)


@app.cell
def _(BLANKS, P_OUTPUT_DIR, pl, run_blank):
    _results = []
    _errors  = []

    for _i, _s in enumerate(BLANKS):
        print(f"[{_i+1:>2}/{len(BLANKS)}] {_s['sample_id']}  "
              f"run={_s['run']}  {_s['dna_ng']} ng", end="  →  ")
        try:
            _row = run_blank(_s)
            _results.append(_row)
            print(f"OK  ({len(_row)} rows)")
            pl.concat(_results, how="diagonal_relaxed").write_parquet(
                f"{P_OUTPUT_DIR}/blank_cohort_raw.parquet"
            )
        except Exception as _e:
            _errors.append({"sample_id": _s["sample_id"], "error": str(_e)})
            print(f"ERROR  {_e}")

    blank_cohort = pl.concat(_results, how="diagonal_relaxed")

    print(f"\n── Summary ─────────────────────────────────────────────────────")
    print(f"Blanks processed  : {len(_results)}/{len(BLANKS)}")
    print(f"Errors            : {len(_errors)}")
    print(f"Total rows        : {len(blank_cohort):,}")
    print(f"Unique organisms  : {blank_cohort['group_taxid'].n_unique():,}")
    if _errors:
        print("\nFailed samples:")
        for _e in _errors:
            print(f"  {_e['sample_id']}: {_e['error']}")

    blank_cohort
    return (blank_cohort,)


@app.cell
def _(blank_cohort, pl):
    blank_cohort.filter(pl.col("kingdom") == "Viruses")["group_name"].sort().unique().sort()
    return


@app.cell
def _(BLANKS, HUMAN_TAXIDS, blank_cohort, pl):
    def _count_dehosted_reads_k1(path: str) -> int:
        """
        Count non-human reads in a Kraken2 output file.
        Uses the globally defined HUMAN_TAXIDS for consistency.
        Returns the effective dehosted library size for RPM normalisation.
        """
        n = 0
        with open(path) as fh:
            for line in fh:
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                try:
                    taxid = int(parts[2])
                except ValueError:
                    continue
                if taxid not in HUMAN_TAXIDS:
                    n += 1
        return n

    lib_sizes = []
    for s in BLANKS:
        total = _count_dehosted_reads_k1(s["p_k1"])
        lib_sizes.append({
            "sample_id":   s["sample_id"],
            "run":         s["run"],
            "dna_ng":      s["dna_ng"],
            "total_reads": total,
        })

    lib_df = pl.DataFrame(lib_sizes).sort("run")
    blank_cohort_norm = (
        blank_cohort
        .join(lib_df.select(["sample_id", "total_reads"]), on="sample_id", how="left")
        .with_columns(
            (pl.col("n_reads_total") / pl.col("total_reads") * 1_000_000).alias("rpm")
        )
    )
    return (blank_cohort_norm,)


@app.cell
def _(blank_cohort, blank_cohort_norm, pl):
    _N = blank_cohort["sample_id"].n_unique()

    _run28_reads = (
        blank_cohort.filter(pl.col("sample_id") == "BLANK_r28")
        .select(["group_taxid", pl.col("n_reads_total").alias("reads_run28")])
    )

    contaminant_db = (
        blank_cohort
        .group_by([
            "group_taxid", "group_name", "kingdom",
            "tax_phylum", "tax_class", "tax_order",
            "tax_family", "tax_genus", "tax_species", "taxonomy",
        ])
        .agg([
            pl.col("sample_id").n_unique().alias("n_blanks_present"),
            pl.col("n_reads_total").mean().alias("mean_reads_blank"),
            pl.col("n_reads_total").max().alias("max_reads_blank"),
            pl.col("n_reads_total").sum().alias("sum_reads_all_blanks"),
            pl.col("best_tier").min().alias("min_tier_blank"),
            (pl.col("run") == 28).any().alias("in_run28_blank"),
            pl.col("run").sort().alias("runs_present"),
        ])
        .with_columns([
            (pl.col("n_blanks_present") / _N).alias("prevalence"),
            (pl.col("sum_reads_all_blanks") / _N).alias("mean_reads_all_blanks"),
        ])
        .join(_run28_reads, on="group_taxid", how="left")
        .with_columns(pl.col("reads_run28").fill_null(0))
        .with_columns(
            pl.when(pl.col("prevalence") >= 0.75)
              .then(pl.lit("UNIVERSAL_KIT"))
              .when(
                  pl.col("in_run28_blank") &
                  (pl.col("reads_run28") > (pl.col("mean_reads_all_blanks") * 10))
              ).then(pl.lit("RUN_BATCH"))
              .when(pl.col("n_blanks_present") >= 3)
              .then(pl.lit("REAGENT"))
              .otherwise(pl.lit("SPORADIC"))
              .alias("contamination_class")
        )
        .sort(["contamination_class", "n_blanks_present", "mean_reads_blank"],
              descending=[False, True, True])
    )

    contaminant_db_final = (
        contaminant_db
        .join(
            blank_cohort_norm
            .group_by(["group_taxid"])
            .agg([
                pl.col("rpm").mean().alias("mean_rpm_blank"),
                pl.col("rpm").max().alias("max_rpm_blank"),
                pl.col("rpm").std().alias("std_rpm_blank"),
            ]),
            on="group_taxid", how="left"
        )
        .with_columns(
            pl.when(pl.col("group_taxid").is_in([1758194, 40324, 40323]))
              .then(pl.lit("27,28,29,30"))
              .when(pl.col("group_taxid").is_in([488447, 32008, 1822464]))
              .then(pl.lit("29,30,31,32,33"))
              .otherwise(pl.lit(""))
              .alias("lot_batch_runs")
        )
    )
    contaminant_db_final
    return (contaminant_db_final,)


@app.cell(hide_code=True)
def _(blank_cohort):
    run_blank_ref = (
        blank_cohort
        .select([
            "run", "sample_id", "group_taxid", "group_name",
            "n_reads_total", "n_reads_tier1", "n_reads_tier2",
            "n_reads_tier3", "best_tier",
        ])
        .rename({
            "n_reads_total": "blank_reads_total",
            "n_reads_tier1": "blank_reads_tier1",
            "n_reads_tier2": "blank_reads_tier2",
            "n_reads_tier3": "blank_reads_tier3",
            "best_tier":     "blank_best_tier",
        })
    )
    return (run_blank_ref,)


@app.cell
def _(run_blank_ref):
    run_blank_ref
    return


@app.cell
def _(P_OUTPUT_DIR, blank_cohort, contaminant_db_final, pl, run_blank_ref):
    _p_raw    = f"{P_OUTPUT_DIR}/blank_cohort_raw.parquet"
    _p_db_pq  = f"{P_OUTPUT_DIR}/contaminant_db.parquet"
    _p_db_tsv = f"{P_OUTPUT_DIR}/contaminant_db.tsv"
    _p_ref    = f"{P_OUTPUT_DIR}/run_blank_ref.parquet"

    blank_cohort.write_parquet(_p_raw)
    contaminant_db_final.write_parquet(_p_db_pq)
    (
        contaminant_db_final
        .with_columns(
            pl.col("runs_present")
              .list.eval(pl.element().cast(pl.String))
              .list.join(",")
              .alias("runs_present")
        )
        .write_csv(_p_db_tsv, separator="\t")
    )
    run_blank_ref.write_parquet(_p_ref)

    print("── Outputs saved ───────────────────────────────────────────────")
    print(f"  {_p_raw}")
    print(f"  {_p_db_pq}")
    print(f"  {_p_db_tsv}")
    print(f"  {_p_ref}")
    print(f"\nTo load in your main notebook:")
    print(f'  contaminant_db = pl.read_parquet("{_p_db_pq}")')
    print(f'  run_blank_ref  = pl.read_parquet("{_p_ref}")')
    return


if __name__ == "__main__":
    app.run()
