# rules/count.smk

# ── helpers ────────────────────────────────────────────────────────────────
# Strip lane suffix from seqkit basename to recover the sample_core
# e.g.  TPOS-ADN_S8_L002_dedup_R1.clean_1.fastq.gz  →  TPOS-ADN
AWK_SUM_BY_CORE = r"""awk 'NR>1 {
    name = $1;
    gsub(/_S[0-9]+_L[0-9]+_dedup_R1\.clean_1\.fastq\.gz$/, "", name);
    sum[name] += $4
} END {
    for (s in sum) print s "," sum[s]
}'"""

AWK_SUM_BY_LANE = r"""awk 'NR>1 {
    name = $1;
    gsub(/_dedup_R1\.fastq\.gz$/, "", name);
    sum[name] += $4
} END {
    for (s in sum) print s "," sum[s]
}'"""


# ── count raw reads (per lane, before QC) ─────────────────────────────────
rule count_raw_reads:
    input:
        r1 = expand(
            "{r1}",                               # raw paths from SAMPLES table
            r1 = SAMPLES["R1"].dropna().tolist()
        )
    output:
        nreads = temp("output/Nreads_raw.txt")
    threads:
        config["threads"]
    shell:
        """
        seqkit stats -j {threads} --basename {input.r1} \
            | awk 'NR>1 {{print $1","$4}}' \
            > {output.nreads}
        """


# ── count deduped reads (per lane, after dedup) ───────────────────────────
rule count_dedup_reads:
    input:
        r1 = expand(
            "output/dedup/{names}_dedup_R1.fastq.gz",
            names = SAMPLES["names"].tolist()
        )
    output:
        nreads = temp("output/Nreads_dedup.txt")
    threads:
        config["threads"]
    params:
        awk = AWK_SUM_BY_LANE
    shell:
        """
        seqkit stats -j {threads} --basename {input.r1} \
            | {params.awk} \
            > {output.nreads}
        """


# ── count dehosted reads (per lane, aggregated by sample_core) ────────────
rule count_dehost_reads:
    input:
        r1 = expand(
            "output/dehost/{names}_dedup_R1.clean_1.fastq.gz",
            names = SAMPLES["names"].tolist()     # <-- same pattern as rule all
        )
    output:
        nreads = temp("output/Nreads_host.txt")
    threads:
        config["threads"]
    params:
        awk = AWK_SUM_BY_CORE
    shell:
        """
        seqkit stats -j {threads} --basename {input.r1} \
            | {params.awk} \
            > {output.nreads}
        """


# ── merge all counts into one CSV ─────────────────────────────────────────
rule merge_counts:
    input:
        raw    = "output/Nreads_raw.txt",
        dedup  = "output/Nreads_dedup.txt",
        dehost = "output/Nreads_host.txt"
    output:
        csv = "output/Nreads.csv"
    run:
        import polars as pl

        def read_counts(path, col_name):
            return (
                pl.read_csv(path, has_header=False,
                            new_columns=["sample", col_name])
            )

        raw    = read_counts(input.raw,    "raw_reads")
        dedup  = read_counts(input.dedup,  "dedup_reads")
        dehost = read_counts(input.dehost, "dehost_reads")

        (
            raw
            .join(dedup,  on="sample", how="left")
            .join(dehost, on="sample", how="left")
            .with_columns([
                (pl.col("dedup_reads")  / pl.col("raw_reads")   * 100).round(1).alias("pct_kept_dedup"),
                (pl.col("dehost_reads") / pl.col("dedup_reads") * 100).round(1).alias("pct_kept_dehost"),
            ])
            .write_csv(output.csv)
        )
