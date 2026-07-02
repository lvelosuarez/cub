# rules/count.smk

AWK_SUM_BY_CORE = r"""awk 'NR>1 {
    name = $1;
    gsub(/_S[0-9]+_L[0-9]+_dedup_R1\.clean_1\.fastq\.gz$/, "", name);
    gsub(/,/, "", $4);
    sum[name] += $4
} END {
    for (s in sum) print s "," sum[s]
}'"""

AWK_SUM_BY_LANE = r"""awk 'NR>1 {
    name = $1;
    gsub(/_S[0-9]+_L[0-9]+_dedup_R1\.fastq\.gz$/, "", name);
    gsub(/,/, "", $4);
    sum[name] += $4
} END {
    for (s in sum) print s "," sum[s]
}'"""

AWK_SUM_RAW = r"""awk 'NR>1 {
    name = $1;
    gsub(/_R1_001\.fastq\.gz$/, "", name);
    gsub(/_R1\.fastq\.gz$/, "", name);
    gsub(/,/, "", $4);
    sum[name] += $4
} END {
    for (s in sum) print s "," sum[s]
}'"""

rule count_raw_reads:
    input:
        r1 = expand(
            "{r1}",
            r1 = SAMPLES["R1"].dropna().tolist()
        )
    output:
        nreads = temp("output/Nreads_raw.txt")
    threads:
        config["threads"]
    params:
        awk = AWK_SUM_RAW
    shell:
        """
        seqkit stats -j {threads} --basename {input.r1} \
            | {params.awk} \
            > {output.nreads}
        """

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

rule count_dehost_reads:
    input:
        r1 = expand(
            "output/dehost/{names}_dedup_R1.clean_1.fastq.gz",
            names = SAMPLES["names"].tolist()
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
            return pl.read_csv(
                path,
                has_header=False,
                new_columns=["sample", col_name]
            )

        raw    = read_counts(input.raw,    "raw_reads")
        dedup  = read_counts(input.dedup,  "dedup_reads")
        dehost = read_counts(input.dehost, "dehost_reads")

        # raw is per-lane (one row per lane), aggregate to sample_core
        # strip _S\d+_L\d+ suffix to get sample_core for joining with dedup/dehost
        raw_agg = (
            raw
            .with_columns(
                pl.col("sample")
                  .str.replace(r"_S\d+_L\d+$", "")
                  .alias("sample_core")
            )
            .group_by("sample_core")
            .agg(pl.col("raw_reads").sum())
            .rename({"sample_core": "sample"})
        )

        (
            raw_agg
            .join(dedup,  on="sample", how="left")
            .join(dehost, on="sample", how="left")
            .with_columns([
                (pl.col("dedup_reads")  / pl.col("raw_reads")   * 100).round(1).alias("pct_kept_dedup"),
                (pl.col("dehost_reads") / pl.col("dedup_reads") * 100).round(1).alias("pct_kept_dehost"),
            ])
            .write_csv(output.csv)
        )
