# rules/count_lib3.smk
#
# Read counting for lib3 (poly-T tagging pipeline).
# Tracks reads at each stage and per stream (tagged / untagged).
#
# Stages:
#   raw       — original input reads
#   qc        — after bbduk adapter trimming
#   tagged    — after split_polyt: R2 IS poly-T (genuine DNA)
#   untagged  — after split_polyt: R2 is NOT poly-T (library-prep contaminant)
#   dedup_tagged / dedup_untagged   — after clumpify per stream
#   dehost_tagged / dehost_untagged — after hostile per stream
#
# Output: output/Nreads.csv
#   sample, raw, qc, tagged, untagged, tag_ratio,
#   dedup_tagged, dedup_untagged, dehost_tagged, dehost_untagged,
#   pct_kept_dedup_tagged, pct_kept_dedup_untagged,
#   pct_kept_dehost_tagged, pct_kept_dehost_untagged


rule count_raw_reads:
    input:
        r1 = expand("{r1}", r1=SAMPLES["R1"].dropna().tolist())
    output:
        temp("output/Nreads_raw.txt")
    threads: config["threads"]
    shell:
        """
        seqkit stats -j {threads} --basename {input.r1} \
            | awk 'NR>1 {{gsub(/_R1\.fastq\.gz$/, "", $1); print $1","$4}}' \
            > {output}
        """


rule count_qc_reads:
    input:
        r1 = expand("output/qc/{sample}_R1.fastq.gz", sample=SAMPLE_CORES)
    output:
        temp("output/Nreads_qc.txt")
    threads: config["threads"]
    shell:
        """
        seqkit stats -j {threads} --basename {input.r1} \
            | awk 'NR>1 {{gsub(/_R1\.fastq\.gz$/, "", $1); print $1","$4}}' \
            > {output}
        """


rule count_split_reads:
    input:
        tagged   = expand("output/split/{sample}_tagged_R1.fastq.gz",   sample=SAMPLE_CORES),
        untagged = expand("output/split/{sample}_untagged_R1.fastq.gz", sample=SAMPLE_CORES),
    output:
        tagged   = temp("output/Nreads_tagged.txt"),
        untagged = temp("output/Nreads_untagged.txt"),
    threads: config["threads"]
    shell:
        """
        seqkit stats -j {threads} --basename {input.tagged} \
            | awk 'NR>1 {{gsub(/_tagged_R1\.fastq\.gz$/, "", $1); print $1","$4}}' \
            > {output.tagged}
        seqkit stats -j {threads} --basename {input.untagged} \
            | awk 'NR>1 {{gsub(/_untagged_R1\.fastq\.gz$/, "", $1); print $1","$4}}' \
            > {output.untagged}
        """


rule count_dedup_reads:
    input:
        tagged   = expand("output/dedup/{sample}_tagged_dedup_R1.fastq.gz",   sample=SAMPLE_CORES),
        untagged = expand("output/dedup/{sample}_untagged_dedup_R1.fastq.gz", sample=SAMPLE_CORES),
    output:
        tagged   = temp("output/Nreads_dedup_tagged.txt"),
        untagged = temp("output/Nreads_dedup_untagged.txt"),
    threads: config["threads"]
    shell:
        """
        seqkit stats -j {threads} --basename {input.tagged} \
            | awk 'NR>1 {{gsub(/_tagged_dedup_R1\.fastq\.gz$/, "", $1); print $1","$4}}' \
            > {output.tagged}
        seqkit stats -j {threads} --basename {input.untagged} \
            | awk 'NR>1 {{gsub(/_untagged_dedup_R1\.fastq\.gz$/, "", $1); print $1","$4}}' \
            > {output.untagged}
        """


rule count_dehost_reads:
    input:
        tagged   = expand("output/dehost/{sample}_tagged_R1.clean_1.fastq.gz",   sample=SAMPLE_CORES),
        untagged = expand("output/dehost/{sample}_untagged_R1.clean_1.fastq.gz", sample=SAMPLE_CORES),
    output:
        tagged   = temp("output/Nreads_dehost_tagged.txt"),
        untagged = temp("output/Nreads_dehost_untagged.txt"),
    threads: config["threads"]
    shell:
        """
        seqkit stats -j {threads} --basename {input.tagged} \
            | awk 'NR>1 {{gsub(/_tagged_R1\.clean_1\.fastq\.gz$/, "", $1); print $1","$4}}' \
            > {output.tagged}
        seqkit stats -j {threads} --basename {input.untagged} \
            | awk 'NR>1 {{gsub(/_untagged_R1\.clean_1\.fastq\.gz$/, "", $1); print $1","$4}}' \
            > {output.untagged}
        """


rule merge_counts:
    input:
        raw            = "output/Nreads_raw.txt",
        qc             = "output/Nreads_qc.txt",
        tagged         = "output/Nreads_tagged.txt",
        untagged       = "output/Nreads_untagged.txt",
        dedup_tagged   = "output/Nreads_dedup_tagged.txt",
        dedup_untagged = "output/Nreads_dedup_untagged.txt",
        dehost_tagged  = "output/Nreads_dehost_tagged.txt",
        dehost_untagged= "output/Nreads_dehost_untagged.txt",
    output:
        csv = "output/Nreads.csv"
    run:
        import polars as pl

        def read_counts(path, col):
            return pl.read_csv(path, has_header=False, new_columns=["sample", col])

        raw            = read_counts(input.raw,             "raw")
        qc             = read_counts(input.qc,              "qc")
        tagged         = read_counts(input.tagged,          "tagged")
        untagged       = read_counts(input.untagged,        "untagged")
        dedup_tagged   = read_counts(input.dedup_tagged,    "dedup_tagged")
        dedup_untagged = read_counts(input.dedup_untagged,  "dedup_untagged")
        dehost_tagged  = read_counts(input.dehost_tagged,   "dehost_tagged")
        dehost_untagged= read_counts(input.dehost_untagged, "dehost_untagged")

        (
            raw
            .join(qc,              on="sample", how="left")
            .join(tagged,          on="sample", how="left")
            .join(untagged,        on="sample", how="left")
            .join(dedup_tagged,    on="sample", how="left")
            .join(dedup_untagged,  on="sample", how="left")
            .join(dehost_tagged,   on="sample", how="left")
            .join(dehost_untagged, on="sample", how="left")
            .with_columns([
                (
                    pl.col("tagged") /
                    (pl.col("tagged") + pl.col("untagged")) * 100
                ).round(1).alias("tag_ratio_pct"),
                (pl.col("dedup_tagged")    / pl.col("tagged")    * 100).round(1).alias("pct_kept_dedup_tagged"),
                (pl.col("dedup_untagged")  / pl.col("untagged")  * 100).round(1).alias("pct_kept_dedup_untagged"),
                (pl.col("dehost_tagged")   / pl.col("dedup_tagged")   * 100).round(1).alias("pct_kept_dehost_tagged"),
                (pl.col("dehost_untagged") / pl.col("dedup_untagged") * 100).round(1).alias("pct_kept_dehost_untagged"),
            ])
            .write_csv(output.csv)
        )
