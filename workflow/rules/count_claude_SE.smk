# rules/count.smk

AWK_SUM_BY_LANE = r"""awk 'NR>1 {
    name = $1;
    gsub(/_S[0-9]+_L[0-9]+_dedup_R1\.fastq\.gz$/, "", name);
    gsub(/,/, "", $4);
    sum[name] += $4
} END {
    for (s in sum) print s "," sum[s]
}'"""

rule count_dedup_reads:
    input:
        r1 = expand(
            "output/dedup/{names}_dedup_R1.fastq.gz",
            names = SAMPLES["names"].tolist()
        )
    output:
        nreads = "output/Nreads.csv"
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


