# Snakemake rules imported in the main Snakefile to count read numbers in raw/quality/cutadapter R1 files
rule count_raw_reads:
    input:
        r1 =  SAMPLES.R1.values
    output:
        nreads= temp("output/Nreads_raw.txt")
    threads:
        config['threads']
    params:
        awk = """awk '!/file/{print $1,$4}'"""
    shell:
        """ 
        seqkit stats -j {threads} --basename {input.r1} | {params.awk} | sed 's/_R1_001.fastq.gz//' | sed 's/,//g'| sed 's/ /,/g' > {output.nreads}
        """ 
rule count_qc_reads:
    input:
        r1 = expand("output/qc/{names}_R1.fastq.gz", names=SAMPLES.names.values)
    output:
        nreads= temp("output/Nreads_quality.txt")
    threads:
        config['threads']
    params:
        awk = """awk '!/file/{print $1,$4}'"""
    shell:
        """ 
        seqkit stats -j {threads} --basename {input.r1} | {params.awk} | sed 's/_R1.fastq.gz//' | sed 's/,//g'| sed 's/ /,/g' > {output.nreads}
        """
rule count_dehost_reads:
    input:
        r1 = lambda wildcards: [
            f
            for sample in SAMPLES["names"].to_list()
            for f in get_dehost_fastq_for_sample(sample)["R1"]
        ]
    output:
        nreads = temp("output/Nreads_host.txt")
    threads:
        config["threads"]
    params:
        # Strip _S{n}_L{nnn}_dedup_R1.clean_1.fastq.gz → sample name, then sum reads per sample
        awk_agg = r"""awk 'NR>1 {
            name = $1;
            gsub(/_S[0-9]+_L[0-9]+_dedup_R1\.clean_1\.fastq\.gz$/, "", name);
            sum[name] += $4
        } END {
            for (s in sum) print s "," sum[s]
        }'"""
    shell:
        """
        seqkit stats -j {threads} --basename {input.r1} \
            | {params.awk_agg} \
            > {output.nreads}
        """
# rule count_dehost_reads:
#     input:
#         r1 = expand("output/dehost/{names}_dedup_R1.clean.fastq.gz",names=SAMPLES.names.values)
#     output:
#         nreads= temp("output/Nreads_host.txt")
#     threads:
#         config['threads']
#     params:
#         awk = """awk '!/file/{print $1,$4}'"""
#     shell:
#         """ 
#         seqkit stats -j {threads} --basename {input.r1} | {params.awk} | sed 's/_dedup_R1.clean.fastq.gz//' | sed 's/,//g'| sed 's/ /,/g' > {output.nreads}
#         """
rule combine_read_counts:
    input:
        'output/Nreads_raw.txt',
        'output/Nreads_quality.txt',
        'output/Nreads_host.txt'
    output:
        report('output/Nreads.csv', category="QC reads"),
    run:
        import pandas as pd
        from natsort import natsorted
        
        D = pd.read_table(input[0],sep=",",names=['nextseq'],index_col=0)                                                                                               
        D = D.join(pd.read_table(input[1],sep=",",names=['quality'],index_col=0))
        D = D.join(pd.read_table(input[2],sep=",",names=['host'],index_col=0))
        D = D[['nextseq','quality','host']]
        sorted_index = natsorted(D.index)
        D = D.reindex(sorted_index)
        D.to_csv(output[0],sep=',')
