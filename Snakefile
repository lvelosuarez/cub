"""
Author: L. Velo Suarez, lourdes.velosuarez@chu-brest.fr, lourdesvelo@gmail.com
Created: February 2022
Last Updated: February 2023
Affiliation: CBAM (Centre Brestois Analyse Microbiota), CRHU Brest
Aim: Snakemake workflow to process paired-end shot gun from Illumina reads (NextSeq)
Run: snakemake -j3 
Run for dag : snakemake --dag | dot -Tsvg > dag.svg
Latest modification: 
- create final rule for reporting
"""
import os
from glob import glob
import pandas as pd 

##### DEFINE PATHS
PlusPF = "/DATA/share/microbio/index_zone/kraken2/PlusPF/"
#### DEFINE samples names
df = pd.read_csv('samples.tsv', sep='\t', index_col=0)
df['sample'] = df.index.to_series().str.split('_').str[0]
sample_id = list(df['sample'].unique())
i =["R1","R2"]
########## USED FUNCTIONS ####
def get_raw_fastq(samples):
    ''' return a dict with the path to the raw fastq files'''
    r1 = list(df[df["sample"] == samples]['R1'])
    r2 = list(df[df["sample"] == samples]['R2'])
    return {'r1': r1, 'r2': r2}
###############################
rule all:
    input:
        expand("results/02_kraken2PlusPF/{sample}.krk", sample = sample_id),
        "results/QC/multiqc_report.html",
        "results/temp/Nreads_raw.txt",
        "results/temp/Nreads_qc.txt"

rule merge_lanes:
    input: 
        unpack(lambda wildcards: get_raw_fastq(wildcards.sample))
    output: 
        r1 = 'results/00_raw/{sample}_R1.fastq.gz',
        r2 = 'results/00_raw/{sample}_R2.fastq.gz'
    shell: 
        """
        zcat {input.r1} | pigz -p 10 > {output.r1}
        zcat {input.r2} | pigz -p 10 > {output.r2}
        """

rule QC:
    input:
        r1 = rules.merge_lanes.output.r1,
        r2 = rules.merge_lanes.output.r2
    output:
        r1 = "results/01_QC/{sample}_R1.fastq.gz",
        r2 = "results/01_QC/{sample}_R2.fastq.gz"
    shell: 
        """
        bbduk.sh in={input.r1} in2={input.r2} \
        ref=adapters,artifacts,phix,lambda,pjet,mtst,kapa  \
        out={output.r1} out2={output.r2} \
        qtrim=rl trimq=20 maq=20  minlen=100       
        """
rule kraken:
    input:
        r1 = rules.QC.output.r1,
        r2 = rules.QC.output.r2
    output:
        out = "results/02_kraken2PlusPF/{sample}.krk",
        report = "results/02_kraken2PlusPF/{sample}.report"
    params:
        ref = PlusPF
    shell:
        """
        kraken2 --db {params.ref} --threads 10 \
        --output {output.out} --report {output.report}  --report-minimizer-data \
        --paired --gzip-compressed {input.r1} {input.r2}
        """
####### Rules QC ######
rule fastqc:
    input: 
        expand("results/00_raw/{sample}_{i}.fastq.gz", sample=sample_id, i=i)
    output: 
        expand("results/QC/{sample}_{i}_fastqc.html", sample=sample_id, i=i)
    shell: 
        "fastqc -o results/QC -t 10 {input}"

rule multiqc:
    input: 
        expand("results/QC/{sample}_{i}_fastqc.html",sample=sample_id,i=i)
    output: 
        "results/QC/multiqc_report.html"
    shell:
        "multiqc results/QC -o results/QC"
## Count reads 
rule count_raw_reads:
    input:
        r1 = expand('results/00_raw/{sample}_R1.fastq.gz', sample = sample_id)
    output:
        nreads= "results/temp/Nreads_raw.txt"
    params:
        awk = """awk '!/file/{print $1,$4}'"""
    threads:
        10 #config['threads']
    shell:
        """ 
        seqkit stats -j {threads} --basename {input.r1} | 
        {params.awk} | 
        sed 's/_R1.fastq.gz//' | 
        sed 's/,//g'| 
        sed 's/ /,/g' > {output.nreads}
        """
rule count_qc_reads:
    input:
        r1 = expand('results/01_QC/{sample}_R1.fastq.gz', sample = sample_id)
    output:
        nreads= "results/temp/Nreads_qc.txt"
    params:
        awk = """awk '!/file/{print $1,$4}'"""
    threads:
        10 #config['threads']
    shell:
        """ 
        seqkit stats -j {threads} --basename {input.r1} | 
        {params.awk} | 
        sed 's/_R1.fastq.gz//' | 
        sed 's/,//g'| 
        sed 's/ /,/g' > {output.nreads}
        """