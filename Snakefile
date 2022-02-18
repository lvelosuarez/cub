"""
Author: L. Velo Suarez, lourdes.velosuarez@chu-brest.fr, lourdesvelo@gmail.com
Created: February 2022
Last Updated: February 2022
Affiliation: CBAM (Centre Brestois Analyse Microbiota), CRHU Brest
Aim: Snakemake workflow to process paired-end shot gun from Illumina reads (NextSeq)
Run: snakemake -j3 
Run for dag : snakemake --dag | dot -Tsvg > dag.svg
Latest modification: 
  - todo
"""
import os
from glob import glob
import pandas as pd 

##### DEFINE PATHS
HUMAN = "/DATA/share/microbio/kraken2_databases/HUMAN"
MICROBIO = "/DATA/share/microbio/kraken2_databases/KRAKEN2"
### Read and write in NAS
nas = "/mnt/R60-11/Bacterio/metagenomique_NextSeq/"
#### DEFINE samples names
df = pd.read_csv('samples.tsv', sep='\t', index_col=0)
df['sample'] = df.index.to_series().str.split('_').str[0]
sample_id = list(df['sample'].unique())
i =["R1","R2"]
########## USED FUNCTIONS ####
def get_raw_fastq(id):
    ''' return a dict with the path to the raw fastq files'''
    r1 = list(df[df["sample"] == id]['R1'])
    r2 = list(df[df["sample"] == id]['R2'])
    return {'r1': r1, 'r2': r2}
###############################
rule all:
    input:
        expand("results/03_not_human/{sample}.report", sample = sample_id),
        "results/QC/multiqc_report.html",
        "results/temp/Nreads_raw.txt",
        "results/temp/Nreads_qc.txt",
        "results/temp/Nreads_nothumain_1.txt",
        "results/temp/Nreads_nothumain_2.txt"

rule merge_lanes:
    input: 
        unpack(lambda wildcards: get_raw_fastq(wildcards.sample))
    output: 
        r1 = 'results/00_raw/{sample}_R1.fastq.gz',
        r2 = 'results/00_raw/{sample}_R2.fastq.gz'
    shell: """
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
        ref=adapters,artifacts,phix,lambda,pjet,mtst,kapa \
        out={output.r1} out2={output.r2} \
        qtrim=rl trimq=20 maq=20  minlen=100       
        """
rule filter_human:
    input:
        r1 = rules.QC.output.r1,
        r2 = rules.QC.output.r2
    output:
        out = "results/02_human/{sample}.out"
    params:
        ref = HUMAN
    shell:
        "kraken2 --db {params.ref} --threads 20 --paired --gzip-compressed {input.r1} {input.r2} --output {output.out} "

rule filter_human_list:
    input:
        kraken = rules.filter_human.output.out
    output:
        lista = "results/02_human/{sample}.list"
    params: 
        awk = """awk '/U/{print $2}'"""
    shell:
        """
       cat {input.kraken} | {params.awk} > {output.lista}
        """  

rule filter_human2:
    input:
        r1 = rules.QC.output.r1,
        r2 = rules.QC.output.r2,
        names = rules.filter_human_list.output.lista
    output:
        r1 = "results/03_not_human/{sample}_R1.fastq.gz",
        r2 = "results/03_not_human/{sample}_R2.fastq.gz",
    shell:
        '''
        filterbyname.sh in={input.r1} in2={input.r2} \
        out={output.r1} out2={output.r2} names={input.names} include=t
        ''' 

rule kraken2_microbio:
    input:
        r1 = rules.filter_human2.output.r1,
        r2 = rules.filter_human2.output.r2,
    output:
        out = "results/03_not_human/{sample}.report"
    params:
        ref = MICROBIO
    shell:
        "kraken2 --db {params.ref} --threads 10 --paired --gzip-compressed {input.r1} {input.r2} --output '-' --report {output.out} "

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
        expand(["results/QC/{sample}_{i}_fastqc.html", 
                "results/QC/{sample}_{i}_fastqc.zip",],sample=sample_id,i=i)
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
        15 #config['threads']
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
        15 #config['threads']
    shell:
        """ 
        seqkit stats -j {threads} --basename {input.r1} | 
        {params.awk} | 
        sed 's/_R1.fastq.gz//' | 
        sed 's/,//g'| 
        sed 's/ /,/g' > {output.nreads}
        """
rule count_nothuman_1_reads:
    input:
        r1 = expand('results/03_not_human/{sample}_R1.fastq.gz', sample = sample_id)
    output:
        nreads= "results/temp/Nreads_nothumain_1.txt"
    params:
        awk = """awk '!/file/{print $1,$4}'"""
    threads:
        15 #config['threads']
    shell:
        """ 
        seqkit stats -j {threads} --basename {input.r1} | 
        {params.awk} | 
        sed 's/_R1.fastq.gz//' | 
        sed 's/,//g'| 
        sed 's/ /,/g' > {output.nreads}
        """
rule count_nothuman_2_reads:
    input:
        r1 = expand('results/03_not_human/{sample}.report', sample = sample_id)
    output:
        nreads = "results/temp/Nreads_nothumain_2.txt"
    params:
        awk = """awk '/Homo sapiens/{print FILENAME , $2}'"""
    shell:
        """ 
        {params.awk} {input.r1} | 
        sed 's/\.\///g' | 
        sed 's/.report//g' | 
        sed 's/results\/03_not_human\///g' |
        sed 's/ /,/g' >  {output.nreads}
        """
#### Create report :
rule report:
    input:
        