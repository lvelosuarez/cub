"""
Author: L. Velo Suarez, lourdes.velosuarez@chu-brest.fr, lourdesvelo@gmail.com
Created: February 2022
Last Updated: March 2022
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
HUMAN = "/DATA/share/microbio/kraken2_databases/HUMAN"
MICROBIO = "/DATA/share/microbio/kraken2_databases/KRAKEN2"
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
        expand("results/03_not_human/report/{sample}.report", sample = sample_id),
        "results/QC/multiqc_report.html",
        "results/temp/Nreads_raw.txt",
        "results/temp/Nreads_qc.txt",
        "results/temp/Nreads_nothumain_1.txt",
        "results/temp/Nreads_nothumain_2.txt",
        'results/' + os.path.basename(os.getcwd()) + '.html'

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
        out = "results/03_not_human/report/{sample}.report"
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
        r1 = expand('results/03_not_human/report/{sample}.report', sample = sample_id)
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
        "results/temp/Nreads_raw.txt",
        "results/temp/Nreads_qc.txt",
        "results/temp/Nreads_nothumain_1.txt",
        "results/temp/Nreads_nothumain_2.txt",
        "results/QC/multiqc_report.html",
        #expand("results/03_not_human/report/{sample}.report",  sample = sample_id))

    output:
        'results/' + os.path.basename(os.getcwd()) + '.html'
    run:
        import sys 
        import os
        import pandas as pd
        import datapane as dp
        import numpy as np
        import plotly.express as px
        import argparse
        from datetime import date
        from collections import defaultdict

        date = date.today() 
        quality = input[4]
        chu_logo = """
        <svg id="Layer_1" data-name="chu_logo" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 768.81 451.59" width="100%" height="100%"  ><defs><style>.cls-1{fill:#1f3f96;}.cls-2{fill:#bd020f;}.cls-3{fill:#224198;}.cls-4{fill:#4fa1ca;}.cls-5{fill:#4fa2cb;}</style></defs><path class="cls-1" d="M134.76,150.6c.62-3.86,1.25-7.73,2.14-13.23l69.95,23.9V210c0,59.93.12,119.87,0,179.8-.13,49.29-38.59,87.78-87.19,87.65-48.79-.14-87.19-38.48-87.33-87.72-.19-71.92,0-143.84,0-215.76V161.13l68.06-23.46c1,4.87,1.89,9,2.75,13l-1.89,40.18Q100,207,98.84,223.15l-3,49.46q-1.35,18.27-2.71,36.55c-2.82,36.36-5.64,72.73-8.74,112.58l64.17-38.13c-2-25.77-3.84-50-5.71-74.27Q142,291,141,272.74q-1.09-24.9-2.21-49.82-.89-16-1.78-32Q135.91,170.75,134.76,150.6Z" transform="translate(-32.19 -25.89)"/><path class="cls-1" d="M103.23,122.39c-.65,3.88-1.3,7.76-2.29,13.72L33.61,113.43C29.45,68.69,70.5,26.54,117.67,25.9c49.48-.67,90.69,39.32,88.1,87l-68.09,23.18-3-13.72c.44-10.36-5.36-16-14.88-16.19C109.79,106,103.1,111.23,103.23,122.39Z" transform="translate(-32.19 -25.89)"/><path class="cls-2" d="M134.76,150.6q1.14,20.15,2.29,40.29H101.17l1.89-40.18Z" transform="translate(-32.19 -25.89)"/><path class="cls-2" d="M103.23,122.39c-.13-11.16,6.56-16.37,16.54-16.19,9.52.18,15.32,5.83,14.88,16.19Z" transform="translate(-32.19 -25.89)"/><path class="cls-2" d="M142.85,309.34c1.87,24.25,3.73,48.5,5.71,74.27L84.39,421.74c3.1-39.85,5.92-76.22,8.74-112.58Z" transform="translate(-32.19 -25.89)"/><path class="cls-2" d="M138.83,222.92l2.22,49.82-45.21-.13,3-49.46Z" transform="translate(-32.19 -25.89)"/><path class="cls-3" d="M480.59,185.86H371.37V153.8H479.89V30.73h36.67v95.15c0,13.73.44,27.12.44,27.12s13.4.33,27.25.33h28.21V30.61h35V259.9H573.12v-73H516.74V259.8H480.59Z" transform="translate(-32.19 -25.89)"/><path class="cls-3" d="M792.59,44.63q.09,75.38,0,150.75c-.1,39.49-26.67,67.6-64,68.05-37.7.46-65.87-28-66.13-67.82-.36-54.66-.09-109.31-.09-165h35.3c.23,4.88.69,10.1.7,15.31q.09,73.88,0,147.76c0,22.54,11,35.58,29.41,35.39,18.09-.19,29.73-13.9,29.77-35.54q.13-73.14,0-146.27V30.38h34" transform="translate(-32.19 -25.89)"/><path class="cls-3" d="M431.24,215.28c0,6.81,0,11.24,0,15.67a105.32,105.32,0,0,1-2.66,23,71.19,71.19,0,0,1-23.74,7.82C342,271,288.18,242,268,187.9c-.51-1.35-.61-2.85-1.68-8.11,20.37,10.59,42.06,10.46,58.13,27.13,24.72,25.66,61.87,30.69,93.87,14.9C421.86,220.09,425.33,218.29,431.24,215.28Z" transform="translate(-32.19 -25.89)"/><path class="cls-3" d="M426.42,73c-23.89-12.38-46.55-17.59-70.78-9.2C333.54,71.48,311.3,78.71,289.08,86c-4.07,1.34-8.3,2.19-15,3.93,13.39-26.44,32.23-43.26,56.71-53.48,25.91-10.83,52.66-13.36,80.38-6.12C417.77,32.07,427,37,427,37a80.91,80.91,0,0,1-.4,13.55C425.56,57.17,426.42,64.08,426.42,73Z" transform="translate(-32.19 -25.89)"/><g id="text87499"><path id="path43022" class="cls-4" d="M296.74,311.64H303a5.64,5.64,0,0,0,4.34-1.93,6.85,6.85,0,0,0,1.83-4.82,6.57,6.57,0,0,0-1.83-4.71,5.57,5.57,0,0,0-4.34-2h-6.27Zm-12,53.42V286.19H303A18.15,18.15,0,0,1,316.09,317a24.51,24.51,0,0,1,9.91,6.1,24.62,24.62,0,0,1-17.41,42Zm12-38.74v26.74h11.85a11.78,11.78,0,0,0,8.89-3.91,13.43,13.43,0,0,0,3.7-9.49,13.14,13.14,0,0,0-3.7-9.37,11.7,11.7,0,0,0-8.89-4Z" transform="translate(-32.19 -25.89)"/></g><g id="text91087"><path id="path43025" class="cls-4" d="M423,348l-.3,14.8h-12V289.08h12v.05c.71,0,1.43-.05,2.14-.05a26.09,26.09,0,0,1,19.08,7.82,26.52,26.52,0,0,1,0,37.71,30,30,0,0,1-4.13,3.4L450,363H438l-9-15h-6Zm-.3-47.84v35a12.29,12.29,0,0,0,2.14.17,14.43,14.43,0,0,0,11.2-5.18,18,18,0,0,0,4.66-12.51A17.84,17.84,0,0,0,436,305.23a14.46,14.46,0,0,0-11.2-5.18c-.71,0-1.43,0-2.14.11Z" transform="translate(-32.19 -25.89)"/></g><g id="text92813"><path id="path43028" class="cls-4" d="M535.14,288.19H567V298H546v30h21v12H546v12h21v10H535.14Z" transform="translate(-32.19 -25.89)"/></g><g id="text95653"><path id="path43031" class="cls-4" d="M666,363a25.59,25.59,0,0,1-8.32-1.24,27.44,27.44,0,0,1-7.68-4.43V343.17a32.48,32.48,0,0,0,3.91,3.86c.92.65,2,1.36,3.33,2.11a23.84,23.84,0,0,0,4.41,1.91,15.51,15.51,0,0,0,4.74.77,10.1,10.1,0,0,0,7.93-3.5,12.28,12.28,0,0,0,3.13-8.4,8.43,8.43,0,0,0-2.64-6.44,42.15,42.15,0,0,0-8.66-5.66c-2.42-1.31-4.45-2.49-6.12-3.56a27.12,27.12,0,0,1-4.79-4,18.6,18.6,0,0,1-3.72-5.62A19.36,19.36,0,0,1,650.1,307a20.28,20.28,0,0,1,2-8,21,21,0,0,1,4.9-6.55,22.43,22.43,0,0,1,15.31-5.92,21.55,21.55,0,0,1,6,1,25,25,0,0,1,6.12,3v11.9a17.31,17.31,0,0,0-2.84-1.75c-.49-.24-1.2-.62-2.15-1.14a14.4,14.4,0,0,0-3.23-1.39,12.84,12.84,0,0,0-3.87-.56,11,11,0,0,0-7.78,3.09,9.26,9.26,0,0,0-3.18,6.95,8,8,0,0,0,2.65,6,40.86,40.86,0,0,0,8.31,5q3.91,1.95,6.37,3.45a25.92,25.92,0,0,1,4.89,3.86,17.4,17.4,0,0,1,3.86,6,23.93,23.93,0,0,1,1.57,8.55,21.23,21.23,0,0,1-1.81,8.65,23.13,23.13,0,0,1-5,7.21,24.08,24.08,0,0,1-7.44,4.95A22.14,22.14,0,0,1,666,363Z" transform="translate(-32.19 -25.89)"/></g><g id="text96839"><path id="path43034" class="cls-5" d="M789,363h-9V297H768v-9h33v9H789Z" transform="translate(-32.19 -25.89)"/></g></svg>
        """
        cub_logo = """
        <svg id="cub_logo" data-name="cub_logo" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 145 73"  width="100%" height="100%" ><defs><style>.cls-1{fill:#0b5ba2;}.cls-2{fill:#333;}.cls-3{font-size:3.6px;fill:#999;font-family:SmoothCirculars, Smooth Circulars;letter-spacing:0.1em;}.cls-4{font-size:4.2px;letter-spacing:0.1em;}</style></defs><path class="cls-1" d="M4.67,40.05A31.38,31.38,0,0,1,25.49,8.43a30.55,30.55,0,0,1,7.74-1.66A31.41,31.41,0,0,1,59,16.65c.28.3.14.46-.07.69s-.39.14-.58-.06a30.35,30.35,0,0,0-17.8-9.41A30.55,30.55,0,0,0,6.62,29.68,29.32,29.32,0,0,0,6.3,45.12,30,30,0,0,0,23.37,65.81a28.84,28.84,0,0,0,11.18,2.71A30.51,30.51,0,0,0,60.85,55.7a.48.48,0,0,1,.87,0,.5.5,0,0,1-.13.43c-1,1.22-1.93,2.48-3,3.61a30.27,30.27,0,0,1-10.81,7.36A31.42,31.42,0,0,1,9,54.18,30.56,30.56,0,0,1,4.75,41.26C4.71,40.86,4.69,40.45,4.67,40.05Z"/><path class="cls-1" d="M10.05,39.71A26.68,26.68,0,0,0,12,47.87,25.83,25.83,0,0,0,29.88,63.08a24.94,24.94,0,0,0,9.27.43,25.31,25.31,0,0,0,13-5.71A25.65,25.65,0,0,0,56.6,53c.12-.16.18-.38.46-.19s.19.36,0,.56a26.62,26.62,0,0,1-10.81,8.78,26,26,0,0,1-9.4,2.18,25.74,25.74,0,0,1-14.27-3.44A26.38,26.38,0,0,1,9.67,34.19a25.49,25.49,0,0,1,2.71-8.32,26.18,26.18,0,0,1,8.47-9.54A25.2,25.2,0,0,1,29,12.59a26.33,26.33,0,0,1,25.9,7.46c.17.18.2.3,0,.47s-.31.13-.47,0a25.38,25.38,0,0,0-11.1-7,25.75,25.75,0,0,0-26,6.59,25.35,25.35,0,0,0-6.54,11.7A26.31,26.31,0,0,0,10.05,39.71Z"/><path class="cls-2" d="M59.76,37.7c.06-.13.13-.26.2-.39l.13,0a6,6,0,0,1,.79-.64c0-.17.08-.23.41-.34s.33-.2.27-.47c-.16-.72-.31-1.46-.45-2.19,0-.27-.17-.34-.43-.26a5.68,5.68,0,0,0-1.91.94A7.52,7.52,0,0,0,56.37,38a26,26,0,0,0-1.23,6.24c-.12,1.13-.24,2.27-.37,3.4a2,2,0,0,1-.13.39,1.19,1.19,0,0,1-.27-.42l-2.1-4.92A1.18,1.18,0,0,0,51,42l.12.32c.9,2.12,1.78,4.26,2.73,6.36a1.17,1.17,0,0,1-.23,1.47l-.12.13c-.25.28-.31.27-.47-.09-1.17-2.7-2.35-5.39-3.52-8.08-.24-.55-.24-.54-.8-.49-.29,0-.37.12-.24.4q1.83,4.2,3.67,8.4c.34.78.32.84-.49,1.11-.23.07-.36,0-.47-.19-.73-1.72-1.48-3.43-2.23-5.14-.55-1.26-1.1-2.52-1.66-3.77,0-.08-.18-.18-.24-.16a.71.71,0,0,0-.42.91q1.32,3,2.65,6.06c.32.74.64,1.48.94,2.22a.52.52,0,0,1,0,.33,3.83,3.83,0,0,1-1.07,0c-.16,0-.27-.35-.36-.56-.93-2.17-1.85-4.35-2.78-6.53,0-.1-.1-.19-.15-.29a2,2,0,0,0-.23,1.78c.72,1.63,1.41,3.28,2.11,4.92.06.14.09.29.16.51a2.24,2.24,0,0,1-.44-.1c-.67-.24-1.34-.48-2-.74-.21-.08-.33-.08-.35.19-.09.8-.19,1.6-.29,2.39a.4.4,0,0,0,.3.5,13,13,0,0,0,4.85.94,7.5,7.5,0,0,0,7.46-5.67A29.45,29.45,0,0,0,58,45.52c.21-1.5.4-3,.61-4.51A4.65,4.65,0,0,1,59.76,37.7Z"/><path class="cls-2" d="M54.42,42.27l-.5-.23a11.23,11.23,0,0,0-4.53-1.31,4.15,4.15,0,0,0-4.34,2.86,11,11,0,0,0-.7,3.77c-.09,1.42-.11,2.84-.22,4.26a9,9,0,0,1-5.34,7.78c-.19.1-.31.11-.4-.13-.27-.67-.55-1.34-.84-2-.09-.21-.11-.35.15-.47a6.5,6.5,0,0,0,3.72-5.5c.15-1.33.07-2.68.11-4A10.58,10.58,0,0,1,44,40a6.64,6.64,0,0,1,5.07-2.12,13.08,13.08,0,0,1,5.55,1.42.42.42,0,0,1,.21.3C54.74,40.49,54.58,41.34,54.42,42.27Z"/><path class="cls-2" d="M40.43,51c-.31-.29-.63-.57-.93-.87a13.08,13.08,0,0,1-3.72-7.51,18.07,18.07,0,0,1,0-2.55c0-.41.12-.46.51-.32l2.1.75a.38.38,0,0,1,.29.44,9.19,9.19,0,0,0,1.61,5.4,1.81,1.81,0,0,1,.36,1.07c0,1.07,0,2.14,0,3.21a1.88,1.88,0,0,1-.05.33Z"/><path class="cls-2" d="M16.91,36.05A3.2,3.2,0,0,1,13.62,33a3.16,3.16,0,0,1,3.31-3,3.2,3.2,0,0,1,3.35,3.1A3.26,3.26,0,0,1,16.91,36.05Z"/><path class="cls-2" d="M33.18,14.4a3,3,0,0,1,3.28,3,3.13,3.13,0,0,1-3.37,3.08,3.25,3.25,0,0,1-3.19-3.08A3,3,0,0,1,33.18,14.4Z"/><path class="cls-2" d="M36.36,28.75a3,3,0,0,1-3.09-3,3.1,3.1,0,0,1,3.16-3,3.22,3.22,0,0,1,3.34,2.95C39.7,27.63,38.45,28.77,36.36,28.75Z"/><path class="cls-2" d="M25.07,23.47a3,3,0,0,1-3.16,2.83,2.85,2.85,0,1,1,.17-5.69A2.88,2.88,0,0,1,25.07,23.47Z"/><path class="cls-2" d="M23,51.3a2.57,2.57,0,0,1-2.82,2.53,2.68,2.68,0,0,1-2.75-2.65,2.78,2.78,0,0,1,2.72-2.52A2.69,2.69,0,0,1,23,51.3Z"/><path class="cls-2" d="M42.1,33a2.32,2.32,0,0,1-2.24,2.27A2.49,2.49,0,0,1,37.31,33a2.63,2.63,0,0,1,2.32-2.21A2.5,2.5,0,0,1,42.1,33Z"/><path class="cls-2" d="M20.21,42.15a2,2,0,0,1-2.05,2,2.13,2.13,0,0,1-2.2-2,2.1,2.1,0,0,1,2.2-2A2,2,0,0,1,20.21,42.15Z"/><path class="cls-2" d="M33.39,48.89a1.82,1.82,0,1,1,.08-3.61,1.76,1.76,0,0,1,1.94,1.79A1.83,1.83,0,0,1,33.39,48.89Z"/><path class="cls-2" d="M27.18,29.44a2.07,2.07,0,0,1-2-1.91A2.11,2.11,0,0,1,27,25.72a2.26,2.26,0,0,1,2.1,2A2,2,0,0,1,27.18,29.44Z"/><path class="cls-2" d="M30,36.21A1.88,1.88,0,0,1,31.89,38,2,2,0,0,1,30,39.76,2.21,2.21,0,0,1,28,38,2,2,0,0,1,30,36.21Z"/><path class="cls-2" d="M23.05,38a1.62,1.62,0,1,1,0,3.23,1.62,1.62,0,1,1,0-3.23Z"/><path class="cls-2" d="M28.08,57.81a1.44,1.44,0,0,1-1.43-1.39A1.41,1.41,0,0,1,28.17,55a1.4,1.4,0,1,1-.09,2.78Z"/><path class="cls-2" d="M24.72,46.72c.44.26.88.53,1.31.81a.57.57,0,0,1,.23.31c0,.77,0,1.54,0,2.3,0,.14-.07.19-.2.19h-.4c-.15,0-.2.06-.2.2s0,.26,0,.39A4.68,4.68,0,0,1,27,50.67a4.88,4.88,0,0,1,1.77.32c0-.17,0-.34,0-.51a.2.2,0,0,0-.13-.14,3.23,3.23,0,0,0-.44,0c-.16,0-.23,0-.23-.22,0-.74,0-1.48,0-2.21a.31.31,0,0,1,.18-.32l1.43-.85a.32.32,0,0,0,.17-.32c0-.9,0-1.8,0-2.69a.32.32,0,0,0-.17-.31L27.27,42A.3.3,0,0,0,27,42l-2.33,1.39a.36.36,0,0,0-.13.24c0,.91,0,1.82,0,2.73A.34.34,0,0,0,24.72,46.72Z"/><path class="cls-2" d="M66.54,50.41a16.42,16.42,0,0,0,23.2,0,15.81,15.81,0,0,0,4.8-11.6V21.23h2.35V38.81a18,18,0,0,1-5.5,13.25,18.72,18.72,0,0,1-26.51,0,18.06,18.06,0,0,1-5.49-13.25V21.23h2.34V38.81A15.78,15.78,0,0,0,66.54,50.41Z"/><path class="cls-2" d="M101.06,29.64a19.46,19.46,0,0,1,3.17-4.1,18.75,18.75,0,0,1,26.5,0,18.07,18.07,0,0,1,5.49,13.26,18.06,18.06,0,0,1-5.49,13.25,18.73,18.73,0,0,1-26.52,0,18.06,18.06,0,0,1-5.49-13.25V14.2h2.34Zm16.41,25.57a16.41,16.41,0,0,0,16.41-16.4,16.4,16.4,0,0,0-28-11.6,15.82,15.82,0,0,0-4.81,11.6,16.41,16.41,0,0,0,16.41,16.4Z"/><text class="cls-3" transform="translate(58.28 64.53)">ngs microbiology analysis<tspan class="cls-4" x="79.88" y="0"> </tspan></text></svg>
        """
        ###
        ## Charge txt files to calculate datatable for raw/humans reads
        #### Read sample sheet from run info
        with open('SampleSheet.csv', 'r') as csvfile:
            lines = csvfile.readlines()
        sections = [lines.index(line) for line in lines if "[" in line]
        SampleSheet = lines[sections[3]+2:]
        header = lines[sections[3]+1].split(sep=",")
        SampleSheet = pd.DataFrame([line.split(",") for line in SampleSheet], columns = header)[['Sample_ID','Description']].set_index('Sample_ID')
        ##### Read info from cub
        df = pd.read_table(input[0], sep=",",names=['raw'],index_col=0)
        df = df.join(pd.read_table(input[1], sep=",",names=['quality'],index_col=0))
        df = df.join(pd.read_table(input[2], sep=",",names=['human1'],index_col=0))
        df = df.join(pd.read_table(input[3], sep=",",names=['human2'],index_col=0))
        df['nonHuman'] = df["human1"] - df["human2"]
        df['percent'] = df['nonHuman'] / df['raw'] * 100
        df = df.drop(labels="Undetermined", axis=0)
        df =df.join(SampleSheet)
        df.index = df.index.astype(int)
        df = df[['Description','raw','quality','nonHuman','percent']].sort_index()
        df['Sample'] = df.index
        df = df[['Sample','Description','raw','quality','nonHuman','percent']]
        df1 = pd.melt(df[['Sample','raw','quality','nonHuman']], id_vars=['Sample'], value_vars=['raw','quality','nonHuman'], var_name='steps', value_name='reads')
        fig = px.line(df1, x="steps", y="reads", color='Sample', markers=True, log_y=True) 
        #### 
        ## Charge all kraken reports and use SPECIES info
        def input_dir(dir):
	        files=os.listdir(dir)
	        path=[]
	        for f in files:
		        fipath=os.path.join(dir,f)
		        path.append(fipath)	
	        return (path)
        filelist=input_dir('results/03_not_human/report/')
        alldf = []
        for file in filelist:
            if 'Undetermined' in file:
                continue
            out = pd.read_csv(file, sep="\t", names=['percent','reads_root','reads_assigned','rank', 'taxon_id','name'])
            out = out[out["rank"] == "S"]
            out["name"] = out["name"].str.strip()
            out["Sample"] = file.split("/")[-1].split(".")[0]
            alldf.append(out)
        superdf = pd.concat(alldf)
        table = superdf[["name","reads_root", "Sample"]]
        table = table.pivot_table(index="name", columns="Sample", values="reads_root")
        table = table.drop(labels="Homo sapiens", axis=0)
        table = table.fillna(0)
        table = table[table.columns.astype(int).sort_values().astype(str)]
        table.loc[:,'Total'] = table.sum(axis=1)
        table = table[table["Total"] > 100]
        #######
        ## Create the report
        ######
        report = dp.Report(
                    dp.Group(
                        dp.HTML(cub_logo),
                        dp.HTML(' <span style="display:inline-block;border:0.01px solid #FFFFFF;width:150px;height:50px;"></span>'),
                        dp.HTML(' <span style="display:inline-block;border:0.01px solid #FFFFFF;width:150px;height:50px;"></span>'),
                        dp.HTML(' <span style="display:inline-block;border:0.01px solid #FFFFFF;width:150px;height:50px;"></span>'),
                        dp.HTML(' <span style="display:inline-block;border:0.01px solid #FFFFFF;width:150px;height:50px;"></span>'),
                        dp.HTML(chu_logo),
                        columns = 6
                    ),
                    '# Run NextSeq ' + os.path.basename(os.getcwd()),
                    dp.Group(
                        dp.Group(
                          dp.BigNumber(heading="Analized on",value = f"{date:%d-%m-%Y}"),
                          dp.Plot(fig, caption = "Read lost per sample"),
                          dp.Attachment(file = quality)
                            ),
                    dp.Table(df.style.format(precision=2, thousands=" ").hide().set_properties(subset=["Description","raw", "quality","nonHuman","percent"], **{'text-align': 'right'}), caption = "Number of reads in raw/qc fastq files and human amount proportion per sample"),
                    columns = 2
                    ),
                    dp.DataTable(table)
                    )                        
        report.save(path = 'results/' + os.path.basename(os.getcwd()) + '.html', formatting=dp.ReportFormatting(text_alignment=dp.TextAlignment.CENTER, width=dp.ReportWidth.FULL))