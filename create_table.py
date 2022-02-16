#! /usr/bin/env python
import os
import logging
import pandas as pd
from collections import defaultdict

def get_sample_files(path, outfile = 'samples.tsv'):
    samples = defaultdict(dict)
    pattern = "results"
    for dir_names, _, files in os.walk(os.path.abspath('DATA')):
        if pattern in dir_names:
            continue
        for fname in files:
            if ".fastq" in fname or ".fq" in fname:
                sample_id = fname.split(".fastq")[0].split(".fq")[0]
                sample_id = sample_id.replace("_R1", "").replace("_R2", "").replace("_001", "").replace(" ", "-")
                fq_path = os.path.join(dir_names, fname)
                if "_R1" in fname:
                    samples[sample_id]['R1'] = fq_path
                else:
                    samples[sample_id]['R2'] = fq_path
    samples= pd.DataFrame(samples).T
    samples.to_csv(os.path.abspath(path) + "/" + outfile,sep='\t')
    
    return samples

if __name__ == '__main__':
    import sys
    get_sample_files(sys.argv[1])