#!/bin/bash

# This file is part of the cub taxonomic sequence classification system.
# Runs the whole bash and snakemake pipeline.  

run="$1"
nas="/mnt/R60-11/Bacterio/metagenomique_NextSeq/"
# run the pipeline
 ./create_table.py $run
snakemake -nfp -d $run