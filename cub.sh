#!/bin/zsh

source /home/lvelo/miniconda3/etc/profile.d/conda.sh 


# This file is part of the cub taxonomic sequence classification system.
# Runs the whole bash and snakemake pipeline.  

RUN_NAME=$1
NAS="/mnt/R60-11/Bacterio/metagenomique_NextSeq/"

RUN_FOLDER=${NAS}/${RUN_NAME}

PROJECT_PATH="/DATA/share/microbio/cub/"

# activate cub environment
conda activate cub
# RUN_FOLDER the pipeline
python $PROJECT_PATH/create_table.py $RUN_FOLDER
snakemake -nfp -d $RUN_FOLDER -s $PROJECT_PATH/Snakefile
#snakemake --use-conda -j 2 -d $RUN_FOLDER
