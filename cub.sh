#!/bin/zsh

# proton path 
# source /home/lvelo/miniconda3/etc/profile.d/conda.sh 

# fermion path
source /data/lourdes/miniconda3/etc/profile.d/conda.sh

# This file is part of the cub taxonomic sequence classification system.
# Runs the whole bash and snakemake pipeline.  

RUN_NAME=$1
# This is the real NAS but I will work sometimes in run located in local NAS
#NAS="/mnt/R60-11/Bacterio/metagenomique_NextSeq"
NAS="/data/lourdes/NAS"
RUN_FOLDER=${NAS}/${RUN_NAME}

# proton path 
#PROJECT_PATH="/DATA/share/microbio/cub"
# fermion path 
PROJECT_PATH="/data/lourdes/metagenomique_clinique"
# activate cub environment
conda activate cub # fermion environment
# RUN_FOLDER the pipeline
echo $RUN_FOLDER
python $PROJECT_PATH/create_table.py $RUN_FOLDER
#snakemake -nfp -d $RUN_FOLDER -s $PROJECT_PATH/Snakefile #--rerun-incomplete
snakemake -pd $RUN_FOLDER -s $PROJECT_PATH/Snakefile -j 5
