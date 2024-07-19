

test: 
	snakemake -n -c 40 --configfile config/config.yaml -p -d /mnt/R60-11/Bacterio/metagenomique_NextSeq/231214_NB501647_0242_AHH5YYAFX5 --rerun-incomplete
