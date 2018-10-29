set -e
snakemake --snakefile Snakefile --configfile=config/test-config.yaml --resources wget_limit=2 -j 20
