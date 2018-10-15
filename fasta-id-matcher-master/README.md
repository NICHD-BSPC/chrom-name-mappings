# fasta-id-matcher

FASTA files from different providers do not always have the same chromosome IDs
(e.g., UCSC vs Ensembl vs NCBI). One option is
https://github.com/dpryan79/ChromosomeMappings, but it's not always clear how
the mappings were created and if a genome is missing it needs to be created and
added to that repository.

The script included here attempts to create a TSV mapping one chromosome naming
scheme to another by randomly sampling the FASTA records to identify matches.

Requires `pyfaidx` package which can be installed via pip or conda. To set up
a fresh environment just for this, use:

```bash
conda create -n fasta-id-matcher pyfaidx -c bioconda -c conda-forge --file requirements.txt
```

and activate the environment:

```bash
source activate fasta-id-matcher
```


Once in an environment with `pyfaidx` installed run to inspect help:

```bash
python fasta-id-matcher.py -h
```

to view the help.

Example usage:

```bash
python fasta-id-matcher.py ucsc.fa ensembl.fa > mappings.tsv
```
