# chrom-name-mappings
Automated workflow for matching chromosome names between genome assemblies.

This worflow attempts to match the chromosomes names from two reference genome fasta files by comparing the full chromosome sequences between
fastas to identify matches. It returns a .tsv table of matched pairs of chromosome names, along with text files containing the lists of filtered-out
chromosome names and non-matching chromosome names, if applicable.

## Software installation

Installation of all dependencies is handled by conda or mamba, ensuring reproducibility, streamlined setup, and no need for root administrator privileges.

Use [bioconda](https://bioconda.github.io/) to automatically install software into the working directory without needing admin rights on the machine.

If you have the Anaconda Python distribution, you already have conda. Otherwise, install [Miniconda](https://conda.io/miniconda.html) or
[Mamba](https://mamba.readthedocs.io/en/latest/index.html).

### 1. Clone the git repo

Clone the repository from github into a new directory and change to that directory.
```
git clone https://github.com/NICHD-BSPC/chrom-name-mappings.git my-project-dir
cd my-project-dir
```

### 2. Create a new conda environment

Create a local environment in your project directory with Snakemake and other requirements. It needs to be activated any time you’ll be working with these workflows.
```
mamba env create -p env --file requirements.yaml
```

Then navigate to your project directory and activate the environment:
```
source activate env/
```

Eventually when you’re done, you can “deactivate”, which removes the environment location from your $PATH until the next time you activate it.
```
source deactivate
```

## Running with test data

Test data are included, run with the following (assuming 20 cores; adjust this value as needed):

```
set -e

snakemake --snakefile Snakefile --configfile=config/test-config.yaml --resources wget_limit=2 -j 20
```


## Running the workflow

Example usage:
```
snakemake --snakefile Snakefile --configfile=config/config.yaml
```

The code can be executed in parralel by setting `-j` to a value greater than 1. In this case, it is recommended to limit the concurent downloads of reference files by setting the number of `--resources wget_limit` to 2.
```
snakemake --snakefile Snakefile --configfile=config/config.yaml -j 20 --resources wget_limit=2
```

Edit config/config.yaml for different genome assemblies and filtering options. The file config.yaml must follow the format:
```
references:
  organism:
    label:
      fastas:
        assembly1:
          url: 'url/to/reference/genome/fasta1.gz'
          filterin: None
          filterout: None
        assembly2:
          url: 'url/to/reference/genome/fasta1.gz'
          filterin: None
          filterout: None
      args:
          from: 'assembly1'
          to: 'assembly2'
```
The code accepts paths to archives .gz or .tar.gz as urls.

If the url points to a local archive, indicate `url: 'file://path/to/archive'`

The argument `filterin` takes a regular expression, a list, or a file containing regular expression, to limit the mapping to chromosome names containing the regular expression(s).

The argument `filterout` takes a regular expression, a list, or a file containing regular expression, to discard from the mapping to chromosome names containing the regular expression(s).

When both `filterin` and `filterout` arguments are indicated, the chromosome names will be included only when there are being both filtered in and not filtered out.

Chromosomes from one assembly that cannot be mapped to the other assembly are indicated in the output file `mappings_{label}_{assembly1}_not_{assembly2}.txt`. It is recommended to add the unmapped chromosomes to the `filterout` parameter before running the workflow again.

The output file `log_number_mismatch` reports mismatches between the number of chromosomes in an assembly and the number of mapped chromosome names + the filtered out chromosomes. This can be due to 2 or more identical chromosomes in an assembly. In this case, only one of the identical chromosomes will be reported in the mapping table.

The DAG of jobs looks like this:

![DAG_chrom-name-mapings](/dag.svg)
