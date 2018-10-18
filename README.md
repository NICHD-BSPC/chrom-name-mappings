# chrom-name-mappings
Automated workflow for matching chromosome names between genome assemblies.

This worflow attempts to match the chromosomes names from two reference genome fasta files by comparing the size of chromosomes and randomly sampling the fasta sequences to identify matches. It returns a .tsv table of matched pairs of chromosome names, along with text files containing the lists of filtered-out chromosome names and non-matching chromosome names, if applicable.

## Software installation

Installation of all dependencies is handled by conda, ensuring reproducibility, streamlined setup, and no need for root administrator privileges.

Use [bioconda](https://bioconda.github.io/) to automatically install software into the working directory without needing admin rights on the machine.

If you have the Anaconda Python distribution, you already have conda. Otherwise, install [Miniconda](https://conda.io/miniconda.html).

### 1. Clone the git repo

Clone the repository from github into a new directory and change to that directory.
```
git clone https://github.com/NICHD-BSPC/chrom-name-mappings.git my-project-dir
cd my-project-dir
```

### 2. Create a new conda environment

Create a top-level environment with Snakemake and other requirements. It needs to be activated any time you’ll be working with these workflows.
```
conda create -n chrom-name-mappings --file requirements.txt --channel bioconda --channel conda-forge
```

Then activate the environment:
```
source activate chrom-name-mappings
```

Eventually when you’re done, you can “deactivate”, which removes the environment location from your $PATH until the next time you activate it.
```
source deactivate
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
        fasta1:
          url: 'url/to/reference/genome/fasta1.gz'
          filterin: None
          filterout: None
        fasta2:
          url: 'url/to/reference/genome/fasta1.gz'
          filterin: None
          filterout: None
      args:
          from: 'fasta1'
          to: 'fasta2'
```
The code accepts paths to archives .gz or .tar.gz as urls.

If the url points to a local archive, indicate `url: 'file://path/to/archive'`

The argument `filterin` takes a regular expression, a list, or a file containing regular expression, to limit the mapping to chromosome names containing the regular expression(s).

The argument `filterout` takes a regular expression, a list, or a file containing regular expression, to discard from the mapping to chromosome names containing the regular expression(s).

When both `filterin` and `filterout` arguments are indicated, the chromosome names will be filtered first according to `filterin` then `filterout`.

The DAG of jobs looks like this:

![DAG_chrom-name-mapings](/dag.svg)
