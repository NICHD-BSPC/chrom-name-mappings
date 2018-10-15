"""
Generate mapping file from Ensembl to UCSC chromosome names
"""

import os
import sys
import gzip
import yaml
import pandas
from snakemake.utils import makedirs


config.update(yaml.load(open('config/config.yaml')))

# reads that config, and for each pair of URLs downloads, does the fasta ID matching, and writes out the mappings file
# print(config['references']['BDGP6_to_UCSC_dm6'])

targets = []
for organism in config['references']:
    for label in config['references'][organism]:
        list_of_fastas = config['references'][organism][label]
        aname = list(list_of_fastas[0].keys())[0]
        bname = list(list_of_fastas[1].keys())[0]
        targets.extend(
            [
                f'mappings/{organism}/{label}/gz/',
                f'mappings/{organism}/{label}/unziped/',
                f'mappings/{organism}/{label}/map/'
            ]
        )

rule all:
     input: targets


def get_fastas(wildcards):
    url_a = list(config['references'][wildcards.organism][wildcards.label][0].values())[0]
    url_b = list(config['references'][wildcards.organism][wildcards.label][1].values())[0]
    return (url_a, url_b)

def get_gzfiles(wildcards):
    gzfile_a = 'mappings/'+wildcards.organism+'/'+wildcards.label+'/gz/'+wildcards.label+'_ref_' + list(config['references'][wildcards.organism][wildcards.label][0].keys())[0] + '.gz'
    gzfile_b = 'mappings/'+wildcards.organism+'/'+wildcards.label+'/gz/'+wildcards.label+'_ref_' + list(config['references'][wildcards.organism][wildcards.label][1].keys())[0] + '.gz'
    return (gzfile_a, gzfile_b)

def get_unzipfiles(wildcards):
    file_a = 'mappings/'+wildcards.organism+'/'+wildcards.label+'/unziped/'+wildcards.label+'_ref_' + list(config['references'][wildcards.organism][wildcards.label][0].keys())[0] + '.fa'
    file_b = 'mappings/'+wildcards.organism+'/'+wildcards.label+'/unziped/'+wildcards.label+'_ref_' + list(config['references'][wildcards.organism][wildcards.label][1].keys())[0] + '.fa'
    return (file_a, file_b)

rule download_refs:
    """ Download urls for each pair of Ensembl / UCSC references
    """
    params:
        get_fastas,
        get_gzfiles
    output:
        #temp('mappings/{organism}/{label}/gz/')
        'mappings/{organism}/{label}/gz/'
    run:
        if params[0][0].endswith("tar.gz") == True:
            shell("wget -qO- {params[0][0]} | tar -xOz | gzip > {params[1][0]}")
        else:
            shell("wget {params[0][0]} -O {params[1][0]}")
        if params[0][1].endswith("tar.gz") == True:
            shell("wget -qO- {params[0][1]} | tar -xOz | gzip > {params[1][1]}")
        else:
            shell("wget {params[0][1]} -O {params[1][1]}")


rule unzip:
    """Unzip files
    """
    input:
        rules.download_refs.output
    params:
        get_gzfiles,
        get_unzipfiles
    output:
        #temp('mappings/{organism}/{label}/unziped/')
        'mappings/{organism}/{label}/unziped/'
    run:
        shell("gunzip -c {params[0][0]} > {params[1][0]}")
        shell("gunzip -c {params[0][1]} > {params[1][1]}")


rule map_ensembl_to_ucsc:
    """ Map chromosome names from Ensembl to UCSC reference
    """
    input:
        rules.unzip.output
    params:
        get_unzipfiles
    output:
        'mappings/{organism}/{label}/map/'
    run:
        mapping = 'mappings/'+wildcards.organism+'/'+wildcards.label+'/map/mappings_'+wildcards.label+'.tsv'
        a_not_b = ('mappings/'+wildcards.organism+'/'+wildcards.label+'/map/'
            +list(config['references'][wildcards.organism][wildcards.label][0].keys())[0]+'_not_'
            +list(config['references'][wildcards.organism][wildcards.label][1].keys())[0]+'.tsv')
        b_not_a = ('mappings/'+wildcards.organism+'/'+wildcards.label+'/map/'
            +list(config['references'][wildcards.organism][wildcards.label][1].keys())[0]+'_not_'
            +list(config['references'][wildcards.organism][wildcards.label][0].keys())[0]+'.tsv')
        filtered = 'mappings/'+wildcards.organism+'/'+wildcards.label+'/map/filtered_out.tsv'
        filterout = config['references'][wildcards.organism][wildcards.label][2]['args']['filterout']
        filterin = config['references'][wildcards.organism][wildcards.label][2]['args']['filterin']

        shell("python3 fasta-id-matcher-master/fasta-id-matcher.py {params[0][0]} {params[0][1]} --filterout {filterout} --filterin {filterin} --output {mapping} --out_a_not_b {a_not_b} --out_b_not_a {b_not_a} --out_filtered {filtered}")
        if os.stat(a_not_b).st_size != 0:
            with open(a_not_b,'rt') as tmp:
                linelist = 'Chromosome(s) "'
                for line in tmp:
                    linelist = linelist + str(line) + '\t'
            raise ValueError(
                '{label}: {linelist}" from {a} not found in {b}'
                .format(label=wildcards.label, linelist=linelist, a=list(config['references'][wildcards.organism][wildcards.label][0].keys())[0],
                b=list(config['references'][wildcards.organism][wildcards.label][1].keys())[0]))
