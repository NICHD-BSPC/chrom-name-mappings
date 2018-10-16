"""
Generate mapping file from fasta1 to fasta2 chromosome names
"""

import os
import sys
import gzip
import yaml
import pandas
from snakemake.utils import makedirs


targets = []
for organism in config['references']:
    for label in config['references'][organism]:
        for name in config['references'][organism][label]['fastas']:
            targets.append(f'mappings/{organism}/{label}/unzipped/{label}_ref_{name}.fa')
            targets.append(f'mappings/{organism}/{label}/mappings_{label}.tsv')
            targets.append(f'mappings/{organism}/{label}/{label}_headers_{name}.txt')


rule all:
     input: targets


rule download_refs:
    """ Download urls for each pair of fasta1 / fasta2 references
    """
    output:
        'mappings/{organism}/{label}/gz/{label}_ref_{name}.fa.gz'
    run:
        url = (
            config['references'][wildcards.organism][wildcards.label]
            ['fastas'][wildcards.name]['url']
        )
        if url.startswith(('http://', 'https://', 'ftp://')):
            if url.endswith('.tar.gz'):
                shell('wget -qO- {url} | tar -xOz | gzip -c > {output}')
            else:
                shell('wget -qO- {url} > {output}')
        else:
            if url.endswith('.tar.gz'):
                shell('tar -xOz -f {url} | gzip > {output}')
            else:
                shell('cp {url} {output}')


rule unzip:
    """Unzip files
    """
    input:
        rules.download_refs.output
    output:
        'mappings/{organism}/{label}/unzipped/{label}_ref_{name}.fa'
    shell:
        'gunzip -c {input} > {output}'


rule grep_headers:
    """ Output headers from  fasta file
    """
    input:
        rules.unzip.output
    output:
        'mappings/{organism}/{label}/{label}_headers_{name}.txt'
    shell:
        'grep ">" {input} > {output}'


def get_unzipfiles(wildcards):
    name_a = config['references'][wildcards.organism][wildcards.label]['args']['from']
    name_b = config['references'][wildcards.organism][wildcards.label]['args']['to']
    file_a = f'mappings/{wildcards.organism}/{wildcards.label}/unzipped/{wildcards.label}_ref_{name_a}.fa'
    file_b = f'mappings/{wildcards.organism}/{wildcards.label}/unzipped/{wildcards.label}_ref_{name_b}.fa'
    return (file_a, file_b)

def get_filter_params(organism, label, name, filter):
    item = (config['references'][organism][label]['fastas'][name][filter])
    if type(item) is str and item.startswith('file://'):
        item = item.replace('file://', '')
        fil = open(item, 'r').readlines()
        fil = [w.replace('\n', '') for w in fil]
    else:
        fil = item
    return fil


rule map_ids:
    """ Map chromosome names from fasta1 to fasta2 reference
    """
    input:
        get_unzipfiles
    output:
        'mappings/{organism}/{label}/mappings_{label}.tsv'
    run:
        base = output[0].replace('.tsv','')
        name_a = config['references'][wildcards.organism][wildcards.label]['args']['from']
        name_b = config['references'][wildcards.organism][wildcards.label]['args']['to']
        a_not_b = f'{base}_{name_a}_not_{name_b}.txt'
        b_not_a = f'{base}_{name_b}_not_{name_a}.txt'
        filtered_a = f'{base}_filtered_out_{name_a}.txt'
        filtered_b = f'{base}_filtered_out_{name_b}.txt'
        organism = wildcards.organism
        label = wildcards.label
        filterin_a = get_filter_params(organism, label, name_a, 'filterin')
        filterout_a = get_filter_params(organism, label, name_a, 'filterout')
        filterin_b = get_filter_params(organism, label, name_b, 'filterin')
        filterout_b = get_filter_params(organism, label, name_b, 'filterout')

        shell(
            "python3 fasta-id-matcher-master/fasta-id-matcher.py "
            "{input[0]} "
            "{input[1]} "
            "--filterin_a {filterin_a} "
            "--filterout_a {filterout_a} "
            "--filterin_b {filterin_b} "
            "--filterout_b {filterout_b} "
            "--output {output} "
            "--out_a_not_b {a_not_b} "
            "--out_b_not_a {b_not_a} "
            "--out_filtered_a {filtered_a} "
            "--out_filtered_b {filtered_b} "
            )


        if os.stat(a_not_b).st_size != 0 or os.stat(b_not_a).st_size != 0:
            linelist_a =[]
            linelist_b = []
            with open(a_not_b,'rt') as tmp:
                for line in tmp:
                    linelist_a.append(line)
            with open(b_not_a,'rt') as tmp:
                for line in tmp:
                    linelist_b.append(line)
            raise ValueError(
                '{label}: from {a} not found in {b}:\n{linelist_a}\nfrom {b} not found in {a}:\n{linelist_b}'
                .format(label=wildcards.label, linelist_a=linelist_a, linelist_b=linelist_b,
                a=config['references'][wildcards.organism][wildcards.label]['args']['from'],
                b=config['references'][wildcards.organism][wildcards.label]['args']['to'])
                )
