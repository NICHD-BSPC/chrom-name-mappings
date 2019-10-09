"""
Generate mapping file from fasta1 to fasta2 chromosome names
"""

import os
import sys
import gzip
import yaml
import pandas
from snakemake.utils import makedirs
from datetime import datetime


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
    resources: wget_limit = 1
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


rule unzip_and_grepheaders:
    """Unzip files and
    Output headers from  fasta file
    """
    input:
        rules.download_refs.output
    output:
        'mappings/{organism}/{label}/unzipped/{label}_ref_{name}.fa',
        'mappings/{organism}/{label}/{label}_headers_{name}.txt'
    run:
        shell('gunzip -c {input} > {output[0]}')
        shell('grep ">" {output[0]} > {output[1]}')


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

def get_names(output, organism, label):
    names = {}
    names['base'] = output.replace('.tsv','')
    names['name_a'] = config['references'][organism][label]['args']['from']
    names['name_b'] = config['references'][organism][label]['args']['to']
    names['a_not_b'] = f'{names["base"]}_{names["name_a"]}_not_{names["name_b"]}.txt'
    names['b_not_a'] = f'{names["base"]}_{names["name_b"]}_not_{names["name_a"]}.txt'
    names['filtered_a'] = f'{names["base"]}_filtered_out_{names["name_a"]}.txt'
    names['filtered_b'] = f'{names["base"]}_filtered_out_{names["name_b"]}.txt'
    names['filterin_a'] = get_filter_params(organism, label, names['name_a'], 'filterin')
    names['filterout_a'] = get_filter_params(organism, label, names['name_a'], 'filterout')
    names['filterin_b'] = get_filter_params(organism, label, names['name_b'], 'filterin')
    names['filterout_b'] = get_filter_params(organism, label, names['name_b'], 'filterout')
    return names

def check_line_sum(organism, label, name, mapping, filtered):
    headers_line = sum(1 for line in open(f'mappings/{organism}/{label}/{label}_headers_{name}.txt'))
    mappings_line = sum(1 for line in open(mapping))
    filtered_line = sum(1 for line in open(filtered))
    return (headers_line - mappings_line - filtered_line)


rule map_ids:
    """ Map chromosome names from fasta1 to fasta2 reference
    """
    input:
        get_unzipfiles
    output:
        'mappings/{organism}/{label}/mappings_{label}.tsv'
    run:
        names = get_names(output[0], wildcards.organism, wildcards.label)
        shell(
            "python3 fasta-id-matcher-master/fasta-id-matcher.py "
            "{input[0]} "
            "{input[1]} "
            "--filterin_a {names[filterin_a]} "
            "--filterout_a {names[filterout_a]} "
            "--filterin_b {names[filterin_b]} "
            "--filterout_b {names[filterout_b]} "
            "--output {output} "
            "--out_a_not_b {names[a_not_b]} "
            "--out_b_not_a {names[b_not_a]} "
            "--out_filtered_a {names[filtered_a]} "
            "--out_filtered_b {names[filtered_b]} "
            )

        # checks for chromosome name in one assembly non-matching or absent from the other assembly
        if os.stat(names['a_not_b']).st_size != 0 or os.stat(names['b_not_a']).st_size != 0:
            linelist_a =[]
            linelist_b = []
            with open(names['a_not_b'],'rt') as tmp:
                for line in tmp:
                    linelist_a.append(line)
            with open(names['b_not_a'],'rt') as tmp:
                for line in tmp:
                    linelist_b.append(line)
            raise ValueError(
                '{label}: from {a} not found in {b}:\n{linelist_a}\nfrom {b} not found in {a}:\n{linelist_b}'
                .format(label=wildcards.label, linelist_a=linelist_a, linelist_b=linelist_b,
                a=config['references'][wildcards.organism][wildcards.label]['args']['from'],
                b=config['references'][wildcards.organism][wildcards.label]['args']['to'])
                )

        # calculate if the number of chromosome names in the fasta files match
        # the number of chromosome names mapped + filtered out
        if check_line_sum(wildcards.organism, wildcards.label, names['name_a'], output[0], names['filtered_a']) != 0:
            with open('log_number_mismatch.txt', "a") as fout:
                fout.write(
                '{date}\nNumber of lines in {name_a} {label} does not match sum of mapped + filtered out; {i} missing\n'
                .format(name_a = names['name_a'], label = wildcards.label,
                i = check_line_sum(wildcards.organism, wildcards.label, names['name_a'], output[0], names['filtered_a']),
                date = datetime.now())
                )

        if check_line_sum(wildcards.organism, wildcards.label, names['name_b'], output[0], names['filtered_b']) != 0:
            with open('log_number_mismatch.txt', "a") as fout:
                fout.write(
                '{date}\nNumber of lines in {name_b} {label} does not match sum of mapped + filtered out; {i} missing\n'
                .format(name_b = names['name_b'], label = wildcards.label,
                i = check_line_sum(wildcards.organism, wildcards.label, names['name_b'], output[0], names['filtered_b']),
                date = datetime.now())
                )
        # cleanout empty output files
        for f in (names['a_not_b'], names['b_not_a'], names['filtered_a'], names['filtered_b']):
            if os.stat(f).st_size == 0:
                shell('rm {f}')
