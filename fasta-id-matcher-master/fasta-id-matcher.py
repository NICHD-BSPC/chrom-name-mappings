#!/usr/bin/env python

import sys
import argparse
import random
from pyfaidx import Fasta
import re
import hashlib

usage = """
Maps chromosome names from one fasta to the chromosome names of another.
Matches are identified by the full chromosome sequence converted
to uppercase. Mapped chromosomes are reported to stdout as TSV.
Any chromosomes unique to one fasta are reported to the two files, "a_not_b"
and "b_not_a", which will appear in the current directory (TODO: expose this
to CLI)
NOTE: this is intended to be used with FASTA files from the same assembly,
ideally with the same exact sequences but with different chromosome labels. It
will not work to match up chromosome IDs across different assemblies (UCSC hg19
and GENCODE hg38, for example).
"""



def construct_dict(fastarec, filterin, filterout):
    """ Build dictionaries where keys are hexdigest of md5sum and values are
    the chromosome ID from which the sequences came that passed 'filternames'.
    If an md5sum already exists, raise an error.
    Parameters
    ----------
    fastarec : pyfaidx.FastaRecord
        Fasta record from assembly, generated with function `Fasta`
    filterin : list
        List of chromosome names to be filtered in
    filterout : list
        List of chromosome names to be filtered out
    Returns
    -------
    Dictionary `fasta_d` where keys are hexdigest of md5sum and values are the
    chromosome IDs, and the list of chromosome names filtered out `out_fasta_head`
    """
    fasta_d = {}
    out_fasta_head = []
    for r in fastarec:
        (keep, out_fasta_head) = filternames(filterin, filterout, r.name, out_fasta_head)
        tmp_key = hashlib.md5(str(r).upper().encode()).hexdigest()
        if keep and tmp_key not in fasta_d.keys():
            fasta_d[tmp_key] = r.name
        elif keep and tmp_key in fasta_d.keys():
            raise ValueError(
                '{rname}\n and \n{fdtmpkey}\n have identical sequences\n'
                .format(rname=r.name, fdtmpkey=fasta_d[tmp_key]
                ))
    return (fasta_d, out_fasta_head)


def filternames(filterin, filterout, name, out_fasta_head):
    """
    Filters in or out according to the list of regular expressions
    in --filterin or --filterout
    Parameters
    ----------
    filterin : args.filterin_a or args_filterin_b
        Regular expression(s) to limit the chromosome names to
    filterout : args.filterout_a or args_filterout_b
        Regular expression(s) to filter out from the chromosome names
    name : str
        Chromosome name ot screen
    out_fasta_head : list
        List of chromosome names filtered out
    Returns
    -------
    True or False for keeping this chromosome name, and the updated list
    of filtered out chromosome names
    """
    if filterin is not None and filterin[0] != 'None':
        regex_in = re.compile(r'|'.join(filterin))
        match = regex_in.search(name)
        if match:
            in_keep = True
        else:
            in_keep = False
            out_fasta_head.append(name)
    else: in_keep = True
    if filterout is not None and filterout[0] != 'None':
        regex_out = re.compile('|'.join(filterout))
        match = regex_out.search(name)
        if match is None:
            out_keep = True
        else:
            out_keep = False
            out_fasta_head.append(name)
    else: out_keep = True

    return(in_keep and out_keep, out_fasta_head)


def difference(f1_d, f2_d, out_f1_d, out_f2_d):
    """
    Figures out the difference between two dictionaries
    and reports the difference
    Parameters
    ----------
    f1_d : dict
        Dictionary for first fasta, chromosome names as values
    f2_d : dict
        Dictionary for second fasta, chromosome names as values
    out_f1_head : list
        Dictionary of filtered out chromosomes from first fasta
    out_f2_head : list
        List of filtered out chromosomes from second fasta
    Returns
    -------
    The dictionaries of chromosome names found in one fasta and not the other
    nor in the filtered out chromosome names
    """
    f1_not_f2 = [
        f1_d[i]
        for i in set(list(f1_d.keys())).difference(list(f2_d.keys()))
    ]
    f1_not_f2_not_out = list(set(f1_not_f2) - set(list(out_f1_head)))

    f2_not_f1 = [
        f2_d[i]
        for i in set(list(f2_d.keys())).difference(list(f1_d.keys()))
    ]
    f2_not_f1_not_out = list(set(f2_not_f1) - set(list(out_f2_head)))
    return(f1_not_f2_not_out, f2_not_f1_not_out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(usage=usage)
    ap.add_argument('f1', help='first fasta')
    ap.add_argument('f2', help='second fasta')
    ap.add_argument('--samplesize', type=int, help='''Length of random sample to take from
                    each chromosome, default is %(default)s''', default=25)
    ap.add_argument('--nsamples', type=int, help='''Number of samples to take from each
                    chromosome, default is %(default)s''', default=3)
    ap.add_argument('--filterout_a', type=str, nargs='+', help='''String for filtering out chromosomes
                    from first fasta with a particular regular expression pattern''')
    ap.add_argument('--filterin_a', type=str, nargs='+', help='''Strings for keeping in chromosomes
                    from first fasta with a particular regular expression pattern''')
    ap.add_argument('--filterout_b', type=str, nargs='+', help='''String for filtering out chromosomes
                    from second fasta with a particular regular expression pattern''')
    ap.add_argument('--filterin_b', type=str, nargs='+', help='''Strings for keeping in chromosomes
                    from second fasta with a particular regular expression pattern''')
    ap.add_argument('--output', help='''Output file to create. Default is to
                    write to stdout.''')
    ap.add_argument('--out_a_not_b', help='''Output file a_not_b to create, default is %(default)s''',
                    default='a_not_b')
    ap.add_argument('--out_b_not_a', help='''Output file b_not_a to create, default is %(default)s''',
                    default='b_not_a')
    ap.add_argument('--out_filtered_a', help='''Output file for filtered-out chromosomes to create, default is %(default)s''',
                    default='filtered_out_a')
    ap.add_argument('--out_filtered_b', help='''Output file for filtered-out chromosomes to create, default is %(default)s''',
                    default='filtered_out_b')
    args = ap.parse_args()

    f1 = Fasta(args.f1, read_long_names=True)
    f2 = Fasta(args.f2, read_long_names=True)


    # Build dictionaries where keys are tuples (containing randomly-sampled
    # sequences) and values are the chromosome ID from which the sequences
    # came that passed 'filternames'.
    (f1_d, out_f1_head) = construct_dict(f1, args.filterin_a, args.filterout_a)
    (f2_d, out_f2_head) = construct_dict(f2, args.filterin_b, args.filterout_b)

    # Removes the > from fasta headers in the f1_d and f2_d dictionaries, coming from empty lines in the fastas
    for i in f1_d.keys():
        f1_d[i] = f1_d[i].lstrip(">")
    for i in f2_d.keys():
        f2_d[i] = f2_d[i].lstrip(">")

    # Same key in both dicts mean randomly-sampled sequences match . . . and we
    # have a mapping between chromosomes.
    if args.output is None:
        out = sys.stdout
    else:
        out = open(args.output, 'w')
    for k in f1_d.keys():
        if k in f2_d:
            print('{0}\t{1}'.format(f1_d[k].split(' ')[0], f2_d[k].split(' ')[0]), file=out)

    # Otherwise, figure out the differences and report to stderr.
    f1_not_f2_not_out, f2_not_f1_not_out = difference(f1_d, f2_d, out_f1_head, out_f2_head)

    with open(args.out_a_not_b, 'w') as fout:
        if f1_not_f2_not_out:
            fout.write('\n'.join(f1_not_f2_not_out) + '\n')

    with open(args.out_b_not_a, 'w') as fout:
        if f2_not_f1_not_out:
            fout.write('\n'.join(f2_not_f1_not_out) + '\n')

    with open(args.out_filtered_a, 'w') as fout:
        if out_f1_head:
            fout.write('\n'.join(out_f1_head) + '\n')

    with open(args.out_filtered_b, 'w') as fout:
        if out_f2_head:
            fout.write('\n'.join(out_f2_head) + '\n')
