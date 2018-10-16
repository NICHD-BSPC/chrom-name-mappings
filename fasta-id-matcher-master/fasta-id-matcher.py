#!/usr/bin/env python

import sys
import argparse
import random
from pyfaidx import Fasta
import re

usage = """
Maps chromosome names from one fasta to the chromosome names of another.

Matches are identified by the combination of length and three 25-mers,
converted to uppercase and randomly sampled from across the sequence. The
number and size of random samples can be adjusted.

Mapped chromosomes are reported to stdout as TSV.

Any chromosomes unique to one fasta are reported to the two files, "a_not_b"
and "b_not_a", which will appear in the current directory (TODO: need better
names, better places, expose this to CLI)

NOTE: this is intended to be used with FASTA files from the same assembly,
ideally with the same exact sequences but with different chromosome labels. It
will not work to match up chromosome IDs across different assemblies (UCSC hg19
and GENCODE hg38, for example).
"""


def sample(rec, sample_size):
    """
    Randomly sample from a record.

    Parameters
    ----------

    rec : pyfaidx.FastaRecord

    sample_size : int
        Randomly sample a sequence of this length from the record

    Returns
    -------
    Randomly-sampled sequence of length `sample_size` as a string.
    """

    s = len(rec)

    # Even for very small seqs or large sample_size, going negative is OK.
    start = random.randint(0, s - sample_size - 1)
    seq = rec[start:start + sample_size]
    seq = seq.seq.upper()
    return seq


def gen_key(rec, nsamples=3, sample_size=25):
    """
    Returns a tuple of randomly-sampled sequences from `rec`.

    Parameters
    ----------
    rec : pyfaidx.FastaRecord

    nsamples : int
        Number of random samples to take

    Returns
    -------
    A tuple of length nsamples + 1 of the form (len(rec), seq1, ... seqN) where
    seq1 thru seqN are sequences of length `sample_size` and there are
    `nsamples` sequences.
    """

    # The first cut for deciding if chroms are equivalent is that they have the
    # same length. If that's true, then we want to sample at similar locations
    # across them (that is, the same sequence provided as input to two
    # different calls of this function). So we set the seed to the chrom
    # length.
    s = len(rec)
    random.seed(s)

    key = [s]
    for _ in range(nsamples):
        key.append(sample(rec, sample_size))

    return tuple(key)


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

    # Build dictionaries where keys are tuples from `gen_key()` (containing
    # randomly-sampled sequences) and values are the chromosome ID from which
    # the sequences came.
    f1_d = {gen_key(r, nsamples=args.nsamples, sample_size=args.samplesize): r.name for r in f1}
    f2_d = {gen_key(r, nsamples=args.nsamples, sample_size=args.samplesize): r.name for r in f2}

    #Filters in or out according to the list of regular expressions in --filterin or --filterout
    out_f1_d = {}
    fil_f1_d = {}
    out_f2_d = {}
    fil_f2_d = {}

    if args.filterin_a is not None and args.filterin_a[0] != 'None':
        regex_in = re.compile(r'|'.join(args.filterin_a))
        for i in f1_d.keys():
            match = regex_in.search(f1_d[i])
            if match:
                fil_f1_d[i] = f1_d[i]
            else:
                out_f1_d[i] = f1_d[i]
        f1_d = fil_f1_d

    fil_f1_d = {}
    if args.filterout_a is not None and args.filterout_a[0] != 'None':
        regex_out = re.compile('|'.join(args.filterout_a))
        for i in f1_d.keys():
            match = regex_out.search(f1_d[i])
            if match is None:
                fil_f1_d[i] = f1_d[i]
            else:
                out_f1_d[i] = f1_d[i]
        f1_d = fil_f1_d

    if args.filterin_b is not None and args.filterin_b[0] != 'None':
        regex_in = re.compile(r'|'.join(args.filterin_b))
        for i in f2_d.keys():
            match = regex_in.search(f2_d[i])
            if match:
                fil_f2_d[i] = f2_d[i]
            else:
                out_f2_d[i] = f2_d[i]
        f2_d = fil_f2_d

    fil_f2_d = {}
    if args.filterout_b is not None and args.filterout_b[0] != 'None':
        regex_out = re.compile('|'.join(args.filterout_b))
        for i in f2_d.keys():
            match = regex_out.search(f2_d[i])
            if match is None:
                fil_f2_d[i] = f2_d[i]
            else:
                out_f2_d[i] = f2_d[i]
        f2_d = fil_f2_d

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
    f1_not_f2 = [
        f1_d[i]
        for i in set(list(f1_d.keys())).difference(list(f2_d.keys()))
    ]
    f1_not_f2_not_out = list(set(f1_not_f2) - set(list(out_f1_d.values())))

    f2_not_f1 = [
        f2_d[i]
        for i in set(list(f2_d.keys())).difference(list(f1_d.keys()))
    ]
    f2_not_f1_not_out = list(set(f2_not_f1) - set(list(out_f2_d.values())))

    with open(args.out_a_not_b, 'w') as fout:
        if f1_not_f2_not_out:
            fout.write('\n'.join(f1_not_f2_not_out) + '\n')

    with open(args.out_b_not_a, 'w') as fout:
        if f2_not_f1_not_out:
            fout.write('\n'.join(f2_not_f1_not_out) + '\n')

    with open(args.out_filtered_a, 'w') as fout:
        if out_f1_d:
            fout.write('\n'.join(out_f1_d.values()) + '\n')

    with open(args.out_filtered_b, 'w') as fout:
        if out_f2_d:
            fout.write('\n'.join(out_f2_d.values()) + '\n')
