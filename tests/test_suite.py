import Bio
from pyfaidx import Fasta
import re
import hashlib
import sys
sys.path.append('./fasta-id-matcher-master')
from fasta_id_matcher import *
import pytest


def test_filternames_output():
    assert filternames(['fa1'], None, 'fa1', []) == (True, [])
    assert filternames(['fa1'], ['fa1'], 'fa1', []) == (False, ['fa1'])
    assert filternames(None, None, 'fa1', []) == (True, [])
    assert filternames(None, ['fa2'], 'fa1', []) == (True, [])
    assert filternames(None, 'fa1', 'fa1', []) == (False, ['fa1'])
    assert filternames(None, 'fa1', 'fa1', ['fa2']) == (False, ['fa2', 'fa1'])

def test_difference():
    f1_d = {'attaga': 'fa1', 'cccc': 'fa2', 'ctgat': 'fa3', 'aaaagtgc': 'fa4'}
    f2_d_all_match = {'attaga': 'FAS1', 'cccc': 'FAS2', 'ctgat': 'FAS3', 'aaaagtgc': 'FAS4'}
    f2_d_2_3_mismatch = {'attaga': 'FAS1', 'cctc': 'FAS2', 'caaat': 'FAS3', 'aaaagtgc': 'FAS4'}
    f2_2_no_match = {'agga': 'FAS1', 'cctcaaa': 'FAS2'}
    out_fa2_fa4 = ['fa2', 'fa4']
    out_FAS2 = ['FAS2']
    no_out = []
    assert difference(f1_d, f2_d_all_match, no_out, no_out) == ([], [])
    assert difference(f1_d, f2_d_all_match, out_fa2_fa4, out_FAS2) == ([], [])
    assert difference(f1_d, f2_d_2_3_mismatch, no_out, no_out) == (['fa2', 'fa3'], ['FAS2', 'FAS3'])
    assert difference(f1_d, f2_d_2_3_mismatch, out_fa2_fa4, out_FAS2) == (['fa3'], ['FAS3'])

def test_remove_char():
    assert remove_char({'atgatgaga': 'fa1'}) == {'atgatgaga': 'fa1'}
    assert remove_char({'atgatgaga': '>fa1'}) == {'atgatgaga': 'fa1'}
    assert remove_char({'atgatgaga': 'f>a1'}) == {'atgatgaga': 'f>a1'}
    assert remove_char({'atgatgaga': '>>fa1'}) == {'atgatgaga': 'fa1'}

def test_Fasta():
    f1 = Fasta('./tests/fasta_test1.fa', read_long_names=True)
    assert str(f1['sequence1 dna:chromosome REF']) == 'GTATTCTTGGACCTAAGCCTGGGAGTCTTTCAAGAAGGGATTTATAATATCTTTTAGTAGCTTAGTGTTACACTGAGCCACATAAAATTAATGTTTTTGAGAGTTAGAGGCAGC'

    with pytest.raises(ValueError) as e:
        Fasta('./tests/fasta_test2.fa', read_long_names=True)
    assert str(e.value) == 'Duplicate key "sequence2 dna:chromosome REF"'

    with pytest.raises(KeyError) as e:
        f1['unk_key']
    assert str(e.value) == "'unk_key not in ./tests/fasta_test1.fa.'"

def test_construct_dict():
    f3 = Fasta('./tests/fasta_test3.fa', read_long_names=True)
    f1 = Fasta('./tests/fasta_test1.fa', read_long_names=True)

    assert construct_dict(f3, None, None) == ({'f16b47f833397dbd6b406ee299ff0992': 'sequence1 dna:chromosome REF',
        'ad767ff4ba776ebd5ec39536dcf32ff1': 'sequence2 dna:chromosome REF'}, [])
    assert construct_dict(f3, None, ['sequence2']) == ({'f16b47f833397dbd6b406ee299ff0992':
        'sequence1 dna:chromosome REF'}, ['sequence2 dna:chromosome REF'])
    assert construct_dict(f3, ['sequence1'], None) == ({'f16b47f833397dbd6b406ee299ff0992':
        'sequence1 dna:chromosome REF'}, ['sequence2 dna:chromosome REF'])

    with pytest.raises(ValueError) as e:
        construct_dict(f1, None, None)
    assert str(e.value) == ">sequence3 dna:chromosome REF\n and \nsequence2 dna:chromosome REF\n have identical sequences\n"
