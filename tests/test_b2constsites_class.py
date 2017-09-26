import pytest
import os
from b2constsites.b2constsites import B2ConstSites


test_folder = os.path.dirname(os.path.realpath(__file__))
test_seq = os.path.join(test_folder, 'test.fa')
test_vcf = os.path.join(test_folder, 'test.vcf')
test_bed = os.path.join(test_folder, 'test.bed')

b2cs_bed = B2ConstSites(seq=test_seq,
                        vcf=test_vcf,
                        maskbed=test_bed,
                        fmt='fasta')


b2cs_nobed = B2ConstSites(seq=test_seq,
                          vcf=test_vcf,
                          fmt='fasta')


def test_load_seq_with_bed():
    b2cs_bed.load_seq()
    assert len(b2cs_bed.seqrec) == 16


def test_load_vcf_with_bed():
    b2cs_bed.load_vcf()
    assert {2, 6} == b2cs_bed.var_pos


def test_load_bed_with_bed():
    b2cs_bed.load_mask()
    assert [range(4, 7)] == b2cs_bed.mask_pos


def test_const_sites_A_with_bed():
    b2cs_bed.const_sites()
    assert 2 == b2cs_bed.tally['A']


def test_const_sites_T_with_bed():
    b2cs_bed.const_sites()
    assert 2 == b2cs_bed.tally['T']


def test_const_sites_C_with_bed():
    b2cs_bed.const_sites()
    assert 4 == b2cs_bed.tally['C']


def test_const_sites_G_with_bed():
    b2cs_bed.const_sites()
    assert 4 == b2cs_bed.tally['G']


def test_output_with_bed(capfd):
    print(b2cs_bed)
    out, err = capfd.readouterr()
    assert out == "<data id='xyz' spec='FilteredAlignment' filter='-'"\
                  " data='@xyzOriginal'"\
                  " constantSiteWeights='2 4"\
                  " 4 2'/>\n"


def test_load_seq_without_bed():
    b2cs_nobed.load_seq()
    assert len(b2cs_nobed.seqrec) == 16


def test_load_vcf_without_bed():
    b2cs_nobed.load_vcf()
    assert {2, 6} == b2cs_nobed.var_pos


def test_load_bed_without_bed():
    b2cs_nobed.load_mask()
    assert b2cs_nobed.mask_pos is None


def test_const_sites_A_without_bed():
    b2cs_nobed.const_sites()
    assert 3 == b2cs_nobed.tally['A']


def test_const_sites_T_without_bed():
    b2cs_nobed.const_sites()
    assert 3 == b2cs_nobed.tally['T']


def test_const_sites_C_without_bed():
    b2cs_nobed.const_sites()
    assert 4 == b2cs_nobed.tally['C']


def test_const_sites_G_without_bed():
    b2cs_nobed.const_sites()
    assert 4 == b2cs_nobed.tally['G']
