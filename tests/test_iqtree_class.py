import pytest
import os
from b2constsites.b2constsites import IQtree


test_folder = os.path.dirname(os.path.realpath(__file__))
test_seq = os.path.join(test_folder, 'test.fa')
test_vcf = os.path.join(test_folder, 'test.vcf')
test_bed = os.path.join(test_folder, 'test.bed')
test_xml = os.path.join(test_folder, 'test.xml')

test_seq_withN = os.path.join(test_folder, 'test_withN.fa')
test_seq_withDash = os.path.join(test_folder, 'test_withdash.fa')

cs_bed = IQtree(seq=test_seq,
                vcf=test_vcf,
                maskbed=test_bed,
                fmt='fasta')

cs_nobed = IQtree(seq=test_seq,
                  vcf=test_vcf,
                  fmt='fasta')

cs_withN = IQtree(seq=test_seq_withN,
                  vcf=test_vcf,
                  fmt='fasta')

cs_withDash = IQtree(seq=test_seq_withDash,
                     vcf=test_vcf,
                     fmt='fasta')


def test_load_seq_with_bed():
    cs_bed.load_seq()
    assert len(cs_bed.seqrec) == 16


def test_load_vcf_with_bed():
    cs_bed.load_vcf()
    assert {2, 6, 10, 11} == cs_bed.var_pos


def test_load_bed_with_bed():
    cs_bed.load_mask()
    assert [range(4, 7)] == cs_bed.mask_pos


def test_const_sites_A_with_bed():
    cs_bed.const_sites()
    assert 2 == cs_bed.tally['A']


def test_const_sites_T_with_bed():
    cs_bed.const_sites()
    assert 2 == cs_bed.tally['T']


def test_const_sites_C_with_bed():
    cs_bed.const_sites()
    assert 2 == cs_bed.tally['C']


def test_const_sites_G_with_bed():
    cs_bed.const_sites()
    assert 4 == cs_bed.tally['G']


def test_output_with_bed(capfd):
    print(cs_bed)
    out, err = capfd.readouterr()
    assert out == 'iqtree <options> -fconst 2,2,4,2\n'


def test_load_seq_without_bed():
    cs_nobed.load_seq()
    assert len(cs_nobed.seqrec) == 16


def test_load_vcf_without_bed():
    cs_nobed.load_vcf()
    assert {2, 6, 10, 11} == cs_nobed.var_pos


def test_load_bed_without_bed():
    cs_nobed.load_mask()
    assert cs_nobed.mask_pos is None


def test_const_sites_A_without_bed():
    cs_nobed.const_sites()
    assert 3 == cs_nobed.tally['A']


def test_const_sites_T_without_bed():
    cs_nobed.const_sites()
    assert 3 == cs_nobed.tally['T']


def test_const_sites_C_without_bed():
    cs_nobed.const_sites()
    assert 2 == cs_nobed.tally['C']


def test_const_sites_G_without_bed():
    cs_nobed.const_sites()
    assert 4 == cs_nobed.tally['G']


def test_output_with_nobed(capfd):
    print(cs_nobed)
    out, err = capfd.readouterr()
    assert out == 'iqtree <options> -fconst 3,2,4,3\n'


def test_class_seq_withN(capfd):
    cs_withN.load_seq()
    cs_withN.load_vcf()
    cs_withN.load_mask()
    cs_withN.const_sites()
    print(cs_withN)
    out, err = capfd.readouterr()
    assert out == 'iqtree <options> -fconst 3,1,4,3\n'


def test_class_seq_withDash(capfd):
    cs_withDash.load_seq()
    cs_withDash.load_vcf()
    cs_withDash.load_mask()
    cs_withDash.const_sites()
    print(cs_withDash)
    out, err = capfd.readouterr()
    assert out == 'iqtree <options> -fconst 3,2,4,3\n'
