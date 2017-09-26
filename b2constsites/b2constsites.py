import vcf
from Bio import SeqIO
import os
import logging
import collections


class B2ConstSites:
    def __init__(self,
                 seq,
                 vcf,
                 maskbed=None,
                 fmt='fasta'):
        self.seq = os.path.realpath(seq)
        self.seq_format = fmt
        self.vcf = os.path.realpath(vcf)
        if maskbed is not None:
            self.mask = os.path.realpath(maskbed)
        else:
            self.mask = None

    def load_seq(self):
        fh = open(self.seq, 'rt')
        self.seqrec = SeqIO.read(fh, format=self.seq_format)
        logging.info(f'Imported sequence {self.seqrec.id}.')
        logging.info(f'it has length {len(self.seqrec)}')
        fh.close()
        pass

    def load_vcf(self):
        fh = open(self.vcf, 'rt')
        self.var_pos = {rec.POS for rec in vcf.Reader(fh)}
        fh.close()
        logging.info(f'Loaded {len(self.var_pos)} SNPs from vcf.')
        pass

    def load_mask(self):
        if self.mask is not None:
            fh = open(self.mask, 'rt')
            mask_pos = [line.strip().split('\t') for line in fh]
            self.mask_pos = [range(int(start)+1, int(end))
                             for chrom, start, end in mask_pos]
        else:
            self.mask_pos = None
        pass

    def const_sites(self):
        tally = collections.Counter({'A': 0,
                                     'C': 0,
                                     'G': 0,
                                     'T': 0})
        for p, b in enumerate(self.seqrec):
            p += 1
            if self.mask_pos is not None\
                    and any([p in r for r in self.mask_pos])\
                    or p in self.var_pos:
                continue
            else:
                tally[b] += 1
        self.tally = tally
        pass

    def __str__(self):
        return f"<data id='xyz' spec='FilteredAlignment' filter='-'"\
               f" data='@xyzOriginal'"\
               f" constantSiteWeights='{self.tally['A']} {self.tally['C']}"\
               f" {self.tally['G']} {self.tally['T']}'/>"
