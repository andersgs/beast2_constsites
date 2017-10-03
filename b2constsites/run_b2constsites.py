import click
from .b2constsites import B2ConstSites


@click.command()
@click.argument('seq')
@click.argument('vcf')
@click.option('-m', '--maskbed', default=None,
              help='A BED file with positions to mask')
@click.option('-f', '--fmt', default='fasta',
              help='seq sequence format', show_default=True)
def run_b2constsites(seq, vcf, maskbed, fmt):
    b2cs = B2ConstSites(seq=seq,
                        vcf=vcf,
                        maskbed=maskbed,
                        fmt=fmt)
    b2cs.load_seq()
    b2cs.load_vcf()
    b2cs.load_mask()
    b2cs.const_sites()
    print(b2cs)
    pass


if __name__ == '__main__':
    run_b2constsites()
