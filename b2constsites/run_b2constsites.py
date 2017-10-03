import click
import logging
from .b2constsites import B2ConstSites


@click.command()
@click.argument('seq')
@click.argument('vcf')
@click.argument('xml')
@click.option('-m', '--maskbed', default=None,
              help='A BED file with positions to mask')
@click.option('-f', '--fmt', default='fasta',
              help='seq sequence format', show_default=True)
def run_b2constsites(seq, vcf, xml, maskbed, fmt):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s B2CONSTSITES %(levelname)-5s %(message)s',
        datefmt='%Y-%m-%d %H:%M')
    b2cs = B2ConstSites(seq=seq,
                        vcf=vcf,
                        xml=xml,
                        maskbed=maskbed,
                        fmt=fmt)
    b2cs.load_seq()
    b2cs.load_xml()
    b2cs.rename_original_data_tag()
    b2cs.load_vcf()
    b2cs.load_mask()
    b2cs.const_sites()
    b2cs.create_new_xml()
    logging.info(f'Created {b2cs.new_xml}.')
    logging.info(f'Use: beast <options> {b2cs.new_xml}')
    logging.info('DONE!')
    pass


if __name__ == '__main__':
    run_b2constsites()
