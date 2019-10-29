import click
import logging
from b2constsites.b2constsites import B2


@click.command()
@click.argument('seq')
@click.argument('vcf')
@click.argument('xml')
@click.option('-m', '--maskbed', default=None,
              help='A BED file with positions to mask')
@click.option('-f', '--fmt', default='fasta',
              help='seq sequence format', show_default=True)
def b2c_beast(seq, vcf, xml, maskbed, fmt):
    '''
    Welcome to Beast2 Constant Sites.

    Example usage:

        b2c-beast myref.fasta myvar.vcf myxml.xml

    If you encounter any issues, please submit them through:

        https://github.com/andersgs/beast2_constsites/issues
    '''
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s B2C-BEAST %(levelname)-5s %(message)s',
        datefmt='%Y-%m-%d %H:%M')
    cs = B2(seq=seq,
            vcf=vcf,
            xml=xml,
            maskbed=maskbed,
            fmt=fmt)
    cs.load_seq()
    cs.load_xml()
    cs.rename_original_data_tag()
    cs.load_vcf()
    cs.load_mask()
    cs.const_sites()
    cs.create_new_xml()
    logging.info(f'Created {cs.new_xml}.')
    logging.info(f'Use: beast <options> {cs.new_xml}')
    logging.info('DONE!')

@click.command()
@click.argument('seq')
@click.argument('vcf')
@click.option('-m', '--maskbed', default=None,
              help='A BED file with positions to mask')
@click.option('-f', '--fmt', default='fasta',
              help='seq sequence format', show_default=True)
def b2c_iqtree():
    """
    Welcome IQTree constant sites
    
    Example usage:
    
        b2c-iqtree iqtree myref.fasta myvar.vcf
        
        If you encounter any issues, please submit them through:

        https://github.com/andersgs/beast2_constsites/issues
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s B2C-IQTREE %(levelname)-5s %(message)s',
        datefmt='%Y-%m-%d %H:%M')

    cs = IQtree(seq=seq,
                vcf=vcf,
                maskbed=maskbed,
                fmt=fmt)
    cs.load_seq()
    cs.load_vcf()
    cs.load_mask()
    cs.const_sites()
    print(cs)
    logging.info('DONE!')


if __name__ == '__main__':
    b2c_beast()
