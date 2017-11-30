import click
import logging
from b2constsites import B2
from b2constsites import IQtree


@click.command()
@click.argument('tool', type=click.Choice(['beast2', 'iqtree']),
                default='beast2')
@click.argument('seq')
@click.argument('vcf')
@click.option('-x', '--xml', default=None,
              help='An XML file produced with BEAUTi')
@click.option('-m', '--maskbed', default=None,
              help='A BED file with positions to mask')
@click.option('-f', '--fmt', default='fasta',
              help='seq sequence format', show_default=True)
def run_b2constsites(tool, seq, vcf, xml, maskbed, fmt):
    '''
    Welcome to Beast2 Constant Sites.

    This version supports output for BEAST2 and IQTree.

    BEAST2 example:

        run_b2cs -x myxml.xml beast2 myref.fasta myvar.vcf

        NOTE: If running BEAST2, an XML file is **REQUIRED**.

    IQTree example:

        run_b2cs iqtree myref.fasta myvar.vcf

    If you encounter any issues, please submit them through:

        https://github.com/andersgs/beast2_constsites/issues
    '''
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s B2CONSTSITES %(levelname)-5s %(message)s',
        datefmt='%Y-%m-%d %H:%M')
    if tool == 'beast2':
        if xml is None:
            raise IOError("When using beast2 please provide an XML file.")
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
        logging.info(f'Created {b2cs.new_xml}.')
        logging.info(f'Use: beast <options> {b2cs.new_xml}')
        logging.info('DONE!')
    elif tool == 'iqtree':
        cs = IQtree(seq=seq,
                    vcf=vcf,
                    maskbed=maskbed,
                    fmt=fmt)
        cs.load_seq()
        cs.load_vcf()
        cs.load_mask()
        cs.const_sites()
        print(cs)
    else:
        logging.warning(f'I don\'t yet support {tool}.')
        logging.warning('Please file an issue on '
                        'https://github.com/andergs/beast2_constsites/issues')
        logging.warning('Exiting without doing anything!')



if __name__ == '__main__':
    run_b2constsites()
