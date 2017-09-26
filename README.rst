beast2_constsites: Create a line for XML with constant sites
============================================================

Introduction
------------

Based on this [suggestion](https://groups.google.com/forum/#!topic/beast-users/QfBHMOqImFE) by Remco,
we can correctly account for constant sites in a BEAST2 analysis by adding the
following ``data`` tag below your current ``data`` tag::

    <data id='xyz' spec='FilteredAlignment' filter='-' data='@xyzOriginal' constantSiteWeights='100 200 300 400'/>

This assumes that your original ``<data>`` tag had ``id=xyz`` and was renamed
to ``id=xyzOriginal``, and that you have 1000 constant sites that were
removed from the alignment, with:

*  100 As
*  200 Cs
*  300 Gs
*  400 Ts

What does this do?
------------------

This script will take a FASTA file with a single DNA sequence (e.g., a
bacterial chromosome), and a VCF file containing the position of
SNPs along the FASTA file (e.g., as outputted from [snippy-core](https://www.github.com/tseemann/snippy))
and will output the ``<data>`` tag ready to copy paste into your XML file.

It will optionally also take a BED file with positions to mask (e.g., positions
of phage).

How to install?
---------------

::

    pip3 install b2constsites
