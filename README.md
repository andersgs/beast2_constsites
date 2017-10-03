# beast2_constsites: Create a line for BEAST2 XML with constant sites

[![Build Status](https://travis-ci.org/andersgs/beast2_constsites.svg?branch=master)](https://travis-ci.org/andersgs/beast2_constsites)

## Introduction

Based on this [suggestion](https://groups.google.com/forum/#!topic/beast-users/QfBHMOqImFE) by Remco,
we can correctly account for constant sites in a BEAST2 analysis by adding the
following ``data`` tag below your current ``data`` tag:

    <data id='xyz' spec='FilteredAlignment' filter='-' data='@xyzOriginal' constantSiteWeights='100 200 300 400'/>

This assumes that your original ``<data>`` tag had ``id=xyz`` and was renamed
to ``id=xyzOriginal``, and that you have 1000 constant sites that were
removed from the alignment, with:

*  100 As
*  200 Cs
*  300 Gs
*  400 Ts

## What does this do?

This script will take a FASTA file with a single DNA sequence (e.g., a
bacterial chromosome), and a VCF file containing the position of
SNPs along the FASTA file (e.g., as outputted from [snippy-core](https://www.github.com/tseemann/snippy))
and will output the ``<data>`` tag ready to copy paste into your XML file.

It will optionally also take a BED file with positions to mask (e.g., positions
of phage).

## How to install?

    pip3 install b2constsites

## How to run it?

If installed correctly, a scripted called `run_b2cs` should be available in
your path:

    run_b2cs --help

At a minimum, you need to supply a sequence file with your reference (assuming
it has a single chromosome entry --- this was designed for bacterial genomics,
but may work with viral too), and a VCF file with variants.

    run_b2cs myref.fasta myvar.vcf


The `data` tag will be printed to screen. The tag then should be cut and pasted
into your XML file just below the original `data` tab. You should take note
to modify the `data` tags as per the Introduction above.

## ASSUMPTIONS

This script will:

*  Only take in to account SNPs and MNPs annotated in the VCF.
Other variant types will be ignored.
*  Will only take into consideration A, C, G, and T bases in your reference
sequence. All other characters will be ignored.

The output will be, therefore, an approximation. However, it should be a close
enough approximation that it will provide a better inference from BEAST2 than
if once uses only variable sites, and then corrects in some *post hoc* manner.


## Issues or Questions

[GitHub Issues](https://github.com/andersgs/beast2_constsites/issues)
