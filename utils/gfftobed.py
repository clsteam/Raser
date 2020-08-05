#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author  : Yao

""" Title

This file(script) can also be imported as a module and contains the following
functions:

    * main - the main function of the script
    * function - returns the column headers of the file
"""

# Standard library
import argparse
import os
import re

# Third party library

# Private module


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("i", help="GTF/GFF file", nargs=1, type=str)
    parser.add_argument("-o", help="output BED file", required=True)
    parser.add_argument("-a", help="gene attribute in GFF/GTF file, which will be writen to BED file, default is gene_id[GTF]/ID[GFF]", default=None)

    args = parser.parse_args()

    gff2bed(args.i[0], args.o, args.a)


def gff2bed(gff, bed, attribute=None):
    fomat = os.path.splitext(gff)[1]
    if "gtf" in fomat:
        attb = attribute if attribute else "gene_id"
        pattern = '{0} "(.+?)";'.format(attb)
    else:
        attb = attribute if attribute else "ID"
        pattern = '{0}=(.+?)(;|\\n)'.format(attb)

    with open(gff, "r") as handle_r, open(bed, "w") as handle_w:
        for line in handle_r:
            if line.startswith("#"):
                continue
            # NW_011499845.1  Gnomon  gene    1474    2247    .       +       .       ID=gene-LOC105127819;Dbxref=GeneID:105127819;Name=LOC105127819;gbkey=Gene;gene=LOC105127819;gene_biotype=protein_coding
            ls = line.split()
            if "gene" in ls:
                res = re.search(pattern, ls[-1])
                if res:
                    gene_id = res.groups()[0]
                    handle_w.write("{0}\t{1}\t{2}\t{3}\n".format(ls[0], ls[3], ls[4], gene_id))


if __name__ == '__main__':
    main()
