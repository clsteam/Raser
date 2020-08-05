#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author  : Yao

import os
import argparse
import logging
import subprocess

from r_script import R_OF_VENN


def conversion_function_symbol(script: str) -> str:
    """Replace curly braces"""

    return script.replace("@<", "{").replace(">@", "}")


def main():
    parser = argparse.ArgumentParser()

    # required

    # optional
    parser.add_argument("i", metavar="input", type=str, nargs=2, help="multiple 'result.csv'(differential expression result) as input")
    parser.add_argument("-n", metavar="name", help="sample name", nargs=2, type=str)
    parser.add_argument("-o", help="output directory, default is ./venn", default="./venn")
    parser.add_argument("-p", "--p", help="padj, default is 0.01", default=0.01)
    parser.add_argument("-lfc", "--lfc", help="log2 fold-change threshold, default is 2", default=2)

    args = parser.parse_args()

    logger = logging.getLogger(__name__)

    params = dict(
        root_path=os.path.abspath(os.path.expanduser(args.o)),
        input_1=os.path.abspath(os.path.expanduser(args.i[0])),
        input_2=os.path.abspath(os.path.expanduser(args.i[1])),
        name1=args.n[0],
        name2=args.n[1],
        padj=args.p,
        lfc=args.lfc
    )

    script_path = os.path.join(params.get("root_path"), "venn.R")

    if not os.path.exists(params.get("root_path")):
        os.makedirs(params.get("root_path"))

    with open(script_path, "w") as handle:
        handle.write(conversion_function_symbol(R_OF_VENN.format(**params)))

    error_code = subprocess.call(" ".join(["Rscript", script_path]), shell=True)

    if error_code:
        raise Exception("R script cannot run correctly, please run the test yourself.")

    logger.info("[venn] success!")


if __name__ == '__main__':
    main()

