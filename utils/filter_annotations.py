#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
trans-spliced isn't handled well by star-fusion, this script can help you remove trans-spliced in genome annotations file.
"""

import re
import os
import argparse
import subprocess


def main():
    # arguments passed
    parser = argparse.ArgumentParser()

    # required
    parser.add_argument("file", help="input GFF annotations file")
    parser.add_argument("--rn_ts", help="rename all trans-splicing", default=True, action="store_true")

    args = parser.parse_args()

    input_gff = os.path.abspath(os.path.expanduser(args.file))
    tmp_gff = os.path.splitext(input_gff)[0] + ".filtered.gff"
    out_gtf = os.path.splitext(input_gff)[0] + ".filtered.gtf"
    with open(args.file, "r") as handle_r, open(tmp_gff, "w") as handle_w:
        for line in handle_r:
            if re.match(".+exception=trans-splicing;.*", line):
                continue
            handle_w.write(line)
    subprocess.call("gffread {0} -T -o {1}".format(tmp_gff, out_gtf), shell=True)
    os.remove(tmp_gff)


if __name__ == '__main__':
    main()
