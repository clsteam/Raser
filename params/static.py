#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from raser.settings import TOOLS_SELECTED

gene_counts = {
    "column": "2",
    "suffix": ".count"
} if TOOLS_SELECTED.get("genecount") == "htseq" else {
    "column": "7",
    "suffix": ".counts"
}

lncrna_counts = {
    "column": "2",
    "suffix": ".lnc.count"
} if TOOLS_SELECTED.get("genecount") == "htseq" else {
    "column": "7",
    "suffix": ".lnc.counts"
}

ballgown_data_dir = "extdata"
