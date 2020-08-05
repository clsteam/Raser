#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" 软件参数限制
例如：
    1> 分析物种对象限制：tophat-fusion，CPAT

"""


class SpeciesLimit(object):
    """
    key: name of allowed species, which user needs to write in configure file
    value: annotations
    """

    tophatfusion = {
        "homo": ("hg18", "hg19", "hg38"),
        "mus musculus": ("mm9",),
        "rattus norvegicus": ("rn4",),
        "canis familiaris": ("canFam2",),
        "gallus gallus": ("galGal3",),
        "oryctolagus cuniculus": ("oryCun2",)
    }
    cpat = {
        "homo": "Human",
        "fly": "Fly",
        "zebrafish": "zebrafish",
        # mouse
        "mus musculus": "Mouse",
        "rattus norvegicus": "Mouse",
        "mouse": "Mouse",
    }
