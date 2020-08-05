#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author  : Yao

""" Title

This file(script) can also be imported as a module and contains the following
functions:

    * main - the main function of the script
    * function - returns the column headers of the file
"""

# Standard library
import os

from state2 import workflow_stateful
from core.abstract import AbstractWorkflow, AbstractTools, AbstractRun

from core.decorator import runshell
from source import featurecounts

from params.tools import ParamsGeneCount


def priority_first(a, b):
    return a if a else b


class FeaturecountsRun(AbstractRun):
    annotations_prefix = "-a"
    strand_specific_prefix = "-s"

    threads_prefix = "-T"

    def __init__(self, *args, **kwargs):
        """
        :param args:
        :param kwargs:
            'out':
            'feature':
            'attribute':
        """
        self.out = priority_first(kwargs.get("out", None), args[1].gene_counts_file)
        self.feature = priority_first(kwargs.get("feature", None), ParamsGeneCount.feature_type)
        self.attribute = priority_first(kwargs.get("attribute", None), ParamsGeneCount.attribute_type)
        super(FeaturecountsRun, self).__init__(*args, **kwargs)

    @property
    def bp(self) -> bool:
        # 对照组改变时，lncrna的差异分析的每个样本lncrna的counts需要重新计算
        if self.out != self.dm.gene_counts_file:
            return True
        if os.path.exists(self.out):
            return False
        return True

    def clear(self):
        pass

    @property
    def se_or_pe(self):
        if self.dm.sequence_format in ("SE",):
            return ""
        return "-p"

    @runshell
    def run(self):
        """
        software: Subread
        :return: cmd
        """
        cmd = " ".join((featurecounts,
                        self.se_or_pe,
                        self.strand_specific,
                        self.threads,
                        self.annotations_gtf,
                        self.dm.bam_file,
                        "-t", self.feature,
                        "-g", self.attribute,
                        "-o", self.out
                        ))
        return cmd

    def void_strand_specific(self, library_type):
        if not library_type:
            return
        if library_type in ("fr-firststrand",):
            return "1"
        elif library_type in ("fr-secondstrand",):
            return "2"
        elif library_type in ("fr-unstranded",):
            return "0"
        else:
            self.logger.warning("[{0}] unkown library type: {1}".format(self.dm.id, library_type))
            return


class HtseqRun(AbstractRun):
    strand_specific_prefix = "-s"

    annotations_prefix = ""

    def __init__(self, *args, **kwargs):
        self.out = priority_first(kwargs.get("out", None), self.dm.gene_counts_file)
        self.feature = priority_first(kwargs.get("feature", None), ParamsGeneCount.feature_type)
        self.attribute = priority_first(kwargs.get("attribute", None), ParamsGeneCount.attribute_type)
        super(HtseqRun, self).__init__(*args, **kwargs)

    @property
    def bp(self) -> bool:
        if os.path.exists(self.out):
            return False
        return True

    def clear(self):
        pass

    @runshell
    def run(self):
        """
        software: Htseq
        :return: cmd
        """
        # if platform.python_version_tuple() > ('3', '7', '0'):
        #     raise ExtensionException(
        #         "Htseq does not support this python version, please see Htseq official documentation")
        # return
        cmd = " ".join(("python3 -m HTSeq.scripts.count",
                        "-f", "bam",
                        self.strand_specific,
                        self.dm.bam_file,
                        self.annotations_gtf,
                        "-t", self.feature,
                        "-i", self.attribute,
                        ">", self.out
                        ))
        return cmd

    def void_strand_specific(self, library_type):
        if not library_type:
            return
        if library_type in ("fr-firststrand",):
            return "yes"
        elif library_type in ("fr-secondstrand",):
            return "reverse"
        elif library_type in ("fr-unstranded",):
            return "no"
        else:
            self.logger.warning("[{0}] unkown library type: {1}".format(self.dm.id, library_type))
            return


@workflow_stateful
class GeneCount(AbstractWorkflow):
    """Gene raw counts"""
    class Featurecounts(AbstractTools):
        class Run(FeaturecountsRun):
            pass

    class Htseq(AbstractTools):
        class Run(HtseqRun):
            pass
