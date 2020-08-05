#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author  : Yao

""" Alignment

"""

# Standard library
import os

# Private module
from state2 import workflow_stateful
from core.abstract import AbstractWorkflow, AbstractTools, AbstractRun

from core.decorator import runshell
from source import tophat2, hisat2, star
from raser.settings import TOOLS_SELECTED, WHETHER_ALIGNMENT_WITH_ANNOTATIONS
from params.tools import ParamsStar


class Tophat2Run(AbstractRun):
    unpaired_prefix = ["", ]
    paired_prefix = ["", ""]

    genome_index_prefix = ""

    annotations_prefix = "-G"

    strand_specific_prefix = "--library-type"

    threads_prefix = "-p"

    def __init__(self, *args, **kwargs):
        self.again = kwargs.get("again", False)
        super(Tophat2Run, self).__init__(*args, **kwargs)

    @property
    def bp(self) -> bool:
        if self.again:
            return True
        if os.path.exists(self.dm.bam_file) or os.path.exists(self.dm.org_bam):
            return False
        return True

    def clear(self):
        pass

    @runshell
    def run(self):
        if os.path.exists(self.dm.org_bam):
            return
            # cmd = " ".join((tophat2,
            #                 "-R", os.path.dirname(self.dm.org_bam)
            #                 ))
            # return cmd

        # (BUG)tophat2 must change current dictionary
        os.chdir(self.dm.home_dir)
        # tophat2: Align the sequences to the reference genome (based on bowtie2(if dont have bowtie2, then use bowtie))
        cmd = " ".join((tophat2,
                        self.strand_specific,
                        self.fac_annotations,
                        self.threads,
                        self.genome_index_bowtie2,
                        self.input_reads,
                        # "-o", os.path.dirname(self.dm.org_bam),
                        ))
        return cmd

    @property
    def fac_annotations(self) -> str:
        if WHETHER_ALIGNMENT_WITH_ANNOTATIONS:
            return self.annotations
        return ""


class Hisat2Run(AbstractRun):
    phred_prefix = "--"

    unpaired_prefix = ["-U", ]
    paired_prefix = ["-1", "-2"]

    genome_index_prefix = "-x"

    annotations_prefix = "--known-splicesite-infile"
    strand_specific_prefix = "--rna-strandness"

    threads_prefix = "-p"

    def __init__(self, *args, **kwargs):
        self.again = kwargs.get("argin", False)
        super(Hisat2Run, self).__init__(*args, **kwargs)

    @property
    def bp(self) -> bool:
        if self.again:
            return True
        if os.path.exists(self.dm.bam_file) or os.path.exists(self.dm.org_bam):
            return False
        return True

    def clear(self):
        pass

    @runshell
    def run(self):
        cmd = " ".join((hisat2,
                        self.phred,
                        self.genome_index_hisat2,
                        self.input_reads,
                        self.strand_specific,
                        self.annotations,
                        self.threads,
                        "-S", self.dm.org_bam
                        ))
        return cmd

    @property
    def annotations_suffix(self) -> str:
        return self.ini.splicesites_txt

    def void_strand_specific(self, library_type):
        if not library_type:
            return
        if library_type in ("fr-firststrand",):
            if self.dm.sequence_format == "SE":
                return "R"
            else:
                return "RF"
        elif library_type in ("fr-secondstrand",):
            if self.dm.sequence_format == "SE":
                return "F"
            else:
                return "FR"
        elif library_type in ("fr-unstranded",):
            return
        else:
            self.logger.warning("[{0}] unkown library type: {1}".format(self.dm.id, library_type))
            return


class StarRun(AbstractRun):
    threads_prefix = "--runThreadN"
    genome_index_prefix = "--genomeDir"
    annotations_prefix = "--sjdbGTFfile"

    unpaired_prefix = ["--readFilesIn", ]
    paired_prefix = ["--readFilesIn", ]

    def __init__(self, *args, **kwargs):
        self.again = kwargs.get("argin", False)
        super(StarRun, self).__init__(*args, **kwargs)

    @property
    def bp(self) -> bool:
        if self.again:
            return True
        if os.path.exists(self.dm.bam_file) or os.path.exists(self.dm.org_bam):
            return False
        return True

    def clear(self):
        pass

    @property
    def add_params(self):
        if TOOLS_SELECTED.get("transcript") == "cufflinks":
            if not self.ini.library_type and self.ini.library_type == "fr-unstranded":
                return "--outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical"
            else:
                return "--outSAMstrandField intronMotif"
        return ""

    @runshell
    def run(self):
        cmd = " ".join((star,
                        self.threads,
                        self.genome_index_star,
                        self.annotations_gtf,
                        "--readFilesCommand", "zcat",
                        self.input_reads,
                        "--outFileNamePrefix", os.path.dirname(self.dm.org_bam),
                        "--outSAMtype BAM SortedByCoordinate",
                        "--quantMode GeneCounts",
                        self.add_params,
                        ParamsStar.common_additional_param,
                        ParamsStar.align_additional_params
                        ))
        return cmd


@workflow_stateful
class Alignment(AbstractWorkflow):
    """Sequence Alignment"""

    class Tophat2(AbstractTools):
        class Run(Tophat2Run):
            pass

    class Hisat2(AbstractTools):
        class Run(Hisat2Run):
            pass

    class Star(AbstractTools):
        class Run(StarRun):
            pass
