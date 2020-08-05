#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author  : Yao

"""
1> sort, remove duplication(PCR) and build index for BAM
2> add @RG

"""

# Standard library
import os

# Private module
from state2 import workflow_stateful
from core.abstract import AbstractWorkflow, AbstractTools, AbstractRun
from core.decorator import runshell
from raser.settings import TOOLS_SELECTED, WHETHER_MARK_DUPLICATES_ONLY, RGPL
from source import picard_jar


class CommonRun(AbstractRun):

    @property
    def bp(self) -> bool:
        return True

    def clear(self):
        pass

    def run(self):
        pass

    @runshell
    def after(self):
        tmp_bam = os.path.join(self.dm.home_dir, "tmp_rg.bam")
        cmd_1 = " ".join(("java",
                          "-jar",
                          picard_jar,
                          "AddOrReplaceReadGroups",
                          "I={0}".format(self.dm.bam_file),
                          "O={0}".format(tmp_bam),
                          "ID=RG{0}".format(self.dm.id),
                          "LB=LB{0}".format(self.dm.id),
                          "PL={0}".format(RGPL),
                          "PU=PU{0}".format(self.dm.id),
                          "SM={0}".format(self.dm.id)
                          ))
        cmd_2 = " ".join(("rm",
                          self.dm.bam_file
                          ))
        cmd_3 = " ".join(("mv",
                          tmp_bam,
                          self.dm.bam_file
                          ))
        cmd_4 = " ".join(("samtools",
                          "index",
                          self.dm.bam_file
                          ))
        return cmd_1, cmd_2, cmd_3, cmd_4


class SamtoolsRun(CommonRun):

    def __init__(self, *args, **kwargs):
        self.tmp_bam = ""
        super(SamtoolsRun, self).__init__(*args, **kwargs)

    @property
    def bp(self) -> bool:
        if os.path.exists(self.dm.bam_file):
            return False
        return True

    def clear(self):
        if os.path.exists(self.dm.org_bam):
            os.remove(self.dm.org_bam)
        if os.path.exists(self.tmp_bam):
            os.remove(self.tmp_bam)

    @runshell
    def before(self):
        """sort BAM|SAM -> sorted.bam"""
        if TOOLS_SELECTED.get("alignment") == "star":
            self.tmp_bam = os.path.join(self.dm.home_dir, "tmp.bam")
            return " ".join(("samtools",
                             "view",
                             "-Su",
                             self.dm.org_bam,
                             "|",
                             "samtools",
                             "sort",
                             "-o",
                             self.tmp_bam
                             ))
        return

    @runshell
    def run(self):
        """remove duplication(PCR)"""
        cmd = " ".join(("samtools",
                        "rmdup",
                        "-s" if self.dm.sequence_format == "SE" else "-S",
                        self.tmp_bam if self.tmp_bam else self.dm.org_bam,
                        self.dm.bam_file
                        ))
        return cmd


class PicardRun(CommonRun):

    @property
    def bp(self) -> bool:
        if os.path.exists(self.dm.bam_file):
            return False
        if os.path.exists(self.dm.org_bam):
            return True
        return False

    def clear(self):
        if os.path.exists(self.dm.org_bam):
            os.remove(self.dm.org_bam)

    @runshell
    def before(self):
        """sort BAM|SAM -> sorted.bam"""
        if TOOLS_SELECTED.get("alignment") == "star":
            self.tmp_bam = os.path.join(self.dm.home_dir, "tmp.bam")
            return " ".join(("samtools",
                             "view",
                             "-Su",
                             self.dm.org_bam,
                             "|",
                             "samtools",
                             "sort",
                             "-o",
                             self.tmp_bam
                             ))
        return

    @runshell
    def run(self):
        cmd = " ".join(("java",
                        "-jar",
                        picard_jar,
                        "MarkDuplicates",
                        "REMOVE_DUPLICATES=true" if WHETHER_MARK_DUPLICATES_ONLY else "",
                        "I={0}".format(self.dm.org_bam),
                        "O={0}".format(self.dm.bam_file),
                        "M={0}".format(os.path.join(self.dm.home_dir, "markdup_metrics.txt"))
                        ))
        return cmd


@workflow_stateful
class Rmdup(AbstractWorkflow):
    """Sequence Alignment"""

    class Samtools(AbstractTools):
        class Run(SamtoolsRun):
            pass

    class Picard(AbstractTools):
        class Run(PicardRun):
            pass
