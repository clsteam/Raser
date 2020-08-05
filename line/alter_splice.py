#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author  : Yao

""" Alter Splice
"""

# Standard library
import os

import pandas as pd

from state2 import workflow_stateful
from core.abstract import AbstractWorkflow, AbstractTools, AbstractRun

from core.decorator import runshell
from source import extract_as, summarize_as, ea_fpkm, gffcompare


class Asprofile1Run(AbstractRun):
    annotations_prefix = "-r"

    threads_prefix = "-T"

    genome_prefix = "-s"
    
    def __init__(self, *args, **kwargs):
        self.org_fpkm = ""
        self.fpkm = ""
        self.as_min = ""
        super(Asprofile1Run, self).__init__(*args, **kwargs)

    @property
    def bp(self) -> bool:
        self.org_fpkm = os.path.join(self.dm.as_dir, self.dm.id + ".xfpkm")
        self.fpkm = os.path.join(self.dm.as_dir, self.dm.id + ".fpkm")
        self.as_min = os.path.join(self.dm.as_dir, self.dm.id + ".mini")
        if not self.is_file_empty(self.fpkm):
            return False
        return True

    def clear(self):
        pass

    @runshell
    def run(self):
        """
        software: ASprofile
        :return: cmd
        """
        cmd_1 = " ".join((gffcompare,
                          self.dm.gtf,
                          self.annotations_gtf,
                          self.genome,
                          "-C",
                          "-o", os.path.join(self.dm.as_dir, "prefix_1")
                          ))
        cmd_2 = " ".join((extract_as,
                          self.dm.gtf,
                          self.ini.hdrs,
                          "-r", os.path.join(self.dm.transcript_dir, "prefix_1.transcripts.gtf.tmap"),
                          self.annotations_gtf.split()[1],
                          ">", os.path.join(self.dm.as_dir, "total.as")
                          ))
        cmd_3 = " ".join((summarize_as,
                          self.dm.gtf,
                          os.path.join(self.dm.as_dir, "total.as"),
                          "-p", os.path.join(self.dm.as_dir, "prefix_2")
                          ))
        cmd_4 = " ".join((ea_fpkm,
                          self.dm.gtf,
                          self.ini.hdrs,
                          os.path.join(self.dm.as_dir, "prefix_2.as.nr"),
                          ">", self.org_fpkm
                          ))
        return cmd_1, cmd_2, cmd_3, cmd_4

    def after(self):
        # remove those variations that allow for some 'wiggle' room at the boundaries of the surrounding features, which prefixed with 'X'
        with open(self.org_fpkm, "r") as handle_r, open(self.fpkm, "w") as handle_w:
            handle_r.readline()
            for line in handle_r:
                if line.split()[1].startswith("X"):
                    continue
                handle_w.write(line)

        # count
        ct = pd.read_table(self.fpkm)
        with open(self.as_min, "w") as handle:
            i = 0
            pre_gene = ct["gene_id"][i]
            start = i

            i += 1
            while i < len(ct):
                if pre_gene != ct["gene_id"][i]:
                    handle.write("{0}\t{1}\n".format(pre_gene, sum(ct["fpkm"][start:i])))

                    pre_gene = ct["gene_id"][i]
                    start = i
                i += 1
            handle.write("{0}\t{1}\n".format(pre_gene, sum(ct["fpkm"][start:i])))
        return


class Asprofile2Run(AbstractRun):

    @property
    def bp(self) -> bool:
        return True

    def clear(self):
        pass

    @runshell
    def run(self):
        pass


@workflow_stateful
class AlterSplice(AbstractWorkflow):
    """alter splice"""

    class Asprofile1(AbstractTools):
        class Run(Asprofile1Run):
            pass

    class Asprofile2(AbstractTools):
        class Run(Asprofile2Run):
            pass
