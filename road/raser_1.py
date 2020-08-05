#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author  : Yao

""" Title

This file(script) can also be imported as a module and contains the following
functions:

    * main - the main function of the script
    * function - returns the column headers of the file
"""


import os
import glob

from state2 import workflow_stateful
from core.abstract import AbstractWorkflow, AbstractTools, AbstractRun
from core.decorator import runshell
from source import fastqc


class FastqcRun(AbstractRun):

    threads_prefix = "-t"

    @property
    def bp(self) -> bool:
        if glob.glob(os.path.join(self.dm.qc_dir, "*clean_fastqc/fastqc_data.txt")):
            return False
        return True

    def clear(self):
        pass

    @runshell
    def run(self):
        if not glob.glob(os.path.join(self.dm.qc_dir, "*/fastqc_data.txt")):
            return self._fastqc(self.dm.sample, "")
        if glob.glob(os.path.join(self.dm.home_dir, "*_clean.fq.gz")):
            return self._fastqc(self.dm.clean_data, "_clean")
        return

    def _fastqc(self, sample_tuple, has_clean) -> tuple:
        """
        :param sample_tuple: orignal fastq file
        :param has_clean: "" or "_clean"
        :return: commands of fastqc
        """
        if self.dm.sequence_format == "SE":
            cmd_1 = " ".join((fastqc,
                              sample_tuple[0],
                              self.threads,
                              "-O", self.dm.qc_dir
                              ))
            qc_zip = "".join((os.path.join(self.dm.qc_dir, self.dm.id), has_clean, "_fastqc.zip"))
            cmd_2 = " ".join(("unzip",
                              qc_zip,
                              "-d", self.dm.qc_dir
                              ))
            return cmd_1, cmd_2
        else:
            cmd_11 = " ".join((fastqc,
                               sample_tuple[0],
                               self.threads,
                               "-O", self.dm.qc_dir
                               ))
            cmd_12 = " ".join((fastqc,
                               sample_tuple[1],
                               self.threads,
                               "-O", self.dm.qc_dir
                               ))
            qc_zip_1 = "".join((os.path.join(self.dm.qc_dir, self.dm.id), "_1", has_clean, "_fastqc.zip"))
            qc_zip_2 = "".join((os.path.join(self.dm.qc_dir, self.dm.id), "_2", has_clean, "_fastqc.zip"))
            cmd_21 = " ".join(("unzip",
                               qc_zip_1,
                               "-d", self.dm.qc_dir
                               ))
            cmd_22 = " ".join(("unzip",
                               qc_zip_2,
                               "-d", self.dm.qc_dir
                               ))
            return cmd_11, cmd_12, cmd_21, cmd_22


@workflow_stateful
class QualityControl(AbstractWorkflow):
    """quality detection"""

    class Fastqc(AbstractTools):

        class Run(FastqcRun):
            pass

