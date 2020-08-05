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
import re
import subprocess
# Private module
from state2 import workflow_stateful
from core.abstract import AbstractWorkflow, AbstractTools, AbstractRun
from raser.settings import STRAND_SPECIFIC_USE_AUTOMATICALLY


class RseqcRun(AbstractRun):

    bed_prefix = "-r"

    @property
    def bp(self) -> bool:
        return True

    def clear(self):
        pass

    def run(self):
        library_type = self._rseqc
        self.logger.info("[0] library type: {1}".format(self.dm.id, library_type))
        if STRAND_SPECIFIC_USE_AUTOMATICALLY:
            self.dm.library_type = library_type

    @property
    def _rseqc(self) -> str:
        cmd = " ".join((
            "infer_experiment.py",
            self.bed,
            "-i", self.dm.bam_file
        ))
        data = subprocess.getoutput(cmd)
        a = float(re.findall(re.compile(".*\s(.*)"), data[-2])[0])
        b = float(re.findall(re.compile(".*\s(.*)"), data[-1])[0])
        if (a - b) > 0.5:
            library_type = "fr-secondstrand"
        elif (b - a) > 0.5:
            library_type = "fr-firststrand"
        else:
            library_type = "fr-unstranded"
        return library_type


@workflow_stateful
class StrandSpecific(AbstractWorkflow):
    """Strand Specific"""

    class Rseqc(AbstractTools):

        class Run(RseqcRun):
            pass
