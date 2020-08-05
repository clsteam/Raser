#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import glob

from state2 import workflow_stateful
from core.abstract import AbstractWorkflow, AbstractTools, AbstractRun
from core.decorator import runshell
from params.tools import ParamsTophatFusion
from source import tophat2, tophat_fusion_post, star_fusion
from core.ilog import creat_dir


class Tophatfusion1Run(AbstractRun):
    """
    parallel
    """
    threads_prefix = "-p"

    genome_index_prefix = ""

    unpaired_prefix = ["", ]
    paired_prefix = ["", ""]

    @property
    def bp(self) -> bool:
        if self.ori_or_done_res:
            return False
        return True

    def clear(self):
        pass

    @property
    def ori_or_done_res(self):
        return glob.glob(os.path.join(self.ini.fusion_dir, "*", "tophat_" + self.dm.id, "accepted_hits.bam"))

    @property
    def origin_res(self):
        return os.path.join(self.ini.fusion_bam, "tophat_" + self.dm.id)

    @runshell
    def run(self):
        """
        software: tophat-fusion
        :return: cmd
        """
        os.chdir(self.ini.fusion_bam)
        cmd = " ".join((tophat2,
                        self.threads,
                        "-o", self.origin_res,
                        "--fusion-search --bowtie1",
                        ParamsTophatFusion.additional_params,
                        self.genome_index_bowtie1,
                        self.input_reads
                        ))
        return cmd


class Tophatfusion2Run(AbstractRun):
    """
    single
    """

    @property
    def bp(self) -> bool:
        if not os.listdir(self.ini.fusion_bam):
            return False
        return True

    def clear(self):
        pass

    def run(self):
        for name, sample in self.ini.treatment_accession:
            absolute_sample = [os.path.join(self.ini.fusion_bam, "tophat_" + x) for x in sample]
            test = map(os.path.exists, absolute_sample)
            if all(test):
                result_dir = os.path.join(self.ini.fusion_dir, name)
                creat_dir(result_dir)
                for per_sample in absolute_sample:
                    shutil.move(per_sample, result_dir)
            if not any(map(os.path.exists, absolute_sample)):
                self.logger.debug("[FUSION BAM] {0}(treatment) had been moved".format(name))
                continue
            self.logger.warning("[FUSION BAM] {0}(treatment) have extra samples".format(name))


class Tophatfusion3Run(AbstractRun):
    """
    parallel
    """
    threads_prefix = "-p"

    genome_index_prefix = ""

    def __init__(self, *args, **kwargs):
        self.name_accession = {x: y[0] for x, y in args[2].treatment_accession}
        self.accession_name = {y: x for x, y in self.name_accession.items()}
        self.result_dir = None
        super(Tophatfusion3Run, self).__init__(*args, **kwargs)

    @property
    def bp(self) -> bool:
        if self.dm.id not in self.name_accession.values():
            return False

        self.result_dir = os.path.join(self.ini.fusion_dir, self.accession_name.get(self.dm.id))
        if self.is_file_empty(os.path.join(self.result_dir, "tophatfusion_out", "result.txt"), 1):
            shutil.rmtree(os.path.join(self.result_dir, "tophatfusion_out"))
            return True
        return False

    def clear(self):
        pass

    @runshell
    def run(self):
        """
        software: tophat-fusion
        :return: cmd
        """
        os.chdir(self.result_dir)
        cmd = " ".join((tophat_fusion_post,
                        self.threads,
                        self.genome_index_bowtie1,
                        "" if self.ini.species == "homo" else "--non-human",

                        ))
        return cmd


class StarfusionRun(AbstractRun):
    unpaired_prefix = ["--left_fq", ]
    paired_prefix = ["--left_fq", "--right_fq"]

    @property
    def bp(self) -> bool:
        pass

    def clear(self):
        pass

    @runshell
    def run(self):
        """
        software: STAR-Fusion
        :return: cmd
        """
        cmd = " ".join((star_fusion,
                        "--genome_lib_dir", self.ini.starfusion_genome_resource_lib,
                        self.input_reads,
                        "--output_dir", os.path.join(self.ini.fusion_dir, self.dm.id)
                        ))
        return cmd


@workflow_stateful
class Fusion(AbstractWorkflow):
    """fusion gene"""

    class Tophatfusion1(AbstractTools):
        class Run(Tophatfusion1Run):
            pass

    class Tophatfusion2(AbstractTools):
        class Run(Tophatfusion2Run):
            pass

    class Tophatfusion3(AbstractTools):
        class Run(Tophatfusion3Run):
            pass

    class Starfusion(AbstractTools):
        class Run(StarfusionRun):
            pass
