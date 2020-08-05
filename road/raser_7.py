#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author  : Yao
import os
import glob
from functools import reduce

from state2 import workflow_stateful
from core.abstract import AbstractWorkflow, AbstractTools, AbstractRun
from core.decorator import runshell
from core.r_rule import conversion_function_symbol, list_to_vector
from raser.settings import TOOLS_SELECTED
from params.static import gene_counts, ballgown_data_dir

from .r_script import R_OF_DESEQ, R_OF_BALLGOWN

def creat_dir(*folder):
    for fd in folder:
        if not os.path.exists(fd):
            os.makedirs(fd)
    return folder[0]


def priority_first(a, b):
    return a if a else b


class Deseq2Run(AbstractRun):

    def __init__(self, *args, **kwargs):
        """
        :param args:
            logger: 多进程主进程日志句柄
            dm: 文件管理实例
            ini: 配置
        :param kwargs:
            'counts_dict':
            'out':
        """
        self.sp_counts = None
        self.treatment_vector = None
        self.treatment_dict = None

        self.out = priority_first(kwargs.get("out", None), args[2].diff_dir)
        self.script = os.path.join(self.out, "deseq2.R")
        self.matrix = os.path.join(self.out, "gene.raw_counts")
        self.counts_dict = kwargs.get("counts_dict", gene_counts)
        super(Deseq2Run, self).__init__(*args, **kwargs)

    def bp(self) -> bool:
        if os.path.exists(self.script):
            return False
        return True

    def clear(self):
        pass

    def run(self):
        if not self.ini.has_comparison:
            self.logger.critical("differential analysis requires 2 group data.")
            return
        self._generate_counts_matrix()
        self._touch_rscript()

    @runshell
    def _generate_counts_matrix(self):
        creat_dir(self.out)
        tmpfile = [os.path.join(self.out, x) for x in ("tmp1", "tmp2")]
        sp1_pat = [os.path.join(x, "*", "*" + self.counts_dict.get("suffix")) for x in self.ini.experiment_group]
        sp2_pat = [os.path.join(x, "*", "*" + self.counts_dict.get("suffix")) for x in self.ini.control_group]
        sp1_file = sorted(reduce(lambda a, b: a + b, [glob.glob(x) for x in sp1_pat]))
        sp2_file = sorted(reduce(lambda a, b: a + b, [glob.glob(x) for x in sp2_pat]))

        # order : Control,Experiment
        list_file = sp1_file + sp2_file
        self.sp_counts = (len(sp1_file), len(sp2_file))
        sp_id_list = [os.path.basename(os.path.dirname(x)) for x in list_file]
        self.treatment_dict = {x: self.ini.treatment.get(x, "None") for x in sp_id_list}
        self.treatment_vector = list_to_vector(sorted(self.treatment_dict.values()))

        shell_merge_cmd_prefix = "awk '{_[FNR]=(_[FNR] OFS $counts_col)}END{for (i=1; i<=FNR; i++) {sub(/^ /,ANCHOR,_[i]); print _[i]}}'".replace('ANCHOR', '""').replace("counts_col", self.counts_dict.get("column"))
        cmd1 = " ".join((shell_merge_cmd_prefix, " ".join(list_file), "|sed '1,2c {0}'".format(" ".join(sp_id_list)) if TOOLS_SELECTED.get("genecount") == "featurecounts" else "|sed '1i {0}'".format(" ".join(sp_id_list)), ">", tmpfile[0]))

        cmd2 = " ".join(("awk '{print $1}'", sp1_file[0], "|sed '1,2d' |sed '1i GeneName' >", tmpfile[1]))
        cmd3 = " ".join(("paste -d' '", tmpfile[1], tmpfile[0], ">", self.matrix))

        cmd4 = " ".join(["rm"] + tmpfile)
        return cmd1, cmd2, cmd3, cmd4

    @runshell
    def _touch_rscript(self):
        params_dict = {
            "root_path": self.out,
            "matrix": os.path.basename(self.matrix),
            "threads": self.ini.ppn,
            "E_name": self.ini.experiment_name,
            "C_name": self.ini.control_name,
            "E_num": self.sp_counts[0],
            "C_num": self.sp_counts[1],
            "treatment_vector": self.treatment_vector,
            "treatment_levels": ", ".join(["".join(["'", x, "'"]) for x in set(self.treatment_dict.values())])
        }
        with open(self.script, "w") as doc:
            doc.write(conversion_function_symbol(R_OF_DESEQ.format(**params_dict)))
        return " ".join(("Rscript", self.script)), 1


class BallgownRun(AbstractRun):

    def __init__(self, *args, **kwargs):
        """
        :param args:
        :param kwargs:
            'out':
        """
        self.sp1 = None
        self.sp2 = None
        super(BallgownRun, self).__init__(*args, **kwargs)

    @property
    def bp(self) -> bool:
        if not self.out:
            self.out = self.ini.diff_dir
        creat_dir(self.out)
        self.script = os.path.join(self.out, "ballgown.R")
        if os.path.exists(self.script):
            return False
        return True
        
    def clear(self):
        pass

    def run(self):
        self.generate_sample_order_dict()
        self._touch_rscript()

    def generate_sample_order_dict(self):
        sp1_pat = reduce(lambda a, b: a + b, [glob.glob(os.path.join(x, "*", "*.bam")) for x in self.ini.experiment_group])
        sp2_pat = reduce(lambda a, b: a + b, [glob.glob(os.path.join(x, "*", "*.bam")) for x in self.ini.control_group])
        self.sp1 = [os.path.basename(os.path.dirname(x)) for x in sp1_pat]
        self.sp2 = [os.path.basename(os.path.dirname(x)) for x in sp2_pat]

    @runshell
    def _touch_rscript(self):
        params_dict = {
            "root_path": self.out,
            "extdata": ballgown_data_dir,
            "sample_vector": list_to_vector(self.sp1 + self.sp2),
            "E_name": self.ini.experiment_name,
            "C_name": self.ini.control_name,
            "E_num": len(self.sp1),
            "C_num": len(self.sp2),
        }
        with open(self.script, "w") as doc:
            doc.write(conversion_function_symbol(R_OF_BALLGOWN.format(**params_dict)))
        return " ".join(("Rscript", self.script)), 1


@workflow_stateful
class DifferentialExpression(AbstractWorkflow):
    """differential expression analysis for RNA-Seq"""

    class Deseq2(AbstractTools):

        class Run(Deseq2Run):
            pass

    class Ballgown(AbstractTools):

        class Run(BallgownRun):
            pass
