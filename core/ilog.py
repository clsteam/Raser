#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" read config file

This file(script) can also be imported as a module and contains the following
functions:

    * creat_dir - function of creat folder
    * Ini - class of return config file handle
"""

import os
import glob
from logging import Logger
from functools import reduce
from configparser import ConfigParser, ExtendedInterpolation

from core.exception import ConfigException
from raser.settings import TOOLS_SELECTED


def creat_dir(*folder):
    for fd in folder:
        if not os.path.exists(fd):
            os.makedirs(fd)
    return folder[0]


def replace_blank_space(s: str) -> str:
    return "_".join(s.split())


class Ini(object):

    category = None  # "ve", "pl"
    treatment = {}

    workflow_dict = {}

    # tmp
    _group_accession = []
    _treatment_accession = []
    _accession_home = {}

    def __init__(self, ini_file, logger=Logger, category="ve"):
        """
        :param ini_file: config.ini
        :param logger: required
        :param category: CNCI params, -m or --model : assign the classification models ("ve" for vertebrate species, "pl" for plat species)
        """
        self.ini_file = ini_file
        self.main_logger = logger
        self.category = category
        self.handle = ConfigParser(interpolation=ExtendedInterpolation())
        self.handle.read_file(open(self.ini_file))
        creat_dir(self.log_path)

    def initialize(self):
        if int(self.ppn) % int(self.pools):
            raise ConfigException("The remainder of 'ppn' and 'pools' must be 0")
        if len(self.group) != 2:
            raise ConfigException("You must enter a name for both sets of data (data can be empty)")
        creat_dir(self.pipe_log_dir, self.log_path, self.fusion_dir, self.fusion_bam, self.ase_result_dir)
        self.check()

    def check(self):
        """treatment
        key: accession
        value: treatment
        """
        treatment_header = self.handle.get("Treatment", "header_name")
        treatment_file = self.handle.get("Treatment", "file")
        if not treatment_header or not treatment_file:
            self.main_logger.warning("treatment file or header not found and each sample entered will be treated as a separate individual")
            self.treatment = {}
            return
        if treatment_file.endswith("csv"):
            sep = ","
        else:
            sep = None
        with open(treatment_file, "r") as handle:
            header = handle.readline().strip().split(sep)
            sample_index, treat_index = map(lambda x: header.index(x), treatment_header.split(","))
            for line in handle:
                ls = line.strip().split(sep)
                self.treatment[ls[sample_index]] = replace_blank_space(ls[treat_index])
        if not self.treatment:
            self.main_logger.error("{0} is empty, and will not use treatment".format(treatment_file))

    def set(self, section, option, value=None):
        self.handle.set(section, option, value)

    def save(self):
        with open(self.ini_file, 'w') as cf:
            self.handle.write(cf)

    @property
    def group(self):
        """
        EXPERIMENT_GROUP_NAME, EXPERIMENT_GROUP
        CONTROL_GROUP_NAME, CONTROL_GROUP
        """
        return self.handle.items("SampleDir")

    @property
    def _sp1_bam(self) -> list:
        if not self.experiment_group:
            return []
        return reduce(lambda a, b: a + b, [glob.glob(os.path.join(x, "*", "*.bam")) for x in self.experiment_group])

    @property
    def _sp2_bam(self) -> list:
        if not self.control_group:
            return []
        return reduce(lambda a, b: a + b, [glob.glob(os.path.join(x, "*", "*.bam")) for x in self.control_group])

    @property
    def accession_home(self) -> dict:
        """
        key: accession (str)
        value: absolute path of sample home dir
        """
        if self._accession_home:
            return self._accession_home
        for x in self._sp1_bam+self._sp2_bam:
            self._accession_home[os.path.basename(os.path.dirname(x))] = os.path.dirname(x)
        return self._accession_home

    @property
    def group_accession(self) -> list:
        """
        [(group name, [accession, ...]), (...)]
        eg:
            [("E", ["", "", ...]), (...)]
        """
        if self._group_accession:
            return self._group_accession
        sp1 = [os.path.basename(os.path.dirname(x)) for x in self._sp1_bam]
        sp2 = [os.path.basename(os.path.dirname(x)) for x in self._sp2_bam]
        self._group_accession = [(self.experiment_name, sp1),
                                 (self.control_name, sp2)
                                 ]
        return self._group_accession

    @property
    def treatment_accession(self) -> list:
        """ 2 layer
        if no treatment input, return:
            [(group name, [accession, ...]), (...)]
        else return:
            [(group name + treatment name, [accession, ...]), (...)]
        """
        if self._treatment_accession:
            return self._treatment_accession
        if not self.treatment:
            return self.group_accession
        from collections import defaultdict
        res = defaultdict(list)
        for group, sp in self.group_accession:
            for accession in sp:
                group_treatment = "_".join([group, self.treatment.get(accession)])
                res[group_treatment].append(accession)
        self._treatment_accession = list(res.items())
        return self._treatment_accession

    @property
    def experiment_group(self) -> list:
        """
        EXPERIMENT_GROUP
        """
        return self.group[0][1].strip().split()

    @property
    def experiment_name(self) -> str:
        """
        EXPERIMENT_GROUP_NAME
        """
        return replace_blank_space(self.group[0][0])

    @property
    def has_comparison(self) -> bool:
        return all([self.experiment_group, self.control_group])

    @property
    def control_group(self) -> list:
        """
        CONTROL_GROUP
        """
        return self.group[1][1].strip().split()

    @property
    def control_name(self) -> str:
        """
        CONTROL_GROUP_NAME
        """
        return replace_blank_space(self.group[1][0])

    @property
    def workflow(self) -> dict:
        """
        :return:{
            "differentialexpression": "deseq2",
            "altersplice": False,
            "fusion": "tophatfusion",
            "lncrna": False,
            "allele": False,
        }
        """
        if self.workflow_dict:
            return self.workflow_dict
        self.workflow_dict = {key: TOOLS_SELECTED.get(key) if eval(value) else "" for key, value in self.handle.items("Workflow")}
        return self.workflow_dict

    @property
    def root_path(self):
        return self.handle.get("Root", "path")

    @property
    def log_path(self):
        return os.path.join(self.root_path, "log")

    @property
    def pipe_log_dir(self):
        return os.path.join(self.log_path, "pipe")

    def pipe_log_file(self, mark):
        class P(object):
            e = os.path.join(self.pipe_log_dir, str(mark) + ".e")
            o = os.path.join(self.pipe_log_dir, str(mark) + ".o")

        class S(object):
            e = os.path.join(self.log_path, "single.e")
            o = os.path.join(self.log_path, "single.o")
        if mark == "global":
            return S
        return P

    @property
    def pools(self):
        return self.handle.get("Resource", "pools")

    @property
    def ppn(self):
        return self.handle.get("Cluster", "ppn")

    @property
    def threads(self):
        return str(int(int(self.ppn) / int(self.pools)))

    @property
    def name(self):
        return self.handle.get("Cluster", "name")

    @property
    def nodes(self):
        return self.handle.get("Cluster", "nodes")

    @property
    def walltime(self):
        return self.handle.get("Cluster", "walltime")

    @property
    def cmd_pbs(self):
        return " ".join(
            ["qsub",
             "-N", self.name,
             "-l", ":".join(["=".join(["nodes", self.nodes]),
                             "=".join(["ppn", self.ppn])
                             ]),
             "-l", "=".join(["walltime", self.walltime], ),
             "-e", self.log_path,
             "-o", self.log_path
             ]
        )

    @property
    def pbs_info(self):
        return {
            "job name": self.name,
            "node": self.nodes,
            "cpu": self.ppn,
            "walltime": self.walltime,
        }

    @property
    def submit_message(self):
        return {
            "threads": self.threads,
            "pools": self.pools,
            "output": self.root_path
        }

    @property
    def species(self):
        return self.handle.get("SampleMessage", "species")

    @property
    def phred(self):
        if self.handle.has_option("SampleMessage", "phred"):
            return self.handle.get("SampleMessage", "phred")
        return

    @property
    def genome(self):
        return self.handle.get("Genome", "genomefile")

    @property
    def annotations(self):
        return self.handle.get("Genome", "annotations")

    @property
    def annotations_gtf(self):
        if self.handle.has_option("Genome", "annotations_gtf"):
            return self.handle.get("Genome", "annotations_gtf")
        return

    @property
    def bowtie1_index(self):
        if self.handle.has_option("Genome", "bowtie1_index"):
            return self.handle.get("Genome", "bowtie1_index")
        return

    @property
    def bowtie2_index(self):
        if self.handle.has_option("Genome", "bowtie2_index"):
            return self.handle.get("Genome", "bowtie2_index")
        return

    @property
    def hisat2_index(self):
        if self.handle.has_option("Genome", "hisat2_index"):
            return self.handle.get("Genome", "hisat2_index")
        return

    @property
    def star_index(self):
        if self.handle.has_option("Genome", "star_index"):
            return self.handle.has_option("Genome", "star_index")
        return

    @property
    def library_type(self):
        if self.handle.has_option("SampleMessage", "library_type"):
            return self.handle.get("SampleMessage", "library_type")
        else:
            return

    @property
    def bed(self):
        if self.handle.has_option("Genome", "bed"):
            return self.handle.get("Genome", "bed")
        return

    @property
    def hdrs(self):
        if self.handle.has_option("Genome", "hdrs"):
            return self.handle.get("Genome", "hdrs")
        return

    @property
    def splicesites_txt(self):
        if self.handle.has_option("Genome", "hisat2_splicesites_txt"):
            return self.handle.get("Genome", "hisat2_splicesites_txt")
        return

    @property
    def dbsnp(self):
        """for workflow： variation"""
        if self.handle.has_option("Allele", "dbsnp"):
            return self.handle.get("Allele", "dbsnp")
        return

    @property
    def hla_bed(self):
        """for workflow： allele"""
        if self.handle.has_option("Allele", "hla_bed"):
            return self.handle.get("Allele", "hla_bed")
        return

    @property
    def haplo_count_bed(self):
        """for workflow： allele"""
        if self.handle.has_option("Allele", "haplo_count_bed"):
            return self.handle.get("Allele", "haplo_count_bed")
        return

    @property
    def known_lncrna_gtf(self):
        if self.handle.has_option("Lncrna", "known_lncrna_gtf"):
            return self.handle.get("Lncrna", "known_lncrna_gtf")
        return

    @property
    def starfusion_genome_resource_lib(self):
        if self.handle.get("Fusion", "starfusion_genome_resource_lib"):
            return self.handle.get("Fusion", "starfusion_genome_resource_lib")
        return

    @property
    def diff_dir(self):
        if TOOLS_SELECTED.get("differentialexpression") == "deseq2":
            return os.path.join(self.root_path, "diff", "deseq2")
        elif TOOLS_SELECTED.get("differentialexpression") == "ballgown":
            return os.path.join(self.root_path, "diff", "ballgown")
        else:
            return os.path.join(self.root_path, "diff", "edger")

    @property
    def fusion_dir(self):
        """Store sample completed fusion-search"""
        return os.path.join(self.root_path, "fusion", TOOLS_SELECTED.get("fusion"))

    @property
    def fusion_bam(self):
        """Store sample unfinished fusion-search"""
        return os.path.join(self.fusion_dir, "bam")

    @property
    def ase_result_dir(self):
        return os.path.join(self.root_path, "allele")

    @property
    def merged_gtf(self):
        return os.path.join(self.root_path, "merged.gtf")

    @property
    def compare_gtf_prefix(self):
        return os.path.join(self.root_path, "origin")

    @property
    def lncrna_gtf(self):
        return os.path.join(self.root_path, "lncrna.gtf")
    # @property
    # def variation_result_dir(self):
    #     """for samtools merge and gatk merge"""
    #     return os.path.join(self.root_path, "variation", TOOLS_SELECTED.get("variation"))
    #
    # @property
    # def variation_bam(self):
    #     return creat_dir(os.path.join(self.variation_result_dir, "bam"))
    #
    # @property
    # def variation_vcf(self):
    #     return creat_dir(os.path.join(self.variation_result_dir, "vcf"))


