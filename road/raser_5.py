#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author  : Yao
""" Transcript

This file(script) can also be imported as a module and contains the following
functions:

    * main - the main function of the script
    * function - returns the column headers of the file
"""

# Standard library
import os
import re
from operator import itemgetter
from state2 import workflow_stateful
from core.abstract import AbstractWorkflow, AbstractTools, AbstractRun
from core.decorator import runshell

from params.bundle import TRANSCRIPT_DIR
from raser.settings import TOOLS_SELECTED, ENFORCE_ASSEMBLY
from source import cufflinks, cuffmerge, cuffcompare, stringtie, gffcompare

from params.static import ballgown_data_dir
from params.tools import ParamsTranscript


class GenerateTpmArray(object):

    def after(self):
        self._extract_tpm()

    def _extract_tpm(self):
        hs = {}
        with open(self.ini.bed, "r") as handle_r:
            for line in handle_r:
                hs[line.strip().split()[3]] = ["0.0", "0.0"]
        pre_gene = ""
        with open(self.dm.gtf, "r") as handle_r:
            for line in handle_r:
                if line.startswith("#"):
                    continue
                if "ref_gene_id" in line:
                    id_mark = "ref_gene_id"
                else:
                    id_mark = "gene_id"
                ls = line.split()
                if ls[2] == "transcript":
                    pres = re.search('{0}\s"(.+?)".*FPKM\s"(.+?)"(.+?TPM\s"(.+?)"|)'.format(id_mark), line)
                    if pres:
                        res = pres.groups()
                        # ('ENSG00000223764', '0.000000', '; TPM "0.000000"', '0.000000')
                        # ('ENSG00000236679', '2.4337356785', '', None)
                        if not res[0].startswith("CUFF") and not res[0].startswith("STRG") and res[0] != pre_gene:
                            pre_gene = res[0]
                            hs[pre_gene] = (res[1], res[3] if res[3] else "0.0")
        with open(self.dm.tpm, "w") as handle_w:
            handle_w.write("{0} {1} {2}\n".format("Name", "FPKM", "TPM"))
            ordered_tuple = sorted(hs.items(), key=itemgetter(0))
            for key, value in ordered_tuple:
                handle_w.write("{0} {1}\n".format(key, " ".join(value)))


class Cufflinks1Run(AbstractRun, GenerateTpmArray):
    """cufflinks"""

    annotations_prefix = "-g"
    strand_specific_prefix = "--library-type"

    threads_prefix = "-p"

    @property
    def bp(self) -> bool:
        if not self.is_file_empty(self.dm.gtf):
            return False
        return True

    def clear(self):
        pass

    @runshell
    def run(self):
        cmd = " ".join((cufflinks,
                        "--no-update-check",
                        self.annotations_gtf,
                        self.strand_specific,
                        self.threads,
                        self.dm.bam_file,
                        "-o", self.dm.transcript_dir
                        ))
        return cmd


class Cufflinks2Run(AbstractRun):
    """Cuffmerge"""

    annotations_prefix = "-g"

    threads_prefix = "-p"

    @property
    def bp(self) -> bool:
        if ENFORCE_ASSEMBLY:
            return True
        if self.ini.has_comparison and self.is_file_empty(self.ini.merged_gtf):
            return True
        return False

    def clear(self):
        pass

    @runshell
    def run(self):
        cmd = " ".join((cuffmerge,
                        self.annotations_gtf,
                        self.threads_all,
                        self.assembly_gtf_list,
                        "-o", os.path.dirname(self.ini.merged_gtf)
                        ))
        return cmd

    def after(self):
        if not self.ini.merged_gtf.endswith("/merged.gtf"):
            os.rename(os.path.join(os.path.dirname(self.ini.merged_gtf), "merged.gtf"), self.ini.merged_gtf)

    @property
    def assembly_gtf_list(self) -> str:
        txt = os.path.join(self.ini.root_path, "assembly_gtf_list.txt")
        gtf_list = [os.path.join(x, TRANSCRIPT_DIR, os.path.basename(self.dm.gtf)) for x in
                    self.ini.accession_home.values()]
        with open(txt, "w") as handle:
            for gtf in gtf_list:
                handle.write(gtf + "\n")
        return txt


class Cufflinks3Run(AbstractRun):
    """cuffcompare"""

    annotations_prefix = "-r"
    genome_prefix = "-s"

    @property
    def bp(self) -> bool:
        if self.is_file_empty(self.ini.merged_gtf) or not self.is_file_empty(
                self.ini.compare_gtf_prefix + ".merged.gtf.tmap"):
            return False
        return True

    def clear(self):
        pass

    @runshell
    def run(self):
        cmd = " ".join((cuffcompare,
                        self.ini.merged_gtf,
                        self.annotations_gtf,
                        self.genome,
                        "-C",
                        "-o", self.ini.compare_gtf_prefix
                        ))
        return cmd


class Stringtie1Run(AbstractRun, GenerateTpmArray):
    annotations_prefix = "-G"
    strand_specific_prefix = ""

    threads_prefix = "-p"

    @property
    def bp(self) -> bool:
        if TOOLS_SELECTED.get("differentialexpression") == "ballgown":
            if not self.is_file_empty(self.dm.e_gtf):
                return False
        if not self.is_file_empty(self.dm.gtf):
            return False
        return True

    def clear(self):
        pass

    @runshell
    def run(self):
        # stringtie ~/code/test/data/test_R/tophat_out/accepted_hits.bam -G ~/genes/genome/hg19.gtf -o new/output.gtf -e -B -C cover.gtf -A abundances.tab
        additional_params = self.dm.e_gtf + " " + ParamsTranscript.stringtie_params.replace("b_param", os.path.join(ballgown_data_dir, self.dm.id))
        cmd = " ".join((stringtie,
                        self.threads,
                        self.annotations,
                        self.strand_specific,
                        self.dm.bam_file,
                        "-o", additional_params if TOOLS_SELECTED.get("differentialexpression") == "ballgown" else self.dm.gtf,
                        ))
        return cmd

    def void_strand_specific(self, library_type):
        if not library_type:
            return
        if library_type in ("fr-firststrand",):
            return "--rf"
        elif library_type in ("fr-secondstrand",):
            return "--fr"
        elif library_type in ("fr-unstranded",):
            return
        else:
            self.logger.warning("[{0}] unkown library type: {1}".format(self.dm.id, library_type))
            return


class Stringtie2Run(AbstractRun):
    annotations_prefix = "-G"

    threads_prefix = "-p"

    @property
    def bp(self) -> bool:
        if ENFORCE_ASSEMBLY:
            return True
        if self.ini.has_comparison and self.is_file_empty(self.ini.merged_gtf):
            return True
        return False

    def clear(self):
        pass

    @runshell
    def run(self):
        cmd = " ".join((stringtie,
                        "--merge",
                        self.threads_all,
                        self.annotations,
                        self.assembly_gtf_list,
                        "-o", self.ini.merged_gtf,
                        ))
        return cmd

    @property
    def assembly_gtf_list(self) -> str:
        txt = os.path.join(self.ini.root_path, "assembly_gtf_list.txt")
        gtf_list = [os.path.join(x, TRANSCRIPT_DIR, "transcripts.gtf") for x in
                    self.ini.accession_home.values()]
        with open(txt, "w") as handle:
            for gtf in gtf_list:
                handle.write(gtf + "\n")
        return txt


class Stringtie3Run(AbstractRun):
    """gffcompare"""
    annotations_prefix = "-r"
    genome_prefix = "-s"

    @property
    def bp(self) -> bool:
        if self.is_file_empty(self.ini.merged_gtf) or not self.is_file_empty(self.ini.compare_gtf_prefix + ".merged.gtf.tmap"):
            return False
        return True

    def clear(self):
        pass

    @runshell
    def run(self):
        cmd = " ".join((gffcompare,
                        self.ini.merged_gtf,
                        self.annotations,
                        self.genome,
                        "-o", self.ini.compare_gtf_prefix,
                        ))
        return cmd


@workflow_stateful
class Transcript(AbstractWorkflow):
    """Transcript construction(conclude FPKM), assembly and expression difference."""

    class Cufflinks1(AbstractTools):

        class Run(Cufflinks1Run):
            pass

    class Cufflinks2(AbstractTools):

        class Run(Cufflinks2Run):
            pass

    class Cufflinks3(AbstractTools):

        class Run(Cufflinks3Run):
            pass

    class Stringtie1(AbstractTools):
        """
        计算转录本
        """
        class Run(Stringtie1Run):
            pass

    class Stringtie2(AbstractTools):
        """
        组装
        """
        class Run(Stringtie2Run):
            pass

    class Stringtie3(AbstractTools):
        """
        compare for downstream of lncrna and AS
        """
        class Run(Stringtie3Run):
            pass
