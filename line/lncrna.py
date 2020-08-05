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
import re
from functools import reduce

from state2 import workflow_stateful
from core.abstract import AbstractWorkflow, AbstractTools, AbstractRun
from core.decorator import runshell
from core.r_rule import conversion_function_symbol

from raser.settings import TOOLS_SELECTED

from params.limit import SpeciesLimit

from source import gtf_to_fa, cpc_dir, cnci_dir, cuffcompare, gffcompare, cpat_dir, pfam_scan, pfam_db

from .r_script import R_OF_LNCFINDER


class CcRun(AbstractRun):
    """Global"""
    annotations_prefix = "-r"
    genome_prefix = "-s"
    threads_prefix = "-p"  # CNCI

    def __init__(self, *args, **kwargs):
        self.lnc_compare_tmap = ""
        self.potential_lnc_gtf = ""
        self.potential_lncid = set()
        self.known_lncid = set()

        self.cpc_res = os.path.join(args[1].nc_dir, "CPC.txt")
        self.cnci_res = os.path.join(args[1].nc_dir, "CNCI.index")
        self.cpat_res = os.path.join(args[1].nc_dir, "CPAT")
        self.lncfinder_res = os.path.join(args[1].nc_dir, "LNCFINDER")

        super(CcRun, self).__init__(*args, **kwargs)

    @property
    def bp(self) -> bool:
        return True

    def clear(self):
        if os.path.exists(self.dm.nc_dir + ".log"):
            os.remove(self.dm.nc_dir + ".log")
            os.removedirs(self.dm.nc_dir + "_Tmp_Dir")

    @runshell
    def compare(self, out_prefix, in_gtf=None, ref_gtf=None):
        """
        cuffcompare/gffcompare
        :param out_prefix:
        :param in_gtf: default: merged.gtf
        :param ref_gtf: default: genome annotations
        :return: cmd
        """
        self.lnc_compare_tmap = out_prefix + ".merged.gtf.tmap"
        if os.path.exists(self.lnc_compare_tmap):
            return
        if TOOLS_SELECTED.get("transcript") == "cufflinks":
            # cufflinks
            cmd = " ".join((cuffcompare,
                            in_gtf if in_gtf else self.ini.merged_gtf,
                            " ".join((self.annotations_prefix, ref_gtf)) if ref_gtf else self.annotations,
                            self.genome,
                            "-C",
                            "-o", out_prefix
                            ))
            return cmd
        else:
            # stringTie
            cmd = " ".join((gffcompare,
                            in_gtf if in_gtf else self.ini.merged_gtf,
                            " ".join((self.annotations_prefix, ref_gtf)) if ref_gtf else self.annotations,
                            self.genome,
                            "-o", out_prefix
                            ))
            return cmd

    def before(self):
        """
        1> filtering
            mrna -> known lnc + others mrna
            mrna -> potential lnc (function: filter_lnc_bytmap())
        :return:
        """
        filter_params = dict(self.ini.handle.items("Lncrna"))
        self.potential_lnc_gtf = os.path.join(os.path.dirname(self.ini.merged_gtf), "lnc_potential.gtf")
        # if os.path.exists(self.potential_lnc_gtf):
        #     return
        if self.ini.known_lncrna_gtf:
            lnc_known_gtf = os.path.join(os.path.dirname(self.ini.merged_gtf), "lnc_known.gtf")
            lnc_unknown_gtf = os.path.join(os.path.dirname(self.ini.merged_gtf), "lnc_unknown.gtf")

            self.compare(os.path.splitext(lnc_known_gtf)[0], ref_gtf=self.ini.known_lncrna_gtf)
            self.known_lncid = self.filter_lnc_bytmap(self.lnc_compare_tmap, lnc_known_gtf, class_code="=c")
            self.filter_lnc_bytmap(self.lnc_compare_tmap, lnc_unknown_gtf, class_code="!=c", **filter_params)

            self.compare(os.path.splitext(self.potential_lnc_gtf)[0], in_gtf=lnc_unknown_gtf)
            self.potential_lncid = self.filter_lnc_bytmap(self.lnc_compare_tmap, self.potential_lnc_gtf, class_code="uxijo", **filter_params)
        else:
            self.compare(self.ini.compare_gtf_prefix)
            self.potential_lncid = self.filter_lnc_bytmap(self.lnc_compare_tmap, self.potential_lnc_gtf, class_code="uxijo", **filter_params)

    @runshell
    def run(self):
        """
        potential lnc -->(identification)--> novel lnc
        2> convert .gtf to .fa
        3> cpc2
        :return:
        """
        # if os.path.exists(self.cnci_res):
        #     return
        lnc_fasta = os.path.join(self.dm.nc_dir, "lnc_potential.fasta")
        cmd_0 = " ".join((gtf_to_fa,
                          self.potential_lnc_gtf,
                          self.genome_suffix,
                          ">", lnc_fasta
                          ))

        cmd_1 = " ".join((os.path.join(cpc_dir, "bin", "CPC2.py"),
                          "-i", lnc_fasta,
                          "-o", os.path.splitext(self.cpc_res)[0]
                          ))

        cmd_2 = " ".join((os.path.join(cnci_dir, "CNCI.py"),
                          self.ini.category,
                          self.threads_all,
                          "-f", lnc_fasta,
                          "-o", self.dm.nc_dir
                          ))
        # [error code == 1]
        cmd_3 = " ".join((os.path.join(cnci_dir, "filter_novel_lincRNA.py"),
                          "-i", self.cnci_res,
                          "-g", self.potential_lnc_gtf,
                          "-o", os.path.join(self.dm.nc_dir, "cnci")
                          ))
        cmd_4 = ""
        if self.ini.species in SpeciesLimit.cpat:
            param_d = os.path.join(cpat_dir, "dat", SpeciesLimit.cpat.get(self.ini.species) + "_logitModel.RData")
            param_x = os.path.join(cpat_dir, "dat", SpeciesLimit.cpat.get(self.ini.species) + "_Hexamer.tsv")
            cmd_4 = " ".join((os.path.join(cpat_dir, "bin", "cpat.py"),
                              "-g", lnc_fasta,
                              "-d", param_d,
                              "-x", param_x,
                              "-o", os.path.join(self.dm.nc_dir, "CPAT")
                              ))

        lncfinder_rscript = os.path.join(self.dm.nc_dir, "lncfinder.R")
        with open(lncfinder_rscript, "w") as handle:
            handle.write(conversion_function_symbol(R_OF_LNCFINDER.format(fasta=lnc_fasta)))
        cmd_5 = " ".join(("Rscript", lncfinder_rscript))

        return cmd_0, cmd_1, cmd_2, cmd_3, cmd_4, cmd_5, 1

    def after(self):
        cpc = self.filter_column(self.cpc_res, -1, "=='noncoding'", start=1)
        cnci = self.filter_column(self.cnci_res, 1, "=='noncoding'", start=1)
        cpat = self.filter_column(self.cpat_res, -1, "<0.5", start=1)
        intersect = reduce(lambda a, b: a & b if b else a, [cnci, cpc, cpat])

        with open(os.path.join(self.dm.nc_dir, "lnc_predict.statistics"), "w") as handle:
            handle.write("software\tpredicted_number\n")
            handle.write("known\t{0}\n".format(len(self.known_lncid)))
            handle.write("potential\t{0}\n".format(len(self.potential_lncid)))
            handle.write("cpc\t{0}\n".format(len(cpc)))
            handle.write("cnci\t{0}\n".format(len(cnci)))
            handle.write("cpat\t{0}\n".format(len(cpat)))
            handle.write("intersect\t{0}\n".format(len(intersect)))

        self.filter_lnc_byid(intersect | self.known_lncid, self.ini.lncrna_gtf)

    def filter_lnc_byid(self, lncid: set, out_gtf):
        with open(self.ini.merged_gtf, "r") as handle_r, open(out_gtf, "w") as handle_w:
            for line in handle_r:
                if line.startswith("#"):
                    handle_w.write(line)
                    continue
                exist_groups = re.search("transcript_id\s\"(.\S+?)\"", line).groups()
                if exist_groups:
                    transcript_id = exist_groups[0]
                    if transcript_id in lncid:
                        handle_w.write(line)

    @runshell
    def pfam(self):
        fasta = os.path.join(self.dm.nc_dir, "lncrna.fasta")
        pfam_fa = os.path.join(self.dm.nc_dir, "lncrna_pfam.fasta")
        cmd_5 = " ".join((gtf_to_fa,
                          self.ini.lncrna_gtf,
                          self.genome_suffix,
                          ">", fasta
                          ))

        cmd_6 = " ".join((pfam_scan,
                          "-fasta", fasta,
                          "-dir", pfam_db,
                          "-as",
                          "-cpu", self.ini.ppn,
                          "-outfile", pfam_fa
                          ))
        return cmd_5, cmd_6

    def filter_lnc_bytmap(self, in_tmap, out_gtf, class_code="!uxijo", min_length=0, min_cov=0, min_fpkm=0, **kwargs) -> set:
        """
        :param class_code: if class code startswith "!", this means exclude those class code, on the contrary it's conclude
        """
        min_length, min_cov, min_fpkm = list(map(int, [min_length, min_cov, min_fpkm]))
        filtered_lnc = []
        ccf = (lambda x: x not in class_code) if class_code.startswith("!") else (lambda x: x in class_code)
        with open(in_tmap, "r") as handle_r1, open(self.ini.merged_gtf, "r") as handle_r2, open(out_gtf, "w") as handle_w:
            handle_r1.readline()
            for line in handle_r1:
                ls = line.strip().split()
                # class_code, cuff_id, fpkm, cov, length
                #      2         4       6    8      9
                if ccf(ls[2]) and float(ls[6]) >= min_fpkm and float(ls[8]) >= min_cov and float(ls[9]) >= min_length:
                    filtered_lnc.append(ls[4])
            if len(filtered_lnc) < 50:
                self.logger.warning("After filtering by condition '{0}', only {1} terms in '{2}' are retained.[]".format(class_code, len(filtered_lnc), in_tmap))
            handle_w.write(
                "# filtered paramters: class_code={0}, fpkm>={1}, cov>={2}, length>={3}\n".format(class_code, min_fpkm, min_cov, min_length))
            for line in handle_r2:
                if line.startswith("#"):
                    handle_w.write(line)
                    continue
                exist_groups = re.search("transcript_id\s\"(.\S+?)\"", line).groups()
                if exist_groups:
                    transcript_id = exist_groups[0]
                    if transcript_id in filtered_lnc:
                        handle_w.write(line)
        return set(filtered_lnc)

    @staticmethod
    def filter_column(doc, column, pattern, start=0, sep=None) -> set:
        """
        :param doc: file
        :param column: column of screening factors
        :param pattern: filter pattern
        :param start: from which line to read
        :param sep: separator
        :return:
        """
        res = set()
        if not os.path.exists(doc):
            return res
        with open(doc, "r") as handle:
            for i in range(start):
                handle.readline()
            for line in handle:
                ls = line.strip().split(sep)
                key = ls[0]
                col = ls[column]
                if "'" in pattern or '"' in pattern:
                    col = "'" + col + "'"
                if eval(col + pattern):
                    res.add(key)
        return res


@workflow_stateful
class Lncrna(AbstractWorkflow):
    """LncRNA"""

    class Cc(AbstractTools):
        class Run(CcRun):
            pass
