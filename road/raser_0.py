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
import itertools
import os
import glob
import re
from functools import reduce
from logging import Logger

from core.ilog import Ini
from core.decorator import runshell
from core.exception import CustomException
from source import fastq_dump


class PickDicts(object):
    """
    Identify the input sample folder.
    """

    def __init__(self, ini):
        self.ini = ini

    @property
    def sp(self) -> list:
        """Allow input of matching characters *
        such as:
            /public/home/xx/lung-seq/normal/*
            /public/home/xx/lung-seq/tumor/*
        """
        sp1 = reduce(lambda x, y: x + y,
                     [glob.glob(d) for d in self.ini.experiment_group]) if self.ini.experiment_group else []
        sp2 = reduce(lambda x, y: x + y,
                     [glob.glob(d) for d in self.ini.control_group]) if self.ini.control_group else []
        return sp1 + sp2

    @property
    def data_hash(self) -> dict:
        """
        :return: one merged data_hash.
            data_hash:
                format:[((file_name,), home_dir),]
        """
        two_layer_list = map(self._redundancy, self.sp)
        return self._determine_the_format(list(itertools.chain(*two_layer_list)))

    @staticmethod
    def _redundancy(dic) -> list:
        """
        Remove the decompressed file of the sra file and folder.
        :param dic: a folder
        :return: a list of original sample file
        """
        dict_list = os.listdir(dic)
        dict_list_copy = dict_list[:]
        for file_name in dict_list_copy:
            # Remove dir,only leave file
            if os.path.isdir(os.path.join(dic, file_name)):
                dict_list.remove(file_name)
            # Remove the decompressed file of the sra file
            if file_name.endswith(".sra"):
                decompressed_file = glob.glob(os.path.join(dic, os.path.splitext(file_name)[0] + "*fastq*"))
                if decompressed_file:
                    for abs_path in decompressed_file:
                        dict_list.remove(os.path.basename(abs_path))

        del dict_list_copy
        return list(os.path.join(dic, name) for name in dict_list)

    def _determine_the_format(self, file_list) -> dict:
        """
        :param file_list: a list of original sample file
        :return:
            data_hash:
                format:[((file_name,), home_dir),]
            {
                "SRA":[(('/PRJ1/SRR1.sra',),'/PRJ1/SRR1')],
                "PE":[(('/PRJ1/SRR21_1.fq.gz', '/PRJ1/SRR21_2.fq.gz'),'/PRJ1/SRR21'),
                (('/PRJ1/SRR31_1.fastq', '/PRJ1/SRR31_2.fastq'),'/PRJ1/SRR31')],
                "SE":[(('/PRJ1/SRR4.fq',),'/PRJ1/SRR4')]
            }
        """
        self.ini.main_logger.spam(
            "determine format of the input file(PE or SE) [SRA or fq or fq.gz or fastq or fastq.gz].")
        data_hash = dict()
        for file_name in file_list:
            base_name = os.path.basename(file_name)
            if base_name.endswith(".sra"):
                home_dir = os.path.splitext(file_name)[0]
                data = (file_name,)
                data_hash.setdefault("SRA", set()).add((data, home_dir))
            elif "_1." in base_name or "_2." in base_name:  # By default, _1.fq.gz and _2.fq.gz is PE
                pat = re.compile('^(.*?)_[1|2](\.f[|ast]q.*)')
                home_dir, suffix = re.search(pat, file_name).groups()
                data = ("".join((home_dir, "_1", suffix)), "".join((home_dir, "_2", suffix)))
                data_hash.setdefault("PE", set()).add((data, home_dir))
            else:
                home_dir = os.path.splitext(file_name)[0]
                data = (file_name,)
                data_hash.setdefault("SE", set()).add((data, home_dir))
        return data_hash


def sequence_format_standardized(logger: Logger, params, ini: Ini):
    """Sra file decompression and determine its sequencing format.
    :param ini:
    :param logger:
    :param params: Parameters of 'sequence_format', 'sample' and 'home_dir'are required.
        eg:
            * params["sequence_format"]
            * params["sample"] --> (file_name,)
            * params["home_dir"]
    """

    @runshell
    def _sra_trans_fq(ini=ini, accession="", name="SRAtoolkit"):
        if not glob.glob(os.path.join(params.get("home_dir"), "*f*q.gz")):
            logger.info("{0} decompressing...".format(params["sample"][0]))
            return " ".join((fastq_dump,
                             "--split-files",
                             "--gzip",
                             params["sample"][0],
                             "-O", os.path.dirname(params["home_dir"])
                             ))
        else:
            return

    if params.get("sequence_format") == "SRA":
        origin_fq = lambda: glob.glob(params.get("home_dir") + "*.fastq.gz")
        clean_fq = lambda: glob.glob(os.path.join(params.get("home_dir"), "*_clean.fq.gz"))
        trans_judgement = [origin_fq(), clean_fq(), glob.glob(os.path.join(params.get("home_dir"), "*.bam"))]
        if not any(trans_judgement):
            _sra_trans_fq(ini=ini, accession=os.path.basename(params.get("home_dir")), name="SRAtoolkit")

        if len(origin_fq()) == 1 or len(clean_fq()) == 1:
            params["sequence_format"] = "SE"
            params["sample"] = (".".join((params.get("home_dir"), "fastq.gz")),)
        elif len(origin_fq()) == 2 or len(clean_fq()) == 2:
            params["sequence_format"] = "PE"
            params["sample"] = (
                "".join((params.get("home_dir"), "_1.fastq.gz")), "".join((params.get("home_dir"), "_2.fastq.gz")))
        else:
            raise CustomException("{0} unable to recognize multiple compressed result files".format(clean_fq))

    return params
