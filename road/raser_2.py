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
import os
import re
import glob

# Private module
from state2 import workflow_stateful
from core.abstract import AbstractWorkflow, AbstractTools, AbstractRun
from core.decorator import runshell
from raser.settings import WHETHER_PE_TO_SE, MINLEN
from source import trimmomatic_jar
from params.tools import ParamsTrim


class TrimmomaticRun(AbstractRun):
    phred_prefix = "-"
    threads_prefix = "-threads"

    @property
    def bp(self) -> bool:
        if os.path.exists(self.dm.bam_file):
            return False
        __clean_fq__ = glob.glob(os.path.join(self.dm.home_dir, "*_clean.fq.gz"))
        if __clean_fq__:
            return False
        return True

    def clear(self):
        for fq in self.dm.sample:
            if os.path.exists(fq):
                os.remove(fq)
        for fq in self.dm.unclean_data:
            if os.path.exists(fq):
                os.remove(fq)

    @runshell
    def run(self):
        """
        :return:commands of trimmomatic
        """
        refactoring_params = self._parameters
        cmd = " ".join(("java",
                        "-jar",
                        trimmomatic_jar,
                        self.dm.sequence_format,
                        self.phred,
                        self.threads,
                        " ".join(self.dm.sample),
                        self.output_data,
                        refactoring_params
                        ))
        return cmd

    @property
    def output_data(self) -> str:
        return " ".join(map(" ".join, zip(self.dm.clean_data, self.dm.unclean_data)))

    @property
    def _parameters(self) -> str:
        """
        :return: parameters
        """
        self._generate_original_adapter()
        params_5 = self._crop_params
        params_1 = ":".join(("ILLUMINACLIP", self.dm.adapter_file, "2:30:10"))
        params_2 = ":".join(("LEADING", ParamsTrim.leading))
        params_3 = ":".join(("TRAILING", ParamsTrim.trailing))
        params_4 = ":".join(("SLIDINGWINDOW", ParamsTrim.slidingwindow))
        return " ".join((params_1, params_2, params_3, params_4, params_5))

    @property
    def _crop_params(self):
        """
        Get HEADCROP, CROP and Overexpression reads.
        :return: string hash
        """
        crop_list = []
        has_orepresented = []
        for fq_txt in self.input_fq_txt:
            with open(self.dm.adapter_file, "a+") as doc_w, open(fq_txt, "rt") as doc_r:
                data = doc_r.readlines()
                self.dm.max_read_length = int(data[8].strip('\n').split('\t')[1].split('-')[-1])
                a = 0
                next_step = 0
                per_index = 0
                for a in range(len(data)):
                    # >>Per base sequence content\tfail\n
                    k = re.match(re.compile('^>>Per base sequence content\t(.*)'), data[a])
                    if k:
                        if k.group(1) == 'pass':
                            break
                        else:
                            next_step = 1
                            per_index = a + 2
                    if next_step == 1 and data[a] == '>>END_MODULE\n':
                        quality_content = data[per_index:a]
                        # Get HEADCROP, CROP
                        postion = []
                        base_quality = ""
                        threshold = ParamsTrim.threshold_per_base_content
                        # base_quality: generate a string to judge which reads should be reserved
                        for i in range(len(quality_content)):
                            atgc_quality = quality_content[i].strip('\n').split('\t')
                            [start, end] = atgc_quality[0].split("-") if re.match(re.compile(".*-.*"),
                                                                                  atgc_quality[0]) else [
                                atgc_quality[0],
                                None]
                            # if the difference between A and T, or G and C is greater than 10% in any position.search
                            # Base	G	A	T	C
                            postion.append((start, end))  # STR convert INT
                            if (abs(float(atgc_quality[1]) - float(atgc_quality[4])) < 100 * threshold) and (
                                    abs(float(atgc_quality[2]) - float(atgc_quality[3])) < 100 * threshold):
                                base_quality += "o"
                            else:
                                base_quality += "x"

                        # print(base_quality)

                        # get CROP,HEADCROP
                        max_str = max(base_quality.split("x"))
                        head_crop = postion[base_quality.find(max_str)]
                        head_crop = head_crop[0] if head_crop[0] else head_crop[1]
                        crop = postion[base_quality.find(max_str) + len(max_str) - 1]
                        crop = crop[1] if crop[1] else crop[0]

                        crop_list.append((int(head_crop), int(crop)))

                        break
                # >>Overrepresented sequences\twarn\n
                next_step = 0
                over_index = 0
                for b in range(a, len(data)):
                    k = re.match(re.compile('^>>Overrepresented sequences\t(.*)'), data[b])
                    if k:
                        if k.group(1) == 'pass':
                            break
                        else:
                            next_step = 1
                            over_index = b + 2
                    if next_step == 1 and data[b] == '>>END_MODULE\n':
                        over_data = data[over_index:b]

                        has_orepresented.append(False)

                        for line in over_data:
                            if line.split('\t')[-1] != 'No Hit\n':
                                has_orepresented[-1] = True
                                doc_w.writelines('>' + line.split('\t')[-1] + line.split('\t')[0] + '\n')
                        break

        if not len(crop_list):
            return ""

        if self.dm.sequence_format == "PE":
            if len(crop_list) == 1 or crop_list[0] != crop_list[1]:
                # 双端质量不一致(len=1 or len=2)
                self.logger.notice("[{0}] QUALITY inconsistent at both pair-ends reads, crop_list={1}".format(self.dm.id, crop_list))

            if len(crop_list) == 2:
                minlen = MINLEN
                is_filter_all_list = [True if (crop - head_crop) < MINLEN else False for head_crop, crop in crop_list]

                if any(is_filter_all_list):
                    if all(is_filter_all_list):
                        minlen = int(ParamsTrim.minlen)
                        self.logger.error("[{0}] both-ends: The length of the longest reads < {1}(MINLEN), and set MINLEN={2}".format(self.dm.id, MINLEN, minlen))
                    else:
                        self.logger.error("[{0}] single-ends: The length of the longest reads < {1}(MINLEN), and you can set WHETHER_PE_TO_SE=True".format(self.dm.id, minlen))
                        if WHETHER_PE_TO_SE:
                            self.logger.notice(
                                "[{0}] single-ends: truncate one end for subsequent single-end analysis".format(self.dm.id))
                            self.dm.sequence_format = "SE"
                            self.dm.clean_data = [
                                "".join((os.path.join(self.dm.home_dir, self.dm.id), "_clean.fq.gz")), ]
                            self.dm.unclean_data = ("",)
                            if is_filter_all_list[0] is True:
                                self.dm.sample = (self.dm.sample[1],)
                                crop_list.pop(0)
                            else:
                                self.dm.sample = (self.dm.sample[0],)
                                crop_list.pop(1)
                        else:
                            minlen = int(ParamsTrim.minlen)
                            self.logger.error(
                                "[{0}] single-ends: The length of the longest reads < {1}(MINLEN), and set MINLEN={2}".format(
                                    self.dm.id, MINLEN, minlen))

        return "HEADCROP:{headcrop} CROP:{crop} MINLEN:{minlen}".format(headcrop=crop_list[0][0], crop=crop_list[0][1], minlen=minlen)

    def _generate_original_adapter(self):
        """Write the standard adapter in ParamsTrim into our file."""
        with open(self.dm.adapter_file, "w") as _:
            pass
        with open(self.dm.adapter_file, "a+") as doc_w:
            # config.get("ParamsTrim", "SE_adapters")
            for trimmomatic_adapter_suffix in getattr(ParamsTrim, "_".join((self.dm.sequence_format.lower(), "adapters"))).split(","):
                adapter_file = os.path.join(os.path.dirname(trimmomatic_jar),
                                            "adapters", trimmomatic_adapter_suffix)
                if os.path.exists(adapter_file):
                    with open(adapter_file, "r") as doc_r:
                        doc_w.write(doc_r.read() + "\n")
                else:
                    self.logger.error("ParamsTrim[adapters] {0} couldn't be found or expired!".format(adapter_file))

    @property
    def input_fq_txt(self) -> list:
        """返回fastq解压结果的txt文件"""
        if self.dm.sequence_format == "SE":
            return [os.path.join(self.dm.qc_dir, self.dm.id + "_fastqc", "fastqc_data.txt")]
        return [os.path.join(self.dm.qc_dir, self.dm.id + "_" + str(i) + "_fastqc", "fastqc_data.txt") for i in
                range(1, 3)]


@workflow_stateful
class Trim(AbstractWorkflow):
    """Trim sequence"""

    def entrance(host):
        host.tool_instance.dm.adapter_file = os.path.join(host.dm.qc_dir, "adapter.fa")
        __clean_fq__ = glob.glob(os.path.join(host.dm.home_dir, "*_clean.fq.gz"))
        if len(__clean_fq__) == 1 and host.dm.sequence_format == "PE":
            host.logger.notice(
                "[{0} {1}:single-ends] truncate one end for subsequent single-end analysis".format(host.dm.id, "PE"))
            host.dm.sequence_format = "SE"
            host.dm.clean_data = ["".join((os.path.join(host.dm.home_dir, host.dm.id), "_clean.fq.gz")), ]
            host.dm.unclean_data = ("",)
            # host.dm.sample = (host.dm.sample[1],)

    class Trimmomatic(AbstractTools):

        class Run(TrimmomaticRun):
            pass
