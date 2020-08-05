#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" 
    SNP INDEL
    1> Indel局部区域重比对（可选），如果后续用GATK HaplotypeCaller，可免去这一步骤
    2> 重新校正碱基质量值（BQSR），矫正系统错误，仅仅GATK可以做，需要下载相关的文件
    3> call SNP
    4> 变异检测质控和过滤（VQSR）
        filter SNP（default QUAL > 20 || DP > 10）
        - QUAL ：as the "Phred-scaled quality score for the assertion made in ALT.
            i.e. −10log10 prob(call in ALT is wrong).
            If ALT is ‘.’ (no variant) then this is −10log10 prob(variant),
            and if ALT is not ‘.’ this is −10log10 prob(no variant).
            If unknown, the missing value should be specified."
        - DP : combined depth across samples, e.g. DP=154
        - QUAL = 20 means 1 out of 100(1/100)
        - DP: 最好是去除大于2倍的平均深度（或者中位数）的SNP
        gatk 步骤较为繁琐，包含VariantRecalibrator， ApplyRecalibration

PS：
    * -DP threshold 受本身测序深度深度影响，低测序深度中样本中snp存在很多误报
    * samtools和GATK的流程区别很大，samtools无法生成gvcf的文件，对于生物学重复，也不可以直接merge VCF，所以首先要merge bam文件，相反GATK每个样本可以生成gvcf
        - samtools
            -bcftools -i 和 -e参数分别是包含和去除，分清楚
        -GATK

有待改进方案：
    1> 每个染色体（独立，除非结构性变异）单独call snp，然后整合在一起[节省时间]
    2> 包含合并bam然后call的类[研究群体变异的时候合并一个群体的vcf]
问题：
    1> 都是先合并然后过滤的
        We don't recommend using the GVCF for final analysis; it is meant to be an intermediate file. You should use GenotypeGVCFs to generate a final VCF which you can use in the filtering step. As for filtering, you have two options: VQSR or hard filtering.

"""
import os

from state2 import workflow_stateful
from core.abstract import AbstractWorkflow, AbstractTools, AbstractRun

from core.decorator import runshell
from params.tools import ParamsVariation
from source import gatk


class SamtoolsRun(AbstractRun):
    genome_prefix = "-f"

    @property
    def bp(self) -> bool:
        if not self.is_file_empty(self.dm.vcf):
            return False
        return True

    def clear(self):
        pass

    @runshell
    def run(self):
        tmp_vcf = os.path.join(self.dm.variation_dir, "org.vcf.gz")
        cmd_1 = " ".join(("samtools",
                          "mpileup -u",
                          self.genome,
                          self.dm.bam_file,
                          "|",
                          "bcftools",
                          "call -mv",
                          "-O z",
                          "-o", tmp_vcf
                          ))
        cmd_2 = " ".join(("bcftools",
                          "filter",
                          tmp_vcf,
                          "-s FAIL",
                          ParamsVariation.bcftools_filter_params,
                          "-O z",
                          "-o", self.dm.vcf
                          ))
        cmd_3 = " ".join(("tabix",
                          "-p vcf",
                          self.dm.vcf
                          ))
        return cmd_1, cmd_2, cmd_3


class GatkRun(AbstractRun):
    genome_prefix = "-R"

    @property
    def bp(self) -> bool:
        if not self.is_file_empty(self.dm.vcf):
            return False
        return True

    def clear(self):
        pass

    @runshell
    def run(self):
        tmp_vcf = os.path.join(self.dm.variation_dir, "org.vcf.gz")
        cmd_1 = " ".join((gatk,
                        "HaplotypeCaller",
                        ParamsVariation.gatk_call_params,
                        "-I", self.dm.bam_file,
                        self.genome,
                        self.ini.dbsnp if self.ini.dbsnp else "",
                        "-O", tmp_vcf
                        ))
        cmd_2 = " ".join((gatk,
                          "VariantFiltration",
                          self.genome,
                          "-V", tmp_vcf,
                          "--filter-name", "PASS",
                          ParamsVariation.gatk_filter_params,
                          "-O", self.dm.vcf
                          ))
        return cmd_1, cmd_2


@workflow_stateful
class Variation(AbstractWorkflow):
    """alter splice"""

    class Samtools(AbstractTools):
        """call snp"""

        class Run(SamtoolsRun):
            pass

    class Gatk(AbstractTools):
        """call snp"""

        class Run(GatkRun):
            pass
