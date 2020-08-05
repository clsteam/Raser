#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author  : Yao

TOOL_ROOT_PATH = "/public/home/yxu/tools"

# The tool is used as a guideline
# All strings must be lowercase
TOOLS_SELECTED = {
    "qualitycontrol": "fastqc",
    "trim": "trimmomatic",
    "alignment": "tophat2",  # tophat2, hisat2, star
    "rmdup": "samtools", # samtools, picard
    "genecount": "featurecounts",  # htseq, featurecounts, star

    "strandspecific": "",

    "transcript": "stringtie",  # cufflinks, stringtie
    "variation": "gatk",  # samtools, gatk

    "differentialexpression": "deseq2",  # ballgown, deseq2, edger
    "altersplice": "asprofile",
    "fusion": "tophatfusion",  # tophatfusion, starfusion
    "lncrna": "cc",  # cc
    "allele": "phaser",
}

# Reads保留的最短长度
MINLEN = 50
# default Read-Group platform (e.g. ILLUMINA, SOLID, LS454, HELICOS and PACBIO)
RGPL = "ILLUMINA"
# 是否将GTF格式作为流程首选，默认 False (GTF兼容性更好，尤其在STAR建立索引的时候)
PRIMARY_GTF_ANNOTATIONS = False
# 双端数据单某一端的质量很差，可以保留质量高的一端进行单端分析
WHETHER_PE_TO_SE = True
# 比对时是否加入参考注释，默认 True
WHETHER_ALIGNMENT_WITH_ANNOTATIONS = True
# 保留仅标记或者去掉PCR重复(仅对picard生效)，默认 True
WHETHER_MARK_DUPLICATES_ONLY = True
# 自动检测链特异性并使用，默认 False (再次比对将花费大量时间)
STRAND_SPECIFIC_USE_AUTOMATICALLY = False
# 即使没有对照组样本，强制组装转录本，默认 False
ENFORCE_ASSEMBLY = False

