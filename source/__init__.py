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

from raser.settings import TOOL_ROOT_PATH

# sratoolkit
fastq_dump = os.path.join(TOOL_ROOT_PATH, "sratoolkit.2.8.2-centos_linux64/bin/fastq-dump")

# FastQC
fastqc = os.path.join(TOOL_ROOT_PATH, "FastQC", "fastqc")

# Trimmomatic
trimmomatic_jar = os.path.join(TOOL_ROOT_PATH, "Trimmomatic-0.36", "trimmomatic-0.36.jar")

# tophat2, hisat2, STAR
bowtie1_build = os.path.join(TOOL_ROOT_PATH, "bowtie-1.0.0", "bowtie-build")
bowtie2_build = os.path.join(TOOL_ROOT_PATH, "bowtie2-2.3.5.1-linux-x86_64", "bowtie2-build")

hisat2_build = os.path.join(TOOL_ROOT_PATH, "hisat2-2.1.0", "hisat2-build")
hista2_splice_sites = os.path.join(TOOL_ROOT_PATH, "hisat2-2.1.0", "hisat2_extract_splice_sites.py")

tophat2 = os.path.join(TOOL_ROOT_PATH, "tophat-2.1.0.Linux_x86_64", "tophat2")
tophat_fusion_post = os.path.join(TOOL_ROOT_PATH, "tophat-2.1.0.Linux_x86_64", "tophat-fusion-post")

hisat2 = os.path.join(TOOL_ROOT_PATH, "hisat2-2.1.0", "hisat2")
star = os.path.join(TOOL_ROOT_PATH, "STAR-2.7.2b/bin/Linux_x86_64_static/STAR")
star_fusion = os.path.join(TOOL_ROOT_PATH, "STAR-Fusion-v1.8.1/STAR-Fusion")

# picard
picard_jar = os.path.join(TOOL_ROOT_PATH, "picard.jar")

# featureCounts
featurecounts = os.path.join(TOOL_ROOT_PATH, "subread-2.0.1-Linux-x86_64/bin/featureCounts")

# GATK
gatk = os.path.join(TOOL_ROOT_PATH, "gatk-4.1.5.0", "gatk")

# cufflinks, stringtie
cufflinks = os.path.join(TOOL_ROOT_PATH, "cufflinks-2.2.1.Linux_x86_64", "cufflinks")
cuffcompare = os.path.join(TOOL_ROOT_PATH, "cufflinks-2.2.1.Linux_x86_64", "cuffcompare")
cuffmerge = os.path.join(TOOL_ROOT_PATH, "cufflinks-2.2.1.Linux_x86_64", "cuffmerge")
stringtie = os.path.join(TOOL_ROOT_PATH, "stringtie-2.0.6.Linux_x86_64", "stringtie")

# gffread, gffcompare
gffread = os.path.join(TOOL_ROOT_PATH, "gffread")
gffcompare = os.path.join(TOOL_ROOT_PATH, "gffcompare", "gffcompare")

#gtfToGenePred, genePredToBed
gtftogenepred = os.path.join(TOOL_ROOT_PATH, "gtfToGenePred")
genepredtobed = os.path.join(TOOL_ROOT_PATH, "genePredToBed")

# allele
phaser_py = os.path.join(TOOL_ROOT_PATH, "phaser/phaser/phaser.py")
phaser_ae_py = os.path.join(TOOL_ROOT_PATH, "phaser/phaser_gene_ae/phaser_gene_ae.py")

# alter splice
extract_as = os.path.join(TOOL_ROOT_PATH, "ASprofile.b-1.0.4", "extract-as")
summarize_as = os.path.join(TOOL_ROOT_PATH, "ASprofile.b-1.0.4", "summarize_as.pl")
ea_fpkm = os.path.join(TOOL_ROOT_PATH, "ASprofile.b-1.0.4", "extract-as-fpkm")

# fusion (tophat)
tophatfusion_library = os.path.join(TOOL_ROOT_PATH, "tophat-fusion")
blast_db = os.path.join(TOOL_ROOT_PATH, "blast")

# fusion (STAR)
starfusion_build_lib_pl = os.path.join(TOOL_ROOT_PATH, "STAR-Fusion-v1.8.1", "ctat-genome-lib-builder", "prep_genome_lib.pl")

# lncRNA
gtf_to_fa = os.path.join(TOOL_ROOT_PATH, "TransDecoder-v5.5.0", "util", "gtf_genome_to_cdna_fasta.pl")
cpc_dir = os.path.join(TOOL_ROOT_PATH, "CPC2-beta")
cnci_dir = os.path.join(TOOL_ROOT_PATH, "CNCI")
cpat_dir = os.path.join(TOOL_ROOT_PATH, "CPAT-1.2.4")

# conda
pfam_scan = "pfam_scan.pl"
pfam_db = os.path.join(TOOL_ROOT_PATH, "pfam_db")