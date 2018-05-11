#!/usr/bin/env python
# -*- coding: UTF-8 -*-
__author__ = 'clsteam'

from bag.configparser import *
from bag.Cluster import *
import subprocess
import os
import sys

class Testing_environment:

    def __init__(self):
        self.devnull=open(os.devnull, 'w')
        self.color=Color()
        configparser_object = configparser()
        self.different_gene_expression=configparser_object.get_bool("Operation-switch","different-gene-expression")
        self.allele_specific_expression=configparser_object.get_bool("Operation-switch","allele-specific-expression")
        self.alternative_splicing=configparser_object.get_bool("Operation-switch","alternative-splicing")
        self.Fusion_gene_analysis=configparser_object.get_bool("Operation-switch","Fusion-gene-analysis")
        self.ncRNA_detection_and_analysis=configparser_object.get_bool("Operation-switch","ncRNA-detection-and-analysis")

    def check_dependency(self,name):
        error_code = subprocess.call("which " + name, shell=True, stdout=self.devnull);
        if error_code == 0:
            return (True);
        else:
            return (False);

    def check(self):
        # check for external dependencies
        if self.check_dependency("java") == False: self.color.fatal_error("External dependency 'java' not installed.");
        if self.check_dependency("fastqc") == False: self.color.fatal_error("External dependency 'fastqc' not installed.");
        if self.check_dependency("samtools") == False: self.color.fatal_error("External dependency 'samtools' not installed.");
        if self.check_dependency("tophat2") == False: self.color.fatal_error("External dependency 'tophat2' not installed.fusion search will not be running");

        if self.different_gene_expression:
            if self.check_dependency("cufflinks") == False: self.color.fatal_error("External dependency 'cufflinks' not installed.");

        if self.allele_specific_expression:
            if self.check_dependency("bcftools") == False: self.color.fatal_error("External dependency 'bcftools' not installed.");
            if self.check_dependency("bgzip") == False: self.color.fatal_error("External dependency 'bgzip' not installed.");
            if self.check_dependency("tabix") == False: self.color.fatal_error("External dependency 'tabix' not installed.");
            if self.check_dependency("bedtools") == False: self.color.fatal_error("External dependency 'bedtools' not installed.");

        if self.alternative_splicing:
            if self.check_dependency("cufflinks") == False: self.color.fatal_error("External dependency 'cufflinks' not installed.");


        if self.Fusion_gene_analysis:
            if self.check_dependency("bowtie") == False: self.color.fatal_error("External dependency 'bowtie' not installed.");

        if self.ncRNA_detection_and_analysis:
            if self.check_dependency("formatdb") == False: self.color.fatal_error("External dependency 'formatdb(blast)' not installed.");
            if self.check_dependency("run_predict.sh") == False: self.color.fatal_error("External dependency 'run_predict.sh(cpc)' not installed.");


