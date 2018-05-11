#!/usr/bin/env python
# -*- coding: UTF-8 -*-
__author__ = 'clsteam'

import os
from bag.Cluster import *
from bag.configparser import *

class Fusion_search:

    def __init__(self,FAD,comparison,mark):
        self.FAD = FAD
        self.mark = mark
        self.config_object = configparser()
        self.fusion_library=self.config_object.get_values("Fusion","fusion_library")
        self.bowtie1_index=""
        self.comparison=comparison
        self.color=Color()
        if comparison:
            self.bam_dir=[self.fusion_library+"/Experimental-Group",self.fusion_library+"/Control-Group",self.fusion_library+"/Experimental-Group-Done",self.fusion_library+"/Control-Group-Done"]
        else:
            self.bam_dir=[self.fusion_library+"/Group",self.fusion_library+"/Group-Done"]

        self.pre_do()
        self.fusion_mapping()

    def fusion_mapping(self):
        if self.comparison:
            if self.Grouping()=="Experimental-Group":
                fusion_bam = self.bam_dir[0] + "/" + self.FAD.accession
            else:
                fusion_bam = self.bam_dir[1]  + "/" + self.FAD.accession
        else:
            fusion_bam = self.bam_dir[0]  + "/" + self.FAD.accession

        if self.Duplicate_detection():
            my_data=self.FAD.clean_data if self.FAD.format == 'SE' else self.FAD.clean_data[0]+' '+self.FAD.clean_data[1]
            shell_cmd='tophat -o '+fusion_bam+' -p '+self.FAD.threads+' --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search '+self.bowtie1_index+' '+my_data
            #shell_cmd='tophat -o '+fusion_bam+' -p '+args.threads+' --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search -r 0 --mate-std-dev 80 --max-intron-length 100000 --fusion-min-dist 100000 --fusion-anchor-length 13 --fusion-ignore-chromosomes chrM '+bowtie1_index+' '+my_data
            run(shell_cmd)

    def pre_do(self):
        for dir in self.bam_dir:
            if not os.path.exists(dir):
                os.makedirs(dir)
        if not self.config_object.has_option("Fusion", "bowtie1_index"):
            self.color.print_warning("No genome index(bowtie1) can be found and build the index now")
            shell_cmd="bowtie-build "+self.config_object.get_values("Genome", "genome")+" "+os.path.dirname(self.config_object.get_values("Genome", "genome"))+"/bowtie1"
            run(shell_cmd)
            self.config_object.add_values("Fusion", "bowtie1_index", os.path.dirname(self.config_object.get_values("Genome", "genome"))+"/bowtie1")

        self.bowtie1_index = self.config_object.get_values("Fusion", "bowtie1_index")

    def Grouping(self):
        if int(self.mark)/2==1:
            return "Experimental-Group"
        else:
            return "Control-Group"

    def Duplicate_detection(self):
        for dir in self.bam_dir:
            if os.path.exists(dir+"/"+self.FAD.accession):
                return False
            else:
                continue
        return True

'''   
#fusion search second step
global sample_class
fusion_sec_step_num=60
fusion_sec_step_threads=24
node=free_nodes(fusion_sec_step_threads)
shell_cmd='python ~/tools/2017/fusion_search_second_step.py '+sample_class+' '+str(fusion_sec_step_num)+' '+str(fusion_sec_step_threads)
shell_cmd='echo "'+shell_cmd+'"|qsub -l walltime=600:00:00 -l nodes='+node+':ppn='+str(fusion_sec_step_threads)+' -e /dev/null -o /dev/null -N fus_sec_tep'
job_status=commands.getoutput("qstat -a |grep fus_sec_tep|awk '{print $10}'")
if job_status != 'R' and job_status != 'Q':
    job=commands.getoutput(shell_cmd)
    print("fusion search second step ...(job id="+job+":node="+node)
'''