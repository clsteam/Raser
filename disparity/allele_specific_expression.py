#!/usr/bin/env python
# -*- coding: UTF-8 -*-
__author__ = 'clsteam'

from bag.configparser import *
from bag.Cluster import *

class Ase:

    def __init__(self,FAD,library_type,mark):
        self.FAD=FAD
        self.library=library_type
        self.mark=mark
        self.config_object = configparser()
        self.genome=self.config_object.get_values("Genome","genome")
        self.annotations=self.config_object.get_values("Genome","annotations")
        self.color=Color()
        self.min_QUAL=self.config_object.get_values("Call_snp","min_QUAL")
        self.min_DP=self.config_object.get_values("Call_snp","min_DP")

    def run(self):
        self.call_snp()
        self.phaser()

    def call_snp(self):
        #call snp
        shell_1="samtools mpileup -uf "+self.genome+" "+self.FAD.bam_file+" |bcftools call -mv >"+self.FAD.snp_dir+"/"+self.FAD.accession+".vcf"
        #filter vcf
        shell_2="bcftools filter -s LowQual -e '%QUAL<"+self.min_QUAL+" || DP<"+self.min_DP+"' "+self.FAD.snp_dir+"/"+self.FAD.accession+".vcf > "+self.FAD.snp_dir+"/"+self.FAD.accession+"_filtered.vcf&&bgzip "+self.FAD.snp_dir+"/"+self.FAD.accession+"_filtered.vcf&&tabix -p vcf "+self.FAD.snp_dir+"/"+self.FAD.accession+"_filtered.vcf.gz"
        shell_cmd=shell_1+"&&"+shell_2
        run(shell_cmd)
        '''
        node=free_nodes(3)
        shell_cmd='echo "'+shell_cmd+'"|qsub -l walltime=600:00:00 -l nodes='+node+':ppn=1 -e /dev/null -o /dev/null -N snp_count'+self.mark
        job=commands.getoutput(shell_cmd)
        self.color.print_debug("@"+self.mark+"call snp,gene count: "+job+";"+node)
        '''

    def raw_gene_count(self):
        if self.FAD.format=='PE':
            mat=' -p '
        else:
            mat=' '
        if self.library=="--library-type fr-unstranded" :
            stranded=" "
        elif self.library=="--library-type fr-firststrand" :
            stranded=" -s 1 "
        else:
            stranded=" -s 2 "
        shell_cmd="featureCounts "+mat+stranded+" -T 1 -t exon -g gene_id -a "+self.annotations+" -o "+self.FAD.home_dir+"/gene_counts.txt "+self.FAD.bam_file
        run(shell_cmd)

    def	phaser(self):
        paired_end = '0' if format == 'SE' else '1'
        phaser_py=self.config_object.get_values("Phaser","phaser_py")
        phaser_gene_ae_py=self.config_object.get_values("Phaser","phaser_gene_ae_py")
        blacklist=self.config_object.get_values("Phaser","blacklist")
        haplo_count_blacklist=self.config_object.get_values("Phaser","haplo_count_blacklist")
        bed=self.config_object.get_values("Genome","bed")
        mapq=self.config_object.get_values("Phaser","mapq")
        baseq=self.config_object.get_values("Phaser","baseq")
        pass_only=self.config_object.get_values("Phaser","pass_only")
        no_gw_phase=self.config_object.get_values("Phaser","no_gw_phase")

        self.color.print_debug("   @"+self.mark+"ASE-1 ...")
        shell_cmd="python "+phaser_py+" --threads "+self.FAD.threads+" --bam "+self.FAD.bam_file+' --mapq '+mapq+' --baseq '+baseq+' --sample '+self.FAD.accession+' --vcf '+self.FAD.vcf_file+' --paired_end '+paired_end+" --pass_only " +pass_only+ " --blacklist "+blacklist+" --haplo_count_blacklist "+haplo_count_blacklist+' --o '+self.FAD.ase_dir+'/ase'
        run(shell_cmd)
        self.color.print_debug("   @"+self.mark+"ASE-2 ...")
        shell_cmd='python '+phaser_gene_ae_py+' --haplotypic_counts '+self.FAD.ase_dir+'/ase.haplotypic_counts.txt --features '+bed+' --no_gw_phase '+no_gw_phase+' --o '+self.FAD.ase_dir+'/ase_count.txt'
        run(shell_cmd)