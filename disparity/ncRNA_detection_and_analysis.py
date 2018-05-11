#!/usr/bin/env python
# -*- coding: UTF-8 -*-
__author__ = 'clsteam'

import numpy as np
from re import *
import sys
from bag.Cluster import *
from bag.configparser import *

class Non_coding_RNA:

	def __init__(self,FAD,mark):
		self.FAD=FAD
		self.mark=mark
		self.color=Color()
		self.config_object = configparser()
		self.annotations=self.config_object.get_values("Genome","annotations")
		self.genome=self.config_object.get_values("Genome","genome")
		self.blast_db=self.config_object.get_values("nc-RNA","blast_db")
		self.cpc_data_db=self.config_object.get_values("nc-RNA","cpc_data_db")
		self.cufflinks_gtf_genome_to_cdna_fasta=self.config_object.get_values("nc-RNA","cufflinks_gtf_genome_to_cdna_fasta")
		self.cnci=self.config_object.get_values("nc-RNA","cnci")
		self.run()

	def run(self):
		self.color.print_debug("   @"+self.mark+"CPC-1...")
		if not os.path.exists(self.FAD.ncRNA_cufflinks_dir+'/transcripts.gtf'):
			#step1:Building transcripts(ncRNA-cpc)
			shell_cmd="cufflinks --no-update-check --library-type  fr-firststrand -g "+self.annotations+" -p "+self.FAD.threads+" -o "+self.FAD.ncRNA_cufflinks_dir+" "+self.FAD.bam_file
			run(shell_cmd)

		self.color.print_debug("   @"+self.mark+"CPC-2...")
		if not os.path.exists(self.FAD.ncRNA_cuffcompare_dir+'/ncRNA.combined.gtf'):
			#step2:comparing transcripts(ncRNA-cpc)
			shell_cmd="cuffcompare -r "+self.annotations+" -s "+self.genome+" -C -o "+self.FAD.ncRNA_cuffcompare_dir+"/ncRNA "+self.FAD.ncRNA_cufflinks_dir+"/transcripts.gtf"
			run(shell_cmd)

		self.color.print_debug("   @"+self.mark+"CPC-3...")
		if not os.path.exists(self.FAD.ncRNA_cuffcompare_dir+"/ncRNA_filter.combined.gtf"):
			#step3:filter transcripts(ncRNA)
			self.fil_go(self.FAD.ncRNA_cuffcompare_dir+"/ncRNA.combined.gtf",self.FAD.ncRNA_cuffcompare_dir+"/ncRNA_filter.combined.gtf")

		self.color.print_debug("   @"+self.mark+"CPC-4...")
		if not os.path.exists(self.FAD.ncRNA_dir+"result.csv"):
			#convert gtf to fa(ncRNA)
			shell_cmd=self.cufflinks_gtf_genome_to_cdna_fasta+' '+self.FAD.ncRNA_cuffcompare_dir+"/ncRNA_filter.combined.gtf "+self.genome+" >"+self.FAD.ncRNA_dir+"/ncRNA.fasta"
			run(shell_cmd)

			#build db(ncRNA)
			shell_cmd="formatdb -i "+self.genome+" -p T -n prot_db"
			run(shell_cmd)

			#generate csv(ncRNA)
			shell_cmd="run_predict.sh "+self.FAD.ncRNA_dir+"/ncRNA.fasta "+self.FAD.ncRNA_dir+"/result.csv "+self.blast_db+" "+self.cpc_data_db
			run(shell_cmd)

		self.color.print_debug("   @"+self.mark+"CPC-5...")
		'''
		if not os.path.exists()(self.FAD.ncRNA_dir+"CNCI_out"):
			os.makedirs(self.FAD.ncRNA_dir+"CNCI_out")
		'''
		if not os.path.exists(self.FAD.ncRNA_dir+"/CNCI_out_Tmp_Dir/CNCI_detil"):
			shell_cmd="python "+self.cnci+" -f "+self.FAD.ncRNA_dir+"/ncRNA.fasta -o "+self.FAD.ncRNA_dir+"/CNCI_out -m ve -p "+self.FAD.threads
			run(shell_cmd)

	def fil_go(in_gtf,out_gtf):
		#Step1 : excert file and create empty vector
		kit = open(in_gtf,'r')
		p = kit.readlines()
		kit.close()

		dict_exon_number = {};
		exon_length = [];
		class_code = [];
		transcript_id = [];
		exon_number = [];

		#Step2 : remind record'class_code\gene_id\exon_number'
		for i in range(len(p)):
			x = p[i]
			head = split(compile('\t'), x)
			exon_length.append(int(head[4]) - int(head[3]))
			ind = head[8]
			h = split(compile('";'), ind)
			gene_idi = h[0];transcrip_idi = h[1];exon_numberi = h[2];gene_namei = h[3]
			for j in range(4,len(h)):
				if h[j][0:11]==' class_code':
					class_codei = h[j];break

			tid = split(compile('"'),transcrip_idi)[-1]
			class_code.append(split(compile('"'),class_codei)[-1])
			transcript_id.append(split(compile('"'),transcrip_idi)[-1])
			exn = int(split(compile('"'),exon_numberi)[-1])
			d={tid : exn}
			dict_exon_number.update(d)


		#Step3 : get real ncRNA in file.gtf
		kit = open(out_gtf,'w')
		for i in range(len(p)):
			x=p[i]
			if  ((class_code[i] in ['x','o','u','l','j']) and (exon_length[i] >= 200)):
				if (dict_exon_number[transcript_id[i]] > 1):
					kit.write(x)
				else:  continue
			else:
				continue
		kit.close()
