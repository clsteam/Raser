#!/usr/bin/env python
# -*- coding: UTF-8 -*-
__author__ = 'clsteam'

import os
class Files_and_directories:

    def __init__(self,home_dir,format,threads):
        self.accession= os.path.basename(home_dir)
        self.format = format
        self.threads = threads
        # Folder
        self.home_dir = home_dir
        self.qc_dir = home_dir + '/fastqc_out'
        self.cufflinks_dir = home_dir + '/DES_cufflinks_out'
        self.snp_dir = home_dir + '/snp_analysis'
        self.ase_dir = home_dir + '/ase_analysis'
        self.alter_splice_dir = home_dir + '/alter_splice'
        self.ncRNA_dir = home_dir + '/ncRNA_identify'
        self.ncRNA_cufflinks_dir = self.ncRNA_dir + '/ncRNA_cufflinks_out'
        self.ncRNA_cuffcompare_dir = self.ncRNA_dir + '/ncRNA_cuffcompare_out'
        self.dir_list=[self.qc_dir,self.cufflinks_dir,self.snp_dir,self.ase_dir,self.ase_dir,self.alter_splice_dir,self.ncRNA_dir,self.ncRNA_cufflinks_dir,self.ncRNA_cuffcompare_dir]
        #deal
        if self.format == 'SE':
            clean_data = home_dir + '/' + self.accession + '_clean.fastq'  # such as ~/data/SRR3291477/SRR3291477_clean.fastq
        else:
            clean_data1 = home_dir + '/' + self.accession + '_1_clean.fastq'  # such as ~/data/SRR3291477/SRR3291477_1_clean.fastq
            clean_data2 = home_dir + '/' + self.accession + '_2_clean.fastq'  # such as ~/data/SRR3291477/SRR3291477_2_clean.fastq
            unclean_data1 = home_dir + '/' + self.accession + '_1_unclean.fastq'  # such as ~/data/SRR3291477/SRR3291477_1_clean.fastq
            unclean_data2 = home_dir + '/' + self.accession + '_2_unclean.fastq'  # such as ~/data/SRR3291477/SRR3291477_2_clean.fastq
            clean_data = [clean_data1, clean_data2]
            unclean_data = [unclean_data1, unclean_data2]
        # file
        self.clean_data = clean_data
        self.unclean_data = unclean_data
        self.bam_file = self.home_dir + "/" + self.accession + '.bam'
        self.original_bam_file = self.home_dir + "/tophat_out/accepted_hits.bam"
        self.vcf_file = self.snp_dir + "/" + self.accession + "_filtered.vcf.gz"
        self.adapter_file = self.qc_dir + '/adapter.fa'

    def creat(self):
        for dir in self.dir_list:
            if not os.path.exists(dir):
                os.makedirs(dir)

class Breakpoint_detection:

    def __init__(self,home_dir,format):
        self.accession = os.path.basename(home_dir)
        self.home_dir=home_dir
        self.format=format
        self.original_bam_file = self.home_dir + "/tophat_out/accepted_hits.bam"
        self.bam_file = self.home_dir + "/" + self.accession + '.bam'
        self.vcf_file = self.home_dir + '/snp_analysis/' + self.accession + "_filtered.vcf.gz"
        self.dse_gtf = self.home_dir + '/DES_cufflinks_out/transcripts.gtf'
        self.raw_gene_counts=self.home_dir+'/gene_counts.txt'
        self.ase_counts=self.home_dir + '/ase_analysis/ase_count.txt'
        self.ncRNA_CNCI_detil=home_dir + '/ncRNA_identify/CNCI_out_Tmp_Dir/CNCI_detil'
        self.Alter_splice=self.home_dir + '/alter_splice/fpkm_'+self.accession+'_stat'

        #Breakpoint
        if self.format == 'SE':
            clean_data = home_dir + '/' +self.accession + '_clean.fastq'  # such as ~/data/SRR3291477/SRR3291477_clean.fastq
            self.ftf = self.file_existed(clean_data)
            fastqc_data_txt = self.home_dir + '/fastqc_out/' + self.accession + '_fastqc/fastqc_data.txt'
            self.Fastqc_1 = self.file_existed(fastqc_data_txt)
            fastqc_data_txt=self.home_dir + '/fastqc_out/'+ self.accession + '_clean_fastqc/fastqc_data.txt'
            self.Fastqc_2=self.file_existed(fastqc_data_txt)
        else:
            clean_data1 = home_dir + '/' +self.accession + '_1_clean.fastq'  # such as ~/data/SRR3291477/SRR3291477_1_clean.fastq
            clean_data2 = home_dir + '/' +self.accession + '_2_clean.fastq'  # such as ~/data/SRR3291477/SRR3291477_2_clean.fastq
            self.ftf = self.file_existed(clean_data1) and self.file_existed(clean_data2)
            fastqc_data_txt1 = self.home_dir + '/fastqc_out/' + self.accession + '1_fastqc/fastqc_data.txt'
            fastqc_data_txt2 = self.home_dir + '/fastqc_out/' + self.accession + '2_fastqc/fastqc_data.txt'
            self.Fastqc_1 = self.file_existed(fastqc_data_txt1) or self.file_existed(fastqc_data_txt2)
            fastqc_data_txt1 = self.home_dir + '/fastqc_out/' + self.accession + '1_clean_fastqc/fastqc_data.txt'
            fastqc_data_txt2 = self.home_dir + '/fastqc_out/' + self.accession + '2_clean_fastqc/fastqc_data.txt'
            self.Fastqc_2 = self.file_existed(fastqc_data_txt1) or self.file_existed(fastqc_data_txt2)

        self.mapping=self.file_existed(self.original_bam_file)
        self.deal_bam=self.file_existed(self.bam_file)
        self.cufflinks=self.file_existed(self.dse_gtf)
        self.call_snp=self.file_existed(self.raw_gene_counts)
        self.phaser=self.file_existed(self.ase_counts)
        self.Non_coding_RNA=self.file_existed(self.ncRNA_CNCI_detil)

    def file_existed(self,file_path):
        if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
            return False
        else:
            return True


