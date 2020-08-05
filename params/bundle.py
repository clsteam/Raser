#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from .static import gene_counts, lncrna_counts
from raser.settings import TOOLS_SELECTED


TRANSCRIPT_DIR = "transcript_out"
QC_DIR = "fastqc_out"
VARIATION_DIR = "variation_out"
ALTER_SPLICE_DIR = "alter_splice_out"
NCRNA_DIR = "ncRNA_out"


def creat_dir(*folder):
    for fd in folder:
        if not os.path.exists(fd):
            os.makedirs(fd)
    return folder[0]


class GlobalDirectoryManagement(object):

    def __setattr__(self, key, value):
        """
        Automatic creation.
        """
        if key.endswith("_dir") and not os.path.exists(value):
            os.makedirs(value)
        else:
            pass
        self.__dict__[key] = value

    def __init__(self, ini):
        self.ini = ini
        self.id = "global"

        self.home_dir = self.ini.root_path
        self.nc_dir = os.path.join(self.home_dir, NCRNA_DIR)


class DirectoryManagement(object):
    """
    Store the output folder for each sample.
    """
    def __setattr__(self, key, value):
        """
        Automatic creation.
        """
        if key.endswith("_dir") and not os.path.exists(value):
            os.makedirs(value)
        else:
            pass
        self.__dict__[key] = value

    def __init__(self, params):
        """
        :param params: Parameters of 'sequence_format', 'home_dir' and 'sample' are required.
        eg:
            * params["sequence_format"]
            * params['sample'] --> (file_name,) == it hasn't SRA in sample
            * params['home_dir']
        """
        self.sequence_format = params.get("sequence_format")
        self.sample = params.get("sample")
        self.home_dir = params.get("home_dir")

        # new vector
        self.id = os.path.basename(self.home_dir)
        self.library_type = ""
        # based on trimmomatic(must had run)
        self.max_read_length = 0
        self.phred = ""

        # new dict
        self.qc_dir = os.path.join(self.home_dir, QC_DIR)
        self.transcript_dir = os.path.join(self.home_dir, TRANSCRIPT_DIR)
        self.variation_dir = os.path.join(self.home_dir, VARIATION_DIR)
        self.as_dir = os.path.join(self.home_dir, ALTER_SPLICE_DIR)

        # new file
        self.clean_data = ["".join((os.path.join(self.home_dir, self.id), "_clean.fq.gz")),
                           ] if params.get("sequence_format") == "SE" else [
            "".join((os.path.join(self.home_dir, self.id), "_", str(i), "_clean.fq.gz")) for i in range(1, 3)]
        # self.unclean_data = ("".join((os.path.join(self.home_dir, self.id), "_unclean.fastq")))
        self.unclean_data = ("", ) if params.get("sequence_format") == "SE" else (
            "".join((os.path.join(self.home_dir, self.id), "_", str(i), "_unclean.fq.gz")) for i in range(1, 3))

        self.bam_file = ".".join((os.path.join(self.home_dir, self.id), "bam"))
        self.vcf = "".join((os.path.join(self.variation_dir, self.id), ".vcf.gz"))
        self.gvcf = "".join((os.path.join(self.variation_dir, self.id), ".g.vcf.gz"))

        self.gene_counts_file = os.path.join(self.home_dir, self.id + gene_counts.get("suffix"))
        self.lncrna_counts_file = os.path.join(self.home_dir, self.id + lncrna_counts.get("suffix"))

        self.tpm = os.path.join(self.home_dir, self.id + ".tpm")
        self.gtf = os.path.join(self.transcript_dir, "transcripts.gtf")
        self.e_gtf = os.path.join(self.transcript_dir, "transcripts_stb.gtf")  # see "-e" in ballgown
        self.lnc_gtf = os.path.join(self.transcript_dir, "pre_lnc.gtf")

        if TOOLS_SELECTED.get("alignment") == "tophat2":
            self.org_bam = os.path.join(self.home_dir, "tophat_out", "accepted_hits.bam")
        elif TOOLS_SELECTED.get("alignment") == "hisat2":
            self.bam_file.replace("bam", "sam")
        else:
            os.path.join(creat_dir(os.path.join(self.home_dir, "star_out/")), "Aligned.sortedByCoord.out.bam")

