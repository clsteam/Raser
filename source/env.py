#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from shutil import which
from logging import Logger

from core.decorator import runshell
from core.exception import DependencyException
from core.ilog import Ini, creat_dir
from raser.settings import TOOLS_SELECTED, WHETHER_ALIGNMENT_WITH_ANNOTATIONS, TOOL_ROOT_PATH

from source import bowtie1_build, bowtie2_build, hisat2_build, star, hista2_splice_sites, gffread, \
    starfusion_build_lib_pl, gatk, gtftogenepred, genepredtobed
from params.tools import ParamsStar
from params.limit import SpeciesLimit

# Tool executable command, it can be multiple (this is used to check tool dependencies)
# All strings must be lowercase

R_PACKAGE = """

"""


class Env(object):
    CHECK = {
        "java": "java",
        "infer_experiment.py": "python model 'RSeQC'",
        "htseq-count": "python model 'HTSeq'",
    }

    def __init__(self, logger: Logger, ini: Ini):
        self.logger = logger
        self.ini = ini

    def initialize(self):
        """main"""
        self.checking()
        self.configuring()

    @property
    def genome_prefix(self):
        """/public/hg19.fa -> /public/hg19"""
        return os.path.splitext(self.ini.genome)[0]

    def checking(self):
        """
        1>检查软件依赖
        2>软件运行检查
        """
        for cmd, mess in self.CHECK.items():
            if not which(cmd):
                raise DependencyException("'{0}' is not installed.".format(mess))

        if self.ini.workflow.get("fusion") == "tophatfusion" and self.ini.species not in SpeciesLimit.tophatfusion:
            self.logger.error("the species tophatfusion can only analyze one of {0}, so this step will be skipped.".format(tuple(SpeciesLimit.tophatfusion.keys())))
            self.ini.workflow["fusion"] = False
        if self.ini.workflow.get("lncrna") and self.ini.species not in SpeciesLimit.cpat:
            self.logger.error("the species CPAT can only analyze one of {0}, so this step doesn't include CPAT.".format(tuple(SpeciesLimit.cpat.keys())))
            # self.ini.workflow["lncrna"] = False


    def configuring(self):
        """2>配置软件依赖"""
        self.configure_genome_index()
        self.configure_genome_relate()
        self.configure_tools_relate()

    def configure_genome_index(self):
        """GENOME INDEX"""
        align_tool = TOOLS_SELECTED.get("alignment").replace("tophat2", "bowtie2")  # tophat2, hisat2, star
        if not getattr(self.ini, align_tool + "_index"):
            self.logger.warning("No {0} index found.".format(align_tool))
            self.logger.info("building genome {0} index...".format(align_tool))
            getattr(self, "build_" + align_tool + "_index")()

        if self.ini.workflow.get("fusion") == "tophatfusion" and not self.ini.bowtie1_index:
            self.build_bowtie1_index()

        if self.ini.workflow.get("allele") and TOOLS_SELECTED.get("variation") == "gatk" and not os.path.exists(
                self.genome_prefix + ".dict"):
            self.build_genome_dict()

        self.ini.save()

    def configure_genome_relate(self):
        """BED, GTF"""
        if not self.ini.annotations.endswith("gtf") and not self.ini.annotations_gtf:
            self.logger.info("convert GFF to GTF...")
            self.convert_gff_to_gtf()

        if not self.ini.bed:
            self.logger.info("convert GFF|GTF BED...")
            self.convert_to_bed()

        if self.ini.workflow.get("fusion") == "starfusion" and not self.ini.starfusion_genome_resource_lib:
            self.build_starfusion_lib()

        self.ini.save()

    def configure_tools_relate(self):
        """HDRS, SPLICESITES.TXT"""
        if self.ini.workflow.get("altersplice") == "asprofile" and not self.ini.hdrs:
            self.generate_hdrs()
        if TOOLS_SELECTED.get(
                "alignment") == "hisat2" and WHETHER_ALIGNMENT_WITH_ANNOTATIONS and not self.ini.splicesites_txt:
            self.generate_splicesites_txt()

        if self.ini.workflow.get("allele"):
            if TOOLS_SELECTED.get("variation") == "gatk" and not self.ini.dbsnp:
                self.logger.warning("GATK HaplotypeCaller need dbsnp file")

            if not self.ini.hla_bed or not self.ini.haplo_count_bed:
                self.logger.error("allele phasing need hla_bed and haplo_count_bed(phaser).")

        self.ini.save()

    def generate_hdrs(self):
        """ASprofile input file"""
        hdrs = self.genome_prefix + ".hdrs"
        res = {}
        with open(self.ini.genome, "r") as handle_r, open(hdrs, "w") as handle_w:
            chrom = ""
            for line in handle_r:
                if line.startswith(">"):
                    if chrom:
                        handle_w.write("{0} /len={1} /nonNlen={2} /org={3}\n".format(chrom, res[chrom][0],
                                                                                     res[chrom][0] - res[chrom][1],
                                                                                     self.ini.species))
                    chrom = line.split()[0]
                    res.setdefault(chrom, [0, 0])
                else:
                    res[chrom][0] += len(line.strip())
                    res[chrom][1] += line.strip().count("N")
            handle_w.write(
                "{0} /len={1} /nonNlen={2} /org={3}\n".format(chrom, res[chrom][0], res[chrom][0] - res[chrom][1],
                                                              self.ini.species))
        self.ini.set("Genome", "hdrs", hdrs)

    def generate_splicesites_txt(self):
        """only support GTF"""
        splicesites_txt = self.genome_prefix + "_splicesites.txt"
        cmd = " ".join((hista2_splice_sites,
                        self.ini.annotations_gtf,
                        ">", splicesites_txt
                        ))
        self.ini.set("Genome", "hisat2_splicesites_txt", splicesites_txt)
        return cmd

    @runshell
    def convert_to_bed(self):
        bed = self.genome_prefix + ".bed"
        if self.ini.annotations.endswith("gtf"):
            annotations = self.ini.annotations
        else:
            annotations = self.ini.annotations_gtf
        cmd_1 = " ".join((gtftogenepred,
                          annotations,
                          self.genome_prefix + ".pred"
                          ))
        cmd_2 = " ".join((genepredtobed,
                          self.genome_prefix + ".pred",
                          bed
                          ))
        # cmd = """cat annotations |awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,$col_genename,$7}}' | tr -d '";' > bed""".replace("annotations", annotations).replace("col_genename", "14").replace("bed", bed)
        self.ini.set("Genome", "bed", bed)
        return cmd_1, cmd_2

    @runshell
    def convert_gff_to_gtf(self):
        """特针对不支持gff3的软件
        such as: cufflinks featurecounts
        """
        annotations_gtf = self.genome_prefix + ".gtf"
        cmd = " ".join((gffread,
                        self.ini.annotations,
                        "-T",
                        "-o",
                        annotations_gtf
                        ))
        self.ini.set("Genome", "annotations_gtf", annotations_gtf)
        return cmd

    @runshell
    def build_bowtie1_index(self):
        index = self.genome_prefix + "_bowtie1"
        self.ini.set("Genome", "bowtie1_index", index)
        return " ".join((bowtie1_build,
                         "-p", self.ini.ppn,
                         self.ini.genome,
                         index
                         ))

    @runshell
    def build_bowtie2_index(self):
        index = self.genome_prefix + "_bowtie2"
        self.ini.set("Genome", "bowtie2_index", index)
        return " ".join((bowtie2_build,
                         "--threads", self.ini.ppn,
                         self.ini.genome,
                         index
                         ))

    @runshell
    def build_hisat2_index(self):
        index = self.genome_prefix + "_hisat2"
        self.ini.set("Genome", "hisat2_index", index)
        return " ".join((hisat2_build,
                         "-p", self.ini.ppn,
                         self.ini.genome,
                         index
                         ))

    @runshell
    def build_star_index(self):
        if not self.ini.annotations.endswith("gtf"):
            self.logger.notice(
                "Star needs more parameters to set up the genome index(especially if the reference genome is in GFF format). If necessary, please add parameters in params/tools.py.")
        index = creat_dir(os.path.join(os.path.dirname(self.genome_prefix), "star_genome_index"))
        self.ini.set("Genome", "star_index", index)
        return " ".join((star,
                         "--runMode genomeGenerate",
                         "--runThreadN", self.ini.ppn,
                         "--genomeFastaFiles", self.ini.genome,
                         "--genomeDir", index,
                         "--sjdbGTFfile", self.ini.annotations,
                         ParamsStar.common_additional_param,
                         ParamsStar.index_additional_params
                         ))

    @runshell
    def build_genome_dict(self):
        return " ".join((gatk,
                         "CreateSequenceDictionary",
                         "-R", self.ini.genome
                         ))

    @runshell
    def build_starfusion_lib(self):
        extra_files_dir = os.path.join(TOOL_ROOT_PATH, "STAR-Fusion-extra-files")
        special_species_dir = creat_dir(os.path.join(extra_files_dir, self.ini.species))
        starfusion_lib = os.path.join(special_species_dir, "ctat_genome_lib_build_dir")
        kwargs = {
            "species_dir": special_species_dir,
            "build_pl": starfusion_build_lib_pl,
            "genome": self.ini.genome,
            "annotations": self.ini.annotations,
            "dfam": os.path.join(special_species_dir, "Dfam.hmm")
        }
        cmd = '''
        cd {species_dir}
        {build_pl} \
        --genome_fa {genome}\
        --gtf {annotations}\
        --dfam_db {dfam}\
        '''.format(**kwargs)

        self.ini.set("Fusion", "starfusion_genome_resource_lib", starfusion_lib)
        return cmd
