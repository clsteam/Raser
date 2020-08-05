#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" 
tools parameters
"""


class ParamsTrim(object):
    threshold_per_base_content = 0.1
    minlen = "20"
    leading = "15"
    trailing = "15"
    slidingwindow = "4:15"

    se_adapters = "TruSeq2-SE.fa,TruSeq3-SE.fa"
    pe_adapters = "NexteraPE-PE.fa,TruSeq2-PE.fa,TruSeq3-PE-2.fa,TruSeq3-PE.fa"


class ParamsGeneCount(object):
    additional_params = "-t exon -g gene_id"
    feature_type = "exon"
    attribute_type = "gene_id"


class ParamsVariation(object):
    min_qual = 30  # QUAL = 20 means 1 out of 100(1/100)
    min_dp = 20

    bcftools_filter_params = "-e '%QUAL<{0} || DP<{1}'".format(min_qual, min_dp)

    gatk_call_params = "-A QualByDepth -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A StrandOddsRatio -A Coverage"
    gatk_filter_params = '--filter-expression "(QD<2.0) || (MQ < 90.0) || (FS > 60.0) || (SOR > 3.0)"'


class ParamsTranscript(object):
    """
    -e: Limits the processing of read alignments to only estimate and output the assembled transcripts matching the reference transcripts given with the -G option (requires -G, recommended for -B/-b). With this option, read bundles with no reference transcripts will be entirely skipped, which may provide a considerable speed boost when the given set of reference transcripts is limited to a set of target genes, for example.
    -b : Just like -B this option enables the output of *.ctab files for Ballgown, but these files will be created in the provided directory <path> instead of the directory specified by the -o option. Note: adding the -e option is recommended with the -B/-b options, unless novel transcripts are still wanted in the StringTie GTF output.
    """
    stringtie_params = "-e -b b_param"


class ParamsStar(object):
    """
    common_additional_params: add parameters in STAR build genome index and alignment.
        such as:
        GTF(default)
            * --sjdbGTFfeatureExon exon
            * --sjdbGTFtagExonParentTranscript transcript id
            * --sjdbGTFtagExonParentGene gene_id
        GFF(should change partly according to the situation, please review the tags in the GFF file yourself)
            * --sjdbGTFfeatureExon exon
            * --sjdbGTFtagExonParentTranscript Parent
            * --sjdbGTFtagExonParentGene gene

    index_additional_params: add parameters in STAR build genome index.
        such as:
        * --genomeSAindexNbases 14  -> min(14, log2(GenomeLength)/2 - 1).
        * --genomeChrBinNbits 15    -> min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)]).
    align_additional_params: add parameters in STAR alignment.
        such as:
        * --twopassMode Basic
        * --genomeLoad LoadAndRemove
    """
    common_additional_param = "--sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene gene"
    index_additional_params = ""
    align_additional_params = ""


class ParamsAllele(object):
    # This should be set to a value that will ensure reads are uniquely mapped.
    mapq_uniquely = 50  # 255,50
    # minimum base quality at the heterozygous SNP for a read to be used.
    baseq_min = 10

    common_params = "--id_separator {0}".format("-")
    phaser_params = "--mapq {0} --baseq {1}".format(mapq_uniquely, baseq_min)


class ParamsTophatFusion(object):
    additional_params = "--keep-fasta-order --no-coverage-search --fusion-ignore-chromosomes chrM"
