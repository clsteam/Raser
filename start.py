#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author  : Yao

""" Raser-main

This file(script) can also be imported as a module and contains the following
functions:

    * head      - Passing command line arguments
    * parallel  - Parallel computing body
    * single    - Serial computing body
    * main      - Main
"""

# Standard library
import logging
import time
import argparse
import sys

# Third party library
import coloredlogs
import verboselogs
from multiprocessing import Pool
from multiprocessing_logging import install_mp_handler

# Private module
from core.ilog import Ini

from raser import __software__, __version__, __author__
from params.bundle import DirectoryManagement, GlobalDirectoryManagement
from params.static import lncrna_counts

from raser.settings import STRAND_SPECIFIC_USE_AUTOMATICALLY, PRIMARY_GTF_ANNOTATIONS

from road.raser_0 import PickDicts, sequence_format_standardized
from road.raser_1 import QualityControl
from road.raser_2 import Trim
from road.raser_3 import Alignment
from road.rmdup import Rmdup
from road.raser_4 import StrandSpecific
from road.raser_5 import Transcript
from road.raser_6 import GeneCount
from road.raser_7 import DifferentialExpression

from line.variation import Variation
from line.allele import Allele
from line.alter_splice import AlterSplice
from line.fusion import Fusion
from line.lncrna import Lncrna

from source.env import Env


def head():
    # arguments passed
    parser = argparse.ArgumentParser()

    # required
    parser.add_argument("--server", help="submit tasks to server compute nodes (PBS)", action="store_true")

    parser.add_argument("--level", help="logger level", choices=["spam", "debug", "verbose", "info", "notice", "warning", "success", "error", "critical"])
    parser.add_argument("--category", choices=["ve", "pl"], help="ve: vertebrate, pl: plant")
    parser.add_argument("--ini_file", help="configuration file, default {RASER_HOME}/config.ini", required=True)
    parser.add_argument("--sim", help="simplified process", choices=[0, 1], type=int)
    parser.add_argument("--test", help="for testing, only run a sample in the main process(0/1)", choices=[0, 1], type=int)

    args = parser.parse_args()

    main(server=args.server, category=args.category, ini_file=args.ini_file, test=args.test, level=args.level, sim=args.sim)


def parallel_1(logger: logging.Logger, params, ini: Ini):
    # standard params key: home_dir, sample, sequence_format
    dm = DirectoryManagement(sequence_format_standardized(logger, params, ini))

    logger.info("[{0}] parallel[1]...".format(dm.id))
    args = [logger, dm, ini]

    qualitycontrol = QualityControl(*args)
    qualitycontrol.initialize()
    qualitycontrol.move()

    trim = Trim(*args)
    trim.initialize()
    trim.move()

    qualitycontrol.move()

    alignment = Alignment(*args)
    alignment.initialize()
    alignment.move()

    # import shutil, os
    # dm.org_bam = os.path.join(dm.home_dir, "tophat_out", "accepted_hits.bam")
    # shutil.move(dm.bam_file, dm.org_bam)
    # print(dm.org_bam, flush=True)

    rmdup = Rmdup(*args)
    rmdup.initialize()
    rmdup.move()

    strandspecific = StrandSpecific(*args)
    strandspecific.initialize()
    strandspecific.move()

    if STRAND_SPECIFIC_USE_AUTOMATICALLY:
        alignment = Alignment(*args, again=True)
        alignment.initialize()
        alignment.move()

        rmdup = Rmdup(*args)
        rmdup.initialize()
        rmdup.move()

    genecount = GeneCount(*args)
    genecount.initialize()
    genecount.move()

    transcript = Transcript(*args)
    transcript.initialize(1)
    transcript.move()

    if ini.workflow.get("allele"):
        variation = Variation(*args)
        variation.initialize()
        variation.move()
        allele = Allele(*args)
        allele.initialize(1)
        allele.move()

    if ini.workflow.get("fusion"):
        fusion = Fusion(*args)
        fusion.initialize(1)
        fusion.move()


def parallel_2(logger: logging.Logger, params, ini: Ini):
    # standard params key: home_dir, sample, sequence_format
    dm = DirectoryManagement(sequence_format_standardized(logger, params, ini))

    logger.info("[{0}] parallel[2]...".format(dm.id))
    args = [logger, dm, ini]

    if ini.workflow.get("fusion"):
        fusion = Fusion(*args)
        fusion.initialize(3)
        fusion.move()

    if ini.workflow.get("altersplice"):
        altersplice = AlterSplice(*args)
        altersplice.initialize(1)
        altersplice.move()

    if ini.workflow.get("lncrna"):
        kwargs = {
            "annotation_file": ini.lncrna_gtf,
            "attribute": "transcript_id",
            "out": dm.lncrna_counts_file,
            "counts_dict": lncrna_counts
        }
        genecount = GeneCount(*args, **kwargs)
        genecount.initialize()
        genecount.move()


def single1(logger, ini: Ini):
    logger.info("single1...")

    args = [logger, GlobalDirectoryManagement(ini), ini]

    transcript = Transcript(*args)
    transcript.initialize(2)
    transcript.move()

    if ini.workflow.get("differentialexpression"):
        differentialexpression = DifferentialExpression(*args)
        differentialexpression.initialize()
        differentialexpression.move()

    if ini.workflow.get("allele"):
        allele = Allele(*args)
        allele.initialize(2)
        allele.move()

    transcript.initialize(3)
    transcript.move()

    if ini.workflow.get("fusion"):
        fusion = Fusion(*args)
        fusion.initialize(2)
        fusion.move()

    if ini.workflow.get("lncrna"):
        lncrna = Lncrna(*args)
        lncrna.initialize()
        lncrna.move()


def single2(logger, ini: Ini):
    logger.info("single2...")
    dm = GlobalDirectoryManagement(ini)

    args = [logger, dm, ini]

    if ini.workflow.get("lncrna"):
        kwargs = {
            "out": dm.nc_dir,
            "counts_dict": lncrna_counts
        }
        differentialexpression = DifferentialExpression(*args, **kwargs)
        differentialexpression.initialize()
        differentialexpression.move()

    if ini.workflow.get("altersplice"):
        altersplice = AlterSplice(*args)
        altersplice.initialize(2)
        altersplice.move()


def main(server=False, category="ve", ini_file="", level="debug", test=0, sim=False):
    test = int(test)

    # 计时开始
    tic = time.perf_counter()

    # logging
    verboselogs.install()
    coloredlogs.install(level=level.upper(), stream=sys.stdout)
    if test:
        coloredlogs.install(level="SPAM", stream=sys.stdout)
    install_mp_handler()
    logger = logging.getLogger(__name__)

    # setting
    if STRAND_SPECIFIC_USE_AUTOMATICALLY:
        logger.warning("[STRAND_SPECIFIC_USE_AUTOMATICALLY] turn on the strand specific automatic detection, it will be align genomen twice, and it will take more time. If you need to change it, please terminate and modify it in setting.py")
    if PRIMARY_GTF_ANNOTATIONS:
        logger.warning("[PRIMARY_GTF_ANNOTATIONS], GTF compatibility is better, but it is best to check if there are redundant retention parameters in params/tools.py, especially the parameters for Star to build the genome index. If you need to change it, please terminate and modify it in setting.py")

    # ini
    ini = Ini(ini_file, logger, category)
    ini.initialize()

    logger.info("#" * 50)
    logger.info("#{0: ^48}#".format("Welcome to the " + __software__))
    logger.info("#{0: ^48}#".format("Version : " + __version__))
    logger.info("#{0: ^48}#".format("Author : " + __author__))
    logger.info("#" * 50)
    # Show the input parameters(samples, the number of threads and process pools)
    logger.info("|{0:-^48}|".format(""))
    logger.info("|{0: ^48}|".format(" : ".join(("Ini file", ini_file))))
    logger.info("|{0: ^48}|".format(" : ".join(("Log level", level))))
    logger.info("|{0: ^48}|".format(" : ".join(("Category", category))))
    logger.info("|{0:-^48}|".format(""))
    for option, value in ini.group:
        logger.info("|{0: ^48}|".format("GROUP: " + option.strip().upper()))
        for single_value in value.strip().split():
            logger.info("|{0: ^48}|".format(single_value))
    if server:
        logger.info("|{0:*<48}|".format("CLUSER"))
        for key, value in ini.pbs_info.items():
            logger.info("|{0: ^48}|".format(" : ".join((key, value))))
    logger.info("|{0:*<48}|".format("MESSAGE"))
    for key, value in ini.submit_message.items():
        logger.info("|{0: ^48}|".format(" : ".join((key, value))))
    logger.info("|{0:-^48}|".format(""))

    if test:
        logger.info("Test mode: only analysis one sample!")

    # 环境检测、建立索引
    env = Env(logger, ini)
    env.initialize()

    def error_callback(x):
        """brief logging callback"""
        pipe_logger = logging.getLogger("PIPE")
        pipe_logger.critical(x)
    # complete logging stack
    debug_log_list = []

    pickdicts = PickDicts(ini).data_hash.items()

    # papal 1
    pool = Pool(processes=int(ini.pools))
    i = 0
    for sequence_format, samples in pickdicts:
        for sample in samples:
            # sequence_format-----sample--------home_dir-
            # ----'SRA'-------('/srr11.sra',)---'/srr1'--
            params = {"sample": sample[0], "home_dir": sample[1], "sequence_format": sequence_format}

            i += 1
            logger_main = logging.getLogger(sequence_format + str(i))  # SRA-11
            if test:
                parallel_1(logger_main, params, ini)
                break
            else:
                x = pool.apply_async(parallel_1, (logger_main, params, ini,), error_callback=error_callback)
                debug_log_list.append(x)

    if not test:
        pool.close()
        if level == "debug":
            for log in debug_log_list:
                log.get()
        pool.join()

    if sim:
        logger.info("The simplified process has been completed, congratulations!")
        return

    single1(logger, ini)

    # papal 2
    debug_log_list = []
    pool = Pool(processes=int(ini.pools))
    i = 0
    for sequence_format, samples in pickdicts:
        for sample in samples:
            # sequence_format-----sample--------home_dir-
            # ----'SRA'-------('/srr11.sra',)---'/srr1'--
            params = {"sample": sample[0], "home_dir": sample[1], "sequence_format": sequence_format}

            i += 1
            logger_main = logging.getLogger(sequence_format + str(i))  # SRA-11
            if test:
                parallel_2(logger_main, params, ini)
                break
            else:
                x = pool.apply_async(parallel_2, (logger_main, params, ini,), error_callback=error_callback)
                debug_log_list.append(x)
    if not test:
        pool.close()
        if level == "debug":
            for log in debug_log_list:
                log.get()
        pool.join()

    single2(logger, ini)

    # 计时结束
    toc = time.perf_counter()
    logger.info("{0:-^25}{1:-^25}".format("USED TIME", "DEAD TIME"))

    def transform_time(t):
        """t: second"""
        if t > 3600:
            return "{0} hour".format(t/3600)
        elif t > 60:
            return "{0} minute".format(t/600)
        else:
            return "{0} second".format(t)
    logger.info("{0: ^25}{1: ^25}".format(transform_time(toc - tic), time.strftime("%Y-%m-%d %X")))


if __name__ == '__main__':
    head()
