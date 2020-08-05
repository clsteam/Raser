#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Standard library
from abc import abstractmethod, ABCMeta
from logging import Logger
import os

from state2 import switch_tools, workflow_stateful

from core.exception import SettingException
from core.ilog import Ini

from raser.settings import TOOLS_SELECTED, PRIMARY_GTF_ANNOTATIONS
from params.bundle import DirectoryManagement

from .guess_encode import GuessEncode


@workflow_stateful
class AbstractWorkflow(metaclass=ABCMeta):
    def __init__(host, logger: Logger, dm: DirectoryManagement, ini: Ini, **kwargs):
        """
        :args
            logger: 多进程主进程日志句柄
            param dm: 文件管理实例
        :param kwargs: 子类工具特有参数
        """
        host.logger = logger
        host.dm = dm
        host.ini = ini
        host.args = (logger, dm, ini)
        host.kwargs = kwargs

        host.workflow_name = host.__class__.__name__.lower()
        host.tool_name = TOOLS_SELECTED.get(host.workflow_name)
        host.is_pass = False
        host.tool_suffix = ""

    @property
    def tool_instance(host):
        """返回Run类的正在运行的最新实例，例如
        self = host.tool_instance
        self.dm 和 host.dm是不同的
        self.dm 可以直接bp属性中修改，直接用于该并行所有，不会对单行造成影响
        host.dm 是全局的，对所有的流程适用,单、并行
        """
        return host.__instance__

    @property
    def tool(host):
        """返回工具类"""
        return getattr(host.__class__, host.tool_name.capitalize() + host.tool_suffix)

    def entrance(host):
        """The entrance to the tool, which can trim the parameters and adjust the process"""
        pass

    # @property
    # def tool_runs(host):
    #     """少用，一般直接在static.py加入参数，直接调用
    #     返回软件的Run类，如此可以在AbstractWorkflow层调用Run类的类属性"""
    #     return host.Run

    def initialize(host, tool_suffix=""):
        """选择工具，并生成Run实例"""

        host.tool_suffix = str(tool_suffix)
        # TOOLS_SELECTED必须包含流程，value=""表示跳过此流程, value=None抛出异常
        if not host.tool_name:
            if host.tool_name is None:
                raise SettingException("[{0}] is not included in TOOLS_SELECTED [setting.py]".format(host.workflow_name))
            else:
                host.logger.verbose("[{0}] {1} skipped".format(host.dm.id, host.workflow_name))
                host.is_pass = True
        else:
            host.logger.spam("[{0}] initialize {1}'s tool {2}".format(host.dm.id, host.workflow_name, host.tool_name))
            switch_tools(host, host.tool)
            host.__initialize__()

    def wf_bp(host) -> bool:
        """workflow bp, different from tool bp"""
        return True

    def move(host):
        """运行"""
        host.entrance()
        if not host.is_pass:
            host.logger.debug("[{0}] running {1}-{2} ...".format(host.dm.id, host.workflow_name, host.tool_name))
            if host.wf_bp():
                host.__work__()
            host.__clear__()


class AbstractRun(metaclass=ABCMeta):

    unpaired_prefix = [None, ]
    paired_prefix = [None, None]

    phred_prefix = None

    genome_prefix = None

    genome_index_prefix = None

    annotations_prefix = None
    strand_specific_prefix = None

    bed_prefix = None

    threads_prefix = None

    def __init__(host, logger: Logger, dm: DirectoryManagement, ini: Ini, **kwargs):
        """
        :args
            logger: 多进程主进程日志句柄
            param dm: 文件管理实例
            ini: 配置
        :param kwargs: 子类工具特有参数
        """
        host.logger = logger
        host.dm = dm
        host.ini = ini
        host.annotation_file = kwargs.get("annotation_file", None)

        host.args = (logger, dm, ini)
        host.kwargs = kwargs

    @staticmethod
    def is_file_empty(f, min_size=500) -> bool:
        if os.path.exists(f) and os.path.getsize(f) > min_size:
            return False
        return True

    @property
    @abstractmethod
    def bp(self) -> bool:
        """ tool bp
        True: 没有中间文件，需要run
        False: 有中间结果
        """
        return False

    @abstractmethod
    def clear(self):
        """for user
        清楚这一步骤产生的冗余数据文件"""
        pass

    @abstractmethod
    def run(self):
        """for user
        一般具体的执行步骤，标准通用tools的运行命令"""
        pass

    def before(self):
        """for user
        可在run之前加入的步骤,用于特定的tools多余步骤"""
        pass

    def after(self):
        """for user
        可在run之前加入的步骤,用于特定的tools多余步骤"""
        pass

    def __run__(self):
        """for admin
        Workflow上层调用"""
        if self.bp:
            self.__resume__()

    def __resume__(self):
        """for admin
        Workflow上层调用"""
        self.before()
        self.run()
        self.after()

    def __clear__(self):
        """for admin
        Workflow上层调用"""
        self.clear()

    @property
    def phred(self) -> str:
        """返回特定的完整输入参数，编码类型，例如 -phred33"""
        if self.dm.phred:
            return self.phred_prefix + self.dm.phred
        if self.ini.phred:
            return self.phred_prefix + self.ini.phred
        phred = GuessEncode(self.dm.sample[0]).judge()
        self.dm.phred = phred
        return self.phred_prefix + phred

    @property
    def genome_suffix(self) -> str:
        """返回特定的完整输入参数，基因组"""
        return self.ini.genome

    @property
    def genome(self) -> str:
        """返回特定的完整输入参数，基因组"""
        return " ".join((self.genome_prefix, self.genome_suffix))

    @property
    def genome_index_bowtie1(self) -> str:
        """返回特定的完整输入参数，bowtie1基因组索引"""
        return " ".join((self.genome_index_prefix, self.ini.bowtie1_index))

    @property
    def genome_index_bowtie2(self) -> str:
        """返回特定的完整输入参数，bowtie2基因组索引"""
        return " ".join((self.genome_index_prefix, self.ini.bowtie2_index))

    @property
    def genome_index_hisat2(self) -> str:
        """返回特定的完整输入参数，hisat2基因组索引"""
        return " ".join((self.genome_index_prefix, self.ini.hisat2_index))

    @property
    def genome_index_star(self) -> str:
        """返回特定的完整输入参数，star基因组索引"""
        return " ".join((self.genome_index_prefix, self.ini.star_index))

    @property
    def annotations_suffix(self) -> str:
        return self.annotation_file if self.annotation_file else self.ini.annotations

    @property
    def annotations(self) -> str:
        """返回特定的完整输入参数，基因组注释"""
        if PRIMARY_GTF_ANNOTATIONS:
            return self.annotations_gtf
        return " ".join((self.annotations_prefix, self.annotations_suffix))

    @property
    def annotations_gtf(self) -> str:
        """返回特定的完整输入参数，基因组注释gtf格式"""
        if self.annotations_suffix.endswith("gtf"):
            return " ".join((self.annotations_prefix, self.annotations_suffix))
        else:
            return " ".join((self.annotations_prefix, self.annotation_file if self.annotation_file else self.ini.annotations_gtf))

    def void_strand_specific(self, library_type):
        """返回特定的部分输入参数，链特异性"""
        if not library_type:
            return
        if library_type in ("fr-firststrand", "fr-secondstrand", "fr-unstranded"):
            return library_type
        else:
            self.logger.warning("[{0}] unkown library type: {1}".format(self.dm.id, library_type))
            return

    @property
    def strand_specific(self) -> str:
        """返回特定的完整输入参数，链特异性"""
        library_type = ""
        if self.ini.library_type:
            library_type = self.ini.library_type
        else:
            if self.dm.library_type:
                library_type = self.dm.library_type

        library_type = self.void_strand_specific(library_type)

        if library_type:
            return " ".join((self.strand_specific_prefix, library_type))
        return ""

    @property
    def bed(self) -> str:
        """返回特定的完整输入参数，BED文件"""
        return " ".join((self.bed_prefix, self.ini.bed))

    @property
    def threads(self) -> str:
        """返回特定的完整输入参数，线程数目"""
        return " ".join((self.threads_prefix, self.ini.threads))

    @property
    def threads_all(self) -> str:
        """返回特定的完整输入参数，总线程数目"""
        return " ".join((self.threads_prefix, self.ini.ppn))

    @property
    def input_reads(self) -> str:
        """返回特定的完整输入参数，reads文件"""
        reads = self.dm.clean_data
        if len(reads) == 1:
            # SE
            return " ".join(self.unpaired_prefix + reads)
        params = self.paired_prefix + reads
        params[1], params[2] = params[2], params[1]
        return " ".join(params)


class AbstractTools(object):

    __instance__ = None

    @classmethod
    def __initialize__(cls, host):
        cls.__instance__ = cls.Run(*host.args, **host.kwargs)

    @classmethod
    def __work__(cls, host):
        cls.__instance__.__run__()

    @classmethod
    def __clear__(cls, host):
        cls.__instance__.__clear__()

    class Run(AbstractRun):

        def bp(self) -> bool:
            return False

        def run(self):
            pass

        def clear(self):
            pass

