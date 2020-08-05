#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

"""

# Standard library
import os

from state2 import workflow_stateful
from core.abstract import AbstractWorkflow, AbstractTools, AbstractRun
from core.decorator import runshell
from core.r_rule import conversion_function_symbol
from source import phaser_py, phaser_ae_py
from params.tools import ParamsAllele
from core.ilog import creat_dir


R_OF_DESEQ = '''
#!/usr/bin/env Rscript
# rm(list = ls())
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!is.installed("DESeq2")) BiocManager::install("DESeq2")


library("DESeq2")


setwd("{root_path}")
countMatrix = read.table("{matrix}",header = T,row.names = 1)




group = c(rep('{E_name}', {E_num}),rep('{C_name}', {C_num}))
          

coldata=data.frame(row.names=colnames(countMatrix),group=group)
coldata$group = factor(x = coldata$group ,levels = c('{E_name}','{C_name}'))
dds = DESeqDataSetFromMatrix(countData = countMatrix, colData = coldata, design = ~ group)

dds <- DESeq(dds)
# 过滤
x <- 3
keep <- rowSums(counts(dds) >= 10) >= x  
dds <- dds[keep,]

res <- results(dds,parallel=TRUE)
summary(res)


resOrdered <- res[order(res$padj),]
write.csv(as.data.frame(resOrdered),file="{res}_diff.csv")

deseq_res <- as.data.frame(resOrdered)
deseq_res$GeneName <- rownames(deseq_res)
for (i in 1:nrow(deseq_res)) @<
  if (abs(deseq_res[i,'log2FoldChange']) >= 1) deseq_res[i,'select_change'] <- 'y' else deseq_res[i,'select_change'] <- 'n'
  if (deseq_res[i,'padj'] %in% NA | abs(deseq_res[i,'padj']) >= 0.01) deseq_res[i,'select_pvalue'] <- 'n' else deseq_res[i,'select_pvalue'] <- 'y'
  deseq_res[i,'select'] <- paste(deseq_res[i,'select_change'], deseq_res[i,'select_pvalue'], sep = '')
>@
deseq_res$select <- factor(deseq_res$select, levels = c('nn', 'ny', 'yn', 'yy'), labels = c('p >= 0.01, FC < 2', 'p < 0.01, FC < 2', 'p >= 0.01, FC >= 2', 'p < 0.01, FC >= 2'))
write.table(deseq_res[c(7, 1:10)], 'DESeq2.csv', row.names = FALSE, sep = ',', quote = FALSE)

#纵轴为显著性 p 值
volcano_plot_pvalue <- ggplot(deseq_res, aes(log2FoldChange, -log(padj, 10))) +
  geom_point(aes(color = select), alpha = 0.6, size=0.5) +
  scale_color_manual(values = c('gray30', 'green4', 'red2', 'blue2')) +
  theme(legend.position = c(0.2, 0.9), legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.5) + 
  geom_hline(yintercept = -log(0.01, 10), color = 'gray', size = 0.5) +
  labs(x = 'log2 Fold Change', y = '-log10 p-value')

#纵轴为 OTU 丰度
volcano_plot_abundance <- ggplot(deseq_res, aes(log2FoldChange, 100 * baseMean / sum(deseq_res$baseMean))) +
  geom_point(aes(color = select), alpha = 0.6, size=0.5, show.legend = FALSE) +
  scale_color_manual(values = c('gray30', 'green4', 'red2', 'blue2')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.5) + 
  labs(x = 'log2 Fold Change', y = 'Abundance (%)')

#组合并输出
library("grid")
png('volcano_plot.png', width = 3000, height = 1600, res = 300, units = 'px')
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(volcano_plot_pvalue, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(volcano_plot_abundance, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
'''


class Phaser1Run(AbstractRun):
    threads_prefix = "--threads"
    bed_prefix = ""

    @property
    def bp(self) -> bool:
        if not self.is_file_empty(os.path.join(self.ini.ase_result_dir, self.dm.id + "_ae.txt")):
            return False
        return True

    def clear(self):
        pass

    @property
    def hla_bed(self):
        if self.ini.hla_bed:
            return " ".join(("--blacklist", self.ini.hla_bed))
        return ""

    @property
    def haplo_count_bed(self):
        if self.ini.haplo_count_bed:
            return " ".join(("--haplo_count_blacklist", self.ini.haplo_count_bed))
        return ""

    @runshell
    def run(self):
        """
        software: phaser
        :return: cmd
        """
        phaser_prefix = os.path.join(self.ini.ase_result_dir, self.dm.id)
        phaser_ae_result = os.path.join(self.ini.ase_result_dir, self.dm.id + "_ae.txt")
        cmd_1 = " ".join(("python2",
                          phaser_py,
                          self.threads,
                          "--bam", self.dm.bam_file,
                          "--sample", self.dm.id,
                          "--vcf", self.dm.vcf,
                          "--paired_end", '0' if self.dm.sequence_format == 'SE' else '1',
                          ParamsAllele.common_params,
                          ParamsAllele.phaser_params,
                          self.hla_bed,
                          self.haplo_count_bed,
                          "--o", phaser_prefix
                          ))
        cmd_2 = " ".join(("python2",
                          phaser_ae_py,
                          "--haplotypic_counts", phaser_prefix + ".haplotypic_counts.txt",
                          "--features", self.bed,
                          ParamsAllele.common_params,
                          "--o", phaser_ae_result
                          ))
        return cmd_1, cmd_2


class Phaser2Run(AbstractRun):
    """single"""

    def __init__(self, *args, **kwargs):
        self.prefix = os.path.join(args[2].ase_result_dir, args[2].experiment_name + "_" + args[2].control_name)
        self.array = self.prefix + ".fc"
        self.res_array = self.prefix + ".diff"
        super(Phaser2Run, self).__init__(*args, **kwargs)

    @property
    def bp(self) -> bool:
        if os.path.exists(self.res_array):
            return False
        return True

    def clear(self):
        pass

    @runshell
    def before(self):
        one_accession = self.ini.treatment_accession[0][1][0]
        cmd = " ".join(("awk",
                        "'{print $4}'",
                        os.path.join(self.ini.ase_result_dir, one_accession + "_ae.txt"),
                        ">",
                        os.path.join(self.ini.ase_result_dir, "GENENAME")
                        ))

        return cmd

    @runshell
    def run(self):
        self.ase_diff_dir = os.path.join(self.ini.ase_result_dir, "diff")
        creat_dir(self.ase_diff_dir)
        for treatment, accession in self.ini.treatment_accession:
            absolute_sample = [os.path.join(self.ini.ase_result_dir, x + "_ae.txt") for x in accession]
            if not all(map(os.path.exists, absolute_sample)):
                self.logger.notice("[ALLELE]The sample corresponding to {0}(treatment) is not found".format(treatment))
                continue
            array = os.path.join(self.ase_diff_dir, treatment + ".counts")
            script = os.path.join(self.ase_diff_dir, treatment + ".R")
            self._merging_array(array, accession, absolute_sample)
            self._touch_rscript(script, array, treatment, len(accession))

    @runshell
    def _merging_array(self, array: str, accession: list, absolute_sample: list):
        shell_merge_cmd_prefix = "awk '{_[FNR]=(_[FNR] OFS $counts_col)}END{for (i=1; i<=FNR; i++) {sub(/^ /,ANCHOR,_[i]); print _[i]}}'".replace(
            'ANCHOR', '""')
        tmp = os.path.join(self.ini.ase_result_dir, "a_counts.tmp")
        cmd_1 = " ".join((shell_merge_cmd_prefix.replace("counts_col", "5"),
                          " ".join(absolute_sample),
                          ">", tmp
                          ))
        cmd_2 = " ".join((shell_merge_cmd_prefix.replace("counts_col", "6"),
                          " ".join(absolute_sample),
                          "|paste -d ' ' {0} {1} -".format(os.path.join(self.ini.ase_result_dir, "GENENAME"), tmp),
                          "|sed '1d'",
                          "|sort |uniq",
                          "|sed '1i {0}'".format("name " + " ".join(["a_" + x for x in accession]+["b_" + x for x in accession])),
                          ">",
                          array
                          ))
        return cmd_1, cmd_2

    @runshell
    def _touch_rscript(self, script, matrix, group_name, group_num):
        """
        script: 生成的R脚本
        """
        params_dict = {
            "root_path": self.ini.ase_result_dir,
            "matrix": matrix,
            "E_name": "a_"+group_name,
            "C_name": "b_"+group_name,
            "E_num": group_num,
            "C_num": group_num,
            "res": group_name
        }
        with open(script, "w") as doc:
            doc.write(conversion_function_symbol(R_OF_DESEQ.format(**params_dict)))
        return " ".join(("Rscript", script)), 1


@workflow_stateful
class Allele(AbstractWorkflow):
    """alter splice"""

    class Phaser1(AbstractTools):
        class Run(Phaser1Run):
            pass

    class Phaser2(AbstractTools):
        class Run(Phaser2Run):
            pass
