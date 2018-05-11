#!/usr/bin/env python
# -*- coding: UTF-8 -*-
__author__ = 'clsteam'

from bag.configparser import *
from bag.Cluster import *
from allele_specific_expression import *
import os

def return_find_exp_sh():
    sh='''
#p_value=0.05
fun_run(){
        awk -v A=$2 '{if($10<-2 && $12<A){print $3"\tdown"}else if($10>2 && $12<A){print $3"\tup"}}' cds_exp.diff >../result/$1/cds
        awk -v A=$2 '{if($10<-2 && $12<A){print $3"\tdown"}else if($10>2 && $12<A){print $3"\tup"}}' isoform_exp.diff >../result/$1/isoform
        awk -v A=$2 '{if($10<-2 && $12<A){print $3"\tdown"}else if($10>2 && $12<A){print $3"\tup"}}' gene_exp.diff >../result/$1/gene
}
dir=../p_value_${p_value}_result
mkdir -p $dir
fun_run $dir $p_value
'''
    return sh

def return_R_script_sh():
    r_script = '''
coldata$treatment = factor(x = coldata$treatment,levels = c("normal","tumor"))
print(coldata)
dds = DESeqDataSetFromMatrix(countData = countMatrix, colData = coldata, design = ~ treatment)
library("BiocParallel")
register(MulticoreParam(4))
dds <- DESeq(dds)
res <- results(dds,parallel=TRUE)
#resLFC <- lfcShrink(dds, coef=2, res=res,parallel=TRUE)
resOrdered <- res[order(res$padj),]
write.csv(as.data.frame(resOrdered),file="condition_Controlreated_results.csv")
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

rld<- rlogTransformation(dds, blind=TRUE)
vsd<-varianceStabilizingTransformation(dds, blind=TRUE)

pdf("deseq2_MAplot.pdf")
plotMA(dds,ylim=c(-2,2),main='DESeq2')
dev.off()

pdf("counts_plot.pdf")
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="treatment",returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=treatment, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
dev.off()

pdf("box_plot.pdf")
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()

pdf("dispersion_plot.pdf")
plotDispEsts(dds)
dev.off()

pdf("curve_plot.pdf")
plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
dev.off()

pdf("PCA_plot.pdf")
plotPCA(rld, intgroup=c("treatment"))
dev.off()

library('RColorBrewer')
library('gplots')
pdf("DESeq2_heatmap1.pdf")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,Rowv = FALSE, Colv = FALSE, scale='none',dendrogram='none', trace='none', margin=c(10,6))
dev.off()

pdf("DESeq2_heatmap2.pdf")
heatmap.2(assay(rld)[select,], col = hmcol,Rowv = FALSE, Colv = FALSE, scale='none',dendrogram='none', trace='none', margin=c(10, 6))
dev.off()

pdf("DESeq2_heatmap3.pdf")
heatmap.2(assay(vsd)[select,], col = hmcol,Rowv = FALSE, Colv = FALSE, scale='none',dendrogram='none', trace='none', margin=c(10, 6))
dev.off()

pdf("clustering.pdf")
distsRL <- dist(t(assay(rld)))
mat<- as.matrix(distsRL)
hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace='none',col = rev(hmcol), margin=c(13, 13))
dev.off()
'''

class Dge:

    def __init__(self,floder_list,threads):
        self.color = Color()
        self.threads=threads
        self.floder_list = floder_list
        self.gtf_list_file=os.path.dirname(self.floder_list[0])+"/gtf_list_file.txt"
        self.cuffmerge_dir=os.path.dirname(self.floder_list[0])+"/merged_asm"
        self.cuffdiff_dir=os.path.dirname(self.floder_list[0])+"/diff_out"
        self.deseq_out_dir=os.path.dirname(self.floder_list[0])+"/deseq_out"
        self.deseq_count_dir=os.path.dirname(self.floder_list[0])+"/deseq_out/others"
        self.deseq_data=os.path.dirname(self.floder_list[0])+"/deseq_out/deseq.data"
        self.config_object = configparser()
        self.annotations=self.config_object.get_values("Genome","annotations")
        self.find_diff_p_value=self.config_object.get_values("Different_gene_expression","find_diff_p_value")
        self.initialize()

    def initialize(self):
        for dir in (self.cuffdiff_dir,self.deseq_out_dir,self.deseq_count_dir):
            if not os.path.exists(dir):
                os.makedirs(dir)

    def deseq2(self):
        n_Experiment = self.converge_file(self.floder_list[0],"Experiment")
        n_Control = self.converge_file(self.floder_list[1],"Control")
        accession_list = self.generate_dataset()
        self.touch_rscript(n_Experiment,n_Control)

    def cuffdiff(self):
        bam=self.generate_cuffmerge_gtf_list_file()
        shell_cmd="cuffmerge -p "+self.threads+" -g "+self.annotations+" -o "+self.cuffmerge_dir+" "+self.gtf_list_file
        run(shell_cmd)
        shell_cmd="cuffdiff -L control,experiment -o "+self.cuffdiff_dir+" -p "+self.threads+" "+self.cuffmerge_dir+"/merged.gtf "+bam
        run(shell_cmd)
        self.extract_cuffdiff_result()

    def touch_rscript(self, n_Experiment, n_Control):
        with open(self.deseq_out_dir + "/script.R", "w") as doc:
            doc.writelines("#!/usr/bin/env Rscript\n")
            doc.writelines('setwd("' + self.deseq_out_dir + '")\n')
            doc.writelines('library("DESeq2")\n')
            doc.writelines('countMatrix = read.table("' + self.deseq_data + '",header = T,row.names = 1)\n')
            # coldata = data.frame(row.names = colnames(countMatrix),group = rep(c("gt1","gt2","gt3","gt4"),2,each = 2),treatment = rep(c("control","treated"),each = 8))
            group = 'c(rep("Experiment",' + str(n_Experiment) + '),rep("Control",' + str(n_Control) + '))'
            treatment = 'c(rep("Experiment",' + str(n_Experiment) + '),rep("Control",' + str(n_Control) + '))'
            doc.writelines('coldata=data.frame(row.names=colnames(countMatrix),group=' + group + ',treatment=' + treatment + ')')
            # design = ~ group + treatment + group:treatment
            doc.writelines(return_R_script_sh())
        os.system("chmod 744 " + self.deseq_out_dir + "/script.R")
        os.system("Rscript " + self.deseq_out_dir + "/script.R")

    def generate_dataset(self):
        accession_list = []
        file_list = []
        os.chdir(self.deseq_count_dir)
        for my_file in os.listdir(self.deseq_count_dir):
            accession_list.append(my_file)
            file_list.append(self.deseq_count_dir + "/" + my_file)

        # merge
        gene_hash = {}
        for i in range(len(file_list)):
            with open(file_list[i], "r") as doc:
                for line in doc.readlines():
                    gene_name = line.strip("\n").split("\t")[0]
                    gene_count = line.strip("\n").split("\t")[1]
                    if gene_hash.has_key(gene_name):
                        gene_hash[gene_name].append(gene_count)
                    else:
                        gene_hash[gene_name] = [gene_count]

        with open(self.deseq_data, "w") as out:
            out.writelines("ID\t" + "\t".join(accession_list) + "\n")
            for key in gene_hash:
                out.writelines(key + "\t" + "\t".join(gene_hash[key]) + "\n")
        return accession_list

    def converge_file(self,in_dir, treatment, n=1000):
        i = 0
        for fr_dir in os.listdir(in_dir):
            for sample in os.listdir(in_dir + "/" + fr_dir):
                if os.path.isdir(in_dir + "/" + fr_dir + "/" + sample):
                    txt = in_dir + "/" + fr_dir + "/" + sample + "/gene_counts.txt"
                    if not os.path.exists(txt):
                        self.color.print_warning(txt+" not existed")
                        continue
                    else:
                        if i < n:
                            os.system(
                                "awk -vOFS='\t' 'NR>2{print $1,$7}' " + txt + " >" + self.deseq_count_dir + "/" + treatment + "_" + sample)
                            i += 1
                        else:
                            break
        return i

    def generate_cuffmerge_gtf_list_file(self):
        bam=""
        with open(self.gtf_list_file,"w") as doc:
            for floder in self.floder_list:
                bam_list=[]
                for element in os.listdir(floder):
                    if os.path.isdir(floder+"/"+element) and os.path.exists(floder+"/"+element+"/DES_cufflinks_out/transcripts.gtf"):
                        doc.writelines(floder+"/"+element+"/DES_cufflinks_out/transcripts.gtf\n")
                        bam_list.append(floder+"/"+element+"/"+element+".bam")
                bam+=",".join(bam_list)+" "
        return bam

    def extract_cuffdiff_result(self):
        with open(self.cuffdiff_dir+"/find_exp.sh","w") as doc:
            doc.writelines("p_value="+self.find_diff_p_value+"\n")
            doc.writelines(return_find_exp_sh())
