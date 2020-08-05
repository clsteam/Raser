#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author  : Yao

"""
    from core.r_rule import conversion_function_symbol

"""

SHELL_OF_CUFFDIFF = '''
#!/usr/bin/env Rscript
p_value=0.01
fun_run(){
        awk -v A=$2 '{if($10<-2 && $12<A){print $3"\tdown"}else if($10>2 && $12<A){print $3"\tup"}}' cds_exp.diff >$1/cds
        awk -v A=$2 '{if($10<-2 && $12<A){print $3"\tdown"}else if($10>2 && $12<A){print $3"\tup"}}' isoform_exp.diff >$1/isoform
        awk -v A=$2 '{if($10<-2 && $12<A){print $3"\tdown"}else if($10>2 && $12<A){print $3"\tup"}}' gene_exp.diff >$1/gene
}
dir=../p_value_${p_value}_result
mkdir -p $dir
fun_run $dir $p_value
'''

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
treatment = {treatment_vector}
          

if ('None' %in% treatment) 
@<
    coldata=data.frame(row.names=colnames(countMatrix),group=group)
    coldata$group = factor(x = coldata$group ,levels = c('{E_name}','{C_name}'))
    dds0 = DESeqDataSetFromMatrix(countData = countMatrix, colData = coldata, design = ~ group)
>@ else
@<
    coldata=data.frame(row.names=colnames(countMatrix),group=group, treatment=treatment)
    coldata$group = factor(x = coldata$group ,levels = c('{E_name}','{C_name}'))
    coldata$treatment = factor(x = coldata$treatment ,levels = c({treatment_levels}))
    dds0 = DESeqDataSetFromMatrix(countData = countMatrix, colData = coldata, design = ~ group + treatment)
>@

# 过滤
x <- 1  # 根据样本多少而定，很多的话x可以定为3
keep <- rowSums(counts(dds0) >= 10) >= x  
dds0 <- dds0[keep,]

library("BiocParallel")
register(MulticoreParam({threads}))

if(file.exists("dds.Rdata"))
@<
    load(file = "dds.Rdata")
>@ else 
@<
    dds <- DESeq(dds0)
    save(dds, file = "dds.Rdata")
>@

res <- results(dds,lfcThreshold=log2(1), parallel=TRUE)
summary(res)


resOrdered <- res[order(res$padj),]
resOrdered <- na.omit(resOrdered)
write.csv(as.data.frame(resOrdered),file="results.csv")

sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

# (res) MA plot
pdf("MAplot.pdf")
plotMA(res,ylim=c(-2,2),main='MAplot-{E_name}|{C_name}')
dev.off()


if(file.exists("resLFC.Rdata"))
@<
    load(file = "resLFC.Rdata")
>@ else 
@<
    resLFC <- lfcShrink(dds, coef=2, res=res, parallel=TRUE)
    save(resLFC, file = "resLFC.Rdata")
>@
summary(resLFC)


if(file.exists("resNorm.Rdata"))
@<
    load(file = "resNorm.Rdata")
>@ else 
@<
    resNorm <- lfcShrink(dds, coef=2, type="normal", parallel=TRUE)
    save(resNorm, file = "resNorm.Rdata")
>@


if(file.exists("resAsh.Rdata"))
@<
    load(file = "resAsh.Rdata")
>@ else 
@<
    resAsh <- lfcShrink(dds, coef=2, type="ashr", parallel=TRUE)
    save(resAsh, file = "resAsh.Rdata")
>@


pdf("muti-MAplot.pdf")
par(mfrow=c(3,1), mar=c(4,8,2,8))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
dev.off()

# (res) curve plot
pdf("curve_plot.pdf")
plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
dev.off()

# (dds) counts box dispersion plot
pdf("counts_plot.pdf")
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="group",returnData=TRUE)
if (!is.installed("ggplot2")) install.packages("ggplot2")
library("ggplot2")
ggplot(d, aes(x=group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
dev.off()

# heatmap
if (!is.installed("pheatmap")) install.packages("pheatmap")
library("pheatmap")
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
if ('None' %in% treatment) 
@<
    df <- as.data.frame(colData(dds)[c("group")])
>@ else
@<
    df <- as.data.frame(colData(dds)[c("group", "treatment")])
>@
pdf("pheatmap-ndt.pdf")
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()

dds <- estimateSizeFactors(dds0)
vsd <- vst(dds, blind=FALSE)
write.csv(assay(vsd),file = "vsd.csv")

pdf("pheatmap-vsd.pdf")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()    

# rld <- rlog(dds, blind=FALSE)
# save(rld, file = "rld.Rdata"
# pdf("pheatmap-rld.pdf")
# pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE, annotation_col=df)
# dev.off()

# cluster PCA plot
if (!is.installed("RColorBrewer")) install.packages("RColorBrewer")
library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
if ('None' %in% treatment) 
@<
    rownames(sampleDistMatrix) <- paste(vsd$group, sep="-") 
>@ else
@<
    rownames(sampleDistMatrix) <- paste(vsd$group, vsd$treatment, sep="-")
>@
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("clustering.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

pdf("PCA.pdf")
plotPCA(vsd, intgroup=c("group", "treatment"))
dev.off()
# pcaData <- plotPCA(vsd, intgroup=c("group", "treatment"), returnData=TRUE)
# ggplot(pcaData, aes(PC1, PC2, color=group, shape=treatment)) +
# geom_point(size=3) +
# xlab(paste0("PC1: ",percentVar[1],"% variance")) +
# ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
# coord_fixed()
# dev.off()

# others
pdf("box_plot.pdf")
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()

pdf("dispersion_plot.pdf")
plotDispEsts(dds)
dev.off()

pdf("MAplot-lfcThreshold.pdf")
par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()
dev.off()


deseq_res <- as.data.frame(resOrdered)
deseq_res$GeneName <- rownames(deseq_res)
for (i in 1:nrow(deseq_res)) @<
  if (abs(deseq_res[i,'log2FoldChange']) >= 1) deseq_res[i,'select_change'] <- 'y' else deseq_res[i,'select_change'] <- 'n'
  if (deseq_res[i,'padj'] %in% NA | abs(deseq_res[i,'padj']) >= 0.01) deseq_res[i,'select_pvalue'] <- 'n' else deseq_res[i,'select_pvalue'] <- 'y'
  deseq_res[i,'select'] <- paste(deseq_res[i,'select_change'], deseq_res[i,'select_pvalue'], sep = '')
>@
deseq_res$select <- factor(deseq_res$select, levels = c('nn', 'ny', 'yn', 'yy'), labels = c('p >= 0.01, FC < 2', 'p < 0.01, FC < 2', 'p >= 0.01, FC >= 2', 'p < 0.01, FC >= 2'))
write.table(deseq_res[c(7, 1:10)], 'DESeq2.csv', row.names = FALSE, sep = ',', quote = FALSE)

# library("ggplot2")


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

R_OF_BALLGOWN = '''
#!/usr/bin/env Rscript
setwd("{root_path}")
library(ballgown)

if(file.exists("bg.rda"))
@<
    load(file = "bg.rda")
>@ else 
@<
    bg = ballgown(dataDir="{extdata}", samplePattern="", meas="all")
    save(bg, file="bg.rda")
>@

# 表型
pData(bg) = data.frame(id={sample_vector}, group = c(rep('{E_name}', {E_num}),rep('{C_name}', {C_num})))

stat_results = stattest(bg, feature='transcript', meas='FPKM', covariate='group')   
write.csv(as.data.frame(stat_results),file="results.csv")

# pdf("trans.pdf")
# plotTranscripts(gene='ENSG00000223972', gown=bg, samples='SRR1706903', 
#     meas='FPKM', colorby='transcript', 
#     main='title')
# dev.off()
'''