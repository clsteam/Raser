# 转换fpkm，TPM
library(GenomicFeatures)
setwd("F:/")
txdb <- makeTxDbFromGFF("hg19.gtf",format="gtf")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})

rawcounts <- read.table("all_array", header = TRUE)
row.names(rawcounts) <- rawcounts[,1]
rawcounts <- rawcounts[,-1]

fpkm <- array(0, dim = c(length(rownames(rawcounts)), length(colnames(rawcounts))))
# rownames(rawcounts) == names(exons_gene_lens)
for (i in 1:length(rownames(rawcounts))) {
  for (j in 1:length(colnames(rawcounts))) {
    fpkm[i, j] <- rawcounts[i, j]/as.integer(exons_gene_lens[names(exons_gene_lens) == rownames(rawcounts)[i]])/sum(rawcounts[,j]) * 10^9
  }
}
colnames(fpkm) <- colnames(rawcounts)
rownames(fpkm) <- rownames(rawcounts)
write.table(fpkm, "F:/fpkm")