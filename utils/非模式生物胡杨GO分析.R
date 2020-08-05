rm(list = ls())
setwd("C:/Users/Xuyao/Desktop/胡杨/data/treatment_0h_24h/SRP028830/allele")
library(AnnotationHub)
hub <- AnnotationHub::AnnotationHub()
# 寻找胡杨属
query(hub, "Populus euphratica")
# download
euphratica.OrgDb <- hub[["AH75896"]]

columns(euphratica.OrgDb)


genename <- read.csv("all_gene.txt", header = FALSE)
genename <- as.vector(genename$V1)
genename <- strsplit(genename,'gene-LOC')
gene <- c()
for (i in 1:length(genename)) {
  gene <- union(gene, genename[[i]][2])
}

library(clusterProfiler)
# 转换ID
keytypes(euphratica.OrgDb)
ls <- c("ACCNUM","ALIAS","ENTREZID","EVIDENCEALL","GENENAME","GID","GOALL","ONTOLOGY","ONTOLOGYALL","PMID","REFSEQ","SYMBOL")
for (x in ls) {
  print(x)
  print(head(keys(euphratica.OrgDb, x)))
}




# "GENENAME","GID", "SYMBOL"
df <- bitr(gene, fromType = "ENTREZID",
                toType = c("GENENAME","GID", "SYMBOL"),
                OrgDb = euphratica.OrgDb)

"105138834" %in% gene
"105138834" %in% keys(euphratica.OrgDb, "ENTREZID")

write.table(df$GENENAME, 'GID.txt', sep = '\n', quote = FALSE, row.names = FALSE, col.names = FALSE)

# go
ego <- enrichGO(gene         = gene,
                 OrgDb         = euphratica.OrgDb,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego))

# 例子
# gene <- sample(keys(euphratica.OrgDb), 100)




sample_test <- enrichGO(genes, OrgDb=euphratica.OrgDb, keyType = 'ENTREZID',pAdjustMethod = "BH",ont = "CC", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
head(summary(sample_test))








