#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author  : Yao

R_OF_VENN = '''
#!/usr/bin/env Rscript
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
if (!is.installed("DESeq2")) install.packages("VennDiagram")
library(VennDiagram)

rm(list = ls())
setwd("{root_path}")
df1 = read.csv("{input_1}", header = TRUE, row.names = 1)
df2 = read.csv("{input_2}", header = TRUE, row.names = 1)
library(dplyr)
set1 <- subset(df1, padj < {padj} & abs(log2FoldChange) > {lfc})
set1_all <- row.names(set1)
set1_up <- row.names(subset(set1, log2FoldChange > {lfc}))
set1_down <- row.names(subset(set1, log2FoldChange < -{lfc}))

set2 <- subset(df2, padj < {padj} & abs(log2FoldChange) > {lfc})
set2_all <- row.names(set2)
set2_up <- row.names(subset(set2, log2FoldChange > {lfc}))
set2_down <- row.names(subset(set2, log2FoldChange < -{lfc}))

overlap <- calculate.overlap(
  x = list(
    x = set1_up,
    y = set2_up
  )
)
overlap <- calculate.overlap(
  x = list(
    x = set1_down,
    y = set2_down
  )
)
overlap = union(overlap1$a3, overlap2$a3)
set1_uniq <- setdiff(set1_all, overlap)
set2_uniq <- setdiff(set2_all, overlap)

write.table(overlap$a3, 'overlap_up_gene.txt', sep = '\n', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(overlap$a3, 'overlap_down_gene.txt', sep = '\n', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(set1_uniq, '{name1}_uniq_gene.txt', sep = '', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(set2_uniq, '{name2}_uniq_gene.txt', sep = '', quote = FALSE, row.names = FALSE, col.names = FALSE)


# sample two-set Venn Diagram
venn_2 <- function(x, y, name) @<
  venn.plot <- venn.diagram(
    x = list(
      {name1} = x,
      {name2} = y
    ),
    filename = name,
    lwd = 0,
    fill = c("cornflowerblue", "darkorchid1"),
    alpha = 0.75,
    label.col = "black",
    cex = 2,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("cornflowerblue", "darkorchid1"),
    cat.cex = 2,
    cat.fontfamily = "serif",
    cat.fontface = "bold",
    cat.dist = c(0.03, 0.03),
    cat.pos = c(-20, 14)
  )
>@

# venn_2 <- function(x, y, name) @<
#   venn.plot <- venn.diagram(
#     x = list(
#       {name1} = x,
#       {name2} = y
#     ),
#     filename = name,
#     cex = 2.5,
#     cat.cex = 2.5,
#     cat.pos = c(-20, 20),
#     ext.line.lty = "dotted",
#     ext.line.lwd = 2,
#     ext.pos = 12,
#     ext.dist = -0.12,
#     ext.length = 0.85
#   )
# >@



venn_2(set1_all, set2_all, "Venn_all.jpeg")
venn_2(set1_up, set2_up, "Venn_up.jpeg")
venn_2(set1_down, set2_down, "Venn_down.jpeg")


# sample four-set Venn Diagram
venn_4 <- function(up1,down1,up2,down2) @<
  venn.plot <- venn.diagram(
    x = list(
      {name1}_up = up1,
      {name1}_down = down1,
      {name2}_up = up2,
      {name2}_down = down2
    ),
    filename = "Venn_4set.jpeg",
    col = "transparent",
    fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
    alpha = 0.50,
    label.col = c("orange", "white", "darkorchid4", "white", 
                  "white", "white", "white", "white", "darkblue", "white", 
                  "white", "white", "white", "darkgreen", "white"),
    cex = 1.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
    cat.cex = 1.5,
    cat.pos = 0,
    cat.dist = 0.07,
    cat.fontfamily = "serif",
    rotation.degree = 270,
    margin = 0.2
  )
>@
venn_4(set1_up, set1_down, set2_up, set2_down)


'''
