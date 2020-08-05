#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author  : Yao

R_OF_LNCFINDER = '''
#!/usr/bin/env Rscript
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

if (!is.installed("LncFinder")) install.packages("LncFinder")

library(Biostrings)
items <- readDNAStringSet("{fasta}")
lnc <- NULL
for (i in 1:length(items))
@<
  res <- LncFinder::lnc_finder(as.character(items[[i]]), SS.features = FALSE, format = "DNA", frequencies.file = "human", svm.model = "human")
    if (res[1]$Pred == "NonCoding") lnc <- c(lnc, names(items[i]))
>@
write.table(lnc, 'LncFinder', row.names = FALSE, col.names=FALSE, quote = FALSE)
'''