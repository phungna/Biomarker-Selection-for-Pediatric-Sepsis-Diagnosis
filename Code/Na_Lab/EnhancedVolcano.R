library(airway)
library(magrittr)
library(EnhancedVolcano)
data('airway')
airway$dex %<>% relevel('untrt')
ens <- rownames(airway)

#Annotate the Ensembl gene IDs to gene symbols:

library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(airway), names(symbols))]
rownames(airway) <- symbols
keep <- !is.na(rownames(airway))
airway <- airway[keep,]

#Conduct differential expression using DESeq2 in order to create 2 sets of results:
library('DESeq2')

dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds, betaPrior=FALSE)
res <- results(dds,
               contrast = c('dex','trt','untrt'))
res <- lfcShrink(dds,
                 contrast = c('dex','trt','untrt'), res=res, type = 'normal')

#Plot the most basic volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')
#https://github.com/kevinblighe/EnhancedVolcano#quick-start

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'N061011 versus N61311',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)





