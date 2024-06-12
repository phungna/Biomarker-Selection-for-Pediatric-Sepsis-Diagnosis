library(affy)
library(GEOquery)
library(tidyverse)


#get supplementary file
#getGEOSuppFiles("GSE26440")

#untar files
untar("D:\\GSE26440\\GSE26440_RAW.tar", exdir = 'data/')


#reading in .cel files
raw.data <- ReadAffy(celfile.path = "data/")

#performing RMA normalization
normalize.data <- rma(raw.data)

#get expression estimates
normalize.expr <- as.data.frame(exprs(normalize.data))

#data.frame()

#map probe IDs to gene symbols
gse <- getGEO("GSE26440", GSEMatrix = TRUE)

#fetch feature data to get ID - gene symbol mapping
feature.data <- gse$GSE26440_series_matrix.txt.gz@featureData@data


#log2 transform
normalize.expr <- log2(normalize.expr)

#write.csv(normalize.expr,"D:\\Na\\data_after_RMA\\diff_GSE26440_afRMA.csv")
#subset

feature.data <- feature.data[,c(1,11)]
feature.data
normalize.expr <- normalize.expr %>%
  rownames_to_column(var = 'Gene Symbol') %>%
  inner_join(., feature.data, by = 'Gene Symbol')

normalize.expr
write.csv(normalize.expr,"D:\\Na\\data_after_RMA\\GSE26440_afRMA.csv")
#write_excel_csv(normalize.expr,"D:\\Na\\data_after_RMA\\GSE26440_afRMA.csv")




