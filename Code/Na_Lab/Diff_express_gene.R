library(DESeq2)
library(tidyverse)
library(airway)
# Step1: Preparing count data

#read in counts data
counts_data <- read.csv("D:\\Na\\data_after_RMA\\GSE26378_afRMA.csv")
counts_data$X <- NULL
counts_data$ID <- NULL
counts_data$Gene.Symbol <- NULL
head(counts_data)

#read in sample info
colData <- read.csv("D:\\Na\\data_after_RMA\\SampleGSE26378.csv")

# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% row.names(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))

#Step 2: construct a DESeqDataSet object...

dds <- DESeqDataSetFromMatrix(counts = counts_data,
                       colData = colData,
                       design = ~ dexamethason)

dds

#pre-filtering: removing rows with low gene counts
#keeping rows that have at least 10 reads total

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

#set the factor level

dds$dexamethasone<- relevel(dds$dexamethasone, ref = "untreated")

#NOTE: collapse technical relicates

#Step 3: Run DESeq
dds <- DESeq(dds)
res <- results(dds)
#Explore Results

summary(res)
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

#contrasts

resultsNames(dds)







