#Limma for differentially expressed genes with Benjamin Hochberg method
library(limma)
library(edgeR)
library(ggrepel)
#read phenotype data
p <- read.csv("D:\\Na\\label_GSE66099.csv",header = TRUE)

# add design matrix
design <- model.matrix(~er, data = p)

#test design matrix
head(design,2)
dim(design)
colSums(design)
table(p[,"er"])

# Construct limma pipeline

#read expression data
x <- read.csv("D:\\Na\\data_GSE66099.csv", header = TRUE)

# fit to linear model
fit <- lmFit(x, design)

# calculate the t-statistics
fit <- eBayes(fit)

# Summarize results
result <- decideTests(fit[,"er"])

#Note: up-redulate gene is gene with expression of non-survival higher than expression of these gene in suvival sample

#export 108 differential expression gene
kq <- topTable(fit,coef = NULL,number = 200, genelist = fit$genes, adjust.method = "BH", sort.by = "M", p.value=0.05, lfc=log2(1.5))
kq 
write.csv(kq,"D:\\Na\\gene_diff_data.csv")
