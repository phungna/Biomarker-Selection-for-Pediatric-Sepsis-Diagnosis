library("limma");
foldChange=1;
padj=0.05;
rawexprSet=read.csv("D:\\Na\\data_after_RMA\\exp.gene.mRNA.GSE26440.csv",header=TRUE,row.names=1,check.names = FALSE);
group <- read.csv("D:\\Na\\data_after_RMA\\label_GSE26440.csv",header=TRUE,row.names=1,check.names = FALSE);
group <- group[,1];
design <- model.matrix(~0+factor(group));
colnames(design)=levels(factor(group));
rownames(design)=colnames(rawexprSet);
fit <- lmFit(rawexprSet,design);
cont.matrix<-makeContrasts(paste0(unique(group),collapse = "-"),levels = design);
fit2=contrasts.fit(fit,cont.matrix);
fit2 <- eBayes(fit2);
tempOutput = topTable(fit2,coef=1,n=Inf,adjust="BH");
nrDEG = na.omit(tempOutput);
allDiff <- nrDEG;
diff=allDiff;
write.csv(diff, "limmaOut.csv")

rawexprSet$X <- NULL
rawexprSet$ID <- NULL
rawexprSet$GeneSymbol <- NULL
