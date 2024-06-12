library(clusterProfiler)
library(KEGGprofile)
library(STRINGdb)

a=read.csv("D:\\Na\\LAB\\Code\\enrich.csv",header=T)

gene=as.character(a[,1])

gene=birt(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")

ego=enrichGO(gene=gene$ENTREZID,OrgDb='org.Hs.eg.db',ont="ALL",pAdjustMethod="BH",pvalueCutoff=0.5,qvalueCutoff=0.5,readable=TRUE)

GO=ego[1:10,c(1,2,3,9,7)]

GO$geneID=str_replace_all(GO$geneID,"/",",")

names(GO)=c("Category","ID","term","Genes","adj_pval")

circ <- circle_dat(GO,a)

chord <- chord_dat(data = circ, genes = a, process = GO$term)

GOChord(data=chord, title="",space = 0.01, gene.order = 'logFC', gene.space = 0.2, gene.size = 3,lfc.col=c('firebrick3','white','royalblue3'),process.label=8)

GOCircle(circ)

kk=enrichKEGG(gene=gene$ENTREZID,organism='hsa',keyType="kegg",pAdjustMethod="BH",pvalueCutoff=0.7,qvalueCutoff=0.7)

barplot(ego,showCategory=20,title="")

barplot(kk,showCategory=20,title="")
