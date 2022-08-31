if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

library("BiocManager")

library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
BiocManager::install("WGCNA")
BiocManager::install("impute",force =TRUE)
BiocManager::install("IRanges",force =TRUE)
BiocManager::install("S4Vectors",force =TRUE)
BiocManager::install("stats",force =TRUE)

library("WGCNA")
BiocManager::install("GSEABase")
BiocManager::install("XML",force =TRUE)
library(GSEABase)
BiocManager::install("GSVA",force = TRUE)
library(GSVA)
library(readxl)

###### mainAll_time.xls
rt=read_excel("mainAll_time.xls")
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp1=rt[,2:ncol(rt)]
exp=exp1[,2:ncol(exp1)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)



#Determine if the original data went to the log
max(rt)
if(max(rt)>30) rt=log2(rt+1)     #If the maximum value of rt is greater than 30, take log

#Use normalizeBetweenArrays for correction, and assign the value of rt1 after correction
rt1=normalizeBetweenArrays(as.matrix(rt))

#Not standardized
cols=rainbow(ncol(rt)) 
par(cex = 0.7)
if(ncol(rt)>40) par(cex = 0.5)  

boxplot(rt,las=2,col =cols ) 
#dev.off()

# standardized
cols=rainbow(ncol(rt1)) 
par(cex = 0.5)
if(ncol(rt1)>40) par(cex = 0.5)  
pdf(file = "nor.pdf",width=5,height = 4.5)
boxplot(rt1,las=2,col =cols )
dev.off()

#save 
rt2=rbind(ID=colnames(rt1),rt1)
write.table(rt2,file="norexp.txt",sep="\t",quote=F,col.names = F)





####DE 
data=rt1

h0Data=data[,as.vector(colnames(data)[1:2])]
h2Data=data[,as.vector(colnames(data)[3:4])]
h8Data=data[,as.vector(colnames(data)[5:6])]
h16Data=data[,as.vector(colnames(data)[7:8])]
rt=cbind(h0Data,h2Data,h8Data,h16Data)

h0Num=ncol(h0Data)
h2Num=ncol(h2Data)
h8Num=ncol(h8Data)
h16Num=ncol(h16Data)

#limma all 
Type=c(rep("h0",h0Num),rep("h2",h2Num),rep("h8",h8Num),rep("h16",h16Num))
#define
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("h0","h2","h8","h16")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(h0-h2,h0-h8,h0-h16,h2-h8,h2-h16,h8-h16,levels=design)
help(make.names)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
head(fit2)
Diff=topTable(fit2,adjust ="fdr",number=length(rownames(data)))

#Save differential results for all genes
DIFFOUT=rbind(id=colnames(Diff),Diff)
write.table(DIFFOUT,file="DIFF_all.xls",sep="\t",quote=F,col.names=F)
diffSig=Diff[with(Diff, (adj.P.Val < 0.1 )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="DIFF_af.xls",sep="\t",quote=F,col.names=F)



#Heat map showing the top 30 most diverse genes
Diff=Diff[order(as.numeric(as.vector(Diff$adj.P.Val))),]
diffGene=as.vector(rownames(Diff))
diffLength=length(diffGene)
afGene=c()
if(diffLength>(100)){
  afGene=diffGene[c(1:30,(diffLength-30+1):diffLength)]
}else{
  afGene=diffGene
}
afExp=rt[afGene,]


#Grouping Tags
Type=c(rep("h0",h0Num),rep("h2",h2Num),rep("h8",h8Num),rep("h16",h16Num))

names(Type)=colnames(rt)
Type=as.data.frame(Type)
#Annotation colors for grouped labels
ann_colors=list(gene_class=c(h0='#CC6666',h2='#3366FF',h8='#FDDCA9',h16="#FF6600"))
pdf(file="DIFF_heatmap.pdf",height=7,width=10)
pheatmap(afExp,                                                                      #热图数据
         annotation=Type,                                                            #分组
         #color = colorRampPalette(c(pal_npg()(2)[2],"white", pal_npg()(1)))(50),     #热图颜色
         color = colorRampPalette(c(pal_npg()(4)[4],"white", pal_npg()(1)))(100),     #热图颜色
         #not sure 
    
         cluster_cols =F,                                                           #不添加列聚类树
         show_colnames = F,                                                         #展示列名
         scale="row", 
         fontsize = 10,
         fontsize_row=6,
         fontsize_col=8,
         annotation_colors=ann_colors
)
dev.off()







####Because of the multiple group limma I did not come up with logFC
#So now use two sets of calculations (0h~16h)

####DE (0h~16h)



data2=rt1
#data=rt
conData=data2[,as.vector(colnames(data2)[1:2])]
treatData=data2[,as.vector(colnames(data2)[7:8])]
rt1=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#limma
Type1=c(rep("con",conNum),rep("treat",treatNum))
design1 <- model.matrix(~0+factor(Type1))
colnames(design1) <- c("con","treat")
fits <- lmFit(rt1,design1)
cont.matrix1<-makeContrasts(treat-con,levels=design1)
fits1 <- contrasts.fit(fits, cont.matrix1)
fits1 <- eBayes(fits1)
Diff=topTable(fits1,adjust='fdr',number=length(rownames(data2)))


## save results
nrDEG = na.omit(Diff) ## 去掉数据中有NA的行或列
diffsig <- nrDEG  
write.csv(diffsig, "all.limmaOut0hVS16h.csv")


###筛选差异基因
#差异基因基因的筛选，我们一般是使用P值和LogFC筛选，常用的筛选标准P<0.05,|LogFC| > 1，这是最常规的筛选标准，如果你的数据差异较大，也可以更改P值和LogFC的大小。
## 我们使用|logFC| > 0.5，padj < 0.05（矫正后P值）
foldChange = 0.5
padj = 0.05
## 筛选出所有差异基因的结果
All_diffSig <- diffsig[(diffsig$adj.P.Val < padj & (diffsig$logFC>foldChange | diffsig$logFC < (-foldChange))),]
#---------------------
dim(All_diffSig)
#[1] 815 6
## 共有815个差异基因

#筛选上调和下调的基因

diffup <-  All_diffSig[(All_diffSig$P.Value < padj & (All_diffSig$logFC > foldChange)),]
write.csv(diffup, "diffup.csv")
#
diffdown <- All_diffSig[(All_diffSig$P.Value < padj & (All_diffSig$logFC < -foldChange)),]
write.csv(diffdown, "diffdown.csv")


#plot to have a view
## 导入R包
library(ggplot2)
library(ggrepel)
##  绘制火山图
## 进行分类别
logFC <- diffsig$logFC
deg.padj <- diffsig$P.Value
data <- data.frame(logFC = logFC, padj = deg.padj)
data$group[(data$padj > 0.05 | data$padj == "NA") | (data$logFC < foldChange) & data$logFC > -foldChange] <- "Not"
data$group[(data$padj <= 0.05 & data$logFC > 1)] <-  "Up"
data$group[(data$padj <= 0.05 & data$logFC < -1)] <- "Down"
x_lim <- max(logFC,-logFC)

# PLOT
##heatmap's
DIFFOUT=rbind(id=colnames(Diff),Diff)
write.table(DIFFOUT,file="DIFF_all_0vs16.xls",sep="\t",quote=F,col.names=F)
diffSig=Diff[with(Diff, (abs(logFC)>1 & adj.P.Val < 0.05 )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="DIFF_af_0vs16.xls",sep="\t",quote=F,col.names=F)

# heat map top 30
Diff=Diff[order(as.numeric(as.vector(Diff$logFC))),]
diffGene=as.vector(rownames(Diff))
diffLength=length(diffGene)
afGene=c()
if(diffLength>(100)){
  afGene=diffGene[c(1:30,(diffLength-30+1):diffLength)]
}else{
  afGene=diffGene
}
afExp=rt[afGene,]
#group lable
Type=c(rep("h0",conNum),rep("h16",treatNum))
names(Type)=colnames(rt)
Type=as.data.frame(Type)
anncolor=list(Type=c(h0=pal_npg()(1),h16=pal_npg()(2)[2]))

pdf(file="DIFF_heatmap.pdf",height=7,width=10)
pheatmap(afExp,                                                                      #热图数据
         annotation=Type,                                                            #分组
         color = colorRampPalette(c(pal_npg()(2)[2],"white", pal_npg()(1)))(50),     #热图颜色
         cluster_cols =F,                                                           #不添加列聚类树
         show_colnames = F,                                                         #展示列名
         scale="row", 
         fontsize = 10,
         fontsize_row=6,
         fontsize_col=8,
         annotation_colors=anncolor
)
dev.off()





################Not SURE 

####GSEA analysis
deg=Diff
logFC_t=0.5
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)

deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort <- DEG %>% 
  arrange(desc(logFC))

geneList = data_all_sort$logFC #Extract the foldchange from the largest to the smallest
names(geneList) <- data_all_sort$ENTREZID #Add the corresponding ENTREZID to the above extracted foldchange
head(geneList)

#enrich
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 10,
               maxGSSize    = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none" )
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af=as.data.frame(kk2@result)
write.table(af,file=paste0("2.","all_GSEA_0hVS16.xls"),sep="\t",quote=F,col.names=T)

#The first 5 and the last 5 of GSEA results are taken after sorting
num=5
pdf(paste0("2.","down_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)])
dev.off()
pdf(paste0("2.","up_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore),num)])
dev.off()
#After sorting, the first 5 and the last 5 are displayed together
num=5
pdf(paste0("2.","all_GSEA.pdf"),width = 10,height = 10)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))])
dev.off()

#save
gseaplot2(kk2,
          title = "name",  
          "hsa04936", 
          color="red",
          base_size = 20, 
          subplots = 1:3, 
          pvalue_table = T) 

#Ridge map, show 10, save 
library(stringr)
kk2@result$Description=gsub("HALLMARK_","",kk2@result$Description)
ridgeplot(kk2,showCategory = 20)
