
#ANOVA

#######import the table 
library(readxl)
#mainAll_time <- read_excel("mainAll_time.xls")
mainAll_time <- read_excel("mainAll_time1.xls")

View(mainAll_time)
a<-data.frame(mainAll_time)
nrow(a)
#[1] 13755


b<-data.frame(t(a))
###
row(b)
row.names(b) <- c("time","0a","0b","2a","2b","8a","8b","16a","16b")
write.csv(b,file = "b1.csv")
##a litle change & write  b1_change
dat2<-data.frame(read.csv("b1_change.csv"))

#
#Then, the time variable is set to the categorical variable.

dat2$time<-factor(dat2$time,levels = c("0","2","8","16"))
dat2$time
#Then.do ANOVA.
#Parametric ANOVA (slow; simple for loop)
baseformula <- " ~ time"
#n Counting, initialization
n<-0
m<-0
#a<-array(data = NA,dim=length(sample))
#b. Store the result of the for loop
b<-array(data = NA,dim=length(sample))
b1<-b
b2<-b
for (i in 2:ncol(dat2)) {
  formula <- paste(colnames(dat2)[i], baseformula, sep="")
  f <-summary(aov(as.formula(formula), data=dat2))[[1]][["F value"]][1]
  p <- summary(aov(as.formula(formula), data=dat2))[[1]][["Pr(>F)"]][1]
  if(f!=0){
    n=n+1
    print(paste(formula, ": f=", f, sep=""))
    b1[n]<-paste(formula, ": f=", f, sep="")
  }
  
  write.table(b1,file="bf.txt")
  
  
  if(p!=0){
    m=m+1
    print(paste(formula, ": p=", p, sep=""))
    b2[m]<-paste(formula, ": p=", p, sep="")
  }
  
  write.table(b2,file="bp.txt")
  
}
print(n)
print(m)
  
###ANOVAP&F anovaP&F.xls
ANOVAPf<-read_excel("anovaP&F.xls")
anova2way_data<-data.frame(ANOVAPf)
##Correcting the p-values with BH (fdr) method
#library(qvalue)
param <- names(anova2way_data)
p <- anova2way_data[,3]
p_fdr <- p.adjust(p, method = "fdr")
#p_fdr<-qvalue(anova2way_data$ANOVAP)
#Storing the corrected p-value back in the ANOVA object, in a new column
colN<-length(anova2way_data)+1
anova2way_data[,colN]<-p_fdr
names(anova2way_data)[colN]<-paste0(param,"FDR")





########can  ignore#######LOG FC 
##########################
#TRY calculate FC (FC for two groups)  
FDRcal<-data.frame(mainAll_time)
a1<-FDRcal


#Pre-generate 2 all-0 vectors of the same length as the number of lines in the input file, which will be used to store the p value and the difference multiplier (log2FC)
log2_FC<-c(rep(0,nrow(a1)))
FC<-c(rep(0,nrow(a1)))


for(i in 1:nrow(a1)){
  
  if(sum(a1[i,2:3])==0&&sum(a1[i,4:9])==0){
    
    log2_FC[i]<- "NA"
    FC[i] <- "NA"
    
  }else{
    
    log2_FC[i]<-log2((mean(as.numeric(a1[i,2:3]))+0.001)/(mean(as.numeric(a1[i,4:9]))+0.001))
    FC[i]<-(mean(as.numeric(a1[i,2:3]))+0.001)/(mean(as.numeric(a1[i,4:9]))+0.001)
    
  }
  
}


# Add log2FC, p value and FDR, in 3 columns, to the end of the original file.

out<-cbind(a1,log2_FC,FC)

write.table(out,file="anova.out.xls",quote=FALSE,sep="\t",row.names=FALSE)
anova2way_data
out
outall<-cbind(out,anova2way_data$ANOVA.F,anova2way_data$ANOVA.P,anova2way_data$GeneNameFDR)
colnames(outall) <- c("GeneName","0a","0b","2a","2b","8a","8b","16a","16b","log2_FC","FC","ANOVAF","ANOVAP","FDR")
write.table(outall,file="anovaoutall.xls",quote=FALSE,sep="\t",row.names=FALSE)


#anova2way_data$GeneNamesFDR<0.02
##CHOSE
install.packages("magrittr") # package installations are only needed the first time you use it
install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%

re1 = outall %>% filter(FDR<0.1)
head(re1)
##5412

#(corrected P < 0.1 plus fold change > 1.5)
#2%FDR
##BUT the author is 3280
data.frame(re1)
re2 = re1 %>% filter(FC>1.5)
  #192
  
DESANO<-re2[,-(2:9)]

###Add change column to mark up or down
#Set the threshold value to  mean+2sd
logFC_cutoff<-with(DESANO,mean(abs(log2_FC))+2*sd(abs(log2_FC)))

k1<-(DESANO$FDR<0.05)&(DESANO$log2_FC<-logFC_cutoff)
k2<-(DESANO$FDR<0.05)&(DESANO$log2_FC>logFC_cutoff)
DESANO$change=ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DESANO$change)

ANOVA_DEG<-DESANO











