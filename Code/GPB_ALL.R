###修稿+文章中画图的程序
###修稿的程序标了序号，正常书写
###文章中的图的程序需要加上的地方，留白的画标了******
###李心慧
###2021年8月15日
##################################
###一.数据处理
###二.Signature lncRNA


######一.数据处理
######1.1 对数据进行批次效应校正
library(dplyr)
library(BiocManager)
BiocManager::install.packages("sva")
library(sva)

setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果一/去批次效应/")
coldata=read.table("colData.txt",sep="\t",header=T,row.names=1)
sample=rownames(coldata)
data1=read.table("去除重复的总的表达谱.txt",sep="\t",header=T,row.names=1)
data2=read.table("pubmed_data.txt",sep="\t",header=T,row.names=1)
alldata=intersect(rownames(data1),rownames(data2))
res=cbind(data1[alldata,],data2[alldata,])
res1=res[,sample]
aa=apply(res1,1,mean)
res2=res1[which(aa!=0),]###删除零表达的基因
write.table(res2,"all_sample_all_RNA_exp.txt",sep="\t",quote=F,row.names=T,col.names=T)

res3=log2(res2+0.001)
write.table(res3,"all_sample_all_RNA_log2exp.txt",sep="\t",quote=F,row.names=T,col.names=T)
modcombat = model.matrix(~1,data=coldata)
batch = coldata$series
#GSE=coldata$series
#platform=coldata$platform
RES=ComBat(dat=res3, batch=batch, mod=modcombat, par.prior=TRUE,mean.only=T)###批次效应校正
#Found 56991 genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.
#Using the 'mean only' version of ComBat
#Found27batches
#Note: one batch has only one sample, setting mean.only=TRUE
#Adjusting for0covariate(s) or covariate level(s)
#Standardizing Data across genes
#Fitting L/S model and finding priors
#Finding parametric adjustments
#Adjusting the Data
write.table(RES,"all_sample_all_RNA_log2exp_ComBat.txt",sep="\t",quote=F,row.names=T,col.names=T)

###验证数据不用进行批次效应的 校正
###计算两种数据的相关性
#按样本计算
A1=read.table("all_sample_all_RNA_log2exp.txt",sep="\t",header=T,row.names=1)
A2=read.table("all_sample_all_RNA_log2exp_ComBat.txt",sep="\t",header=T,row.names=1)
res_cor=c()
for(i in 1:218){
  cc=cor(A1[,i],A2[,i])
  res_cor=rbind(res_cor,cc)
}
rownames(res_cor)=colnames(A1)
write.table(res_cor,"all_RNA_cor_ComBat_res.txt",sep="\t",quote=F)###相关系数都是在0.99以上
###也是因为数据中0值（在同一批次里面一样的数）很多，所以ComBat得到的结果与原结果相差不多
###进一步检验数据，说明数据的批次效应影响很小，所以不做批次效应的校正，用的方法是聚类，加上样本的series的信息和platform的信息
###用所有RNA(文章里面的图是用的lncRNA画的)聚类了批次效应去除前后的相关性聚类图

A=read.table("all_sample_all_RNA_log2exp.txt",sep="\t",header=T,row.names=1)
B=matrix(0,nrow=218,ncol=218)
colnames(B)=colnames(A)
rownames(B)=colnames(A)
for(i in 1:218){
for(j in 1:218){
B[i,j]=cor(A[,i],A[,j],method="spearman")
}
}
write.table(B,"cor_218sample_allRNA_spearman.txt",quote=FALSE,sep="\t",row.names = T,col.names = T)

##层次聚类
### 导入需要的包
library(pheatmap)
library(gplots)
library(vegan)
library(permute)
library(lattice)
coldata=read.table("colData.txt",sep="\t",header=T,row.names=1)

cor=as.matrix(B)
annotation_col=data.frame(  cells = factor(coldata[,1]),
                            platform= factor(coldata[,2]),series=factor(coldata[,3]))

annotation_colors=list(cells = c("B"= "#DA70D6","TCD4"= "#32CD32","TCD3"= "#00FA9A","TCD8"= "#5F9EA0",
                                "Treg"= "#90EE90","Mono"= "#FFF8DC","M0"= "#EED5D2","M1"= "#EEE0E5","M2"= "#EED2EE",
                                "NK"= "#FFFFE0","Neu"= "#6A5ACD","macro"= "#FFC0C8","DC"= "#778899","Granulocytes"="#F5FFFA"),

                       platform = c("Illumina Genome Analyzer IIx"= "#FAFAD2", "Illumina Hiseq 2000" = "#F8F8FF", 
                               "Illumina Hiseq 2500" = "#6A5ACD", "Illumina Genome Analyzer" = "#48D1CC", "AB SOLiD 4 System" = "#EECBAD", "Illumina Genome Analyzer II" = "#EEEE00",
                               "Ion Torrent Proton" = "#2E8B57", "Illumina HiScanSQ" = "#EED2EE"),

series=c("GSE26530"="#F5F5DC","GSE45734"="#FFFFE0","GSE45982"="#FFFF00","GSE55536"="#FFD700","GSE57494"="#DAA520","GSE58596"="#B8860B","GSE59846"="#F5DEB3",
         "GSE60482"="#FFA500","GSE62408"="#FFDEAD","GSE64182"="#D2B48C","GSE64655"="#D2691E","GSE64713"="#8B4513","GSE66117"="#FFA07A","GSE66385"="#FF7F50",
         "GSE66763"="#FF4500","GSE66895"="#E9967A","GSE68482"="#FF6347","GSE68795"="#FF44E1","GSE72502"="#FA8072","GSE30811"="#FFFAFA","GSE33772"="#BC8F8F",
         "GSE34260"="#CD5C5C","GSE36952"="#FF0000","GSE40131"="#A52A2A","GSE40548"="#B22222","GSE40718"="#FFF8DC","E-MTAB-2319"="80000"))

colnames(cor)=rownames(annotation_col)  

pheatmap(cor, show_colnames= FALSE, show_rownames= F, scale= "none", fontsize= 6.5,
         clustering_method ="complete",
         annotation_col= annotation_col, 
         annotation_colors= annotation_colors, 
         main ="Heatmap-sperman-All",
         col = colorRampPalette(c("dodgerblue","white","red3"))(1000),
filename ="pheatmap_218samples_all_RNA_sperman1.pdf",width = 10,height = 10)

#####2021年9月8日
#####统计校正前后的表达谱的相关性大于0.9的比例，在讨论里面可以说一下
###将校正前后的数据用Tsne展示，可以分别用细胞类型，GSE和platform信息上色
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果一/去批次效应/")
cor_res=read.table("all_RNA_cor_ComBat_res.txt",sep="\t",header=T)
ratio=length(which(cor_res[,1]>0.9))/length(cor_res[,1])
ratio
[1] 0.9633028
> length(which(cor_res[,1]>0.9))
[1] 210
> length(cor_res[,1])
[1] 218
#####这些数值可以在讨论里面说一下

library(Rtsne)
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果一/去批次效应/")
coldata=read.table("colData.txt",sep="\t",header=T,row.names=1)
A1=read.table("all_sample_all_RNA_log2exp.txt",sep="\t",header=T,row.names=1)
A2=read.table("all_sample_all_RNA_log2exp_ComBat.txt",sep="\t",header=T,row.names=1)


colors = c("#DA70D6","#00FA9A","#EED5D2","#EEE0E5","#EED2EE","#FFF8DC","#32CD32",
           "#FFFFE0","#5F9EA0","#F5FFFA","#778899","#6A5ACD","#90EE90","#FFC0C8")

names(colors) = unique(coldata$cells)
head(colors)

a11=A2
A11=as.matrix(t(a11))
# 使用tsne函数进行tSNE降维分析
#tsne_iris = tsne(A11,k=2,perplexity=50)

iris_unique <- unique(A11)
tsne_out <- Rtsne(iris_unique,pca=FALSE,dims=2,
                  perplexity=30,theta=0.0) # Run TSNE

plot(tsne_out$Y,col=colors, asp=1,pch=20,
     xlab = "tSNE_1",ylab = "tSNE_2",main = "tSNE plot")
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")
# 添加图例
legend("topright",title = "Cell type",inset = 0.001,
       legend = unique(unique(coldata$cells)),pch=16,
       col = colors,text.width=3,cex = 0.5)




***************
Fig 1A 是样本数量的pie()的结果手动合并的
***************
Fig 1C
类似上面的程序，输入是所有lncRNA计算的相关
##提出所有的lncRNA对应的表达谱
****************
Fig 1B 14种细胞的lncRNA的表达累计概率分布
**************************

###二.计算lncRNA在不同细胞类型中的方差、或者差异性P值，FC值
### 整理一下15444个lncRNA对应所有样本的表达谱，lncRNA_exp_all_sample.txt
### 整理一下mRNA对应所有样本的表达谱，mRNA_exp_all_sample.txt
###1.1 计算lncRNA的方差
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/")
RNA=read.table("all_sample_all_RNA_log2exp.txt",sep="\t",header=T,row.names=1)
lnc=read.table("lncRNA_exp_all_sample.txt",header=T,sep="\t",row.names=1)##没有进行log转换
coldata=read.table("colData.txt",sep="\t",header=T,row.names=1)
group=coldata[,1]
LNC=as.matrix(RNA[row.names(lnc),])
write.table(LNC,"lncRNA_log2exp_all_sample.txt",quote=F,sep="\t")
b=c()
res=c()
for(i in 1:15444){
a=anova(lm(LNC[i,]~group))
p=a[[5]][1]
b=cbind(as.character(rownames(LNC)[i]),p)
res=rbind(res,b)
}
write.table(res,file="anova_p_15444Lnc.txt",quote=F,sep="\t",row.names=F,col.names=F)

sig_union=read.table("union_All.txt",sep="\t",header=T,row.names=1)
sig_p=res[which(res[,1]%in%rownames(sig_union)),]
write.table(sig_p,file="anova_p_all_sig_Lnc.txt",quote=F,sep="\t",row.names=F,col.names=F)

###p的阈值取0.05
###每个特征集中去掉p>0.05的lncRNA，然后计算剩下的sig 
####计算14cell的均值表达谱
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/")
matrix1=read.table("lncRNA_exp_all_sample.txt",header=T,sep="\t",row.names=1)
coldata=read.table("colData.txt",sep="\t",header=T,row.names=1)
group=coldata[,1]
sig=unique(coldata[,1])
matrix_14cell=matrix(0,ncol=14,nrow=length(matrix1[,1]))
colnames(matrix_14cell)=sig
rownames(matrix_14cell)=rownames(matrix1)
for(i in c(1:9,11:length(sig))){
cell=matrix1[,which(coldata[,1]==sig[i])]
matrix_14cell[,i]=apply(cell,1,mean)
}

i=10
cell=matrix1[,which(coldata[,1]==sig[i])]
matrix_14cell[,i]=cell
write.table(matrix_14cell,"lncRNA_exp_14cell.txt",sep="\t",quote=F)

###对均值进行log转换，作为CIBERSORT的输入
matrix_14cell=read.table("lncRNA_exp_14cell.txt",sep="\t",header=T,row.names=1)
matrix_14cell_log=log2(matrix_14cell+0.001)
write.table(matrix_14cell_log,"lncRNA_log2exp_14cell.txt",sep="\t",quote=F)

####取出每种细胞的signature lncRNA的表达谱（14cell)
###并且对应算出其log转换的表达谱，作为cIBERSORT输入
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/sig")
for(i in 1:length(sig)){
a=read.table(paste0(sig[i],".txt"),sep="\t",header=T,row.names=1)
sig_a=matrix_14cell[rownames(a),]
write.table(sig_a,paste0(sig[i],"_1.txt"),sep="\t",quote=F)

sig_log=log2(sig_a+0.001)
write.table(sig_log,paste0(sig[i],"_log2.txt"),sep="\t",quote=F)
}

####筛选出在14中cell中的方差显著差异的sig lncRNA的log后的表达谱
var_res=read.table("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/anova_p_all_sig_Lnc.txt",sep="\t",header=T,row.names=1)
for(i in 1:length(sig)){
logexp=read.table(paste0(sig[i],"_log2.txt"),sep="\t",header=T,row.names=1)
a_p=var_res[rownames(logexp),1]
length(a_p)
sig_p=logexp[which(a_p<0.05),]
ll=length(logexp[,1])-length(sig_p[,1])
print(ll)
write.table(sig_p,paste0(sig[i],"_log2_anova_P_0.05.txt"),sep="\t",quote=F)
}
[1] 73
[1] 45
[1] 9
[1] 20
[1] 55
[1] 20
[1] 53
[1] 45
[1] 75
[1] 13
[1] 26
[1] 6
[1] 24
[1] 12

####分别将方差筛选前后的sig set作为输入，查看预测值是否因为方差筛选而提高
BiocManager::install("e1071")
library(e1071)
source('D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/CIBERSORT.R')
for(i in 1:length(sig)){
results <- CIBERSORT("lncRNA_log2exp_14cell.txt",paste0(sig[i],"_log2_anova_P_0.05.txt"), 100, TRUE)
write.table(results,file=paste0(sig[i],"_anova_P_CIBERSORT_res.txt"),sep="\t",quote=F,row.names=T,col.names=T)
results1 <- CIBERSORT("lncRNA_log2exp_14cell.txt",paste0(sig[i],"_log2.txt"), 100, TRUE)
write.table(results1,file=paste0(sig[i],"_CIBERSORT_res.txt"),sep="\t",quote=F,row.names=T,col.names=T)
}




####2021-8-20，用CIBERSORT中的方法检验sig的差异性，双侧T检验
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/")
lnc=read.table("lncRNA_log2exp_all_sample.txt",header=TRUE,sep="\t",row.names=1)##没有进行log转换
coldata=read.table("colData.txt",sep="\t",header=TRUE,row.names=1)
group=unique(coldata[,1])
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/sig/")
for( i in c(1:9,11:14)){
sig=read.table(paste0(group[i],"_log2.txt"),sep="\t",header=TRUE,row.names=1)
sig_lnc=lnc[which(rownames(lnc)%in%rownames(sig)),]
res=matrix(1,ncol=3,nrow=length(sig_lnc[,1]))
rownames(res)=rownames(sig_lnc)
colnames(res)=c("var.test.p","T.test.p","FDR")
for(j in 1:length(sig_lnc[,1])){
    xa=sig_lnc[j,which(coldata[,1]==group[i])]
    xb=sig_lnc[j,which(coldata[,1]!=group[i])]
    qi=var.test(as.numeric(xa),as.numeric(xb),alternative = "two.sided")
    qi_p=qi[3][[1]]
    res[j,1]=qi_p
    T=t.test(xa,xb,alternative="two.sided",paired=F,var.equal=F,conf.level=0.95)
    T_p=T[3][[1]]
    res[j,2]=T_p
}
fdr=p.adjust(res[,2],method="fdr",n=length(res[,1]))
res[,3]=fdr
write.table(res,file=paste0(group[i],"_var.test_T.test_Lnc_P.txt"),quote=F,sep="\t",row.names=TRUE,col.names=TRUE)
}


#####
####筛选出在14中cell中的T检验显著差异的sig lncRNA的log后的表达谱
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/")
coldata=read.table("colData.txt",sep="\t",header=TRUE,row.names=1)
group=unique(coldata[,1])
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/sig/")
for( i in c(1:9,11:14)){
pp=read.table(paste0(group[i],"_var.test_T.test_Lnc_P.txt"),sep="\t",header=TRUE,row.names=1)
logexp=read.table(paste0(group[i],"_log2.txt"),sep="\t",header=TRUE,row.names=1)
sig_p=logexp[rownames(pp)[which(pp[,3]<0.3)],]
ll=length(logexp[,1])-length(sig_p[,1])
print(ll)
write.table(sig_p,paste0(group[i],"_log2_T_FDR_0.3.txt"),sep="\t",quote=F)
}
[1] 552
[1] 184
[1] 101
[1] 198
[1] 340
[1] 133
[1] 171
[1] 39
[1] 81
[1] 252
[1] 24
[1] 69
[1] 28



####分别将方差筛选前后的sig set作为输入，查看预测值是否因为方差筛选而提高
BiocManager::install("e1071")
library(e1071)
source('D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/CIBERSORT.R')
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/")
coldata=read.table("colData.txt",sep="\t",header=TRUE,row.names=1)
group=unique(coldata[,1])
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/sig/")
for(i in c(1:9,11:14)){
results <- CIBERSORT("lncRNA_log2exp_14cell.txt",paste0(group[i],"_log2_T_FDR_0.3.txt"), 100, TRUE)
write.table(results,file=paste0(group[i],"_T_FDR_CIBERSORT_res.txt"),sep="\t",quote=F,row.names=T,col.names=T)
results1 <- CIBERSORT("lncRNA_log2exp_14cell.txt",paste0(group[i],"_log2.txt"), 100, TRUE)
write.table(results1,file=paste0(group[i],"_CIBERSORT_res11.txt"),sep="\t",quote=F,row.names=T,col.names=T)
}



####最后用的是这个####方差分析的p值校正，按0.3筛一下看看
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/")
p=read.table("anova_p_all_sig_Lnc.txt",header=T,sep="\t",row.names=1)
PP=p.adjust(p[,1],method="fdr",n=length(p[,1]))
#PP[which(PP>0.3)]
PP=as.matrix(PP)
rownames(PP)=rownames(p)
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/")
coldata=read.table("colData.txt",sep="\t",header=TRUE,row.names=1)
group=unique(coldata[,1])
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/sig/")
for( i in 1:14){
logexp=read.table(paste0(group[i],"_log2.txt"),sep="\t",header=TRUE,row.names=1)
pp_sig=PP[rownames(logexp),]
sig_p=logexp[which(pp_sig<0.3),]
ll=length(logexp[,1])-length(sig_p[,1])
print(ll)
write.table(sig_p,paste0(group[i],"_log2_annova_FDR_0.3.txt"),sep="\t",quote=F)
}
[1] 19
[1] 18
[1] 3
[1] 7
[1] 15
[1] 8
[1] 12
[1] 15
[1] 17
[1] 6
[1] 8
[1] 1
[1] 7
[1] 4

library(e1071)
source('D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/CIBERSORT.R')
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/")
coldata=read.table("colData.txt",sep="\t",header=TRUE,row.names=1)
group=unique(coldata[,1])
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果二/sig/")
for(i in 1:14){
results <- CIBERSORT("lncRNA_log2exp_14cell.txt",paste0(group[i],"_log2_annova_FDR_0.3.txt"), 100, TRUE)
write.table(results,file=paste0(group[i],"_annova_FDR_0.3_CIBERSORT_res.txt"),sep="\t",quote=F,row.names=T,col.names=T)
}



*******************************
Fig 2 包含了一些排序和计算百分比的程序
然后是循环迭代调用Cibersort
***********************
********
#Fig 3B upset
****************
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果三（更新b-e图)/Fig 3B/sig-new/")
listinput=list()
cell=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD3 T","CD4 T","CD8 T","Treg")
for(i in 1:14){
A=read.table(paste0(cell[i],"_log2_annova_FDR_0.3.txt"),sep="\t",header=T,stringsAsFactor=F)
listinput[[i]] <- A[,1]
}
library(UpSetR) 
a=fromList(listinput)
names(a)=cell
upset(a, nsets = 14,nintersects = 31,order.by = "freq",matrix.color = "IndianRed", 
main.bar.color = "Salmon",sets.bar.color = "lavender",
point.size = 3.5, line.size = 1.0, mb.ratio = c(0.7, 0.3)) 

************
#Fig 3C jaccard 
****************
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果三（更新b-e图)/Fig 3C/sig-new/")
cell=c("DC","M0","M1","M2","Macrophage","CD3 T","CD4 T","CD8 T","B","Treg","NK","Monocyte","Granulocytes","Neutrophil")

res=matrix(0,nrow=14,ncol=14)
colnames(res)=cell
rownames(res)=cell
for(i in 1:14){
 A1=read.table(paste0(cell[i],"_log2_annova_FDR_0.3.txt"),sep="\t",header=T,stringsAsFactor=F)
  for(j in 1:14){
   A2=read.table(paste0(cell[j],"_log2_annova_FDR_0.3.txt"),sep="\t",header=T,stringsAsFactor=F)
   res[i,j]=length(intersect(A1[,1],A2[,1]))/length(union(A1[,1],A2[,1]))
  }
}
library(pheatmap)
pheatmap(res, show_colnames= T, show_rownames= F, scale= "none", fontsize= 6.5,
         #clustering_method ="complete",
         cluster_rows = F,
         cluster_cols = F,
         breaks=c(seq(0,0.3,length.out=100),seq(0.3001,0.6,length.out=90),seq(0.6001,1,length.out=10)),
         col = colorRampPalette(c("navy", "white", "firebrick3"))(200),
         width = 10,height = 10)

******************************
Fig 3D 3E几个预测方法预测结果比较
********************
Fig 3D
*****
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果三（更新b-e图)/Fig 3D/")
mrna=read.table("mRNA.txt",sep="\t",header=T,row.names=1)
mrna_log=log2(mrna+0.001)
write.table(mrna_log,"mRNA_log2_14cell.txt",sep="\t",quote=F)

res=c(0.973456089427049,0.949649294982554,0.941250676814501,0.955263729452141,0.947135721071536,
      0.941225221683006,0.931854475199465,
      0.951322340439403,0.916862461987048,0.946206081066131,0.971471668313657,0.949966659615177,
      0.926348909142735,0.914278524797303)
res=as.matrix(t(res))
colnames(res)=c("B","DC","Granulocytes","M0","M1","M2","Macrophage",
                "Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg")
corrplot(res, method = "pie")


###用的是网页版预测结果：
###xcell(邮箱），有TCGA的结果
###CIBERSORT找一下硬盘里面的LM22.txt

***
Fig 3E
********
##整理7种免疫细胞的lncRNA和mRNA的表达谱
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果三（更新b-e图)/Fig 3E/修稿修改/数据")

cell=c("B_CELL_NAIVE","CD4_NAIVE","CD4_STIM","CD8_NAIVE","CD8_STIM","M2","MONOCYTES","NK","TFH","TH1","TH17","TH2","THSTAR","TREG_MEM","TREG_NAIVE")
res=c()
for(i in 1:15){
  A=read.table(paste0(cell[i],"_TPM.csv"),header=T,sep=",")
  dim(A)
  all=A[,4:length(A[1,])]
  all_mean=apply(all,1,mean)
  res=cbind(res,all_mean)
  }
  colnames(res)=cell
  rownames(res)=A[,1]
  write.table(res,"15cell-mean-exp.txt",sep="\t",quote=F)


##手动整合成7种细胞类型7cell-mean-exp.txt
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果三（更新b-e图)/Fig 3E/修稿修改")

exp=read.table("7cell-mean-exp.txt",sep="\t",header=T,row.names=1)
mm=apply(exp,1,mean)
exp0=exp[which(mm!=0),]
exp_log=log2(exp0+0.001)
write.table(exp_log,"7cell-mean-exp-all-RNA-log2.txt",sep="\t",quote=F)

lnc=read.table("all_lncRNA_ID.txt",header=F,sep="\t")
mRNA=read.table("mRNA_log2_14cell.txt",header=T,sep="\t")
exp_log_lnc=exp_log[lnc[,1],]
exp_log_m=exp_log[mRNA[,1],]

write.table(exp_log_lnc,"7cell-mean-exp-lnc-log2.txt",sep="\t",quote=F)###这种提取方法造成很多NA,手动删除一下
write.table(exp_log_m,"7cell-mean-exp-mRNA-log2.txt",sep="\t",quote=F)###这种提取方法造成很多NA,手动删除一下
###注意：mRNA和lncRNA的样本顺序不一样，但是这里没用mRNA的谱，用的时候不影响结果（其他几种方法应用的时候，但是结果不用改）


####用新的sig，预测这七个免疫细胞的预测值结果
###其他几种已知方法的结果需要重新网页版跑一下，也可以直接改图，保留原来后面几列的扇形图
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果三（更新b-e图)/Fig 3E/修稿修改/")
source('CIBERSORT.R')
cell_log=read.table("7cell-mean-exp-lnc-log2.txt",sep="\t",header=T,row.names=1)
cell_7=c("B","CD4 T","CD8 T","M2","Monocyte","NK","Treg")
for(i in 1:7){
  sig=read.table(paste0("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果三（更新b-e图)/Fig 3E/修稿修改/sig-new/",cell_7[i],"_log2_annova_FDR_0.3.txt"),sep="\t",header=T,row.names=1)
  inter=intersect(rownames(sig),rownames(cell_log))
   sig_matrix=cell_log[inter,]
   write.table(sig_matrix,paste0("DICE_sig_matrix_",cell_7[i],"_log2.txt"),sep="\t",quote=F)
  }
   for(i in 1:7){

   results <- CIBERSORT(paste0("DICE_sig_matrix_",cell_7[i],"_log2.txt"),"7cell-mean-exp-lnc-log2.txt", 100, TRUE)
   write.table(results,paste0("DICE_CIBEERSORT_result_",cell_7[i],".txt"),sep="\t",quote=F)
 }


res=c(0.713247068669213,0.820168471541082,0.931447396305705,0.770023971263164,0.630047149556392,
       0.915508428775188,0.536841587585105)
res=as.matrix(t(res))
colnames(res)=c("B","M2","Monocyte","NK","CD8 T","CD4 T","Treg")
corrplot(res, method = "pie")

***********************************************************************************************
Fig 4
*****************
Fig 4A
************
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果四(更新a-d图）/Fig 4A/修稿修改/")
###重新计算lnc和mRNA的相关性，根据之前的阈值找出共表达gene
###注意：需要更改细胞的顺序
mRNA=read.table("mRNA.txt",sep="\t",header=T,row.names=1)
mm=apply(mRNA,1,mean)
ll=which(mm>1e-6)
mRNA=log2(mRNA[ll,]+0.001)
cell=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg")
for(n in 1:14){
lncRNA=read.table(paste0(cell[n],"_log2_annova_FDR_0.3.txt"),sep="\t",header=T,row.names=1)

for(i in 1:length(mRNA[,1])){
  count=0
  for(j in 1:length(lncRNA[,1])){
    cc=cor(as.numeric(mRNA[i,]),as.numeric(lncRNA[j,]))
    if(abs(cc)>=0.6){
      count=count+1
    }
  }
  ratio=count/length(lncRNA[,1])
  if(ratio>=0.1){
    write.table(mRNA[i,],paste0(cell[n],"_new_co-expression-mRNA-exp-1.txt"),sep="\t",quote=FALSE,append=T,col.names = FALSE)
  }
}
}


####将共表达基因集合进行功能富集分析
library(TCGAbiolinks)
cell=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg")
ID=read.table("DICE_57820_RNA_annotation.txt",sep="\t",header=T,row.names=1)
ID=ID[which(ID[,3]=="protein_coding"),]
for(n in 1:14){
Gene <- read.table(paste0(cell[n],"_new_co-expression-mRNA-exp-1.txt"),header=T,sep="\t")
Genelist=ID[Gene[,1],2]
Genelist <- toupper(as.character(Genelist))
print(length(Genelist))
system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist))
write.table(ansEA$ResBP,file=paste0(cell[n],"_bp_res-1.txt"),sep="\n",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(ansEA$ResPat,file=paste0(cell[n],"_pathway_res-1.txt"),sep="\n",row.names=FALSE,col.names=FALSE,quote=FALSE)
}
9713,9757,10924,11093,10165,9960,16074,15163, 10605,9300, 9667,9459,10891,11689


###把多个细胞的FDR结果整合到一起，合成对那种，ggplot画图用
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果四(更新a-d图）/Fig 4A/修稿修改/筛免疫相关term/")
cell=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg")
BP_res=c()
P_res=c()
for(n in 1:14){
BP=read.table(paste0(cell[n],"_bp_res-1.txt"),sep=";",header=FALSE)
Pathway=read.table(paste0(cell[n],"_pathway_res-1.txt"),sep=";",header=FALSE)
BP[,2]=substring(BP[,2],7)
BP[,3]=substring(BP[,3],6,(nchar(BP[,3])-1))
BP[,4]=substring(BP[,4],11,(nchar(BP[,4])-1))
#nchar(BP[,3])
BP1=cbind(cell[n],as.character(BP[,1]),BP[,2:4])
BP_res=rbind(BP_res,BP1)
Pathway[,2]=substring(Pathway[,2],7)
Pathway[,3]=substring(Pathway[,3],6,(nchar(Pathway[,3])-1))
Pathway[,4]=substring(Pathway[,4],11,(nchar(Pathway[,4])-1))
P1=cbind(cell[n],as.character(Pathway[,1]),Pathway[,2:4])
P_res=rbind(P_res,P1)
}
colnames(BP_res)=c("Cell","Term","FDR","ng","ncommon")
colnames(P_res)=c("Cell","Term","FDR","ng","ncommon")
write.table(BP_res,"BP_res_14cell.txt",sep="\t",quote=FALSE,col.names=T,row.names=FALSE)
write.table(P_res,"Pathway_res_14cell.txt",sep="\t",quote=FALSE,col.names=T,row.names=FALSE)
###筛选出免疫相关的功能进行展示
BP_old=read.table("BP_imm_cell-old.txt",sep="\t",header=FALSE)
BP=read.table("BP_res_14cell.txt",sep="\t",header=T)
Pa_old=read.table("Pathway_imm_cell-old.txt",sep="\t",header=T)
Pa=read.table("Pathway_res_14cell.txt",sep="\t",header=T)
for(i in 1:length(BP[,1])){
  if(BP[i,2]%in%BP_old[,2]){
    write.table(BP[i,],"BP_immune_res_last_SUP.txt",sep="\t",quote=FALSE,append=T,row.names=FALSE,col.names=FALSE)
  }
}
for(i in 1:length(Pa[,1])){
  if(Pa[i,2]%in%Pa_old[,2]){
    write.table(Pa[i,],"Pathway_immune_res_last__SUP.txt",sep="\t",quote=FALSE,append=T,row.names=FALSE,col.names=FALSE)
  }
}
  

###另外，附表里面有一个富集具体结果的表，需要根据这个结果整理一下
library(ggplot2)
data=read.table("BP_immune_res_last.txt",sep="\t",header=T)
#data$Cell=factor(data$Cell,levels=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD3 T","CD4 T","CD8 T","Treg"), ordered=TRUE)
plot1=ggplot(data,aes(Cell,Term))+geom_point(aes(size=-log10(FDR),color=as.factor(Cell)))
plot1+theme(axis.text.x=element_text(angle=45,hjust=1))

Pathway=read.table("Pathway_immune_res_last.txt",sep="\t",header=T)
Pathway$Cell=factor(Pathway$Cell,levels=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD3 T","CD4 T","CD8 T","Treg"), ordered=TRUE)
plot2=ggplot(Pathway,aes(Cell,Term))+geom_point(aes(size=-log10(FDR),color=as.factor(Cell)))+theme(axis.text.x=element_text(angle=45,hjust=1))+theme(axis.text.y=element_text(size=7))






*****************
***Fig 4C****
###先计算marker和共表达基因集合的交集占marker中的比例
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果四(更新a-d图）/Fig 4C/修稿修改/")
mRNA=read.table("mRNA.txt",sep="\t",header=T,row.names=1)
mRNA=log2(mRNA+0.001)
lncRNA=read.table("lncRNA_exp_14cell.txt",sep="\t",header=T,row.names=1)
lncRNA=log2(lncRNA+0.001)
###有几个marker不分亚型的那种，需要合成并集
##M(M0,M1,M2)
Mcell=c("M0","M1","M2")
M0=read.table("M0_new_co-expression-mRNA-exp-1.txt",sep="\t",header=T,row.names=1)
M1=read.table("M1_new_co-expression-mRNA-exp-1.txt",sep="\t",header=T,row.names=1)
M2=read.table("M2_new_co-expression-mRNA-exp-1.txt",sep="\t",header=T,row.names=1)
mrna_M=union(rownames(M2),union(rownames(M0),rownames(M1)))
  mrna_M_exp=mRNA[mrna_M,]
  write.table(mrna_M_exp,"M_new_co-expression-mRNA-exp-1.txt",sep="\t",quote=FALSE)
Tcell=c("CD8 T","CD4 T","CD3 T","Treg")
CD8=read.table("CD8 T_new_co-expression-mRNA-exp-1.txt",sep="\t",header=T,row.names=1)
CD4=read.table("CD4 T_new_co-expression-mRNA-exp-1.txt",sep="\t",header=T,row.names=1)
CD3=read.table("CD3 T_new_co-expression-mRNA-exp-1.txt",sep="\t",header=T,row.names=1)
Treg=read.table("Treg_new_co-expression-mRNA-exp-1.txt",sep="\t",header=T,row.names=1)
mrna_T=union(rownames(Treg),union(rownames(CD8),union(rownames(CD4),rownames(CD3))))
mrna_T_exp=mRNA[mrna_T,]
write.table(mrna_T_exp,"T_new_co-expression-mRNA-exp-1.txt",sep="\t",quote=FALSE)

other=c("B","DC","Granulocytes","M","Macrophage","Monocyte","Neutrophil","NK","T")
ratio_res=c()
for(i in 1:9){
  marker=read.table(paste0(other[i],"_human_marker_ensg_result111.txt"),sep="\t",header=T,row.names=1)
  mm=read.table(paste0(other[i],"_new_co-expression-mRNA-exp-1.txt"),sep="\t",header=T,row.names=1)
  inte=length(intersect(rownames(marker),rownames(mm)))
  aa=c(length(rownames(marker)),length(mm[,1]),inte)
  print(aa)
  ratio=inte/length(rownames(marker))
  ratio_res=c(ratio_res,ratio)
}
c("B","DC","Granulocytes","M","Macrophage","Monocyte","Neutrophil","NK","T")
c(length(rownames(marker)),length(mm[,1]),inte)
[1]  819 9713  247
[1]  849 9757  507
[1]    34 10924    17
[1]    37 11664    26
[1]   122 16074    83
[1]  1398 15163   775
[1]    94 10605    73
[1]   84 9300   12
[1]  4379 17179  2043

> ratio_res

0.3015873 0.5971731 0.5000000 0.7027027 0.6803279 0.5543634 0.7765957 0.1428571 0.4665449

M0=read.table("M0_log2_annova_FDR_0.3.txt",sep="\t",header=T,row.names=1)
M1=read.table("M1_log2_annova_FDR_0.3.txt",sep="\t",header=T,row.names=1)
M2=read.table("M2_log2_annova_FDR_0.3.txt",sep="\t",header=T,row.names=1)
sig_M=union(rownames(M2),union(rownames(M0),rownames(M1)))
  sig_M_exp=lncRNA[sig_M,]
  write.table(sig_M_exp,"M_log2_annova_FDR_0.3.txt",sep="\t",quote=F)
Tcell=c("CD8 T","CD4 T","CD3 T","Treg")
CD8=read.table("CD8 T_log2_annova_FDR_0.3.txt",sep="\t",header=T,row.names=1)
CD4=read.table("CD4 T_log2_annova_FDR_0.3.txt",sep="\t",header=T,row.names=1)
CD3=read.table("CD3 T_log2_annova_FDR_0.3.txt",sep="\t",header=T,row.names=1)
Treg=read.table("Treg_log2_annova_FDR_0.3.txt",sep="\t",header=T,row.names=1)
sig_T=union(rownames(Treg),union(rownames(CD8),union(rownames(CD4),rownames(CD3))))
sig_T_exp=lncRNA[sig_T,]
write.table(sig_T_exp,"T_log2_annova_FDR_0.3.txt",sep="\t",quote=F)

cell=c("B","DC","Granulocytes","M","Macrophage","Monocyte","Neutrophil","NK","T")
result=c()
for(i in 1:9){
  marker=read.table(paste0(other[i],"_human_marker_ensg_result111.txt"),sep="\t",header=T,row.names=1)
  sig=read.table(paste0(cell[i],"_log2_annova_FDR_0.3.txt"),sep="\t",header=T)
  all_sig=lncRNA[setdiff(rownames(lncRNA),rownames(sig)),]##去掉特征集后的所有lncRNA
  m=length(sig[,1])
  Bi=c()
for(j in 1:500){##随机500次
print(j)
sig_new=all_sig[sample(1:length(all_sig[,1]),m),]###随机取的新的lncRNA集合
M=c()
cor=cor(t(sig_new),t(mRNA))
AA=colnames(cor)
for(n in 1:length(cor[1,])){
if(length(which(abs(cor[,n])>0.6))>m/10){
M=rbind(M,AA[n])
}
}
bizhi=length(intersect(M[,1],rownames(marker)))/length(rownames(marker))##计算比值
Bi=rbind(Bi,bizhi)
j=j+1
}
result=cbind(result,Bi)
}
colnames(result)=cell
write.table(result,"9cell_500bizhi.txt",sep="\t",row.names=FALSE,col.names=T,quote=FALSE)

suiji=read.table("9cell_500bizhi.txt",header=T,sep="\t")
boxplot(suiji,las=2,outline=FALSE,boxwex = 0.35,col="#B4EEB4",ylim=c(0,0.8))
text(1:9,ratio_res,"*",col="#FFAEB9",cex=3)




***************************************
**********Fig 4D******
  
  #######原来的程序在G:\1-1-LncRNACIBERSORT\待整理\整理好的\2019年3月26日指标聚类-富集免疫通路-sig数量\7sig_marker相关
#####计算所有lncRNA与marker的相关性，取所有lncRNA对应的均值？
  library(fgsea)
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果四(更新a-d图）/Fig 4D/修稿修改/")
RNA=read.table("all_sample_all_RNA_exp.txt",sep="\t",header=TRUE,row.names=1)
lnc=read.table("lncRNA_exp_14cell.txt",sep="\t",header=TRUE,row.names=1)
lncRNA=RNA[rownames(lnc),]
ID=read.table("DICE_57820_RNA_annotation.txt",sep="\t",header=TRUE,row.names=1)
ID_pro=ID[which(ID[,3]=="protein_coding"),]
mRNA=RNA[rownames(ID_pro),]

cell=c("B","DC","Granulocytes","M","Macrophage","Monocyte","Neutrophil","NK","T")
result=c()
coldata=read.table("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果一/去批次效应/colData.txt",sep="\t",header=TRUE,row.names=1)
coldata[which(coldata[,1]%in%c("M0","M1","M2")),1]="M"
coldata[which(coldata[,1]%in%c("TCD3","TCD4","TCD8","Treg")),1]="T"
coldata[which(coldata[,1]=="macro"),1]="Macrophage"
coldata[which(coldata[,1]=="Mono"),1]="Monocyte"
coldata[which(coldata[,1]=="Neu"),1]="Neutrophil"


for(i in c(1:2,4:9)){
  marker=read.table(paste0(cell[i],"_human_marker_ensg_result111.txt"),sep="\t",header=TRUE,row.names=1)
  sig=read.table(paste0(cell[i],"_log2_annova_FDR_0.3.txt"),sep="\t",header=TRUE)
  inter=intersect(rownames(marker),rownames(mRNA))
  marker_exp=mRNA[inter,which(coldata[,1]==cell[i])]
  mm1=apply(marker_exp,1,mean)
  marker_exp1=marker_exp[which(mm1!=0),]
  lnc_exp=lncRNA[,which(coldata[,1]==cell[i])]
  mm2=apply(lnc_exp,1,mean)
  lnc_exp1=lnc_exp[which(mm2!=0),]
  lnc_cor=cor(t(lnc_exp1),t(marker_exp1))
  #lnc_cor=cor(lncRNA[1,],marker[1,])
  m_lnc_cor=apply(lnc_cor,1,median)
  rank_score=cbind(rownames(lnc_exp1),m_lnc_cor)
  write.table(rank_score,paste0(cell[i],"_alllnc_rank_result_marker_median.txt"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
 
  }


#####函数：将结果中的最后一个list格式的转换成data.frame,以便输出
listToframe <- function(x){
  toFrame = as.data.frame(matrix(nrow=0,ncol=1),stringsAsFactors=F)
  for(n in 1:lengths(x[, "leadingEdge"])){
    toChar = paste(as.character(unlist(x[, "leadingEdge"][n])), sep="", collapse="\t")
    toFrame[n,1] = toChar
  }
  
  x = cbind(x[,-"leadingEdge"], toFrame)
  colnames(x)[dim(x)[2]] = "leadingEdge"
  return(x) 
}



#all_gene=read.table(paste0(cell[i],"_rank_result_mRNA_mean.txt"),header=FALSE,sep="\t")
#all <- all_gene[order(all_gene[,2],decreasing=T),]
#rank <- all[,2]
#names(rank) <- all[,1]
#marker_list=list(rownames(sig))

#fgseaRes2 <- fgsea(marker_list,rank,nperm=10000)#, nperm=10000,minSize=10,maxSize=100
#fgseaRes=listToframe(fgseaRes2)
#write.table(fgseaRes,file=paste0(cell[i],"_cor_marker_GSEA_res.txt"),sep="\t",quote=F,row.names=F)
#pdf(paste0(cell[i],"_cor_marker_GSEA_res.pdf"))
#plotEnrichment(marker_list[[1]],rank,gseaParam = 1, ticksSize = 0.2)
#dev.off()
for(i in c(1:2,4:9)){
all_gene=read.table(paste0(cell[i],"_alllnc_rank_result_marker_median.txt"),header=FALSE,sep="\t")
aa=which(table(all_gene[,2])>1)
for(j in 1:length(aa)){
  num=length(which(all_gene[,2]%in%names(aa)[j]))
  for(n in 1:num){
all_gene[which(all_gene[,2]%in%names(aa)[j]),][num,2]=all_gene[which(all_gene[,2]%in%names(aa)[j]),][num-1,2]+0.0000001*(num-n)

}
all <- all_gene[order(all_gene[,2],decreasing=TRUE),]
rank <- all[,2]
names(rank) <- all[,1]


cell=c("B","T","Monocyte","Macrophage","NK","Neutrophil","DC")
for(i in 1:7){
sig=read.table(paste0(cell[i],"_log2_annova_FDR_0.3.txt"),sep="\t",header=T,row.names=1)
all_gene=read.table(paste0(cell[i],"_alllnc_rank_results_marker_median_new-abs-new1.txt"),sep="\t",header=FALSE)
all <- all_gene[order(all_gene[,2],decreasing=TRUE),]
rank <- all[,2]
names(rank) <- all[,1]
sam=list(as.character(rownames(sig)))
names(sam)=paste0(cell[i],"_sig")###不命名会报错
length(rank)
  fgseaRes1 <- fgsea(sam,rank)
  print(fgseaRes1)#, nperm=10000[1:12332]
  #fgseaRes2=fgseaMultilevel(sam,rank[1:13683])
  fgseaRes=listToframe(fgseaRes1)
  write.table(fgseaRes,file=paste0("GSEA_res_",cell[i],".txt"),sep="\t",quote=FALSE,row.names=FALSE)
  
  pdf(paste0(cell[i],"__GSEA_result_abs.pdf"))
  plotEnrichment(sam[[1]], rank)
  dev.off()
}
}

##i=1
pathway  pval  padj log2err        ES     NES size
B_sig 1e-10 1e-10      NA 0.6324892 2.27324  771
T_sig 1e-10 1e-10      NA 0.6195184 1.841248  907
Monocyte_sig    1    1       0 0.2395061 0.7644501  397
Macrophage_sig 1e-10 1e-10      NA 0.6901901 2.235413  426
NK_sig 1e-10 1e-10      NA 0.3651022 1.714825  895
Neutrophil_sig 0.999001 0.999001 0.001442695 0.2085976 0.7132656  189
DC_sig 1e-10 1e-10      NA 0.5385696 2.453345  612
setwd("F:/修稿文件汇总/改稿进展/据_结果_程序_整理/结果四(更新a-d图）/Fig 4D/修稿修改/")
cell=c("NK","Neutrophil","DC")
for(i in 1:3){
  all_gene=read.table(paste0(cell[i],"_alllnc_rank_results_marker_median_new-abs-new.txt"),header=T,sep="\t",row.names = 1)
  all_gene[,2]=abs(all_gene[,2])
  bb=which(table(all_gene[,2])<2)
  bb=as.data.frame(bb)
  bbb=all_gene[which(all_gene[,2] %in% rownames(bb)),]
  aa=which(table(all_gene[,2])>1)
  aaaa=c()
  for(j in 1:length(aa)){
    num=length(which(all_gene[,2]==names(aa)[j]))
    for(n in 1:num){
      aaa=cbind(all_gene[which(all_gene[,2]==names(aa)[j]),][n,1],all_gene[which(all_gene[,2]==names(aa)[j]),][n,2]*(1-0.001*n))
      aaaa=rbind(aaaa,aaa)
    }
  }

  write.table(bbb,paste0(cell[i],"_alllnc_rank_results_marker_median_new-abs-new1.txt"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(aaaa,paste0(cell[i],"_alllnc_rank_results_marker_median_new-abs-new1.txt"),row.names=FALSE,col.names=FALSE,sep="\t",append=T,quote=FALSE)
}
  




fgseaRes2 <- fgsea(sam, rank[1:10], maxSize=500)
data(examplePathways)
data(exampleRanks)
fgseaRes <- fgsea(examplePathways, exampleRanks, maxSize=500)
# Testing only one pathway is implemented in a more efficient manner
fgseaRes1 <- fgsea(examplePathways[1], exampleRanks)
write.table(examplePathways[1],"examplePathways.txt",quote=F,sep="\t")
write.table(exampleRanks,"exampleRanks.txt",quote=F,sep="\t")













****************************************************************************
Fig 5
*****************
****Fig 5A****注意！！！！这个数据里面好像包含了正常样本，在作图和后续分析中记得把癌症样本和正常样本分开

setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5A/")
source('CIBERSORT.R')
lncRNA=read.table("lncRNA_exp_14cell.txt",sep="\t",header=T,row.names=1)
cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
for(i in 1:33){
  cancer_all=read.table(paste0("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5A/cancer-new/TCGA-",cs[i],"-htseq_fpkm-uq.txt"),sep="\t",header=T,row.names=1)
  cancer_lnc=cancer_all[rownames(lncRNA),]
  cancer_log=log2(cancer_lnc+0.001)
  write.table(cancer_log,paste0(cs[i],"_lnc_log2_exp.txt"),sep="\t",quote=F)
}
cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
#cs=c("HNSC","KIRC","LGG","LUAD","LUSC","PRAD","THCA","UCEC")
cell=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg")
for(i in 1:33){
   for(j in 1:14){

   results <- CIBERSORT(paste0(cell[j],"_log2_annova_FDR_0.3.txt"),paste0(cs[i],"_lnc_log2_exp.txt"), 100, TRUE)
   res=results[,cell[j]]
   write.table(t(res),paste0(cs[i],"_14cell_pre_CIBEERSORT_result.txt"),sep="\t",quote=F,append=T,col.names = F,row.names = F)
 }
}


###整合每种癌症的均值为每种癌症的免疫细胞浸润水平,画雷达图
###！！！！！要挑出癌症样本出来再画图
cell=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg")

cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
res=c()
for(i in 1:33){
wd1="F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5A/cnacer_LNCCIBERSORT_results/"
wd2="F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5A/33_cancer_label/"
data=read.table(paste0(wd1,cs[i],"_14cell_pre_CIBEERSORT_result.txt"),sep="\t",header=T,row.names=1)
label=read.table(paste0(wd2,"TCGA-",cs[i],".txt"),sep="\t",header=T)
number_cancer=gsub("-",".",label[which(label[,2]=="Tumor"),1])
number_normal=gsub("-",".",label[which(label[,2]=="Normal"),1])
cancer=data[cell,number_cancer]
cancer1=apply(cancer,1,mean)
res=cbind(res,cancer1)
}
colnames(res)=cs
write.table(res,"33cancer_14cell_CIBER_res_mean.txt",sep="\t",quote=F)


library(fmsb)
maxmin=matrix(rep(c(0.5,0),each=33),byrow=T,ncol=33)
colnames(maxmin)=names(res)
Data=rbind(maxmin,res)
Data=as.data.frame(Data)
rownames(Data)
radarchart(Data, seg=6,plty=1,pty=32,plwd=3,cglty=2,pcol=c("B"="PaleTurquoise1","DC"="red","Granulocytes"="grey","M0"="Goldenrod1","M1"="Sienna3","M2"="Sienna4","Macrophage"="Burlywood1",
                                 "Monocyte"="PaleGreen","Neutrophil"="Plum2","NK"="RosyBrown1","CD8.T"="DarkOliveGreen3","CD4.T"="MediumPurple2","CD3.T"="Sienna1","Treg"="black"))
text(2,c(0.10,0.12,0.2,0.28,0.36,0.44,0.52,0.60,0.68,0.76,0.84,0.92,0.96,1),c("—— B","—— DC","—— Granulocytes","—— M0","—— M1","—— M2","—— Macrophage","—— Monocyte","—— Neutrophil","—— Nk","—— CD8 T","—— CD4 T","—— CD3 T","—— Treg"),
col=c("PaleTurquoise1","red","grey","Goldenrod1","Sienna3","Sienna4","Burlywood1","PaleGreen","Plum2","RosyBrown1","DarkOliveGreen3","MediumPurple2","Sienna1","black"))
max(res)
[1] 0.4605035
****************************************
********Fig 5B***
####将免疫细胞浸润水平和肿瘤纯度计算相关性，气泡图展示
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5B/")
tur=read.table("tum_pur.txt",sep="\t",header=T,row.names=1)
cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
nn=gsub("-",".",substring(rownames(tur),1,16))
#write.table(nn,"tur_sample_names.txt",sep="\t",quote=F)
tur=read.table("tum_pur.txt",sep="\t",header=T,row.names=1)
RES_cancer=c()
res2=c()
for(i in 1:33){
  wd1="F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5A/cnacer_LNCCIBERSORT_results/"
  wd2="F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5A/33_cancer_label/"
  data=read.table(paste0(wd1,cs[i],"_14cell_pre_CIBEERSORT_result.txt"),sep="\t",header=T,row.names=1)
  label=read.table(paste0(wd2,"TCGA-",cs[i],".txt"),sep="\t",header=T)
  number_cancer=gsub("-",".",label[which(label[,2]=="Tumor"),1])
  number_normal=gsub("-",".",label[which(label[,2]=="Normal"),1])
  cancer=data[cell,number_cancer]
  res_cell=c()
 
   for(j in 1:14){
    cell_cancer_cor=cor.test(as.numeric(cancer[cell[j],]),as.numeric(tur[colnames(cancer),4]))
    res_cell=cbind(res_cell,cell_cancer_cor$estimate[[1]])
    plot_res=c(cs[i],cell[j],cell_cancer_cor$estimate[[1]],cell_cancer_cor$p.value)
    res2=rbind(res2,plot_res)

  }
  RES_cancer=rbind(RES_cancer,res_cell)
}
  rownames(RES_cancer)=cs
  colnames(RES_cancer)=cell
  write.table(RES_cancer,"cor_cell_cancer.txt",sep="\t",quote=F,row.names=T,col.names=T)
  write.table(res2,"cor_cell_cancer_plot_res.txt",sep="\t",quote=F,row.names=F,col.names=F)

  res2=read.table("cor_cell_cancer_plot_res.txt",sep="\t",header=F)
library(ggplot2)
ggplot(res2,aes(x=res2[,3],y=res2[,1],size=-log10(res2[,4]),col=res2[,2]))+geom_point()+theme_bw()






**********************
***Fig 5C*******
####相似性聚类
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5C/")
data=read.table("33cancer_14cell_CIBER_res_mean.txt",sep="\t",header=T,row.names=1)
res=cor(data,data)##可能要转置
library(pheatmap)
Data=read.table("33cancer_cor.txt",sep="\t",header=T,row.names=1)
pheatmap(res, show_colnames= T, show_rownames= T, scale= "none", fontsize= 6.5,
         
         clustering_method ="complete",
         cluster_cols=T,cluster_rows=T,
         breaks = c(seq(0.51,0.8,length.out=50),seq(0.81,1,length.out=50)),
         
         col = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#c("navy", "white", "firebrick3")
#c("#A4D3EE","white","#D8BFD8")

####加上那几个癌症的相关性散点图
library(ggplot2)
data=read.table("33cancer_14cell_CIBER_res_mean.txt",sep="\t",header=T,row.names=1)
cs=c("STAD","ESCA","LUAD","LUSC","UCS","UCEC","COAD","READ")
for(i in c(1,3,5,7)){
res=cor.test(data[,cs[i]],data[,cs[i+1]])
res_p=res[3][[1]]
res_cor=res[4][[1]]
p <- ggplot(data,aes(x=data[,cs[i]],y=data[,cs[i+1]])) + geom_point(colour="red",size=3)
p1=p + geom_smooth(method='lm',colour='black') +
  theme(axis.text=element_text(colour = 'black',size = 12)) +
  geom_text(aes(x = 0.12, y = 0.25, label = paste0("p=",res_p,"   r=",res_cor)))+
  labs(x=cs[i],y=cs[i+1])+theme_bw()
ggsave(p1,file=paste0(cs[i],"_",cs[i+1],"_cor_plot.pdf"))
}



*******************************
*****FIg 5D***
####正常样本和癌症样本的浸润水平的差异性
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5D/修稿修改/")
#cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
cs=c("BRCA","CHOL","COAD","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")

#res2=c()
cell=c("B","DC","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T")
for(i in 1:14){
  wd1="F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5D/修稿修改/cnacer_LNCCIBERSORT_results/"
  wd2="F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5D/修稿修改/33_cancer_label/"
  data=read.table(paste0(wd1,cs[i],"_14cell_pre_CIBEERSORT_result.txt"),sep="\t",header=T,row.names=1)
  label=read.table(paste0(wd2,"TCGA-",cs[i],".txt"),sep="\t",header=T)
  number_cancer=gsub("-",".",label[which(label[,2]=="Tumor"),1])
  number_normal=gsub("-",".",label[which(label[,2]=="Normal"),1])
  cancer=data[cell,number_cancer]
  normal=data[cell,number_normal]
  

    for(j in 1:12){
        res=c()
        TT=t.test(as.numeric(cancer[j,]),as.numeric(normal[j,]))
        TP=TT$p.value
        WT=wilcox.test(as.numeric(cancer[j,]),as.numeric(normal[j,]))
        WP=WT$p.value
        FC=mean(as.numeric(cancer[j,]))/mean(as.numeric(normal[j,]))
        res=cbind(cs[i],cell[j],TP,FC)
        #res=cbind(cs[i],cell[j],WP,FC)
        write.table(res,"14cancer_T_P_FC.txt",sep="\t",row.names=F,col.names=F,quote=F,append=T)
      }
    }
  

  ###改进res111.pdf
  library(ggplot2)
  res=read.table("14cancer_T_P_FC1.txt",sep="\t",header=T)
  res[which(res[,4]>2),4]=2
  res[intersect(which(res[,4]>1),which(res[,4]<2)),4]=1.5
  res[which(res[,4]<0.5),4]=0.3
  res[intersect(which(res[,4]>0.5),which(res[,4]<1)),4]=0.8
  ggplot(as.data.frame(res),aes(cell,Cancer))+geom_point(aes(size=-log(Pvalue,10),color=as.factor(FC)))+theme_bw()
  
  
*****************************************8
  ********************8
  ####Fig 5E#####
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5E/")

  cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  
  #res2=c()
  cell=c("B","DC","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T")
  for(i in 1:33){
    wd1="F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5E/cnacer_LNCCIBERSORT_results/"
    wd2="F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5E/33_cancer_label/"
    data=read.table(paste0(wd1,cs[i],"_14cell_pre_CIBEERSORT_result.txt"),sep="\t",header=T,row.names=1)
    label=read.table(paste0(wd2,"TCGA-",cs[i],".txt"),sep="\t",header=T)
    number_cancer=gsub("-",".",label[which(label[,2]=="Tumor"),1])
    #number_normal=gsub("-",".",label[which(label[,2]=="Normal"),1])
    cancer=data[cell,number_cancer]
    #normal=data[cell,number_normal]
    wd3="F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5E/TCGA_new_survival/"
    sur=read.table(paste0(wd3,"TCGA-",cs[i],".survival.tsv"),sep="\t",header=T,row.names = 1)
    rownames(sur)=gsub("-",".",rownames(sur))
    cancer_sample=intersect(colnames(cancer),rownames(sur))
    days=sur[cancer_sample,3]
    status=sur[cancer_sample,1]
    y=Surv(days,status)
    Cancer_data=cancer[,cancer_sample]
    for(j in 1:12){
      res=c()
      c1=coxph(y~as.numeric(Cancer_data[j,]),data=Cancer_data)
      cc1=summary(c1)
      res=cbind(cs[i],cell[j],cc1[7][[1]][2],cc1[7][[1]][5])
      write.table(res,"33cancer_HR_P.txt",sep="\t",row.names=F,col.names=F,quote=F,append=T)
    }
}


  setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5E/")
HR_P=read.table("33cancer_HR_P.txt",sep="\t",header=T)
res=HR_P[which(HR_P[,4]<0.1),]
library(ggplot2)
##第三列是HR第四列是P
res[which(res[,3]>1),5]=2
res[intersect(which(res[,3]>1),which(res[,4]<0.05)),5]=1.5
res[which(res[,3]<1),5]=0.3
res[intersect(which(res[,3]<1),which(res[,4]<0.05)),5]=0.8
res[which(res[,3]<1),6]=0.5
res[intersect(which(res[,3]>=1),which(res[,3]<10)),6]=5
res[intersect(which(res[,3]>=10),which(res[,3]<1000)),6]=500
res[intersect(which(res[,3]>=1000),which(res[,3]<10000)),6]=5000
res[which(res[,3]>=10000),6]=50000
ggplot(as.data.frame(res),aes(x=Cell,y=Cancer))+
  geom_point(aes(size=as.factor(res[,6]),color=as.factor(res[,5])))+theme_bw()


*****************************************8
********************8
####Fig 5E#####
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5E/")

cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

#res2=c()
cell=c("B","DC","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T")
for(i in 1:33){
  wd1="F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5E/cnacer_LNCCIBERSORT_results/"
  wd2="F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5E/33_cancer_label/"
  data=read.table(paste0(wd1,cs[i],"_14cell_pre_CIBEERSORT_result.txt"),sep="\t",header=T,row.names=1)
  label=read.table(paste0(wd2,"TCGA-",cs[i],".txt"),sep="\t",header=T)
  number_cancer=gsub("-",".",label[which(label[,2]=="Tumor"),1])
  #number_normal=gsub("-",".",label[which(label[,2]=="Normal"),1])
  cancer=data[cell,number_cancer]
  #normal=data[cell,number_normal]
  wd3="F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5E/TCGA_new_survival/"
  sur=read.table(paste0(wd3,"TCGA-",cs[i],".survival.tsv"),sep="\t",header=T,row.names = 1)
  rownames(sur)=gsub("-",".",rownames(sur))
  cancer_sample=intersect(colnames(cancer),rownames(sur))
  days=sur[cancer_sample,3]
  status=sur[cancer_sample,1]
  y=Surv(days,status)
  Cancer_data=cancer[,cancer_sample]
  for(j in 1:12){
    res=c()
    c1=coxph(y~as.numeric(Cancer_data[j,]),data=Cancer_data)
    cc1=summary(c1)
    res=cbind(cs[i],cell[j],cc1[7][[1]][2],cc1[7][[1]][5])
    write.table(res,"33cancer_HR_P.txt",sep="\t",row.names=F,col.names=F,quote=F,append=T)
  }
}


setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5E/")
HR_P=read.table("33cancer_HR_P.txt",sep="\t",header=T)
res=HR_P[which(HR_P[,4]<0.1),]
library(ggplot2)
##第三列是HR第四列是P
res[which(res[,3]>1),5]=2
res[intersect(which(res[,3]>1),which(res[,4]<0.05)),5]=1.5
res[which(res[,3]<1),5]=0.3
res[intersect(which(res[,3]<1),which(res[,4]<0.05)),5]=0.8
res[which(res[,3]<1),6]=0.5
res[intersect(which(res[,3]>=1),which(res[,3]<10)),6]=5
res[intersect(which(res[,3]>=10),which(res[,3]<1000)),6]=500
res[intersect(which(res[,3]>=1000),which(res[,3]<10000)),6]=5000
res[which(res[,3]>=10000),6]=50000
ggplot(as.data.frame(res),aes(x=Cell,y=Cancer))+
  geom_point(aes(size=as.factor(res[,6]),color=as.factor(res[,5])))+theme_bw()






******************************************************************************************************
  ***********Fig 6********
  *************************************************************
  ****Fig 6***
  #####Fig 6A
  library(preprocessCore)##分位数标准化需要的包
library(pheatmap)
library(gplots)
library(vegan)
library(permute)
library(lattice)
BiocManager::install("vegan")
library(pheatmap)
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果六/Fig 6A/修稿修改/")
data=read.table("Tumorsample_Sur_lncRNACibersort_immunefiltration_merge.txt",sep="\t",header=T,row.names = 1)
exp_data=as.matrix(data[,c(2,4:7,13,12,11,1,10,8,9)])
#png("cancers_ciber_all33_new.png")
cluster=data[,18]
annotation_col=data.frame(sample=factor(as.numeric(cluster)))
unique(annotation_col)
annotation_colors =list(sample=c("1"="red","2"="MediumPurple2","3"="PaleGreen","4"="Plum2","5"="DarkOliveGreen3","6"="Gold2"))
rownames(exp_data)=rownames(annotation_col)
pheatmap(exp_data, show_colnames= T, show_rownames= F, scale= "none", fontsize= 6.5,
         clustering_method ="complete",
         cluster_cols=F,cluster_rows=T,
         annotation_row= annotation_col, 
         annotation_colors= annotation_colors, 
         breaks = c(seq(0,0.1,length.out=50),seq(0.11,0.6,length.out=50)),
         col = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#dev.off()
> table(cluster)
cluster
1    2    3    4    5    6 
3210 4161 1778  535  207  230 

list=pheatmap(exp_data,cluster_cols=F,cluster_rows=T,scale = "none")
row_cluster=cutree(list$tree_row,k=6)
write.table(row_cluster,"row_cluster_K6.txt",row.names=T,col.names=F,quote=F,sep="\t")
##聚类图中的顺序
newOrder=exp_data[list$tree_row$order,]
newOrder[,ncol(newOrder)+1]=row_cluster[match(rownames(newOrder),names(row_cluster))]
colnames(newOrder)[ncol(newOrder)]="Cluster"
write.table(newOrder,"row_cluster_order_K6.txt",row.names=T,col.names=F,quote=F,sep="\t")


###Xcell和Cibersort的结果也聚类一下,更新结果：用他们预测出来的所有细胞类型聚类
library(pheatmap)
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/加-改的附图/CIBERSORT_xcell/更新Fig6部分")
data=read.table("xcell_all.txt",sep="\t",header=T,row.names = 1)
data=data[which(data[,68]=="Tumor"),]
exp_data=as.matrix(data[,1:64])
list=pheatmap(exp_data,cluster_cols=F,cluster_rows=T,scale = "none")
row_cluster=cutree(list$tree_row,k=10)
write.table(cbind(exp_data,row_cluster),"xcell_row_cluster_K6.txt",row.names=T,col.names=T,quote=F,sep="\t")
##聚类图中的顺序
table(row_cluster)
row_cluster
1    2    3    4    5    6 
9892   11  277   78   33   36 
data=read.table("xcell_row_cluster_K6.txt",sep="\t",header=T,row.names = 1)
exp_data=as.matrix(data[,1:64])

#png("cancers_ciber_all33_new.png")
cluster=data[,65]
annotation_col=data.frame(sample=factor(as.numeric(cluster)))
unique(annotation_col)
annotation_colors =list(sample=c("1"="red","2"="MediumPurple2","3"="PaleGreen","4"="Plum2","5"="DarkOliveGreen3","6"="Gold2"))
rownames(exp_data)=rownames(annotation_col)
pheatmap(exp_data, show_colnames= T, show_rownames= F, scale= "none", fontsize= 6.5,
         clustering_method ="complete",
         cluster_cols=T,cluster_rows=T,
         annotation_row= annotation_col, 
         annotation_colors= annotation_colors, 
         breaks = c(seq(0,0.05,length.out=50),seq(0.051,0.3,length.out=50)),
         col = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#dev.off()
> table(cluster)
cluster
1    2    3    4    5    6 
9509   93  245  457   20    3 
> row_cluster=cutree(list$tree_row,k=15)
> table(row_cluster)
row_cluster
row_cluster
1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
9416  150    8  268   59  255   39   32   22    3    7    9   11   12   36 
###cibersort
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/加-改的附图/CIBERSORT_xcell/更新Fig6部分")
data=read.table("Cibersort_all.txt",sep="\t",header=T,row.names = 1)
data=data[which(data[,26]=="Tumor"),]
exp_data=as.matrix(data[,1:22])
list=pheatmap(exp_data,cluster_cols=F,cluster_rows=T,scale = "none")
row_cluster=cutree(list$tree_row,k=6)
write.table(cbind(exp_data,row_cluster),"cibersort_row_cluster_K6.txt",row.names=T,col.names=T,quote=F,sep="\t")
##聚类图中的顺序
table(row_cluster)
> table(row_cluster)
row_cluster
1    2    3    4    5    6 
2162 4614 1167 1153  974  257 
data=read.table("cibersort_row_cluster_K6.txt",sep="\t",header=T,row.names = 1)
exp_data=as.matrix(data[,1:22])


#png("cancers_ciber_all33_new.png")
cluster=data[,23]
annotation_col=data.frame(sample=factor(as.numeric(cluster)))
unique(annotation_col)
annotation_colors =list(sample=c("1"="red","2"="MediumPurple2","3"="PaleGreen","4"="Plum2","5"="DarkOliveGreen3","6"="Gold2"))
rownames(exp_data)=rownames(annotation_col)
pheatmap(exp_data, show_colnames= T, show_rownames= F, scale= "none", fontsize= 6.5,
         clustering_method ="complete",
         cluster_cols=F,cluster_rows=T,
         annotation_row= annotation_col, 
         annotation_colors= annotation_colors, 
         breaks = c(seq(0,0.05,length.out=50),seq(0.051,0.9,length.out=50)),
         col = colorRampPalette(c("navy", "white", "firebrick3"))(100))




#*******************************************************
############################
##Fig 6B
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果六/Fig 6B/")
data=read.table("Tumorsample_Sur_lncRNACibersort_immunefiltration_merge.txt",sep="\t",header=T,row.names=1)
clu_ratio=matrix(0,ncol=7,nrow=33)
cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
cluster=unique(data[,18])
for(i in 1:33){
  cancer=data[which(data[,17]==cs[i]),18]
  cancer_clu=table(cancer)
  clu_ratio[i,7]=length(cancer)
  for(j in names(cancer_clu)){
    n=which(names(cancer_clu)==j)
    clu_ratio[i,as.numeric(j)]=cancer_clu[[n]]/length(cancer)
  }
}

colnames(clu_ratio)=c("C1","C2","C3","C4","c5","C6","samples")
rownames(clu_ratio)=cs

barplot(t(as.matrix(clu_ratio[,1:6])),space=0.5,col=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2"),width=clu_ratio[,7]/100)
#barplot(t(as.matrix(clu_ratio[,1:6])),yaxt="n",space=0.5,col=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2"),width=clu_ratio[,7]/100)


###xcell和Cibersort结果
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/加-改的附图/CIBERSORT_xcell/")
data=read.table("xcell_row_cluster_K6.txt",sep="\t",header=T,row.names=1)
data1=read.table("xcell.txt",sep="\t",header=T,row.names=1)
data1=data1[which(data1[,15]=="Tumor"),]
clu_ratio=matrix(0,ncol=7,nrow=33)
cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
cluster=unique(data[,12])
for(i in 1:33){
  cancer=data[which(data1[,12]==cs[i]),12]
  cancer_clu=table(cancer)
  clu_ratio[i,7]=length(cancer)
  for(j in names(cancer_clu)){
    n=which(names(cancer_clu)==j)
    clu_ratio[i,as.numeric(j)]=cancer_clu[[n]]/length(cancer)
  }
}

colnames(clu_ratio)=c("C1","C2","C3","C4","c5","C6","samples")
rownames(clu_ratio)=cs

barplot(t(as.matrix(clu_ratio[,1:6])),space=0.5,col=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2"),width=clu_ratio[,7]/100)


data=read.table("cibersort_row_cluster_K6.txt",sep="\t",header=T,row.names=1)
data1=read.table("Cibersort.txt",sep="\t",header=T,row.names=1)
data1=data1[which(data1[,11]=="Tumor"),]
clu_ratio=matrix(0,ncol=7,nrow=33)
cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
#cluster=unique(data[,12])
for(i in 1:33){
  cancer=data[which(data1[,8]==cs[i]),8]
  cancer_clu=table(cancer)
  clu_ratio[i,7]=length(cancer)
  for(j in names(cancer_clu)){
    n=which(names(cancer_clu)==j)
    clu_ratio[i,as.numeric(j)]=cancer_clu[[n]]/length(cancer)
  }
}

colnames(clu_ratio)=c("C1","C2","C3","C4","c5","C6","samples")
rownames(clu_ratio)=cs

barplot(t(as.matrix(clu_ratio[,1:6])),space=0.5,col=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2"),width=clu_ratio[,7]/100)

*********************************************
  #################################333
  ###Fig 6C
  setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果六/Fig 6C/修稿修改/")
cancer=read.table("Tumorsample_Sur_lncRNACibersort_immunefiltration_merge.txt",sep="\t",header=T,row.names=1)
res1=substring(rownames(cancer),1,12)
write.table(res1,"names1.txt",sep="\t",quote=F)
score=read.table("HDR_new.txt",header = T,sep="\t",row.names=1)
rownames(score)=gsub("-",".",rownames(score))
cancer=read.table("Tumorsample_Sur_lncRNACibersort_immunefiltration_merge.txt",sep="\t",header=T,row.names=1)
inter=intersect(rownames(cancer),rownames(score))
cancer1=cancer[inter,18]
score1=score[inter,5]
boxplot(score1~cancer1,xlab="Cluster",ylab="HRD",outline=FALSE,boxwex = 0.6,col=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2"))

score=read.table("Tumor_CYT.txt",header = T,sep="\t",row.names=1)
rownames(score)=gsub("-",".",rownames(score))
cancer=read.table("Tumorsample_Sur_lncRNACibersort_immunefiltration_merge.txt",sep="\t",header=T,row.names=1)
inter=intersect(rownames(cancer),rownames(score))
cancer1=cancer[inter,18]
score1=score[inter,3]
boxplot(score1~cancer1,xlab="Cluster",ylab="CYT",outline=FALSE,boxwex = 0.6,col=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2"))

score=read.table("Tumor_score.txt",header = T,sep="\t",row.names=1)
rownames(score)=gsub("-",".",rownames(score))
cancer=read.table("Tumorsample_Sur_lncRNACibersort_immunefiltration_merge.txt",sep="\t",header=T,row.names=1)
inter=intersect(rownames(cancer),rownames(score))
cancer1=cancer[inter,18]
score1=score[inter,3]
boxplot(score1~cancer1,xlab="Cluster",ylab="Immune score",outline=FALSE,boxwex = 0.6,col=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2"))



*********************************************
  #####Fig 6D
  #####
library(survival)
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果六/Fig 6D/修稿修改/")
data=read.table("Tumorsample_Sur_lncRNACibersort_immunefiltration_merge.txt",sep="\t",header=T,row.names = 1)
cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
y=Surv(data$OS.time,data$OS_event)
kmfit2<-survfit(y~data$Cluster, data=data) 
plot(kmfit2, main="Kaplan-Meier curves of six groups", xlab="Time(days)",  ylab="Overall survival", mark=3,mark.time=T,col=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2"))
text(10000,c(0.8,0.84,0.88,0.92,0.96,1),c("—— 1","—— 2","—— 3","—— 4","—— 5","—— 6"),col=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2"))
text(6000,0.8,"p= <2e-16") 
##计算P值
r<-survdiff(formula = y~data$Cluster,  rho = 0,data=data) 
p<-1 - pchisq(r$chisq, length(r$n) - 1)

###分别画33个癌症的生存结果，去掉样本数少于10个的类
color=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2")
for(i in 1:33){
  data1=data[which(data[,17]==cs[i]),]
  a=unique(data1[,18])
  aa=c()
  for(n in 1:length(a)){
    cc=length(which(data1[,18]==a[n]))
    if(cc>=10){
      aa=c(aa,a[n])
    }
  }
  if(length(aa)>1){
    data1=data1[which(data1[,18]%in%aa),]
    y=Surv(data1$OS.time,data1$OS_event)
    kmfit2<-survfit(y~data1$Cluster, data=data1) 
    main_name<-paste("module_",cs[i],sep="")
    data.survdiff=survdiff(y~data1$Cluster,data = data1)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    pdf(paste(main_name,"_survival.pdf",sep=''))
    plot(kmfit2,main=main_name,xlab="survival time",ylab="survival probability",mark=3,mark.time=T,cex.lab=1.5,cex.main=1.5,col=color[aa],lwd=5)
    legend("bottomright",col=color[aa],legend=aa,lty=1,lwd=3)
    text(max(data1$OS.time)-1000,0.8,paste0("P=",p.val))
    dev.off()
  }
}



  
  
  
  
  


****************************************************************
****Fig S2****
####先计算出每个特征集中的lncRNA在其他几个集合中出现的次数和比例
####先计算出每个特征集中的lncRNA在其他几个集合中出现的次数和比例
setwd("D:/0-课题/1-1-LncRNACIBERSORT/修稿文件汇总/改稿进展/数据_结果_程序_整理/加-改的附图/Fig S2/")
cell=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg")
RES=c()
for(i in 1:14){
sig=read.table(paste0(cell[i],"_log2_annova_FDR_0.3.txt"),sep="\t",header=T,row.names=1)
res=rownames(sig)
RES=c(RES,res)
}
RES=unique(RES)
###2187
Data_res=c()
for(j in 1:length(RES)){
  lnc1=RES[j]
  count=0
for(n in 1:14){
  sig=read.table(paste0(cell[n],"_log2_annova_FDR_0.3.txt"),sep="\t",header=T,row.names=1)
  res=rownames(sig)
  if(lnc1%in%res){
    count=count+1
    }
  }
Data=c(lnc1,count)
Data_res=rbind(Data_res,Data)
}

Num_Res=c()
Ratio_Res=c()
for(n in 1:14){
  sig=read.table(paste0(cell[n],"_log2_annova_FDR_0.3.txt"),sep="\t",header=T,row.names=1)
  sig_num=Data_res[which(Data_res[,1]%in%rownames(sig)),]
  num_res=cbind(cell[n],table(sig_num[,2]))
  ratio_res=cbind(cell[n],table(sig_num[,2])/length(sig[,1]))

  ###比例
  Num_Res=rbind(Num_Res,num_res)
  Ratio_Res=rbind(Ratio_Res,ratio_res)

}
write.table(Ratio_Res,"sig_lnc_number_ratio_res.txt",sep="\t",quote=F)
write.table(Num_Res,"sig_lnc_number_res.txt",sep="\t",quote=F)

A=read.table("sig_lnc_number_ratio_res.txt",sep="\t",header=T)
AA=matrix(rep(0,14*14),ncol=14,nrow=14)
col=c(1:14)
rownames(AA)=cell
colnames(AA)=col

for(i in 1:length(A[,1])){
a=which(cell==A[i,1])
b=which(col==A[i,2])
AA[a,b]=A[i,3]
}
barplot(t(AA),space=0.5, angle =45,col=c("1"="PaleTurquoise1","2"="HotPink1","3"="MediumOrchid1","4"="Chartreuse1","5"="DarkSeaGreen1","6"="OliveDrab1",
                                 "7"="OliveDrab","8"="DarkTurquoise","9"="Brown1","10"="Khaki1","11"="LightYellow1","12"="Sienna1",
                                 "13"="Gold1","14"="Goldenrod2"))
text(20,c(400,425,450,475,500,525,550,575,600,625,650,675,700,725),
c("——1","——2","——3","——4","——5","——6","——7","——8","——9","——10","——11","——12","——13","——14"),
col=c("PaleTurquoise1","HotPink1","MediumOrchid1","Chartreuse1","DarkSeaGreen1","OliveDrab1",
                                 "OliveDrab","DarkTurquoise","Brown1","Khaki1","LightYellow1","Sienna1",
                                 "Gold1","Goldenrod2"))


*********************************
****Fig S3****
###计算几个特征集之间的交集，提出交集的表达谱，画连线的箱式图
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/加-改的附图/Fig S3/")
sig_M0=read.table("M0_log2_annova_FDR_0.3.txt",sep="\t",header=TRUE,row.names=1)
sig_M01=rownames(sig_M0)
sig_M1=read.table("M1_log2_annova_FDR_0.3.txt",sep="\t",header=TRUE,row.names=1)
sig_M11=rownames(sig_M1)
sig_M2=read.table("M2_log2_annova_FDR_0.3.txt",sep="\t",header=TRUE,row.names=1)
sig_M21=rownames(sig_M2)
sig_CD4T=read.table("CD4 T_log2_annova_FDR_0.3.txt",sep="\t",header=TRUE,row.names=1)
sig_CD4T1=rownames(sig_CD4T)
sig_CD8T=read.table("CD8 T_log2_annova_FDR_0.3.txt",sep="\t",header=TRUE,row.names=1)
sig_CD8T1=rownames(sig_CD8T)

M=intersect(sig_M01,intersect(sig_M11,sig_M21))
T=intersect(sig_CD4T1,sig_CD8T1)

M_exp=sig_M0[M,3:5]
for(i in 1:length(M_exp[,1])){
for(j in 1:3){
res=c(rownames(M_exp)[i],colnames(M_exp)[j],M_exp[i,j])
write.table(t(res),"M_366lnc_exp_res.txt",sep="\t",quote=F,append=T,row.names=F,col.names=F)
}
}

T_exp=sig_CD4T[T,c(7,9)]
for(i in 1:length(T_exp[,1])){
for(j in 1:2){
res=c(rownames(T_exp)[i],colnames(T_exp)[j],T_exp[i,j])
write.table(t(res),"T_245lnc_exp_res.txt",sep="\t",quote=F,append=TRUE,row.names=F,col.names=F)
}
}
write.table(M_exp,"M_co-lnc_exp_res.txt",sep="\t",quote=F,row.names=TRUE,col.names=TRUE)
write.table(T_exp,"T_co-lnc_exp_res.txt",sep="\t",quote=F,row.names=TRUE,col.names=TRUE)

MM=read.table("M_366lnc_exp_res.txt",sep="\t",header=F)
colnames(MM)=c("gene","cell","exp")
TT=read.table("T_245lnc_exp_res.txt",sep="\t",header=F)
colnames(TT)=c("gene","cell","exp")

wp1=wilcox.test(M_exp[,1],M_exp[,2])
p1=wp1$p.value
wp2=wilcox.test(M_exp[,1],M_exp[,3])
p2=wp2$p.value
wp3=wilcox.test(M_exp[,3],M_exp[,2])
p3=wp3$p.value
wp4=wilcox.test(T_exp[,1],T_exp[,2])
p4=wp4$p.value
> p1(M0-M1)
[1] 0.8576989
> p2(M0-M2)
[1] 0.398203
> p3(M1-M2)
[1] 0.3024667
> p4(TCD4-CD8)
[1] 0.0001860221


M_exp=read.table("M_co-lnc_exp_res.txt",sep="\t",header=T,row.names=1)
T_exp=read.table("T_co-lnc_exp_res.txt",sep="\t",header=T,row.names=1)
library(ggpubr)
library(ggplot2)
ggpaired(T_exp, cond1 = "TCD4", cond2 = "TCD8",
fill = "condition",color="condition" ,line.color = "gray", line.size = 0.4)







************************************
***Fig S4****
####用共表达的基因集合富集到免疫相关功能
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/加-改的附图/Fig S4/")
cell=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg")
ID=read.table("DICE_57820_RNA_annotation.txt",sep="\t",header=T,row.names=1)
ID=ID[which(ID[,3]=="protein_coding"),]
target=read.table("Gene-from_ImmPort.txt",header=T,sep="\t")
term=unique(target[,2])
for(i in 1:length(term)){
   tar_gene=target[which(target[,2]==term[i]),1]
for(n in 1:14){
  Gene <- read.table(paste0(cell[n],"_new_co-expression-mRNA-exp-1.txt"),header=T,sep="\t")
  Genelist=ID[Gene[,1],2]
  inter=length(intersect(Genelist,tar_gene))
  #P=1-phyper(inter,length(hall[,2]),(64914-length(hall[,2])),length(tar_gene))
  P=phyper(inter-1,length(Genelist),(64914-length(Genelist)),length(tar_gene),lower.tail=FALSE)
  print(i)
  res=c(term[i],cell[n],P)
  write.table(t(res),"P_value_new-1.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE,append=T)
}
}

P=read.table("P_value_new-1.txt",sep="\t",header=FALSE)
P[,3]=p.adjust(P[,3],method="fdr",n=length(P[,3]))
P=P[-which(P[,3]>0.05),]
colnames(P)=c("term","cell","FDR")
ggplot(P,aes(x=cell,y=term,color=cell)) +geom_point(aes(size=-log10(FDR)))+
scale_size_continuous(breaks = c(2,5,10,20), guide = guide_legend())+
theme(axis.text.x = element_text(angle = 45, hjust = 0.5))+theme_bw()

###########
*****Fig S5 
###S5在F:\修稿文件汇总\改稿进展\数据_结果_程序_整理\加-改的附图\Fig S5\波-circos中
#文件见/base
#程序是circos.conf

******************************
  ***Fig S6-S7
###14种免疫细胞在每种癌症中的平均浸润情况
library(plyr)
setwd("C:/Users/Administrator/Desktop/心慧/cnacer_LNCCIBERSORT_results/")
cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")#res2=c()
cell=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg")
for(i in 1:14){
  list1<-list() 
  M=c()
  for(j in 1:length(cs)){
    Data=read.table(paste0(cs[j],"_14cell_pre_CIBEERSORT_result.txt"),sep="\t",header=T,row.names=1) 
    Data=t(Data)
    list1[[j]]=as.data.frame(t(as.data.frame(Data[,i])))
    m=median(Data[,i])
    M=rbind(M,m)   
  }
  A=do.call(rbind.fill,list1)
  rownames(A)=cs
  B=arrange(A,M)  
  cs_new=arrange(as.data.frame(cs),M)
  pdf(paste("C://Users/Administrator/Desktop/心慧/Cancer_immunefiltratuion_boxplot/",cell[i],"_33cancer_boxplot.pdf",sep=''))
  boxplot(t(B),outline=F,names=t(cs_new),las=2)
  dev.off()
}

###Fig S7
###肿瘤样本的免疫浸润提取
cs=c("BRCA","CHOL","COAD","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")
cell=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg")
for(j in 1:14){
  wd1="C://Users/Administrator/Desktop/心慧/cnacer_LNCCIBERSORT_results/"
  wd2="C://Users/Administrator/Desktop/心慧/33_cancer_label/"
  wd3="C://Users/Administrator/Desktop/心慧/Sample_immunefiltratuion_group/"
  data=read.table(paste0(wd1,cs[j],"_14cell_pre_CIBEERSORT_result.txt"),sep="\t",header=T,row.names=1)
  label=read.table(paste0(wd2,"TCGA-",cs[j],".txt"),sep="\t",header=T)
  number_cancer=gsub("-",".",label[which(label[,2]=="Tumor"),1])
  number_normal=gsub("-",".",label[which(label[,2]=="Normal"),1])
  cancer=data[cell,number_cancer]
  normal=data[cell,number_normal]
  write.table(cancer,paste0(wd3,"TCGA_",cs[j],"_tumor_sample_TILs.txt"),sep="\t",quote=F,row.names=T)
  write.table(normal,paste0(wd3,"TCGA_",cs[j],"_normal_sample_TILs.txt"),sep="\t",quote=F,row.names=T)
}

##画图
setwd("C://Users/Administrator/Desktop/心慧/Sample_immunefiltratuion_group/")
cs=c("BRCA","CHOL","COAD","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")
cell=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg")
for(i in 1:14){
  list1<-list() 
  list2<-list()
  M=c()
  MM=c()
  for(j in 1:length(cs)){
    Data=read.table(paste0("TCGA_",cs[j],"_tumor_sample_TILs.txt"),sep="\t",header=T,row.names=1) 
    Data=t(Data)
    list1[[j]]=as.data.frame(t(as.data.frame(Data[,i])))
    m=median(Data[,i])
    M=rbind(M,m) 
    
    Data1=read.table(paste0("TCGA_",cs[j],"_normal_sample_TILs.txt"),sep="\t",header=T,row.names=1) 
    Data1=t(Data1)
    list2[[j]]=as.data.frame(t(as.data.frame(Data1[,i])))
    mm=median(Data1[,i])
    MM=rbind(MM,mm)
    
  }
  A=do.call(rbind.fill,list1)
  rownames(A)=cs
  B=arrange(A,M)  
  cs_new=arrange(as.data.frame(cs),M)
  
  AA=do.call(rbind.fill,list2)
  rownames(AA)=cs
  BB=arrange(AA,MM)  
  cs_new_1=arrange(as.data.frame(cs),MM)
  
  pdf(paste("C://Users/Administrator/Desktop/心慧/Sample_immunefiltratuion_group/",cell[i],"_14_cancer(normal_cancer)_boxplot.pdf",sep=''))
  boxplot(t(B),outline=F,names=t(cs_new),boxwex=0.25,col="#e08080",at=1:14-0.15,xaxt="n",yaxt="n")
  boxplot(t(BB),outline=F,names=t(cs_new_1),boxwex=0.25,add=T,at=1:14+0.15,las=2)
  dev.off()
}




#############################################################
*****************************************************
  ###修稿中要加的图 
  ***SS1
  *********1.加在fig2后面，附图，回答的是Reviewer #2: 的第三个问题
##计算sig在免疫细胞中的方差值

setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/加-改的附图/SS1/")
lnc=read.table("lncRNA_log2exp_all_sample.txt",header=TRUE,sep="\t",row.names=1)##没有进行log转换
coldata=read.table("colData.txt",sep="\t",header=TRUE,row.names=1)
cell=unique(coldata[,1])
lnc_cell=read.table("lncRNA_log2exp_14cell.txt",sep="\t",header=T,row.names=1)
uu=c()
for(i in 1:14){
  sig=read.table(paste0(cell[i],"_log2_annova_FDR_0.3.txt"),sep="\t",header=T,row.names=1)
  uu=c(uu,rownames(sig))
}
other_res=c()
union_res=unique(uu)
other=rownames(lnc_cell)[-which(rownames(lnc_cell)%in%union_res)]


RES_var=c()
RES_var_other=c()
annova_P_RES=c()
other_annova_P_RES=c()
for(i in 1:14){
  sig=read.table(paste0(cell[i],"_log2_annova_FDR_0.3.txt"),sep="\t",header=T,row.names=1)
  sig_group=lnc[rownames(sig),] 
  var_res=apply(sig_group,1,var)
  var_res1=cbind(rep(cell[i],length(sig[,1])),var_res)
  RES_var=rbind(RES_var,var_res1)
  other_lnc=lnc[sample(other,length(sig[,1])),]
  var_other=apply(other_lnc,1,var)
  var_other1=cbind(rep(paste0("NOT",cell[i]),length(other_lnc[,1])),var_other)
  RES_var_other=rbind(RES_var_other,var_other1)
  
  for (j in 1:length(sig[,1])) {
    annova_test=anova(lm(as.numeric(sig_group[j,])~coldata[,1]))
    annova_P=annova_test[5][[1]][1]
    annova_P_res=cbind(cell[i],annova_P)
    annova_P_RES=rbind(annova_P_RES,annova_P_res)
    
    other_annova_test=anova(lm(as.numeric(other_lnc[j,])~coldata[,1]))
    other_annova_P=other_annova_test[5][[1]][1]
    other_annova_P_res=cbind(paste0("NOT",cell[i]),other_annova_P)
    other_annova_P_RES=rbind(other_annova_P_RES,other_annova_P_res)
  }
}
annova_P_RES[,2]=p.adjust(annova_P_RES[,2],method="fdr")
other_annova_P_RES[,2]=p.adjust(other_annova_P_RES[,2],method="fdr")
write.table(RES_var,"sig_var_res.txt",sep="\t",quote=F,col.names = F,row.names = F)
write.table(RES_var_other,"Not_sig_var_res.txt",sep="\t",quote=F,col.names = F,row.names = F)
write.table(annova_P_RES,"sig_annova_P_RES.txt",sep="\t",quote=F,col.names = F,row.names = F)
write.table(other_annova_P_RES,"Not_sig_annova_P_RES.txt",sep="\t",quote=F,col.names = F,row.names = F)


RES=read.table("sig_var_res.txt",sep="\t",header=F)
RES1=read.table("Not_sig_var_res.txt",sep="\t",header=F)
RES_all=rbind(RES,RES1)
colnames(RES_all)=c("Cell","Var")
ggplot(RES_all,aes(x=factor(RES_all[,1],level=c(cell[1],paste0("NOT",cell[1]),
      cell[2],paste0("NOT",cell[2]),cell[3],paste0("NOT",cell[3]),cell[4],
      paste0("NOT",cell[4]),cell[5],paste0("NOT",cell[5]),cell[6],
      paste0("NOT",cell[6]),cell[7],paste0("NOT",cell[7]),cell[8],
      paste0("NOT",cell[8]),cell[9],paste0("NOT",cell[9]),cell[10],
      paste0("NOT",cell[10]),cell[11],paste0("NOT",cell[11]),cell[12],
      paste0("NOT",cell[12]),cell[13],paste0("NOT",cell[13]),cell[14],
      paste0("NOT",cell[14]))),y=RES_all[,2],fill=Cell))+geom_boxplot(outlier.colour = NA)+
      coord_cartesian(ylim=c(0,55))+
      labs(title="Var of signature and No-signature",x="Cell", y = "Var")+theme_bw()
      ##+geom_jitter(shape=16, position=position_jitter(0.2))##添加点

RES=read.table("sig_annova_P_RES.txt",sep="\t",header=F)
RES1=read.table("Not_sig_annova_P_RES.txt",sep="\t",header=F)
RES_all=rbind(RES,RES1)
colnames(RES_all)=c("Cell","VarP")
ggplot(RES_all,aes(x=factor(RES_all[,1],level=c(cell[1],paste0("NOT",cell[1]),
        cell[2],paste0("NOT",cell[2]),cell[3],paste0("NOT",cell[3]),cell[4],
        paste0("NOT",cell[4]),cell[5],paste0("NOT",cell[5]),cell[6],
        paste0("NOT",cell[6]),cell[7],paste0("NOT",cell[7]),cell[8],
        paste0("NOT",cell[8]),cell[9],paste0("NOT",cell[9]),cell[10],
        paste0("NOT",cell[10]),cell[11],paste0("NOT",cell[11]),cell[12],
        paste0("NOT",cell[12]),cell[13],paste0("NOT",cell[13]),cell[14],
        paste0("NOT",cell[14]))),y=-log10(VarP),fill=Cell))+geom_boxplot(outlier.colour = NA)+
  labs(title="annova of signature and No-signature",x="Cell", y = "annova_p")+
  coord_cartesian(ylim=c(0,25))+theme_bw()
##+geom_jitter(shape=16, position=position_jitter(0.2))##添加点


*************************************
  ******SS2
#回答的是Reviewer #3:4.
###首先加了一个LNCRNACIBER和CIBERSORT结果的相关性，合到FigS8里面，程序在S8那里
###对于FIg5种的结果，加上一些xcell和CIBERSORT自己的预测结果
###雷达图，癌症样本和正常样本的差异性，预测值对预后的影响
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/加-改的附图/CIBERSORT_xcell/")
cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
CIB=read.table("Cibersort.txt",sep="\t",header=T,row.names=1)
xcell=read.table("xcell.txt",sep="\t",header=T,row.names = 1)
for(i in 1:33){
  CIB1=CIB[which(CIB[,8]==cs[i]),]
  a=CIB1[which(CIB1[,11]=="Tumor"),1:7]
  CIB11=apply(CIB1[which(CIB1[,11]=="Tumor"),1:7],2,mean)
  xcell1=xcell[which(xcell[,12]==cs[i]),]
  xcell11=apply(xcell1[which(xcell1[,15]=="Tumor"),1:11],2,mean)
  write.table(t(CIB11),"CIBERSORT_result_6cell_mean.txt",sep="\t",quote=F,append=T,row.names=F,col.names=F)
  write.table(t(xcell11),"Xcell_result_11cell_mean.txt",sep="\t",quote=F,append=T,row.names=F,col.names=F)
}
#write.table(cs,"cs.txt",sep="\t",quote=F,row.names=F,col.names=F)

res=read.table("CIBERSORT_result_6cell_mean.txt",sep="\t",header=T,row.names = 1)
library(fmsb)
maxmin=matrix(rep(c(0.5,0),each=33),byrow=T,ncol=33)
colnames(maxmin)=rownames(res)
Data=rbind(maxmin,t(res))
Data=as.data.frame(Data)
rownames(Data)

radarchart(Data, seg=6,plty=1,pty=32,plwd=3,
           cglty=2,pcol=c("M0"="Goldenrod1","M1"="Sienna3","M2"="Sienna4",
          "Monocyte"="PaleGreen","Neutrophil"="Plum2",
          "CD8.T"="DarkOliveGreen3","Treg"="black"))
text(2,c(0.10,0.2,0.3,0.4,0.5,0.6,0.7),c("—— M0","—— M1","—— M2",
    "—— Monocyte","—— Neutrophil","—— CD8 T","—— Treg"),
     col=c("Goldenrod1","Sienna3","Sienna4","PaleGreen",
           "Plum2","DarkOliveGreen3","black"))
max(res)
#[1] 0.3849654,CIBERSORT_max_0.39_figA.pdf



res=read.table("Xcell_result_11cell_mean.txt",sep="\t",header=T,row.names = 1)
library(fmsb)
maxmin=matrix(rep(c(0.3,0),each=33),byrow=T,ncol=33)
colnames(maxmin)=rownames(res)
Data=rbind(maxmin,t(res))
Data=as.data.frame(Data)
rownames(Data)

radarchart(Data, seg=6,plty=1,pty=32,plwd=3,cglty=2,pcol=c("B"="PaleTurquoise1","DC"="red","M1"="Sienna3","M2"="Sienna4","Macrophage"="Burlywood1",
                                                           "Monocyte"="PaleGreen","Neutrophil"="Plum2","NK"="RosyBrown1","CD8.T"="DarkOliveGreen3",
                                                           "CD4.T"="MediumPurple2","Treg"="black"))
text(2,c(0.10,0.12,0.2,0.28,0.36,0.44,0.52,0.60,0.68,0.76,0.84,0.92,0.96,1),c("—— B","—— DC","—— M1","—— M2","—— Macrophage","—— Monocyte",
                                                                              "—— Neutrophil","—— Nk","—— CD8 T","—— CD4 T","—— Treg"),
     col=c("PaleTurquoise1","red","Sienna3","Sienna4","Burlywood1","PaleGreen","Plum2","RosyBrown1","DarkOliveGreen3","MediumPurple2","black"))

max(res)

###tumor/normal差异检验气泡图

setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5D/修稿修改/")
#cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
cs=c("BRCA","CHOL","COAD","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")

#res2=c()
cell=c("B","DC","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T")
for(i in 1:14){
  wd1="F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5D/修稿修改/cnacer_LNCCIBERSORT_results/"
  wd2="F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5D/修稿修改/33_cancer_label/"
  data=read.table(paste0(wd1,cs[i],"_14cell_pre_CIBEERSORT_result.txt"),sep="\t",header=T,row.names=1)
  label=read.table(paste0(wd2,"TCGA-",cs[i],".txt"),sep="\t",header=T)
  number_cancer=gsub("-",".",label[which(label[,2]=="Tumor"),1])
  number_normal=gsub("-",".",label[which(label[,2]=="Normal"),1])
  cancer=data[cell,number_cancer]
  normal=data[cell,number_normal]
  
  
  for(j in 1:12){
    res=c()
    TT=t.test(as.numeric(cancer[j,]),as.numeric(normal[j,]))
    TP=TT$p.value
    WT=wilcox.test(as.numeric(cancer[j,]),as.numeric(normal[j,]))
    WP=WT$p.value
    FC=mean(as.numeric(cancer[j,]))/mean(as.numeric(normal[j,]))
    res=cbind(cs[i],cell[j],TP,FC)
    #res=cbind(cs[i],cell[j],WP,FC)
    write.table(res,"14cancer_T_P_FC.txt",sep="\t",row.names=F,col.names=F,quote=F,append=T)
  }
}


###改进res111.pdf
library(ggplot2)
res=read.table("14cancer_T_P_FC1.txt",sep="\t",header=T)
res[which(res[,4]>2),4]=2
res[intersect(which(res[,4]>1),which(res[,4]<2)),4]=1.5
res[which(res[,4]<0.5),4]=0.3
res[intersect(which(res[,4]>0.5),which(res[,4]<1)),4]=0.8
ggplot(as.data.frame(res),aes(cell,Cancer))+geom_point(aes(size=-log(Pvalue,10),color=as.factor(FC)))+theme_bw()



#########################################################3
###Xcell和Cibersort的生存分析
##xcell
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/加-改的附图/CIBERSORT_xcell/更新Fig6部分/")
data1=read.table("xcell_all.txt",sep="\t",header=T,row.names = 1)
data2=read.table("xcell_row_cluster_K6.txt",sep="\t",header=T,row.names = 1)
data1=data1[which(data1[,68]=="Tumor"),65:67]
data=cbind(data1,data2[,65])
colnames(data)=c("Cancer_lable","OS","OS.time","Cluster")
cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
y=Surv(data$OS.time,data$OS)
kmfit2<-survfit(y~data$Cluster, data=data) 
plot(kmfit2, main="Kaplan-Meier curves of six groups", xlab="Time(days)",  ylab="Overall survival", mark=3,mark.time=T,col=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2"))
text(10000,c(0.8,0.84,0.88,0.92,0.96,1),c("—— 1","—— 2","—— 3","—— 4","—— 5","—— 6"),col=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2"))
text(6000,0.8,paste0("p=",p) )
##计算P值
r<-survdiff(formula = y~data$Cluster,  rho = 0,data=data) 
p<-1 - pchisq(r$chisq, length(r$n) - 1)

###分别画33个癌症的生存结果，去掉样本数少于10个的类
color=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2")
for(i in 1:33){
  data1=data[which(data[,1]==cs[i]),]
  a=unique(data1[,4])
  aa=c()
  for(n in 1:length(a)){
    cc=length(which(data1[,4]==a[n]))
    if(cc>=10){
      aa=c(aa,a[n])
    }
  }
  if(length(aa)>1){
    data1=data1[which(data1[,4]%in%aa),]
    y=Surv(data1$OS.time,data1$OS)
    kmfit2<-survfit(y~data1$Cluster, data=data1) 
    main_name<-paste("module_",cs[i],sep="")
    data.survdiff=survdiff(y~data1$Cluster,data = data1)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    if(p.val<0.05){
    pdf(paste(main_name,"_xcell_survival.pdf",sep=''))
    plot(kmfit2,main=main_name,xlab="survival time",ylab="survival probability",mark=3,mark.time=T,cex.lab=1.5,cex.main=1.5,col=color[aa],lwd=5)
    legend("bottomright",col=color[aa],legend=aa,lty=1,lwd=3)
    legend("bottomleft",legend=paste0("P=",p.val))
    dev.off()
  }
}
}

###Cibersort

setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/加-改的附图/CIBERSORT_xcell/更新Fig6部分/")
data1=read.table("Cibersort_all.txt",sep="\t",header=T,row.names = 1)
data2=read.table("Cibersort_row_cluster_K6.txt",sep="\t",header=T,row.names = 1)
data1=data1[which(data1[,26]=="Tumor"),23:25]
data=cbind(data1,data2[,23])
colnames(data)=c("Cancer_lable","OS","OS.time","Cluster")
y=Surv(data$OS.time,data$OS)
kmfit2<-survfit(y~data$Cluster, data=data) 
plot(kmfit2, main="Kaplan-Meier curves of six groups", xlab="Time(days)",  ylab="Overall survival", mark=3,mark.time=T,col=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2"))
text(10000,c(0.8,0.84,0.88,0.92,0.96,1),c("—— 1","—— 2","—— 3","—— 4","—— 5","—— 6"),col=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2"))
text(6000,0.8,paste0("p=",p) )
table(data$Cluster)
##计算P值
r<-survdiff(formula = y~data$Cluster,  rho = 0,data=data) 
p<-1 - pchisq(r$chisq, length(r$n) - 1)

###分别画33个癌症的生存结果，去掉样本数少于10个的类
color=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2")
for(i in 1:33){
  data1=data[which(data[,1]==cs[i]),]
  a=unique(data1[,4])
  aa=c()
  for(n in 1:length(a)){
    cc=length(which(data1[,4]==a[n]))
    if(cc>=10){
      aa=c(aa,a[n])
    }
  }
  if(length(aa)>1){
    data1=data1[which(data1[,4]%in%aa),]
    y=Surv(data1$OS.time,data1$OS)
    kmfit2<-survfit(y~data1$Cluster, data=data1) 
    main_name<-paste("module_",cs[i],sep="")
    data.survdiff=survdiff(y~data1$Cluster,data = data1)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    if(p.val<0.05){
    pdf(paste(main_name,"_Cibersort_survival.pdf",sep=''))
    plot(kmfit2,main=main_name,xlab="survival time",ylab="survival probability",mark=3,mark.time=T,cex.lab=1.5,cex.main=1.5,col=color[aa],lwd=5)
    legend("bottomright",col=color[aa],legend=aa,lty=1,lwd=3)
    legend("bottomleft",legend=paste0("P=",p.val))
    dev.off()
  }
}
}






#############计算lncRNA预测结果的免疫亚型聚类的类和Cibersort/xcell的聚类结果的交集
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/加-改的附图/CIBERSORT_xcell/更新Fig6部分/和lnc计算交集/")
xcell=read.table("xcell_row_cluster_K6.txt",sep="\t",header=T,row.names = 1)
ciber=read.table("cibersort_row_cluster_K6.txt",sep="\t",header=T,row.names = 1)
lnc=read.table("Tumorsample_Sur_lncRNACibersort_immunefiltration_merge.txt",sep="\t",header=T,row.names = 1)
rownames(xcell)=gsub("-",".",rownames(xcell))
rownames(ciber)=gsub("-",".",rownames(ciber))
inter1=matrix(0,ncol=6,nrow=6)
inter2=matrix(0,ncol=6,nrow=6)
for(i in 1:6){
  for(j in 1:6){
  lnc1=lnc[which(lnc[,18]==i),]
  ciber1=ciber[which(ciber[,23]==j),]
  xcell1=xcell[which(xcell[,65]==j),]
  inter1[i,j]=length(intersect(rownames(lnc1),rownames(xcell1)))
  inter2[i,j]=length(intersect(rownames(lnc1),rownames(ciber1)))
  }
}
table(lnc[,18])
table(xcell[,65])
table(ciber[,23])

>table(lnc[,18])

1    2    3    4    5    6 
3210 4161 1778  535  207  230 
> table(xcell[,65])

1    2    3    4    5    6 
9509   93  245  457   20    3 
> table(ciber[,23])

1    2    3    4    5    6 
2162 4614 1167 1153  974  257 

lnc_ciber=lnc[rownames(ciber),]
lnc_xcell=lnc[rownames(ciber),]
ciber[,24]=lnc_ciber[,18]
xcell[,66]=lnc_xcell[,18]


exp_data=as.matrix(data[,1:22])
#png("cancers_ciber_all33_new.png")
cluster=ciber[,23]
annotation_col=data.frame(Cluster_CIBERSORT=factor(as.numeric(cluster)),Cluster_LNCRNA=factor(as.numeric(ciber[,24])))

annotation_colors =list(Cluster_CIBERSORT=c("1"="red","2"="MediumPurple2","3"="PaleGreen","4"="Plum2","5"="DarkOliveGreen3","6"="Gold2"),
                        Cluster_LNCRNA=c("1"="red","2"="MediumPurple2","3"="PaleGreen","4"="Plum2","5"="DarkOliveGreen3","6"="Gold2"))
rownames(exp_data)=rownames(annotation_col)
pheatmap(exp_data, show_colnames= T, show_rownames= F, scale= "none", fontsize= 6.5,
         clustering_method ="complete",
         cluster_cols=F,cluster_rows=T,
         annotation_row= annotation_col, 
         annotation_colors= annotation_colors, 
         breaks = c(seq(0,0.05,length.out=50),seq(0.051,0.9,length.out=50)),
         col = colorRampPalette(c("navy", "white", "firebrick3"))(100))

exp_data=as.matrix(xcell[,1:64])
#png("cancers_ciber_all33_new.png")
cluster=xcell[,65]
annotation_col=data.frame(Cluster_xcell=factor(as.numeric(cluster)),Cluster_LNCRNA=factor(as.numeric(xcell[,66])))

annotation_colors =list(Cluster_xcell=c("1"="red","2"="MediumPurple2","3"="PaleGreen","4"="Plum2","5"="DarkOliveGreen3","6"="Gold2"),
                        Cluster_LNCRNA=c("1"="red","2"="MediumPurple2","3"="PaleGreen","4"="Plum2","5"="DarkOliveGreen3","6"="Gold2"))
rownames(exp_data)=rownames(annotation_col)
pheatmap(exp_data, show_colnames= T, show_rownames= F, scale= "none", fontsize= 6.5,
         clustering_method ="complete",
         cluster_cols=F,cluster_rows=T,
         annotation_row= annotation_col, 
         annotation_colors= annotation_colors, 
         breaks = c(seq(0,0.05,length.out=50),seq(0.051,1,length.out=50)),
         col = colorRampPalette(c("navy", "white", "firebrick3"))(100))


*************************************************************************************
  






####################################################################################3
############Fig S9
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/加-改的附图/Fig S9/修稿修改/")
A=read.table("TT_lnc_exp_new.txt",sep="\t",header=T,row.names=1)
#A=log2(A+0.001)
#write.table(A,"TT_lnc_exp_new_log.txt",sep="\t",quote=F)
source('CIBERSORT.R')
cell=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","Treg","CD3 T","CD4 T","CD8 T")
for(i in 1:14){
result=CIBERSORT("TT_lnc_exp_new_log.txt",paste0(cell[i],"_log2_annova_FDR_0.3.txt"),0,TRUE)
write.table(result,paste0(cell[i],"LGG_lncRNA_cibersort_result.txt"),sep="\t",quote=F,row.names=T,col.names=T)
}
###手动整理了结果文件



##画箱式图
res=read.table("LGG_lncCiber_result.txt",sep="\t",header=T,row.names=1)
res=t(res)
LGG_grade<-read.table("TT_sample.txt",header=T,sep="\t",row.names=1,fill=T)

###把肿瘤分级数据加到最后一列
inter=intersect(rownames(res),rownames(LGG_grade))
res1=res[inter,]
grade1=LGG_grade[inter,5]
data=cbind(res1,grade1)
#数据前14列是预测结果，最后一列是癌症恶性程度
library(cowplot)
BiocManager::install("tidyverse")
library(tidyverse)
library(ggplot2)
BiocManager::install("ggsci")
library(ggsci)
BiocManager::install("ggpubr")
library(ggpubr)

cell=c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","Treg","CD3 T","CD4 T","CD8 T")
for(i in 1:14){
  F=aov(data[,i]~data[,15])
  P=summary(F)
  PP=unlist(P)
  print(c(cell[i],PP[9][[1]]))
}
class(data)
data=as.data.frame(data)
colnames(data)
data("ToothGrowth")
df <- ToothGrowth
for(i in 1:14){
my_comparisons <- list( c("2", "3"),c("3","4"), c("2", "4"))
ggboxplot(data,x="grade1",y=cell[i], color = "grade1",palette = "jco", add = "jitter",xlab="Grade")+
          stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test")
#画箱式图并根据T检验得出P值进行标注显著性
}



************************************************************************
#########Fig S11
####计算几个免疫检查点基因在不同分类中的差异性

setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/加-改的附图/Fig S11/修稿修改/")
cs=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
gene=c("CD80","CD86","CGAS","CTLA4","PD-1","STING")

ID=c("ENSG00000121594","ENSG00000114013","ENSG00000164430","ENSG00000163599","ENSG00000188389","ENSG00000184584")

  for(i in 1:33){
    Cancer=read.table(paste0("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/结果五/Fig 5A/cancer-new/TCGA-",cs[i],"-htseq_fpkm-uq.txt"),sep="\t",header=T,row.names=1)
    Cancer_res=t(Cancer[ID,])
    write.table(Cancer_res,"6gene_33cancer_exp.txt",quote=F,sep="\t",append=T)
  }

pre=read.table("Tumorsample_Sur_lncRNACibersort_immunefiltration_merge.txt",sep="\t",header=T,row.names = 1)
geneexp=read.table("6gene_33cancer_exp.txt",sep="\t",header=T,row.names=1)
inter=intersect(rownames(pre),rownames(geneexp))
res=cbind(pre[inter,18],geneexp[inter,])
res[1,]
boxplot(res[,7]~res[,1],outline=FALSE,main="STING",xlab="Cluster",ylab="Expression",col=c("red","MediumPurple2","PaleGreen","Plum2","DarkOliveGreen3","Gold2"))
text(6,22,"P<2e-16 ***")
##方差分析
F=aov(res[,7]~res[,1])
F_p=summary(F)
F_p











#*****************************************************************8
#####创建R
包
###整理需要得数据
#lncRNA_signature_matrix
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/其他需要加的材料/LNCRNACIBERSORT/")
lnc=read.table("lncRNA_log2exp_14cell.txt",sep="\t",header=T,row.names=1)
cell<- c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg")
res=c()
for(i in 1:14){
sig=read.table(paste0(cell[i],"_log2_annova_FDR_0.3.txt"),sep="\t",header=T,row.names=1)
res=c(res,rownames(sig))
}
RES=lnc[unique(res),]
write.table(RES,"lncRNA_signature_matrix.txt",sep="\t",quote=F)

#
setwd("F:/修稿文件汇总/改稿进展/数据_结果_程序_整理/其他需要加的材料/LNCRNACIBERSORT/")
lnc=read.table("lncRNA_signature_matrix.txt",sep="\t",header=T,row.names=1)
cell<- c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg")
res=matrix(0,ncol=14,nrow=length(RES[,1]))
rownames(res)=rownames(RES)
colnames(res)=cell
for(i in 1:14){
  sig=read.table(paste0(cell[i],"_log2_annova_FDR_0.3.txt"),sep="\t",header=T,row.names=1)
  res[rownames(sig),i]=1
}

write.table(res,"lncRNA_signature_label.txt",sep="\t",quote=F)

length(sig[,1])
sum(res[,14])
source("LncRNACIBERSORT.R")
y=ce
perm=0
x="lncRNA_log2exp_14cell.txt"
LncRNACIBERSORT <- function(x,y,perm=0){
  lnc_sig=read.table("lncRNA_signature_matrix.txt",sep="\t",header=T,row.names=1)
  sig_label=read.table("lncRNA_signature_label.txt",sep="\t",header=T,row.names=1)
  cell<- c("B","DC","Granulocytes","M0","M1","M2","Macrophage","Monocyte","Neutrophil","NK","CD8 T","CD4 T","CD3 T","Treg")
  cellnumber <- which(cell %in% y)
  for(n in cellnumber){
    lnc_sig_cell<-lnc_sig[rownames(sig_label)[which(sig_label[,n]==1)],]
    write.table(lnc_sig_cell,paste0(cell[n],"lncRNA_sig.txt"),sep="\t",quote=F)
    pre <- CIBERSORT(paste0(cell[n],"lncRNA_sig.txt"),x, perm, TRUE)
    #return matrix object containing all results
    precell <- rbind(precell,pre[cell[n],])
    
  }
  #save results
  rownames(precell)=cell[cellnumber]
  write.table(precell,file=paste0("LNCRNACIBERSORT_",x),sep="\t",quote=F,append=T,row.names=T,col.names=T)
  precell
}
source("LncRNACIBERSORT.R")
ce=c("B","DC")
result=LncRNACIBERSORT("lncRNA_log2exp_14cell.txt",ce,0)











