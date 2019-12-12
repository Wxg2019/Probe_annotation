# Probe_annotation
Simplified methods of probe annotaiton, referred from Jimmy(jmzeng1314/AnnoProbe)【https://github.com/jmzeng1314/annoprobe】. Here just the notes from that.


#有趣的是，因为这些包存储在GitHub，而且每个包自带的数据是40~50M，对很多在中国大陆的朋友来说， 
#几乎是不可能完成，所以我把这4个包整合成为了一个GitHub包（AnnoProbe）！总共不到5M，相信大家使用
#起来应该是很方便啦！


library(devtools)
install_github("jmzeng1314/AnnoProbe")
library(GEOmirror)

rm(list = ls())
library(AnnoProbe) 
suppressPackageStartupMessages(library(GEOquery)) 
gset=AnnoProbe::geoChina('GSE1009')
gset


# check the ExpressionSet
eSet=gset[[1]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)
head(probes_expr[,1:4])
boxplot(probes_expr,las=2)
probes_expr=log2(probes_expr+1)
boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])

#然后对表达芯片的探针进行基因注释
#这一个步骤也是在线下载我们的芯片注释信息
## check GPL and annotate the probes to genes.
(gpl=eSet@annotation)
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene=idmap(gpl)
head(probe2gene)
genes_expr <- filterEM(probes_expr,probe2gene )
head(genes_expr)


#走limma的经典2组差异分析
# do DEG
## define the group
group_list=factor(c(rep('Control',3),rep('Diabetes',3)))
table(group_list)
library(limma)
design=model.matrix(~factor(group_list))
design
fit=lmFit(genes_expr,design)
fit=eBayes(fit)
DEG=topTable(fit,coef=2,n=Inf)
head(DEG)


#对差异分析结果进行一些检验
## visualization
need_deg=data.frame(symbols=rownames(DEG), logFC=DEG$logFC, p=DEG$P.Value)
deg_volcano(need_deg,1)
deg_volcano(need_deg,2)

deg_heatmap(DEG,genes_expr,group_list)
deg_heatmap(DEG,genes_expr,group_list,30)

check_diff_genes('PLCE1',genes_expr,group_list)
check_diff_genes('MPP6',genes_expr,group_list)


#如果你做了GO/KEGG注释后也可以挑选基因集进行可视化
# 假设我这里对hsa03410感兴趣
library(KEGGREST)
cg <- KEGGREST::keggGet("hsa03410")[[1]]$GENE
cg=as.character(sapply(cg[seq(2,length(cg),by=2)], function(x) strsplit(x,';')[[1]][1]))
check_diff_genes( cg ,genes_expr,group_list)



#上面的代码你可以套用到任何一个表达芯片数据集
#当然了，你需要有一点R语言基础知识啦，不然，上面的代码你不知道应该修改哪里。
#马上试试看下面的数据集吧，我觉得蛮有意义的。 GSE1462/GSE20950/GSE21785/GSE26526/GSE32575
#GSE43837/GSE474/GSE58979/GSE60291/GSE62832/GSE70529/GSE72158


##三个万能芯片探针注释


#第一个万能芯片探针ID注释平台R包,from jimmy
rm(list=ls())

library(devtools)
#install_github("jmzeng1314/idmap1")
library(idmap1)

#获取同样的GPL570    hgu133plus2 [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array注释信息，一行代码就搞定！

ids=getIDs('gpl570')
head(ids)


#实战演练
#结合我们昨天的发布的 GEO数据库中国区镜像横空出世，随时随地方便下载GEO数据集，并且进行ID转换！
rm(list=ls())
#install_github("jmzeng1314/GEOmirror")
#source('http://raw.githubusercontent.com/jmzeng1314/GEOmirror/master/R/geoChina.R') 
library(devtools)
library(GEOmirror)
library(idmap1) 
geoChina('gse6222')
load('GSE6222_eSet.Rdata')
gset
a=exprs(gset[[1]])
a[1:4,1:4]
gset[[1]]@annotation
b=getIDs("GPL570")
head(b)

#还有一个小功能
#就是对基因名字进行注释
IDs <- c("DDX11L1", "MIR6859-1", "OR4G4P", "OR4F5")
ID_type = "SYMBOL"
annoGene(IDs, ID_type)

#同样的，我们前面注释好的探针基因对应关系，也可以进行进一步注释：
ids=getIDs('gpl570')
head(ids)
annoGene(ids$symbol,'SYMBOL')



#第二个万能芯片探针ID注释平台R包
#一定要跟我们的 idmap1 区分开来哦，那个idmap1是把bioconductor里面有的37个芯片平台整合
#了一下，而我们的这个idmap2包不得了啦，有122个GPL之多，它们的6G多的soft信息被我下载整
#理成为了不到40M的R包
rm(list=ls())
library(devtools)  
install_github("jmzeng1314/idmap2")  
#安装失败，另想办法，
#如下方：https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247491955&idx=1&sn=1ab7d8dfeeca67221aff28e84df2b4d5&scene=21#wechat_redirect
#devtools::install("jmzeng1314/idmap2") #该方法尝试失败
library(idmap2)

#同样的获取同样的GPL570    hgu133plus2 [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array注释信息，一行代码就搞定！
library(idmap1)
ids=getIDs('gpl570')
head(ids)
library(idmap2)
ids=get_soft_IDs('gpl570')
head(ids)

#当然了，第二个包有一百多个平台，而第一个只有37个，注意哦，比如
idmap1::getIDs('GPL13912') # 失败
idmap2::get_soft_IDs('GPL13912') # 成功

#这个 GPL13912 平台，就存在于第二个包，但是不在第一个哈。
#你想知道我们支持哪些平台吗，当然是可以看的：
data(gpl_list)
gpl_list[,1:4]

#一个芯片数据挖掘实战
#结合我们发布的 GEO数据库中国区镜像横空出世，随时随地方便下载GEO数据集，并且进行ID转换！
library(GEOmirror)
library(idmap1) 
library(idmap2) 
gset=geoChina('GSE31731') 
gset
a=exprs(gset[[1]])
a[1:4,1:4]
gpl=gset[[1]]@annotation
b=idmap2::get_soft_IDs(gpl)
head(b)


#第三个万能芯片探针ID注释平台R包
#idmap1解决了bioconductor包下载困难的问题，idmap2解决了GPL平台的soft文件下载困难，而这个idmap3解决了那些并不提供探针的注释信息的平台。
library(devtools)
install_github("jmzeng1314/idmap3")
library(idmap3)

library(idmap3)
ids=idmap3::get_pipe_IDs('GPL21827')
head(ids) 

#你想知道我们支持哪些平台吗，当然是可以看的：
data(gpl_list)
gpl_list[,1:4]
