library(WGCNA)
library(Seurat)
library(tidyverse)
library(reshape2)
library(stringr)
# 多线程聚类
library(fastcluster)
# 多线程计算距离
library(parallelDist)
library(org.Hs.eg.db)
library(clusterProfiler)
load("../AST_PFC_wgcna.rdata")
df = AggregateExpression(object = AST_PFC, group.by = "Individual")
df=as.data.frame(df)
cpm <- apply(df ,2, function(x) { x/sum(x)*1000000 })
df = log(cpm+1)
df = as.data.frame(df)
allDEGs = read.csv("../AllUpDown.csv")
uniqueDegs = read.csv("../uniqueDegs.csv",header = 1)
# 所有的
x = df[uniqueDegs$degs,]
x = t(x)
# ASD的
# x = df[allDEGs$ASD,]
# x = na.omit(x)
# x = t(x)
powers = c(1:50)
sft = pickSoftThreshold(x, powerVector = powers, verbose = 3)
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

cor=WGCNA::cor
net = blockwiseModules(x, power = 30,
                       TOMType = "unsigned", minModuleSize = 10,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       verbose = 5)
moduleColors = labels2colors(net$colors)
table(moduleColors)

plotDendroAndColors(net$dendrograms[[1]],moduleColors[net$blockGenes[[1]]],"Module colors",dendroLabels =F ,hang =0.03,addGuide=T,guideHang=0.05)


names(x)=colnames(x)
adjacency = adjacency(x, power = 30)
TOM = TOMsimilarity(adjacency)

# 选出要的模块的TOM矩阵
probes = names(x)
modules = "turquoise"
inModule = moduleColors=="turquoise" | moduleColors=="blue"
modProbes = probes[inModule]
modTOM = TOM[inModule, inModule]

# 导出网络 前500条边
cyt = exportNetworkToCytoscape(modTOM,edgeFile = "edge_blue_turquoise.txt", nodeFile = "node_blue_turquoise.txt",weighted = TRUE, threshold = 0.627067069909246,nodeNames = modProbes, nodeAttr = moduleColors[inModule])


# 找hub基因
module = "turquoise"
MEList = moduleEigengenes(x, colors = moduleColors)
# KME值接近0,说明这个基因不是该模块的成员：KME接近1或者－1,说明这个基因与该模块密切相关（正相关或者负相关）
MEs = MEList$eigengenes
# 基因和模块的kme
datKME=signedKME(x, MEs, outputColumnName="kME_MM.")
# modNames就是MEs的列名只保留颜色,从kmeblue变成blue
modNames = substring(names(MEs), 3)
column = match(module,modNames)
moduleGenes = moduleColors==module
turquoise_module = as.data.frame(dimnames(data.frame(x))[[2]][moduleGenes])
names(turquoise_module)="genename"

turquoise_KME = as.data.frame(datKME[moduleGenes,column]) 
names(turquoise_KME)="KME"
rownames(turquoise_KME)=turquoise_module$genename
# 自己指定阈值
FilterGenes = abs(turquoise_KME$KME) > 0.99
turquoise_hub = turquoise_module[FilterGenes,"genename"]



module = "blue"
MEList = moduleEigengenes(x, colors = moduleColors)
# KME值接近0,说明这个基因不是该模块的成员：KME接近1或者－1,说明这个基因与该模块密切相关（正相关或者负相关）
MEs = MEList$eigengenes
# 基因和模块的kme
datKME=signedKME(x, MEs, outputColumnName="kME_MM.")
# modNames就是MEs的列名只保留颜色,从kmeblue变成blue
modNames = substring(names(MEs), 3)
column = match(module,modNames)
moduleGenes = moduleColors==module
blue_module = as.data.frame(dimnames(data.frame(x))[[2]][moduleGenes])
names(blue_module)="genename"

blue_KME = as.data.frame(datKME[moduleGenes,column]) 
names(blue_KME)="KME"
rownames(blue_KME)=blue_module$genename
# 自己指定阈值
FilterGenes = abs(blue_KME$KME) > 0.99
blue_hub = blue_module[FilterGenes,"genename"]


