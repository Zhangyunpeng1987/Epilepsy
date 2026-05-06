###查看海马体兴奋性神经元与颞叶皮层兴奋性神经元之间的同源性与差异----------
library("Seurat")
library("ggplot2")
hippo <- readRDS(file = "./Hippo\\subcluster\\Excitatory\\anno\\Excitatory_sub.rds")
hippo_clusters<-hippo@meta.data$Anno_Idents

table(hippo$orig.ident) # 10
table(hippo@meta.data$Anno_Idents)
temp <- readRDS(file = "./temporal_lobe\\subcluster\\Excitatory\\Excitatory_anno_cluster_1.rds")
table(temp@meta.data$Anno_Idents)
table(temp[["idents"]])
temp_clusters<-temp@meta.data$Anno_Idents
temp<-NormalizeData(temp, normalization.method="LogNormalize", scale.factor=10000)
temp <- FindVariableFeatures(temp, selection.method = "vst", nfeatures = 2000)
var.genes<-intersect(VariableFeatures(hippo),VariableFeatures(temp))

###斯皮尔曼----------
hippo_exp <-AverageExpression(hippo,
                              group.by = "Anno_Idents",
                              assays = "RNA")
hippo_exp=hippo_exp[[1]]
head(hippo_exp)
hippo_exp <- hippo_exp[var.genes,]

temp_exp <-AverageExpression(temp,
                             group.by = "Anno_Idents",
                             assays = "RNA")
temp_exp=temp_exp[[1]]
head(temp_exp)
temp_exp <- temp_exp[var.genes,]

exp <- cbind(hippo_exp,temp_exp)
pheatmap::pheatmap(cor(exp,method = 'spearman')) #
View(cor(exp,method = 'spearman'))
library(RColorBrewer)
library(corrplot)
M <-cor(exp,method = 'spearman')
color <- c("#3A6963","#80B1D3","#FDB462","#E59CC4", "#BC80BD")
corrplot(M,method = 'square',
         is.corr = F,order = "AOE",addCoef.col = "black",number.cex = 0.8,
         tl.cex=0.8,tl.col = "black",col = color,
         cl.length = 5,tl.srt = 45)
####海马体兴奋性神经元与颞叶兴奋性神经元之间存在同源性，CA3  CA1与L6相似  接下来是DG5,4,3,2，1


###查看海马体抑制性神经元与颞叶皮层抑制性神经元之间的同源性与差异----------
library("Seurat")
library("ggplot2")
Inhibitory <-readRDS(file = "./subcluster\\Inhibitory\\anno\\Inhibitory_sub.rds")
DimPlot(Inhibitory, reduction='umap',group.by="Anno_Idents", pt.size=3,label=T,label.size = 5,raster = F,
        cols = cols)

hippo_clusters<-Inhibitory@meta.data$Anno_Idents

table(Inhibitory$orig.ident) # 10
table(Inhibitory@meta.data$Anno_Idents)
# PVALB   SST   VIP   CCK LAMP5 
# 1467   964  1191  1345  1308
temp <- readRDS(file = ".\\temporal_lobe\\subcluster\\Inhibitory\\Inhibitory_anno_cluster.rds")
table(temp@meta.data$Anno_Idents)
# In_Vip   In_Sst In_Pvalb In_Lamp5 
# 14632    11543     7536     5901
table(temp[["idents"]])
temp_clusters<-temp@meta.data$Anno_Idents
temp<-NormalizeData(temp, normalization.method="LogNormalize", scale.factor=10000)
temp <- FindVariableFeatures(temp, selection.method = "vst", nfeatures = 2000)
var.genes<-intersect(VariableFeatures(Inhibitory),VariableFeatures(temp))
###斯皮尔曼----------
hippo_exp <-AverageExpression(Inhibitory,
                              group.by = "Anno_Idents",
                              assays = "RNA")
hippo_exp=hippo_exp[[1]]
head(hippo_exp)
hippo_exp <- hippo_exp[var.genes,]

temp_exp <-AverageExpression(temp,
                             group.by = "Anno_Idents",
                             assays = "RNA")
temp_exp=temp_exp[[1]]
head(temp_exp)

temp_exp <- temp_exp[var.genes,]

exp <- cbind(hippo_exp,temp_exp)
pheatmap::pheatmap(cor(exp,method = 'spearman')) #默认是Pearson
View(cor(exp,method = 'spearman'))
library(RColorBrewer)
library(corrplot)
M <-cor(exp,method = 'spearman')
color <- c("#3A6963","#80B1D3","#FDB462","#E59CC4", "#BC80BD")
corrplot(M,method = 'square',
         is.corr = F,order = "AOE",addCoef.col = "black",number.cex = 0.8,
         tl.cex=0.8,tl.col = "black",col = color,
         cl.length = 5,tl.srt = 45)














