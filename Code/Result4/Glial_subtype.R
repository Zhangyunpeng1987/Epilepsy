##*星形胶质细胞---------
scRNA <- readRDS(file = "./Hippo\\scRNA_anno_cluster.rds")
table(scRNA$Anno_Idents)
Astrocyte <- subset(scRNA,Anno_Idents %in% "Astrocyte")
table(Astrocyte$group)
metadata <- Astrocyte@meta.data
dim(metadata)
metadata <- metadata[,-c(6:18,20)]
Astrocyte <- CreateSeuratObject(counts = Astrocyte@assays$RNA@counts,
                                meta.data = metadata) 
dim(Astrocyte)  
####
####
####
####**seurat流程--------
Astrocyte <- NormalizeData(Astrocyte, normalization.method = "LogNormalize", scale.factor = 10000)
Astrocyte <- FindVariableFeatures(Astrocyte, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Astrocyte)
Astrocyte <- ScaleData(Astrocyte, features = all.genes)
Astrocyte <- RunPCA(Astrocyte, features = VariableFeatures(object = Astrocyte))
DimPlot(Astrocyte, reduction = "pca")


ElbowPlot(Astrocyte,ndims = 50)
Astrocyte <- FindNeighbors(Astrocyte, dims = 1:20)
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  Astrocyte <- FindClusters(Astrocyte, resolution = res)
}

Astrocyte <- RunUMAP(Astrocyte, dims = 1:20)
Astrocyte <- RunTSNE(Astrocyte, dims = 1:20)
#####**手动注释---------
Astrocyte.markers <- FindAllMarkers(Astrocyte, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Astrocyte.markers <- read.table("./Hippo\\subcluster\\Astrocyte\\Astrocyte.markers.txt",sep = "\t")
Astrocyte.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(Astrocyte, features = top5$gene) + scale_fill_gradientn(colors = c("#006699","#FFFFCC","#FF6600"))
####参考文献 1. 参考A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation   2. Resolving cellular and molecular diversity along the hippocampal anterior-to-posterior axis in humans
library(dplyr) 
d <- dist(Astrocyte@reductions[["pca"]]@cell.embeddings, method = "euclidean")
sample_cor <- cor(Matrix::t(Astrocyte@reductions[["pca"]]@cell.embeddings))
sample_cor <- (1 - sample_cor)/2
d2 <- as.dist(sample_cor)
h_euclidean <- hclust(d, method = "ward.D2")
h_correlation <- hclust(d2, method = "ward.D2")
Astrocyte$hc_euclidean_5 <- cutree(h_euclidean,k = 5)

Astrocyte$hc_corelation_5 <- cutree(h_correlation,k = 5)
###采用hc_corelation_5--------------
Idents(Astrocyte) <- 'hc_corelation_5'  #5 cluster
saveRDS(Astrocyte,"./Hippo\\subcluster\\Astrocyte\\anno\\Astrocyte.rds")

Astrocyte.markers <- FindAllMarkers(Astrocyte, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Astrocyte.markers,"./Hippo\\subcluster\\Astrocyte\\anno\\Astrocyte.markers_hc_corelation_5.txt",sep = "\t",quote = F)
Astrocyte.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
Astrocyte.markers[which(Astrocyte.markers$gene %in% "EGFR"),]
DotPlot(Astrocyte, features = c("TNFSF4","IL33","NFKB1","JAK1","STAT3","IKBKB","CABLES1","GFAP","SOX2","EGFR","NR4A2","PDE1C","FGF12","MAP1B","MEF2C"),
        group.by = "group",cols = c("lightgrey", "#CC3333"))+coord_flip()+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size = 11,colour = "black",angle = 30),
        axis.text.y = element_text(size = 11, colour = "black"))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
sce.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DotPlot(Astrocyte, features  =top5$gene,
        group.by = "Anno_Idents",cols = c("lightgrey", "#CC3333"))+coord_flip()+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size = 11,colour = "black",angle = 30),
        axis.text.y = element_text(size = 11, colour = "black"))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))

library(dplyr)
sce.markers <- read.table("./Hippo\\subcluster\\Astrocyte\\anno\\Astrocyte.markers_hc_corelation_5.txt",sep = "\t")
###### step2:富集分析 ###### 
library(ClusterGVis)
library(org.Hs.eg.db)
scRNA.markers <- sce.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)

head(scRNA.markers)

scRNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> markGenes

L6.enrich.go <- enrichGO(gene = scRNA.markers$gene,  #基因列表文件中的基因名称
                         OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                         keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                         ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                         pAdjustMethod = 'fdr',  #指定 p 值校正方法
                         pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                         qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                         readable = FALSE)
L6.enrich.go <-summary(L6.enrich.go)  #79
write.table(L6.enrich.go,"./Hippo\\subcluster\\Astrocyte\\Astrocyte_sub_top20.enrich.go.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

GO <- L6.enrich.go[L6.enrich.go$Description %in% enrich$Description[enrich$group %in% unique(enrich$group)[1]],3]
for (i in 1:length(GO)) {
  genes <- L6.enrich.go[L6.enrich.go$Description %in% GO[i],9]
  genes <- str_split(genes,"/")[[1]]
  genes
  aaa <- intersect(genes,scRNA.markers$gene[scRNA.markers$cluster %in% unique(scRNA.markers$cluster)[9]])
  enrich$num[enrich$Description %in% GO[i]] <- length(aaa)
  genesymbol <- paste(aaa,collapse=",")
  enrich$gene[enrich$Description %in% GO[i]] <- genesymbol
}

st.data <- prepareDataFromscRNA(object = Astrocyte,
                                diffData = scRNA.markers,
                                showAverage = F,
                                keep.uniqGene = FALSE,
                                sep = "_")
str(st.data)
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 5201314)
write.table(enrich,"./Hippo\\subcluster\\Astrocyte\\Astrocyte_sub_top20.visCluster.go.txt",sep = '\t',col.names = T,row.names = T,quote = FALSE,na='')
Astrocyte_sub <- read.table("./Hippo\\subcluster\\Astrocyte\\Astrocyte_sub_top20.visCluster.go.txt",header = T,sep = "\t")
head(enrich)
visCluster(object = st.data,
           plot.type = "line")
pdf('./Hippo\\subcluster\\Astrocyte\\Astrocyte细胞类型间差异基因功能富集.pdf',height = 10,width = 20,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 45,
           markGenes = unique(markGenes$gene),
           markGenes.side = "left",
           annoTerm.data = enrich,
           genes.gp = c('italic',fontsize = 12,col = "black"),
           show_column_names = F,
           line.side = "left",
           cluster.order = c(1:5),
           add.bar = T,
           #sample.cell.order = rev(Anno_Idents),
           sample.col = c("#6699CC","#339999","#CCCC99","#FFCC00","#FF9900"))
dev.off()
####**细胞注释完成-------------
# 1: GFAP   # 2:GFAP         #3:GFAP    # 4:EGFR,CABLES1,NR4A2  # 5: EGFR,CABLES1
levels(Astrocyte)
new.cluster.ids <- c("AST1","AST2","AST3","AST4","AST5")
unique(new.cluster.ids)
names(new.cluster.ids) <- levels(Astrocyte)
Astrocyte <- RenameIdents(Astrocyte, new.cluster.ids)
Astrocyte <- StashIdent(Astrocyte, save.name = 'Anno_Idents')
cols<- c("#3A6963","#BC80BD","#E59CC4","#cabbe9","#58A4C3")

DimPlot(Astrocyte, reduction='umap',group.by="Anno_Idents", pt.size=1.5,label=F,label.size = 3,raster = F,
        label.box = F,cols = cols)
ggsave("./Hippo\\subcluster\\Astrocyte\\anno\\Anno_Idents.pdf",width = 13,height = 11)  
jjDotPlot(object = Astrocyte,
          gene = c("GFAP","HSPB1","CD44","EGFR","CABLES1","WIF1","ETNPPL","NR4A2","NFKB1","JAK1","STAT3"),
          id = 'Anno_Idents',
          xtree = F,
          ytree = F,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1,
          point.shape = 22)
ggsave("./Hippo\\subcluster\\Astrocyte\\anno\\anno_dotplot.pdf",width = 8,height = 8)  

###细胞比例----------
names(table(Astrocyte$Anno_Idents))
all.an <- prop.table(table(Astrocyte$Anno_Idents[Astrocyte@meta.data$group %in% c("A")]))
All.po <- prop.table(table(Astrocyte$Anno_Idents[Astrocyte@meta.data$group %in% c("P")]))
df_prop <- cbind(all.an,All.po)
df_prop <- t(apply(df_prop, 1, function(x) log(x / x[1])))
df_prop <- as.data.frame(df_prop)
df_prop$celltype <- rownames(df_prop)
df_prop <- df_prop[,-1]
class(df_prop)
colnames(df_prop) <- c("log","celltype")
ggplot(df_prop, aes(x=celltype,y=log))+ 
  geom_bar(position="dodge",stat="identity",width=0.8)+
  theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size =8),
        panel.background = element_rect(color = 'black', fill = 'transparent'))
ggsave("./Hippo\\subcluster\\Astrocyte\\anno\\logPA.pdf",height = 5,width = 5)

######细胞类型在前后端比例----
AST1 <- prop.table(table(Astrocyte$group[Astrocyte@meta.data$Anno_Idents %in% c("AST1")]))
AST2 <- prop.table(table(Astrocyte$group[Astrocyte@meta.data$Anno_Idents %in% c("AST2")]))
AST3 <- prop.table(table(Astrocyte$group[Astrocyte@meta.data$Anno_Idents %in% c("AST3")]))
AST4 <- prop.table(table(Astrocyte$group[Astrocyte@meta.data$Anno_Idents %in% c("AST4")]))
AST5 <- prop.table(table(Astrocyte$group[Astrocyte@meta.data$Anno_Idents %in% c("AST5")]))

cluster <- c(sort(rep(names(table(Astrocyte$Anno_Idents)),2)))
pos <- c(rep(rep(names(table(Astrocyte$group))),5))
library("ggplot2")
cell.prop<-as.data.frame(c(AST1,AST2,AST3,AST4,AST5),pos)
cell.prop$pos <- rownames(cell.prop) 
cell.prop$cluster <- cluster
colnames(cell.prop)<-c("proportion","pos","cluster")
library('ggplot2')
library('reshape2')
ggplot(cell.prop,aes(cluster,proportion,fill=pos))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+scale_fill_manual(values = c("#99CCCC","#FFCC99"))+
  theme(axis.text.x = element_text(angle=30))
ggsave("./Hippo\\subcluster\\Astrocyte\\anno\\细胞类型比例.pdf",height = 5,width = 5)

##*小胶质细胞---------
scRNA <- readRDS(file = "./Hippo\\scRNA_anno_cluster.rds")
table(scRNA$Anno_Idents)
Microglial <- subset(scRNA,Anno_Idents %in% "Microglial")
table(Microglial$group)
metadata <- Microglial@meta.data
dim(metadata)
metadata <- metadata[,-c(6:18,20)]
Microglial <- CreateSeuratObject(counts = Microglial@assays$RNA@counts,
                                 meta.data = metadata) 
dim(Microglial)  
####**seurat流程--------
Microglial <- NormalizeData(Microglial, normalization.method = "LogNormalize", scale.factor = 10000)
Microglial <- FindVariableFeatures(Microglial, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Microglial)
Microglial <- ScaleData(Microglial, features = all.genes)
Microglial <- RunPCA(Microglial, features = VariableFeatures(object = Microglial))
DimPlot(Microglial, reduction = "pca")
ElbowPlot(Microglial,ndims = 50)
Microglial <- FindNeighbors(Microglial, dims = 1:10)
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  Microglial <- FindClusters(Microglial, resolution = res)
}
Microglial <- RunUMAP(Microglial, dims = 1:10)
Microglial <- RunTSNE(Microglial, dims = 1:10)
library(clustree)
library(patchwork)
clustree(Microglial, prefix = 'RNA_snn_res.') + coord_flip()
ggsave("./Hippo\\subcluster\\Microglial\\clustree.pdf",width = 13,height = 11)
DimPlot(Microglial, group.by = 'RNA_snn_res.0.3', pt.size = 1,label = T,reduction = "umap",label.box=T,cols = mycolors)
Idents(Microglial) <- 'RNA_snn_res.0.3'  
ggsave("./Hippo\\subcluster\\Microglial\\RNA_snn_res.0.3.pdf",width = 13,height = 11)
#####**手动注释---------
Microglial.markers <- FindAllMarkers(Microglial, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Microglial.markers,"./Hippo\\subcluster\\Microglial\\Microglial.markers.txt",sep = "\t",quote = F)
Microglial.markers <- read.table("./Hippo\\subcluster\\Microglial\\Microglial.markers.txt",sep = "\t")
Microglial.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(Microglial, features = top5$gene) + scale_fill_gradientn(colors = c("#006699","#FFFFCC","#FF6600"))

library(ClusterGVis)
library(org.Hs.eg.db)
table(scRNA$group)
table(scRNA$Anno_Idents)
Idents(scRNA) <- scRNA$group
scRNA.markers <- Microglial.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)

head(scRNA.markers)

scRNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> markGenes

st.data <- prepareDataFromscRNA(object = Microglial,
                                diffData = scRNA.markers,
                                showAverage = F,
                                keep.uniqGene = FALSE,
                                sep = "_")
str(st.data)
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 5201314)
write.table(enrich,"./Hippo\\subcluster\\Microglial\\Microglial_sub_top20.visCluster.go.txt",sep = '\t',col.names = T,row.names = T,quote = FALSE,na='')

head(enrich)
visCluster(object = st.data,
           plot.type = "line")
# heatmap plot
pdf('./Hippo\\subcluster\\Microglial\\Anno\\小胶质细胞类型间差异基因功能富集.pdf',height = 10,width = 20,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 45,
           markGenes = unique(markGenes$gene),
           markGenes.side = "left",
           annoTerm.data = enrich,
           genes.gp = c('italic',fontsize = 12,col = "black"),
           show_column_names = F,
           line.side = "left",
           cluster.order = c(1:9),
           add.bar = T,
           #sample.cell.order = rev(Anno_Idents),
           sample.col = c("#CCCCFF","#FFCC00","#0099CC","#99CC99","#CC6699","#FF9966","#CC9999","blue","red"))
dev.off()
###**综合功能富集和基因表达，手动注释完成------------
#0,6  KANK1,DOCK5  1
#1  CCL4,CCL3L3,HIF1A,SPP1  2
#4 ST6GAL1,DOCK8     2 
#7 "CD74","CD86"     2
#2,3  TGFBR2,CSF1R  P2RY12  3  
#5  "DPP10","SERPINE1"    4
#8  "FYN","CBLB","SKAP1"  5

DotPlot(Microglial, features = unique(c("KANK1","DOCK5","SPP1","CCL4","HIF1A","ST6GAL1","DOCK8",
                                        "CD74","CD86","CD163",
                                        "P2RY12","TGFBR2","CSF1R","DPP10","SERPINE1","FYN","CBLB","SKAP1")),
        group.by = "RNA_snn_res.0.3",cols = c("lightgrey", "#CC3333"))+coord_flip()+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size = 11,colour = "black",angle = 30),
        axis.text.y = element_text(size = 11, colour = "black"))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
levels(Microglial)
new.cluster.ids <- c("Micro1","Micro2","Micro3","Micro3","Micro2","Micro4",
                     "Micro1","Micro2","Micro5")
unique(new.cluster.ids)
names(new.cluster.ids) <- levels(Microglial)
Microglial <- RenameIdents(Microglial, new.cluster.ids)
Microglial <- StashIdent(Microglial, save.name = 'Anno_Idents')
cols<- c("#8F797E","#FFC2B5","#FFE3CC","#646C8F","#DCC3A1")
DimPlot(Microglial, reduction='umap',group.by="Anno_Idents", pt.size=1.5,label=F,label.size = 3,raster = F,
        label.box = F,cols = cols)
saveRDS(Microglial,file = "./Hippo\\subcluster\\Microglial\\anno\\Microglial_sub.rds")
DotPlot(Microglial, features = unique(c("KANK1","DOCK5","SPP1","CCL4","HIF1A","ST6GAL1","DOCK8","CD74","CD86",
                                        "P2RY12","TGFBR2","CSF1R","DPP10","SERPINE1","FYN","CBLB","SKAP1")),
        group.by = "Anno_Idents",cols = c("lightgrey", "#CC3333"))+coord_flip()+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size = 11,colour = "black",angle = 30),
        axis.text.y = element_text(size = 11, colour = "black"))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
ggsave('./Hippo\\subcluster\\Microglial\\anno\\dotplot_anno_idents.pdf',height = 11,width = 13)
#0,6  KANK1,DOCK5  1
#1  CCL4,CCL3L3,HIF1A,SPP1  2
#4 ST6GAL1,DOCK8     2 
#7 "CD74","CD86"     2
#2,3  TGFBR2,CSF1R  P2RY12  3  
#5  "DPP10","SERPINE1"    4
#8  "FYN","CBLB","SKAP1"  5

jjDotPlot(object = Microglial,
          gene = c("KANK1","DOCK5","SPP1","CCL4","ST6GAL1","DOCK8","CD74","CD86",
                   "P2RY12","CSF1R","TMEM119","SALL1","DPP10","SERPINE1","FYN","CBLB"),
          id = 'Anno_Idents',
          xtree = F,
          ytree = F,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1,
          point.shape = 22)
ggsave("./Hippo\\subcluster\\Microglial\\anno\\anno_dotplot.pdf",width = 8,height = 8)  

###细胞比例----------
Microglial <- readRDS(file = "./Hippo\\subcluster\\Microglial\\anno\\Microglial_sub.rds")
cols<- c("#8F797E","#FFC2B5","#FFE3CC","#646C8F","#DCC3A1")
DimPlot(Microglial, reduction='umap',group.by="Anno_Idents", pt.size=1.5,label=T,label.size = 3,raster = F,
        label.box = T,cols = cols)

names(table(Microglial$Anno_Idents))
all.an <- prop.table(table(Microglial$Anno_Idents[Microglial@meta.data$group %in% c("A")]))
All.po <- prop.table(table(Microglial$Anno_Idents[Microglial@meta.data$group %in% c("P")]))
df_prop <- cbind(all.an,All.po)
df_prop <- t(apply(df_prop, 1, function(x) log(x / x[1])))
df_prop <- as.data.frame(df_prop)
df_prop$celltype <- rownames(df_prop)
df_prop <- df_prop[,-1]
class(df_prop)
colnames(df_prop) <- c("log","celltype")
p1 <- ggplot(df_prop, aes(x=celltype,y=log))+ 
  geom_bar(position="dodge",stat="identity",width=0.8)+
  theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size =8),
        panel.background = element_rect(color = 'black', fill = 'transparent'))
ggsave("./Hippo\\subcluster\\Microglial\\anno\\logPA.pdf",height = 5,width = 5)


######细胞类型在前后端比例----
Micro1 <- prop.table(table(Microglial$group[Microglial@meta.data$Anno_Idents %in% c("Micro1")]))
Micro2 <- prop.table(table(Microglial$group[Microglial@meta.data$Anno_Idents %in% c("Micro2")]))
Micro3 <- prop.table(table(Microglial$group[Microglial@meta.data$Anno_Idents %in% c("Micro3")]))
Micro4 <- prop.table(table(Microglial$group[Microglial@meta.data$Anno_Idents %in% c("Micro4")]))
Micro5 <- prop.table(table(Microglial$group[Microglial@meta.data$Anno_Idents %in% c("Micro5")]))

cluster <- c(sort(rep(names(table(Microglial$Anno_Idents)),2)))
pos <- c(rep(rep(names(table(Microglial$group))),5))
library("ggplot2")
cell.prop<-as.data.frame(c(Micro1,Micro2,Micro3,Micro4,Micro5),pos)
cell.prop$pos <- rownames(cell.prop) 
cell.prop$cluster <- cluster
colnames(cell.prop)<-c("proportion","pos","cluster")
library('ggplot2')
library('reshape2')

p2<- ggplot(cell.prop,aes(cluster,proportion,fill=pos))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+scale_fill_manual(values = c("#99CCCC","#FFCC99"))+
  theme(axis.text.x = element_text(angle=30))
ggsave("./Hippo\\subcluster\\Microglial\\anno\\细胞类型比例.pdf",p2,height = 5,width = 5)

##*少突胶质细胞-------------
scRNA <- readRDS(file = "./Hippo\\scRNA_anno_cluster.rds")
table(scRNA$Anno_Idents)
Oligodendrocyte <- subset(scRNA,Anno_Idents %in% c("Oligodendrocyte","Oligodendrocyte progenitor cell"))
table(scRNA$Anno_Idents)
table(Oligodendrocyte$group)
dim(Oligodendrocyte)  #16956 77544
metadata <- Oligodendrocyte@meta.data
dim(metadata)
metadata <- metadata[,-c(6:18,20)]
Oligodendrocyte <- CreateSeuratObject(counts = Oligodendrocyte@assays$RNA@counts,
                                      meta.data = metadata) 
dim(Oligodendrocyte) 

####**seurat流程--------
Oligodendrocyte <- NormalizeData(Oligodendrocyte, normalization.method = "LogNormalize", scale.factor = 10000)
Oligodendrocyte <- FindVariableFeatures(Oligodendrocyte, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Oligodendrocyte)
Oligodendrocyte <- ScaleData(Oligodendrocyte, features = all.genes)
Oligodendrocyte <- RunPCA(Oligodendrocyte, features = VariableFeatures(object = Oligodendrocyte))
DimPlot(Oligodendrocyte, reduction = "pca")

ElbowPlot(Oligodendrocyte,ndims = 50)
Oligodendrocyte <- FindNeighbors(Oligodendrocyte, dims = 1:10)
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  Oligodendrocyte <- FindClusters(Oligodendrocyte, resolution = res)
}
Oligodendrocyte <- RunUMAP(Oligodendrocyte, dims = 1:10)
Oligodendrocyte <- RunTSNE(Oligodendrocyte, dims = 1:10)
library(clustree)
library(patchwork)
clustree(Oligodendrocyte, prefix = 'RNA_snn_res.') + coord_flip()
ggsave("./Hippo\\subcluster\\Oligodendrocyte\\clustree.pdf",width = 13,height = 11)
DimPlot(Oligodendrocyte, group.by = 'RNA_snn_res.0.4', pt.size = 1,label = T,reduction = "umap",label.box=T,cols = mycolors)
Idents(Oligodendrocyte) <- 'RNA_snn_res.0.4'  #13 cluster
levels(Oligodendrocyte)  
ggsave("./Hippo\\subcluster\\Oligodendrocyte\\RNA_snn_res.0.4.pdf",width = 13,height = 11)
levels(Oligodendrocyte)
DimPlot(Oligodendrocyte, group.by = 'RNA_snn_res.0.4', pt.size = 1,label = T,reduction = "umap",label.box=T,cols = mycolors)

#####**手动注释---------
Oligodendrocyte.markers <- FindAllMarkers(Oligodendrocyte, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Oligodendrocyte.markers,"./Hippo\\subcluster\\Oligodendrocyte\\Oligodendrocyte.markers.txt",sep = "\t",quote = F)
Oligodendrocyte.markers <- read.table("./Hippo\\subcluster\\Oligodendrocyte\\Oligodendrocyte.markers.txt",sep = "\t")
Oligodendrocyte.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(Oligodendrocyte, features = top5$gene) + scale_fill_gradientn(colors = c("#006699","#FFFFCC","#FF6600"))
####参考文献 1. 参考Single-cell transcriptomic reveals molecular diversity and developmental heterogeneity of human stem cell-derived oligodendrocyte lineage cells
##A single-cell transcriptome atlas of glial diversity in the human hippocampus across the postnatal lifespan
c("PCDH15","SOX6","FGF12")  #祖细胞基因 4,7,11,12
c("OPALIN","CNP") #中间少突细胞  1.9 神经胶质细胞发育（CNP、CD9）
c("SERINC3") #凋亡信号   3  终末状态
c("MAG","MOG","ZNF488","KLK6","MOBP")#成熟少突胶质
c("STMN2","NREP","MAP1B,SOX11") #神经发生和神经元分化  
c("SOX6","SIRT2","SOX10","NREP")#神经胶质生成（SOX6、SIRT2、SOX10）12
c("HEY1","TRO","HIS1")  #神经发生（HEY1、TRO、HIS1）#未表达
c("IFI6","ISG15","IFIT1","HLA-A","HLA-B","HLA-C")#细胞因子 11,12
# 神经胶质细胞发育（PLP1、CNP、CD9）和凋亡信号（SEPTIN4、SERINC3）的转录本。  PMID35381189
c("MOBP","KLK6","OPALIN","BCAS1","SOX6","PCDH15","CSPG4","PDGFRA","CD74") #祖细胞、未成熟、成熟
FeaturePlot(Oligodendrocyte,features = c("CD74","BCAS1"),cols = c("grey","#B40F20"))
DotPlot(Oligodendrocyte, features = c("MOBP","KLK6","OPALIN","BCAS1","OLIG2","SOX6","PCDH15","CSPG4","PDGFRA","CD74"),
        group.by = "RNA_snn_res.0.4",cols = c("lightgrey", "#CC3333"))+coord_flip()+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size = 11,colour = "black",angle = 30),
        axis.text.y = element_text(size = 11, colour = "black"))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))

library(ClusterGVis)
library(org.Hs.eg.db)
scRNA.markers <- Oligodendrocyte.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)

head(scRNA.markers)

scRNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> markGenes

L6.enrich.go <- enrichGO(gene = scRNA.markers$gene,  #基因列表文件中的基因名称
                         OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                         keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                         ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                         pAdjustMethod = 'fdr',  #指定 p 值校正方法
                         pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                         qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                         readable = FALSE)
L6.enrich.go <-summary(L6.enrich.go)  #79
write.table(L6.enrich.go,"./Hippo\\subcluster\\Oligodendrocyte\\Oligodendrocyte_sub_top20.enrich.go.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

GO <- L6.enrich.go[L6.enrich.go$Description %in% enrich$Description[enrich$group %in% c("C1")],3]
length(GO)


GO <- L6.enrich.go[L6.enrich.go$Description %in% enrich$Description[enrich$group %in% unique(enrich$group)[13]],3]
for (i in 1:length(GO)) {
  genes <- L6.enrich.go[L6.enrich.go$Description %in% GO[i],9]
  genes <- str_split(genes,"/")[[1]]
  genes
  aaa <- intersect(genes,scRNA.markers$gene[scRNA.markers$cluster %in% unique(scRNA.markers$cluster)[9]])
  enrich$num[enrich$Description %in% GO[i]] <- length(aaa)
  genesymbol <- paste(aaa,collapse=",")
  enrich$gene[enrich$Description %in% GO[i]] <- genesymbol
}

st.data <- prepareDataFromscRNA(object = Oligodendrocyte,
                                diffData = scRNA.markers,
                                showAverage = F,
                                keep.uniqGene = FALSE,
                                sep = "_")
str(st.data)
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 5201314)
write.table(enrich,"./Hippo\\subcluster\\Oligodendrocyte\\Oligodendrocyte_sub_top20.visCluster.go.txt",sep = '\t',col.names = T,row.names = T,quote = FALSE,na='')
enrich<-read.table("./Hippo\\subcluster\\Oligodendrocyte\\Oligodendrocyte_sub_top20.visCluster.go.txt",sep = '\t',header = T)

head(enrich)
visCluster(object = st.data,
           plot.type = "line")
pdf('./Hippo\\subcluster\\Oligodendrocyte\\anno\\少突胶质细胞类型间差异基因功能富集.pdf',height = 10,width = 20,onefile = F)
aa <- visCluster(object = st.data,
                 plot.type = "both",
                 column_title_rot = 45,
                 markGenes = unique(markGenes$gene),
                 markGenes.side = "left",
                 annoTerm.data = enrich,
                 genes.gp = c('italic',fontsize = 12,col = "black"),
                 show_column_names = F,
                 line.side = "left",
                 cluster.order = c(1:13),
                 add.bar = T,
                 #sample.cell.order = rev(Anno_Idents),
                 sample.col = mycolors)
dev.off()
###**综合功能富集和基因表达，手动注释完成------------
Idents(Oligodendrocyte) <- "RNA_snn_res.0.4"
new.cluster.ids <- c("C0","C1","C2","C3","C4","C5",
                     "C6","C7","C8","C9","C10",
                     "C11","C12")
# unique(new.cluster.ids)
names(new.cluster.ids) <- levels(Oligodendrocyte)
Oligodendrocyte <- RenameIdents(Oligodendrocyte, new.cluster.ids)
Oligodendrocyte <- StashIdent(Oligodendrocyte, save.name = 'celltype')
cols<- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C')
DimPlot(Oligodendrocyte, reduction='umap',group.by="celltype", pt.size=1.5,label=T,label.size = 3,raster = F,
        label.box = T,cols = cols)

jjDotPlot(object = Oligodendrocyte,
          gene = c("PCDH15","SOX6","FGF12",#祖细胞 4.7.11.12
                   "BCAS1",#神经发生、神经分化12   未成熟细胞
                   "IFI6","HLA-A","HLA-B","HLA-C",#细胞因子 11,12
                   "MOBP","MAG","MOG",#成熟胶质035
                   "OPALIN",#中间成熟胶质  1，9
                   "ZNF565",#2 &9 中间
                   "LINGO1","CIRBP","DNAJB2",#mRNA加工  6.10
                   "CNTN5","LINGO2"#突触相关  8
          ),
          id = 'Anno_Idents',
          xtree = F,
          ytree = T,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1,
          point.shape = 22)
ggsave("./Hippo\\subcluster\\Oligodendrocyte\\anno\\anno_dotplot.pdf",width = 8,height = 8)  

saveRDS(Oligodendrocyte,file = "./Hippo\\subcluster\\Oligodendrocyte\\anno\\Oligodendrocyte_sub.rds")

###细胞比例----------
all.an <- prop.table(table(Oligodendrocyte$Anno_Idents[Oligodendrocyte@meta.data$group %in% c("A")]))
All.po <- prop.table(table(Oligodendrocyte$Anno_Idents[Oligodendrocyte@meta.data$group %in% c("P")]))
df_prop <- cbind(all.an,All.po)
df_prop <- t(apply(df_prop, 1, function(x) log(x / x[1])))
df_prop <- as.data.frame(df_prop)
df_prop$celltype <- rownames(df_prop)
df_prop <- df_prop[,-1]
class(df_prop)
colnames(df_prop) <- c("log","celltype")
p1 <- ggplot(df_prop, aes(x=celltype,y=log))+ 
  geom_bar(position="dodge",stat="identity",width=0.8)+
  theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size =8),
        panel.background = element_rect(color = 'black', fill = 'transparent'))
ggsave("./Hippo\\subcluster\\Oligodendrocyte\\anno\\logPA.pdf",p1,height = 5,width = 5)


######细胞类型在前后端比例----
mOli <- prop.table(table(Oligodendrocyte$group[Oligodendrocyte@meta.data$Anno_Idents %in% c("mOli")]))
imOli <- prop.table(table(Oligodendrocyte$group[Oligodendrocyte@meta.data$Anno_Idents %in% c("imOli")]))
OPCs <- prop.table(table(Oligodendrocyte$group[Oligodendrocyte@meta.data$Anno_Idents %in% c("OPCs")]))
Oli1 <- prop.table(table(Oligodendrocyte$group[Oligodendrocyte@meta.data$Anno_Idents %in% c("Oli1")]))
Oli2 <- prop.table(table(Oligodendrocyte$group[Oligodendrocyte@meta.data$Anno_Idents %in% c("Oli2")]))

cluster <- c(sort(rep(names(table(Oligodendrocyte$Anno_Idents)),2)))
pos <- c(rep(rep(names(table(Oligodendrocyte$group))),5))
library("ggplot2")
cell.prop<-as.data.frame(c(mOli,imOli,OPCs,Oli1,Oli2),pos)
cell.prop$pos <- rownames(cell.prop) 
cell.prop$cluster <- cluster
colnames(cell.prop)<-c("proportion","pos","cluster")
library('ggplot2')
library('reshape2')

p2<- ggplot(cell.prop,aes(cluster,proportion,fill=pos))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+scale_fill_manual(values = c("#99CCCC","#FFCC99"))+
  theme(axis.text.x = element_text(angle=30))
ggsave("./Hippo\\subcluster\\Oligodendrocyte\\anno\\细胞类型比例.pdf",p2,height = 5,width = 5)
