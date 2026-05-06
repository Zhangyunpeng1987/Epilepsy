###亚型分析-------------
###*兴奋性神经元----------
scRNA <- readRDS(file = "./scRNA_anno_cluster.rds")
metadata <- scRNA@meta.data
##meta.data添加信息
scRNA@meta.data$color[metadata$Anno_Idents %in% "Oligodendrocyte"] <- "#CCCCFF"
scRNA@meta.data$color[metadata$Anno_Idents %in% "Oligodendrocyte progenitor cell"] <- "#FFCC00"
scRNA@meta.data$color[metadata$Anno_Idents %in% "Microglial"] <- "#0099CC"
scRNA@meta.data$color[metadata$Anno_Idents %in% "Astrocyte"] <- "#99CC99"
scRNA@meta.data$color[metadata$Anno_Idents %in% "Excitatory neuron"] <- "#CC6699"
scRNA@meta.data$color[metadata$Anno_Idents %in% "Inhibitory neuron"] <- "#FF9966"
scRNA@meta.data$color[metadata$Anno_Idents %in% "Endothelial"] <- "#CC9999"
saveRDS(scRNA,"./scRNA_anno_cluster.rds")
DimPlot(scRNA, reduction='umap', group.by="Anno_Idents", pt.size=0.3,raster=FALSE,
        label=T,label.size = 5,cols = c("#CCCCFF","#FFCC00","#0099CC","#99CC99","#CC6699","#FF9966","#CC9999"))
Excitatory <- subset(scRNA,Anno_Idents %in% "Excitatory neuron")

metadata <- Excitatory@meta.data
dim(metadata)
metadata <- metadata[,-c(6:18,20)]
Excitatory <- CreateSeuratObject(counts = Excitatory@assays$RNA@counts,
                                 meta.data = metadata) 
dim(Excitatory)  # 16956 16638
####
####
####
####**seurat流程--------
Excitatory <- NormalizeData(Excitatory, normalization.method = "LogNormalize", scale.factor = 10000)
Excitatory <- FindVariableFeatures(Excitatory, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Excitatory)
Excitatory <- ScaleData(Excitatory, features = all.genes)
Excitatory <- RunPCA(Excitatory, features = VariableFeatures(object = Excitatory))
DimPlot(Excitatory, reduction = "pca")

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Excitatory <- JackStraw(Excitatory, num.replicate = 100)
Excitatory <- ScoreJackStraw(Excitatory, dims = 1:20)
#JackStrawPlot(Excitatory, dims = 1:15)
ElbowPlot(Excitatory,ndims = 50)
Excitatory <- FindNeighbors(Excitatory, dims = 1:30)
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  Excitatory <- FindClusters(Excitatory, resolution = res)
}
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
Excitatory <- RunUMAP(Excitatory, dims = 1:30)
Excitatory <- RunTSNE(Excitatory, dims = 1:30)
library(clustree)
library(patchwork)
clustree(Excitatory, prefix = 'RNA_snn_res.') + coord_flip()
DimPlot(Excitatory, group.by = 'RNA_snn_res.0.6', pt.size = 1,label = T,reduction = "umap",label.box=T,cols = mycolors)
Idents(Excitatory) <- 'RNA_snn_res.0.6'  #25 cluster


# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
saveRDS(Excitatory,"./subcluster\\Excitatory\\Excitatory-umap_tsne.rds")

Excitatory <- readRDS("./subcluster\\Excitatory\\Excitatory-umap_tsne.rds")
BrainMarker <- read.table("./CellMarker_brain_embryo_marker.txt",sep = "\t",header = T)
names(table(BrainMarker$cell_name))

####
#**小提琴图展示findallmarker基因与文献中marker的交集基因
celltype.list <- c("Broad excitatory neuron","Deep layer excitatory neuron",
                   "Excitatory neuron", 
                   "Glutamatergic neuron","Layer 4 cell","Deep layer cell","Lower layer neuron","Neuron","Neuronal precursor cell",
                   "Neuronal progenitor cell",
                   "Superficial layer excitatory neuron","Upper layer cell")
setwd("./subcluster\\Excitatory")
# pdf("E:\\original\\epilepsy_Data\\Epilepsy19-master\\FIRST/celltype_marker_scRna.pdf")
a <- c()
for(i in 1:length(celltype.list)){
  #type <- paste(celltype.list[[i]], collapse = "_")
  marker <- intersect(top5$gene,BrainMarker[BrainMarker$cell_name %in% celltype.list[[i]],2])
  if(length(marker)==0){
    name <- strsplit(celltype.list[[i]]," ")
    type <- paste(name[[1]], collapse = "_")
    print(paste(type, "Have no marker",sep = ":"))
    next;
  }
  name <- strsplit(celltype.list[[i]]," ")
  type <- paste(name[[1]], collapse = "_")
  p <- print(paste(type, marker,sep = ":"))
  a <- c(a,p)
  write.table(a,"cell_type_marker.txt",sep = "\t",quote = F,row.names = F)
  plots <- StackedVlnPlot(Excitatory, marker, pt.size=0, cols=my36colors)
  ggsave(plots,file=paste(type,".pdf",sep="_"),height = length(marker),width = 20,limitsize = FALSE)
  #print(plots)
}
DimPlot(Excitatory, reduction = "umap",group.by = "RNA_snn_res.0.6",label=T,cols = mycolors,label.box = T)
####参考文献 1. Decoding the development of the human hippocampus   2. Resolving cellular and molecular diversity along the hippocampal anterior-to-posterior axis in humans
#首先区分DG   CA
#DG: "MAML2","SEMA5A"
#CA："SV2B"
##CA1: "PID1"
##CA1/3 : "TYRO3"(不采用)
##CA2/3 : "PFKP"(不采用)
##CA3:"HS3ST4","TYRO3"
###DG1:0.3.7.14.10.21:  "SLC47A1","COLEC12",#"SLC4A4","COL6A3"
###DG2: 13 "SERPINE1"
###DG3: 15  "SLC14A1","TNC"
###DG4: 20 "LHFPL3","TMEM132C"
###DG5: 1.4.6 "MAML2","SEMA5A"
StackedVlnPlot(Excitatory, c("MAML2","SEMA5A","SV2B"),pt.size=0, cols=my36colors)
StackedVlnPlot(Excitatory, c("PID1","TYRO3","HS3ST4"),pt.size=0, cols=my36colors)
StackedVlnPlot(Excitatory, c("SLC47A1","COLEC12","SERPINE1","SLC14A1","TNC","LHFPL3","TMEM132C","MAML2","SEMA5A"),pt.size=0, cols=my36colors)
#FeaturePlot(Excitatory,features = c("PID1","STAB2"),cols = c("grey","#B40F20"))
#
# p <- DotPlot(Excitatory,features=c("SLC47A1","COLEC12","SERPINE1","SLC14A1","TNC","LHFPL3","TMEM132C","MAML2","SEMA5A"))
# p+ theme(axis.text.x = element_text(angle = 45))


####**细胞注释完成-------------
# 0:DG_Ex1   # 1:DG_Ex2    # 2:CA1_Ex         #3:DG_Ex1    # 4:DG_Ex2  # 5:CA1_Ex 
# 6:DG_Ex2   # 7:DG_Ex1    # 8:CA_Ex_GAPDH    #9:CA3_Ex    #10:DG_Ex1
# 11:CA1_Ex  # 12:CA3_Ex   # 13:DG_Ex3   # 14:DG_Ex1  # 15:DG_Ex5
# 16:CA3_Ex  # 17:CA1_Ex   # 18:CA1_Ex   # 19:CA3_Ex  # 20:DG_Ex4
# 21:DG_Ex1  # 22:CA3_Ex   # 23:CA1_Ex   # 24:DG_Ex1
#*
levels(Excitatory)
#齿状回dentate gyrus (DG)  
#CA
new.cluster.ids <- c("DG_Ex1","DG_Ex2","CA1_Ex","DG_Ex1","DG_Ex2","CA1_Ex",
                     "DG_Ex2","DG_Ex1","CA_Ex_GAPDH","CA3_Ex","DG_Ex1",
                     "CA1_Ex","CA3_Ex","DG_Ex3","DG_Ex1","DG_Ex5",
                     "CA3_Ex","CA1_Ex", "CA1_Ex","CA3_Ex","DG_Ex4",
                     "DG_Ex1","CA3_Ex", "CA1_Ex","DG_Ex1")
unique(new.cluster.ids)
names(new.cluster.ids) <- levels(Excitatory)
Excitatory <- RenameIdents(Excitatory, new.cluster.ids)
Excitatory <- StashIdent(Excitatory, save.name = 'Anno_Idents')
# cols <- c(rgb(173,206,215,150,maxColorValue = 255),#ex1
#           rgb(161,169,208,150,maxColorValue = 255), #ex3
#           rgb(246,202,229,150,maxColorValue = 255),#ex2
#           rgb(207,234,241,150,maxColorValue = 255),#ex4
#           rgb(196,165,222,150,maxColorValue = 255),#ex6
#           rgb(199,109,162,150,maxColorValue = 255),#ex5
#           rgb(129,184,223,150,maxColorValue = 255))#ex7
cols<- c("#AAD0E3","#277AB4","#FABF74","#FB7D1A","#CAB4D6","#B5DF90","#3AA12F","#693C9A")
DimPlot(Excitatory, reduction='umap',group.by="Anno_Idents", pt.size=1.5,label=T,label.size = 3,raster = F,
        label.box = T,cols = cols)
##meta.data添加信息
Excitatory@meta.data$color[Excitatory@meta.data$Anno_Idents %in% "DG_Ex1"] <- "#AAD0E3"
Excitatory@meta.data$color[Excitatory@meta.data$Anno_Idents %in% "DG_Ex2"] <- "#277AB4"
Excitatory@meta.data$color[Excitatory@meta.data$Anno_Idents %in% "CA1_Ex"] <- "#FABF74"
Excitatory@meta.data$color[Excitatory@meta.data$Anno_Idents %in% "CA_Ex_GAPDH"] <- "#FB7D1A"
Excitatory@meta.data$color[Excitatory@meta.data$Anno_Idents %in% "CA3_Ex"] <- "#CAB4D6"
Excitatory@meta.data$color[Excitatory@meta.data$Anno_Idents %in% "DG_Ex3"] <- "#B5DF90"
Excitatory@meta.data$color[Excitatory@meta.data$Anno_Idents %in% "DG_Ex5"] <- "#3AA12F"
Excitatory@meta.data$color[Excitatory@meta.data$Anno_Idents %in% "DG_Ex4"] <- "#693C9A"

saveRDS(Excitatory,file = "./subcluster\\Excitatory\\anno\\Excitatory_sub.rds")
####*****marker展示----------------
#DG: "MAML2","SEMA5A"
#CA："SV2B"
##CA1: "PID1"
##CA1/3 : "TYRO3"
##CA3:"HS3ST4","SULF2"
###DG1:0.3.7.14.10.21:  "SLC47A1","COLEC12",#"SLC4A4","COL6A3"
###DG3: 13 "SERPINE1"
###DG5: 15  "SLC14A1","TNC"
###DG4: 20 "LHFPL3","TMEM132C"
###DG2: 1.4.6 "MAML2","SEMA5A"
Excitatory_sub <- readRDS(file = "./subcluster\\Excitatory\\anno\\Excitatory_sub.rds")
cols<- c("#AAD0E3","#277AB4","#FABF74","#FB7D1A","#CAB4D6","#B5DF90","#3AA12F","#693C9A")
table(Excitatory_sub$color)
#c("#AAD0E3","#277AB4","#FABF74","#FB7D1A","#CAB4D6","#B5DF90","#3AA12F","#693C9A")
#c("DG_Ex1" ,"DG_Ex2" ,"CA1_Ex","CA_Ex_GAPDH","CA3_Ex" , "DG_Ex3","DG_Ex5","DG_Ex4")
DimPlot(Excitatory_sub, reduction='umap',group.by="Anno_Idents", pt.size=1.5,label=T,label.size = 3,raster = F,
        label.box = T,cols = cols)
table(Excitatory_sub$group)
# A     P 
# 5572 11066

DimPlot(Excitatory_sub, reduction='umap',group.by="Anno_Idents", pt.size=1,label=T,label.size = 3,raster = F,
        label.box = T,cols = cols)

marker <- c("MAML2","SEMA5A","SV2B","SLC47A1","COLEC12","SERPINE1","LHFPL3","TMEM132C","SLC14A1","TNC","PID1","TYRO3","HS3ST4","GAPDH")
library(ggplot2)
DotPlot(Excitatory_sub, features = marker,group.by = "Anno_Idents",cols = c("lightgrey", "#CC3333"))+coord_flip()+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size = 11,colour = "black",angle = 30),
        axis.text.y = element_text(size = 11, colour = "black"))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
ggsave("./subcluster\\Excitatory\\anno\\marker_dotplot.pdf")

a <- subset(Excitatory_sub,group %in% "A")
DimPlot(a, reduction='umap',group.by="Anno_Idents", pt.size=1,label=T,label.size = 3,raster = F,
        label.box = T,cols = cols)
#                A    P
# DG_Ex1        24 5559
# DG_Ex2      2031 1586
# CA1_Ex      1614 2088
# CA_Ex_GAPDH  612  281
# CA3_Ex       959  929
# DG_Ex3        89  322
# DG_Ex5       200  108
# DG_Ex4        43  193
###细胞比例----------
######相关性柱状图   纵轴log (percentage P/percentage A)----
names(table(Excitatory_sub$Anno_Idents))
Excitatory_sub$Anno_Idents <- factor(Excitatory_sub$Anno_Idents,levels = c("CA1_Ex","CA_Ex_GAPDH","CA3_Ex","DG_Ex1","DG_Ex2","DG_Ex3","DG_Ex4","DG_Ex5"))
# DG_Ex1      DG_Ex2      CA1_Ex CA_Ex_GAPDH      CA3_Ex      DG_Ex3      DG_Ex5      DG_Ex4 
# 5583        3617        3702         893        1888         411         308         236
all.an <- prop.table(table(Excitatory_sub$Anno_Idents[Excitatory_sub@meta.data$group %in% c("A")]))
All.po <- prop.table(table(Excitatory_sub$Anno_Idents[Excitatory_sub@meta.data$group %in% c("P")]))
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
ggsave("./subcluster\\Excitatory\\anno\\logPA.pdf",p1,height = 5,width = 5)


######细胞类型在前后端比例----
DG_Ex1 <- prop.table(table(Excitatory_sub$group[Excitatory_sub@meta.data$Anno_Idents %in% c("DG_Ex1")]))
DG_Ex2 <- prop.table(table(Excitatory_sub$group[Excitatory_sub@meta.data$Anno_Idents %in% c("DG_Ex2")]))
DG_Ex3 <- prop.table(table(Excitatory_sub$group[Excitatory_sub@meta.data$Anno_Idents %in% c("DG_Ex3")]))
DG_Ex4 <- prop.table(table(Excitatory_sub$group[Excitatory_sub@meta.data$Anno_Idents %in% c("DG_Ex4")]))
DG_Ex5 <- prop.table(table(Excitatory_sub$group[Excitatory_sub@meta.data$Anno_Idents %in% c("DG_Ex5")]))
CA1_Ex <- prop.table(table(Excitatory_sub$group[Excitatory_sub@meta.data$Anno_Idents %in% c("CA1_Ex")]))
CA3_Ex <- prop.table(table(Excitatory_sub$group[Excitatory_sub@meta.data$Anno_Idents %in% c("CA3_Ex")]))
CA_Ex_GAPDH <- prop.table(table(Excitatory_sub$group[Excitatory_sub@meta.data$Anno_Idents %in% c("CA_Ex_GAPDH")]))


cluster <- c(sort(rep(names(table(Excitatory_sub$Anno_Idents)),2)))
#Samples_1 <- factor(Samples_1,levels = names(table(scRNA@meta.data$idents)))
pos <- c(rep(rep(names(table(Excitatory_sub$group))),8))
library("ggplot2")
cell.prop<-as.data.frame(c(CA_Ex_GAPDH,CA1_Ex,CA3_Ex,DG_Ex1,DG_Ex2,DG_Ex3,DG_Ex4,DG_Ex5),pos)
cell.prop$pos <- rownames(cell.prop) 
cell.prop$cluster <- cluster
colnames(cell.prop)<-c("proportion","pos","cluster")
library('ggplot2')
library('reshape2')

p2 <- ggplot(cell.prop,aes(cluster,proportion,fill=pos))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+scale_fill_manual(values = c("#99CCCC","#FFCC99"))+
  theme(axis.text.x = element_text(angle=30))
ggsave("./subcluster\\Excitatory\\anno\\细胞类型比例.pdf",p2,height = 5,width = 5)

###统计细胞类型在前后端数目差异---------
DG_Ex1 <- prop.table(table(Excitatory_sub$orig.ident[Excitatory_sub@meta.data$Anno_Idents %in% c("DG_Ex1")]))
DG_Ex2 <- prop.table(table(Excitatory_sub$orig.ident[Excitatory_sub@meta.data$Anno_Idents %in% c("DG_Ex2")]))
DG_Ex3 <- prop.table(table(Excitatory_sub$orig.ident[Excitatory_sub@meta.data$Anno_Idents %in% c("DG_Ex3")]))
DG_Ex4 <- prop.table(table(Excitatory_sub$orig.ident[Excitatory_sub@meta.data$Anno_Idents %in% c("DG_Ex4")]))
DG_Ex5 <- prop.table(table(Excitatory_sub$orig.ident[Excitatory_sub@meta.data$Anno_Idents %in% c("DG_Ex5")]))
CA1_Ex <- prop.table(table(Excitatory_sub$orig.ident[Excitatory_sub@meta.data$Anno_Idents %in% c("CA1_Ex")]))
CA3_Ex <- prop.table(table(Excitatory_sub$orig.ident[Excitatory_sub@meta.data$Anno_Idents %in% c("CA3_Ex")]))
CA_Ex_GAPDH <- prop.table(table(Excitatory_sub$orig.ident[Excitatory_sub@meta.data$Anno_Idents %in% c("CA_Ex_GAPDH")]))
cluster <- c(sort(rep(names(table(Excitatory_sub$Anno_Idents)),10)))
#Samples_1 <- factor(Samples_1,levels = names(table(scRNA@meta.data$idents)))
pos <- rep(c(sort(rep(rep(names(table(Excitatory_sub$group))),5))),8)
library("ggplot2")
cell.prop<-as.data.frame(c(CA_Ex_GAPDH,CA1_Ex,CA3_Ex,DG_Ex1,DG_Ex2,DG_Ex3,DG_Ex4,DG_Ex5),pos)
cell.prop$pos <- rownames(cell.prop) 
cell.prop$cluster <- cluster
colnames(cell.prop)<-c("proportion","pos","cluster")
library('ggplot2')
library('reshape2')
ggplot(cell.prop, aes(fill=pos, y=proportion, x=cluster))+ 
  geom_bar(position="dodge",stat="identity",width=0.8)+
  theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size =8),
        panel.background = element_rect(color = 'black', fill = 'transparent'))
library('ggplot2')
library('ggpubr')
p3 <- ggplot(cell.prop, aes(fill=pos, y=proportion, x=cluster))+
  geom_bar(alpha=0.3,width=0.45,position=position_dodge(width=0.8),stat="identity",size=2)+
  theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size =8),
        panel.background = element_rect(color = 'black', fill = 'transparent'))+
  theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size =8),
        panel.background = element_rect(color = 'black', fill = 'transparent'))+
  theme(legend.direction = "horizontal", legend.position = "top")+
  labs(title = "", y="proportion", x = "")+
  theme(axis.text.x = element_text(size = 12,angle = 15,hjust=1, vjust=1))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title = element_text(size = 14))+
  geom_boxplot(width=0.45,
               position=position_dodge(width=0.8),
               size=0.2,outlier.colour = NA)+
  stat_compare_means(aes(group=pos),method = "t.test",
                     label="p.signif")+
  geom_jitter(data = cell.prop, aes(fill=pos, y=proportion, x=cluster),
              position = position_jitterdodge(0.2),size = 2, shape = 21, show.legend = FALSE)+
  scale_fill_manual(values = c("#99CCCC","#FFCC99"))

ggsave("./subcluster\\Excitatory\\anno\\箱线图.pdf",p3,height = 5,width = 5)
ggsave("./subcluster\\Excitatory\\anno\\组合图.pdf",p1|p2|p3,height = 5,width = 10)

##*抑制性神经元---------
scRNA <- readRDS(file = "./scRNA_anno_cluster.rds")
Inhibitory <- subset(scRNA,Anno_Idents %in% "Inhibitory neuron")

table(scRNA$group)
# A     P 
# 64895 66201
metadata <- Inhibitory@meta.data
dim(metadata)
metadata <- metadata[,-c(6:18,20)]
Inhibitory <- CreateSeuratObject(counts = Inhibitory@assays$RNA@counts,
                                 meta.data = metadata) 
dim(Inhibitory)  # 16956  6275
####
####
####
####**seurat流程--------
Inhibitory <- NormalizeData(Inhibitory, normalization.method = "LogNormalize", scale.factor = 10000)
Inhibitory <- FindVariableFeatures(Inhibitory, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Inhibitory)
Inhibitory <- ScaleData(Inhibitory, features = all.genes)
Inhibitory <- RunPCA(Inhibitory, features = VariableFeatures(object = Inhibitory))
DimPlot(Inhibitory, reduction = "pca")

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
#JackStrawPlot(Excitatory, dims = 1:15)
ElbowPlot(Inhibitory,ndims = 50)
Inhibitory <- FindNeighbors(Inhibitory, dims = 1:30)
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  Inhibitory <- FindClusters(Inhibitory, resolution = res)
}
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
Inhibitory <- RunUMAP(Inhibitory, dims = 1:30)
Inhibitory <- RunTSNE(Inhibitory, dims = 1:30)
library(clustree)
library(patchwork)
clustree(Inhibitory, prefix = 'hc_euclidean_') + coord_flip()
ggsave("./subcluster\\Inhibitory\\clustree.pdf",width = 13,height = 11)
d <- dist(Inhibitory@reductions[["pca"]]@cell.embeddings, method = "euclidean")
# Compute sample correlations
# 计算细胞之间的相关性
sample_cor <- cor(Matrix::t(Inhibitory@reductions[["pca"]]@cell.embeddings))

# Transform the scale from correlations
sample_cor <- (1 - sample_cor)/2

# Convert it to a distance object
d2 <- as.dist(sample_cor)
# euclidean
h_euclidean <- hclust(d, method = "ward.D2")

# correlation
h_correlation <- hclust(d2, method = "ward.D2")
#euclidean distance
Inhibitory$hc_euclidean_5 <- cutree(h_euclidean,k = 5)
Inhibitory$hc_euclidean_10 <- cutree(h_euclidean,k = 10)
Inhibitory$hc_euclidean_15 <- cutree(h_euclidean,k = 15)

#correlation distance
Inhibitory$hc_corelation_5 <- cutree(h_correlation,k = 5)
Inhibitory$hc_corelation_10 <- cutree(h_correlation,k = 10)
Inhibitory$hc_corelation_15 <- cutree(h_correlation,k = 15)
library(cowplot)
plot_grid(ncol = 3,
          DimPlot(Inhibitory, reduction = "umap", group.by = "hc_euclidean_5",cols = mycolors)+ggtitle("hc_euc_5"),
          DimPlot(Inhibitory, reduction = "umap", group.by = "hc_euclidean_10")+ggtitle("hc_euc_10"),
          DimPlot(Inhibitory, reduction = "umap", group.by = "hc_euclidean_15")+ggtitle("hc_euc_15"),
          
          DimPlot(Inhibitory, reduction = "umap", group.by = "hc_corelation_5",cols = mycolors)+ggtitle("hc_cor_5"),
          DimPlot(Inhibitory, reduction = "umap", group.by = "hc_corelation_10")+ggtitle("hc_cor_10"),
          DimPlot(Inhibitory, reduction = "umap", group.by = "hc_corelation_15")+ggtitle("hc_cor_15")
)
###不采用RNA_snn_res，难以注释，采用层次聚类注释---------
#DimPlot(Inhibitory, group.by = 'RNA_snn_res.0.8', pt.size = 1,label = T,reduction = "umap",label.box=T,cols = mycolors)
#Idents(Inhibitory) <- 'RNA_snn_res.0.8'  #25 cluster
Idents(Inhibitory) <- 'hc_corelation_5'  #5 cluster
DimPlot(Inhibitory, reduction = "umap", group.by = "hc_corelation_5",cols = mycolors,label.box=T,pt.size = 1,label = T)
ggsave("./subcluster\\Inhibitory\\hc_corelation_5.pdf",width = 13,height = 11)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
saveRDS(Inhibitory,"./subcluster\\Inhibitory\\Inhibitory-umap_tsne.rds")

Inhibitory <- readRDS("./subcluster\\Inhibitory\\Inhibitory-umap_tsne.rds")
levels(Inhibitory)
# DimPlot(Inhibitory, reduction='umap',group.by="RNA_snn_res.0.8", pt.size=3,label=T,label.size = 5,
#         cols = mycolors,label.box = T)
#ggsave("./subcluster\\Inhibitory\\RNA_snn_res.0.8.pdf",width = 13,height = 11)
#####**手动注释---------
Inhibitory.markers <- FindAllMarkers(Inhibitory, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Inhibitory.markers,"./subcluster\\Inhibitory\\Inhibitory.markers.txt",sep = "\t",quote = F)
Inhibitory.markers <- read.table("./subcluster\\Inhibitory\\Inhibitory.markers.txt",sep = "\t")
Inhibitory.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(Inhibitory, features = top5$gene) + scale_fill_gradientn(colors = c("#006699","#FFFFCC","#FF6600"))
BrainMarker <- read.table("./CellMarker_brain_embryo_marker.txt",sep = "\t",header = T)
names(table(BrainMarker$cell_name))
####参考文献 1. 参考A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation   2. Resolving cellular and molecular diversity along the hippocampal anterior-to-posterior axis in humans
PVALB <- c("ST18") #ST18 参考Evolutionarily conservative and non-conservative regulatory networks during primate interneuron development revealed by single-cell RNA and ATAC sequencing
SST <- 	c("SST")  
VIP <- c("VIP","CALB2") 
CCK <- c("CXCL14","CNR1")
LAMP5 <- c("LAMP5","SV2C") #采用
marker <- unique(c(PVALB,SST,VIP,CCK,LAMP5))
DotPlot(Inhibitory, features = c(marker),group.by = "hc_corelation_5",cols = c("lightgrey", "#CC3333"))+coord_flip()+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size = 11,colour = "black",angle = 30),
        axis.text.y = element_text(size = 11, colour = "black"))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
ggsave("./subcluster\\Inhibitory\\DotPlot_hc_corelation_5.pdf",width = 13,height = 11)

####**细胞注释完成-------------
# 1: PVALB   # 2:SST         #3:VIP    # 4:CCK  # 5: LAMP5

#*
levels(Inhibitory)
#齿状回dentate gyrus (DG)  
#CA
new.cluster.ids <- c("PVALB","SST","VIP","CCK","LAMP5")
unique(new.cluster.ids)
names(new.cluster.ids) <- levels(Inhibitory)
Inhibitory <- RenameIdents(Inhibitory, new.cluster.ids)
Inhibitory <- StashIdent(Inhibitory, save.name = 'Anno_Idents')
#cols<- c("#6699CC","#339999","#CCCC99","#FFCC00","#FF9900")
DimPlot(Inhibitory, reduction='umap',group.by="Anno_Idents", pt.size=1.5,label=T,label.size = 3,raster = F,
        label.box = T,cols = cols)
##meta.data添加信息
Inhibitory@meta.data$color[Inhibitory@meta.data$Anno_Idents %in% "PVALB"] <- "#6699CC"
Inhibitory@meta.data$color[Inhibitory@meta.data$Anno_Idents %in% "SST"] <- "#339999"
Inhibitory@meta.data$color[Inhibitory@meta.data$Anno_Idents %in% "VIP"] <- "#CCCC99"
Inhibitory@meta.data$color[Inhibitory@meta.data$Anno_Idents %in% "CCK"] <- "#FFCC00"
Inhibitory@meta.data$color[Inhibitory@meta.data$Anno_Idents %in% "LAMP5"] <- "#FF9900"

saveRDS(Inhibitory,file = "./subcluster\\Inhibitory\\anno\\Inhibitory_sub.rds")
####**marker展示----------------
PVALB <- c("ST18") #ST18 参考Evolutionarily conservative and non-conservative regulatory networks during primate interneuron development revealed by single-cell RNA and ATAC sequencing
SST <- 	c("SST")  
VIP <- c("VIP","CALB2") 
CCK <- c("CXCL14","CNR1")
LAMP5 <- c("LAMP5","SV2C") #采用
marker <- unique(c(PVALB,SST,VIP,CCK,LAMP5))
DotPlot(Inhibitory, features = c(marker),group.by = "Anno_Idents",cols = c("lightgrey", "#CC3333"))+coord_flip()+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size = 11,colour = "black",angle = 30),
        axis.text.y = element_text(size = 11, colour = "black"))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
ggsave("./subcluster\\Inhibitory\\anno\\marker_dotplot.pdf",width = 13,height = 11)

###细胞比例----------
######相关性柱状图   纵轴log (percentage P/percentage A)----
Inhibitory <-readRDS(file = "./subcluster\\Inhibitory\\anno\\Inhibitory_sub.rds")

names(table(Inhibitory$Anno_Idents))
all.an <- prop.table(table(Inhibitory$Anno_Idents[Inhibitory@meta.data$group %in% c("A")]))
All.po <- prop.table(table(Inhibitory$Anno_Idents[Inhibitory@meta.data$group %in% c("P")]))
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
ggsave("./subcluster\\Inhibitory\\anno\\logPA.pdf",height = 5,width = 5)


######细胞类型在前后端比例----
PVALB <- prop.table(table(Inhibitory$group[Inhibitory@meta.data$Anno_Idents %in% c("PVALB")]))
SST <- prop.table(table(Inhibitory$group[Inhibitory@meta.data$Anno_Idents %in% c("SST")]))
VIP <- prop.table(table(Inhibitory$group[Inhibitory@meta.data$Anno_Idents %in% c("VIP")]))
CCK <- prop.table(table(Inhibitory$group[Inhibitory@meta.data$Anno_Idents %in% c("CCK")]))
LAMP5 <- prop.table(table(Inhibitory$group[Inhibitory@meta.data$Anno_Idents %in% c("LAMP5")]))


cluster <- c(sort(rep(names(table(Inhibitory$Anno_Idents)),2)))
#Samples_1 <- factor(Samples_1,levels = names(table(scRNA@meta.data$idents)))
pos <- c(rep(rep(names(table(Inhibitory$group))),5))
library("ggplot2")
cell.prop<-as.data.frame(c(PVALB,SST,VIP,CCK,LAMP5),pos)
cell.prop$pos <- rownames(cell.prop) 
cell.prop$cluster <- cluster
colnames(cell.prop)<-c("proportion","pos","cluster")
library('ggplot2')
library('reshape2')

p2 <- ggplot(cell.prop,aes(cluster,proportion,fill=pos))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+scale_fill_manual(values = c("#99CCCC","#FFCC99"))+
  theme(axis.text.x = element_text(angle=30))
ggsave("./subcluster\\Inhibitory\\anno\\细胞类型比例.pdf",height = 5,width = 5)

###统计细胞类型在前后端数目差异---------
PVALB <- prop.table(table(Inhibitory$orig.ident[Inhibitory@meta.data$Anno_Idents %in% c("PVALB")]))
SST <- prop.table(table(Inhibitory$orig.ident[Inhibitory@meta.data$Anno_Idents %in% c("SST")]))
VIP <- prop.table(table(Inhibitory$orig.ident[Inhibitory@meta.data$Anno_Idents %in% c("VIP")]))
CCK <- prop.table(table(Inhibitory$orig.ident[Inhibitory@meta.data$Anno_Idents %in% c("CCK")]))
LAMP5 <- prop.table(table(Inhibitory$orig.ident[Inhibitory@meta.data$Anno_Idents %in% c("LAMP5")]))

cluster <- c(sort(rep(names(table(Inhibitory$Anno_Idents)),10)))
#Samples_1 <- factor(Samples_1,levels = names(table(scRNA@meta.data$idents)))
pos <- rep(c(sort(rep(rep(names(table(Inhibitory$group))),5))),5)
library("ggplot2")
cell.prop<-as.data.frame(c(PVALB,SST,VIP,CCK,LAMP5),pos)
cell.prop$pos <- rownames(cell.prop) 
cell.prop$cluster <- cluster
colnames(cell.prop)<-c("proportion","pos","cluster")
library('ggplot2')
library('reshape2')
ggplot(cell.prop, aes(fill=pos, y=proportion, x=cluster))+ 
  geom_bar(position="dodge",stat="identity",width=0.8)+
  theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size =8),
        panel.background = element_rect(color = 'black', fill = 'transparent'))


p3 <- ggplot(cell.prop, aes(fill=pos, y=proportion, x=cluster))+
  geom_bar(alpha=0.3,width=0.45,position=position_dodge(width=0.8),stat="identity",size=2)+
  theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size =8),
        panel.background = element_rect(color = 'black', fill = 'transparent'))+
  theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size =8),
        panel.background = element_rect(color = 'black', fill = 'transparent'))+
  theme(legend.direction = "horizontal", legend.position = "top")+
  labs(title = "", y="proportion", x = "")+
  theme(axis.text.x = element_text(size = 12,angle = 15,hjust=1, vjust=1))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title = element_text(size = 14))+
  geom_boxplot(width=0.45,
               position=position_dodge(width=0.8),
               size=0.2,outlier.colour = NA)+
  stat_compare_means(aes(group=pos),method = "t.test",
                     label="p.signif")+
  geom_jitter(data = cell.prop, aes(fill=pos, y=proportion, x=cluster),
              position = position_jitterdodge(0.2),size = 2, shape = 21, show.legend = FALSE)+
  scale_fill_manual(values = c("#99CCCC","#FFCC99"))

ggsave("./subcluster\\Inhibitory\\anno\\箱线图.pdf",height = 5,width = 5)
ggsave("./subcluster\\Inhibitory\\anno\\组合图.pdf",p1|p2|p3,height = 5,width = 10)

