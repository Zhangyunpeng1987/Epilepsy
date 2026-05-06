####*兴奋性神经元分析---------------------------
####****兴奋性神经元分亚型---------
####*****加载颜色-------
cols <- c(rgb(196,165,222,150,maxColorValue = 255),#L23
          rgb(129,184,223,150,maxColorValue = 255),#456
          rgb(246,202,229,150,maxColorValue = 255)) #L6 
####*****提取数据--------
Excitatory <- readRDS(file = "./temporal_lobe\\scRNA_anno_cluster.rds")
Excitatory <- subset(Excitatory, group %in% c("Nor","Ep"))
table(Excitatory$Anno_Idents)
Excitatory <- subset(Excitatory, Anno_Idents=="Excitatory neuron")  #67823
saveRDS(Excitatory,"./temporal_lobe\\subcluster\\Excitatory\\Excitatory.rds")

Excitatory <- readRDS("./temporal_lobe\\subcluster\\Excitatory\\Excitatory.rds")
metadata <- Excitatory@meta.data
metadata <- metadata[,-c(8:19)]
Excitatory <- CreateSeuratObject(counts = Excitatory@assays$RNA@counts,
                                 meta.data = metadata) 
dim(Excitatory)  #28306 67823
####
####
####
####*****seurat流程--------
Excitatory <- NormalizeData(Excitatory, normalization.method = "LogNormalize", scale.factor = 10000)
Excitatory <- FindVariableFeatures(Excitatory, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Excitatory), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Excitatory)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
all.genes <- rownames(Excitatory)
Excitatory <- ScaleData(Excitatory, features = all.genes)
Excitatory <- RunPCA(Excitatory, features = VariableFeatures(object = Excitatory))
DimPlot(Excitatory, reduction = "pca")
DimHeatmap(Excitatory, dims = 1, cells = 500, balanced = TRUE)
saveRDS(Excitatory,"./temporal_lobe\\subcluster\\Excitatory\\Excitatory-pca.rds")

Excitatory <- JackStraw(Excitatory, num.replicate = 100)
Excitatory <- ScoreJackStraw(Excitatory, dims = 1:20)
JackStrawPlot(Excitatory, dims = 1:15)
ElbowPlot(Excitatory,ndims = 50)
Excitatory <- FindNeighbors(Excitatory, dims = 1:20)
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  Excitatory <- FindClusters(Excitatory, resolution = res)
}

library(clustree)
library(patchwork)
p1 <- clustree(Excitatory, prefix = 'RNA_snn_res.') + coord_flip()
p2 <- DimPlot(Excitatory, group.by = 'RNA_snn_res.0.8', label = T)
p1 + p2 + plot_layout(widths = c(3, 1))
Idents(Excitatory) <- 'RNA_snn_res.0.8'
saveRDS(Excitatory,"./temporal_lobe\\subcluster\\Excitatory\\Excitatory-umap.rds")

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
Excitatory <- RunUMAP(Excitatory, dims = 1:30)
Excitatory <- RunTSNE(Excitatory, dims = 1:30)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
saveRDS(Excitatory,"./temporal_lobe\\subcluster\\Excitatory\\Excitatory-umap_tsne.rds")

Excitatory <- readRDS("./temporal_lobe\\subcluster\\Excitatory\\Excitatory-umap_tsne.rds")
levels(Excitatory)
DimPlot(Excitatory, reduction = "umap")
DimPlot(Excitatory, reduction='umap',group.by="RNA_snn_res.0.8", pt.size=1,label=T,label.size = 5,
        cols = mycolors)

DimPlot(Excitatory, reduction='tsne',group.by="RNA_snn_res.0.6", pt.size=1,label=T,label.size = 5,
        cols = mycolors)

L2_3_Cux2 <- c("LAMP5","CUX2")  #*可用
L4_L5_6_Rorb <- c("RORB","COBLL1")  #*可用
L5_6 <- c('TOX',"RXFP1")#*可用
L6 <- c("TLE4","NTNG2")  #*可用
all <- c("LAMP5","CUX2","RORB","FOXP2","RXFP1","TLE4","NTNG2")
StackedVlnPlot(Excitatory,all, pt.size=0)
FeaturePlot(Excitatory,all,cols = c("grey","#B40F20"))
#####****手动注释---------
Excitatory.markers <- FindAllMarkers(Excitatory, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Excitatory.markers,"./temporal_lobe\\subcluster\\Excitatory/resolution8\\Excitatory.markers.txt",sep = "\t",quote = F)
Excitatory.markers <- read.table("./temporal_lobe\\subcluster\\Excitatory/resolution8\\result\\Excitatory.markers.txt",sep = "\t")
Excitatory.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(Excitatory, features = top5$gene) +scale_fill_gradientn(colors = c("#006699","#FFFFCC","#FF6600"))
#####****将doheatmap数据提取出来用heatmap对簇进行聚类，观察是否有多个簇存在相同的趋势（存在即多个簇归为一类）--------
subobj <- subset(Excitatory, downsample = 300)
cluster <- data.frame(subobj$RNA_snn_res.0.8)
cluster$cell <- rownames(cluster)
colnames(cluster) <- c("cluster","cellid")
cluster <- cluster[order(cluster$cluster),]
raw.data <- as.matrix(subobj@assays$RNA@scale.data)
raw.data[1:4,1:4]
aa <- intersect(top5$gene,rownames(raw.data))
bb <- intersect( cluster$cellid,colnames(raw.data))
diff_exp <- raw.data[aa,bb]
diff_exp[1:4,1:4]
dim(diff_exp)  #151 67823
library(pheatmap)
p1 <- pheatmap(t(diff_exp),show_rownames = F, show_colnames = T,
               angle_col="90",annotation_names_row = F, annotation_names_col = TRUE,
               cluster_row = F,cluster_col = T,annotation_row = cluster,
               color = c("#006699","#FFFFCC","#FF6600"))
ggsave("./temporal_lobe\\subcluster\\Excitatory\\monocle3\\p1.pdf", plot = p1, width = 40, height = 40)

####
####
####
#####*****查看每个cluster之间是否存在共同表达趋势---------
setwd("./temporal_lobe\\subcluster\\Excitatory\\resolution8")
a <- c()
celltype.list <- as.character(unique(top5$cluster))
pdf("Excitatory_cluster_top5.pdf",width = 30)
for(i in 1:length(celltype.list)){
  marker <- top5$gene[top5$cluster %in% celltype.list[i]]
  name <- celltype.list[i]
  p <- print(paste(name, marker,sep = ":"))
  plots <- StackedVlnPlot(Excitatory, marker, pt.size=0, cols=c(my36colors,color_cols))
  print(plots)
}
dev.off()
DimPlot(Excitatory, reduction = "umap",group.by = "RNA_snn_res.0.8",label=T,cols = mycolors,label.box = T)
setwd("./temporal_lobe\\subcluster\\Excitatory\\resolution8")
a <- c()
celltype.list <- as.character(unique(top30$cluster))
for(i in 1:length(celltype.list)){
  #type <- paste(celltype.list[[i]], collapse = "_")
  marker <- top30$gene[top30$cluster %in% celltype.list[i]]
  if(length(marker)==0){
    name <- celltype.list[[i]]
    print(paste(name, "Have no marker",sep = ":"))
    next;
  }
  name <- celltype.list[i]
  p <- print(paste(name, marker,sep = ":"))
  plots <- StackedVlnPlot(Excitatory, marker, pt.size=0, cols=c(my36colors,color_cols,MYCOLOR))
  ggsave(plots,file=paste(name,".pdf",sep="_"),height = length(marker),width = 20,limitsize = FALSE)
  #print(plots)
}

####*****细胞注释完成-------------
# 0:Ex1 (L2_3_CUX2)  # 1:Ex1   #2:Ex1     # 3:Ex1 # 4:Ex1 # 5:Ex3(L5_6_RXFP1) 
# 6:Ex1   # 7:Ex1    # 8:Ex1   #9:Ex2(L4_RORB)   #10:Ex1
# 11:Ex1  # 12:Ex3  # 13:Ex1   # 14:Ex3   # 15:Ex4(L4_5_6)
# 16:Ex4  # 17:Ex1   # 18:Ex1  # 19:Ex3 # 20:Ex3
# 21:Ex6(L6_TLE4)    # 22:Ex1  # 23:Ex3   # 24:Ex6 # 25:Ex2
# 26:Ex6  # 27:Ex2   # 28:Ex5(L2_3_4_5_6)   # 29:Ex7(L6_NTNG2) # 30:Ex4
# 31:Ex3  # 32:Ex3   # 33:Ex2   # 34:Ex2    # 35:Ex2
# 36:Ex7  # 37:Ex3   # 38:Ex4   # 39:Ex4  # 40:Ex2
# 41:Ex6  # 42:Ex3   # 43:Ex6   # 44:Ex3 # 45:Ex4
# 46:Ex6  # 47:Ex3   # 48:Ex1   # 49:Ex3 # 50:Ex4
# 51:Ex1  # 52:Ex3   
#*
Idents(Excitatory) <- 'RNA_snn_res.0.8'
levels(Excitatory)
new.cluster.ids <- c("L2_3_Cux2","L2_3_Cux2","L2_3_Cux2","L2_3_Cux2","L2_3_Cux2","L4_L5_6_Rorb",
                     "L2_3_Cux2","L2_3_Cux2","L2_3_Cux2","L4_L5_6_Rorb","L2_3_Cux2",
                     "L2_3_Cux2","L4_L5_6_Rorb","L2_3_Cux2","L4_L5_6_Rorb","L4_L5_6_Rorb",
                     "L4_L5_6_Rorb","L2_3_Cux2", "L2_3_Cux2","L4_L5_6_Rorb","L4_L5_6_Rorb",
                     "L6","L2_3_Cux2", "L4_L5_6_Rorb","L6","L4_L5_6_Rorb",
                     "L6","L4_L5_6_Rorb", "L4_L5_6_Rorb","L6","L4_L5_6_Rorb",
                     "L4_L5_6_Rorb","L4_L5_6_Rorb", "L4_L5_6_Rorb","L4_L5_6_Rorb","L4_L5_6_Rorb",
                     "L6","L4_L5_6_Rorb", "L4_L5_6_Rorb","L4_L5_6_Rorb","L4_L5_6_Rorb",
                     "L6","L4_L5_6_Rorb", "L6","L4_L5_6_Rorb","L4_L5_6_Rorb",
                     "L6","L4_L5_6_Rorb", "L2_3_Cux2","L4_L5_6_Rorb","L4_L5_6_Rorb",
                     "L2_3_Cux2","L4_L5_6_Rorb")
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

DimPlot(Excitatory, reduction='umap',group.by="Anno_Idents", pt.size=3,label=T,label.size = 5,raster = F,
        cols = cols)
saveRDS(Excitatory,file = "./temporal_lobe\\subcluster\\Excitatory\\resolution8\\Excitatory_anno_cluster_1.rds")
####*****marker展示----------------
Ex1 <- c("LAMP5","CUX2")  #*可用L2_3_Cux2
Ex2 <- c("RORB","COBLL1")  #*可用L4_Rorb
Ex3 <- c('TOX',"RXFP1")#*可用L5_6
Ex4 <- c("RORB","COBLL1","TOX","RXFP1")#*可用L4_5_6
Ex5 <- c("LAMP5","CUX2","RORB","COBLL1","TOX","RXFP1")#*可用L2_3_4_5_6
Ex6_7 <- c("TLE4","NTNG2")  #*可用L6
all <- c("LAMP5","CUX2","RORB","COBLL1","TOX","RXFP1","TLE4","NTNG2")
##可以将通过rgb颜色改个透明度   去掉框顺便改个颜色
#print(FeaturePlot(object = scRNA, features = Excitatory_neuron,cols = c(rgb(220,212,213,180,maxColorValue = 255),rgb(174,27,52,50, maxColorValue = 255)),label.size = 6,pt.size = 1.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()))
FeaturePlot(Excitatory,features = all,cols = c("grey","#B40F20"))

###############
table(Excitatory$Anno_Idents)
# L2_3_Cux2 L4_L5_6_Rorb           L6 
# 34223        27219         6381 
table(Excitatory$group)
# Ep  NeuN   Nor 
# 27797     0 40026


##meta.data添加信息
metadata <- Excitatory@meta.data
group.id <- as.data.frame(metadata[,4])
rownames(group.id) <- rownames(metadata)
group.id$group <- t(as.data.frame(stringr::str_extract_all(group.id$`metadata[, 4]`, '\\D+')))
group <- factor(group.id[,2])
class(group)
Excitatory <- AddMetaData(Excitatory, group,col.name = "group")
dim(Excitatory)  #28306 67823

table(Excitatory$group)
# Ep   Nor 
# 27797 40026
saveRDS(Excitatory,file = "./temporal_lobe\\subcluster\\Excitatory\\resolution8\\Excitatory_anno_cluster_1.rds")
###细胞比例----------
######相关性柱状图   纵轴log (percentage P/percentage A)----
Excitatory <-readRDS(file = "./temporal_lobe\\subcluster\\Excitatory\\resolution8\\Excitatory_anno_cluster_1.rds")

names(table(Excitatory$Anno_Idents))
All.ep <- prop.table(table(Excitatory$Anno_Idents[Excitatory@meta.data$group %in% c("Ep")]))
All.Nor <- prop.table(table(Excitatory$Anno_Idents[Excitatory@meta.data$group %in% c("Nor")]))
df_prop <- cbind(All.ep,All.Nor)
df_prop <- t(apply(df_prop, 1, function(x) log(x / x[1])))
df_prop <- as.data.frame(df_prop)
df_prop$celltype <- rownames(df_prop)
df_prop <- df_prop[,-1]
class(df_prop)
colnames(df_prop) <- c("log","celltype")
df_prop$celltype <- factor(df_prop$celltype,levels = c(names(table(Excitatory$Anno_Idents))))

p1 <- ggplot(df_prop, aes(x=celltype,y=log))+ 
  geom_bar(position="dodge",stat="identity",width=0.8)+
  theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size =8),
        panel.background = element_rect(color = 'black', fill = 'transparent'))
ggsave("./temporal_lobe\\subcluster\\Excitatory\\resolution8\\result\\logPA.pdf",height = 5,width = 5)


######细胞类型在前后端比例----
L2_3_Cux2 <- prop.table(table(Excitatory$group[Excitatory@meta.data$Anno_Idents %in% c("L2_3_Cux2")]))
L4_L5_6_Rorb <- prop.table(table(Excitatory$group[Excitatory@meta.data$Anno_Idents %in% c("L4_L5_6_Rorb")]))
L6 <- prop.table(table(Excitatory$group[Excitatory@meta.data$Anno_Idents %in% c("L6")]))

cluster <- c(sort(rep(names(table(Excitatory$Anno_Idents)),2)))
#Samples_1 <- factor(Samples_1,levels = names(table(scRNA@meta.data$idents)))
pos <- c(rep(rep(names(table(Excitatory$group))),3))
library("ggplot2")
cell.prop<-as.data.frame(c(L2_3_Cux2,L4_L5_6_Rorb,L6),pos)
cell.prop$pos <- rownames(cell.prop) 
cell.prop$cluster <- cluster
colnames(cell.prop)<-c("proportion","pos","cluster")
library('ggplot2')
library('reshape2')
p2 <- ggplot(cell.prop,aes(cluster,proportion,fill=pos))+
  geom_bar(stat="identity",position="fill",alpha=0.7)+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+scale_fill_manual(values = c('#c82423','#2878b5'))+
  theme(axis.text.x = element_text(angle=30))

###统计细胞类型在前后端数目差异---------
L2_3_Cux2 <- prop.table(table(Excitatory$idents[Excitatory@meta.data$Anno_Idents %in% c("L2_3_Cux2")]))
L4_L5_6_Rorb <- prop.table(table(Excitatory$idents[Excitatory@meta.data$Anno_Idents %in% c("L4_L5_6_Rorb")]))
L6 <- prop.table(table(Excitatory$idents[Excitatory@meta.data$Anno_Idents %in% c("L6")]))

cluster <- c(sort(rep(names(table(Excitatory$Anno_Idents))[1:2],19)),rep(names(table(Excitatory$Anno_Idents))[3],18))

pos <- c(rep(c(rep(names(table(Excitatory$group))[1],9),rep(names(table(Excitatory$group))[2],10)),2),
         rep(names(table(Excitatory$group))[1],9),rep(names(table(Excitatory$group))[2],9))
library("ggplot2")
cell.prop<-as.data.frame(c(L2_3_Cux2,L4_L5_6_Rorb,L6),pos)
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
  geom_boxplot(width=0.45,alpha=0.2,
               position=position_dodge(width=0.8),
               size=0.2,outlier.colour = NA)+
  stat_compare_means(aes(group=pos),method = "t.test",
                     label="p.signif")+
  geom_jitter(data = cell.prop, aes(fill=pos, y=proportion, x=cluster),
              position = position_jitterdodge(0.2),size = 2, shape = 21, show.legend = FALSE)+
  scale_fill_manual(values = c('#c82423','#2878b5'))

ggsave("./temporal_lobe\\subcluster\\Excitatory\\resolution8\\result\\boxplot.pdf",height = 5,width = 5)
ggsave("./temporal_lobe\\subcluster\\Excitatory\\resolution8\\result\\combind.pdf",p1|p2|p3,height = 5,width = 10)
####*抑制性神经元亚型分析------------
###****颜色设置-----------
cols_2<-c("#FFBE7A","#82B0D2","#8ECFC9","#D76364")  #抑制性神经元亚型颜色
cols_1 <- c("#CC9966","#CCCC66","#669999")  #####抑制性神经元monocle2状态颜色
####****抑制性神经元分亚型---------
####*****提取数据--------
Inhibitory <- readRDS(file = "./temporal_lobe\\scRNA_anno_cluster.rds")
Inhibitory <- subset(Inhibitory, group %in% c("Nor","Ep"))
table(Inhibitory$Anno_Idents)
Inhibitory <- subset(Inhibitory, Anno_Idents=="Inhibitory neuron")  #39612
saveRDS(Inhibitory,"./temporal_lobe\\subcluster\\Inhibitory\\Inhibitory.rds")

Inhibitory <- readRDS("./temporal_lobe\\subcluster\\Inhibitory\\Inhibitory.rds")
metadata <- Inhibitory@meta.data
metadata <- metadata[,-c(8:19,21)]
Inhibitory <- CreateSeuratObject(counts = Inhibitory@assays$RNA@counts,
                                 meta.data = metadata) 

##meta.data添加信息
group.id <- as.data.frame(metadata[,4])
rownames(group.id) <- rownames(metadata)
group.id$group <- t(as.data.frame(stringr::str_extract_all(group.id$`metadata[, 4]`, '\\D+')))
group <- factor(group.id[,2])
class(group)
Inhibitory <- AddMetaData(Inhibitory, group,col.name = "group")
dim(Inhibitory)  #28306 39612
####
####
####
####*****seurat流程--------
Inhibitory <- NormalizeData(Inhibitory, normalization.method = "LogNormalize", scale.factor = 10000)
Inhibitory <- FindVariableFeatures(Inhibitory, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Inhibitory), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Inhibitory)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
all.genes <- rownames(Inhibitory)
Inhibitory <- ScaleData(Inhibitory, features = all.genes)
Inhibitory <- RunPCA(Inhibitory, features = VariableFeatures(object = Inhibitory))
DimPlot(Inhibitory, reduction = "pca")
DimHeatmap(Inhibitory, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Inhibitory <- JackStraw(Inhibitory, num.replicate = 100)
Inhibitory <- ScoreJackStraw(Inhibitory, dims = 1:20)
JackStrawPlot(Inhibitory, dims = 1:15)
ElbowPlot(Inhibitory,ndims = 50)
Inhibitory <- FindNeighbors(Inhibitory, dims = 1:30)
#Inhibitory <- FindClusters(Inhibitory, resolution = 0.5)
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  Inhibitory <- FindClusters(Inhibitory, resolution = res)
}

library(clustree)
library(patchwork)
p1 <- clustree(Inhibitory, prefix = 'RNA_snn_res.') + coord_flip()
p2 <- DimPlot(Inhibitory, group.by = 'RNA_snn_res.0.4', label = T)
p1 + p2 + plot_layout(widths = c(3, 1))

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
Inhibitory <- RunUMAP(Inhibitory, dims = 1:30)
Inhibitory <- RunTSNE(Inhibitory, dims = 1:30)
Idents(Inhibitory) <- 'RNA_snn_res.0.4'

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
saveRDS(Inhibitory,"./temporal_lobe\\subcluster\\Inhibitory\\Inhibitory-umap.rds")
Inhibitory <- readRDS("./temporal_lobe\\subcluster\\Inhibitory\\Inhibitory-umap.rds")
levels(Inhibitory)
DimPlot(Inhibitory, reduction = "umap")
DimPlot(Inhibitory, reduction='umap',group.by="RNA_snn_res.0.4", pt.size=5,label=T,label.size = 5,raster = F,
        cols = mycolors)

DimPlot(Inhibitory, reduction='tsne',group.by="RNA_snn_res.0.4", pt.size=5,label=T,label.size = 5,raster = F,
        cols = mycolors)

#####****手动注释---------
Inhibitory.markers <- FindAllMarkers(Inhibitory, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Inhibitory.markers,"./temporal_lobe\\subcluster\\Inhibitory/Inhibitory.markers.txt",sep = "\t",quote = F)
# Inhibitory.markers <- read.table("./temporal_lobe\\subcluster\\Inhibitory/Inhibitory.markers.txt",sep = "\t")
Inhibitory.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(Inhibitory, features = top5$gene) +scale_fill_gradientn(colors = c("#006699","#FFFFCC","#FF6600"))




BrainMarker <- read.table("E:\\original\\epilepsy_Data\\CellMarker_brain_embryo_marker.txt",sep = "\t",header = T)
names(table(BrainMarker$cell_name))
Inhibitory.markers <- read.table("./temporal_lobe\\subcluster\\Inhibitory/Inhibitory.markers.txt",sep = "\t")

####
####
####
Inhibitory.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top50

marker <- intersect(Inhibitory.markers$gene[Inhibitory.markers$cluster %in% c("3")],Inhibitory.markers$gene[Inhibitory.markers$cluster %in% c("7")])
StackedVlnPlot(Inhibitory, pv, pt.size=0, cols=my36colors)
setwd("./temporal_lobe\\subcluster\\Inhibitory")
# pdf("E:\\original\\epilepsy_Data\\Epilepsy19-master\\FIRST/celltype_marker_scRna.pdf")
a <- c()
celltype.list <- as.character(unique(top50$cluster))
for(i in 1:length(celltype.list)){
  #type <- paste(celltype.list[[i]], collapse = "_")
  marker <- top50$gene[top50$cluster %in% celltype.list[i]]
  if(length(marker)==0){
    name <- celltype.list[[i]]
    print(paste(name, "Have no marker",sep = ":"))
    next;
  }
  name <- celltype.list[i]
  p <- print(paste(name, marker,sep = ":"))
  plots <- StackedVlnPlot(Inhibitory, marker, pt.size=0, cols=my36colors)
  ggsave(plots,file=paste(name,".pdf",sep="_"),height = length(marker),width = 20,limitsize = FALSE)
  #print(plots)
}
DimPlot(Inhibitory, reduction = "umap",group.by = "RNA_snn_res.0.4",label=T,cols = mycolors,label.box = T)
####*****细胞注释完成-------------
# 0:In_Vip  # 1:In_Sst   #2:In_Vip # 3:In_Pvalb # 4:In_Vip # 5:In_Sst # 6:In_Sst # 7:In_Pvalb
# 8:In_Vip  #9:In_Lamp5   #10:In_Sst
# 11:In_Lamp5# 12:In_Lamp5# 13:In_Pvalb # 14:In_Vip # 15:In_Lamp5
# 16:In_Pvalb# 17:In_Lamp5# 18:In_Vip# 19:In_Sst# 20:In_Vip
#*pv <- c("PVALB","NOS1","SULF1","LHX6",'KCNS3',"CRH","PLEKHH2",'LGR5') "PVALB""SULF1""KCNS3" "RHOBTB3""GHR""SLC4A4""PTPN13""RUNX2" "PPARGC1A" "EPB41" "DPP10-AS3" 7,13,16  "CNTNAP3B"   "CNTNAP3"
#*Sst <- 	c("SST","NOS1","SEMA6A","FAM89A","LHX6","GRIK1")  1,5,6,10,19
#*Vip <- c("VIP","TAC3","CALB2","LAMA3","FAM19A1","NPR3")  0,2,4,8,14,18,20 ARPP21 CNR1  CHRNA7 
#*Lamp5 <- c("LAMP5","SV2C")  9.11.12.15.17
#*
Idents(Inhibitory) <- 'RNA_snn_res.0.4'
levels(Inhibitory)
new.cluster.ids <- c("In_Vip","In_Sst","In_Vip","In_Pvalb","In_Vip","In_Sst",
                     "In_Sst","In_Pvalb","In_Vip","In_Lamp5","In_Sst",
                     "In_Lamp5","In_Lamp5","In_Pvalb","In_Vip","In_Lamp5",
                     "In_Pvalb","In_Lamp5", "In_Vip","In_Sst","In_Vip")
names(new.cluster.ids) <- levels(Inhibitory)
Inhibitory <- RenameIdents(Inhibitory, new.cluster.ids)
Inhibitory <- StashIdent(Inhibitory, save.name = 'Anno_Idents')
cols <- c(rgb(227,174,60,150,maxColorValue = 255),#CC99CC
          rgb(77,105,112,150,maxColorValue = 255), #99CCCC
          rgb(74,125,144,150,maxColorValue = 255),#669999
          rgb(194,79,26,150,maxColorValue = 255))#336699

DimPlot(Inhibitory, reduction='umap',group.by="Anno_Idents", pt.size=5,label=T,label.size = 5,raster = F,
        cols = cols_2)
#c(color_cols
saveRDS(Inhibitory,file = "./temporal_lobe\\subcluster\\Inhibitory\\Inhibitory_anno_cluster.rds")
Inhibitory <- readRDS(file = "./temporal_lobe\\subcluster\\Inhibitory\\Inhibitory_anno_cluster.rds")

####*****marker展示----------------
In_Lamp5 <- c("LAMP5","SV2C")
Vip <- c("VIP","ARPP21") 
Sst <- 	c("SST","PDE1A")
pv <- c("PVALB","CNTNAP3") 
##可以将通过rgb颜色改个透明度   去掉框顺便改个颜色
#print(FeaturePlot(object = scRNA, features = Excitatory_neuron,cols = c(rgb(220,212,213,180,maxColorValue = 255),rgb(174,27,52,50, maxColorValue = 255)),label.size = 6,pt.size = 1.5) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()))
FeaturePlot(Inhibitory,features = pv,cols = c("grey","#B40F20"))
###细胞比例----------
######相关性柱状图   纵轴log (percentage P/percentage A)----
Inhibitory <- readRDS(file = "./temporal_lobe\\subcluster\\Inhibitory\\Inhibitory_anno_cluster.rds")

ep <- subset(Inhibitory,group %in% "Ep")
nor <- subset(Inhibitory,group %in% "Nor")
table(ep$Anno_Idents)

table(nor$Anno_Idents)

names(table(Inhibitory$Anno_Idents))
All.ep <- prop.table(table(Inhibitory$Anno_Idents[Inhibitory@meta.data$group %in% c("Ep")]))
All.Nor <- prop.table(table(Inhibitory$Anno_Idents[Inhibitory@meta.data$group %in% c("Nor")]))
df_prop <- cbind(All.ep,All.Nor)
df_prop <- t(apply(df_prop, 1, function(x) log(x / x[1])))
df_prop <- as.data.frame(df_prop)
df_prop$celltype <- rownames(df_prop)
df_prop <- df_prop[,-1]
class(df_prop)
colnames(df_prop) <- c("log","celltype")
df_prop$celltype <- factor(df_prop$celltype,levels = c("In_Lamp5","In_Pvalb","In_Sst","In_Vip"))

p1 <- ggplot(df_prop, aes(x=celltype,y=log))+ 
  geom_bar(position="dodge",stat="identity",width=0.8)+
  theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size =8),
        panel.background = element_rect(color = 'black', fill = 'transparent'))
ggsave("./temporal_lobe\\subcluster\\Inhibitory\\logPA.pdf",height = 5,width = 5)


######细胞类型在前后端比例----
LAMP5 <- prop.table(table(Inhibitory$group[Inhibitory@meta.data$Anno_Idents %in% c("In_Lamp5")]))
PVALB <- prop.table(table(Inhibitory$group[Inhibitory@meta.data$Anno_Idents %in% c("In_Pvalb")]))
SST <- prop.table(table(Inhibitory$group[Inhibitory@meta.data$Anno_Idents %in% c("In_Sst")]))
VIP <- prop.table(table(Inhibitory$group[Inhibitory@meta.data$Anno_Idents %in% c("In_Vip")]))


cluster <- c(sort(rep(names(table(Inhibitory$Anno_Idents)),2)))
#Samples_1 <- factor(Samples_1,levels = names(table(scRNA@meta.data$idents)))
pos <- c(rep(rep(names(table(Inhibitory$group))),4))
library("ggplot2")
cell.prop<-as.data.frame(c(LAMP5,PVALB,SST,VIP),pos)
cell.prop$pos <- rownames(cell.prop) 
cell.prop$cluster <- cluster
colnames(cell.prop)<-c("proportion","pos","cluster")
library('ggplot2')
library('reshape2')

p2 <- ggplot(cell.prop,aes(cluster,proportion,fill=pos))+
  geom_bar(stat="identity",position="fill",alpha=0.7)+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+scale_fill_manual(values = c('#c82423','#2878b5'))+
  theme(axis.text.x = element_text(angle=30))
ggsave("./temporal_lobe\\subcluster\\Inhibitory\\细胞类型比例.pdf",height = 5,width = 5)

###统计细胞类型在前后端数目差异---------
LAMP5 <- prop.table(table(Inhibitory$idents[Inhibitory@meta.data$Anno_Idents %in% c("In_Lamp5")]))
PVALB <- prop.table(table(Inhibitory$idents[Inhibitory@meta.data$Anno_Idents %in% c("In_Pvalb")]))
SST <- prop.table(table(Inhibitory$idents[Inhibitory@meta.data$Anno_Idents %in% c("In_Sst")]))
VIP <- prop.table(table(Inhibitory$idents[Inhibitory@meta.data$Anno_Idents %in% c("In_Vip")]))

cluster <- c(sort(rep(names(table(Inhibitory$Anno_Idents)),19)))
#Samples_1 <- factor(Samples_1,levels = names(table(scRNA@meta.data$idents)))
pos <- c(rep(c(rep(names(table(Inhibitory$group))[1],9),rep(names(table(Inhibitory$group))[2],10)),4))
library("ggplot2")
cell.prop<-as.data.frame(c(LAMP5,PVALB,SST,VIP),pos)
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
  geom_boxplot(width=0.45,alpha=0.2,
               position=position_dodge(width=0.8),
               size=0.2,outlier.colour = NA)+
  stat_compare_means(aes(group=pos),method = "t.test",
                     label="p.signif")+
  geom_jitter(data = cell.prop, aes(fill=pos, y=proportion, x=cluster),
              position = position_jitterdodge(0.2),size = 2, shape = 21, show.legend = FALSE)+
  scale_fill_manual(values = c('#c82423','#2878b5'))

ggsave("./temporal_lobe\\subcluster\\Inhibitory\\boxplot.pdf",height = 5,width = 5)
ggsave("./temporal_lobe\\subcluster\\Inhibitory\\combind.pdf",p1|p2|p3,height = 5,width = 10)
