###*GSE160189-----
###**读取数据------------
library(data.table)
Hippo_snRNA<-fread(input="./GSE160189_Hippo_Counts.csv", stringsAsFactors = T,header = T)
dim(Hippo_snRNA)  #17180 131326
class(Hippo_snRNA)
head(Hippo_snRNA[1:4,1:4])
sample_1 <- str_split(colnames(Hippo_snRNA),'_',simplify = T)[,1]
sample_1 <- unique(sample_1)
sample_1
Hippo_snRNA <- as.data.frame(Hippo_snRNA)
rownames(Hippo_snRNA) <- Hippo_snRNA$gene
Hippo_snRNA <- Hippo_snRNA[,-1]
head(Hippo_snRNA[1:4,1:4])
###**A--前海马体   P--后海马体  各五个样本---------
Hippo_snRNA_1 <- CreateSeuratObject(counts = Hippo_snRNA,min.cells = 3, min.features = 200)
###**计算线粒体和红细胞基因比例----------
Hippo_snRNA_1[["percent.mt"]] <- PercentageFeatureSet(Hippo_snRNA_1, pattern = "^MT-")
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(Hippo_snRNA_1@assays$RNA)) 
HB.genes <- rownames(Hippo_snRNA_1@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
Hippo_snRNA_1[["percent.HB"]]<-PercentageFeatureSet(Hippo_snRNA_1, features=HB.genes) 

###**质控标准------------
###质控小提琴图
##设置可能用到的主题
theme.set2 <- theme(axis.title.x = element_blank())
##设置绘图元素
plot.features <- c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB")
#质控前小提琴图
plots <- list()
for (i in 1:length(plot.features)) {
  plots[[i]] <- VlnPlot(Hippo_snRNA_1, group.by = "orig.ident", pt.size = 0.1,
                        features = plot.features[i])+theme.set2+NoLegend()
}
violin <- wrap_plots(plots = plots,nrow = 2)
ggsave("./vlnplot_before_qc.pdf", plot = violin, width = 12, height = 12) 
plot2 <- FeatureScatter(Hippo_snRNA_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("./pearplot_before_qc.pdf", plot = plot2, width = 5, height = 5) 

All.scRNA <- subset(Hippo_snRNA_1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
dim(All.scRNA)  #16956 131096
col.num <- length(levels(as.factor(All.scRNA@meta.data$orig.ident)))
violin <-VlnPlot(All.scRNA, group.by = "orig.ident",
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("./vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 

####**查看批次效应------------
All.scRNA <- NormalizeData(All.scRNA) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData()
scRNA <- RunPCA(All.scRNA, verbose = F)
ElbowPlot(scRNA,ndims = 50)
DimPlot(scRNA, reduction = "pca",cols = mycolors,group.by = "orig.ident",raster = F)
scRNA <- scRNA %>% RunUMAP(dims = 1:25)
scRNA <- FindNeighbors(scRNA, dims = 1:25) %>% FindClusters(resolution = 0.6)
pdf("./batch_before_cluster.pdf")
for (i in c(0.2, 0.3, 0.4, 0.5,0.6,0.7, 0.8, 1, 1.2, 1.5, 2)) {
  scRNA <- FindClusters(scRNA, resolution = i)
  p1 <- DimPlot(scRNA, reduction = "umap",group.by = "orig.ident",split.by = "orig.ident",ncol = 2) + labs(title = paste0("resolution: ", i))
  p2 <- DimPlot(scRNA, reduction = "tsne",group.by="orig.ident", pt.size=0.3,label=T) + labs(title = paste0("resolution: ", i))
  p3 <- DimPlot(scRNA, reduction = "umap",group.by = "orig.ident",pt.size=0.3,label=T)
  p4 <- DimPlot(scRNA, reduction = "tsne",label=T)
  p5 <- DimPlot(scRNA, reduction = "umap",label=T)
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
}
dev.off()
saveRDS(scRNA,"./GSE160189.scRNA_cluster.rds")
####不需要重新去除批次效应---
#####
#####
#####
####**手动注释----------
scRNA <- readRDS(file = "./GSE160189.scRNA_cluster.rds")
levels(scRNA)
library(clustree)
library(patchwork)
clustree(scRNA, prefix = 'RNA_snn_res.') + coord_flip()
DimPlot(scRNA, group.by = 'RNA_snn_res.0.6', label = T)

scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(scRNA.markers,"./scRNA.markers.txt",sep = "\t",quote = F)
scRNA.markers <- read.table("./scRNA.markers.txt",sep = "\t")
scRNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(scRNA, features = top5$gene) +scale_fill_gradientn(colors = c("#006699","#FFFFCC","#FF6600"))

BrainMarker <- read.table("./CellMarker_brain_embryo_marker.txt",sep = "\t",header = T)
names(table(BrainMarker$cell_name))
####
####
####
#**小提琴图展示findallmarker基因与文献中marker的交集基因
celltype.list <- c("A1 astrocyte","A2 astrocyte","B cell" , "Astrocyte","Adult neuronal stem cell","Broad excitatory neuron","Cajal–Retzius cell","Deep layer excitatory neuron",
                   "Dendritic cell","Dopaminergic neuron","Endothelial cell","Epithelial cell","Excitatory neuron", "GABAergic neuron","Glial cell","Glial neural progenitor cell",
                   "Glutamatergic neuron","Immature neuron" ,"Inhibitory neuron","Intermediate progenitor cell","Interneuron","Interneuron precursor",
                   "Late neuronal progenitor" ,"Layer 4 cell","Deep layer cell","Lower layer neuron","M1 macrophage","M1 microglial cell","M2 macrophage" ,
                   "Macrophage" ,"Mast cell","Mature oligodendrocyte","Microglial cell","Natural killer cell","Neural stem cell","Neuron","Neuronal precursor cell",
                   "Neuronal progenitor cell","Oligodendrocyte","Oligodendrocyte precursor cell","Oligodendrocyte progenitor cell","Outer radial glia","Pericyte",
                   "Proneural cell","Radial glial cell","Sst interneuron","Superficial layer excitatory neuron","T cell","T follicular helper(Tfh) cell","Upper layer cell",
                   "von Economo neuron(VEN)")
setwd(".\\cellmarker_marker")
a <- c()
for(i in 1:length(celltype.list)){
  #type <- paste(celltype.list[[i]], collapse = "_")
  marker <- intersect(scRNA.markers$gene,BrainMarker[BrainMarker$cell_name %in% celltype.list[[i]],2])
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
  plots <- StackedVlnPlot(scRNA, marker, pt.size=0, cols=my36colors)
  ggsave(plots,file=paste(type,".pdf",sep="_"),height = length(marker),width = 20,limitsize = FALSE)
  #print(plots)
}
DimPlot(scRNA, reduction = "umap",label=T,cols = mycolors[36:63],label.box = T)


####**细胞注释完成-------------
# 0:Oli  # 1:Oli   # 2:Oli  # 3:Oli  # 4:opc   # 5:Oli# 6:micro   # 7:astro# 8:micro# 9:astro# 10:ex
# 11:opc# 12:ex# 13:ex# 14:In# 15:micro# 16:micro# 17:ex# 18:In# 19:endo# 20:ex
# 21:ex # 22:In# 23:ex# 24:opc# 25:astro# 26:astro
#Idents(scRNA) <- 'RNA_snn_res.0.6'

levels(scRNA)
new.cluster.ids <- c("Oligodendrocyte","Oligodendrocyte","Oligodendrocyte","Oligodendrocyte","Oligodendrocyte progenitor cell","Oligodendrocyte",
                     "Microglial","Astrocyte","Microglial","Astrocyte","Excitatory neuron",
                     "Oligodendrocyte progenitor cell","Excitatory neuron","Excitatory neuron","Inhibitory neuron","Microglial",
                     "Microglial","Excitatory neuron", "Inhibitory neuron","Endothelial",
                     "Excitatory neuron","Excitatory neuron","Inhibitory neuron","Excitatory neuron","Oligodendrocyte progenitor cell",
                     "Astrocyte","Astrocyte")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
scRNA <- StashIdent(scRNA, save.name = 'Anno_Idents')
DimPlot(scRNA, reduction='umap', group.by="Anno_Idents", pt.size=0.3,label=T,label.size = 5,cols = c("#CCCCFF","#FFCC00","#0099CC","#99CC99","#CC6699","#FF9966","#CC9999"))

#c("#CCCCFF",         "#FFCC00",                        "#0099CC",        "#99CC99",  "#CC6699",  "#FF9966", "#CC9999")
# Oligodendrocyte    Oligodendrocyte_progenitor_cell    Microglial_cell  Astrocyte     ex           in        endo     
saveRDS(scRNA,file = "E./scRNA_anno_cluster.rds")

a <- subset(scRNA,group %in% "A")
DimPlot(a, reduction='umap', group.by="group", raster = F,pt.size=0.3,label=F,label.size = 5,cols = c("#99CCCC"))
ggsave("./A-umap1.pdf",height = 11,width = 11)

p<- subset(scRNA,group %in% "P")
DimPlot(p, reduction='umap', group.by="group", pt.size=0.3,raster = F,label=F,label.size = 5,cols = c("#FFCC99"))

ggsave("./P_umap1.pdf",height = 11,width = 11)


#**展示注释所用marker-------
#*
#*
#*
scRNA <- readRDS(file = "./scRNA_anno_cluster.rds")
scRNA.markers <- read.table("//scRNA.markers.txt",sep = "\t")
celltype.list <- unique(c("Oligodendrocyte","Oligodendrocyte","Oligodendrocyte","Oligodendrocyte","Oligodendrocyte progenitor cell","Oligodendrocyte",
                          "Microglial","Astrocyte","Microglial","Astrocyte","Excitatory neuron",
                          "Oligodendrocyte progenitor cell","Excitatory neuron","Excitatory neuron","Inhibitory neuron","Microglial",
                          "Microglial","Excitatory neuron", "Inhibitory neuron","Endothelial",
                          "Excitatory neuron","Excitatory neuron","Inhibitory neuron","Excitatory neuron","Oligodendrocyte progenitor cell",
                          "Astrocyte","Astrocyte"))
marker <- intersect(scRNA.markers$gene,BrainMarker[BrainMarker$cell_name %in% "Endothelial cell",2])
marker
StackedVlnPlot(scRNA, In_Sst, pt.size=0, cols=my36colors)
Astrocyte <- c("GJA1","AQP4","SLC14A1","RFX4")
Endothelial<- c("ABCB1","CLDN5","EBF1","APOLD1")
Excitatory_neuron<- c("STMN2","NRGN","SLC17A7")
Inhibitory_neuron<- c("GAD1","RELN","CCK","GAD2")
Oligodendrocyte_progenitor_cell<- c("BCAN","CA10","COL9A1","LIMA1")
Microglial_cell<- c("CD74","APBB1IP","TBXAS1","DOCK8")
Oligodendrocyte<- c("MOG","CNP","MAG","CNDP1")
FeaturePlot(scRNA,features = Microglial_cell,cols = c("grey","#B40F20"))
###**marker散点图------
####Marker expression--------------------
marker <- c("STMN2","NRGN","SLC17A7",#Excitatory_neuron
            "GAD1","RELN","CCK","GAD2",#Inhibitory_neuron
            "MOG","CNP","MAG","CNDP1",#Oligodendrocyte 
            "GJA1","AQP4","SLC14A1","RFX4",#Astrocyte
            "CD74","APBB1IP","TBXAS1","DOCK8",#Microglial_cell
            "ABCB1","CLDN5","EBF1","APOLD1",#Endothelial
            "BCAN","CA10","COL9A1","LIMA1"#Oligodendrocyte_progenitor_cell
)
marker <- c("NRGN","SLC17A7",#Excitatory_neuron
            "GAD1","GAD2",#Inhibitory_neuron
            "MOG","CNDP1",#Oligodendrocyte 
            "GJA1","AQP4",#Astrocyte
            "APBB1IP","DOCK8",#Microglial_cell
            "ABCB1","CLDN5",#Endothelial
            "BCAN","CA10"#Oligodendrocyte_progenitor_cell
)
DotPlot(scRNA, features = marker)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("grey","#B40F20"))
ggsave("marker_final.pdf")

####**celltype  porpotion------------------
###细胞比例----------
######相关性柱状图   纵轴log (percentage ep/percentage nor)----
scRNA.hippo <- readRDS(file = "./scRNA_anno_cluster.rds")

names(table(scRNA$Anno_Idents))
scRNA$Anno_Idents <- factor(scRNA$Anno_Idents,levels = c("Excitatory neuron","Inhibitory neuron","Astrocyte","Microglial",
                                                         "Oligodendrocyte" ,"Oligodendrocyte progenitor cell",
                                                         "Endothelial"))
All.an <- prop.table(table(scRNA$Anno_Idents[scRNA@meta.data$group %in% c("A")]))
All.po <- prop.table(table(scRNA$Anno_Idents[scRNA@meta.data$group %in% c("P")]))
df_prop <- cbind(All.an,All.po)
df_prop <- t(apply(df_prop, 1, function(x) log(x / x[1])))
df_prop <- as.data.frame(df_prop)
df_prop$celltype <- rownames(df_prop)
df_prop <- df_prop[,-1]
class(df_prop)
colnames(df_prop) <- c("log","celltype")
df_prop$celltype <- factor(df_prop$celltype,levels = c(names(table(scRNA$Anno_Idents))))

p1 <- ggplot(df_prop, aes(x=celltype,y=log))+ 
  geom_bar(position="dodge",stat="identity",width=0.8)+
  theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size =8),
        panel.background = element_rect(color = 'black', fill = 'transparent'))
ggsave("./logPA.pdf",p1,height = 5,width = 5)


######细胞类型在前后端比例----
Excitatory <- prop.table(table(scRNA$group[scRNA@meta.data$Anno_Idents %in% c("Excitatory neuron")]))
Inhibitory <- prop.table(table(scRNA$group[scRNA@meta.data$Anno_Idents %in% c("Inhibitory neuron")]))
Oligodendrocyte <- prop.table(table(scRNA$group[scRNA@meta.data$Anno_Idents %in% c("Oligodendrocyte")]))
opc <- prop.table(table(scRNA$group[scRNA@meta.data$Anno_Idents %in% c("Oligodendrocyte progenitor cell")]))
Microglial <- prop.table(table(scRNA$group[scRNA@meta.data$Anno_Idents %in% c("Microglial")]))
Astrocyte <- prop.table(table(scRNA$group[scRNA@meta.data$Anno_Idents %in% c("Astrocyte")]))
Endothelial <- prop.table(table(scRNA$group[scRNA@meta.data$Anno_Idents %in% c("Endothelial")]))


cluster <- c(rep(names(table(scRNA$Anno_Idents))[1],2),
             rep(names(table(scRNA$Anno_Idents))[2],2),rep(names(table(scRNA$Anno_Idents))[3],2),
             rep(names(table(scRNA$Anno_Idents))[4],2),rep(names(table(scRNA$Anno_Idents))[5],2),
             rep(names(table(scRNA$Anno_Idents))[6],2),rep(names(table(scRNA$Anno_Idents))[7],2))
#Samples_1 <- factor(Samples_1,levels = names(table(scRNA@meta.data$idents)))
pos <- c(rep(rep(names(table(scRNA$group))),7))
library("ggplot2")
cell.prop<-as.data.frame(c(Excitatory,Inhibitory,Astrocyte,Microglial,Oligodendrocyte,opc,Endothelial),pos)
cell.prop$pos <- rownames(cell.prop) 
cell.prop$cluster <- cluster
colnames(cell.prop)<-c("proportion","pos","cluster")
library('ggplot2')
library('reshape2')
cell.prop$cluster <- factor(cell.prop$cluster,levels = c(names(table(scRNA$Anno_Idents))))
p2 <- ggplot(cell.prop,aes(cluster,proportion,fill=pos))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+scale_fill_manual(values = c("#99CCCC","#FFCC99"))+
  theme(axis.text.x = element_text(angle=30))

###统计细胞类型在前后端数目差异---------
Excitatory <- prop.table(table(scRNA$orig.ident[scRNA@meta.data$Anno_Idents %in% c("Excitatory neuron")]))
Inhibitory <- prop.table(table(scRNA$orig.ident[scRNA@meta.data$Anno_Idents %in% c("Inhibitory neuron")]))
Oligodendrocyte <- prop.table(table(scRNA$orig.ident[scRNA@meta.data$Anno_Idents %in% c("Oligodendrocyte")]))
opc <- prop.table(table(scRNA$orig.ident[scRNA@meta.data$Anno_Idents %in% c("Oligodendrocyte progenitor cell")]))
Microglial <- prop.table(table(scRNA$orig.ident[scRNA@meta.data$Anno_Idents %in% c("Microglial")]))
Astrocyte <- prop.table(table(scRNA$orig.ident[scRNA@meta.data$Anno_Idents %in% c("Astrocyte")]))
Endothelial <- prop.table(table(scRNA$orig.ident[scRNA@meta.data$Anno_Idents %in% c("Endothelial")]))

cluster <- c(rep(names(table(scRNA$Anno_Idents))[1],10),
             rep(names(table(scRNA$Anno_Idents))[2],10),rep(names(table(scRNA$Anno_Idents))[3],10),
             rep(names(table(scRNA$Anno_Idents))[4],10),rep(names(table(scRNA$Anno_Idents))[5],10),
             rep(names(table(scRNA$Anno_Idents))[6],10),rep(names(table(scRNA$Anno_Idents))[7],10))
pos <- c(rep(names(table(scRNA$group))[1],5),rep(names(table(scRNA$group))[2],5),
         rep(names(table(scRNA$group))[1],5),rep(names(table(scRNA$group))[2],5),
         rep(names(table(scRNA$group))[1],5),rep(names(table(scRNA$group))[2],5),
         rep(names(table(scRNA$group))[1],5),rep(names(table(scRNA$group))[2],5),
         rep(names(table(scRNA$group))[1],5),rep(names(table(scRNA$group))[2],5),
         rep(names(table(scRNA$group))[1],5),rep(names(table(scRNA$group))[2],5),
         rep(names(table(scRNA$group))[1],5),rep(names(table(scRNA$group))[2],5))
library("ggplot2")
cell.prop<-as.data.frame(c(Excitatory,Inhibitory,Astrocyte,Microglial,Oligodendrocyte,opc,Endothelial),pos)
cell.prop$pos <- rownames(cell.prop) 
cell.prop$cluster <- cluster
colnames(cell.prop)<-c("proportion","pos","cluster")
library('ggplot2')
library('reshape2')
cell.prop$cluster <- factor(cell.prop$cluster,levels = c(names(table(scRNA$Anno_Idents))))
ggplot(cell.prop, aes(fill=pos, y=proportion, x=cluster))+ 
  geom_bar(width=0.8,position=position_dodge(width=0.8),stat="identity",size=2)+
  theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size =8),
        panel.background = element_rect(color = 'black', fill = 'transparent'))+
  scale_fill_manual(values = c("#99CCCC","#FFCC99"))

library(ggpubr)
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

###统计细胞类型在样本中差异---------
samplelist <- names(table(scRNA@meta.data$orig.ident))
cell.prop <- data.frame()
for (i in 1:length(samplelist)) {
  a <- prop.table(table(scRNA$Anno_Idents[scRNA@meta.data$orig.ident %in% samplelist[i]]))
  a <- as.data.frame(a)
  a$sample <- rep(samplelist[i],7)
  cell.prop <- rbind(cell.prop,a)
}
All.A <- prop.table(table(scRNA$Anno_Idents[scRNA@meta.data$group %in% c("A")]))
All.P <- prop.table(table(scRNA$Anno_Idents[scRNA@meta.data$group %in% c("P")]))
All.A <- as.data.frame(All.A)
All.A$sample <- rep("Anterior",7)

All.P <- as.data.frame(All.P)
All.P$sample <- rep("Posterior",7)

cell.prop <- rbind(cell.prop,All.A)
cell.prop <- rbind(cell.prop,All.P)

colnames(cell.prop)<-c("celltype","proportion","sample")

library('ggplot2')
library('reshape2')
cell.prop$cluster <- factor(cell.prop$celltype,levels = c(names(table(scRNA$Anno_Idents))))
cell.prop$Samples <- factor(cell.prop$sample,levels = c(samplelist,"Anterior","Posterior"))

ggplot(cell.prop,aes(Samples,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+scale_fill_manual(values = c("#CCCCFF","#FFCC00","#1590BF","#99CC99","#CC6699","#FF9966","#CC9999"))+
  theme(axis.text.x = element_text(angle=30))
