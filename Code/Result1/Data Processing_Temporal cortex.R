library(patchwork)
library(data.table)
library(ggplot2)
library(stringr)
library(SeuratObject, lib.loc = "/data3/software/R/4.4.1/lib64/R/library")
library(Seurat, lib.loc = "/data3/software/R/4.4.1/lib64/R/library")
library(Matrix)
library(dplyr)
library(tidyr)
library(stats)
library(readr)
GSE160189_celltype_color <- c( "Oligodendrocyte" = "#CCCCFF",
                               "Oligodendrocyte progenitor cell" = "#FFCC00",
                               "Microglial" = "#1590BF","Astrocyte" = "#99CC99",
                               "Excitatory neuron" = "#CC6699","Inhibitory neuron" = "#FF9966",
                               "Endothelial" = "#CC9999")
Epilepsy19_celltype_color <- c( "Oligodendrocyte" = "#669999",
                                "Oligodendrocyte precursor cell" = "#B40F20",
                                "Microglial" = "#E58601","Astrocyte" = "#6699CC",
                                "Excitatory neuron" = "#CC99CC","Inhibitory neuron" = "#99CCCC",
                                "Endothelial" = "#71C08F")
###*读取数据-----
###**Epilesy-19数据集（全颞叶）----
##****过滤后的count计数-----
count_matrices <- readRDS("./count_matrices.rds")
##All.sample.list <- names(count_matrices)[c(1,14:17,2:4,12,13,18:20,5:10,11)]
###****批量创建seurat对象,计算线粒体和红细胞基因比例------
###****全部样本---------
All.sample.list <- names(count_matrices)[c(1,14:17,2:4,12,13,18:20,5:10,11)]
list.order <- c(paste("Nor", 1:10, sep = ""),paste("Ep", 1:9, sep = ""),"NeuN")
All.scRNAlist <- list()
for (i in 1:length(All.sample.list)) {
  All.scRNAlist[[i]] <- CreateSeuratObject(count_matrices[[All.sample.list[i]]], min.cells = 3, min.features = 200,project = All.sample.list[i])
  print(head(x=All.scRNAlist[[i]][[]]))
  All.scRNAlist[[i]]$orig.ident <- NULL
  ##meta.data添加信息
  orig.ident <- data.frame(orig.ident=rep(All.sample.list[[i]],ncol(All.scRNAlist[[i]])))
  rownames(orig.ident) <- row.names(All.scRNAlist[[i]]@meta.data)
  All.scRNAlist[[i]] <- AddMetaData(All.scRNAlist[[i]], orig.ident)
  print("a")
  #查看active.ident是否是"one“
  if(levels(All.scRNAlist[[i]])=="one"){
    print(levels(All.scRNAlist[[i]]))
    ##重新设置active.ident
    names(All.scRNAlist)
    All.scRNAlist[[i]] <- RenameIdents(All.scRNAlist[[i]],"one" = list.order[i])
    print("d")
    ##将active.ident设置为matadata中'idents'
    All.scRNAlist[[i]] <- StashIdent(All.scRNAlist[[i]], save.name = 'idents')
    print("e")
    newIdent <- colnames(x = All.scRNAlist[[i]][["RNA"]])
    All.scRNAlist[[i]] <- RenameCells(All.scRNAlist[[i]],new.names = gsub("one",list.order[i],newIdent))
    print("b")
  }else{
    Idents(All.scRNAlist[[i]])
    
    levels(All.scRNAlist[[i]])#0-20
    
    head(All.scRNAlist[[i]]@meta.data)
    
    new.cluster.ids <- list.order[i]
    
    names(new.cluster.ids) <- levels(All.scRNAlist[[i]])
    
    All.scRNAlist[[i]] <- RenameIdents(All.scRNAlist[[i]], new.cluster.ids)
    
    
    Idents(All.scRNAlist[[i]])
    
    levels(All.scRNAlist[[i]])
    print("z")
    All.scRNAlist[[i]] <- StashIdent(All.scRNAlist[[i]], save.name = 'idents')
    All.scRNAlist[[i]] <- RenameCells(All.scRNAlist[[i]],add.cell.id = list.order[i])
    print("c")
  }
  
  if(T){
    print("oo")
    All.scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(All.scRNAlist[[i]], pattern = "^MT-")
  }
  
  if(T){
    print("KK")
    All.scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(All.scRNAlist[[i]], pattern = "^RP[SL]")
  }
  
  if(T){
    print("nn")
    HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    HB.genes <- CaseMatch(HB.genes,rownames(All.scRNAlist[[i]]))
    All.scRNAlist[[i]][["percent.HB"]] <- PercentageFeatureSet(All.scRNAlist[[i]], features = HB.genes)
  }
}
names(All.scRNAlist) <- All.sample.list

####
####
####
All.scRNA <- merge(All.scRNAlist[[1]], y = All.scRNAlist[2:length(All.scRNAlist)])
All.scRNA
dim(All.scRNA)  
table(All.scRNA@meta.data$idents)  

###
###****质控标准------------
###质控小提琴图
##设置可能用到的主题
theme.set2 <- theme(axis.title.x = element_blank())
##设置绘图元素
plot.features <- c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB","percent.rb")
#质控前小提琴图
plots <- list()
for (i in 1:length(plot.features)) {
  plots[[i]] <- VlnPlot(All.scRNA, group.by = "idents", pt.size = 0.1,
                        features = plot.features[i])+theme.set2+NoLegend()
}
violin <- wrap_plots(plots = plots,nrow = 2)
ggsave("./vlnplot_before_qc.pdf", plot = violin, width = 12, height = 12) 
plot2 <- FeatureScatter(All.scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("./pearplot_before_qc.pdf", plot = plot2, width = 5, height = 5) 

All.scRNA <- subset(All.scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
dim(All.scRNA)  #28306 112671
col.num <- length(levels(as.factor(All.scRNA@meta.data$idents)))
violin <-VlnPlot(All.scRNA, group.by = "idents",
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("./vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 
saveRDS(All.scRNA,"./All.scRNA_afterQC.rds")

###****批次效应--------
All.scRNA <- readRDS(".\\All.scRNA_afterQC.rds")
All.scRNA <- NormalizeData(All.scRNA) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData()
scRNA <- RunPCA(All.scRNA, verbose = F)
ElbowPlot(scRNA,ndims = 50)
DimPlot(scRNA, reduction = "pca",cols = mycolors,group.by = "idents")
scRNA <- scRNA %>% RunUMAP(dims = 1:30, verbose = T) %>% RunTSNE(dims = 1:30)
scRNA <- FindNeighbors(scRNA, dims = 1:30)
pdf("./batch_before_cluster.pdf")
for (i in c(0.2, 0.3, 0.4, 0.5,0.6,0.7, 0.8, 1, 1.2, 1.5, 2)) {
  scRNA <- FindClusters(scRNA, resolution = i)
  p1 <- DimPlot(scRNA, reduction = "umap",group.by = "idents",split.by = "idents",ncol = 2) + labs(title = paste0("resolution: ", i))
  p2 <- DimPlot(scRNA, reduction = "tsne",group.by="idents", pt.size=0.3,label=T) + labs(title = paste0("resolution: ", i))
  p3 <- DimPlot(scRNA, reduction = "umap",group.by = "idents",pt.size=0.3,label=T)
  p4 <- DimPlot(scRNA, reduction = "tsne",label=T)
  p5 <- DimPlot(scRNA, reduction = "umap",label=T)
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
}
dev.off()
library(clustree)
library(patchwork)
p1 <- clustree(scRNA, prefix = 'RNA_snn_res.') + coord_flip()
p2 <- DimPlot(scRNA, group.by = 'RNA_snn_res.0.6', label = T)
dim(scRNA) #28306 112671
colnames(scRNA@meta.data)
unique(scRNA@meta.data$idents)
table(scRNA@meta.data$idents)

scRNA@meta.data$Group[scRNA@meta.data$idents %in% c("Nor1",  "Nor2",  "Nor3",  "Nor4",  "Nor5",  
                                                    "Nor6",  "Nor7",  "Nor8",  "Nor9",  "Nor10")] <- "Nor"
scRNA@meta.data$Group[scRNA@meta.data$idents %in% c("Ep1",   "Ep2",   "Ep3",   "Ep4",   "Ep5",   
                                                    "Ep6",   "Ep7",   "Ep8",   "Ep9")] <- "Ep"
scRNA@meta.data$Group[scRNA@meta.data$idents %in% c("NeuN")] <- "NeuN"
scRNA <- subset(scRNA, Group %in% c("Nor","Ep"))
dim(scRNA) #28306 110169
colnames(scRNA@meta.data)
p2 <- DimPlot(scRNA, group.by = 'idents', label = F,cols = mycolors,raster = F)
p3 <- DimPlot(scRNA, group.by = 'Group', label = F,cols = c("Nor"="#2878B5","Ep"="#C82423"),raster = F)
ggsave("./Epilepsy19_orig.ident_batch.pdf",p2,width = 10,height = 8) 
ggsave("./Epilepsy19_Group_batch.pdf",p3,width = 10,height = 8) 
###
####
###****经观察，不存在批次效应
#将resolution=0.6结果作为降维聚类结果
saveRDS(scRNA,file = "./scRNA_cluster.rds")

#####****手动注释---------
scRNA <- readRDS(file = "./scRNA_cluster.rds")
Idents(scRNA) <- 'RNA_snn_res.0.6'
levels(scRNA)
scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(scRNA.markers,"./scRNA.markers.txt",sep = "\t",quote = F)
# scRNA.markers <- read.table("./scRNA.markers.txt",sep = "\t")
scRNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(scRNA, features = top5$gene,slot = ) +scale_fill_gradientn(colors = c("#006699","#FFFFCC","#FF6600"))

BrainMarker <- read.table("./CellMarker_brain_embryo_marker.txt",sep = "\t",header = T)
names(table(BrainMarker$cell_name))
scRNA.markers <- read.table("./scRNA.markers.txt",sep = "\t")

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
setwd("./cellmarker_marker")
# pdf("E:\\original\\epilepsy_Data\\Epilepsy19-master\\FIRST/celltype_marker_scRna.pdf")
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
DimPlot(scRNA, reduction = "umap",group.by = "RNA_snn_res.0.6",label=T,cols = mycolors,label.box = T)

####
####
####
####*****细胞注释完成-------------
# 0:Excitatory neuron  # 1:Inhibitory neuron   # 2:In# 3:ex# 4:In  sst# 5:ex# 6:In   sst# 7:ex# 8:ex# 9:ex# 10:ex
# 11:In# 12:In# 13:ex# 14:In# 15:ex# 16:ex# 17:ex# 18:Oli# 19:ex# 20:ex
# 21:ex# 22:ex# 23:ex# 24:ex# 25:In# 26:ex# 27:astro# 28:ex# 29:ex# 30:ex
# 31:In# 32:micro# 33:In# 34:astro# 35:endo/pericyte# 36:ex
Idents(scRNA) <- 'RNA_snn_res.0.6'
levels(scRNA)
new.cluster.ids <- c("Excitatory neuron","Inhibitory neuron","Excitatory neuron","Inhibitory neuron","Excitatory neuron","Inhibitory neuron",
                     "Excitatory neuron","Inhibitory neuron","Excitatory neuron","Excitatory neuron","Inhibitory neuron",
                     "Excitatory neuron","Excitatory neuron","Inhibitory neuron","Excitatory neuron","Excitatory neuron",
                     "Excitatory neuron","Oligodendrocyte", "Excitatory neuron","Excitatory neuron","Excitatory neuron",
                     "Excitatory neuron","Inhibitory neuron","Astrocyte","Excitatory neuron","Inhibitory neuron",
                     "Excitatory neuron","Excitatory neuron","Inhibitory neuron","Microglial","Oligodendrocyte precursor cell",
                     "Endothelial","Excitatory neuron")
unique(new.cluster.ids)
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
scRNA <- StashIdent(scRNA, save.name = 'Anno_Idents')

DimPlot(scRNA, reduction='umap',group.by="Anno_Idents", pt.size=5,label=T,label.size = 5,raster = F,
        cols = c(color_cols))
saveRDS(scRNA,file = ".\\scRNA_anno_cluster.rds")


#*****展示注释所用marker-------
library(here)
scRNA <- readRDS(file = ".\\scRNA_anno_cluster.rds")

celltype.list <- unique(c("Excitatory neuron","Inhibitory neuron","Excitatory neuron","Inhibitory neuron","Excitatory neuron","Inhibitory neuron",
                          "Excitatory neuron","Inhibitory neuron","Excitatory neuron","Excitatory neuron","Inhibitory neuron",
                          "Excitatory neuron","Excitatory neuron","Inhibitory neuron","Excitatory neuron","Excitatory neuron",
                          "Excitatory neuron","Oligodendrocyte", "Excitatory neuron","Excitatory neuron","Excitatory neuron",
                          "Excitatory neuron","Inhibitory neuron","Astrocyte","Excitatory neuron","Inhibitory neuron",
                          "Excitatory neuron","Excitatory neuron","Inhibitory neuron","Microglial","Oligodendrocyte precursor cell",
                          "Endothelial","Excitatory neuron"))
marker <- intersect(scRNA.markers$gene,BrainMarker[BrainMarker$cell_name %in% celltype.list[[7]],2])
marker
StackedVlnPlot(scRNA, Excitatory_neuron, pt.size=0, cols=my36colors)
Astrocyte <- c("GJA1","AQP4")
Endothelial<- c("VWF","FLT1")
Excitatory_neuron<- c("SATB2","SLC17A7")
Oligodendrocyte_precursor_cell<- c("VCAN","PDGFRA")
Inhibitory_neuron<- c("GAD2","DLX6-AS1")
Microglial_cell<- c("APBB1IP","C1QB")
Oligodendrocyte<- c("CLDN11","MOBP")
FeaturePlot(scRNA,features = Astrocyte,cols = c("grey","#B40F20"))
###*****marker散点图------
####Marker expression--------------------
library(here)
marker <- c("SATB2","SLC17A7",#Excitatory_neuron
            "GAD2","DLX6-AS1",#Inhibitory_neuron
            "VCAN","PDGFRA",#Oligodendrocyte_precursor_cell
            "CLDN11","MOBP",#Oligodendrocyte 
            "GJA1","AQP4",#Astrocyte
            "APBB1IP","C1QB",#Microglial_cell
            "VWF","FLT1"#Endothelial
)
DotPlot(scRNA, features = marker)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
###*****marker散点图+小提琴图------
sce <- as.SingleCellExperiment(scRNA)
plotExpressionCustom <- function(sce, features, features_name, anno_name = "ident",
                                 point_alpha=0.5, point_size=4, ncol=2, xlab = NULL,
                                 exprs_values = "logcounts", scales = "fixed"){
  scater::plotExpression(sce, 
                         exprs_values = exprs_values, 
                         features = features,
                         x = anno_name, 
                         colour_by = anno_name,
                         ncol = ncol,
                         xlab = xlab,
                         point_alpha = point_alpha, 
                         point_size = point_size,
                         add_legend = F,
                         scales = scales) +
    stat_summary(fun = median, 
                 fun.min = median, 
                 fun.max = median,
                 geom = "crossbar", 
                 width = 0.3) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(face = "italic")) +  
    ggtitle(label=paste0(features_name, " markers"))
}

markers.mathys.custom = list(
  'Excitatory_neuron' = c("SATB2","SLC17A7"),
  'Inhibitory_neuron' = c("GAD2","DLX6-AS1"),
  'Oligodendrocyte_precursor_cell' = c("VCAN","PDGFRA"),
  'Oligodendrocyte' = c("CLDN11","MOBP"),
  'Astrocyte' = c("GJA1","AQP4"),
  'Microglial_cell' = c("APBB1IP","C1QB"),
  'Endothelial' = c("VWF","FLT1")
)
pdf(".\\marker1.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  
  p <- plotExpressionCustom(sce = sce,
                            features = markers.mathys.custom[[i]], 
                            features_name = names(markers.mathys.custom)[[i]], 
                            anno_name = "ident") +
    scale_color_manual(values = color_cols)
  print(p)
}
dev.off()
pdf("./plots_markers.pdf")
for (i in 0:36) {
  p <- StackedVlnPlot(use_tsne, top10_tsne$gene[which(top10_tsne$cluster==i)], pt.size=0, cols=mycolors) 
  print(p)
  
}
dev.off()

####*****sample celltype  porpotion------------------
###
names(table(scRNA@meta.data$idents))
############比例图:
Ep1 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Ep1")]))
Ep2 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Ep2")]))
Ep3 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Ep3")]))
Ep4 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Ep4")]))
Ep5 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Ep5")]))
Ep6 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Ep6")]))
Ep7 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Ep7")]))
Ep8 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Ep8")]))
Ep9 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Ep9")]))
Nor1 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Nor1")]))
Nor2 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Nor2")]))
Nor3 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Nor3")]))
Nor4 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Nor4")]))
Nor5 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Nor5")]))
Nor6 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Nor6")]))
Nor7 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Nor7")]))
Nor8 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Nor8")]))
Nor9 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Nor9")]))
Nor10 <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Nor10")]))
NeuN <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("NeuN")]))

Samples <- c(rep("Ep1",7),rep("Ep2",7),rep("Ep3",7),rep("Ep4",7),rep("Ep5",7),rep("Ep6",7),
             rep("Ep7",7),rep("Ep8",7),rep("Ep9",7),rep("Nor1",7),rep("Nor2",7),rep("Nor3",7),
             rep("Nor4",7),rep("Nor5",7),rep("Nor6",7),rep("Nor7",7),rep("Nor8",7),rep("Nor9",7),rep("Nor10",7),rep("NeuN",7))
Samples <- factor(Samples,levels = names(table(scRNA@meta.data$idents)))
cluster <- rep(names(table(scRNA@active.ident)),20)
## 绘制堆叠条形图
library("ggplot2")
cell.prop<-as.data.frame(c(Ep1,Ep2,Ep3,Ep4,Ep5,Ep6,Ep7,Ep8,Ep9,Nor1,
                           Nor2,Nor3,Nor4,Nor5,Nor6,Nor7,Nor8,Nor9,Nor10,NeuN),cluster)
cell.prop$cluster <- rownames(cell.prop) 
cell.prop$Samples <- Samples
colnames(cell.prop)<-c("proportion","Cluster","Samples")
cell.prop <- cell.prop[-which(cell.prop$Samples=="NeuN"),]
library('ggplot2')
library('reshape2')
color_cols_1 <- c(rgb(102,153,204,250,maxColorValue = 255),#6699CC
                  rgb(108,211,152,250,maxColorValue = 255),#336699
                  rgb(204,153,204,250,maxColorValue = 255),#CC99CC
                  rgb(153,204,204,250,maxColorValue = 255), #99CCCC
                  rgb(229,134,1,250,maxColorValue = 255),#E58601
                  rgb(102,153,153,250,maxColorValue = 255),#669999 
                  rgb(180,15,32,250,maxColorValue = 255))#B40F20)

ggplot(cell.prop,aes(Samples,proportion,fill=Cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+scale_fill_manual(values = color_cols_1)+
  theme(axis.text.x = element_text(angle=30))
# "#6699CC","#336699", "#B40F20",'#99CCCC', "#E58601","#669999","#CC99CC"
# Astrocyte  endo   ex  in  Microglial_cell   Oligodendrocyte   Oligodendrocyte_progenitor_cell


###****计算疾病-正常细胞类型比例图------------
All.ep <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Ep1","Ep2","Ep3","Ep4","Ep5","Ep6","Ep7","Ep8",'Ep9')]))
All.Nor <- prop.table(table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Nor1","Nor2","Nor3","Nor4","Nor5","Nor6","Nor7","Nor8",'Nor9',"Nor10")]))
Samples_1 <- c(rep("All.ep",7),rep("All.Nor",7),rep("NeuN",7))
#Samples_1 <- factor(Samples_1,levels = names(table(scRNA@meta.data$idents)))
cluster <- c(rep(rep(names(table(scRNA@active.ident))),3))

library("ggplot2")
cell.prop<-as.data.frame(c(All.ep,All.Nor,NeuN),cluster)
cell.prop$cluster <- rownames(cell.prop) 
cell.prop$Samples <- Samples_1
colnames(cell.prop)<-c("proportion","Cluster","Samples")
cell.prop <- cell.prop[-which(cell.prop$Samples=="NeuN"),]
library('ggplot2')
library('reshape2')

ggplot(cell.prop,aes(Samples,proportion,fill=Cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+scale_fill_manual(values = color_cols_1)+
  theme(axis.text.x = element_text(angle=30))
##计算卡fang检验，二维列联表
table(Idents(scRNA)[scRNA@meta.data$idents %in% c("Ep1","Ep2","Ep3","Ep4","Ep5","Ep6","Ep7","Ep8",'Ep9')])
table(scRNA@active.ident[scRNA@meta.data$idents%in% c("Nor1","Nor2","Nor3","Nor4","Nor5","Nor6","Nor7","Nor8",'Nor9',"Nor10")])

Nor<-c(40026,21310,377,735,9,17,8)
Ep<-c(27797,18302,662,413,223,216,74)
data<-data.frame(Nor,Ep)
data<-t(data)
chisq.test(data)   # p-value < 2.2e-16






