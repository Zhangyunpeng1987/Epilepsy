####差异分析以及功能富集--------------
####*每个细胞类型在海马体前后差异分析--------
scRNA <- readRDS(file = "./scRNA_anno_cluster.rds")
metadata <- scRNA@meta.data
##meta.data添加信息
group.id <- as.data.frame(metadata[,1])
rownames(group.id) <- rownames(metadata)
group.id$group <- t(as.data.frame(stringr::str_extract_all(group.id$`metadata[, 1]`, '\\D+')))
group <- factor(group.id[,2])
class(group)
scRNA <- AddMetaData(scRNA, group,col.name = "group")
dim(scRNA)  #16956 131096
####*Astrocyte----------------
Astrocyte <- subset(scRNA, Anno_Idents=="Astrocyte")
diff_Astrocyte <- FindMarkers(Astrocyte, min.pct = 0.1, 
                              logfc.threshold = 0.25,
                              test.use = "wilcox",
                              group.by = "group",
                              ident.1 ="P",
                              ident.2="A")
diff_Astrocyte<-cbind(rownames(diff_Astrocyte),diff_Astrocyte)  ##对数据增加一列
colnames(diff_Astrocyte)<- c('gene_id',"pvalue" ,"avglog2FoldChange","pct.1" ,"pct.2","padj" )
##对列名进行修改  1600
write.table(diff_Astrocyte,"./wilcox\\case-control-astro-wilcox.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
diff_Astrocyte <- read.table("./case-control-astro-wilcox.txt",sep = '\t',header = T)
diff_Astrocyte <- diff_Astrocyte[which(diff_Astrocyte$pvalue<0.01),]
astro.gene <- diff_Astrocyte[,1]


#对于加载的注释库的使用，以上述为例，就直接在 OrgDb 中指定人（org.Hs.eg.db）或绵羊（sheep）
astro.enrich.go <- enrichGO(gene = astro.gene,  #基因列表文件中的基因名称
                            OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                            keyType = 'SYMBOL',  
                            ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                            pAdjustMethod = 'fdr',  #指定 p 值校正方法
                            pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                            qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                            readable = FALSE)
astro.enrich.go <-summary(astro.enrich.go)  #179
write.table(astro.enrich.go,"./enrich\\astro.enrich.go.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
#其他细胞类型上同

####******富集热图+点图----------
#BiocManager::install("BoutrosLab.plotting.general")
library(BoutrosLab.plotting.general)
library(ggplot2)
##根据GO term 的富集程度即P值画热图，行为GO term，列为细胞类型，值为logP
##ENDo
Endo.enrich.go <- read.table("./enrich\\Endo.enrich.go.txt",sep = '\t',header = T)
Endo.BP <- Endo.enrich.go[order(Endo.enrich.go$pvalue),]
Endo.BP <- Endo.BP[1:10,c(2,5,9)]
rownames(Endo.BP) <- Endo.BP$Description
colnames(Endo.BP) <- c("Description","Endo_p","Endo_count")

####
Astro.BP <- Astro.enrich.go[order(Astro.enrich.go$pvalue),]
Astro.BP <- Astro.BP[1:10,c(2,5,9)]
rownames(Astro.BP) <- Astro.BP$Description
colnames(Astro.BP) <- c("Description","Astro_p","Astro_count")

####
Microglial.BP <- Microglial.enrich.go[order(Microglial.enrich.go$pvalue),]
Microglial.BP <- Microglial.BP[1:10,c(2,5,9)]
rownames(Microglial.BP) <- Microglial.BP$Description
colnames(Microglial.BP) <- c("Description","Micro_p","Micro_count")

####
OPC.BP <- OPC.enrich.go[order(OPC.enrich.go$pvalue),]
OPC.BP <- OPC.BP[1:10,c(2,5,9)]
rownames(OPC.BP) <- OPC.BP$Description
colnames(OPC.BP) <- c("Description","OPC_p","OPC_count")

####
Oligodendrocyte.BP <- Oligodendrocyte.enrich.go[order(Oligodendrocyte.enrich.go$pvalue),]
Oligodendrocyte.BP <- Oligodendrocyte.BP[1:10,c(2,5,9)]
rownames(Oligodendrocyte.BP) <- Oligodendrocyte.BP$Description
colnames(Oligodendrocyte.BP) <- c("Description","Oli_p","Oli_count")

####
Excitatory.BP <- Excitatory_neuron.enrich.go[order(Excitatory_neuron.enrich.go$pvalue),]
Excitatory.BP <- Excitatory.BP[1:10,c(2,5,9)]
rownames(Excitatory.BP) <- Excitatory.BP$Description
colnames(Excitatory.BP) <- c("Description","Ex_p","Ex_count")

####
Inhbitory.BP <- Inhbitory_neuron.enrich.go[order(Inhbitory_neuron.enrich.go$pvalue),]
Inhbitory.BP <- Inhbitory.BP[1:10,c(2,5,9)]
rownames(Inhbitory.BP) <- Inhbitory.BP$Description
colnames(Inhbitory.BP) <- c("Description","In_p","In_count")

####构建矩阵--------
EX_IN <- merge(Excitatory.BP, Inhbitory.BP, by='Description', all=TRUE)
EX_IN_Oli <- merge(EX_IN, Oligodendrocyte.BP, by='Description', all=TRUE)
EX_IN_Oli_Astro <- merge(EX_IN_Oli, Astro.BP, by='Description', all=TRUE)
EX_IN_Oli_Astro_Opcs <- merge(EX_IN_Oli_Astro, OPC.BP, by='Description', all=TRUE)
EX_IN_Oli_Astro_Opcs_micro <- merge(EX_IN_Oli_Astro_Opcs, Microglial.BP, by='Description', all=TRUE)

EX_IN_Oli_Astro_Opcs_micro_Endo <- merge(EX_IN_Oli_Astro_Opcs_micro, Endo.BP, by='Description', all=TRUE)
dim(EX_IN_Oli_Astro_Opcs_micro_Endo)
#33 15
exp <- EX_IN_Oli_Astro_Opcs_micro_Endo[,c(2:ncol(EX_IN_Oli_Astro_Opcs_micro_Endo))]
rownames(exp) <- EX_IN_Oli_Astro_Opcs_micro_Endo$Description

spot.size.function<-function(x){
  x= (-log2(x))/15
}
#c("#CCCCFF",         "#FFCC00",                        "#0099CC",        "#99CC99",  "#CC6699",  "#FF9966", "#CC9999")
# Oligodendrocyte    Oligodendrocyte_progenitor_cell    Microglial_cell  Astrocyte     ex           in        endo     

spot.colour.function <- function(x) {
  colours <- rep("white", length(x));
  colours[x==exp$Ex_p] <- "#CC6699"; 
  colours[x==exp$In_p] <- "#FF9966"; 
  colours[x==exp$Oli_p] <- "#CCCCFF";
  colours[x==exp$Astro_p] <- "#99CC99";
  colours[x==exp$OPC_p] <- "#FFCC00";
  colours[x==exp$Micro_p] <- "#0099CC";
  colours[x==exp$Endo_p] <- "#CC9999";
  return(colours);
}
color_cols <- c("#CC6699", "#FF9966", "#CCCCFF","#99CC99","#FFCC00","#0099CC","#CC9999")
#加图例
key.trans <- list(title = "cell type",space = "right",columns = 1,
                  points = list(pch = c(20,20,20,20,20,20,20),
                                col = color_cols,
                                cex=c(1,1,1,1,1,1,1)),
                  text = list(c("Excitatory neuron","Inhibitory neuron","Oligodendrocyte","Astrocyte",
                                "OPCs","Microglial","Endothelial")),
                  cex.title = 1,cex = .9)
create.dotmap(
  exp[,c(1,3,5,7,9,11,13)],bg.data = exp[,c(2,4,6,8,10,12,14)],
  pch = 20,na.spot.size=2,spot.size.function=spot.size.function,
  spot.colour.function=spot.colour.function,
  ylab.cex=1,yaxis.cex=0.5,xaxis.cex=0.5,
  colour.scheme=colorRampPalette(c("#80B1D3","#FDB462","#E59CC4", "#BC80BD"))(100),
  key = key.trans,total.colours = 10,colourkey = T
)


####**细胞类型间差异分析--------
devtools::install_github("junjunlab/jjAnno")
devtools::install_github("junjunlab/ClusterGVis")
library(ClusterGVis)
library(org.Hs.eg.db)
scRNA <- readRDS(file = "./scRNA_anno_cluster.rds")
Idents(scRNA) <- "Anno_Idents"
diff_celltype <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(diff_celltype,".\\wilcox\\diff_celltype.txt",sep = "\t")
####******富集分析----------
diff_celltype<- read.table(".\\wilcox\\diff_celltype.txt",sep = "\t")

table(diff_celltype$cluster)
# Oligodendrocyte Oligodendrocyte progenitor cell                      Microglial                       Astrocyte 
# 505                             374                             346                             537 
# Excitatory neuron               Inhibitory neuron                     Endothelial 
# 956                             793                             187
Excitatory.gene <- diff_celltype[diff_celltype$cluster=="Excitatory neuron",7]
Excitatory <- enrichGO(gene = Excitatory.gene,  #基因列表文件中的基因名称
                       OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                       keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                       ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                       pAdjustMethod = 'fdr',  #指定 p 值校正方法
                       pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                       qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                       readable = FALSE)
Excitatory <-as.data.frame(Excitatory)  #439
write.table(Excitatory,"./enrich\\celltype_Excitatory.go.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
#same other cell types

######******细胞类型间差异基因-----------------
Endo.enrich.go <- read.table("./enrich\\celltype_Endothelial.go.txt",sep = '\t',header = T)
Endo.BP <- Endo.enrich.go[order(Endo.enrich.go$pvalue),]

####
Astro.BP <- Astro.enrich.go[order(Astro.enrich.go$pvalue),]

####
Microglial.BP <- Microglial.enrich.go[order(Microglial.enrich.go$pvalue),]

####
OPC.BP <- OPC.enrich.go[order(OPC.enrich.go$pvalue),]

####
Oligodendrocyte.BP <- Oligodendrocyte.enrich.go[order(Oligodendrocyte.enrich.go$pvalue),]

####
Excitatory.BP <- Excitatory_neuron.enrich.go[order(Excitatory_neuron.enrich.go$pvalue),]

####
Inhbitory.BP <- Inhbitory_neuron.enrich.go[order(Inhbitory_neuron.enrich.go$pvalue),]

hippo_go <- rbind(Inhbitory.BP,Excitatory.BP,Oligodendrocyte.BP,OPC.BP,Microglial.BP,Astro.BP,Endo.BP)

####******构建矩阵--------
Astro.BP <- Astro.BP[1:10,c(2,5,9)]
rownames(Astro.BP) <- Astro.BP$Description
colnames(Astro.BP) <- c("Description","Astro_p","Astro_count")

Oligodendrocyte.BP <- Oligodendrocyte.BP[1:10,c(2,5,9)]
rownames(Oligodendrocyte.BP) <- Oligodendrocyte.BP$Description
colnames(Oligodendrocyte.BP) <- c("Description","Oli_p","Oli_count")


OPC.BP <- OPC.BP[1:10,c(2,5,9)]
rownames(OPC.BP) <- OPC.BP$Description
colnames(OPC.BP) <- c("Description","OPC_p","OPC_count")


Endo.BP <- Endo.BP[1:10,c(2,5,9)]
rownames(Endo.BP) <- Endo.BP$Description
colnames(Endo.BP) <- c("Description","Endo_p","Endo_count")


Microglial.BP <- Microglial.BP[1:10,c(2,5,9)]
rownames(Microglial.BP) <- Microglial.BP$Description
colnames(Microglial.BP) <- c("Description","Microglial_p","Microglial_count")


Excitatory.BP <- Excitatory.BP[1:10,c(2,5,9)]
rownames(Excitatory.BP) <- Excitatory.BP$Description
colnames(Excitatory.BP) <- c("Description","Excitatory_p","Excitatory_count")


Inhbitory.BP <- Inhbitory.BP[1:10,c(2,5,9)]
rownames(Inhbitory.BP) <- Inhbitory.BP$Description
colnames(Inhbitory.BP) <- c("Description","Inhbitory_p","Inhbitory_count")



EX_IN <- merge(Excitatory.BP, Inhbitory.BP, by='Description', all=TRUE)
EX_IN_Astro <- merge(EX_IN, Astro.BP, by='Description', all=TRUE)
EX_IN_Astro_Micro <- merge(EX_IN_Astro, Microglial.BP, by='Description', all=TRUE)
EX_IN_Astro_Micro_Oli <- merge(EX_IN_Astro_Micro, Oligodendrocyte.BP, by='Description', all=TRUE)
EX_IN_Astro_Micro_Oli_Opcs <- merge(EX_IN_Astro_Micro_Oli, OPC.BP, by='Description', all=TRUE)
EX_IN_Astro_Micro_Oli_Opcs_Endo <- merge(EX_IN_Astro_Micro_Oli_Opcs, Endo.BP, by='Description', all=TRUE)
dim(EX_IN_Astro_Micro_Oli_Opcs_Endo)
#46 15
write.table(EX_IN_Astro_Micro_Oli_Opcs_Endo,"./enrich\\celltype_EX_IN_Astro_Micro_Oli_Opcs_Endo.go.txt",sep = '\t',col.names = T,row.names = F,na='')

exp <- EX_IN_Astro_Micro_Oli_Opcs_Endo[,c(2:ncol(EX_IN_Astro_Micro_Oli_Opcs_Endo))]
rownames(exp) <- EX_IN_Astro_Micro_Oli_Opcs_Endo$Description

spot.size.function<-function(x){
  x= (-log2(x))/50
}
color_cols <- c("#CC6699", "#FF9966", "#99CC99","#0099CC","#CCCCFF","#FFCC00","#CC9999")

spot.colour.function <- function(x) {
  colours <- rep("white", length(x));
  colours[x==exp$Excitatory_p] <- "#CC6699"; 
  colours[x==exp$Inhbitory_p] <- "#FF9966"; 
  colours[x==exp$Astro_p] <- "#99CC99";
  colours[x==exp$Microglial_p] <- "#0099CC";
  colours[x==exp$Oli_p] <- "#CCCCFF";
  colours[x==exp$OPC_p] <- "#FFCC00";
  colours[x==exp$Endo_p] <- "#CC9999";
  return(colours);
}
#加图例
key.trans <- list(title = "cell type",space = "right",columns = 1,
                  points = list(pch = c(20,20,20,20,20,20),
                                col = color_cols,
                                cex=c(1,1,1,1)),
                  text = list(c("Excitatory neuron","Inhibitory neuron","Astrocyte",
                                "Microglial","Oligodendrocyte","OPCs","Endothelial")),
                  
                  #lines = list(col = colors,lty = lines),
                  cex.title = 1,cex = .9)
create.dotmap(
  exp[,c(1,3,5,7,9,11,13)],bg.data = exp[,c(2,4,6,8,10,12,14)],
  pch = 20,na.spot.size=2,spot.size.function=spot.size.function,
  spot.colour.function=spot.colour.function,
  ylab.cex=1,yaxis.cex=0.5,xaxis.cex=0.5,
  colour.scheme=colorRampPalette(c("#80B1D3","#FDB462","#E59CC4", "#BC80BD"))(100),
  key = key.trans,total.colours = 10,colourkey = T,xaxis.rot = 45
)