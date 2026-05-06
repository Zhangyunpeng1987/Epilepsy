library(SCENIC)
packageVersion("SCENIC")
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(zoo)
library(mixtools)
library(rbokeh)
library(DT)
library(NMF)
library(pheatmap)
library(R2HTML)
library(Rtsne)
library(doMC)
library(doRNG)
library(scRNAseq)
library(SCopeLoomR)
library(SCENIC)
library(Seurat)
###数据准备-----------
cols_2<-c("#FFBE7A","#82B0D2","#8ECFC9","#D76364")  #抑制性神经元亚型颜色
Inhibitory <- readRDS(file = "./temporal_lobe\\subcluster\\Inhibitory\\Inhibitory_anno_cluster.rds")
dim(Inhibitory)
DimPlot(Inhibitory, reduction='umap',group.by="Anno_Idents", pt.size=5,label=T,label.size = 5,raster = F,
        cols = cols_2)
table(Inhibitory$Anno_Idents)
table(Inhibitory$group)
Inhibitory_Nor <- subset(Inhibitory, group %in% c("Nor"))  #Ep: 18302   Nor:21310
count <- Inhibitory_Nor@assays$RNA@counts
write.csv(t(as.matrix(count)),file = "./for.scenic.data.csv")

##可视化---------------------
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
###提取细胞类型特异的转录调控因子--------------
###兴奋性神经元-------------
Excitatory <- readRDS(file = "./temporal_lobe\\subcluster\\Excitatory\\Excitatory_anno_cluster_1.rds")
table(Excitatory$Anno_Idents)
 
table(Excitatory$group)

Excitatory_Nor <- subset(Excitatory, group %in% c("Nor"))  
Excitatory_Ep <- subset(Excitatory, group %in% c("Ep"))  

loom <- open_loom('./temporal_lobe\\subcluster\\Excitatory\\pyscenic\\Ex_Nor\\out_SCENIC.loom') 
loom[["matrix"]]
loom[["col_attrs/CellID"]]
loom[["row_attrs/Gene"]]
loom[["col_graphs"]]
full.matrix <- loom[["matrix"]][,]
gene.names <- loom[["row_attrs/Gene"]][]
head(gene.names)
get.attribute.df()
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])
embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)
library(dplyr)
library(tibble)
regulons.df <- regulons %>% unlist() %>% as.data.frame()
regulons.df <- tibble::rownames_to_column(regulons.df)
colnames(regulons.df) <- c("TF", "TargetGene")
TF <-  str_split(regulons.df$TF,'\\(',simplify = T)[,1]
regulons.df$TF <- TF

seurat.data = Excitatory_Ep
sub_regulonAUC <- regulonAUC[,match(colnames(seurat.data),colnames(regulonAUC))]
regulonAUC[1:4,1:4]
dim(sub_regulonAUC) 
seurat.data
identical(colnames(sub_regulonAUC), colnames(seurat.data))

cellClusters <- data.frame(row.names = colnames(seurat.data), 
                           seurat_clusters = as.character(seurat.data$Anno_Idents))
cellTypes <- data.frame(row.names = colnames(seurat.data), 
                        celltype = seurat.data$Anno_Idents)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 


selectedResolution <- "celltype" 
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])

sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)

regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

Heatmap(
  regulonActivity_byGroup_Scaled,
  name                         = "z-score",
  col                          = colorRampPalette(c('#FFFFFF','#6699CC'))(100),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

###rss查看特异TF
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusters[colnames(regulonAUC), ])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
spec.rss <- rssPlot[["plot"]][["data"]]
spec.rss %>%
  group_by(cellType) %>%
  top_n(n = 5, wt = RSS) -> spec.rss
TF <- str_split(spec.rss$Topic,'\\(',simplify = T)[,1]
spec.rss$TF <- TF
spec_target <- merge(spec.rss,regulons.df,by="TF")  
write.table(spec_target[,-2],"./temporal_lobe\\subcluster\\Excitatory\\pyscenic\\Ex_Nor\\spec_target.txt",row.names = F,sep = "\t",quote = F)

###构建网络---------
###靶基因筛选----------------
###Ex_Ep--------------
Ex_Ep <- read.table("./temporal_lobe\\subcluster\\Excitatory\\pyscenic\\Ex_Ep\\spec_target.txt",sep = "\t",header = T)
unique(Ex_Ep$TF)
dim(Ex_Ep) 
length(unique(Ex_Ep$TF))
length(unique(Ex_Ep$TargetGene)) 
Tf<- unique(Ex_Ep$TF)
Ex_Ep_new <- data.frame(matrix(ncol = 5, nrow = 0))
for (i in 1:length(Tf)) {
  TF_module <- enrichGO(gene = c(Tf[i],Ex_Ep$TargetGene[Ex_Ep$TF %in% Tf[i]]),  #基因列表文件中的基因名称
                        OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                        keyType = 'SYMBOL',  
                        ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                        pAdjustMethod = 'BH',  #指定 p 值校正方法
                        pvalueCutoff = 1,  #指定 p 值阈值，不显著的值将不显示在结果中
                        qvalueCutoff = 1,  #指定 q 值阈值，不显著的值将不显示在结果中
                        readable = FALSE)
  TF_module <-as.data.frame(TF_module)  #
  TF_module <- TF_module[TF_module$pvalue < 0.01,]
  setwd("./temporal_lobe\\subcluster\\Excitatory\\pyscenic\\Ex_Ep")
  name=paste(Tf[i],"Ex_Ep.enrich.txt",sep = "_")
  write.table(TF_module,name,sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
  print(dim(TF_module))
  if(dim(TF_module)[1] >20){
    TF_module <- TF_module[order(TF_module$pvalue),][1:20,]
    print(dim(TF_module))
    targetgene <- TF_module$geneID
    print(length(targetgene))
    library(stringr)
    genes <- str_split(targetgene,"/")
    genes <- unique(unlist(genes))
    Ex_Ep_1 <- Ex_Ep[Ex_Ep$TF %in% Tf[i] & Ex_Ep$TargetGene %in% genes,]
    Ex_Ep_new <- rbind(Ex_Ep_new,Ex_Ep_1)
    print(dim(Ex_Ep_new))
  }else{
    targetgene <- TF_module$geneID
    print(length(targetgene))
    library(stringr)
    genes <- str_split(targetgene,"/")
    genes <- unique(unlist(genes))
    Ex_Ep_1 <- Ex_Ep[Ex_Ep$TF %in% Tf[i] & Ex_Ep$TargetGene %in% genes,]
    Ex_Ep_new <- rbind(Ex_Ep_new,Ex_Ep_1)
    print(dim(Ex_Ep_new))
  }
}
write.table(Ex_Ep_new,"./temporal_lobe\\subcluster\\Excitatory\\pyscenic\\Ex_Ep\\Ex_Ep_TF.txt",row.names = F,sep = "\t")
length(intersect(Ex_Ep_new$TF,Ex_Ep$TF))
length(unique(Ex_Ep_new$TargetGene))  

####转录因子筛选-----------
Ex_Ep <- read.table("./temporal_lobe\\subcluster\\Excitatory\\pyscenic\\Ex_Ep\\Ex_Ep_TF.txt",header = T,sep = "\t")
Ex_Nor <- read.table("./temporal_lobe\\subcluster\\Excitatory\\pyscenic\\Ex_Nor\\Ex_Nor_TF.txt",header = T,sep = "\t")
In_Nor <- read.table("./temporal_lobe\\subcluster\\Inhibitory\\pyscenic\\IN_NOR\\In_Nor_TF.txt",header = T,sep = "\t")
In_Ep <- read.table("./temporal_lobe\\subcluster\\Inhibitory\\pyscenic\\In_Ep\\In_Ep_TF.txt",header = T,sep = "\t")
df <- rbind(Ex_Ep,Ex_Nor,In_Nor,In_Ep)
df$Group <- c(rep("Ex_Ep",dim(Ex_Ep)[1]),
              rep("Ex_Nor",dim(Ex_Nor)[1]),
              rep("In_Nor",dim(In_Nor)[1]),
              rep("In_Ep",dim(In_Ep)[1]))
df <- df[,c(1,2,5,6)]
df <- unique(df)
###转录因子筛选：经文献验证---------------
####http://humantfs.ccbr.utoronto.ca/download.php  29425488  The Human Transcription Factors  Cell
Human_TF  <- read.csv("./Hippo\\subcluster\\Human_TF.csv",header = T,sep = ",")
Human_TF <- Human_TF[,-1]
Tf_Cell <- unique(Human_TF$HGNC.symbol[Human_TF$TF.assessment %in% "Known motif" & Human_TF$Is.TF. %in% "Yes" & Human_TF$TF.tested.by.HT.SELEX. %in% c("DBD", "DBD and Full","Full") & Human_TF$TF.tested.by.PBM. %in% "Yes"]) #106
Tf_Cell <- intersect(Tf_Cell,df$TF)
length(Tf_Cell) #40
##http://bioinfo.life.hust.edu.cn/AnimalTFDB4/#/Download
Homo_sapiens_TF  <- read.table("./Hippo\\subcluster\\Homo_sapiens_TF.txt",header = T,sep = "\t")  #1659
Tf_AnimalTFDB4 <- unique(Homo_sapiens_TF$Symbol) #1638
Tf_AnimalTFDB4 <- intersect(Tf_AnimalTFDB4,Tf_Cell)  #98 
length(Tf_AnimalTFDB4) #98

###可参考A comprehensive library of human transcription factors for cell fate  engineering
TF_2020  <- read.table("./Hippo\\subcluster\\A comprehensive library of human transcription factors for cell fate  engineering.txt",header = T,sep = "\t")  #1659
TF_2020_TF <- unique(TF_2020$TF) 
TF_2020_TF <- intersect(TF_2020_TF,Tf_AnimalTFDB4)  
length(TF_2020_TF)  #37

df <- df[df$TF %in% TF_2020_TF,]  

write.table(df,"./temporal_lobe\\subcluster\\all_spec_target.txt",row.names = F,sep = "\t",quote = F)
library(tidyr)
aa <- as.data.frame(table(df$TF,df$Group))
aa <- spread(aa,Var2,Freq)
rownames(aa) <- aa$Var1
aa <- aa[,-1]
rowSums(aa==0)
aa[aa != 0] <- 1
aa$TF <- rownames(aa)
write.table(aa,"./temporal_lobe\\subcluster\\node.txt",row.names = F,sep = "\t",quote = F)

intersect(df$TF,df$TargetGene) #22
node1 <- df[,1]  #1940
length(unique(node1)) #37
node1 <- c(node1,df$TargetGene)
node1 <- as.data.frame(node1)
node1$Type <- c(rep("TF",length(df[,1])),rep("TargetGene",length(df[,3])))
node1$Type[node1$node1 %in% intersect(df$TF,df$TargetGene)] <- "TF"
node1 <- unique(node1)  #1248
table(node1$Type)
# TargetGene         TF 
# 1211         37 
write.table(node1,"./temporal_lobe\\subcluster\\TF_TargetGENE.txt",row.names = F,sep = "\t",quote = F)

#cytoscape细胞类型特异的用细胞类型的颜色，两种细胞类型相同的用不同颜色，疾病正常状态下相同的用不同颜色






