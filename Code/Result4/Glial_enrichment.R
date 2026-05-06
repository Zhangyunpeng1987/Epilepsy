####*胶质细胞富集功能总结-------------
###***小胶质细胞每个亚群在海马体前后端差异基因--------
Microglial <- readRDS(file = "./Hippo\\subcluster\\Microglial\\anno\\Microglial_sub.rds")
celltype <- names(table(Microglial$Anno_Idents))
diff_micro_all <- data.frame(matrix(ncol = 5, nrow = 0))
for (i in 1:length(celltype)) {
  data <- subset(Microglial,Anno_Idents %in% celltype[i])
  diff_celltype <- FindMarkers(data, min.pct = 0.1, 
                               logfc.threshold = 0.25,
                               test.use = "wilcox",
                               group.by = "group",
                               ident.1 ="P",
                               ident.2="A")
  
  diff_celltype$celltype <- rep(celltype[i],dim(diff_celltype)[1])
  diff_micro_all <- rbind.data.frame(diff_micro_all,diff_celltype)
}
write.table(diff_micro_all,"./Hippo\\subcluster\\Microglial\\diff\\diff_micro_all_marker.txt",sep = "\t")

celltype <- unique(diff_micro_all$celltype)
micro.GO <- data.frame(matrix(ncol = 10, nrow = 0))
for (i in 1:length(celltype)) {
  ##在P端上调
  up_exp <- diff_micro_all[diff_micro_all$celltype %in% celltype[i] & diff_micro_all$p_val<0.05 & diff_micro_all$avg_log2FC>0.25,]
  up <- rownames(up_exp)
  up.enrich.go <- enrichGO(gene = up,  #基因列表文件中的基因名称
                           OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                           keyType = 'SYMBOL',  
                           ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                           pAdjustMethod = 'BH',  #指定 p 值校正方法
                           pvalueCutoff = 1,  #指定 p 值阈值，不显著的值将不显示在结果中
                           qvalueCutoff = 1,  #指定 q 值阈值，不显著的值将不显示在结果中
                           readable = FALSE)
  up.enrich.go <-as.data.frame(up.enrich.go)  #
  up.enrich.go <- up.enrich.go[up.enrich.go$pvalue <0.01,]
  up.enrich.go$Group <- rep("Micro_P",dim(up.enrich.go)[1])
  up.enrich.go$Celltype <- rep(celltype[i],dim(up.enrich.go)[1])
  micro.GO <- rbind.data.frame(micro.GO,up.enrich.go)
  ##在A端上调
  down_exp <- diff_micro_all[diff_micro_all$celltype %in% celltype[i] & diff_micro_all$p_val<0.05 & diff_micro_all$avg_log2FC<(-0.25),]
  down <- rownames(down_exp)
  down.enrich.go <- enrichGO(gene = down,  #基因列表文件中的基因名称
                             OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                             keyType = 'SYMBOL',  
                             ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                             pAdjustMethod = 'BH',  #指定 p 值校正方法
                             pvalueCutoff = 1,  #指定 p 值阈值，不显著的值将不显示在结果中
                             qvalueCutoff = 1,  #指定 q 值阈值，不显著的值将不显示在结果中
                             readable = FALSE)
  down.enrich.go <-as.data.frame(down.enrich.go)  #
  down.enrich.go <- down.enrich.go[down.enrich.go$pvalue <0.01,]
  down.enrich.go$Group <- rep("Micro_A",dim(down.enrich.go)[1])
  down.enrich.go$Celltype <- rep(celltype[i],dim(down.enrich.go)[1])
  micro.GO <- rbind(micro.GO,down.enrich.go)
  
}
write.table(micro.GO,"./Hippo\\subcluster\\Microglial\\diff\\micro.GO.txt",sep = "\t")
table(micro.GO$Group,micro.GO$Celltype)

micro5_p <- micro.GO[micro.GO$Group %in% "Micro_P" & micro.GO$Celltype %in% "Micro5",]
micro5_a <- micro.GO[micro.GO$Group %in% "Micro_A" & micro.GO$Celltype %in% "Micro5",]


###***星形胶质细胞每个亚群间差异基因--------
Astrocyte <- readRDS(file = "./Hippo\\subcluster\\Astrocyte\\anno\\Astrocyte_sub.rds")
celltype <- names(table(Astrocyte$Anno_Idents))
diff_Astrocyte_all <- data.frame(matrix(ncol = 5, nrow = 0))
for (i in 1:length(celltype)) {
  data <- subset(Astrocyte,Anno_Idents %in% celltype[i])
  diff_celltype <- FindMarkers(data, min.pct = 0.1, 
                               logfc.threshold = 0.25,
                               test.use = "wilcox",
                               group.by = "group",
                               ident.1 ="P",
                               ident.2="A")
  
  diff_celltype$celltype <- rep(celltype[i],dim(diff_celltype)[1])
  diff_Astrocyte_all <- rbind.data.frame(diff_Astrocyte_all,diff_celltype)
}
write.table(diff_Astrocyte_all,"./Hippo\\subcluster\\Astrocyte\\diff\\diff_Astrocyte_all_marker.txt",sep = "\t")

celltype <- unique(diff_Astrocyte_all$celltype)
Astrocyte.GO <- data.frame(matrix(ncol = 10, nrow = 0))
for (i in 1:length(celltype)) {
  ##在P端上调
  up_exp <- diff_Astrocyte_all[diff_Astrocyte_all$celltype %in% celltype[i] & diff_Astrocyte_all$p_val<0.05 & diff_Astrocyte_all$avg_log2FC>0.25,]
  up <- rownames(up_exp)
  up.enrich.go <- enrichGO(gene = up,  #基因列表文件中的基因名称
                           OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                           keyType = 'SYMBOL', 
                           ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                           pAdjustMethod = 'BH',  #指定 p 值校正方法
                           pvalueCutoff = 1,  #指定 p 值阈值，不显著的值将不显示在结果中
                           qvalueCutoff = 1,  #指定 q 值阈值，不显著的值将不显示在结果中
                           readable = FALSE)
  up.enrich.go <-as.data.frame(up.enrich.go)  #
  up.enrich.go <- up.enrich.go[up.enrich.go$pvalue <0.01,]
  up.enrich.go$Group <- rep("Astro_P",dim(up.enrich.go)[1])
  up.enrich.go$Celltype <- rep(celltype[i],dim(up.enrich.go)[1])
  Astrocyte.GO <- rbind.data.frame(Astrocyte.GO,up.enrich.go)
  ##在A端上调
  down_exp <- diff_Astrocyte_all[diff_Astrocyte_all$celltype %in% celltype[i] & diff_Astrocyte_all$p_val<0.05 & diff_Astrocyte_all$avg_log2FC<(-0.25),]
  down <- rownames(down_exp)
  down.enrich.go <- enrichGO(gene = down,  #基因列表文件中的基因名称
                             OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                             keyType = 'SYMBOL',  
                             ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                             pAdjustMethod = 'BH',  #指定 p 值校正方法
                             pvalueCutoff = 1,  #指定 p 值阈值，不显著的值将不显示在结果中
                             qvalueCutoff = 1,  #指定 q 值阈值，不显著的值将不显示在结果中
                             readable = FALSE)
  down.enrich.go <-as.data.frame(down.enrich.go)  #
  down.enrich.go <- down.enrich.go[down.enrich.go$pvalue <0.01,]
  down.enrich.go$Group <- rep("Astro_A",dim(down.enrich.go)[1])
  down.enrich.go$Celltype <- rep(celltype[i],dim(down.enrich.go)[1])
  Astrocyte.GO <- rbind(Astrocyte.GO,down.enrich.go)
  
}
write.table(Astrocyte.GO,"./Hippo\\subcluster\\Astrocyte\\diff\\Astrocyte.GO.txt",sep = "\t")

table(Astrocyte.GO$Group,Astrocyte.GO$Celltype)

AST5_p <- Astrocyte.GO[Astrocyte.GO$Group %in% "Astro_P" & Astrocyte.GO$Celltype %in% "AST5",]
AST5_a <- Astrocyte.GO[Astrocyte.GO$Group %in% "Astro_A" & Astrocyte.GO$Celltype %in% "AST5",]



######***少突胶质细胞每个亚群间差异基因--------
Oligodendrocyte <- readRDS(file = "./Hippo\\subcluster\\Oligodendrocyte\\anno\\Oligodendrocyte_sub.rds")
celltype <- names(table(Oligodendrocyte$Anno_Idents))
diff_Oligodendrocyte_all <- data.frame(matrix(ncol = 5, nrow = 0))
for (i in 1:length(celltype)) {
  data <- subset(Oligodendrocyte,Anno_Idents %in% celltype[i])
  diff_celltype <- FindMarkers(data, min.pct = 0.1, 
                               logfc.threshold = 0.25,
                               test.use = "wilcox",
                               group.by = "group",
                               ident.1 ="P",
                               ident.2="A")
  
  diff_celltype$celltype <- rep(celltype[i],dim(diff_celltype)[1])
  diff_Oligodendrocyte_all <- rbind.data.frame(diff_Oligodendrocyte_all,diff_celltype)
}
write.table(diff_Oligodendrocyte_all,"./Hippo\\subcluster\\Oligodendrocyte\\diff\\diff_Oligodendrocyte_all_marker.txt",sep = "\t")

celltype <- unique(diff_Oligodendrocyte_all$celltype)
Oligodendrocyte.GO <- data.frame(matrix(ncol = 10, nrow = 0))
for (i in 1:length(celltype)) {
  ##在P端上调
  up_exp <- diff_Oligodendrocyte_all[diff_Oligodendrocyte_all$celltype %in% celltype[i] & diff_Oligodendrocyte_all$p_val<0.05 & diff_Oligodendrocyte_all$avg_log2FC>0.25,]
  up <- rownames(up_exp)
  up.enrich.go <- enrichGO(gene = up,  #基因列表文件中的基因名称
                           OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                           keyType = 'SYMBOL',  
                           ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                           pAdjustMethod = 'BH',  #指定 p 值校正方法
                           pvalueCutoff = 1,  #指定 p 值阈值，不显著的值将不显示在结果中
                           qvalueCutoff = 1,  #指定 q 值阈值，不显著的值将不显示在结果中
                           readable = FALSE)
  up.enrich.go <-as.data.frame(up.enrich.go)  #
  up.enrich.go <- up.enrich.go[up.enrich.go$pvalue <0.01,]
  up.enrich.go$Group <- rep("Oli_P",dim(up.enrich.go)[1])
  up.enrich.go$Celltype <- rep(celltype[i],dim(up.enrich.go)[1])
  Oligodendrocyte.GO <- rbind.data.frame(Oligodendrocyte.GO,up.enrich.go)
  ##在A端上调
  down_exp <- diff_Oligodendrocyte_all[diff_Oligodendrocyte_all$celltype %in% celltype[i] & diff_Oligodendrocyte_all$p_val<0.05 & diff_Oligodendrocyte_all$avg_log2FC<(-0.25),]
  down <- rownames(down_exp)
  down.enrich.go <- enrichGO(gene = down,  #基因列表文件中的基因名称
                             OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                             keyType = 'SYMBOL',  
                             ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                             pAdjustMethod = 'BH',  #指定 p 值校正方法
                             pvalueCutoff = 1,  #指定 p 值阈值，不显著的值将不显示在结果中
                             qvalueCutoff = 1,  #指定 q 值阈值，不显著的值将不显示在结果中
                             readable = FALSE)
  down.enrich.go <-as.data.frame(down.enrich.go)  #
  down.enrich.go <- down.enrich.go[down.enrich.go$pvalue <0.01,]
  down.enrich.go$Group <- rep("Oli_A",dim(down.enrich.go)[1])
  down.enrich.go$Celltype <- rep(celltype[i],dim(down.enrich.go)[1])
  Oligodendrocyte.GO <- rbind(Oligodendrocyte.GO,down.enrich.go)
  
}
write.table(Oligodendrocyte.GO,"./Hippo\\subcluster\\Oligodendrocyte\\diff\\Oligodendrocyte.GO.txt",sep = "\t")

table(Oligodendrocyte.GO$Group,Oligodendrocyte.GO$Celltype)

OPCs_p <- Oligodendrocyte.GO[Oligodendrocyte.GO$Group %in% "Oli_P" & Oligodendrocyte.GO$Celltype %in% "OPCs",]
OPCs_a <- Oligodendrocyte.GO[Oligodendrocyte.GO$Group %in% "Oli_A" & Oligodendrocyte.GO$Celltype %in% "OPCs",]

####每个亚群top5功能合并---------
astro_cols <- c("#3A6963","#BC80BD","#E59CC4","#cabbe9","#58A4C3")
micro_cols<- c("#8F797E","#FFC2B5","#FFE3CC","#646C8F","#DCC3A1",)
oli_cols<- c('#DB9A56','#8DBEB2','#C8B985','#B5D673','#803736')

Oligodendrocyte.GO<- read.table("./Hippo\\subcluster\\Oligodendrocyte\\diff\\Oligodendrocyte.GO.txt",sep = "\t",header = T)
micro.GO<- read.table("./Hippo\\subcluster\\Microglial\\diff\\micro.GO.txt",sep = "\t",header = T)
Astrocyte.GO<- read.table("./Hippo\\subcluster\\Astrocyte\\diff\\Astrocyte.GO.txt",sep = "\t",header = T)

A_Oligodendrocyte.GO <- Oligodendrocyte.GO[Oligodendrocyte.GO$Group %in% "Oli_A",]
A_Oligodendrocyte.GO$logp<- -log10(A_Oligodendrocyte.GO$pvalue)
A_Oligodendrocyte.GO <- A_Oligodendrocyte.GO %>%
  dplyr::top_n(n = 5, wt = logp)


A_micro.GO <- micro.GO[micro.GO$Group %in% "Micro_A",]
A_micro.GO$logp<- -log10(A_micro.GO$pvalue)
A_micro.GO <- A_micro.GO %>%
  dplyr::top_n(n = 5, wt = logp)

A_Astrocyte.GO <- Astrocyte.GO[Astrocyte.GO$Group %in% "Astro_A",]
A_Astrocyte.GO$logp<- -log10(A_Astrocyte.GO$pvalue)
A_Astrocyte.GO <- A_Astrocyte.GO %>%
  dplyr::top_n(n = 5, wt = logp)


Oli_micro <- intersect(A_Oligodendrocyte.GO$Description,A_micro.GO$Description)
Oli_micro_astro <- intersect(Oli_micro,A_Astrocyte.GO$Description)


###***在海马体前端共同的功能----------
A_Oligodendrocyte.GO <- Oligodendrocyte.GO[Oligodendrocyte.GO$Group %in% "Oli_A",]
table(A_Oligodendrocyte.GO$Celltype)
A_Oligodendrocyte.GO$logp<- -log10(A_Oligodendrocyte.GO$pvalue)
A_Oligodendrocyte.GO <- A_Oligodendrocyte.GO %>%
  dplyr::group_by(Celltype) %>%
  dplyr::top_n(n = 5, wt = logp)




A_micro.GO <- micro.GO[micro.GO$Group %in% "Micro_A",]
table(A_micro.GO$Celltype)
A_micro.GO$logp<- -log10(A_micro.GO$pvalue)
A_micro.GO <- A_micro.GO %>%
  dplyr::group_by(Celltype) %>%
  dplyr::top_n(n = 5, wt = logp)




A_Astrocyte.GO <- Astrocyte.GO[Astrocyte.GO$Group %in% "Astro_A",]
table(A_Astrocyte.GO$Celltype)
A_Astrocyte.GO$logp<- -log10(A_Astrocyte.GO$pvalue)
A_Astrocyte.GO <- A_Astrocyte.GO %>%
  dplyr::group_by(Celltype) %>%
  dplyr::top_n(n = 5, wt = logp)


Oli_micro <- intersect(A_Oligodendrocyte.GO$Description,A_micro.GO$Description)
Oli_micro_astro <- intersect(Oli_micro,A_Astrocyte.GO$Description)
#"cell junction assembly"

A_go.all <- rbind.data.frame(A_Astrocyte.GO,A_micro.GO,A_Oligodendrocyte.GO)
length(unique(A_go.all$Description))
#64 
write.table(A_go.all,"./Hippo\\subcluster\\enrich\\A_go.all_top5.go.txt",sep = '\t',col.names = T,row.names = F,na='')
A_go.all$ID[A_go.all$Description %in% "cell junction assembly"]

A_go.gial <- A_go.all %>%
  dplyr::group_by(Celltype) %>%
  dplyr::top_n(n = 5, wt = logp)
###***在海马体后端共同的功能----------
P_Oligodendrocyte.GO <- Oligodendrocyte.GO[Oligodendrocyte.GO$Group %in% "Oli_P",]
table(P_Oligodendrocyte.GO$Celltype)
P_Oligodendrocyte.GO$logp<- -log10(P_Oligodendrocyte.GO$pvalue)
P_Oligodendrocyte.GO <- P_Oligodendrocyte.GO %>%
  dplyr::group_by(Celltype) %>%
  dplyr::top_n(n = 5, wt = logp)




P_micro.GO <- micro.GO[micro.GO$Group %in% "Micro_P",]
table(P_micro.GO$Celltype)
P_micro.GO$logp<- -log10(P_micro.GO$pvalue)
P_micro.GO <- P_micro.GO %>%
  dplyr::group_by(Celltype) %>%
  dplyr::top_n(n = 5, wt = logp)




P_Astrocyte.GO <- Astrocyte.GO[Astrocyte.GO$Group %in% "Astro_P",]
table(P_Astrocyte.GO$Celltype)
P_Astrocyte.GO$logp<- -log10(P_Astrocyte.GO$pvalue)
P_Astrocyte.GO <- P_Astrocyte.GO %>%
  dplyr::group_by(Celltype) %>%
  dplyr::top_n(n = 5, wt = logp)


Oli_micro_P <- intersect(P_Oligodendrocyte.GO$Description,P_micro.GO$Description)
Oli_micro_astro_P <- intersect(Oli_micro_P,P_Astrocyte.GO$Description)
P_go.all$ID[P_go.all$Description %in% "regulation of trans-synaptic signaling"]
#[1] "GO:0050808":"synapse organization"                         "GO:0050804":"modulation of chemical synaptic transmission"
#[3] "GO:0099177":"regulation of trans-synaptic signaling"   

P_go.all <- rbind.data.frame(P_Astrocyte.GO,P_micro.GO,P_Oligodendrocyte.GO)
length(unique(P_go.all$Description))

write.table(P_go.all,"./Hippo\\subcluster\\enrich\\P_go.all_top5.go.txt",sep = '\t',col.names = T,row.names = F,na='')

synapse_organization <- P_go.all[P_go.all$Description %in% "synapse organization",]

####海马体胶质细胞细胞类型海马体前后端差异基因功能富集点图+热图-----------------
astro_cols <- c("#3A6963","#BC80BD","#E59CC4","#cabbe9","#58A4C3")
micro_cols<- c("#8F797E","#FFC2B5","#FFE3CC","#646C8F","#DCC3A1",)
oli_cols<- c('#DB9A56','#8DBEB2','#C8B985','#B5D673','#803736')
####******构建矩阵--------
A_go.all <- read.table("./Hippo\\subcluster\\enrich\\A_go.all_top5.go.txt",sep = '\t',header = T)
table(A_go.all$Celltype)
Celltype <- unique(A_go.all$Celltype)
exp <- A_go.all[A_go.all$Celltype %in% Celltype[1],c(2,9,12)]
colnames(exp) <- c("Description",paste(Celltype[1],"A_count",sep = "_"),paste(Celltype[1],"A_logP",sep = "_"))
library(plyr)
for (i in 2:length(Celltype)) {
  class_go <- A_go.all[A_go.all$Celltype %in% Celltype[i],c(2,9,12)]
  colnames(class_go) <- c("Description",paste(Celltype[i],"A_count",sep = "_"),paste(Celltype[i],"A_logP",sep = "_"))
  exp <- merge(exp,class_go,by="Description",all=TRUE)
}

P_go.all <- read.table("./Hippo\\subcluster\\enrich\\P_go.all_top5.go.txt",sep = '\t',header = T)
Celltype <- unique(P_go.all$Celltype)
for (i in 1:length(Celltype)) {
  class_go <- P_go.all[P_go.all$Celltype %in% Celltype[i],c(2,9,12)]
  colnames(class_go) <- c("Description",paste(Celltype[i],"P_count",sep = "_"),paste(Celltype[i],"P_logP",sep = "_"))
  exp <- merge(exp,class_go,by="Description",all=TRUE)
}

exp_1 <- exp[,c(2:ncol(exp))]
rownames(exp_1) <- exp$Description
exp <- exp_1
spot.size.function<-function(x){
  x= x/4
}
astro_cols <- c("#3A6963","#BC80BD","#E59CC4","#cabbe9","#58A4C3")
micro_cols<- c("#8F797E","#FFC2B5","#FFE3CC","#646C8F","#DCC3A1")
oli_cols<- c('#DB9A56','#8DBEB2','#C8B985','#B5D673','#803736')

spot.colour.function <- function(x) {
  colours <- rep("white", length(x));
  colours[x==exp$AST1_A_logP] <- "#3A6963"; 
  colours[x==exp$AST2_A_logP] <- "#BC80BD"; 
  colours[x==exp$AST3_A_logP] <- "#E59CC4";
  colours[x==exp$AST4_A_logP] <- "#cabbe9";
  colours[x==exp$AST5_A_logP] <- "#58A4C3";
  
  colours[x==exp$Micro1_A_logP] <- "#8F797E"; 
  colours[x==exp$Micro2_A_logP] <- "#FFC2B5"; 
  colours[x==exp$Micro3_A_logP] <- "#FFE3CC";
  colours[x==exp$Micro4_A_logP] <- "#646C8F";
  colours[x==exp$Micro5_A_logP] <- "#DCC3A1";
  
  colours[x==exp$mOli_A_logP] <- "#DB9A56"; 
  colours[x==exp$imOli_A_logP] <- "#8DBEB2"; 
  colours[x==exp$OPCs_A_logP] <- "#C8B985";
  colours[x==exp$Oli1_A_logP] <- "#B5D673";
  colours[x==exp$Oli2_A_logP] <- "#803736";
  
  colours[x==exp$AST1_P_logP] <- "#3A6963"; 
  colours[x==exp$AST2_P_logP] <- "#BC80BD"; 
  colours[x==exp$AST3_P_logP] <- "#E59CC4";
  colours[x==exp$AST4_P_logP] <- "#cabbe9";
  colours[x==exp$AST5_P_logP] <- "#58A4C3";
  
  colours[x==exp$Micro1_P_logP] <- "#8F797E"; 
  colours[x==exp$Micro2_P_logP] <- "#FFC2B5"; 
  colours[x==exp$Micro3_P_logP] <- "#FFE3CC";
  colours[x==exp$Micro4_P_logP] <- "#646C8F";
  colours[x==exp$Micro5_P_logP] <- "#DCC3A1";
  
  colours[x==exp$mOli_P_logP] <- "#DB9A56"; 
  colours[x==exp$imOli_P_logP] <- "#8DBEB2"; 
  colours[x==exp$OPCs_P_logP] <- "#C8B985";
  colours[x==exp$Oli1_P_logP] <- "#B5D673";
  colours[x==exp$Oli2_P_logP] <- "#803736";
  return(colours);
}
#加图例
key.trans <- list(title = "cell type",space = "right",columns = 1,
                  points = list(pch = c(rep(20,15)),
                                col = c(astro_cols,micro_cols,oli_cols),
                                cex=c(1,1,1,1)),
                  text = list(c(Celltype)),
                  
                  #lines = list(col = colors,lty = lines),
                  cex.title = 1,cex = .9)
library(BoutrosLab.plotting.general)
library(ggplot2)
create.dotmap(
  exp[,c(seq(from=2, to=60, by=2))],bg.data = exp[,c(seq(from=1, to=60, by=2))],
  pch = 20,na.spot.size=0,spot.size.function=spot.size.function,
  spot.colour.function=spot.colour.function,
  ylab.cex=1,yaxis.cex=0.5,xaxis.cex=0.5,
  colour.scheme=colorRampPalette(c("#80B1D3","#FDB462","#E59CC4", "#BC80BD"))(100),
  key = key.trans,total.colours = 10,colourkey = T,xaxis.rot = 45)

#########top5功能个数展示upset--------
A_go.all <- read.table("./Hippo\\subcluster\\enrich\\A_go.all_top5.go.txt",sep = '\t',header = T)
aa <- as.data.frame(table(A_go.all$Description,A_go.all$Group))
aa <- spread(aa,Var2,Freq)
rownames(aa) <- aa$Var1
aa <- aa[,-1]
aa[aa != 0] <- 1
GOTerms <- colnames(aa)
#devtools::install_github("krassowski/complex-upset")
library(ComplexUpset)
upset(aa,GOTerms,width_ratio = 0.1,
      base_annotations = list(
        "intersection" = intersection_size(
          counts = T,
          mapping = aes(fill="bars_color")
        ) 
        + scale_fill_manual(values = c("bars_color"="#99CCCC"),guide="none")
      )
)
ggsave("./Hippo\\subcluster\\enrich\\top5_A端三类胶质细胞GO交集upset.pdf",height = 5,width = 5)
overlap <- intersect(A_go.all$Description[A_go.all$Group %in% "Micro_A"],A_go.all$Description[A_go.all$Group %in% "Astro_A"])
overlap <- intersect(overlap,A_go.all$Description[A_go.all$Group %in% "Oli_A"])
#"cell junction assembly"

P_go.all <- read.table("./Hippo\\subcluster\\enrich\\P_go.all_top5.go.txt",sep = '\t',header = T)
aa <- as.data.frame(table(P_go.all$Description,P_go.all$Group))
aa <- spread(aa,Var2,Freq)
rownames(aa) <- aa$Var1
aa <- aa[,-1]
aa[aa != 0] <- 1
GOTerms <- colnames(aa)
#devtools::install_github("krassowski/complex-upset")
library(ComplexUpset)
upset(aa,GOTerms,width_ratio = 0.1,
      base_annotations = list(
        "intersection" = intersection_size(
          counts = T,
          mapping = aes(fill="bars_color")
        ) 
        + scale_fill_manual(values = c("bars_color"="#FFCC99"),guide="none")
      )
)
ggsave("./Hippo\\subcluster\\enrich\\top5_P端三类胶质细胞GO交集upset.pdf",height = 5,width = 5)
overlap <- intersect(P_go.all$Description[P_go.all$Group %in% "Micro_P"],P_go.all$Description[P_go.all$Group %in% "Astro_P"])
overlap <- intersect(overlap,P_go.all$Description[P_go.all$Group %in% "Oli_P"])
# [1] "synapse organization"                         "modulation of chemical synaptic transmission"
# [3] "regulation of trans-synaptic signaling" 
#########将这四个功能与基因关系展示桑基图-----------
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cowplot)
library(cols4all)

df <- read.table("./Hippo\\subcluster\\all_spec_target.txt",header = T,sep = "\t")
A_go.all <- read.table("./Hippo\\subcluster\\enrich\\A_go.all_top5.go.txt",sep = '\t',header = T)
P_go.all <- read.table("./Hippo\\subcluster\\enrich\\P_go.all_top5.go.txt",sep = '\t',header = T)
A_cell_junction_assembly <- A_go.all[A_go.all$Description %in% "cell junction assembly",]
exp <- data.frame(matrix(ncol = 2,nrow = 0))
colnames(exp) <- c("genes","Go")
for (i in 1:length(A_cell_junction_assembly$Description)) {
  genes <- A_cell_junction_assembly[i,8]
  genes <- str_split(genes,"/")[[1]]
  genes <- as.data.frame(genes)
  genes$Go <- A_cell_junction_assembly[i,2]
  genes$Celltype <- A_cell_junction_assembly[i,11]
  exp <- rbind(exp,genes)
}
intersect(exp$genes,df$TargetGene[df$Group %in% c("Astro_A","Micro_A","Oligodendrocyte_A")])
unique(df$TF[df$TargetGene %in% intersect(exp$genes,df$TargetGene[df$Group %in% c("Astro_A","Micro_A","Oligodendrocyte_A")]) & df$Group %in% c("Astro_A","Micro_A","Oligodendrocyte_A")])
a <- as.data.frame(table(exp$genes))
exp_A <- exp[exp$genes %in% a$Var1[a$Freq >1],]
table(a$Freq)
P_selected <- P_go.all[P_go.all$Description %in% c("synapse organization",
                                                   "modulation of chemical synaptic transmission",
                                                   "regulation of trans-synaptic signaling"),]
exp_p <- data.frame(matrix(ncol = 2,nrow = 0))
colnames(exp_p) <- c("genes","Go")
for (i in 1:length(P_selected$Description)) {
  genes <- P_selected[i,8]
  genes <- str_split(genes,"/")[[1]]
  genes <- as.data.frame(genes)
  genes$Go <- P_selected[i,2]
  genes$Celltype <- P_selected[i,11]
  exp_p <- rbind(exp_p,genes)
}
a <- as.data.frame(table(exp_p$genes))
table(a$Freq)
exp_P <- exp_p[exp_p$genes %in% a$Var1[a$Freq >7],]

exp_all <- rbind(exp_A,exp_P)
write.table(exp_all,"./Hippo\\subcluster\\enrich\\筛选出的四个功能-基因.txt",row.names = T,sep = "\t")
exp_all <-read.table("./Hippo\\subcluster\\enrich\\筛选出的四个功能-基因.txt",header = T,sep = "\t")

###桑葚图----
df_a <- exp_all %>%
  make_long(genes,Go, Celltype)

head(df_a)
df_a$node <- factor(df_a$node,levels = c(exp_all$Celltype %>% unique()%>% rev(),
                                         exp_all$genes %>% unique() %>% rev(),
                                         exp_all$Go %>% unique() %>% rev()))
#自定义配色：
#绘图：
ggplot(df_a, aes(x = x,
                 next_x = next_x,
                 node = node,
                 next_node = next_node,
                 fill = node,
                 label = node)) +scale_fill_manual(values = mycol)+
  geom_sankey(flow.alpha = 0.5,
              
              smooth = 8,
              width = 0.08) +scale_fill_manual(values = mycol)+
  geom_sankey_text(size = 3.2,
                   color = "black")+
  theme_void() +
  theme(legend.position = 'none')
#筛选出的四个功能-基因-细胞亚簇桑基图
########全部功能个数展示upset--------------
Oligodendrocyte.GO<- read.table("./Hippo\\subcluster\\Oligodendrocyte\\diff\\Oligodendrocyte.GO.txt",sep = "\t",header = T)
micro.GO<- read.table("./Hippo\\subcluster\\Microglial\\diff\\micro.GO.txt",sep = "\t",header = T)
Astrocyte.GO<- read.table("./Hippo\\subcluster\\Astrocyte\\diff\\Astrocyte.GO.txt",sep = "\t",header = T)
###A端------------
A_Oligodendrocyte.GO <- Oligodendrocyte.GO[Oligodendrocyte.GO$Group %in% "Oli_A",]
A_micro.GO <- micro.GO[micro.GO$Group %in% "Micro_A",]
A_Astrocyte.GO <- Astrocyte.GO[Astrocyte.GO$Group %in% "Astro_A",]
library(plyr)
aaa <- rbind.fill(A_Astrocyte.GO,A_micro.GO,A_Oligodendrocyte.GO)
aa <- as.data.frame(table(aaa$Description,aaa$Group))
aa <- spread(aa,Var2,Freq)
rownames(aa) <- aa$Var1
aa <- aa[,-1]
aa[aa != 0] <- 1
GOTerms <- colnames(aa)
#devtools::install_github("krassowski/complex-upset")
library(ComplexUpset)
upset(aa,GOTerms,width_ratio = 0.1,
      base_annotations = list(
        "intersection" = intersection_size(
          counts = T,
          mapping = aes(fill="bars_color")
        ) 
        + scale_fill_manual(values = c("bars_color"="#99CCCC"),guide="none")
      )
)
ggsave("./Hippo\\subcluster\\enrich\\A端三类胶质细胞GO交集upset.pdf",height = 5,width = 5)
overlap1 <- intersect(A_Astrocyte.GO$Description,A_Oligodendrocyte.GO$Description)
overlap1 <- intersect(overlap1,A_micro.GO$Description)
####细胞类型分开------
celltype <- unique(aaa$Group)
for(i in 1:length(celltype)){
  aa <- as.data.frame(table(aaa[aaa$Group %in% celltype[i],2],aaa[aaa$Group %in% celltype[i],11]))
  aa <- spread(aa,Var2,Freq)
  rownames(aa) <- aa$Var1
  aa <- aa[,-1]
  aa[aa != 0] <- 1
  GOTerms <- colnames(aa)
  #devtools::install_github("krassowski/complex-upset")
  library(ComplexUpset)
  upset(aa,GOTerms,width_ratio = 0.1,
        base_annotations = list(
          "intersection" = intersection_size(
            counts = T,
            mapping = aes(fill="bars_color")
          ) 
          + scale_fill_manual(values = c("bars_color"="#99CCCC"),guide="none")
        )
  )
  setwd("./Hippo\\subcluster\\enrich")
  names <- paste(celltype[i],"亚簇GO交集upset.pdf",sep = "")
  ggsave(names,height = 5,width = 5)
  
}
###P端------------
P_Oligodendrocyte.GO <- Oligodendrocyte.GO[Oligodendrocyte.GO$Group %in% "Oli_P",]
P_micro.GO <- micro.GO[micro.GO$Group %in% "Micro_P",]
P_Astrocyte.GO <- Astrocyte.GO[Astrocyte.GO$Group %in% "Astro_P",]

aaa <- rbind.fill(P_Astrocyte.GO,P_micro.GO,P_Oligodendrocyte.GO)
aa <- as.data.frame(table(aaa$Description,aaa$Celltype))
aa <- spread(aa,Var2,Freq)
rownames(aa) <- aa$Var1
aa <- aa[,-1]
aa[aa != 0] <- 1
GOTerms <- colnames(aa)
#devtools::install_github("krassowski/complex-upset")
library(ComplexUpset)
upset(aa,GOTerms,width_ratio = 0.1,
      base_annotations = list(
        "intersection" = intersection_size(
          counts = T,
          mapping = aes(fill="bars_color")
        ) 
        + scale_fill_manual(values = c("bars_color"="#FFCC99"),guide="none")
      )
)
ggsave("./Hippo\\subcluster\\enrich\\P端三类胶质细胞亚簇GO交集upset.pdf",height = 10,width = 20)
overlap <- intersect(P_Astrocyte.GO$Description,P_Oligodendrocyte.GO$Description)
overlap <- intersect(overlap,P_micro.GO$Description) #221
####AP端都富集到的功能--------------
all <- intersect(overlap,overlap1)
####细胞类型分开------
celltype <- unique(aaa$Group)
for(i in 1:length(celltype)){
  aa <- as.data.frame(table(aaa[aaa$Group %in% celltype[i],2],aaa[aaa$Group %in% celltype[i],11]))
  aa <- spread(aa,Var2,Freq)
  rownames(aa) <- aa$Var1
  aa <- aa[,-1]
  aa[aa != 0] <- 1
  GOTerms <- colnames(aa)
  #devtools::install_github("krassowski/complex-upset")
  library(ComplexUpset)
  upset(aa,GOTerms,width_ratio = 0.1,
        base_annotations = list(
          "intersection" = intersection_size(
            counts = T,
            mapping = aes(fill="bars_color")
          ) 
          + scale_fill_manual(values = c("bars_color"="#FFCC99"),guide="none")
        )
  )
  setwd("./Hippo\\subcluster\\enrich")
  names <- paste(celltype[i],"亚簇GO交集upset.pdf",sep = "")
  ggsave(names,height = 5,width = 5)
  
}
