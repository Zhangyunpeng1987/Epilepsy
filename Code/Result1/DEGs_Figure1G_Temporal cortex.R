####****差异分析以及功能富集--------------
####*****每个细胞类型在不同样本之间差异分析--------
scRNA <- readRDS(file = ".\\scRNA_anno_cluster.rds")
metadata <- scRNA@meta.data
##meta.data添加信息
group.id <- as.data.frame(metadata[,4])
rownames(group.id) <- rownames(metadata)
group.id$group <- t(as.data.frame(stringr::str_extract_all(group.id$`metadata[, 4]`, '\\D+')))
group <- factor(group.id[,2])
class(group)
scRNA <- AddMetaData(scRNA, group,col.name = "group")
dim(scRNA)  #28306 112671
saveRDS(scRNA,file = ".\\scRNA_anno_cluster.rds")
####******Astrocyte----------------
scRNA <- readRDS(file = ".\\scRNA_anno_cluster.rds")
Astrocyte <- subset(scRNA, Anno_Idents=="Astrocyte")
diff_Astrocyte <- FindMarkers(Astrocyte, min.pct = 0.1, 
                              logfc.threshold = 0.25,
                              test.use = "wilcox",
                              group.by = "group",
                              ident.1 ="Ep",
                              ident.2="Nor")
diff_Astrocyte<-cbind(rownames(diff_Astrocyte),diff_Astrocyte)  ##对数据增加一列
colnames(diff_Astrocyte)<- c('gene_id',"pvalue" ,"avglog2FoldChange","pct.1" ,"pct.2","padj" )
##对列名进行修改  1600
write.table(diff_Astrocyte,".\\wilcox\\case-control-astro-wilcox.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

####******Microglial----------------
levels(scRNA)
Microglial <- subset(scRNA, Anno_Idents=="Microglial")
diff_Microglial <- FindMarkers(Microglial, min.pct = 0.1, 
                               logfc.threshold = 0.25,
                               test.use = "wilcox",
                               group.by = "group",
                               ident.1 ="Ep",
                               ident.2="Nor")
diff_Microglial<-cbind(rownames(diff_Microglial),diff_Microglial)  ##对数据增加一列
colnames(diff_Microglial)<- c('gene_id',"pvalue" ,"avglog2FoldChange","pct.1" ,"pct.2","padj" )
##对列名进行修改   5342
write.table(diff_Microglial,".\\wilcox\\case-control-Microglial-wilcox.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')


####******Oligodendrocyte precursor cell----------------
levels(scRNA)
Oligodendrocyte_precursor_cell <- subset(scRNA, Anno_Idents=="Oligodendrocyte precursor cell")
diff_Oligodendrocyte_precursor_cell <- FindMarkers(Oligodendrocyte_precursor_cell, min.pct = 0.1, 
                                                   logfc.threshold = 0.25,
                                                   test.use = "wilcox",
                                                   group.by = "group",
                                                   ident.1 ="Ep",
                                                   ident.2="Nor")
diff_Oligodendrocyte_precursor_cell<-cbind(rownames(diff_Oligodendrocyte_precursor_cell),diff_Oligodendrocyte_precursor_cell)  ##对数据增加一列
colnames(diff_Oligodendrocyte_precursor_cell)<- c('gene_id',"pvalue" ,"avglog2FoldChange","pct.1" ,"pct.2","padj" )
##对列名进行修改  4894
write.table(diff_Oligodendrocyte_precursor_cell,".\\wilcox\\case-control-Oligodendrocyte_precursor_cell-wilcox.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

####******Endothelial----------------
levels(scRNA)
Endothelial <- subset(scRNA, Anno_Idents=="Endothelial")
diff_Endothelial <- FindMarkers(Endothelial, min.pct = 0.1,
                                logfc.threshold = 0.25,
                                test.use = "wilcox",
                                group.by = "group",
                                ident.1 ="Ep",
                                ident.2="Nor")
diff_Endothelial<-cbind(rownames(diff_Endothelial),diff_Endothelial)  ##对数据增加一列
colnames(diff_Endothelial)<- c('gene_id',"pvalue" ,"avglog2FoldChange","pct.1" ,"pct.2","padj" )
##对列名进行修改  5557
write.table(diff_Endothelial,".\\wilcox\\case-control-Endothelial-wilcox.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

####******Oligodendrocyte-----
levels(scRNA)
Oligodendrocyte <- subset(scRNA, Anno_Idents=="Oligodendrocyte")
diff_Oligodendrocyte <- FindMarkers(Oligodendrocyte, min.pct = 0.1, 
                                    logfc.threshold = 0.25,
                                    test.use = "wilcox",
                                    group.by = "group",
                                    ident.1 ="Ep",
                                    ident.2="Nor")
diff_Oligodendrocyte<-cbind(rownames(diff_Oligodendrocyte),diff_Oligodendrocyte)  ##对数据增加一列
colnames(diff_Oligodendrocyte)<- c('gene_id',"pvalue" ,"avglog2FoldChange","pct.1" ,"pct.2","padj" )
##对列名进行修改  1808
write.table(diff_Oligodendrocyte,".\\wilcox\\case-control-Oligodendrocyte-wilcox.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

####******Inhibitory neuron----------------
levels(scRNA)
Inhibitory_neuron <- subset(scRNA, Anno_Idents=="Inhibitory neuron")
diff_Inhibitory_neuron <- FindMarkers(Inhibitory_neuron, min.pct = 0.1, 
                                      logfc.threshold = 0.25,
                                      test.use = "wilcox",
                                      group.by = "group",
                                      ident.1 ="Ep",
                                      ident.2="Nor")
diff_Inhibitory_neuron<-cbind(rownames(diff_Inhibitory_neuron),diff_Inhibitory_neuron)  ##对数据增加一列
colnames(diff_Inhibitory_neuron)<- c('gene_id',"pvalue" ,"avglog2FoldChange","pct.1" ,"pct.2","padj" )
##对列名进行修改
write.table(diff_Inhibitory_neuron,".\\wilcox\\case-control-Inhibitory_neuron-wilcox.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

####******Excitatory neuron----------------
levels(scRNA)
Excitatory_neuron <- subset(scRNA, Anno_Idents=="Excitatory neuron")
diff_Excitatory_neuron <- FindMarkers(Excitatory_neuron, min.pct = 0, 
                                      logfc.threshold = 0,
                                      test.use = "wilcox",
                                      group.by = "group",
                                      ident.1 ="Ep",
                                      ident.2="Nor")
diff_Excitatory_neuron<-cbind(rownames(diff_Excitatory_neuron),diff_Excitatory_neuron)  ##对数据增加一列
colnames(diff_Excitatory_neuron)<- c('gene_id',"pvalue" ,"avglog2FoldChange","pct.1" ,"pct.2","padj" )
##对列名进行修改
write.table(diff_Excitatory_neuron,".\\wilcox\\case-control-Excitatory_neuron-wilcox.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

####*****功能富集------------
####******Endo圈图------------
diff_Endothelial <- read.table(".\\wilcox\\case-control-Endothelial-wilcox.txt",sep = '\t',header = T)
diff_Endothelial <- diff_Endothelial[which(diff_Endothelial$pvalue<0.01),]
Endo.gene <- diff_Endothelial[,1]

#222
#对于加载的注释库的使用，以上述为例，就直接在 OrgDb 中指定人（org.Hs.eg.db）或绵羊（sheep）
Endo.enrich.go <- enrichGO(gene = Endo.gene,  #基因列表文件中的基因名称
                           OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                           keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                           ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                           pAdjustMethod = 'fdr',  #指定 p 值校正方法
                           pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                           qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                           readable = FALSE)
Endo.enrich.go <-summary(Endo.enrich.go)  #11
write.table(Endo.enrich.go,".\\enrich\\Endo.enrich.go.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')


DE.Endo <- diff_Endothelial
genelist <- data.frame(ID = DE.Endo$gene_id, logFC = DE.Endo$avglog2FoldChange)
row.names(genelist)=genelist[,1]

Endo.enrich.go <- read.table(".\\enrich\\Endo.enrich.go.txt",sep = '\t',header = T)
library(GOplot)
Endo.go=data.frame(Category = Endo.enrich.go$ONTOLOGY,ID = Endo.enrich.go$ID,Term = Endo.enrich.go$Description, Genes = gsub("/", ", ", Endo.enrich.go$geneID), adj_pval = Endo.enrich.go$p.adjust)
#读取基因的logFC文件
circ <- circle_dat(Endo.go, genelist)
termNum = 10                                     #限定term数目
geneNum = nrow(genelist)                        #限定基因数目
chord <- chord_dat(circ, genelist[1:geneNum,], Endo.go$Term[1:termNum])
color_cols <- c(rgb(204,153,204,150,maxColorValue = 255),#CC99CC
                rgb(153,204,204,150,maxColorValue = 255), #99CCCC
                rgb(102,153,153,150,maxColorValue = 255),#669999
                rgb(102,153,204,150,maxColorValue = 255),#6699CC
                rgb(229,134,1,150,maxColorValue = 255),#E58601
                rgb(180,15,32,150,maxColorValue = 255),#B40F20
                rgb(108,211,152,150,maxColorValue = 255),#336699
                rgb(255,181,197,150,maxColorValue = 255),
                rgb(46,139,87,150,maxColorValue = 255),
                rgb(0,76,153,150,maxColorValue = 255))
pdf(file=".\\enrich\\Endo_circ.pdf",width = 30,height = 30)
GOChord(chord, 
        space = 0.001,           #基因之间的间距
        gene.order = 'logFC',    #按照logFC值对基因排序
        gene.space = 0.25,       #基因名跟圆圈的相对距离
        gene.size = 5,           #基因名字体大小 
        border.size = 0.1,       #线条粗细
        process.label = 8,
        ribbon.col=color_cols,
        lfc.col=c("#CC3333","#FFCCCC","#99CCFF"))       #term字体大小
dev.off()
termCol <- color_cols
####******Endo聚类圈图------------
pdf(file=".\\enrich\\Endo_cluster.pdf",width = 40,height = 10)
GOCluster(circ.gsym, 
          Endo.go$Term[1:termNum], 
          lfc.space = 0.5,                   #倍数跟树间的空隙大小
          lfc.width = 1,
          #clust.by = 'logFC',#变化倍数的圆圈宽度
          term.col = termCol[1:termNum],     #自定义term的颜色
          term.space = 0.2,                  #倍数跟term间的空隙大小
          term.width = 1)                    #富集term的圆圈宽度
dev.off()   


#
####******Astro圈图------------
diff_Astro <- read.table(".\\wilcox\\case-control-astro-wilcox.txt",sep = '\t',header = T)
diff_Astro <- diff_Astro[which(diff_Astro$pvalue<0.01),] #1082
Astro.gene <- diff_Astro[,1] #
#对于加载的注释库的使用，以上述为例，就直接在 OrgDb 中指定人（org.Hs.eg.db）或绵羊（sheep）
Astro.enrich.go <- enrichGO(gene = Astro.gene,  #基因列表文件中的基因名称
                            OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                            keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                            ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                            pAdjustMethod = 'fdr',  #指定 p 值校正方法
                            pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                            qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                            readable = FALSE)
Astro.enrich.go <-summary(Astro.enrich.go) #561
write.table(Astro.enrich.go,".\\enrich\\Astro.enrich.go.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
#astro.ego <- Astro.enrich.go[Astro.enrich.go$p.adjust<0.01,]  #23

DE.Astro <- diff_Astro
genelist <- data.frame(ID = DE.Astro$gene_id, logFC = DE.Astro$avglog2FoldChange)
row.names(genelist)=genelist[,1]

Astro.enrich.go <- read.table(".\\enrich\\Astro.enrich.go.txt",sep = '\t',header = T)
library(GOplot)
Astro.go=data.frame(Category = Astro.enrich.go$ONTOLOGY,ID = Astro.enrich.go$ID,Term = Astro.enrich.go$Description, Genes = gsub("/", ", ", Astro.enrich.go$geneID), adj_pval = Astro.enrich.go$p.adjust)
#读取基因的logFC文件
circ <- circle_dat(Astro.go, genelist)
termNum = 10                                     #限定term数目
geneNum = nrow(genelist)                        #限定基因数目
chord <- chord_dat(circ, genelist[1:geneNum,], Astro.go$Term[1:termNum])
color_cols <- c(rgb(204,153,204,150,maxColorValue = 255),#CC99CC
                rgb(153,204,204,150,maxColorValue = 255), #99CCCC
                rgb(102,153,153,150,maxColorValue = 255),#669999
                rgb(102,153,204,150,maxColorValue = 255),#6699CC
                rgb(229,134,1,150,maxColorValue = 255),#E58601
                rgb(180,15,32,150,maxColorValue = 255),#B40F20
                rgb(108,211,152,150,maxColorValue = 255),#336699
                rgb(255,181,197,150,maxColorValue = 255),
                rgb(46,139,87,150,maxColorValue = 255),
                rgb(0,76,153,150,maxColorValue = 255))
pdf(file=".\\enrich\\Astro_circ.pdf",width = 30,height = 30)
GOChord(chord, 
        space = 0.001,           #基因之间的间距
        gene.order = 'logFC',    #按照logFC值对基因排序
        gene.space = 0.25,       #基因名跟圆圈的相对距离
        gene.size = 3,           #基因名字体大小 
        border.size = 0.1,       #线条粗细
        process.label = 8,
        ribbon.col=color_cols,
        lfc.col=c("#CC3333","#FFCCCC","#99CCFF"))       #term字体大小
dev.off()
termCol <- color_cols
####******Astro聚类圈图------------
pdf(file=".\\enrich\\Astro_cluster.pdf",width = 40,height = 10)
GOCluster(circ.gsym, 
          Astro.go$Term[1:termNum], 
          lfc.space = 0.5,                   #倍数跟树间的空隙大小
          lfc.width = 1,
          #clust.by = 'logFC',#变化倍数的圆圈宽度
          term.col = termCol[1:termNum],     #自定义term的颜色
          term.space = 0.2,                  #倍数跟term间的空隙大小
          term.width = 1)                    #富集term的圆圈宽度
dev.off()   
####******Microglial圈图------------
diff_Microglial <- read.table(".\\wilcox\\case-control-Microglial-wilcox.txt",sep = '\t',header = T)
diff_Microglial <- diff_Microglial[which(diff_Microglial$pvalue<0.01),]
Microglial.gene <- diff_Microglial[,1] #309


Microglial.enrich.go <- enrichGO(gene = Microglial.gene,  #基因列表文件中的基因名称
                                 OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                                 keyType = 'SYMBOL',  
                                 ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                                 pAdjustMethod = 'fdr',  #指定 p 值校正方法
                                 pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                                 qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                                 readable = FALSE)
Microglial.enrich.go <-summary(Microglial.enrich.go)
write.table(Microglial.enrich.go,".\\enrich\\Microglial.enrich.go.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
#1

DE.Microglial <- diff_Microglial
genelist <- data.frame(ID = DE.Microglial$gene_id, logFC = DE.Microglial$avglog2FoldChange)
row.names(genelist)=genelist[,1]

Microglial.enrich.go <- read.table(".\\enrich\\Microglial.enrich.go.txt",sep = '\t',header = T)
library(GOplot)
Microglial.go=data.frame(Category = Microglial.enrich.go$ONTOLOGY,ID = Microglial.enrich.go$ID,Term = Microglial.enrich.go$Description, Genes = gsub("/", ", ", Microglial.enrich.go$geneID), adj_pval = Microglial.enrich.go$p.adjust)
#读取基因的logFC文件
circ <- circle_dat(Microglial.go, genelist)
termNum = 8                                     #限定term数目
geneNum = nrow(genelist)                        #限定基因数目
chord <- chord_dat(circ, genelist[1:geneNum,], Microglial.go$Term[1:termNum])
color_cols <- c(rgb(204,153,204,150,maxColorValue = 255),#CC99CC
                rgb(153,204,204,150,maxColorValue = 255), #99CCCC
                rgb(102,153,153,150,maxColorValue = 255),#669999
                rgb(102,153,204,150,maxColorValue = 255),#6699CC
                rgb(229,134,1,150,maxColorValue = 255),#E58601
                rgb(180,15,32,150,maxColorValue = 255),#B40F20
                rgb(108,211,152,150,maxColorValue = 255),#336699
                rgb(255,181,197,150,maxColorValue = 255),
                rgb(46,139,87,150,maxColorValue = 255),
                rgb(0,76,153,150,maxColorValue = 255))
pdf(file=".\\enrich\\Microglial_circ.pdf",width = 30,height = 30)
GOChord(chord, 
        space = 0.001,           #基因之间的间距
        gene.order = 'logFC',    #按照logFC值对基因排序
        gene.space = 0.25,       #基因名跟圆圈的相对距离
        gene.size = 5,           #基因名字体大小 
        border.size = 0.1,       #线条粗细
        process.label = 8,
        ribbon.col=color_cols[1:termNum],
        lfc.col=c("#CC3333","#FFCCCC","#99CCFF"))       #term字体大小
dev.off()
termCol <- color_cols
####******Microglial聚类圈图------------
pdf(file=".\\enrich\\Microglial_cluster.pdf",width = 40,height = 10)
GOCluster(circ.gsym, 
          Microglial.go$Term[1:termNum], 
          lfc.space = 0.5,                   #倍数跟树间的空隙大小
          lfc.width = 1,
          #clust.by = 'logFC',#变化倍数的圆圈宽度
          term.col = termCol[1:termNum],     #自定义term的颜色
          term.space = 0.2,                  #倍数跟term间的空隙大小
          term.width = 1)                    #富集term的圆圈宽度
dev.off()   

####******Oligodendrocyte_precursor_cell圈图------------
diff_OPC <- read.table(".\\wilcox\\case-control-Oligodendrocyte_precursor_cell-wilcox.txt",sep = '\t',header = T)
diff_OPC <- diff_OPC[which(diff_OPC$pvalue<0.01),]  #213
OPC.gene <- diff_OPC[,1]  

#213
#对于加载的注释库的使用，以上述为例，就直接在 OrgDb 中指定人（org.Hs.eg.db）或绵羊（sheep）
OPC.enrich.go <- enrichGO(gene = OPC.gene,  #基因列表文件中的基因名称
                          OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                          keyType = 'SYMBOL',  
                          ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                          pAdjustMethod = 'fdr',  #指定 p 值校正方法
                          pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                          qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                          readable = FALSE)
OPC.enrich.go <-summary(OPC.enrich.go)
write.table(OPC.enrich.go,".\\enrich\\OPC.enrich.go.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
#342

DE.OPC <- diff_OPC
genelist <- data.frame(ID = DE.OPC$gene_id, logFC = DE.OPC$avglog2FoldChange)
row.names(genelist)=genelist[,1]

OPC.enrich.go <- read.table(".\\enrich\\OPC.enrich.go.txt",sep = '\t',header = T)
library(GOplot)
OPC.go=data.frame(Category = OPC.enrich.go$ONTOLOGY,ID = OPC.enrich.go$ID,Term = OPC.enrich.go$Description, Genes = gsub("/", ", ", OPC.enrich.go$geneID), adj_pval = OPC.enrich.go$p.adjust)
#读取基因的logFC文件
circ <- circle_dat(OPC.go, genelist)
termNum = 10                                     #限定term数目
geneNum = nrow(genelist)                        #限定基因数目
chord <- chord_dat(circ, genelist[1:geneNum,], OPC.go$Term[1:termNum])
color_cols <- c(rgb(204,153,204,150,maxColorValue = 255),#CC99CC
                rgb(153,204,204,150,maxColorValue = 255), #99CCCC
                rgb(102,153,153,150,maxColorValue = 255),#669999
                rgb(102,153,204,150,maxColorValue = 255),#6699CC
                rgb(229,134,1,150,maxColorValue = 255),#E58601
                rgb(180,15,32,150,maxColorValue = 255),#B40F20
                rgb(108,211,152,150,maxColorValue = 255),#336699
                rgb(255,181,197,150,maxColorValue = 255),
                rgb(46,139,87,150,maxColorValue = 255),
                rgb(0,76,153,150,maxColorValue = 255))
pdf(file=".\\enrich\\OPC_circ.pdf",width = 30,height = 30)
GOChord(chord, 
        space = 0.001,           #基因之间的间距
        gene.order = 'logFC',    #按照logFC值对基因排序
        gene.space = 0.25,       #基因名跟圆圈的相对距离
        gene.size = 5,           #基因名字体大小 
        border.size = 0.1,       #线条粗细
        process.label = 8,
        ribbon.col=color_cols,
        lfc.col=c("#CC3333","#FFCCCC","#99CCFF"))       #term字体大小
dev.off()
termCol <- color_cols
####******OPC聚类圈图------------
pdf(file=".\\enrich\\OPC_cluster.pdf",width = 40,height = 10)
GOCluster(circ.gsym, 
          OPC.go$Term[1:termNum], 
          lfc.space = 0.5,                   #倍数跟树间的空隙大小
          lfc.width = 1,
          #clust.by = 'logFC',#变化倍数的圆圈宽度
          term.col = termCol[1:termNum],     #自定义term的颜色
          term.space = 0.2,                  #倍数跟term间的空隙大小
          term.width = 1)                    #富集term的圆圈宽度
dev.off()   
####******Oligodendrocyte圈图------------
diff_Oligodendrocyte <- read.table(".\\wilcox\\case-control-Oligodendrocyte-wilcox.txt",sep = '\t',header = T)
diff_Oligodendrocyte <- diff_Oligodendrocyte[which(diff_Oligodendrocyte$pvalue<0.01),]
Oligodendrocyte.gene <- diff_Oligodendrocyte[,1]  #1211

#对于加载的注释库的使用，以上述为例，就直接在 OrgDb 中指定人（org.Hs.eg.db）或绵羊（sheep）
Oligodendrocyte.enrich.go <- enrichGO(gene = Oligodendrocyte.gene,  #基因列表文件中的基因名称
                                      OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                                      keyType = 'SYMBOL',  
                                      ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                                      pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                                      qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                                      readable = FALSE)
Oligodendrocyte.enrich.go <-summary(Oligodendrocyte.enrich.go) #823
write.table(Oligodendrocyte.enrich.go,".\\enrich\\Oligodendrocyte.enrich.go.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')


DE.Oligodendrocyte <- diff_Oligodendrocyte
genelist <- data.frame(ID = DE.Oligodendrocyte$gene_id, logFC = DE.Oligodendrocyte$avglog2FoldChange)
row.names(genelist)=genelist[,1]

Oligodendrocyte.enrich.go <- read.table(".\\enrich\\Oligodendrocyte.enrich.go.txt",sep = '\t',header = T)
library(GOplot)
Oligodendrocyte.go=data.frame(Category = Oligodendrocyte.enrich.go$ONTOLOGY,ID = Oligodendrocyte.enrich.go$ID,Term = Oligodendrocyte.enrich.go$Description, Genes = gsub("/", ", ", Oligodendrocyte.enrich.go$geneID), adj_pval = Oligodendrocyte.enrich.go$p.adjust)
#读取基因的logFC文件
circ <- circle_dat(Oligodendrocyte.go, genelist)
termNum = 10                                     #限定term数目
geneNum = nrow(genelist)                        #限定基因数目
chord <- chord_dat(circ, genelist[1:geneNum,], Oligodendrocyte.go$Term[1:termNum])
color_cols <- c(rgb(204,153,204,150,maxColorValue = 255),#CC99CC
                rgb(153,204,204,150,maxColorValue = 255), #99CCCC
                rgb(102,153,153,150,maxColorValue = 255),#669999
                rgb(102,153,204,150,maxColorValue = 255),#6699CC
                rgb(229,134,1,150,maxColorValue = 255),#E58601
                rgb(180,15,32,150,maxColorValue = 255),#B40F20
                rgb(108,211,152,150,maxColorValue = 255),#336699
                rgb(255,181,197,150,maxColorValue = 255),
                rgb(46,139,87,150,maxColorValue = 255),
                rgb(0,76,153,150,maxColorValue = 255))
pdf(file=".\\enrich\\Oligodendrocyte_circ.pdf",width = 30,height = 30)
GOChord(chord, 
        space = 0.001,           #基因之间的间距
        gene.order = 'logFC',    #按照logFC值对基因排序
        gene.space = 0.25,       #基因名跟圆圈的相对距离
        gene.size = 3,           #基因名字体大小 
        border.size = 0.1,       #线条粗细
        process.label = 8,
        ribbon.col=color_cols,
        lfc.col=c("#CC3333","#FFCCCC","#99CCFF"))       #term字体大小
dev.off()
termCol <- color_cols
####******Oligodendrocyte聚类圈图------------
pdf(file=".\\enrich\\Oligodendrocyte_cluster.pdf",width = 40,height = 10)
GOCluster(circ.gsym, 
          Oligodendrocyte.go$Term[1:termNum], 
          lfc.space = 0.5,                   #倍数跟树间的空隙大小
          lfc.width = 1,
          #clust.by = 'logFC',#变化倍数的圆圈宽度
          term.col = termCol[1:termNum],     #自定义term的颜色
          term.space = 0.2,                  #倍数跟term间的空隙大小
          term.width = 1)                    #富集term的圆圈宽度
dev.off()
####******Excitatory neuron圈图--------
diff_Excitatory_neuron <- read.table(".\\wilcox\\case-control-Excitatory_neuron-wilcox.txt",sep = '\t',header = T)
diff_Excitatory_neuron <- diff_Excitatory_neuron[which(diff_Excitatory_neuron$padj<0.01),]
Excitatory_neuron.gene <- diff_Excitatory_neuron[,1]
#1693
#对于加载的注释库的使用，以上述为例，就直接在 OrgDb 中指定人（org.Hs.eg.db）或绵羊（sheep）
Excitatory_neuron.enrich.go <- enrichGO(gene = Excitatory_neuron.gene,  #基因列表文件中的基因名称
                                        OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                                        keyType = 'SYMBOL',  
                                        ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                                        pAdjustMethod = 'fdr',  #指定 p 值校正方法
                                        pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                                        qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                                        readable = FALSE)
Excitatory_neuron.enrich.go <-summary(Excitatory_neuron.enrich.go)
write.table(Excitatory_neuron.enrich.go,".\\enrich\\Excitatory_neuron.enrich.go.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
#727

DE.Excitatory_neuron <- diff_Excitatory_neuron
genelist <- data.frame(ID = DE.Excitatory_neuron$gene_id, logFC = DE.Excitatory_neuron$avglog2FoldChange)
row.names(genelist)=genelist[,1]

Excitatory_neuron.enrich.go <- read.table(".\\enrich\\Excitatory_neuron.enrich.go.txt",sep = '\t',header = T)
library(GOplot)
Excitatory_neuron.go=data.frame(Category = Excitatory_neuron.enrich.go$ONTOLOGY,ID = Excitatory_neuron.enrich.go$ID,Term = Excitatory_neuron.enrich.go$Description, Genes = gsub("/", ", ", Excitatory_neuron.enrich.go$geneID), adj_pval = Excitatory_neuron.enrich.go$p.adjust)
#读取基因的logFC文件
circ <- circle_dat(Excitatory_neuron.go, genelist)
termNum = 10                                     #限定term数目
geneNum = nrow(genelist)                        #限定基因数目
chord <- chord_dat(circ, genelist[1:geneNum,], Excitatory_neuron.go$Term[1:termNum])
color_cols <- c(rgb(204,153,204,150,maxColorValue = 255),#CC99CC
                rgb(153,204,204,150,maxColorValue = 255), #99CCCC
                rgb(102,153,153,150,maxColorValue = 255),#669999
                rgb(102,153,204,150,maxColorValue = 255),#6699CC
                rgb(229,134,1,150,maxColorValue = 255),#E58601
                rgb(180,15,32,150,maxColorValue = 255),#B40F20
                rgb(108,211,152,150,maxColorValue = 255),#336699
                rgb(255,181,197,150,maxColorValue = 255),
                rgb(46,139,87,150,maxColorValue = 255),
                rgb(0,76,153,150,maxColorValue = 255))
pdf(file=".\\enrich\\Excitatory_neuron_circ.pdf",width = 30,height = 30)
GOChord(chord, 
        space = 0.001,           #基因之间的间距
        gene.order = 'logFC',    #按照logFC值对基因排序
        gene.space = 0.25,       #基因名跟圆圈的相对距离
        gene.size = 2,           #基因名字体大小 
        border.size = 0.1,       #线条粗细
        process.label = 8,
        ribbon.col=color_cols,
        lfc.col=c("#CC3333","#FFCCCC","#99CCFF"))       #term字体大小
dev.off()
termCol <- color_cols
####******Excitatory_neuron聚类圈图------------
pdf(file=".\\enrich\\Excitatory_neuron_cluster.pdf",width = 40,height = 10)
GOCluster(circ.gsym, 
          Excitatory_neuron.go$Term[1:termNum], 
          lfc.space = 0.5,                   #倍数跟树间的空隙大小
          lfc.width = 1,
          #clust.by = 'logFC',#变化倍数的圆圈宽度
          term.col = termCol[1:termNum],     #自定义term的颜色
          term.space = 0.2,                  #倍数跟term间的空隙大小
          term.width = 1)                    #富集term的圆圈宽度
dev.off()

####******Inhibitory_neuron圈图---------
diff_Inhibitory_neuron <- read.table(".\\wilcox\\case-control-Inhibitory_neuron-wilcox.txt",sep = '\t',header = T)
diff_Inhibitory_neuron <- diff_Inhibitory_neuron[which(diff_Inhibitory_neuron$padj<0.01),]
Inhibitory_neuron.gene <- diff_Inhibitory_neuron[which(diff_Inhibitory_neuron$padj<0.01),1]

#1475

#对于加载的注释库的使用，以上述为例，就直接在 OrgDb 中指定人（org.Hs.eg.db）或绵羊（sheep）
Inhibitory_neuron.enrich.go <- enrichGO(gene = Inhibitory_neuron.gene,  #基因列表文件中的基因名称
                                        OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                                        keyType = 'SYMBOL',  
                                        ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                                        pAdjustMethod = 'fdr',  #指定 p 值校正方法
                                        pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                                        qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                                        readable = FALSE)
Inhibitory_neuron.enrich.go <-summary(Inhibitory_neuron.enrich.go)
# aa <- grep("glu",Inhibitory_neuron.enrich.go$Description)
# bb <- Inhibitory_neuron.enrich.go[aa,]
write.table(Inhibitory_neuron.enrich.go,".\\enrich\\Inhibitory_neuron.enrich.go.txt",sep = '\t',col.names = T,row.names = F)
#726

DE.Inhibitory_neuron <- diff_Inhibitory_neuron
genelist <- data.frame(ID = DE.Inhibitory_neuron$gene_id, logFC = DE.Inhibitory_neuron$avglog2FoldChange)
row.names(genelist)=genelist[,1]

Inhbitory_neuron.enrich.go <- read.table(".\\enrich\\Inhibitory_neuron.enrich.go.txt",sep = '\t',header = T)

library(GOplot)
Inhibitory_neuron.go=data.frame(Category = Inhbitory_neuron.enrich.go$ONTOLOGY,ID = Inhbitory_neuron.enrich.go$ID,Term = Inhbitory_neuron.enrich.go$Description, Genes = gsub("/", ", ", Inhbitory_neuron.enrich.go$geneID), adj_pval = Inhbitory_neuron.enrich.go$p.adjust)
#读取基因的logFC文件
circ <- circle_dat(Inhibitory_neuron.go, genelist)
termNum = 10                                     #限定term数目
geneNum = nrow(genelist)                        #限定基因数目
chord <- chord_dat(circ, genelist[1:geneNum,], Inhibitory_neuron.go$Term[1:termNum])
color_cols <- c(rgb(204,153,204,150,maxColorValue = 255),#CC99CC
                rgb(153,204,204,150,maxColorValue = 255), #99CCCC
                rgb(102,153,153,150,maxColorValue = 255),#669999
                rgb(102,153,204,150,maxColorValue = 255),#6699CC
                rgb(229,134,1,150,maxColorValue = 255),#E58601
                rgb(180,15,32,150,maxColorValue = 255),#B40F20
                rgb(108,211,152,150,maxColorValue = 255),#336699
                rgb(255,181,197,150,maxColorValue = 255),
                rgb(46,139,87,150,maxColorValue = 255),
                rgb(0,76,153,150,maxColorValue = 255))
pdf(file=".\\enrich\\Inhibitory_neuron_circ.pdf",width = 30,height = 30)
GOChord(chord, 
        space = 0.001,           #基因之间的间距
        gene.order = 'logFC',    #按照logFC值对基因排序
        gene.space = 0.25,       #基因名跟圆圈的相对距离
        gene.size = 3,           #基因名字体大小 
        border.size = 0.1,       #线条粗细
        process.label = 8,
        ribbon.col=color_cols,
        lfc.col=c("#CC3333","#FFCCCC","#99CCFF"))       #term字体大小
dev.off()
termCol <- color_cols

####******Inhibitory_neuron聚类圈图------------
pdf(file=".\\enrich\\Inhibitory_neuron_cluster.pdf",width = 40,height = 10)
GOCluster(circ.gsym, 
          Inhibitory_neuron.go$Term[1:termNum], 
          lfc.space = 0.5,                   #倍数跟树间的空隙大小
          lfc.width = 1,
          #clust.by = 'logFC',#变化倍数的圆圈宽度
          term.col = termCol[1:termNum],     #自定义term的颜色
          term.space = 0.2,                  #倍数跟term间的空隙大小
          term.width = 1)                    #富集term的圆圈宽度
dev.off()


#####
#####
#####
####******富集热图+点图----------
BiocManager::install("BoutrosLab.plotting.general")
library(BoutrosLab.plotting.general)
library(ggplot2)
####癫痫与正常细胞间差异基因--------------
##根据GO term 的富集程度即P值画热图，行为GO term，列为细胞类型，值为logP
##ENDo
Endo.enrich.go <- read.table(".\\enrich\\Endo.enrich.go.txt",sep = '\t',header = T)
Endo.BP <- Endo.enrich.go[Endo.enrich.go$ONTOLOGY=="BP",]
Endo.BP <- Endo.BP[order(Endo.BP$pvalue),]
####
Astro.enrich.go <- read.table(".\\enrich\\Astro.enrich.go.txt",sep = '\t',header = T)
Astro.BP <- Astro.enrich.go[Astro.enrich.go$ONTOLOGY=="BP",]
Astro.BP <- Astro.BP[order(Astro.BP$pvalue),]

####
Microglial.enrich.go <- read.table(".\\enrich\\Microglial.enrich.go.txt",sep = '\t',header = T)
Microglial.BP <- Microglial.enrich.go[Microglial.enrich.go$ONTOLOGY=="BP",]
Microglial.BP <- Microglial.BP[order(Microglial.BP$pvalue),]

####
OPC.enrich.go <- read.table(".\\enrich\\OPC.enrich.go.txt",sep = '\t',header = T)
OPC.BP <- OPC.enrich.go[OPC.enrich.go$ONTOLOGY=="BP",]
OPC.BP <- OPC.BP[order(OPC.BP$pvalue),]

####
Oligodendrocyte.enrich.go <- read.table(".\\enrich\\Oligodendrocyte.enrich.go.txt",sep = '\t',header = T)
Oligodendrocyte.BP <- Oligodendrocyte.enrich.go[Oligodendrocyte.enrich.go$ONTOLOGY=="BP",]
Oligodendrocyte.BP <- Oligodendrocyte.BP[order(Oligodendrocyte.BP$pvalue),]

####
Excitatory_neuron.enrich.go <- read.table(".\\enrich\\Excitatory_neuron.enrich.go.txt",sep = '\t',header = T)
Excitatory.BP <- Excitatory_neuron.enrich.go[Excitatory_neuron.enrich.go$ONTOLOGY=="BP",]
Excitatory.BP <- Excitatory.BP[order(Excitatory.BP$pvalue),]

####
Inhbitory_neuron.enrich.go <- read.table(".\\enrich\\Inhibitory_neuron.enrich.go.txt",sep = '\t',header = T)
Inhbitory.BP <- Inhbitory_neuron.enrich.go[Inhbitory_neuron.enrich.go$ONTOLOGY=="BP",]
Inhbitory.BP <- Inhbitory.BP[order(Inhbitory.BP$pvalue),]
####******构建矩阵--------
###选取前50BP画图。细胞类型之间富集的功能肯定不一致，若没有重合，值为1
Astro.BP <- Astro.BP[1:10,c(3,6,10)]
rownames(Astro.BP) <- Astro.BP$Description
colnames(Astro.BP) <- c("Description","Astro_p","Astro_count")

Oligodendrocyte.BP <- Oligodendrocyte.BP[1:10,c(3,6,10)]
rownames(Oligodendrocyte.BP) <- Oligodendrocyte.BP$Description
colnames(Oligodendrocyte.BP) <- c("Description","Oli_p","Oli_count")


OPC.BP <- OPC.BP[1:10,c(3,6,10)]
rownames(OPC.BP) <- OPC.BP$Description
colnames(OPC.BP) <- c("Description","OPC_p","OPC_count")


Endo.BP <- Endo.BP[,c(3,6,10)]
rownames(Endo.BP) <- Endo.BP$Description
colnames(Endo.BP) <- c("Description","Endo_p","Endo_count")

# Microglial.BP <- Microglial.enrich.go[,c(3,6)]
# rownames(Microglial.BP) <- Microglial.BP$Description

Excitatory.BP <- Excitatory.BP[1:10,c(3,6,10)]
rownames(Excitatory.BP) <- Excitatory.BP$Description
colnames(Excitatory.BP) <- c("Description","Excitatory_p","Excitatory_count")


Inhbitory.BP <- Inhbitory.BP[1:10,c(3,6,10)]
rownames(Inhbitory.BP) <- Inhbitory.BP$Description
colnames(Inhbitory.BP) <- c("Description","Inhbitory_p","Inhbitory_count")


EX_IN <- merge(Excitatory.BP, Inhbitory.BP, by='Description', all=TRUE)
EX_IN_Oli <- merge(EX_IN, Oligodendrocyte.BP, by='Description', all=TRUE)
EX_IN_Oli_Astro <- merge(EX_IN_Oli, Astro.BP, by='Description', all=TRUE)
EX_IN_Oli_Astro_Opcs <- merge(EX_IN_Oli_Astro, OPC.BP, by='Description', all=TRUE)
EX_IN_Oli_Astro_Opcs_Endo <- merge(EX_IN_Oli_Astro_Opcs, Endo.BP, by='Description', all=TRUE)
dim(EX_IN_Oli_Astro_Opcs_Endo)
#24 13
exp <- EX_IN_Oli_Astro_Opcs_Endo[,c(2:ncol(EX_IN_Oli_Astro_Opcs_Endo))]
rownames(exp) <- EX_IN_Oli_Astro_Opcs_Endo$Description

spot.size.function<-function(x){
  x= (-log2(x))/15
}

spot.colour.function <- function(x) {
  colours <- rep("white", length(x));
  colours[x==exp$Excitatory_p] <- "#CC99CC"; 
  colours[x==exp$Inhbitory_p] <- "#99CCCC"; 
  colours[x==exp$Oli_p] <- "#669999";
  colours[x==exp$Astro_p] <- "#6699CC";
  colours[x==exp$OPC_p] <- "#B40F20";
  colours[x==exp$Endo_p] <- "#336699";
  return(colours);
}
color_cols <- c(rgb(204,153,204,150,maxColorValue = 255),#CC99CC
                rgb(153,204,204,150,maxColorValue = 255), #99CCCC
                rgb(102,153,153,150,maxColorValue = 255),#669999
                rgb(102,153,204,150,maxColorValue = 255),#6699CC
                rgb(180,15,32,150,maxColorValue = 255),#B40F20
                rgb(108,211,152,150,maxColorValue = 255))#336699
#加图例
key.trans <- list(title = "cell type",space = "right",columns = 2,
                  points = list(pch = c(20,20,20,20,20,20),
                                col = c(color_cols,color_cols),
                                cex=c(1,1,1,1,1,1,1,2,3,4,5,6)),
                  text = list(c("Excitatory neuron","Inhibitory neuron","Oligodendrocyte","Astrocyte",
                                "OPCs","Endothelial",
                                "Excitatory neuron","Inhibitory neuron","Oligodendrocyte","Astrocyte",
                                "OPCs","Endothelial")),
                  #lines = list(col = colors,lty = lines),
                  cex.title = 1,cex = .9)

create.dotmap(
  exp[,c(1,3,5,7,9,11)],bg.data = exp[,c(2,4,6,8,10,12)],
  pch = 20,na.spot.size=2,spot.size.function=spot.size.function,
  spot.colour.function=spot.colour.function,
  ylab.cex=1,yaxis.cex=0.5,xaxis.cex=0.5,
  colour.scheme=colorRampPalette(c('#FFFFFF','#6699CC'))(100),
  key = key.trans,
  colourkey = T
)



######******细胞类型间差异基因-----------------
Endo.enrich.go <- read.table(".\\enrich\\celltype_Endothelial.go.txt",sep = '\t',header = T)
Endo.BP <- Endo.enrich.go[Endo.enrich.go$ONTOLOGY=="BP",]
Endo.BP <- Endo.BP[order(Endo.BP$pvalue),]
####
Astro.enrich.go <- read.table(".\\enrich\\celltype_Astrocyte.go.txt",sep = '\t',header = T)
Astro.BP <- Astro.enrich.go[Astro.enrich.go$ONTOLOGY=="BP",]
Astro.BP <- Astro.BP[order(Astro.BP$pvalue),]

####
Microglial.enrich.go <- read.table(".\\enrich\\celltype_Microglial.go.txt",sep = '\t',header = T)
Microglial.BP <- Microglial.enrich.go[Microglial.enrich.go$ONTOLOGY=="BP",]
Microglial.BP <- Microglial.BP[order(Microglial.BP$pvalue),]

####
OPC.enrich.go <- read.table(".\\enrich\\celltype_OPC.go.txt",sep = '\t',header = T)
OPC.BP <- OPC.enrich.go[OPC.enrich.go$ONTOLOGY=="BP",]
OPC.BP <- OPC.BP[order(OPC.BP$pvalue),]

####
Oligodendrocyte.enrich.go <- read.table(".\\enrich\\celltype_Oligodendrocyte.go.txt",sep = '\t',header = T)
Oligodendrocyte.BP <- Oligodendrocyte.enrich.go[Oligodendrocyte.enrich.go$ONTOLOGY=="BP",]
Oligodendrocyte.BP <- Oligodendrocyte.BP[order(Oligodendrocyte.BP$pvalue),]

####
Excitatory_neuron.enrich.go <- read.table(".\\enrich\\celltype_Excitatory.go.txt",sep = '\t',header = T)
Excitatory.BP <- Excitatory_neuron.enrich.go[Excitatory_neuron.enrich.go$ONTOLOGY=="BP",]
Excitatory.BP <- Excitatory.BP[order(Excitatory.BP$pvalue),]

####
Inhbitory_neuron.enrich.go <- read.table(".\\enrich\\celltype_Inhibitory.go.txt",sep = '\t',header = T)
Inhbitory.BP <- Inhbitory_neuron.enrich.go[Inhbitory_neuron.enrich.go$ONTOLOGY=="BP",]
Inhbitory.BP <- Inhbitory.BP[order(Inhbitory.BP$pvalue),]

####******构建矩阵--------
###选取前50BP画图。细胞类型之间富集的功能肯定不一致，若没有重合，值为1
Astro.BP <- Astro.BP[1:10,c(3,6,10)]
rownames(Astro.BP) <- Astro.BP$Description
colnames(Astro.BP) <- c("Description","Astro_p","Astro_count")

Oligodendrocyte.BP <- Oligodendrocyte.BP[1:10,c(3,6,10)]
rownames(Oligodendrocyte.BP) <- Oligodendrocyte.BP$Description
colnames(Oligodendrocyte.BP) <- c("Description","Oli_p","Oli_count")


OPC.BP <- OPC.BP[1:10,c(3,6,10)]
rownames(OPC.BP) <- OPC.BP$Description
colnames(OPC.BP) <- c("Description","OPC_p","OPC_count")


Endo.BP <- Endo.BP[1:10,c(3,6,10)]
rownames(Endo.BP) <- Endo.BP$Description
colnames(Endo.BP) <- c("Description","Endo_p","Endo_count")


Microglial.BP <- Microglial.BP[1:10,c(3,6,10)]
rownames(Microglial.BP) <- Microglial.BP$Description
colnames(Microglial.BP) <- c("Description","Microglial_p","Microglial_count")


Excitatory.BP <- Excitatory.BP[1:10,c(3,6,10)]
rownames(Excitatory.BP) <- Excitatory.BP$Description
colnames(Excitatory.BP) <- c("Description","Excitatory_p","Excitatory_count")


Inhbitory.BP <- Inhbitory.BP[1:10,c(3,6,10)]
rownames(Inhbitory.BP) <- Inhbitory.BP$Description
colnames(Inhbitory.BP) <- c("Description","Inhbitory_p","Inhbitory_count")






EX_IN <- merge(Excitatory.BP, Inhbitory.BP, by='Description', all=TRUE)
EX_IN_Oli <- merge(EX_IN, Oligodendrocyte.BP, by='Description', all=TRUE)
EX_IN_Oli_Astro <- merge(EX_IN_Oli, Astro.BP, by='Description', all=TRUE)
EX_IN_Oli_Astro_Micro <- merge(EX_IN_Oli_Astro, Microglial.BP, by='Description', all=TRUE)
EX_IN_Oli_Astro_Micro_Opcs <- merge(EX_IN_Oli_Astro_Micro, OPC.BP, by='Description', all=TRUE)
EX_IN_Oli_Astro_Micro_Opcs_Endo <- merge(EX_IN_Oli_Astro_Micro_Opcs, Endo.BP, by='Description', all=TRUE)
dim(EX_IN_Oli_Astro_Micro_Opcs_Endo)
#46 15
exp <- EX_IN_Oli_Astro_Micro_Opcs_Endo[,c(2:ncol(EX_IN_Oli_Astro_Micro_Opcs_Endo))]
rownames(exp) <- EX_IN_Oli_Astro_Micro_Opcs_Endo$Description

spot.size.function<-function(x){
  x= (-log2(x))/30
}

spot.colour.function <- function(x) {
  colours <- rep("white", length(x));
  colours[x==exp$Excitatory_p] <- "#CC99CC"; 
  colours[x==exp$Inhbitory_p] <- "#99CCCC"; 
  colours[x==exp$Oli_p] <- "#669999";
  colours[x==exp$Astro_p] <- "#6699CC";
  colours[x==exp$Microglial_p] <- "#E58601";
  colours[x==exp$OPC_p] <- "#B40F20";
  colours[x==exp$Endo_p] <- "#336699";
  return(colours);
}
color_cols <- c(rgb(204,153,204,150,maxColorValue = 255),#CC99CC
                rgb(153,204,204,150,maxColorValue = 255), #99CCCC
                rgb(102,153,153,150,maxColorValue = 255),#669999
                rgb(102,153,204,150,maxColorValue = 255),#6699CC
                rgb(229,134,1,150,maxColorValue = 255),#E58601
                rgb(180,15,32,150,maxColorValue = 255),#B40F20
                rgb(108,211,152,150,maxColorValue = 255))#336699
#加图例
key.trans <- list(title = "cell type",space = "right",columns = 1,
                  points = list(pch = c(20,20,20,20,20,20),
                                col = color_cols,
                                cex=c(1,1,1,1)),
                  text = list(c("Excitatory neuron","Inhibitory neuron","Oligodendrocyte","Astrocyte",
                                "Microglial","OPCs","Endothelial")),
                  
                  #lines = list(col = colors,lty = lines),
                  cex.title = 1,cex = .9)
create.dotmap(
  exp[,c(1,3,5,7,9,11,13)],bg.data = exp[,c(2,4,6,8,10,12,14)],
  pch = 20,na.spot.size=1,spot.size.function=spot.size.function,
  spot.colour.function=spot.colour.function,
  ylab.cex=1,yaxis.cex=0.5,xaxis.cex=0.5,
  colour.scheme=colorRampPalette(c('#FFFFFF','#6699CC'))(100),
  key = key.trans
)