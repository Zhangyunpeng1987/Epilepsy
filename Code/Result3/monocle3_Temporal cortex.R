###*****monocle3--------
Excitatory <- readRDS(file = "./temporal_lobe\\subcluster\\Excitatory\\resolution8\\Excitatory_anno_cluster_1.rds")
library(Seurat)
library(monocle3)
data <- GetAssayData(Excitatory, assay = "RNA", slot = "counts")  
cell_metadata <- Excitatory@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)    
plot_pc_variance_explained(cds)   
cds <- reduce_dimension(cds,preprocess_method = "PCA")
cds <- reduce_dimension(cds, reduction_method="UMAP")
plot_cells(cds)
cds <- cluster_cells(cds, resolution=1e-5)
plot_cells(cds)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition")
plot_cells(cds, color_cells_by="group")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(Excitatory, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed   

cds <- learn_graph(cds,use_partition = FALSE) 
plot_cells(cds,
           color_cells_by = "Anno_Idents",
           label_groups_by_cluster=F,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           group_label_size=4,
           cell_size=1.5)
get_earliest_principal_node  <- function(cds, time_bin="Nor"){
  cell_ids <- which(colData(cds)[, "group"] == time_bin)
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))    #很多人这里一直在哭诉error，那是都没有理解这一步在干嘛，很多无脑运行别人的代码，别人代码选择XX细胞作为根，你的数据集里又没有，所以报错说没有node啦
library(ggplot2)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1)
plot_cells(cds,
           reduction_method = "UMAP",
           #trajectory_graph_color="blue",
           color_cells_by = "group",
           label_groups_by_cluster=F,
           label_cell_groups = TRUE,
           label_leaves=F,
           label_branch_points=TRUE,
           group_label_size=6,
           cell_size=1)+scale_colour_manual(values = c("#FFCCCC","#6699CC"))                                   #这个图就可以看出细胞直接的分化轨迹了
plot_cells(cds,
           reduction_method = "UMAP",
           #trajectory_graph_color="blue",
           color_cells_by = "Anno_Idents",
           label_groups_by_cluster=T,
           label_cell_groups = TRUE,
           label_leaves=T,
           label_branch_points=TRUE,
           group_label_size=5,
           cell_size=1)+scale_colour_manual(values = cols) 
saveRDS(cds,file="./temporal_lobe\\subcluster\\Excitatory\\monocle3\\Excitatory_reault-cds_1.rds")

####******差异分析------------
cds <- readRDS(file="./temporal_lobe\\subcluster\\Excitatory\\monocle3\\Nor\\Ex_Nor_1-cds.rds")
pr_graph_test_res = graph_test(cds, neighbor_graph="principal_graph", cores=8)
write.table(pr_graph_test_res,file="./temporal_lobe\\subcluster\\Excitatory\\monocle3\\Nor\\diff\\pr_graph_test_res.txt",col.names=T,row.names=T,sep="\t",quote=F)
pr_graph_test_res <- read.table(file="./temporal_lobe\\subcluster\\Excitatory\\monocle3\\Nor\\diff\\pr_graph_test_res.txt",header = T,sep="\t")
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
length(pr_deg_ids)  
###********关键基因展示----------
Track_genes_sig <- pr_graph_test_res %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="Anno_Idents", 
                         min_expr=0.5, ncol = 2)+
  scale_colour_manual(values = cols)

###*******将monocle3差异基因提取出来，构建monocle2差异基因伪时间热图-------
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
cds_NOR <- readRDS(file="./temporal_lobe\\subcluster\\Excitatory\\monocle3\\Nor\\Ex_Nor_1-cds.rds")
modulated_genes_NOR <- read.table(file="./temporal_lobe\\subcluster\\Excitatory\\monocle3\\Nor\\diff\\pr_graph_test_res.txt",header = T,sep="\t")

genes <- row.names(subset(modulated_genes_NOR, q_value < 0.05 & morans_I > 0.5))
genes 
pt.matrix <- exprs(cds_NOR)[match(genes,rownames(rowData(cds_NOR))),order(pseudotime(cds_NOR))]

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,km=4)
print(hthc)
#dev.off()
###*******伪时间热图不同簇基因富集分析-------
dhm <- draw(hthc)
dhm1 <-as.data.frame(dhm)
aa <- unlist(row_order(dhm)) %>% as.data.frame()
###*******簇1-------
C1 <- rownames(pt.matrix)[row_order(dhm)[[1]]]
library(clusterProfiler)
library(org.Hs.eg.db)
C1.go <- enrichGO(gene = C1,  #
                  OrgDb = 'org.Hs.eg.db',  
                  keyType = 'SYMBOL', 
                  ont = 'ALL',  
                  pAdjustMethod = 'fdr',  
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.2,  
                  readable = FALSE)  
C1.go <-summary(C1.go) 
write.table(C1.go,"./temporal_lobe\\subcluster\\Excitatory\\monocle3\\Nor\\diff\\C1.go.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
C1.go_top=C1.go %>% group_by(ONTOLOGY) %>% top_n(n=(-5),wt=pvalue)

C1.go_top <- C1.go_top[c(1:5),]
C1.go_top$Description<- factor(C1.go_top$Description,levels = C1.go_top$Description)
C1.go_top$Group <- rep("C1",dim(C1.go_top)[1])

library(ggplot2)
library(RColorBrewer)
color <- brewer.pal(3,"Dark2")
colorl <- rep(color,each=10)
ggplot(C1.go_top, aes(p.adjust, Description)) +
  geom_point(aes(color=p.adjust, size=Count))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text.y = element_text(colour = "black"))+
  scale_color_gradient(low='#6699CC',high='#CC3333')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
###*******簇2-------
C2 <- rownames(pt.matrix)[row_order(dhm)[[2]]]
library(clusterProfiler)
library(org.Hs.eg.db)
C2.go <- enrichGO(gene = C2,  #基因列表文件中的基因名称
                  OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                  keyType = 'SYMBOL',  
                  ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                  pAdjustMethod = 'fdr',  #指定 p 值校正方法
                  pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                  qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                  readable = FALSE)  #
C2.go <-summary(C2.go)  #8
write.table(C2.go,"./temporal_lobe\\subcluster\\Excitatory\\monocle3\\Nor\\diff\\C2.go.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
C2.go_top=C2.go %>% group_by(ONTOLOGY) %>% top_n(n=(-5),wt=p.adjust)

C2.go_top <- C2.go_top[c(1:4),]
C2.go_top$Description<- factor(C2.go_top$Description,levels = C2.go_top$Description)
C2.go_top$Group <- rep("C2",dim(C2.go_top)[1])
ggplot(C2.go_top, aes(p.adjust, Description)) +
  geom_point(aes(color=p.adjust, size=Count))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text.y = element_text(colour = "black"))+
  scale_color_gradient(low='#6699CC',high='#CC3333')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))

###*******簇3-------
C3 <- rownames(pt.matrix)[row_order(dhm)[[3]]]
library(clusterProfiler)
library(org.Hs.eg.db)
C3.go <- enrichGO(gene = C3,  #基因列表文件中的基因名称
                  OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                  keyType = 'SYMBOL',  
                  ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                  pAdjustMethod = 'fdr',  #指定 p 值校正方法
                  pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                  qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                  readable = FALSE)  #
C3.go <-summary(C3.go) 
write.table(C3.go,"./temporal_lobe\\subcluster\\Excitatory\\monocle3\\Nor\\diff\\C3.go.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
C3.go_top=C3.go %>% group_by(ONTOLOGY) %>% top_n(n=(-5),wt=p.adjust)

C3.go_top <- C3.go_top[c(1:5),]
C3.go_top$Description<- factor(C3.go_top$Description,levels = C3.go_top$Description)
C3.go_top$Group <- rep("C3",dim(C3.go_top)[1])
###*******簇4-------
C4 <- rownames(pt.matrix)[row_order(dhm)[[4]]]
library(clusterProfiler)
library(org.Hs.eg.db)
C4.go <- enrichGO(gene = C4,  #基因列表文件中的基因名称
                  OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                  keyType = 'SYMBOL',  
                  ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                  pAdjustMethod = 'fdr',  #指定 p 值校正方法
                  pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                  qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                  readable = FALSE)  #
C4.go <-summary(C4.go)  
write.table(C4.go,"./temporal_lobe\\subcluster\\Excitatory\\monocle3\\Nor\\diff\\C4.go.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
C4.go_top=C4.go %>% group_by(ONTOLOGY) %>% top_n(n=(-5),wt=p.adjust)

C4.go_top <- C4.go[C4.go$Description %in%c("regulation of excitatory synapse assembly","regulation of neurogenesis","excitatory synapse assembly","central nervous system projection neuron axonogenesis","regulation of dendrite development"),]
C4.go_top$Description<- factor(C4.go_top$Description,levels = C4.go_top$Description)
C4.go_top$Group <- rep("C4",dim(C4.go_top)[1])
exp <- rbind(C1.go_top,C2.go_top,C3.go_top,C4.go_top)
ggplot(exp, aes(Group, Description)) +
  geom_point(aes(color=p.adjust, size=Count))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low='#6699CC',high='#CC3333')+guides(size=guide_legend(order=1))
write.table(exp,"./temporal_lobe\\subcluster\\Excitatory\\monocle3\\Nor\\diff\\Exp_c1234.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
