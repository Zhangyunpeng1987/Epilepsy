##***伪时间轨迹分析------------
Excitatory_sub <- readRDS(file = "./Hippo\\subcluster\\Excitatory\\anno\\Excitatory_sub.rds")
cols<- c("#AAD0E3","#277AB4","#FABF74","#FB7D1A","#CAB4D6","#B5DF90","#3AA12F","#693C9A")
table(Excitatory_sub$group)
A_sub <- subset(Excitatory_sub,group %in% "A")
table(A_sub$Anno_Idents)
table(Excitatory_sub$group)
##monocle3-------------
##A-sub-----------
library(Seurat)
library(monocle3)
DimPlot(A_sub, reduction='umap',group.by="Anno_Idents", pt.size=1.5,label=F,label.size = 5,raster = F,
        cols = cols)
data <- GetAssayData(A_sub, assay = "RNA", slot = "counts")  
cell_metadata <- A_sub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)     
# #pca降维
cds <- reduce_dimension(cds,preprocess_method = "PCA") 
# #umap降维
cds <- reduce_dimension(cds, reduction_method="UMAP")
# #聚类
cds <- cluster_cells(cds, resolution=1e-5)
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(A_sub, reduction = "umap")
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
get_earliest_principal_node  <- function(cds, time_bin="CA3_Ex"){
  cell_ids <- which(colData(cds)[, "Anno_Idents"] == time_bin)
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
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
           label_cell_groups = F,
           label_leaves=F,
           label_branch_points=F,
           group_label_size=1)+scale_colour_manual(values = c("#8ECFC9","#FFBE7A"))                                #这个图就可以看出细胞直接的分化轨迹了
plot_cells(cds,
           reduction_method = "UMAP",
           #trajectory_graph_color="blue",
           color_cells_by = "Anno_Idents",
           label_groups_by_cluster=T,
           label_cell_groups = F,
           label_leaves=F,
           label_branch_points=F,
           group_label_size=1)+scale_colour_manual(values = cols) 
saveRDS(cds,file="./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\Excitatory_monocle3_cds.rds")
cds <- readRDS(file="./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\Excitatory_monocle3_cds.rds")

####*差异分析------------
pr_graph_test_res = graph_test(cds, neighbor_graph="principal_graph", cores=1)
write.table(pr_graph_test_res,file="./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\pr_graph_test_res.txt",col.names=T,row.names=T,sep="\t",quote=F)
pr_graph_test_res <- read.table(file="./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\pr_graph_test_res.txt",header = T,sep="\t")
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
length(pr_deg_ids)  
###*关键基因展示----------
Track_genes_sig <- pr_graph_test_res %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
cols<- c("#AAD0E3","#277AB4","#FABF74","#FB7D1A","#CAB4D6","#B5DF90","#3AA12F","#693C9A")
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="Anno_Idents", 
                         min_expr=0.5, ncol = 2)+
  scale_colour_manual(values = cols)
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\Track_genes_sig.pdf",height = 5,width = 5)

####*找模块-----------
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-03)
write.table(gene_module_df,file="./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\gene_module_df.txt",col.names=T,row.names=F,sep="\t",quote=F)
gene_module_df <- read.table(file="./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\gene_module_df.txt",header=T,sep="\t")

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$Anno_Idents)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
p <- pheatmap::pheatmap(agg_mat,
                        #color = colorRampPalette(brewer.pal(n = 7, name = "RdPu"))(100),
                        color = colorRampPalette(c('#FFFFFF','#8ECFC9'))(100),
                        treeheight_row=0,
                        treeheight_col=0,
                        scale="column", clustering_method="ward.D2",
                        cellwidth = 20, cellheight = 15)
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\Module_pheatmap.pdf",p,height = 5,width = 5)

library("RColorBrewer")
#####*细胞类型相关模块中基因提取以及功能富集-----------
#####*
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(aPEAR)
install.packages("ggupset")
library(ggupset)
names(table(A_sub$Anno_Idents))
go.all <- data.frame(matrix(ncol = 9, nrow = 0))
for (i in 1:length(names(table(A_sub$Anno_Idents)))) {
  gene_module <- rownames(agg_mat)[order(agg_mat[,colnames(agg_mat) %in% names(table(A_sub$Anno_Idents))[8]],decreasing = T)[1:3]]
  gene_module <- as.data.frame(str_split(gene_module," ",simplify = F))[2,]
  genes_exp <- gene_module_df %>% filter(module %in% gene_module)
  genes <- genes_exp$id
  enrich.go <- enrichGO(gene = genes,  #基因列表文件中的基因名称
                        OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                        keyType = 'SYMBOL',  
                        ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                        pAdjustMethod = 'fdr',  #指定 p 值校正方法
                        pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                        qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                        readable = FALSE)
  enrich.go <-as.data.frame(enrich.go)  #
  setwd("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub")
  name=paste(names(table(A_sub$Anno_Idents)),"enrich.txt",sep = "_")
  write.table(enrich.go,name[i],sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
  setwd("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub")
  enrich.go$Description<- factor(enrich.go$Description,levels = enrich.go$Description)
  enrich.go$Group <- rep(names(table(A_sub$Anno_Idents))[i],dim(enrich.go)[1])
  go.all <- rbind(enrich.go,go.all)
  
}
write.table(go.all,"./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\go.all.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
go.all <- read.table("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\go.all.txt",sep = '\t',header = T)
table(go.all$Group)

###CA_Ex_GAPDH----------
CA_Ex_GAPDH <- go.all[go.all$Group %in% "CA_Ex_GAPDH",]
CA_Ex_GAPDH <- CA_Ex_GAPDH %>% slice_min(order_by = p.adjust, n = 50)
enrichmentNetwork(CA_Ex_GAPDH, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#FB7D1A",'#944102'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\CA_Ex_GAPDH_go_top50.pdf",height = 6,width = 6)

###CA1_Ex----------
CA1_Ex <- go.all[go.all$Group %in% "CA1_Ex",]
CA1_Ex <- CA1_Ex %>% slice_min(order_by = p.adjust, n = 50)
enrichmentNetwork(CA1_Ex, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#FABF74",'#C3A279'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\CA1_Ex_go_top50.pdf",height = 6,width = 6)


###CA3_Ex----------
CA3_Ex <- go.all[go.all$Group %in% "CA3_Ex",] ##7
enrichmentNetwork(CA3_Ex, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#CAB4D6",'#b4d6b9'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\CA3_Ex_go_top50.pdf",height = 6,width = 6)


table(go.all$Group)
###DG_Ex1----------
DG_Ex1 <- go.all[go.all$Group %in% "DG_Ex1",] ##14
#DG_Ex1 <- DG_Ex1 %>% slice_min(order_by = p.adjust, n = 50)
enrichmentNetwork(DG_Ex1, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#AAD0E3",'#367FA3'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\DG_Ex1_go_top50.pdf",height = 6,width = 6)

table(go.all$Group)
###DG_Ex2----------
DG_Ex2 <- go.all[go.all$Group %in% "DG_Ex2",] ##35
enrichmentNetwork(DG_Ex2, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#277AB4",'#587D98'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\DG_Ex2_go_top50.pdf",height = 6,width = 6)

table(go.all$Group)
###DG_Ex3----------
DG_Ex3 <- go.all[go.all$Group %in% "DG_Ex3",] ##29
enrichmentNetwork(DG_Ex3, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#B5DF90",'#4DA3A5'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\DG_Ex3_go_top50.pdf",height = 6,width = 6)

table(go.all$Group)
###DG_Ex4----------
DG_Ex4 <- go.all[go.all$Group %in% "DG_Ex4",] ##18
enrichmentNetwork(DG_Ex4, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#693C9A",'#BB9ED9'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\DG_Ex4_go_top50.pdf",height = 6,width = 6)
table(go.all$Group)
###DG_Ex5----------
DG_Ex5 <- go.all[go.all$Group %in% "DG_Ex5",] ##238
DG_Ex5 <- DG_Ex5 %>% slice_min(order_by = p.adjust, n = 50)
enrichmentNetwork(DG_Ex5, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#3AA12F",'#4DA3A5'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\DG_Ex5_go_top50.pdf",height = 6,width = 6)

##P-sub-----------
Excitatory_sub <- readRDS(file = "./Hippo\\subcluster\\Excitatory\\anno\\Excitatory_sub.rds")
cols<- c("#AAD0E3","#277AB4","#FABF74","#FB7D1A","#CAB4D6","#B5DF90","#3AA12F","#693C9A")
table(Excitatory_sub$group)
P_sub <- subset(Excitatory_sub,group %in% "P")
table(P_sub$Anno_Idents)
table(Excitatory_sub$group)
library(Seurat)
library(monocle3)
DimPlot(P_sub, reduction='umap',group.by="Anno_Idents", pt.size=1,label=F,label.size = 5,raster = F,
        cols = cols)
data <- GetAssayData(P_sub, assay = "RNA", slot = "counts")  
cell_metadata <- P_sub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)  

cds <- reduce_dimension(cds,preprocess_method = "PCA") 

cds <- reduce_dimension(cds, reduction_method="UMAP")

cds <- cluster_cells(cds, resolution=1e-5)
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(P_sub, reduction = "umap")
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
get_earliest_principal_node  <- function(cds, time_bin="CA3_Ex"){
  cell_ids <- which(colData(cds)[, "Anno_Idents"] == time_bin)
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))  
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
           label_cell_groups = F,
           label_leaves=F,
           label_branch_points=F,
           group_label_size=1)+scale_colour_manual(values = c("#8ECFC9","#FFBE7A"))                                #这个图就可以看出细胞直接的分化轨迹了
plot_cells(cds,
           reduction_method = "UMAP",
           #trajectory_graph_color="blue",
           color_cells_by = "Anno_Idents",
           label_groups_by_cluster=T,
           label_cell_groups = F,
           label_leaves=F,
           label_branch_points=F,
           group_label_size=1)+scale_colour_manual(values = cols) 
saveRDS(cds,file="./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\Excitatory_monocle3_cds.rds")
cds <- readRDS(file="./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\Excitatory_monocle3_cds.rds")

####*差异分析------------
pr_graph_test_res = graph_test(cds, neighbor_graph="principal_graph", cores=1)
write.table(pr_graph_test_res,file="./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\pr_graph_test_res.txt",col.names=T,row.names=T,sep="\t",quote=F)
pr_graph_test_res <- read.table(file="./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\pr_graph_test_res.txt",header = T,sep="\t")
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05 & morans_I>0.5))
length(pr_deg_ids)  # 10749
###*关键基因展示----------
Track_genes_sig <- pr_graph_test_res %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#基因表达趋势图
cols<- c("#AAD0E3","#277AB4","#FABF74","#FB7D1A","#CAB4D6","#B5DF90","#3AA12F","#693C9A")
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="Anno_Idents", 
                         min_expr=0.5, ncol = 2)+
  scale_colour_manual(values = cols)
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\Track_genes_sig.pdf",height = 5,width = 5)

####*找模块-----------
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-03)
write.table(gene_module_df,file="./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\gene_module_df.txt",col.names=T,row.names=F,sep="\t",quote=F)
gene_module_df <- read.table(file="./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\gene_module_df.txt",header=T,sep="\t")

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$Anno_Idents)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
p <- pheatmap::pheatmap(agg_mat,
                        #color = colorRampPalette(brewer.pal(n = 7, name = "RdPu"))(100),
                        color = colorRampPalette(c('#FFFFFF','#FFBE7A'))(100),
                        treeheight_row=0,
                        treeheight_col=0,
                        scale="column", clustering_method="ward.D2",
                        cellwidth = 20, cellheight = 15)
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\Module_pheatmap.pdf",p,height = 5,width = 5)

library("RColorBrewer")
#####*细胞类型相关模块中基因提取以及功能富集-----------
#####*
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(aPEAR)
library(ggupset)
names(table(P_sub$Anno_Idents))
go.all <- data.frame(matrix(ncol = 9, nrow = 0))
for (i in 1:length(names(table(P_sub$Anno_Idents)))) {
  gene_module <- rownames(agg_mat)[order(agg_mat[,colnames(agg_mat) %in% names(table(P_sub$Anno_Idents))[i]],decreasing = T)[1:3]]
  gene_module <- as.data.frame(str_split(gene_module," ",simplify = F))[2,]
  genes_exp <- gene_module_df %>% filter(module %in% gene_module)
  genes <- genes_exp$id
  enrich.go <- enrichGO(gene = genes,  #基因列表文件中的基因名称
                        OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                        keyType = 'SYMBOL',  
                        ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                        pAdjustMethod = 'fdr',  #指定 p 值校正方法
                        pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                        qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                        readable = FALSE)
  enrich.go <-as.data.frame(enrich.go)  #
  setwd("./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub1")
  name=paste(names(table(P_sub$Anno_Idents)),"enrich.txt",sep = "_")
  write.table(enrich.go,name[i],sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
  setwd("./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub1")
  enrich.go$Description<- factor(enrich.go$Description,levels = enrich.go$Description)
  enrich.go$Group <- rep(names(table(P_sub$Anno_Idents))[i],dim(enrich.go)[1])
  go.all <- rbind(enrich.go,go.all)
  
}
write.table(go.all,"./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub1\\go.all.txt",sep = '\t',col.names = T,row.names = F)
go.all <- read.table("./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub1\\go.all.txt",sep = '\t',header = T)
table(go.all$Group)  #1472
###CA_Ex_GAPDH----------
CA_Ex_GAPDH <- go.all[go.all$Group %in% "CA_Ex_GAPDH",]  #595
CA_Ex_GAPDH <- CA_Ex_GAPDH %>% slice_min(order_by = p.adjust, n = 50)
enrichmentNetwork(CA_Ex_GAPDH, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#FB7D1A",'#944102'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\CA_Ex_GAPDH_go_top50.pdf",height = 6,width = 6)

###CA1_Ex----------
CA1_Ex <- go.all[go.all$Group %in% "CA1_Ex",]#181
CA1_Ex <- CA1_Ex %>% slice_min(order_by = p.adjust, n = 50)
enrichmentNetwork(CA1_Ex, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#FABF74",'#C3A279'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\CA1_Ex_go_top50.pdf",height = 6,width = 6)


###CA3_Ex----------
CA3_Ex <- go.all[go.all$Group %in% "CA3_Ex",] ##155
CA3_Ex <- CA3_Ex %>% slice_min(order_by = p.adjust, n = 50)
enrichmentNetwork(CA3_Ex, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#CAB4D6",'#b4d6b9'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\CA3_Ex_go_top50.pdf",height = 6,width = 6)


table(go.all$Group)
###DG_Ex1----------
DG_Ex1 <- go.all[go.all$Group %in% "DG_Ex1",] ##39
enrichmentNetwork(DG_Ex1, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#AAD0E3",'#367FA3'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\DG_Ex1_go_top50.pdf",height = 6,width = 6)
exp <- data.frame(matrix(ncol = 2))
colnames(exp) <- c("genes","Go")
for (i in 1:length(DG_Ex1$Description)) {
  genes <- DG_Ex1[i,8]
  genes <- str_split(genes,"/")[[1]]
  genes <- as.data.frame(genes)
  genes$Go <- DG_Ex1[i,2]
  exp <- rbind(exp,genes)
}
exp$Freq <- "1"
exp <- exp[-1,]
head(exp)
gene1 <- unique(exp$genes[exp$Go %in% "macroautophagy"])
intersect(gene1,P_up_DEG)
P_up_DEG <- clusterdeg[clusterdeg$avg_log2FC >0.25 & clusterdeg$p_val_adj <0.01,7]

table(go.all$Group)
###DG_Ex2(未富集到）----------
DG_Ex2 <- go.all[go.all$Group %in% "DG_Ex2",] ##0
enrichmentNetwork(DG_Ex2, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#277AB4",'#587D98'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\DG_Ex2_go_top50.pdf",height = 6,width = 6)

table(go.all$Group)
###DG_Ex3----------
DG_Ex3 <- go.all[go.all$Group %in% "DG_Ex3",] ##43
enrichmentNetwork(DG_Ex3, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#B5DF90",'#4DA3A5'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\DG_Ex3_go_top50.pdf",height = 6,width = 6)

DG_Ex3



table(go.all$Group)
###DG_Ex4----------
DG_Ex4 <- go.all[go.all$Group %in% "DG_Ex4",] ##140
DG_Ex4 <- DG_Ex4 %>% slice_min(order_by = p.adjust, n = 50)
enrichmentNetwork(DG_Ex4, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#693C9A",'#BB9ED9'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\DG_Ex4_go_top50.pdf",height = 6,width = 6)
table(go.all$Group)
###DG_Ex5----------
DG_Ex5 <- go.all[go.all$Group %in% "DG_Ex5",] ##319
DG_Ex5 <- DG_Ex5 %>% slice_min(order_by = p.adjust, n = 50)
enrichmentNetwork(DG_Ex5, 
                  simMethod = 'jaccard',
                  clustMethod = 'markov',
                  colorBy = 'pvalue', 
                  colorType = 'pval',
                  nodeSize = 'Count',
                  #innerCutoff = 0.05,
                  fontSize = 3,
                  drawEllipses = F,
                  repelLabels =F,
                  plotOnly =F
)+
  scale_color_gradientn(colours = c("#3AA12F",'#4DA3A5'),
                        name = "logPval")
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\DG_Ex5_go_top50.pdf",height = 6,width = 6)

###兴奋性神经元在前后端差异明显，求每个细胞亚型在前后端的差异表达基因---------------
devtools::install_github("sajuukLyu/ggunchull", type = "source")
devtools::install_github('junjunlab/scRNAtoolVis')
library(jjAnno)
library(scRNAtoolVis)
library(devtools)
library(BiocManager)
library(ComplexHeatmap)
library(ggunchull)
library(jjAnno)
library(scRNAtoolVis)
Excitatory_sub <- readRDS(file = "./Hippo\\subcluster\\Excitatory\\anno\\Excitatory_sub.rds")
cols<- c("#AAD0E3","#277AB4","#FABF74","#FB7D1A","#CAB4D6","#B5DF90","#3AA12F","#693C9A")
table(Excitatory_sub$group)
table(Excitatory_sub$Anno_Idents)
table(Excitatory_sub$Anno_Idents,Excitatory_sub$group)

clusterdeg <- data.frame()
for (i in 1:length(names(table(Excitatory_sub$Anno_Idents)))) {
  cell_exp <- subset(Excitatory_sub, Anno_Idents==names(table(Excitatory_sub$Anno_Idents))[i])
  print(names(table(Excitatory_sub$Anno_Idents))[i])
  print(table(cell_exp$group))
  diff_cell_exp <- FindMarkers(cell_exp, min.pct = 0.1, 
                               logfc.threshold = 0.25,
                               test.use = "wilcox",
                               group.by = "group",
                               ident.1 ="P", ##上调
                               ident.2="A")##下调
  diff_cell_exp$cluster=names(table(Excitatory_sub$Anno_Idents))[i]
  clusterdeg=rbind(diff_cell_exp,clusterdeg)
  #setwd("./Hippo\\subcluster\\Excitatory\\diff")
  #name=paste(names(table(Excitatory_sub$Anno_Idents)),"P-A-wilcox.txt",sep = "_")
  #write.table(diff_cell_exp,name[i],sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
}
clusterdeg$gene=rownames(clusterdeg)  ##对数据增加一列
write.table(clusterdeg,"./Hippo\\subcluster\\Excitatory\\diff\\clusterdeg.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
clusterdeg <- read.table("./Hippo\\subcluster\\Excitatory\\diff\\clusterdeg.txt",sep = '\t',header = T)
table(clusterdeg$cluster)
###显著差异的基因-------------
clusterdeg <- clusterdeg[clusterdeg$p_val_adj <0.05,]
table(clusterdeg$cluster)
cols<- c("#FB7D1A","#FABF74","#CAB4D6","#AAD0E3","#277AB4","#B5DF90","#693C9A","#3AA12F")
P_sub_pr_graph_test_res <- read.table(file="./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\pr_graph_test_res.txt",header = T,sep="\t")
P_sub_pr_deg_ids <- row.names(subset(P_sub_pr_graph_test_res, q_value < 0.05 & morans_I>0.5))
A_sub_pr_graph_test_res <- read.table(file="./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\pr_graph_test_res.txt",header = T,sep="\t")
A_sub_pr_deg_ids <- row.names(subset(A_sub_pr_graph_test_res, q_value < 0.05 & morans_I>0.5))

P_up_DEG <- clusterdeg[clusterdeg$avg_log2FC >0.25 & clusterdeg$p_val_adj <0.01,7]
A_up_DEG <- clusterdeg[clusterdeg$avg_log2FC < (-0.25) & clusterdeg$p_val_adj <0.01,7]
p_common <- intersect(P_up_DEG,P_sub_pr_deg_ids)
p_common  #"GRM1" "NBEA" "TLL1"
A_common <- intersect(A_up_DEG,A_sub_pr_deg_ids)
A_common  #"VWC2L"  "TSHZ2"  "COBLL1"
#"TSHZ2"
library(jjAnno)
library(scRNAtoolVis)
library(devtools)
library(BiocManager)
library(ComplexHeatmap)
library(ggunchull)
library(jjAnno)
library(scRNAtoolVis)

jjDotPlot(object = Excitatory_sub,
          gene = c(p_common,A_common),
          id = 'group',
          xtree = F,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1,
          point.shape = 22)
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\group_ap_轨迹基因vs细胞类型在ap差异基因交集.pdf",height = 5,width = 15)
jjDotPlot(object = Excitatory_sub,
          gene = c(p_common,A_common),
          id = 'Anno_Idents',
          xtree = F,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1,
          point.shape = 22)
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\Anno_ap_轨迹基因vs细胞类型在ap差异基因交集.pdf",height = 5,width = 15)
cds <- readRDS(file="./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\Excitatory_monocle3_cds.rds")
cols<- c("#AAD0E3","#277AB4","#FABF74","#FB7D1A","#CAB4D6","#B5DF90","#3AA12F","#693C9A")
plot_genes_in_pseudotime(cds[c("TLL1"),], color_cells_by="Anno_Idents", 
                         min_expr=0.5, ncol = 2)+
  scale_colour_manual(values = cols)
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\P_sub\\TLL1_ap_轨迹基因vs细胞类型在ap差异基因交集.pdf",height = 5,width = 5)

cds <- readRDS(file="./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\Excitatory_monocle3_cds.rds")
cols<- c("#AAD0E3","#277AB4","#FABF74","#FB7D1A","#CAB4D6","#B5DF90","#3AA12F","#693C9A")
plot_genes_in_pseudotime(cds[c("TSHZ2"),], color_cells_by="Anno_Idents", 
                         min_expr=0.5, ncol = 2)+
  scale_colour_manual(values = cols)
ggsave("./Hippo\\subcluster\\Excitatory\\monocle3\\A_sub\\TSHZ2_ap_轨迹基因vs细胞类型在ap差异基因交集.pdf",height = 5,width = 5)
