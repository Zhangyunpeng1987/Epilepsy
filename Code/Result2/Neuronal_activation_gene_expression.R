###细胞活性基因在细胞中的表达------------------
Neuronal_activation <- toupper(c("Arc", "Btg2", "Coq10b", "Crem", "Dusp1", "Dusp5", "Egr1", "Egr3", "Fbxo33", "Fos",
                                 "Fosl2", "Gadd45g", "Gmeb2", "Grasp", "Junb", "Nr4a1", "Nr4a2", "Nr4a3", "Per1", "Rgs2",
                                 "Sertad1","Tiparp"))

aver_dt <- AverageExpression(all.scRNA,
                             features = Neuronal_activation,
                             group.by = 'Anno_Idents',
                             slot = 'data')
aver_dt <- as.data.frame(aver_dt$RNA)
aver_dt[1:6,1:6]

aver_dt <- aver_dt[,c("CA_Ex_GAPDH","CA1_Ex","CA3_Ex","DG_Ex1","DG_Ex2","DG_Ex3","DG_Ex4","DG_Ex5",
                      "L2_3_Cux2","L4_L5_6_Rorb","L6",
                      "CCK","LAMP5","PVALB","SST","VIP",
                      "In_Lamp5","In_Pvalb","In_Sst","In_Vip")]
#列注释：celltype
col_anno <- data.frame(cell_anno = colnames(aver_dt),
                       row.names = colnames(aver_dt))
col_anno$cell_anno <- factor(col_anno$cell_anno,levels = c("CA_Ex_GAPDH","CA1_Ex","CA3_Ex","DG_Ex1","DG_Ex2","DG_Ex3","DG_Ex4","DG_Ex5",
                                                           "L2_3_Cux2","L4_L5_6_Rorb","L6",
                                                           "CCK","LAMP5","PVALB","SST","VIP",
                                                           "In_Lamp5","In_Pvalb","In_Sst","In_Vip"))

head(col_anno)

#热图美化：
#热图配色自定义：
mycol <- colorRampPalette(c("#108dc7", "white", "#ef8e38"))(50)

#行列注释配色自定义：

celltype_color <- c("L2_3_Cux2"="#C4A5DE96","L4_L5_6_Rorb"="#81B8DF96","L6"="#F6CAE596",
                    "In_Lamp5"="#D76364","In_Pvalb"="#8ECFC9","In_Sst"="#82B0D2","In_Vip"="#FFBE7A",
                    "CA_Ex_GAPDH"="#FB7D1A","CA1_Ex"="#FABF74","CA3_Ex"="#CAB4D6","DG_Ex1"="#AAD0E3",
                    "DG_Ex2"="#277AB4","DG_Ex3"="#B5DF90","DG_Ex4"="#693C9A","DG_Ex5"="#3AA12F",
                    "CCK"="#FFCC00","LAMP5"="#FF9900","PVALB"="#6699CC","SST"="#339999","VIP"="#CCCC99"
)

anno_col <- list(cell_anno = celltype_color
                 #gene_anno = celltype_col
)
anno_col
library(Seurat)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
pheatmap(scale(as.matrix(aver_dt)),
         scale = "row",
         cluster_rows = T,
         cluster_cols = F,
         gaps_col = c(8, 11,16),
         #cellwidth  = 1,cellwidth  = 1,
         annotation_col = col_anno,
         show_rownames = T,
         cellwidth = 20, cellheight = 20,
         #annotation_row = row_anno,
         annotation_colors = anno_col, #注释配色
         color = mycol, #热图配色
         border_color = 'white') #描边颜色
#Neuronal_activation_expression_hippo_temp.pdf
