###基因集打分总-----------------
###活性基因集--------------
Neuronal_activation <- toupper(c("Arc", "Btg2", "Coq10b", "Crem", "Dusp1", "Dusp5", "Egr1", "Egr3", "Fbxo33", "Fos",
                                 "Fosl2", "Gadd45g", "Gmeb2", "Grasp", "Junb", "Nr4a1", "Nr4a2", "Nr4a3", "Per1", "Rgs2",
                                 "Sertad1","Tiparp"))

###功能集合----------
temp_GO <- read.table("./temp/temp_go.all_top5.csv",header = T,sep = ",")
hippo_GO <- read.table("./hippo/hippo_go.all_top5.csv",header = T,sep = ",")
all_GO <- rbind(temp_GO[,c(1,2,12)],hippo_GO[,c(1,2,12)])
all_GO <- unique(all_GO) #45

unique(all_GO$BP)
# [1] "differentiation"            "signaling"                  "metabolic process"          "transport"                 
# [5] "cell junction organization" "system process"             "cell adhesion"              "developmental process" 

library(msigdbr)
GO_df = msigdbr(species = "Homo sapiens",category = "C5") %>% 
  dplyr::select(gene_symbol,gs_name,gs_subcat,gs_url)
dim(GO_df)
## [1] 1424471       3
table(GO_df$gs_subcat)
## 
##  GO:BP  GO:CC  GO:MF    HPO 
## 721379 115769 122674 464649
#GO_df = GO_df[GO_df$gs_subcat!="HPO",]
GO_id <- GO_df$gs_url %>%
  str_split("/",simplify = T) %>%
  .[,6]
GO_df$GO_id <- GO_id

save(GO_df,file = "./msigdbr_GO.RData")

load(file = "./msigdbr_GO.RData")

genelist <- list()
for (i in 1:8) {
  a <-all_GO[all_GO$BP %in% unique(all_GO$BP)[i],]
  print(length(intersect(GO_df$GO_id,a$ID)))
  print(setdiff(unique(a$ID),intersect(GO_df$GO_id,a$ID)))
  gene <- unique(GO_df[GO_df$GO_id %in% a$ID,1]$gene_symbol)
  #gene <- intersect(gene,rownames(all.scRNA))
  genelist[[i]] <- gene
  names(genelist)[i] <- unique(all_GO$BP)[i]
}

genelist[9] <- list(Neuronal_activation)
names(genelist)[9] <- "Neuronal_activation"
for (i in 1:length(genelist)) {
  all.scRNA <- AddModuleScore(object = all.scRNA, features = genelist[i],name = names(genelist)[i])
  
}
save(all.scRNA,file = "./all.scRNA.RData")
save(genelist,file = "./function_genelist.RData")
colnames(all.scRNA@meta.data)
score <- all.scRNA@meta.data[,c(6,18,28:38)]
table(score$group)
score$tissue[score$group %in% c("A","P")] <- "Hippo"
score$tissue[score$group %in% c("Nor","Ep")] <- "Temp"

score$celltype[score$Anno_Idents %in% c("CA_Ex_GAPDH","CA1_Ex","CA3_Ex","DG_Ex1","DG_Ex2","DG_Ex3","DG_Ex4","DG_Ex5",
                                        "L2_3_Cux2","L4_L5_6_Rorb","L6")] <- "Excitatory neurons"
score$celltype[score$Anno_Idents %in% c("CCK","LAMP5","PVALB","SST","VIP",
                                        "In_Lamp5","In_Pvalb","In_Sst","In_Vip")] <- "Inhibitory neurons"

hippo_score <- score[score$tissue %in% "Hippo",]

temp_score <- score[score$tissue %in% "Temp",]

###颞叶皮层-----------
temp_score <- aggregate(x = temp_score[,colnames(temp_score) != c("tissue","celltype")],
                        # Mean by group
                        by = list(temp_score$celltype),
                        FUN = mean)


table(temp_score$celltype)


library(ggradar) 
library(ggplot2)
temp_score <- temp_score[,-c(13:14)]
temp_score <- temp_score[,1:10]
#"#8EA9A9","#A98EA9"  in  ex
max(temp_score[,2:12])
ggradar(temp_score,
        base.size=18,
        grid.min = 0, #网格线最小值
        grid.mid = max(temp_score[,2:10])/2, #网格线均值
        grid.max = max(temp_score[,2:10]), #网格线最大值
        values.radar = c(0,max(temp_score[,2:10])/2,max(temp_score[,2:10])), #轴标签显示
        axis.label.size = 6,
        group.colours = c("#A98EA9","#8EA9A9"),
        group.point.size = 4,#分组点大小
        group.line.width = 1.2, #线宽
        background.circle.colour = 'grey', #背景填充色
        background.circle.transparency = 0, #背景填充不透明度(这里改为0可去掉背景填充)
        legend.position = 'right', #图例位置
        legend.text.size = 16, #图例标签大小
        fill = T, #各分组是否填充色
        centre.y = grid.min - ((1/2) * (grid.max - grid.min)),
        fill.alpha = 0.5 #分组填充色不透明度
)
#c("#CC99CC","#99CCCC")  ex in
#temp_ggradar.pdf

###海马体-----------
hippo_score <- aggregate(x = hippo_score[,colnames(hippo_score) != c("tissue","celltype")],
                         # Mean by group
                         by = list(hippo_score$celltype),
                         FUN = mean)


table(hippo_score$celltype)


library(ggradar) 
library(ggplot2)
hippo_score <- hippo_score[,-c(13:14)]
dt <- hippo_score[,1:10]
#"#8EA9A9","#A98EA9"  in  ex
hippo_score <- hippo_score[,1:10]
ggradar(hippo_score,
        base.size=18,
        grid.min = 0, #网格线最小值
        grid.mid = max(hippo_score[,2:10])/2, #网格线均值
        grid.max = max(hippo_score[,2:10]), #网格线最大值
        values.radar = c(0,max(hippo_score[,2:10])/2,max(hippo_score[,2:10])), #轴标签显示
        axis.label.size = 6,
        group.colours = c("#A9728E","#C58E72"),
        group.point.size = 4,#分组点大小
        group.line.width = 1.2, #线宽
        background.circle.colour = 'grey', #背景填充色
        background.circle.transparency = 0, #背景填充不透明度(这里改为0可去掉背景填充)
        legend.position = 'right', #图例位置
        legend.text.size = 16, #图例标签大小
        fill = T, #各分组是否填充色
        centre.y = grid.min - ((1/2) * (grid.max - grid.min)),
        fill.alpha = 0.5 #分组填充色不透明度
)
#c("#A9728E","#C58E72") ex in
#hippo_ggradar.pdf
