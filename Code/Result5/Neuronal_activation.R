####神经元活化------------------
target <- read.table("/home/data/t120425/epilepsy/neuron/temp_hippo_node.txt",sep = "\t",header = T)
unique(target$Group)
load(file = "/home/data/t120425/epilepsy/function_genelist.RData")
for(i in 1:11){
  genelist[[i]] <- intersect(genelist[[i]],target$TF)
}
genelist[["Neuronal_activation"]]
#"EGR1"  "EGR3"  "FOS"   "FOSL2"
table(target[target$TF %in% genelist[["Neuronal_activation"]],4])
# Astro_A           Astro_P              EX_A              IN_A              IN_P           Micro_A          Micro_Ep 
# 68                41                93                16                 6                11                21 
# Oligodendrocyte_A 
# 16  
loom <- open_loom('./Hippo\\subcluster\\Excitatory\\EX_A\\out_SCENIC.loom') 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')

ex_Nor_auc <- as.data.frame(regulonAUC@assays@data$AUC)
rownames(ex_Nor_auc) <- str_split(rownames(ex_Nor_auc),'\\(',simplify = T)[,1]

ex_nor_TF <- unique(target[target$Group %in% "EX_A",1])
ex_nor_TF
ex_nor_TF_auc <- ex_Nor_auc[genelist[["Neuronal_activation"]],]
ex_nor_TF_auc <- t(ex_nor_TF_auc)
library(monocle3)
cds <- readRDS(file="./Excitatory_monocle3_cds.rds")
pseudotime <- pseudotime(cds)
pseudotime <- as.data.frame(pseudotime)
ex_nor_dat <- cbind(ex_nor_TF_auc,pseudotime)
p_AM <- ggplot() + 
  geom_smooth(data = ex_nor_dat, aes(x = pseudotime, y = EGR1, color = '#E59CC4'),  method = "loess", se = T) +
  ggpubr::stat_cor(data = ex_nor_dat, aes(x = pseudotime, y = EGR1, color = '#E59CC4'),size=3.5,label.y.npc = 0.1)+
  
  geom_smooth(data = ex_nor_dat, aes(x = pseudotime, y = EGR3, color = "#476D87"),  method = "loess", se = T) +
  ggpubr::stat_cor(data = ex_nor_dat, aes(x = pseudotime, y = EGR3, color = "#476D87"),size=3.5,label.y.npc = 0.12)+
  
  geom_smooth(data = ex_nor_dat, aes(x = pseudotime, y = FOS, color = '#B53E2B'),  method = "loess", se = T) +
  ggpubr::stat_cor(data = ex_nor_dat, aes(x = pseudotime, y = FOS, color = '#B53E2B'),size=3.5,label.y.npc = 0.14)+
  
  geom_smooth(data = ex_nor_dat, aes(x = pseudotime, y = FOSL2, color = '#68A180'),  method = "loess", se = T) +
  ggpubr::stat_cor(data = ex_nor_dat, aes(x = pseudotime, y = FOSL2, color = '#68A180'),size=3.5,label.y.npc = 0.18)+
  scale_color_manual(values = c('#E59CC4',"#476D87",'#B53E2B','#68A180')) +  
  coord_cartesian(ylim = c(0, 0.15))+
  theme(panel.background = element_rect(fill = "white", colour = NA),  
        panel.border = element_rect(fill = NA,  colour = "black"),
        panel.grid.minor = element_line(linewidth = rel(1)),
        legend.key = element_rect(fill = "black", 
                                  colour = NA))+
  labs(title = "", x = "Pseudotime", y = "TF activatity", color = "TF")+theme(legend.position = "none") 

print(p_AM)

####靶基因功能富集-----------------
###海马体-----------
unique(target$Group)
hippo <- target[target$TF %in% genelist[["Neuronal_activation"]],]
hippo <- hippo[hippo$Group %in% c("EX_A","EX_P","IN_A","IN_P","Astro_A","Astro_P",
                                  "Micro_A","Micro_P","Oligodendrocyte_A","Oligodendrocyte_P"),]
table(hippo$Group)
# Astro_A           Astro_P              EX_A              IN_A              IN_P           Micro_A Oligodendrocyte_A 
# 68                41                93                16                 6                11                16 

hippo_gene <- unique(c(genelist[["Neuronal_activation"]],hippo$TargetGene))

enrich.go <- enrichGO(gene = hippo_gene,  #基因列表文件中的基因名称
                      OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                      keyType = 'SYMBOL', 
                      ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                      qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                      readable = FALSE)
enrich.go <-as.data.frame(enrich.go)  #
enrich.go$Description<- factor(enrich.go$Description,levels = enrich.go$Description)
hippo_go <- enrich.go
save(hippo_go,file = "/home/data/t120425/epilepsy/neuron/hippo_go.TF.RData")




load(file ="/home/data/t120425/epilepsy/neuron/hippo_go.TF.RData")
hippo_go$logPvalue <- -log10(hippo_go$pvalue)
hippo_go <- hippo_go %>% top_n(10,logPvalue)
head(hippo_go)
library(ggplot2)
library(ggsci)
write.table(hippo_go,file = "/home/data/t120425/epilepsy/neuron/hippo_go.TF.txt",sep = "\t")


#功能富集花瓣图---------------

mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  plot.title = element_text(size = 14,
                            hjust = 0.5,
                            face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5,
                       r = 10,
                       l = 5.5,
                       b = 5.5)
)
mytheme2 <- mytheme + theme(axis.text.y = element_blank())

hippo_go <- hippo_go[order(hippo_go$p.adjust, decreasing = T), ]
hippo_go$Description <- factor(hippo_go$Description,levels = hippo_go$Description)
p1 <- ggplot(data = hippo_go, aes(x = -log10(pvalue), y = Description, fill = -log10(pvalue))) +
  scale_fill_gradient(low = "#F5E3DE",high = "#831A1F") +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.7) +
  geom_col(colour="black",just = 0.5,width = 0.5)+
  labs(x = "-log10(pvalue)", y = "", title = "Hippo enrichment") +
  geom_text(aes(x = 0.03, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic()
p1

epilepsy_risk_gene <- read.csv("/home/data/t120425/epilepsy/neuron/epilepsy_risk_gene.csv",sep = ",",header = T)
epilepsy_risk_gene <- unique(epilepsy_risk_gene$Approved.Symbol)

gene <- intersect(GO_hippo$geneID,epilepsy_risk_gene)
#"GABRA1" "NTRK2"  "GABRG2" "GRIN1"  "PLCB1" 
intersect(FOSL2_target,c("GABRA1","NTRK2","GABRG2","GRIN1","PLCB1"))
#"GABRA1" "NTRK2"  "GABRG2" "GRIN1"  "PLCB1"
intersect(c("EGR1","EGR3","FOS","FOSL2"),epilepsy_risk_gene)

hippo_go_gene <- c(intersect(GO_hippo$geneID,epilepsy_risk_gene),
                   intersect(GO_hippo$geneID,FOS_target),
                   intersect(GO_hippo$geneID,EGR1_target),intersect(GO_hippo$geneID,EGR3_target))

###颞叶皮层-----------
unique(target$Group)
temp <- target[target$TF %in% genelist[["Neuronal_activation"]],]
temp <- temp[temp$Group %in% c("Ex_Ep","In_Ep","Micro_Ep","Astro_Ep","ODCs_Ep" ,"Oligodendrocyte_Ep"),]
table(temp$Group)
# Micro_Ep 
# 21
temp_gene <- unique(c(genelist[["Neuronal_activation"]],temp$TargetGene)) #25

enrich.go <- enrichGO(gene = temp_gene,  #基因列表文件中的基因名称
                      OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                      keyType = 'SYMBOL',  
                      ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                      qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                      readable = FALSE)
enrich.go <-as.data.frame(enrich.go)  #
enrich.go$Description<- factor(enrich.go$Description,levels = enrich.go$Description)
temp_go <- enrich.go
save(temp_go,file = "/home/data/t120425/epilepsy/neuron/temp_go.TF.RData")




load(file ="/home/data/t120425/epilepsy/neuron/temp_go.TF.RData")
temp_go$logPvalue <- -log10(temp_go$pvalue)
temp_go <- temp_go %>% top_n(10,logPvalue)
head(temp_go)
library(ggplot2)
library(ggsci)
write.table(temp_go,file = "/home/data/t120425/epilepsy/neuron/temp_go.TF.txt",sep = "\t")

mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  plot.title = element_text(size = 14,
                            hjust = 0.5,
                            face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5,
                       r = 10,
                       l = 5.5,
                       b = 5.5)
)
mytheme2 <- mytheme + theme(axis.text.y = element_blank())

temp_go <- temp_go[order(temp_go$p.adjust, decreasing = T), ]
temp_go$Description <- factor(temp_go$Description,levels = temp_go$Description)
p1 <- ggplot(data = temp_go, aes(x = -log10(pvalue), y = Description, fill = -log10(pvalue))) +
  scale_fill_gradient(low = "#F5E3DE",high = "#831A1F") +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.7) +
  geom_col(colour="black",just = 0.5,width = 0.5)+
  labs(x = "-log10(pvalue)", y = "", title = "Temp enrichment") +
  geom_text(aes(x = 0.03, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic()
p1
