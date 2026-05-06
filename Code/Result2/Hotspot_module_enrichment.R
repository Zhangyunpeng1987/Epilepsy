###temp_hotspot_module_GO-----------
library(sceasy)
library(reticulate)
Sys.setenv(RETICULATE_MINICONDA_PATH="~/miniconda")
reticulate::use_condaenv("hotspotenv")
library(SeuratDisk)
SaveH5Seurat(Temp_scRNA, filename = "./temp/temp.h5Seurat")

Excitatory <- readRDS(file = "./temp/Excitatory/Excitatory_anno_cluster_1.rds")
names(table(Excitatory$Anno_Idents))
Inhibitory <- readRDS(file = "./temp/Inhibitory/Inhibitory_anno_cluster.rds")
Temp_scRNA <- merge(x=Excitatory,y=Inhibitory)

Temp_genemodule <- read.csv("./temp/temp_genemodule.csv")
Temp_genemodule <- Temp_genemodule[!is.na(Temp_genemodule$Module),]
Temp_genemodule <- Temp_genemodule[order(Temp_genemodule$Module),]
Temp_genemodule <- Temp_genemodule[-which(Temp_genemodule$Module ==(-1)),]

Temp_genemodule$Module_defined[Temp_genemodule$Module %in% c(2,9)] <- "M1"
Temp_genemodule$Module_defined[Temp_genemodule$Module %in% c(1,3,11)] <- "M2"
Temp_genemodule$Module_defined[Temp_genemodule$Module %in% c(5)] <- "M3"
Temp_genemodule$Module_defined[Temp_genemodule$Module %in% c(6)] <- "M4"
Temp_genemodule$Module_defined[Temp_genemodule$Module %in% c(12)] <- "M5"
Temp_genemodule$Module_defined[Temp_genemodule$Module %in% c(4,10)] <- "M6"
Temp_genemodule$Module_defined[Temp_genemodule$Module %in% c(8)] <- "M7"
Temp_genemodule$Module_defined[Temp_genemodule$Module %in% c(7)] <- "M8"
Module_defined <- names(table(Temp_genemodule$Module_defined))

unique(Temp_genemodule$Module)
Temp_genemodule$Module_defined <- factor(Temp_genemodule$Module_defined,
                                         levels = c("M1"))
Module_defined <- names(table(Temp_genemodule$Module_defined))
library(clusterProfiler)
library(org.Hs.eg.db)

go.all <- data.frame(matrix(ncol = 5, nrow = 0))
for (i in 1:length(Module_defined)) {
  gene <- Temp_genemodule$Gene[Temp_genemodule$Module_defined %in% Module_defined[i]]
  enrich.go <- enrichGO(gene = gene,  #基因列表文件中的基因名称
                        OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                        keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                        ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                        pAdjustMethod = 'fdr',  #指定 p 值校正方法
                        pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                        qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                        readable = FALSE)
  enrich.go <-as.data.frame(enrich.go)  #
  enrich.go$Description<- factor(enrich.go$Description,levels = enrich.go$Description)
  enrich.go$Module <- rep(Module_defined[i],dim(enrich.go)[1])
  go.all <- rbind(enrich.go,go.all)
}
save(go.all,file = "./temp/temp_go.all.RData")

load(file ="./temp/temp_go.all.RData")
go.all$logPvalue <- -log10(go.all$pvalue)
temp_go.all_top5 <- go.all %>% group_by(Module) %>% top_n(5,logPvalue)
head(temp_go.all_top5)
library(ggplot2)
library(ggsci)
write.table(temp_go.all_top5,file = "./temp/temp_go.all_top5.txt",sep = "\t")

gg_final <- ggplot(temp_go.all_top5, aes(y=Module, x=`Description`, 
                                         fill=`logPvalue`, size=Count)) + 
  geom_point(colour="black", pch=21) + theme_light() +
  scale_fill_material("grey") + 
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_discrete(expand=expansion(mult = 0, add = 0.75)) +
  scale_size(range = c(0, 5)) +
  theme(plot.margin=unit(c(5.5, 40, 5.5, 5.5), "points"),
        axis.text.x = element_text(size=7)) 
gg_final
#temp_module_GO.pdf
gg_final <- ggplot(temp_go.all_top5, aes(y=Module, x=`Description`,
                                         fill=`logPvalue`, size=logPvalue)) +
  geom_point(pch=21) + theme_light() +
  scale_fill_gradient(low = "#F5E3DE",high = "#831A1F")  +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_discrete(expand=expansion(mult = 0, add = 0.75)) +
  scale_size(range = c(0, 5)) +
  theme(plot.margin=unit(c(5.5, 40, 5.5, 5.5), "points"),
        axis.text.x = element_text(size=7))
gg_final
#temp_module_GO_color.pdf
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

temp_go.all_top5 <- temp_go.all_top5[order(temp_go.all_top5$logPvalue, decreasing = F), ]
temp_go.all_top5$Description <- paste(temp_go.all_top5$Module,temp_go.all_top5$Description,sep="_")


temp_go.all_top5$Description <- factor(temp_go.all_top5$Description,levels = temp_go.all_top5$Description)
p1 <- ggplot(data = temp_go.all_top5, aes(x = -log10(pvalue), y = Description, fill = -log10(pvalue))) +
  scale_fill_gradient(low = "#F5E3DE",high = "#831A1F") +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.7) +
  geom_col(colour="black",just = 0.5,width = 0.6)+
  labs(x = "-log10(pvalue)", y = "", title = "Hippo enrichment") +
  geom_text(aes(x = 0.03, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic()
p1
#temp_module_GO_barplot.pdf


###hippo_hotspot_module_GO-----------
unique(hippo_genemodule$Module_defined)
hippo_genemodule$Module_defined <- factor(hippo_genemodule$Module_defined,
                                          levels = c("M1","M2","M3","M4","M5","M6" ))
Module_defined <- names(table(hippo_genemodule$Module_defined))
library(clusterProfiler)
library(org.Hs.eg.db)

hippo_genemodule


go.all <- data.frame(matrix(ncol = 5, nrow = 0))
for (i in 1:length(Module_defined)) {
  gene <- hippo_genemodule$Gene[hippo_genemodule$Module_defined %in% Module_defined[i]]
  enrich.go <- enrichGO(gene = gene,  #基因列表文件中的基因名称
                        OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                        keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                        ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                        pAdjustMethod = 'fdr',  #指定 p 值校正方法
                        pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                        qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                        readable = FALSE)
  enrich.go <-as.data.frame(enrich.go)  #
  enrich.go$Description<- factor(enrich.go$Description,levels = enrich.go$Description)
  enrich.go$Module <- rep(Module_defined[i],dim(enrich.go)[1])
  go.all <- rbind(enrich.go,go.all)
}
save(go.all,file = "./hippo/hippo_go.all.RData")

load(file = "./hippo/hippo_go.all.RData")
go.all$logPvalue <- -log10(go.all$pvalue)
hippo_go.all_top5 <- go.all %>% group_by(Module) %>% top_n(5,logPvalue)
head(hippo_go.all_top5)
library(ggplot2)
library(ggsci)
save(hippo_go.all_top5,file = "./hippo/hippo_go.all_top5.RData")
write.table(hippo_go.all_top5,file = "./hippo/hippo_go.all_top5.csv",sep = "\t")


load(file = "./hippo/hippo_go.all_top5.RData")

hippo_go.all_top5 <- hippo_go.all_top5[order(hippo_go.all_top5$p.adjust, decreasing = T), ]
hippo_go.all_top5$Description <- factor(hippo_go.all_top5$Description,levels = hippo_go.all_top5$Description)
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



gg_final <- ggplot(hippo_go.all_top5, aes(y=Module, x=`Description`,
                                          fill=`logPvalue`, size=logPvalue)) +
  geom_point(pch=21) + theme_light() +
  scale_fill_gradient(low = "#F5E3DE",high = "#831A1F")  +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_discrete(expand=expansion(mult = 0, add = 0.75)) +
  scale_size(range = c(0, 5)) +
  theme(plot.margin=unit(c(5.5, 40, 5.5, 5.5), "points"),
        axis.text.x = element_text(size=7))
gg_final
# #hippo_module_GO_color.pdf

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

hippo_go.all_top5 <- hippo_go.all_top5[order(hippo_go.all_top5$logPvalue, decreasing = F), ]
hippo_go.all_top5$Description <- paste(hippo_go.all_top5$Module,hippo_go.all_top5$Description,sep="_")


hippo_go.all_top5$Description <- factor(hippo_go.all_top5$Description,levels = hippo_go.all_top5$Description)
p1 <- ggplot(data = hippo_go.all_top5, aes(x = -log10(pvalue), y = Description, fill = -log10(pvalue))) +
  scale_fill_gradient(low = "#F5E3DE",high = "#831A1F") +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.7) +
  geom_col(colour="black",just = 0.5,width = 0.6)+
  labs(x = "-log10(pvalue)", y = "", title = "Hippo enrichment") +
  geom_text(aes(x = 0.03, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic()
p1
#hippo_module_GO_barplot.pdf

