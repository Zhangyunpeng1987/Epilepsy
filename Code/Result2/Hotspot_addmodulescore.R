###hotspots-----------
library(Seurat)
###海马体-----------
Excitatory <- readRDS(file = "./hippo/Excitatory_sub.rds")
names(table(Excitatory$Anno_Idents))
Inhibitory <- readRDS(file = "./hippo/Inhibitory_sub.rds")
Hippo_scRNA <- merge(x=Excitatory,y=Inhibitory)
#devtools::install_github("cellgeni/sceasy")
# library(sceasy)
# library(reticulate)
# Sys.setenv(RETICULATE_MINICONDA_PATH="~/miniconda")
# reticulate::use_condaenv("hotspotenv")
# sceasy::convertFormat(Hippo_scRNA, from="seurat", to="anndata", assay = "RNA", main_layer = "count",
#                       outFile='./hippo/hippo.h5ad')

# library(SeuratDisk) #采用
# SaveH5Seurat(Hippo_scRNA, filename = "./hippo/hippo.h5Seurat")
# Convert("./hippo/hippo.h5Seurat", dest = "h5ad")
#
#
# Hippo_scRNA[["RNA"]] <- as(Hippo_scRNA[["RNA"]], "Assay")
# sceasy::convertFormat(Hippo_scRNA, from="seurat", to="anndata", outFile='./hippo/hippo.h5ad')
hippo_genemodule <- read.csv("./hippo/hippo_genemodule.csv")
hippo_genemodule <- hippo_genemodule[!is.na(hippo_genemodule$Module),]
hippo_genemodule <- hippo_genemodule[order(hippo_genemodule$Module),]
hippo_genemodule <- hippo_genemodule[-which(hippo_genemodule$Module ==(-1)),]

hippo_genemodule$Module_defined[hippo_genemodule$Module %in% c(1,5,10)] <- "M1"
hippo_genemodule$Module_defined[hippo_genemodule$Module %in% c(2,8)] <- "M2"
hippo_genemodule$Module_defined[hippo_genemodule$Module %in% c(6)] <- "M3"
hippo_genemodule$Module_defined[hippo_genemodule$Module %in% c(4,7)] <- "M4"
hippo_genemodule$Module_defined[hippo_genemodule$Module %in% c(13)] <- "M5"
hippo_genemodule$Module_defined[hippo_genemodule$Module %in% c(3,9,11,12)] <- "M6"
###模块基因对细胞打分---------------
Module_defined <- names(table(hippo_genemodule$Module_defined))

for (i in 1:length(Module_defined)) {
  genelist <- list(hippo_genemodule$Gene[hippo_genemodule$Module_defined %in% Module_defined[i]])
  Hippo_scRNA <- AddModuleScore(object = Hippo_scRNA, features = genelist,name = Module_defined[i])
  
}
save(Hippo_scRNA,file = "./neuron/Hippo_scRNA_AddModuleScore.RData")

###棒棒糖图------------
df1 <- Hippo_scRNA@meta.data[,c(18,26:31)]
df1$Anno_Idents <- factor(df1$Anno_Idents,levels = c("CA_Ex_GAPDH","CA1_Ex","CA3_Ex","DG_Ex1","DG_Ex2","DG_Ex3","DG_Ex4","DG_Ex5",
                                                     "CCK","LAMP5","PVALB","SST","VIP"))
# df1$celltype[df1$Anno_Idents %in% c("CA_Ex_GAPDH","CA1_Ex","CA3_Ex","DG_Ex1","DG_Ex2","DG_Ex3","DG_Ex4","DG_Ex5")] <- "Excitatory"
# df1$celltype[df1$Anno_Idents %in% c("CCK","LAMP5","PVALB","SST","VIP")] <- "Inhibitory"


names(table(df1$Anno_Idents))
library(ggplot2)
library(ggpubr)
library(tidyr)
celltype_color <- c("CA_Ex_GAPDH"="#FB7D1A","CA1_Ex"="#FABF74","CA3_Ex"="#CAB4D6","DG_Ex1"="#AAD0E3",
                    "DG_Ex2"="#277AB4","DG_Ex3"="#B5DF90","DG_Ex4"="#693C9A","DG_Ex5"="#3AA12F",
                    "CCK"="#FFCC00","LAMP5"="#FF9900","PVALB"="#6699CC","SST"="#339999","VIP"="#CCCC99")

df1 <- melt(df1)
colnames(df1) <- c("Anno_Idents","Module","Score" )
df1 <- aggregate(x = df1$Score,by = list(df1$Anno_Idents,df1$Module),
                 FUN = mean)

colnames(df1) <- c("Anno_Idents","Module","Score" )

library(ggh4x)
library(ggtext)
ggplot(df1,aes(x = Anno_Idents, y = Score)) +
  geom_segment(aes(x = Anno_Idents, xend = Anno_Idents,
                   y = 0, yend = Score,
                   #linetype = Anno_Idents,
                   color = Anno_Idents)) + 
  geom_point(aes(color = Anno_Idents),size = 3.5) + #添加散点/气泡
  scale_fill_manual(values = celltype_color) +
  scale_color_manual(values = celltype_color)+
  # geom_hline(yintercept = 1.30103,colour="grey", linetype="dashed")+
  # geom_hline(yintercept = -1.30103,colour="grey", linetype="dashed")+
  theme_bw()+
  facet_nested(~ Module,scales = "free_x",space = "free") +
  theme(legend.title = element_markdown(color="black",face="bold"),
        legend.text = element_text(color="black",face= "bold.italic"),
        axis.ticks = element_blank(),
        axis.text = element_text(color="black"),
        axis.title.x = element_markdown(color="black",face="bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin=unit(c(2,0.3,0.3,0.3),unit="cm"),
        panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

#hippo_module_bangbangtang.pdf  5*12


###颞叶皮层-----------
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
##细胞打分---------
for (i in 1:length(Module_defined)) {
  genelist <- list(Temp_genemodule$Gene[Temp_genemodule$Module_defined %in% Module_defined[i]])
  Temp_scRNA <- AddModuleScore(object = Temp_scRNA, features = genelist,name = Module_defined[i])
  
}

###打分采用----------------
colnames(Temp_scRNA@meta.data)
df1 <- Temp_scRNA@meta.data[,c(8,21:28)]
###棒棒糖图------------
df1$Anno_Idents <- factor(df1$Anno_Idents,levels = c("L2_3_Cux2","L4_L5_6_Rorb","L6",
                                                     "In_Lamp5","In_Pvalb","In_Sst","In_Vip"))
names(table(df1$Anno_Idents))
library(ggplot2)
library(ggpubr)
library(tidyr)
celltype_color <- c("L2_3_Cux2"="#C4A5DE","L4_L5_6_Rorb"="#81B8DF","L6"="#F6CAE5",
                    "In_Lamp5"="#D76364","In_Pvalb"="#8ECFC9","In_Sst"="#82B0D2","In_Vip"="#FFBE7A")

df1 <- melt(df1)
colnames(df1) <- c("Anno_Idents","Module","Score" )
df1 <- aggregate(x = df1$Score,by = list(df1$Anno_Idents,df1$Module),
                 FUN = mean)

colnames(df1) <- c("Anno_Idents","Module","Score" )

library(ggh4x)
library(ggtext)
ggplot(df1,aes(x = Anno_Idents, y = Score)) +
  geom_segment(aes(x = Anno_Idents, xend = Anno_Idents,
                   y = 0, yend = Score,
                   #linetype = Anno_Idents,
                   color = Anno_Idents)) + 
  geom_point(aes(color = Anno_Idents),size = 3.5) + #添加散点/气泡
  scale_fill_manual(values = celltype_color) +
  scale_color_manual(values = celltype_color)+
  # geom_hline(yintercept = 1.30103,colour="grey", linetype="dashed")+
  # geom_hline(yintercept = -1.30103,colour="grey", linetype="dashed")+
  theme_bw()+
  facet_nested(~ Module,scales = "free_x",space = "free") +
  theme(legend.title = element_markdown(color="black",face="bold"),
        legend.text = element_text(color="black",face= "bold.italic"),
        axis.ticks = element_blank(),
        axis.text = element_text(color="black"),
        axis.title.x = element_markdown(color="black",face="bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin=unit(c(2,0.3,0.3,0.3),unit="cm"),
        panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

#temp_module_bangbangtang.pdf  5*12











