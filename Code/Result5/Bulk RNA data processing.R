##BULK RNA------------
# 读取数据
#GSE256068--------------
library(data.table)
GSE256068_raw_data <- fread("GSE256068_raw_data.csv.gz")
class(GSE256068_raw_data)
GSE256068_raw_data <- as.data.frame(GSE256068_raw_data)
library(stringr)
GSE256068_raw_data[1:4,1:4]
library(org.Hs.eg.db)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = GSE256068_raw_data$V1,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
GSE256068_raw_data$gene_symbol <- gene_symbols[GSE256068_raw_data$V1]
GSE256068_raw_data[1:4,1:4]
GSE256068_raw_data[1:4,(ncol(GSE256068_raw_data)-4):ncol(GSE256068_raw_data)]
length(GSE256068_raw_data$gene_symbol)
length(unique(GSE256068_raw_data$gene_symbol))

GSE256068_raw_data<-na.omit(GSE256068_raw_data)
table(duplicated(GSE256068_raw_data$gene_symbol))
library(limma)
GSE256068_raw_data[1:4,c(1,(ncol(GSE256068_raw_data)-4):ncol(GSE256068_raw_data))]
exp_unique<-avereps(GSE256068_raw_data[,-c(1,ncol(GSE256068_raw_data))],ID=GSE256068_raw_data$gene_symbol)
exp_unique[1:4,1:4]
dim(exp_unique) 
exp_unique[1:4,c(1,(ncol(exp_unique)-4):ncol(exp_unique))]

GSE256068_exp <- exp_unique
save(GSE256068_exp,file = "./GSE256068_exp.RData")

SampleFile <- fread("./GSE256068_SampleFile.txt")
unique(SampleFile$Tissue)
SampleFile <- SampleFile[SampleFile$Tissue %in% c("Temporal","Hippocampus"),]  #96*5
table(SampleFile$Tissue)
# Hippocampus    Temporal 
# 13          83 
colnames(SampleFile)[4] <- "Disease"
table(SampleFile$Disease)
table(SampleFile$Tissue,SampleFile$Disease)


#因子型
SampleFile <- SampleFile[SampleFile$Sample_ID %in% colnames(GSE256068_exp),]
SampleFile <- SampleFile[order(SampleFile$Disease),]
table(SampleFile$Disease)
save(SampleFile,file = "./SampleFile_GSE256068.RData")

GSE256068_exp[1:4,1:4]
dim(GSE256068_exp) #20714   162
#sample <- intersect(SampleFile$Sample,colnames(GSE256068_exp))
GSE256068_exp <- GSE256068_exp[,SampleFile$Sample_ID]
identical(colnames(GSE256068_exp),SampleFile$Sample_ID)
colnames(GSE256068_exp) <- SampleFile$Sample
boxplot(GSE256068_exp,outline=FALSE, notch=T,las=2)
colnames(GSE256068_exp)
#GSE256068_exp_batch_boxplot.pdf
range(GSE256068_exp)# -5.933978 14.931358 在20以内的范围就是已经log2了
dim(GSE256068_exp) #20714    96
save(GSE256068_exp,file = "./GSE256068_exp.RData")

# #1.差异表达分析--------------
#limma-------------
#构建分组矩阵--design
load(file = "./GSE256068_exp.RData")
load(file = "./SampleFile_GSE256068.RData")
# table(SampleFile$Disease)
# # ControlCortex ControlHippocampus              FCD2a              FCD2b             TLE_HS                TSC 
# # 4                 13                  3                  6                 64                  6 
# SampleFile$Group[SampleFile$Disease %in% c("ControlCortex","ControlHippocampus")] <- "Control"
# SampleFile$Group[SampleFile$Disease %in% c("FCD2a" ,"FCD2b","TLE_HS" ,"TSC")] <- "Epilepsy"
table(SampleFile$Group)
# Control Epilepsy 
# 17       79
#save(SampleFile,file = "./SampleFile_GSE256068.RData")

design <- model.matrix(~0+factor(SampleFile$Group))
colnames(design) <- levels(factor(SampleFile$Group))

rownames(design) <- colnames(GSE256068_exp)
#构建比较矩阵——contrast
contrast.matrix <- makeContrasts(Epilepsy-Control,levels = design)
#limma DEG
fit <- lmFit(GSE256068_exp,design)##线性拟合模型构建
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2, coef = 1,n = Inf)

DEG$type <- ifelse(DEG$adj.P.Val > 0.05, "no_change",
                   ifelse(DEG$logFC > 1, "up",
                          ifelse(DEG$logFC < -1, "down", "no_change")))
DEG  <- dplyr::filter(DEG,  !is.na(DEG$type))

save(DEG,file = "./DEGs.RData")
dif <- DEG[DEG$adj.P.Val<0.05&abs(DEG$logFC)>1,]
dif <- dif[order(dif$logFC),]
save(dif,file = "./DEGs_filter.RData")
dim(dif)
##1916    7
##火山图-------------
library(ggplot2)
load(file = "./DEGs.RData")
head(DEG)
table(DEG$type)
# down no_change        up 
# 927     18798       989 
GSE256068_exp[1:4,1:4]
ggplot(DEG, aes(logFC, -log10(adj.P.Val)))+   #读取差异表达结果，X轴为logFC,y轴为-log10(P.Value)
  
  geom_point(aes(col=type))+    #设置点数据的来源
  
  scale_color_manual(values=c("#0072B5","grey","#BC3C28"))+  #这里可以对途中，上调，不变和下调的点的颜色进行设置
  
  labs(x="log2(FoldChange)",y="-log10(adj.P.Val)")+  #设置x轴和y轴的标签
  
  geom_vline(xintercept=c(-0.5,0.5), colour="grey", linetype="dashed")+  #设置x轴的分界线
  
  geom_hline(yintercept = -log10(0.05),colour="grey", linetype="dashed") #设置y轴的分界线

ggsave("./diff_gene_volcano.pdf",height = 5,width = 5)

#2.功能富集分析------------
library(clusterProfiler)
library(org.Hs.eg.db)
##GO------------
load(file = "./DEGs_filter.RData")
head(dif) #
##全部---------------
enrich.go <- enrichGO(gene = rownames(dif),  #基因列表文件中的基因名称
                      OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                      keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                      qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                      readable = FALSE)
enrich.go <-as.data.frame(enrich.go)  #794
save(enrich.go,file = "./DEGs_GO_enrichment.RData")
load(file = "./DEGs_GO_enrichment.RData")
##可视化---------
dat <- enrich.go[1:10,]
dat <- dat[order(dat$Count,decreasing = F),]
dat$Description <- factor(dat$Description, levels = dat$Description)

#柱形图，横坐标 p 值的对数转换，纵坐标是 GO Term，颜色按 Category 着色
p2 <- ggplot(dat, aes(Description, Count)) +
  geom_col(aes(fill = -log10(pvalue)), width = 0.5) +
  scale_fill_gradient(low = "#F5E3DE",high = "#831A1F") +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  coord_flip() +
  labs(y = 'Count',x = "", title = "DEGs GO enrichment")+
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14,
                              hjust = 0.5,vjust = 0.5,
                              face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )
p2
#上调------------------
dif$logFC
up_genes <- dif %>%
  filter(adj.P.Val < 0.05 & logFC > 1) %>%
  rownames()  #
down_genes <- dif %>%
  filter(adj.P.Val < 0.05 & logFC < -1) %>%
  rownames()  #

enrich.go <- enrichGO(gene = up_genes,  #基因列表文件中的基因名称
                      OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                      keyType = 'SYMBOL',  
                      ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                      qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                      readable = FALSE)
enrich.go <-as.data.frame(enrich.go)  #826
save(enrich.go,file = "./up_genes_go.RData")
load(file = "./up_genes_go.RData")
##可视化---------
dat <- enrich.go[1:10,]
dat <- dat[order(dat$Count,decreasing = F),]
dat$Description <- factor(dat$Description, levels = dat$Description)
p2 <- ggplot(dat, aes(Description, Count)) +
  geom_col(aes(fill = -log10(pvalue)), width = 0.5) +
  scale_fill_gradient(low = "#F5E3DE",high = "#831A1F") +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  coord_flip() +
  labs(y = 'Count',x = "", title = "DEGs GO enrichment")+
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14,
                              hjust = 0.5,vjust = 0.5,
                              face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )
p2

#下调------------------
enrich.go <- enrichGO(gene = down_genes,  #基因列表文件中的基因名称
                      OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                      keyType = 'SYMBOL',  
                      ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                      qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                      readable = FALSE)
enrich.go <-as.data.frame(enrich.go)  #94
save(enrich.go,file = "./down_genes_go.RData")
load(file = "./down_genes_go.RData")
##可视化---------
dat <- enrich.go[1:10,]
dat <- dat[order(dat$Count,decreasing = F),]
dat$Description <- factor(dat$Description, levels = dat$Description)

p2 <- ggplot(dat, aes(Description, Count)) +
  geom_col(aes(fill = -log10(pvalue)), width = 0.5) +
  scale_fill_gradient(low = "#F5E3DE",high = "#831A1F") +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  coord_flip() +
  labs(y = 'Count',x = "", title = "DEGs GO enrichment")+
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14,
                              hjust = 0.5,vjust = 0.5,
                              face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )
p2

##KEGG------------
load(file = "./DEGs_filter.RData")
head(dif) #
dif$gene <- rownames(dif)
genelist <- bitr(dif$gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
library(dplyr)
#inner_join() 函数要基于 DEG 数据框的 "Gene" 列和 genelist 数据框的 "SYMBOL" 列进行连接。
dif <- inner_join(dif,genelist,by=c("gene"="SYMBOL"))
enrich.kegg <- enrichKEGG(gene = dif$ENTREZID,  #基因列表文件中的基因名称
                          organism = 'hsa',  #指定物种的基因数据库
                          pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                          qvalueCutoff = 0.2  #指定 q 值阈值，不显著的值将不显示在结果中
)
enrich.kegg <-as.data.frame(enrich.kegg)  #45
save(enrich.kegg,file = "./dif_gene_kegg.RData")

##可视化---------
dat <- enrich.kegg[1:10,]
dat <- dat[order(dat$Count,decreasing = F),]
dat$Description <- factor(dat$Description, levels = dat$Description)
p2 <- ggplot(dat, aes(Description, Count)) +
  geom_col(aes(fill = -log10(pvalue)), width = 0.5) +
  scale_fill_gradient(low = "#F5E3DE",high = "#831A1F") +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  coord_flip() +
  labs(y = 'Count',x = "", title = "DEGs KEGG enrichment")+
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14,
                              hjust = 0.5,vjust = 0.5,
                              face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )
p2
##与单细胞神经元激活相关的转录因子和靶基因取交集-----------
#GSE256068--------------------
load(file = "./DEGs_filter.RData")
dim(dif) #

library(readr)
Neuronal_activation_subnetwork <- read_delim("./Neuronal_activation_subnetwork.txt", 
                                             delim = "\t", escape_double = FALSE, 
                                             trim_ws = TRUE)
TF <- unique(c(Neuronal_activation_subnetwork$TF,Neuronal_activation_subnetwork$TargetGene))
Neuronal_activation_subnetwork_gene <- intersect(unique(c(Neuronal_activation_subnetwork$TF,Neuronal_activation_subnetwork$TargetGene)),
                                                 rownames(dif)
)

library(VennDiagram)
library(grid)
# 创建韦恩图对象
venn.plot <- draw.pairwise.venn(
  area1 = length(rownames(dif)),       # 第一个集合大小
  area2 = length(TF),            # 第二个集合大小
  cross.area = length(Neuronal_activation_subnetwork_gene),              # 交集大小
  category = c("DEGs", "Neuronal activation gene"), # 集合名称
  fill = c("#c86f5e", "skyblue"),           # 颜色填充
  alpha = 0.7,                           # 透明度
  cat.pos = c(0, 0),                     # 标签位置
  cat.dist = 0.05,                       # 标签距离
  ext.text = FALSE
)

# 显示图形
grid.newpage()
grid.draw(venn.plot)
# [1] "FOS"     "EGR3"    "EGR1"    "HAVCR2"  "CD74"    "PTPRC"   "CCL3"    "CX3CR1"  "IL1B"    "CCL4"    "SOCS6"  
# [12] "CCL2"    "HOMER1"  "SLC7A11" "BTG2"    "PTGS2"   "COL4A1"  "EGR2"    "FOSB"    "TPPP3"   "GFAP"    "JUNB"   
# [23] "ZFP36"   "HSPA1A" 
gene_expr <- GSE256068_exp[rownames(GSE256068_exp) %in% Neuronal_activation_subnetwork_gene, ]      # 选出目标基因表达
gene_expr <- t(gene_expr)  # 转置为样本×基因

# group 向量，样本分组，长度与样本数一致
group <- SampleFile$Group  # 0=normal, 1=epilepsy
group[group %in% "Control"] <- "0"
group[group %in% "Epilepsy"] <- "1"

group <- as.factor(group)  # 0=normal, 1=epilepsy
# LASSO回归筛选特征基因---------
library(glmnet)
# 训练LASSO模型
set.seed(123)
cvfit <- cv.glmnet(
  as.matrix(gene_expr),
  as.numeric(as.character(group)),   # 注意要转为0/1数值
  family = "binomial",               # 二分类
  alpha = 1                          # LASSO
)
plot(cvfit$glmnet.fit, xvar="lambda", label=TRUE, main="LASSO coefficient paths")
#lasso_lujing.pdf  8*8
plot(cvfit, main="LASSO cross-validation curve")
abline(v=log(cvfit$lambda.min), col="red", lty=2)
#lasso_CV.pdf  8*8

# 提取最佳lambda时的基因
coef_lasso <- coef(cvfit, s = "lambda.min")
lasso_genes <- rownames(coef_lasso)[which(coef_lasso != 0)][-1]
print(lasso_genes)
#"IL1B"    "GFAP"    "SLC7A11" "SOCS6"   "COL4A1"  "HSPA1A" 
save(lasso_genes,file = "./lasso_genes.RData")

#SVM递归特征消除（RFE）----------
load(file = "./GSE256068_exp.RData")
load(file = "./SampleFile_GSE256068.RData")

library(caret)
library(e1071)
gene_expr <- GSE256068_exp[rownames(GSE256068_exp) %in% Neuronal_activation_subnetwork_gene, ]  
gene_expr <- t(gene_expr)  

dat <- data.frame(gene_expr)
dat$group <- group
library(caret)
library(e1071)

x <- dat[, !(names(dat) %in% "group")]
y <- as.factor(dat$group)


stopifnot(!any(is.na(x)))
stopifnot(!any(is.na(y)))
set.seed(123)
ctrl <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
svmProfile <- rfe(
  x,
  y,
  sizes = c(1:10, 15, 20),
  rfeControl = ctrl,
  method = "rf"
)
plot(svmProfile, type = c("g", "o"), main = "RFE Accuracy vs. Number of Features")
print(svmProfile)
svmProfile$optVariables
predictors(svmProfile)
#"CCL3"   "IL1B"   "CX3CR1" "CCL4"   "COL4A1" "EGR2"   "SOCS6"  "FOSB"
varImp(svmProfile)
svm_genes <- predictors(svmProfile)
print(svm_genes)
save(svm_genes,file = "./svm_genes.RData")

load("./lasso_genes.RData")
load(file = "./svm_genes.RData")
common_genes <- intersect(lasso_genes, svm_genes)
print(common_genes)
#"IL1B"   "SOCS6"  "COL4A1"


library(VennDiagram)
library(grid)
# 创建韦恩图对象
venn.plot <- draw.pairwise.venn(
  area1 = length(lasso_genes),       # 第一个集合大小
  area2 = length(svm_genes),            # 第二个集合大小
  cross.area = length(common_genes),              # 交集大小
  category = c("Lasso genes", "RF genes"), # 集合名称
  fill = c("#c86f5e", "skyblue"),           # 颜色填充
  alpha = 0.7,                           # 透明度
  cat.pos = c(0, 0),                     # 标签位置
  cat.dist = 0.05,                       # 标签距离
  ext.text = FALSE
)

# 显示图形
grid.newpage()
grid.draw(venn.plot)


library(caret)
library(pROC)

X <- t(GSE256068_exp[common_genes, ])
X <- as.data.frame(X)
y <- as.factor(group) 

stopifnot(nrow(X) == length(y))

model <- glm(y ~ ., data = X, family = binomial)

summary(model)


pred_prob <- predict(model, type = "response")

roc_obj <- roc(y, pred_prob)
auc_score <- auc(roc_obj)
cat("模型AUC:", auc_score, "\n")
plot(roc_obj, col="red", lwd=2,main=paste("Logistic Regression ROC, AUC =", round(auc_score, 3)))
coef(summary(model))

library(pROC)
library(RColorBrewer)

y <- as.factor(group)
if (is.character(y) || is.factor(y)) y <- as.numeric(as.character(y))

cols <- brewer.pal(min(6, length(common_genes)), "Set1")

plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab="1-Specificity", ylab="Sensitivity", main="Single Gene ROC Curves")
abline(0,1,lty=2,col="gray")

legend_text <- c()
cols <- c("#99CCCC","#336699","#996699")
for (i in seq_along(common_genes)) {
  gene <- common_genes[i]
  gene_expr <- as.numeric(GSE256068_exp[gene, ])
  roc_obj <- roc(y, gene_expr, quiet=TRUE)
  lines(1 - roc_obj$specificities, roc_obj$sensitivities, col=cols[i], lwd=2)
  legend_text <- c(legend_text, paste0(gene, " (AUC=", round(auc(roc_obj), 3), ")"))
}

legend("bottomright", legend=legend_text, col=cols, lwd=2, cex=0.9)

##GSE139914  8癫痫  39正常--------------
library(data.table)
GSE139914_raw_data <- fread("GSE139914_Within_Subject_RawCounts.txt.gz")

GSE139914_exp<-na.omit(GSE139914_raw_data)
table(duplicated(GSE139914_exp$V1))
library(limma)
exp_unique<-avereps(GSE139914_exp[,-1],ID=GSE139914_exp$V1)
exp_unique[1:4,1:4]
dim(exp_unique) 
exp_unique[1:4,c(1,(ncol(exp_unique)-4):ncol(exp_unique))]

GSE139914_exp <- exp_unique
range(GSE139914_exp) #0 798857
GSE139914_exp <- log2(GSE139914_exp+1) #0.00000 19.60758
save(GSE139914_exp,file = "./GSE139914_exp.RData")

load(file = "./GSE139914_exp.RData")
GSE139914_SampleFile <- fread("GSE139914_SampleFile.txt")

group_new <- GSE139914_SampleFile$Group
group_new <- factor(group_new, levels = c("Control", "Epilepsy"))

library(pROC)
gene_list <- c("IL1B","SOCS6","COL4A1") # 

plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab="1-Specificity", ylab="Sensitivity", main="Single Gene ROC")
abline(0,1,lty=2,col="gray")
cols <- c("#99CCCC","#336699","#996699")[1:length(gene_list)]
legend_text <- c()

for(i in seq_along(gene_list)){
  gene <- gene_list[i]
  gene_expr <- as.numeric(GSE139914_exp[gene, ])
  roc_obj <- roc(group_new, gene_expr, quiet=TRUE, levels=rev(levels(group_new)))
  lines(1-roc_obj$specificities, roc_obj$sensitivities, col=cols[i], lwd=2)
  legend_text <- c(legend_text, paste0(gene, " (AUC=", round(auc(roc_obj),3), ")"))
}
legend("bottomright", legend=legend_text, col=cols, lwd=2)

X_new <- t(GSE139914_exp[gene_list, ]) 
X_new <- as.data.frame(X_new)

model_new <- glm(group_new ~ ., data=X_new, family=binomial)
pred_prob_new <- predict(model_new, type="response")
roc_model <- roc(group_new, pred_prob_new, levels=rev(levels(group_new)))
auc_model <- auc(roc_model)
plot(roc_model, col="red", lwd=2, main=paste("3-Gene Model ROC, AUC =", round(auc_model, 3)))

#GSE140393  9正常  12癫痫----------------
library(data.table)
GSE140393_raw_data <- fread("GSE140393_rldnormalized_EGFR.txt.gz")
GSE140393_raw_data <- GSE140393_raw_data[c(1,86:nrow(GSE140393_raw_data))]
colnames(GSE140393_raw_data)[2:4] <- as.character(GSE140393_raw_data[1,2:4])
GSE140393_raw_data <- GSE140393_raw_data[-1,-c(1,6,7)]
GSE140393_exp<-na.omit(GSE140393_raw_data)
table(duplicated(GSE140393_exp$hgnc_symbol))
GSE140393_exp <- GSE140393_exp[!(is.na(GSE140393_exp$hgnc_symbol) & GSE140393_exp$hgnc_symbol == ""), ]

row_data <- GSE140393_exp[ , -4]
non_zero_row <- rowSums(row_data != 0) > 0
GSE140393_exp <- GSE140393_exp[non_zero_row, ]

dim(GSE140393_exp)
table(duplicated(GSE140393_exp$hgnc_symbol))

GSE140393_exp <-  as.data.frame(GSE140393_exp)
library(limma)
gene_ids <- as.character(GSE140393_exp$hgnc_symbol)
exp_unique <- avereps(GSE140393_exp[,-4], ID = gene_ids)
if (is.null(rownames(exp_unique))) {
  rownames(exp_unique) <- unique(gene_ids)
}
rownames(exp_unique)
exp_unique[1:10,1:3]
exp_unique <- as.matrix(exp_unique)
dim(exp_unique) #20726     3

GSE140393_exp <- exp_unique
range(GSE140393_exp) #"-0.001887921" "9.999648117" 
save(GSE140393_exp,file = "./GSE140393_exp.RData")

GSE140393_raw_data_1 <- fread("GSE140393_rldnormalized_nuclei.txt.gz")
GSE140393_raw_data <- GSE140393_raw_data_1[c(1,86:nrow(GSE140393_raw_data_1))]
colnames(GSE140393_raw_data)[2:19] <- as.character(GSE140393_raw_data[1,2:19])
GSE140393_raw_data <- GSE140393_raw_data[-1,-c(1,21,22)]
GSE140393_exp<-na.omit(GSE140393_raw_data)
table(duplicated(GSE140393_exp$hgnc_symbol))

GSE140393_exp <- GSE140393_exp[!(is.na(GSE140393_exp$hgnc_symbol) & GSE140393_exp$hgnc_symbol == ""), ]

row_data <- GSE140393_exp[ , -19]
non_zero_row <- rowSums(row_data != 0) > 0
GSE140393_exp <- GSE140393_exp[non_zero_row, ]

dim(GSE140393_exp)
table(duplicated(GSE140393_exp$hgnc_symbol))
library(limma)
gene_ids <- as.character(GSE140393_exp$hgnc_symbol)
exp_unique <- avereps(GSE140393_exp[,-19], ID = gene_ids)
if (is.null(rownames(exp_unique))) {
  rownames(exp_unique) <- unique(gene_ids)
}
rownames(exp_unique)
exp_unique[1:10,1:3]
dim(exp_unique) 

GSE140393_exp <- exp_unique
range(GSE140393_exp) #"-0.000123268" "9.999995611"
GSE140393_exp_1 <- GSE140393_exp
save(GSE140393_exp_1,file = "./GSE140393_exp_1.RData")


load("./GSE140393_exp.RData")
load("./GSE140393_exp_1.RData")
GSE140393_exp <- as.data.frame(GSE140393_exp)
GSE140393_exp_1 <- as.data.frame(GSE140393_exp_1)

gene <- intersect(rownames(GSE140393_exp),rownames(GSE140393_exp_1))
GSE140393_exp <- GSE140393_exp[gene,]
GSE140393_exp_1 <- GSE140393_exp_1[gene,]

GSE140393_exp <- cbind(GSE140393_exp,GSE140393_exp_1)
GSE140393_exp_1 <- apply(GSE140393_exp, 2, as.numeric)
rownames(GSE140393_exp_1) <- rownames(GSE140393_exp) 

SampleFile <- fread("GSE140393_SampleFile.txt")
SampleFile$Group
GSE140393_exp <-GSE140393_exp_1[,SampleFile$Patient]
colnames(GSE140393_exp) <- SampleFile$Sample

boxplot(GSE140393_exp,outline=FALSE, notch=T,las=2)
#
GSE140393_exp<-na.omit(GSE140393_exp)

range(GSE140393_exp)#-2.359978 21.928633 


save(GSE140393_exp,file = "./GSE140393_exp_final.RData")

load(file = "./GSE140393_exp_final.RData")
gene_list <- c("IL1B","SOCS6","COL4A1") # 

GSE140393_SampleFile <- fread("GSE140393_SampleFile.txt")

group_new <- GSE140393_SampleFile$Group
group_new <- factor(group_new, levels = c("Control", "Epilepsy"))

library(pROC)

plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab="1-Specificity", ylab="Sensitivity", main="Single Gene ROC")
abline(0,1,lty=2,col="gray")
cols <- c("#99CCCC","#336699","#996699")[1:length(gene_list)]
legend_text <- c()

for(i in seq_along(gene_list)){
  gene <- gene_list[i]
  gene_expr <- as.numeric(GSE140393_exp[gene, ])
  roc_obj <- roc(group_new, gene_expr, quiet=TRUE, levels=rev(levels(group_new)))
  lines(1-roc_obj$specificities, roc_obj$sensitivities, col=cols[i], lwd=2)
  legend_text <- c(legend_text, paste0(gene, " (AUC=", round(auc(roc_obj),3), ")"))
}
legend("bottomright", legend=legend_text, col=cols, lwd=2)

X_new <- t(GSE140393_exp[gene_list, ]) 
X_new <- as.data.frame(X_new)

model_new <- glm(group_new ~ ., data=X_new, family=binomial)
pred_prob_new <- predict(model_new, type="response")
roc_model <- roc(group_new, pred_prob_new, levels=rev(levels(group_new)))
auc_model <- auc(roc_model)
plot(roc_model, col="red", lwd=2, main=paste("3-Gene Model ROC, AUC =", round(auc_model, 3)))
