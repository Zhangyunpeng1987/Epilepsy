library(devtools)
library(lisi)
library(kBET)
library(cluster)
library(factoextra)
library(dplyr)
library(ggplot2)
library(ggpubr)#
###数据准备###
set.seed(123) # 确保可重复性
#颞叶皮层----------------
load(file = "./Epilepsy19_scRNA_anno_cluster.RData")
umap_emb <- Embeddings(Epilepsy19,"umap")
### 批次混合评估 ####
#kBET评估 #
kbet_result <- kBET(
  df = umap_emb,
  batch = Epilepsy19$idents,
  plot = FALSE,
  heuristic = TRUE, # 内存优化
  k0 = 100,# 增加邻居数
  n_repeat = 10 # 增加重复次数
)
kbet_acceptance <- round(kbet_result$stats$kBET.observed*100, 1)
cat(paste("kBET acceptance rate:", kbet_acceptance,"%\n"))
# 计算整体kBET接受率（取各批次平均值）
overall_kbet_acceptance <- round(mean(kbet_result$stats$kBET.observed) * 100, 1)
# 计算整体kBET接受率（取各批次平均值）
# 计算kBET效应大小
kbet_effect_size <- (mean(kbet_result$stats$kBET.observed) -
                       kbet_result$stats$kBET.expected[1]) /
  sd(kbet_result$stats$kBET.observed)
# kBET结果显著性检验
kbet_pvalue <- t.test(
  x = kbet_result$stats$kBET.observed,
  mu = kbet_result$stats$kBET.expected[1],
  alternative = "greater")$p.value
# kBET置信区间
kbet_ci <- t.test(kbet_result$stats$kBET.observed)$conf.int
save(kbet_result,file = "./kbet_result_Temp.RData")

### 可视化 ###
# kBET图 #
# 提取kBET统计结果
kbet_df <- as.data.frame(kbet_result$stats) %>%
  mutate(category = rownames(.))  #保留批次名称
p <- ggplot(kbet_df, aes(x = category, y = kBET.observed)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = paste0("kBET Batch Acceptance (Overall: ", kbet_acceptance, "%)"),
       x = "Batch",
       y = "Acceptance Rate") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    panel.grid.major.x = element_blank() ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  geom_text(aes(label = round(kBET.observed, 2)),  # 添加文本标签
            vjust = -0.5, size = 3.5, color = "darkblue")
p
ggsave("kBET_validation_Temp.pdf", plot = p, width = 10, height = 6) #保存图形
