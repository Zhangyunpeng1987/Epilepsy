# sccoda 数据准备=========================
Epilepsy19_cell_counts <- as.data.frame(table(Epilepsy19@meta.data$idents, Epilepsy19@meta.data$Anno_Idents))
library(dplyr)
library(tidyr)
#pandas.Dataframe：第一列为样本名，其余每列各代表一种细胞类型，值表示细胞数量
Epilepsy19_cell_counts <- Epilepsy19@meta.data %>%
  count(idents, Anno_Idents, name = "n") %>%
  pivot_wider(
    names_from = Anno_Idents,
    values_from = n,
    values_fill = 0
  ) %>%
  rename(sample = idents)
write.table(Epilepsy19_cell_counts,file = "./Epilepsy19_cell_counts.txt",sep = "\t",quote = F)
#接着是python运行
