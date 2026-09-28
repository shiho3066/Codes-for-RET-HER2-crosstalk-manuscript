######### boxplot ##############
library(ggpubr)
######### whole dataset ##########
setdiff(colnames(gsva), rownames(clinical_her2))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical_her2), colnames(gsva))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva), rownames(clinical_her2))

# 重新筛选数据集
gsva_common <- gsva[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical_her2[common_samples, ]  # 筛选临床数据

# RET expression by HER2 status
plot_data <- data.frame(
  Expression = as.numeric(gsva_common["RET", ]),  # RET expression
  HER2_Status = clinical_common$HER2  # HER2 status
)

ggplot(plot_data, aes(x = HER2_Status, y = Expression, fill = HER2_Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("Negative", "Positive")), method = "t.test") +  # 只比较这两组
  labs(title = "TCGA: RET expression by HER2 Status",
       x = "HER2 Status",
       y = "RET expression") +
  theme_minimal()

# RET signature by HER2 status
plot_data <- data.frame(
  Expression = as.numeric(gsva_common["GDNF", ]),  # RET signature
  HER2_Status = clinical_common$HER2  # HER2 status
)

ggplot(plot_data, aes(x = HER2_Status, y = Expression, fill = HER2_Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("Negative", "Positive")), method = "t.test") +  # 只比较这两组
  labs(title = "TCGA: RET signature by HER2 Status",
       x = "HER2 Status",
       y = "RET signature") +
  theme_minimal()

######### ERpos ##########
setdiff(colnames(gsva_erpos), rownames(clinical_her2))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical_her2), colnames(gsva_erpos))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples_erpos <- intersect(colnames(gsva_erpos), rownames(clinical_her2))

# 重新筛选数据集
gsva_erpos_common <- gsva_erpos[, common_samples_erpos]  # 筛选 RNA-seq 数据
clinical_erpos_common <- clinical_her2[common_samples_erpos, ]  # 筛选临床数据

# RET expression by HER2 status
plot_data_erpos <- data.frame(
  Expression = as.numeric(gsva_erpos_common["RET", ]),  # RET expression
  HER2_Status = clinical_erpos_common$HER2  # HER2 状态
)

ggplot(plot_data_erpos, aes(x = HER2_Status, y = Expression, fill = HER2_Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("Negative", "Positive")), method = "t.test") +  # 只比较这两组
  labs(title = "TCGA_ERpos: RET expression by HER2 Status",
       x = "HER2 Status",
       y = "RET expression") +
  theme_minimal()

# RET signature by HER2 status
plot_data_erpos <- data.frame(
  Expression = as.numeric(gsva_erpos_common["GDNF", ]),  # RET signature
  HER2_Status = clinical_erpos_common$HER2  # HER2 状态
)

ggplot(plot_data_erpos, aes(x = HER2_Status, y = Expression, fill = HER2_Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("Negative", "Positive")), method = "t.test") +  # 只比较这两组
  labs(title = "TCGA_ERpos: RET signature by HER2 Status",
       x = "HER2 Status",
       y = "RET signature") +
  theme_minimal()

######### ERneg ##########
setdiff(colnames(gsva_erneg), rownames(clinical_her2))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical_her2), colnames(gsva_erneg))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples_erneg <- intersect(colnames(gsva_erneg), rownames(clinical_her2))

# 重新筛选数据集
gsva_erneg_common <- gsva_erneg[, common_samples_erneg]  # 筛选 RNA-seq 数据
clinical_erneg_common <- clinical_her2[common_samples_erneg, ]  # 筛选临床数据

# RET expression by HER2 status
plot_data_erneg <- data.frame(
  Expression = as.numeric(gsva_erneg_common["RET", ]),  # RET expression
  HER2_Status = clinical_erneg_common$HER2  # HER2 status
)

ggplot(plot_data_erneg, aes(x = HER2_Status, y = Expression, fill = HER2_Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("Negative", "Positive")), method = "t.test") +  # 只比较这两组
  labs(title = "TCGA_ERneg: RET expression by HER2 Status",
       x = "HER2 Status",
       y = "RET expression") +
  theme_minimal()

# RET signature by HER2 status
plot_data_erneg <- data.frame(
  Expression = as.numeric(gsva_erneg_common["GDNF", ]),  # RET signature
  HER2_Status = clinical_erneg_common$HER2  # HER2 status
)

ggplot(plot_data_erneg, aes(x = HER2_Status, y = Expression, fill = HER2_Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("Negative", "Positive")), method = "t.test") +  # 只比较这两组
  labs(title = "TCGA_ERneg: RET signature by HER2 Status",
       x = "HER2 Status",
       y = "RET signature") +
  theme_minimal()
