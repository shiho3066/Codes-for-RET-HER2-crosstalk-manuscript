install.packages("GSVA")
# Load necessary libraries
library(GSVA)
library(Biobase)


# Convert the data frame to a matrix
expression_matrix <- as.matrix(rnaseq)

# Define the gene signatures
gene_signatures <- list(
  GDNF_Up = c("F2RL1", "PHLDA1", "SHC4", "EGLN3", "GAL", "CEMIP", "MALL", "PHLDA2", "ACYP1", "AHNAK2", "CENPK", "CPE", "IFIT1", "ISG15", "MTMR11", "MSMO1", "IFI27", "IFI6", "COL4A5", "DDX60", "EIF2AK2", "HERC5", "HERC6", "IFIH1", "IFIT2", "IFITM1", "IFITM3", "IRF9", "OAS1", "OAS3", "PARP9", "RAB29", "SAMD9", "SP110", "UBE2L6"),
  GDNF_Down = c("AMTN", "CIRBP", "FOSB", "H3-3B", "SGK1", "RET", "ABHD11", "ANXA2", "APOD", "ADIRF", "COLEC12", "NECAB1", "EFHD1", "EGR1", "FAM107B", "FGD3", "GCLM", "GCNT1", "HMOX1", "KLK11", "KRT15", "NECAB1", "NQO1", "PDLIM1", "RBM3", "SNORD13", "SUSD2", "VSIG2"),
  Novel = c("ERBB2", "SPDEF", "TFAP2B", "CD24", "SERHL2", "CNTNAP2", "RPL19", "CAPN13", "RPL23", "LRRC26", "PRODH", "GPRC5C", "GGCT", "CLCA2", "KDM5B", "SPP1", "PHLDA1", "C15orf48", "SUSD3", "SERPINA1"),
  Established = c("ERBB2", "POLD1", "PSMD3", "PNMT", "GSDMB", "CASC3", "LASP1", "WIPF2", "EPN3", "PHB1", "CLCA2", "ORMDL2", "RAP1GAP", "CUEDC1", "HOXC11", "CYP2J2", "HGD", "ABCA12", "ATP2C2", "ITGA3", "CEACAM5", "ANO10", "NR1D1", "SNX7", "FJX1", "KCTD9", "CDK18", "CREG1")
)

gene_signatures_old <- list(
  GDNF_Up = c("F2RL1", "PHLDA1", "SHC4", "EGLN3", "GAL", "KIAA1199", "MALL", "PHLDA2", "ACYP1", "C14orf78", "CENPK", "CPE", "IFIT1", "ISG15", "MTMR11", "SC4MOL", "IFI27", "IFI6", "COL4A5", "DDX60", "EIF2AK2", "HERC5", "HERC6", "IFIH1", "IFIT2", "IFITM1", "IFITM3", "IRF9", "OAS1", "OAS3", "PARP9", "RAB7L1", "SAMD9", "SP110", "UBE2L6"),
  GDNF_Down = c("LOC441763", "AMTN", "CIRBP", "FOSB", "H3F3B", "SGK", "RET", "ABHD11", "ANXA2", "APOD", "C10orf116", "COLEC12", "EFCBP1", "EFHD1", "EGR1", "FAM107B", "FGD3", "GCLM", "GCNT1", "HMOX1", "KLK11", "KRT15", "NECAB1", "NQO1", "PDLIM1", "RBM3", "SNORD13", "SUSD2", "VSIG2"),
  Novel = c("ERBB2", "SPDEF", "TFAP2B", "CD24", "SERHL2", "CNTNAP2", "RPL19", "CAPN13", "RPL23", "LRRC26", "PRODH", "GPRC5C", "GGCT", "CLCA2", "KDM5B", "SPP1", "PHLDA1", "C15orf48", "SUSD3", "SERPINA1"),
  Established = c("ERBB2", "POLD1", "PSMD3", "PNMT", "GSDMB", "CASC3", "LASP1", "WIPF2", "EPN3", "PHB", "CLCA2", "ORMDL2", "RAP1GAP", "CUEDC1", "HOXC11", "CYP2J2", "HGD", "ABCA12", "ATP2C2", "ITGA3", "CEACAM5", "ANO10", "NR1D1", "SNX7", "FJX1", "KCTD9", "CDK18", "CREG1")
)

gene_signatures_tcga <- list(
  GDNF_Up = c("F2RL1", "PHLDA1", "SHC4", "EGLN3", "GAL", "KIAA1199", "MALL", "PHLDA2", "ACYP1", "AHNAK2", "CENPK", "CPE", "IFIT1", "ISG15", "MTMR11", "MSMO1", "IFI27", "IFI6", "COL4A5", "DDX60", "EIF2AK2", "HERC5", "HERC6", "IFIH1", "IFIT2", "IFITM1", "IFITM3", "IRF9", "OAS1", "OAS3", "PARP9", "RAB7L1", "SAMD9", "SP110", "UBE2L6"),
  GDNF_Down = c("LOC441763", "AMTN", "CIRBP", "FOSB", "H3F3B", "SGK1", "RET", "ABHD11", "ANXA2", "APOD", "C10orf116", "COLEC12", "NECAB1", "EFHD1", "EGR1", "FAM107B", "FGD3", "GCLM", "GCNT1", "HMOX1", "KLK11", "KRT15", "NECAB1", "NQO1", "PDLIM1", "RBM3", "SNORD13", "SUSD2", "VSIG2"),
  Novel = c("ERBB2", "SPDEF", "TFAP2B", "CD24", "SERHL2", "CNTNAP2", "RPL19", "CAPN13", "RPL23", "LRRC26", "PRODH", "GPRC5C", "GGCT", "CLCA2", "KDM5B", "SPP1", "PHLDA1", "C15orf48", "SUSD3", "SERPINA1"),
  Established = c("ERBB2", "POLD1", "PSMD3", "PNMT", "GSDMB", "CASC3", "LASP1", "WIPF2", "EPN3", "PHB", "CLCA2", "ORMDL2", "RAP1GAP", "CUEDC1", "HOXC11", "CYP2J2", "HGD", "ABCA12", "ATP2C2", "ITGA3", "CEACAM5", "ANO10", "NR1D1", "SNX7", "FJX1", "KCTD9", "CDK18", "CREG1")
)
# 直接融合成一个向量
gene_list <- unlist(gene_signatures_tcga)

# 去重（如果有重复的基因）
gene_list <- unique(gene_list)

# 查看结果
print(gene_list)

# 获取 RNA-seq 表中的基因名
rnaseq_genes <- rownames(rnaseq)

# 找出 gene_list 中 **不在** rnaseq_genes 里的基因
missing_genes <- setdiff(gene_list, rnaseq_genes)

# 查看缺失的基因
print(missing_genes)


# Perform ssGSEA
gsvapar <- ssgseaParam(expression_matrix, gene_signatures)
gsva_es <- gsva(gsvapar)

# Convert the results to a data frame for easier handling
ssgsea_results_df <- as.data.frame(gsva_es)
ssgsea_results_df_t <- t(ssgsea_results_df)
ssgsea_results_df_t <- as.data.frame(ssgsea_results_df_t)

GDNF <- ssgsea_results_df_t$GDNF_Up - ssgsea_results_df_t$GDNF_Down
gsva <- rbind(ssgsea_results_df, GDNF)
rownames(gsva)[5] <- "GDNF"
View(gsva)
gsva <- as.data.frame(gsva)
row_index_1 <- which(rownames(rnaseq) == "RET")
row_index_2 <- which(rownames(rnaseq) == "ERBB2")
gsva <- rbind(gsva, rnaseq[row_index_1, ])
rownames(gsva)[6] <- "RET"
gsva <- rbind(gsva, rnaseq[row_index_2, ])
rownames(gsva)[7] <- "ERBB2"

gsva_t <- as.data.frame(t(gsva))

# 检查 HER2 Status 是否正确
table(clinical$ER)
table(clinical$HER2)
table(clinical$histology)

clinical$sampleID <- rownames(clinical)

#分别提取erpos和erneg的样本
valid_erpos <- clinical %>%
  filter(ER == "Positive") %>% 
  pull(sampleID)  # 提取样本 ID
gsva_erpos <- gsva %>%
  select(all_of(intersect(valid_erpos, colnames(gsva))))

valid_erneg <- clinical %>%
  filter(ER == "Negative") %>% 
  pull(sampleID)  # 提取样本 ID
gsva_erneg <- gsva %>%
  select(all_of(intersect(valid_erneg, colnames(gsva))))

# 将 HER2 Status 转换为因子
clinical$HER2 <- as.factor(clinical$HER2)

library(dplyr)

# 删除 HER2 列为空的样本
#clinical_her2 <- clinical %>% filter(HER2_STATUS != "" & !is.na(HER2_STATUS))
#clinical_her2neg <- clinical %>% filter(Her2 == "Negative")
clinical_her2 <- subset(clinical, HER2 %in% c("Positive", "Negative"))

################################## RET signature ################################
######### whole dataset ##########
setdiff(colnames(gsva), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva["GDNF", ]),  # RET signature
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC:RET signature by primary/metastasis",
       x = "Tumor",
       y = "RET signature") +
  theme_minimal()

######### ERpos ##########
setdiff(colnames(gsva_erpos), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva_erpos))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva_erpos), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva_erpos[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva_erpos["GDNF", ]),  # RET signature
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC_ERpos:RET signature by primary/metastasis",
       x = "Tumor",
       y = "RET signature") +
  theme_minimal()


######### ERneg ##########
setdiff(colnames(gsva_erneg), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva_erneg))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva_erneg), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva_erneg[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva_erneg["GDNF", ]),  # RET signature
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC_erneg:RET signature by primary/metastasis",
       x = "Tumor",
       y = "RET signature") +
  theme_minimal()


################################## RET expression ################################
######### whole dataset ##########
setdiff(colnames(gsva), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva["RET", ]),  # RET expression
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC:RET expression by primary/metastasis",
       x = "Tumor",
       y = "RET expression") +
  theme_minimal()

######### ERpos ##########
setdiff(colnames(gsva_erpos), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva_erpos))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva_erpos), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva_erpos[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva_erpos["RET", ]),  # RET signature
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC_ERpos:RET expression by primary/metastasis",
       x = "Tumor",
       y = "RET expression") +
  theme_minimal()


######### ERneg ##########
setdiff(colnames(gsva_erneg), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva_erneg))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva_erneg), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva_erneg[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva_erneg["RET", ]),  # RET expression
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC_erneg:RET expression by primary/metastasis",
       x = "Tumor",
       y = "RET expression") +
  theme_minimal()

################################## HER2 expression ################################
######### whole dataset ##########
setdiff(colnames(gsva), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva["ERBB2", ]),  
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC:HER2 expression by primary/metastasis",
       x = "Tumor",
       y = "HER2 expression") +
  theme_minimal()

######### ERpos ##########
setdiff(colnames(gsva_erpos), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva_erpos))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva_erpos), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva_erpos[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva_erpos["ERBB2", ]),  # RET signature
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC_ERpos:HER2 expression by primary/metastasis",
       x = "Tumor",
       y = "HER2 expression") +
  theme_minimal()


######### ERneg ##########
setdiff(colnames(gsva_erneg), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva_erneg))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva_erneg), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva_erneg[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva_erneg["ERBB2", ]),  # RET expression
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC_erneg:HER2 expression by primary/metastasis",
       x = "Tumor",
       y = "HER2 expression") +
  theme_minimal()


############################## HER2 signature_Novel ################################
######### whole dataset ##########
setdiff(colnames(gsva), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva["Novel", ]),  
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC:HER2 signature_Novel by primary/metastasis",
       x = "Tumor",
       y = "HER2 signature_Novel") +
  theme_minimal()

######### ERpos ##########
setdiff(colnames(gsva_erpos), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva_erpos))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva_erpos), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva_erpos[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva_erpos["Novel", ]),  
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC_ERpos:HER2 signature_Novel by primary/metastasis",
       x = "Tumor",
       y = "HER2 signature_Novel") +
  theme_minimal()


######### ERneg ##########
setdiff(colnames(gsva_erneg), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva_erneg))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva_erneg), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva_erneg[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva_erneg["Novel", ]),  
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC_erneg:HER2 signature_Novel by primary/metastasis",
       x = "Tumor",
       y = "HER2 signature_Novel") +
  theme_minimal()


############################## HER2 signature_Established ################################
######### whole dataset ##########
setdiff(colnames(gsva), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva["Established", ]),  
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC:HER2 signature_Established by primary/metastasis",
       x = "Tumor",
       y = "HER2 signature_Established") +
  theme_minimal()

######### ERpos ##########
setdiff(colnames(gsva_erpos), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva_erpos))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva_erpos), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva_erpos[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva_erpos["Established", ]),  
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC_ERpos:HER2 signature_Established by primary/metastasis",
       x = "Tumor",
       y = "HER2 signature_Established") +
  theme_minimal()


######### ERneg ##########
setdiff(colnames(gsva_erneg), rownames(clinical))  # RNA-seq 中有，但临床数据中没有
setdiff(rownames(clinical), colnames(gsva_erneg))  # 临床数据中有，但 RNA-seq 没有
# 找到 RNA-seq 和临床数据的交集样本
common_samples <- intersect(colnames(gsva_erneg), rownames(clinical))

# 重新筛选数据集
gsva_common <- gsva_erneg[, common_samples]  # 筛选 RNA-seq 数据
clinical_common <- clinical[common_samples, ]  # 筛选临床数据

# whole cohort:RET signature by primary/metastasis
plot_data <- data.frame(
  Expression = as.numeric(gsva_erneg["Established", ]),  
  Status = clinical_common$tumourID
)

ggplot(plot_data, aes(x = Status, y = Expression, fill = Status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("M", "P")), method = "t.test") +  # 只比较这两组
  labs(title = "WCRC_erneg:HER2 signature_Established",
       x = "Tumor",
       y = "HER2 signature_Established") +
  theme_minimal()
