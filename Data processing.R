library(tidyverse)
library(edgeR)

rnaseq_raw <- read.csv("~/Desktop/scanb/GSE202203_RawCounts_gene_3207.csv")
View(rnaseq_raw)
clinical <- read.delim("~/Desktop/scanb/Clinical features SCANB.tsv")
View(clinical)

library(stringr)
clinical$SAMPLE <- str_sub(clinical$SAMPLE_ID, 4, 7)

#data processing
merged_rnaseq <- merged_rnaseq[, c("external_gene_name", setdiff(names(merged_rnaseq), "external_gene_name"))]

anyDuplicated(merged_rnaseq$external_gene_name)

duplicated_rows <- merged_rnaseq[duplicated(merged_rnaseq[, 1]), 1]  # Looking for replicated values
print(duplicated_rows)  # view duplicated values

merged_rnaseq[, 1] <- make.unique(as.character(merged_rnaseq[, 1]))  # 让重复的值变唯一
rownames(merged_rnaseq) <- merged_rnaseq[, 1]  # 设定行名

merged_rnaseq <- merged_rnaseq[!is.na(merged_rnaseq$external_gene_name), ]

merged_rnaseq <- merged_rnaseq[, -1]  # 删除第一列

merged_rnaseq <- merged_rnaseq[, 1:(ncol(merged_rnaseq) - 5)]


#row name
library(tibble)

colnames(rnaseq) <- gsub("\\.", "-", substr(colnames(rnaseq), 1, 7))
clinical$SAMPLE <- gsub("\\.", " ", substr(clinical$SAMPLE_ID, 1, 7))

# 去除 rnaseq 数据中的 Ensembl 版本号
rnaseq$sample <- sub("\\..*$", "", rnaseq$sample)
colnames(rnaseq)[which(colnames(rnaseq) == "sample")] <- "ensembl_gene_id"

# 合并 rnaseq 数据与 gene_map
merged_rnaseq <- merge(rnaseq, geneMap, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)

# 用 external gene name 替换原来的 gene 名（可选）
merged_df$gene <- ifelse(merged_df$external_gene_name == "" | is.na(merged_df$external_gene_name),
                         merged_df$ensembl_clean,
                         merged_df$external_gene_name)

# 删除中间变量列（可选）
merged_df$ensembl_clean <- NULL
merged_df$external_gene_name <- NULL


