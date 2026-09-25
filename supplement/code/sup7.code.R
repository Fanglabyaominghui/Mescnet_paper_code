library(Seurat)
library(Scissor)


print(paste("⏳ 正在处理:", cancer_type, "，请耐心等待..."))

# 1. 查询数据
query <- GDCquery(
  project = cancer_type,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# 2. 下载数据 (如果你的文件夹里已经下载好了，这步会自动识别并跳过)
GDCdownload(query, files.per.chunk = 50, directory = base_dir)

# 3. 整合数据为 SummarizedExperiment 对象 (这步会读取你截图里的那些散落的文件夹)
se_data <- GDCprepare(query, directory = base_dir)



library(TCGAbiolinks)
library(SummarizedExperiment)

# 指定你下载数据存放的根目录 (包含 TCGA-LUAD 和 TCGA-LUSC 文件夹的上一级目录)
base_dir <- "E:/deskup/spatial_luad/"
setwd(base_dir)

print("1. 正在从本地解析 LUAD 数据...")
query_luad <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
# 直接 prepare，它会自动读取本地存在的文件夹
data_luad <- GDCprepare(query_luad, directory = base_dir)

print("2. 正在从本地解析 LUSC 数据...")
query_lusc <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
data_lusc <- GDCprepare(query_lusc, directory = base_dir)

# ==========================================
# 3. 提取 TPM 并严格过滤肿瘤样本 (< 10)
# ==========================================
get_tumor_tpm <- function(se_obj) {
  mat <- assay(se_obj, "tpm_unstrand")
  # 提取第14-15位数字
  sample_types <- as.numeric(substr(colnames(mat), 14, 15))
  # < 10 为肿瘤 (01A 等)
  return(mat[, sample_types < 10])
}

tpm_luad <- get_tumor_tpm(data_luad)
tpm_lusc <- get_tumor_tpm(data_lusc)

print(paste("过滤后，LUAD 肿瘤样本数:", ncol(tpm_luad)))
print(paste("过滤后，LUSC 肿瘤样本数:", ncol(tpm_lusc)))

# ==========================================
# 4. 合并与 $log_2(TPM + 1)$ 转换
# ==========================================
common_genes <- intersect(rownames(tpm_luad), rownames(tpm_lusc))
bulk_expr <- cbind(tpm_luad[common_genes, ], tpm_lusc[common_genes, ])

# 核心转换
bulk_expr <- log2(bulk_expr + 1)


# 构建表型标签 (LUAD = 0, LUSC = 1)
phenotype <- c(rep(0, ncol(tpm_luad)), rep(1, ncol(tpm_lusc)))
tag <- c("LUAD", "LUSC")


setwd('E:/deskup/spatial_luad/data/')







# 载入所需的包
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(tibble)

# 1. 提取行名并去除 Ensembl ID 的版本号后缀
ensembl_ids_with_version <- rownames(bulk_expr)
clean_ensembl <- gsub("\\..*", "", ensembl_ids_with_version)

# 2. 获取 Ensembl 到 Symbol 的映射字典
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = clean_ensembl,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# 3. 将原矩阵转换为数据框，并加入 Symbol 和行平均表达量信息
expr_df <- as.data.frame(bulk_expr)
expr_df$Symbol <- gene_symbols
expr_df$mean_expr <- rowMeans(bulk_expr)

# 4. 清洗数据：去除没有匹配到 Symbol 的行 (NA)
expr_df_clean <- expr_df[!is.na(expr_df$Symbol), ]

# 5. 去重：如果有多个 Ensembl ID 对应同一个 Symbol，保留该 Symbol 平均表达量最大的一行
expr_df_unique <- expr_df_clean %>%
  group_by(Symbol) %>%
  slice_max(order_by = mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  as.data.frame()

# 6. 将整理好的数据框还原为标准的纯数字表达矩阵
# 提取原本的样本列 (去除最后两列：Symbol 和 mean_expr)
sample_cols <- 1:ncol(bulk_expr)
bulk_expr_symbol <- as.matrix(expr_df_unique[, sample_cols])

# 将唯一的 Gene Symbol 赋给行名
rownames(bulk_expr_symbol) <- expr_df_unique$Symbol




str(bulk_expr_symbol)
str(phenotype)
str(tag)

max(bulk_expr_symbol)
# 保存为 RData，准备上传服务器

setwd('E:/deskup/spatial_luad/data/')
save(bulk_expr_symbol, phenotype, tag, file = "TCGA_LUAD_LUSC_Scissor_Ref.RData")
print("🎉 本地数据处理完成，请将 TCGA_LUAD_LUSC_Scissor_Ref.RData 上传至服务器！")




Scissor<-function (bulk_dataset, sc_dataset, phenotype, tag = NULL, alpha = NULL, 
                   cutoff = 0.2, family = c("gaussian", "binomial", "cox"), 
                   Save_file = "Scissor_inputs.RData", Load_file = NULL) 
{
  library(Seurat)
  library(Matrix)
  library(limma) # 替换 preprocessCore
  if (is.null(Load_file)) {
    common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
    if (length(common) == 0) {
      stop("There is no common genes between the given single-cell and bulk samples.")
    }
    if (class(sc_dataset) == "Seurat") {
      # 替换为兼容 Seurat V5 的图层提取函数
      sc_exprs <- as.matrix(GetAssayData(sc_dataset, assay = "RNA", layer = "data"))
      network <- as.matrix(sc_dataset@graphs$RNA_snn)
    }
    else {
      sc_exprs <- as.matrix(sc_dataset)
      Seurat_tmp <- CreateSeuratObject(sc_dataset)
      Seurat_tmp <- FindVariableFeatures(Seurat_tmp, selection.method = "vst", 
                                         verbose = F)
      Seurat_tmp <- ScaleData(Seurat_tmp, verbose = F)
      Seurat_tmp <- RunPCA(Seurat_tmp, features = VariableFeatures(Seurat_tmp), 
                           verbose = F)
      Seurat_tmp <- FindNeighbors(Seurat_tmp, dims = 1:10, 
                                  verbose = F)
      network <- as.matrix(Seurat_tmp@graphs$RNA_snn)
    }
    diag(network) <- 0
    network[which(network != 0)] <- 1
    dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, 
    ])
    dataset1 <- limma::normalizeQuantiles(as.matrix(dataset0))
    rownames(dataset1) <- rownames(dataset0)
    colnames(dataset1) <- colnames(dataset0)
    Expression_bulk <- dataset1[, 1:ncol(bulk_dataset)]
    Expression_cell <- dataset1[, (ncol(bulk_dataset) + 
                                     1):ncol(dataset1)]
    X <- cor(Expression_bulk, Expression_cell)
    quality_check <- quantile(X)
    print("|**************************************************|")
    print("Performing quality-check for the correlations")
    print("The five-number summary of correlations:")
    print(quality_check)
    print("|**************************************************|")
    if (quality_check[3] < 0.01) {
      warning("The median correlation between the single-cell and bulk samples is relatively low.")
    }
    if (family == "binomial") {
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)) {
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
      }
      else {
        print(sprintf("Current phenotype contains %d %s and %d %s samples.", 
                      z[1], tag[1], z[2], tag[2]))
        print("Perform logistic regression on the given phenotypes:")
      }
    }
    if (family == "gaussian") {
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)) {
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
      }
      else {
        tmp <- paste(z, tag)
        print(paste0("Current phenotype contains ", 
                     paste(tmp[1:(length(z) - 1)], collapse = ", "), 
                     ", and ", tmp[length(z)], " samples."))
        print("Perform linear regression on the given phenotypes:")
      }
    }
    if (family == "cox") {
      Y <- as.matrix(phenotype)
      if (ncol(Y) != 2) {
        stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")
      }
      else {
        print("Perform cox regression on the given clinical outcomes:")
      }
    }
    save(X, Y, network, Expression_bulk, Expression_cell, 
         file = Save_file)
  }
  else {
    load(Load_file)
  }
  if (is.null(alpha)) {
    alpha <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
               0.6, 0.7, 0.8, 0.9)
  }
  for (i in 1:length(alpha)) {
    set.seed(123)
    fit0 <- APML1(X, Y, family = family, penalty = "Net", 
                  alpha = alpha[i], Omega = network, nlambda = 100, 
                  nfolds = min(10, nrow(X)))
    fit1 <- APML1(X, Y, family = family, penalty = "Net", 
                  alpha = alpha[i], Omega = network, lambda = fit0$lambda.min)
    if (family == "binomial") {
      Coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])
    }
    else {
      Coefs <- as.numeric(fit1$Beta)
    }
    Cell1 <- colnames(X)[which(Coefs > 0)]
    Cell2 <- colnames(X)[which(Coefs < 0)]
    percentage <- (length(Cell1) + length(Cell2))/ncol(X)
    print(sprintf("alpha = %s", alpha[i]))
    print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", 
                  length(Cell1), length(Cell2)))
    print(sprintf("The percentage of selected cell is: %s%%", 
                  formatC(percentage * 100, format = "f", digits = 3)))
    if (percentage < cutoff) {
      break
    }
    cat("\n")
  }
  print("|**************************************************|")
  return(list(para = list(alpha = alpha[i], lambda = fit0$lambda.min, 
                          family = family), Coefs = Coefs, Scissor_pos = Cell1, 
              Scissor_neg = Cell2))
}






library(Seurat)
library(Scissor)

# 1. 设置服务器根路径并加载刚才上传的 TCGA 参考数据
server_base <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/"
setwd(server_base)
load("TCGA_LUAD_LUSC_Scissor_Ref.RData")

# 2. 定义自动化处理空间数据并跑 Scissor 的函数
run_spatial_scissor <- function(sample_dir, out_csv_path, rdata_save_name, alpha_vec, cutoff_val) {
  print(sprintf("⏳ 正在加载空转数据: %s", sample_dir))
  
  # A. 加载高分辨率空转多边形数据
  obj <- Load10X_Spatial(
    data.dir = sample_dir, 
    bin.size = "polygons", 
    image.name = "tissue_hires_image.png"
  )
  
  # B. 【终极防错修复】：利用 do.call 动态抓取原生名字，自适应重命名为 RNA
  # 不管叫 Spatial.Polygons、Spatial.polygons 还是 Spatial，都能完美兼容
  native_assay <- DefaultAssay(obj)
  print(sprintf("  -> 成功抓取到当前原生 Assay 名字为: %s", native_assay))
  
  if (native_assay != "RNA") {
    rename_args <- list(obj)
    rename_args[[native_assay]] <- "RNA"
    obj <- do.call(RenameAssays, rename_args)
    DefaultAssay(obj) <- "RNA"
    print("  -> 已经自动完成 Assay 重命名转换：原生图层 -> RNA")
  }
  
  # C. 计算单细胞网络 (RNA_snn)，这是 Scissor alpha 惩罚项必需的输入
  print("⏳ 正在计算 PCA 和空间邻接网络 (RNA_snn)...")
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)
  
  # D. 运行 Scissor
  print(sprintf("⏳ 正在运行 Scissor (Cutoff=%.2f)...", cutoff_val))
  infos <- Scissor(
    bulk_dataset = bulk_expr_symbol, 
    sc_dataset = obj, 
    phenotype = phenotype, 
    tag = tag, 
    alpha = alpha_vec, 
    cutoff = cutoff_val, 
    family = "binomial",
    Save_file = rdata_save_name
  )
  
  # E. 提取细胞 ID，初始化标签
  cell_ids <- colnames(obj)
  scissor_labels <- rep("Background", length(cell_ids))
  names(scissor_labels) <- cell_ids
  
  # F. 严格映射标签：TCGA 表型中 LUAD=0, LUSC=1
  scissor_labels[infos$Scissor_pos] <- "LUSC_like"
  scissor_labels[infos$Scissor_neg] <- "LUAD_like"
  
  # G. 组装并导出 CSV
  df_export <- data.frame(
    cell_id = names(scissor_labels),
    Scissor_Subtype = scissor_labels,
    stringsAsFactors = FALSE
  )
  
  print(sprintf("💾 正在保存分类结果至: %s", out_csv_path))
  write.csv(df_export, file = out_csv_path, row.names = FALSE)
  
  print("📊 当前样本分类统计:")
  print(table(df_export$Scissor_Subtype))
  print("---------------------------------------------------")
}


# ==========================================
# 执行 AXB_1913 (背景较干净，参数适中)
# ==========================================
dir_1913 <- paste0(server_base, "AXB-1913-0017-A1/AXB_1913_Cellbin_Result/outs/")
csv_1913 <- paste0(server_base, "AXB-1913-0017-A1/AXB_1913_Scissor_labels.csv")
rdata_1913 <- paste0(server_base, "AXB-1913-0017-A1/Scissor_1913_result.RData")

run_spatial_scissor(
  sample_dir = dir_1913,
  out_csv_path = csv_1913,
  rdata_save_name = rdata_1913,
  alpha_vec = c(0.005, 0.01, 0.05),
  cutoff_val = 0.40
)

# ==========================================
# 执行 AXB_6123 (弥漫性强，参数激进)
# ==========================================
dir_6123 <- paste0(server_base, "AXB-6123-2431-6976-A1/AXB_6123_Cellbin_Result/outs/")
csv_6123 <- paste0(server_base, "AXB-6123-2431-6976-A1/AXB_6123_Scissor_labels.csv")
rdata_6123 <- paste0(server_base, "AXB-6123-2431-6976-A1/Scissor_6123_result.RData")

run_spatial_scissor(
  sample_dir = dir_6123,
  out_csv_path = csv_6123,
  rdata_save_name = rdata_6123,
  alpha_vec = c(0.001, 0.005, 0.01),
  cutoff_val = 0.60
)


# 载入必要的绘图包
library(ggplot2)

print("1. 正在提取高变基因并计算 PCA...")

# 计算每个基因的方差，取 Top 3000 加速 PCA 计算且降噪
gene_vars <- apply(bulk_expr_symbol, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:3000])
pca_input <- bulk_expr_symbol[top_genes, ]

# 运行主成分分析 (转置矩阵：PCA要求行为样本，列为基因)
pca_res <- prcomp(t(pca_input), scale. = TRUE)

# 提取前两个主成分
pca_df <- as.data.frame(pca_res$x[, 1:2])

# 加入表型标签 (phenotype: 0 = LUAD, 1 = LUSC)
pca_df$Phenotype <- ifelse(phenotype == 0, "LUAD", "LUSC")

# 计算 PC1 和 PC2 的方差解释百分比
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
pc1_label <- sprintf("PC1 (%.1f%%)", var_explained[1] * 100)
pc2_label <- sprintf("PC2 (%.1f%%)", var_explained[2] * 100)

print("2. 正在生成 PCA 矢量图 (严格锁定 Arial 7pt, 纯黑)...")

# 生成 ggplot 图像
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size = 1.2, alpha = 0.8, stroke = 0) +
  # Scissor 经典红蓝配色
  scale_color_manual(values = c("LUAD" = "#3C5488", "LUSC" = "#E64B35")) +
  labs(x = pc1_label, y = pc2_label) +
  theme_bw() +
  theme(
    # 严格锁定 7pt，纯黑，常规粗细
    text = element_text(family = "sans", size = 7, color = "black", face = "plain"),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 7, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  )

# 导出到服务器目录
out_pdf <- paste0(server_base, "TCGA_LUAD_vs_LUSC_PCA_7pt.pdf")

ggsave(
  filename = out_pdf,
  plot = p,
  width = 4,
  height = 2.5,
  units = "in",
  device = "pdf",
  useDingbats = FALSE
)

print(paste0("✅ PCA 绘图完成！纯矢量 PDF 已保存至: ", out_pdf))








# ==========================================
# 1. 数据准备与名称大清洗
# ==========================================
cat("⏳ 1. 正在读取外部打分并提取共有细胞...\n")

base_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/"

rctd_df <- read.csv(paste0(base_dir, "RCTD_weights_percentage.csv"), row.names = 1)
disease_df <- read.csv(paste0(base_dir, "AXB_6123_UCell_ADC_SCC_scores.csv"), row.names = 1)

valid_cells <- intersect(rownames(seurat_obj@meta.data), rownames(rctd_df))
valid_cells <- intersect(valid_cells, rownames(disease_df))

# 提取特征
mod_cols <- grep("Module_[0-9]+_UCell", colnames(seurat_obj@meta.data), value = TRUE)
raw_modules <- seurat_obj@meta.data[valid_cells, mod_cols]
raw_cells <- rctd_df[valid_cells, ]
raw_diseases <- disease_df[valid_cells, c("ADC_score", "SCC_score")]

# 🚨 名称大清洗！
# 1. 将 Module_1_UCell 变为 M1
colnames(raw_modules) <- gsub("Module_([0-9]+)_UCell", "M\\1", colnames(raw_modules))

# 2. 将 ADC/SCC 改为 LUAD/LUSC
colnames(raw_diseases) <- gsub("ADC_score", "LUAD", colnames(raw_diseases))
colnames(raw_diseases) <- gsub("SCC_score", "LUSC", colnames(raw_diseases))

# 3. 将细胞名字里的下划线替换为空格
colnames(raw_cells) <- gsub("_", " ", colnames(raw_cells))

cat("⚖️ 2. 正在执行各自特征内部的 Scale 标准化并横向拼合...\n")

mat_modules <- scale(as.matrix(raw_modules))
mat_cells <- scale(as.matrix(raw_cells))
mat_diseases <- scale(as.matrix(raw_diseases))

# 把细胞和疾病在列方向拼在一起作为 X 轴
mat_x_axis <- cbind(mat_cells, mat_diseases)

# ==========================================
# 2. 计算相关性 r 与 P 值矩阵
# ==========================================
cat("⏳ 3. 正在计算大合并矩阵的皮尔逊相关性...\n")

cor_mat <- cor(mat_modules, mat_x_axis, method = "pearson", use = "pairwise.complete.obs")

p_mat <- matrix(1, nrow = ncol(mat_modules), ncol = ncol(mat_x_axis))

for(i in 1:ncol(mat_modules)) {
  for(j in 1:ncol(mat_x_axis)) {
    try({
      test <- cor.test(mat_modules[, i], mat_x_axis[, j], method = "pearson")
      p_mat[i, j] <- test$p.value
    }, silent = TRUE)
  }
}

rownames(cor_mat) <- colnames(mat_modules)
colnames(cor_mat) <- colnames(mat_x_axis)
rownames(p_mat) <- colnames(mat_modules)
colnames(p_mat) <- colnames(mat_x_axis)

cor_mat[is.na(cor_mat)] <- 0
p_mat[is.na(p_mat)] <- 1

# ==========================================
# 3. 生成 7pt 纯矢量大热图 (宽6，高1.6)
# ==========================================
cat("🎨 4. 正在生成名称清洗后的纯矢量 7pt 相关性热图...\n")

my_palette <- colorRampPalette(c("#4393C3", "#FFFFFF", "#D6604D"))(200)

out_path <- paste0(base_dir, "Figure_Combined_Cor_Heatmap.pdf")

# 画布锁定：宽 6，高 1.6
cairo_pdf(out_path, width = 8, height = 2, family = "Arial")

# 排版锁定：Arial 7号字体，黑色，不加粗
par(ps = 7, family = "Arial", font = 1, col.axis = "black")

corrplot::corrplot(
  cor_mat,
  method = "circle",
  order = "original",
  tl.pos = "lt",
  tl.col = "black",
  tl.cex = 1,          # 继承 ps=7
  tl.srt = 45,
  cl.pos = "r",
  cl.cex = 1,          # 继承 ps=7
  cl.offset = 0.5,
  col = my_palette,
  p.mat = p_mat,
  sig.level = 0.05,
  insig = "blank",
  mar = c(1, 1, 1, 1)
)

dev.off()

# 保存清洗后的矩阵数据
save(cor_mat, p_mat, file = paste0(base_dir, "Figure_Combined_Cor_Matrix_Clean.RData"))

cat("🎉 完美收工！最终清洗版热图已保存至:", out_path, "\n")

# ==========================================
# 0. 载入必备包与排版设定 (严格 Arial 7pt, 黑色)
# ==========================================
if (!requireNamespace("ggrastr", quietly = TRUE)) install.packages("ggrastr")

library(uwot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggrastr)

# 强制 Arial 7pt 纯黑不加粗
family_use <- "Arial"
font_style_7pt <- element_text(size = 7, color = "black", face = "plain", family = family_use)

# ==========================================
# 1. 提取底层数据并清洗名称 (M1 - M7)
# ==========================================
cat("⏳ 1. 正在从 seurat_obj 提取 TOM 矩阵与 kME...\n")

datExpr <- seurat_obj@misc$MEscnet$datExpr
moduleColors <- seurat_obj@misc$MEscnet$moduleColors
tom <- seurat_obj@misc$MEscnet$bayes_cor1

if(is.null(rownames(tom))) rownames(tom) <- colnames(datExpr)
if(is.null(colnames(tom))) colnames(tom) <- colnames(datExpr)

# 提取特征值与 kME
MEs0 <- WGCNA::moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- WGCNA::orderMEs(MEs0)
datKME <- WGCNA::signedKME(datExpr, MEs, outputColumnName="kME_M")

# 清洗模块名为 M1, M2...
mod_nums <- gsub("[^0-9]", "", as.character(moduleColors))
gene_modules <- paste0("M", mod_nums)
names(gene_modules) <- colnames(datExpr)

old_cols <- colnames(datKME)
kme_nums <- gsub("[^0-9]", "", old_cols)
colnames(datKME) <- paste0("kME_M", kme_nums)

gene_info <- data.frame(
  gene = names(gene_modules),
  module = unname(gene_modules),
  stringsAsFactors = FALSE
) %>%
  filter(!module %in% c("M", "M0", "Grey", "grey"))

# 匹配 kME
gene_info$kME <- sapply(1:nrow(gene_info), function(i) {
  g <- gene_info$gene[i]
  m <- gene_info$module[i]
  col_name <- paste0("kME_", m)
  if (g %in% rownames(datKME) && col_name %in% colnames(datKME)) {
    return(datKME[g, col_name])
  } else {
    return(0)
  }
})

# ==========================================
# 2. 执行无监督 UMAP 降维
# ==========================================
cat("🗺️ 2. 正在执行无监督 UMAP 降维...\n")

n_hubs <- 10
n_genes_per_mod <- 200

# 提取 Hub 基因
hub_genes <- gene_info %>%
  group_by(module) %>%
  slice_max(order_by = kME, n = n_hubs, with_ties = FALSE) %>%
  pull(gene)

selected_genes <- gene_info %>%
  group_by(module) %>%
  slice_max(order_by = kME, n = n_genes_per_mod, with_ties = FALSE)

valid_selected <- intersect(selected_genes$gene, rownames(tom))
valid_hubs <- intersect(hub_genes, colnames(tom))

feature_mat <- tom[valid_selected, valid_hubs]

set.seed(42)
hub_umap <- uwot::umap(
  X = feature_mat,
  n_neighbors = 15,
  min_dist = 0.1,
  spread = 1,
  metric = "cosine"
)

umap_df <- as.data.frame(hub_umap)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$gene <- valid_selected
umap_df <- left_join(umap_df, selected_genes, by = "gene")

# 对每个模块内部的 kME 进行缩放
umap_df <- umap_df %>%
  group_by(module) %>%
  mutate(kME_scaled = (kME - min(kME)) / (max(kME) - min(kME) + 1e-6)) %>%
  ungroup()

# ==========================================
# 3. 绘制 7pt 顶刊格式图表 (自定义配色)
# ==========================================
cat("🎨 3. 正在生成定制配色的栅格化 UMAP 图...\n")

# 自定义颜色字典
mod_palette <- c(
  "M1" = "#D94848",    # 红色
  "M2" = "#7F76D7",    # Slate purple
  "M3" = "#3C5488",    # 蓝色
  "M4" = "#F59929",    # Orange
  "M5" = "#41C4C5",    # Teal/Cyan
  "M6" = "#8491B4",    # 柔和灰蓝色
  "M7" = "#91D1C2"     # 薄荷绿
)

hubs_to_label <- umap_df %>%
  group_by(module) %>%
  slice_max(order_by = kME, n = 3, with_ties = FALSE)

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  # 开启 600 DPI 散点栅格化
  geom_point_rast(
    aes(color = module, size = kME_scaled * 2),
    alpha = 0.9,
    raster.dpi = 600,
    stroke = 0
  ) +
  geom_text_repel(
    data = hubs_to_label,
    aes(label = gene),
    size = 7 / .pt,        # 严格 7pt
    fontface = "italic",   # 基因名标准斜体
    family = family_use,
    color = "black",
    bg.color = "white",
    bg.r = 0.15,
    max.overlaps = Inf,
    segment.color = "grey50",
    segment.size = 0.3
  ) +
  scale_color_manual(values = mod_palette) +
  scale_size_identity() +
  theme_void() +
  theme(
    text = font_style_7pt,
    legend.text = font_style_7pt,
    legend.title = element_blank(),
    legend.position = "right",
    aspect.ratio = 1
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# ==========================================
# 4. 导出 PDF
# ==========================================
out_pdf <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/Figure_Unsupervised_Module_UMAP_7pt.pdf"

cairo_pdf(out_pdf, width = 2.7, height = 1.8, family = family_use)
print(p_umap)
dev.off()

cat("✅ 搞定！定制配色的 UMAP 已经生成，散点栅格化，且文字均锁定为 Arial 7pt！\n")


#fig6
# ==========================================
# 0. 载入必需包
# ==========================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

cat("⏳ [Curve] 正在提取 TP63 表达量并计算一阶导数...\n")

local_csv_path <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/AXB_6123_Trajectory_Export_for_Local_2.csv"
pt_data <- read.csv(local_csv_path)

meta_data <- seurat_obj@meta.data
meta_data$cell_id <- rownames(meta_data)
tp63_expr <- FetchData(seurat_obj, vars = c("TP63"))
tp63_expr$cell_id <- rownames(tp63_expr)

meta_combined <- merge(meta_data, tp63_expr, by = "cell_id")
pt_data <- merge(pt_data, meta_combined, by = "cell_id", all.x = TRUE)
pt_data <- pt_data[order(pt_data$ptime), ]

# ==========================================
# 1. 计算 AST 分数与一阶导数
# ==========================================
pt_data$M1_Scaled <- scale(pt_data$Module_1_UCell)[, 1]
pt_data$M3_Scaled <- scale(pt_data$Module_3_UCell)[, 1]
pt_data$AST_Score <- pt_data$M1_Scaled - pt_data$M3_Scaled

fit_ast <- loess(AST_Score ~ ptime, data = pt_data, span = 0.75)
seq_pt <- seq(min(pt_data$ptime, na.rm = TRUE), max(pt_data$ptime, na.rm = TRUE), length.out = 1000)
pred_ast <- predict(fit_ast, newdata = data.frame(ptime = seq_pt))

deriv_ast <- diff(pred_ast) / diff(seq_pt)
deriv_ast <- c(deriv_ast[1], deriv_ast)

deriv_df <- data.frame(ptime = seq_pt, AST_Derivative = deriv_ast)

# 🌟 核心修改：转折点定义为一阶导数最大值
turn_ptime <- seq_pt[which.max(deriv_ast)]

# ==========================================
# 2. 严格排版规范 (绝对 Arial 7pt, 黑色)
# ==========================================
state_colors <- c("LUAD+" = "#3C5488", "LUSC+" = "#E64B35", "Transitional" = "#B0B0B0")

font_settings <- theme(
  text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
  axis.text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
  axis.title = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
  legend.text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
  legend.title = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
  legend.position = "right",
  panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
  panel.grid.major = element_line(color = "#E5E5E5", linetype = "dotted", linewidth = 0.3),
  panel.grid.minor = element_blank(),
  legend.key.size = unit(0.3, "cm")
)

# ==========================================
# 3. 绘图 - 上半部分：TP63 表达量散点与演化 (隐藏 X 轴)
# ==========================================
p_top <- ggplot(pt_data, aes(x = ptime, y = TP63)) +
  geom_vline(xintercept = turn_ptime, color = "grey50", linetype = "dotted", linewidth = 0.5) +
  geom_point(aes(color = Tumor_Label), alpha = 0.5, size = 0.5, stroke = 0) +
  scale_color_manual(values = state_colors, name = "Cell State") +
  geom_smooth(method = "loess", span = 0.75, color = "black", linewidth = 0.8, se = FALSE) +
  labs(x = NULL, y = "TP63 Expression") +
  theme_classic() + font_settings +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) # 隐藏上图 X 轴

# ==========================================
# 4. 绘图 - 下半部分：AST 一阶导数
# ==========================================
p_bottom <- ggplot(deriv_df, aes(x = ptime, y = AST_Derivative)) +
  geom_hline(yintercept = 0, color = "#E64B35", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = turn_ptime, color = "grey50", linetype = "dotted", linewidth = 0.5) +
  geom_line(color = "black", linewidth = 0.8) +
  labs(x = "Pseudotime", y = "AST 1st Derivative\n(Rate of Change)") +
  theme_classic() + font_settings

# ==========================================
# 5. 拼接与导出 PDF
# ==========================================
# 上面 TP63 (占 2 份高度)，下面 导数 (占 1 份高度)
combined_plot <- p_top / p_bottom + plot_layout(heights = c(2, 1), guides = 'collect')

out_pdf <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/TP63_and_AST_Deriv_7pt.pdf"
cairo_pdf(out_pdf, width = 3.5, height = 5, family = "Arial")
print(combined_plot)
dev.off()

cat("✅ [Curve] 完美！TP63在上，一阶导在下的拼接图 (Arial 7pt) 已保存为 TP63_and_AST_Deriv_7pt.pdf\n")




# ==========================================
# 0. 载入必需包
# ==========================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)

cat("⏳ [Boxplot] 正在提取 CNV 肿瘤细胞的 TP63 表达数据...\n")

base_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1"

# 1. 读取 InferCNV 结果文件
cnv_files <- Sys.glob(file.path(base_dir, "InferCNV_Results/*Cancer_Cells_Only_With_UCell.csv"))
if (length(cnv_files) == 0) {
  stop("❌ 未找到 InferCNV 结果文件，请检查路径！")
}
cnv_df <- read.csv(cnv_files[1])

# 2. 从 Seurat 对象中提取 metadata 与 TP63 表达量
meta_data <- seurat_obj@meta.data
meta_data$cell_id <- rownames(meta_data)

tp63_expr <- FetchData(seurat_obj, vars = c("TP63"))
tp63_expr$cell_id <- rownames(tp63_expr)

seurat_combined <- merge(meta_data, tp63_expr, by = "cell_id")

# 3. 筛选 CNV 肿瘤细胞
cnv_df$join_id <- gsub("[^0-9]", "", cnv_df[, 1])
seurat_combined$join_id <- gsub("[^0-9]", "", seurat_combined$cell_id)

cancer_join_ids <- cnv_df %>%
  filter(malignant_status == "Cancer") %>%
  pull(join_id)

plot_cells <- seurat_combined %>%
  filter(join_id %in% cancer_join_ids)

# 4. 智能匹配细胞状态标签 (Tumor_Label)
if (!"Tumor_Label" %in% colnames(plot_cells)) {
  scissor_path <- file.path(base_dir, "AXB_6123_Scissor_labels.csv")
  if (file.exists(scissor_path)) {
    scissor_df <- read.csv(scissor_path)
    scissor_df$join_id <- gsub("[^0-9]", "", scissor_df[, 1])
    colnames(scissor_df)[2] <- "Scissor_Subtype"
    
    plot_cells <- merge(plot_cells, scissor_df[, c("join_id", "Scissor_Subtype")], 
                        by = "join_id", all.x = TRUE)
    
    plot_cells$Tumor_Label <- "Transitional"
    plot_cells$Tumor_Label[plot_cells$Scissor_Subtype == "LUAD_like"] <- "LUAD+"
    plot_cells$Tumor_Label[plot_cells$Scissor_Subtype == "LUSC_like"] <- "LUSC+"
  } else {
    plot_cells$Tumor_Label <- "Transitional"
  }
}

plot_df <- plot_cells %>%
  select(cell_id, Tumor_Label, TP63) %>%
  na.omit()

plot_df$Tumor_Label <- factor(plot_df$Tumor_Label, 
                              levels = c("LUAD+", "Transitional", "LUSC+"))

cat(paste0("📌 成功锁定 ", nrow(plot_df), " 个肿瘤细胞，准备绘制纯箱线图...\n"))

# ==========================================
# 5. 绘制纯箱线图 (去散点，加 p 值)
# ==========================================
my_comparisons <- list(
  c("LUAD+", "Transitional"),
  c("Transitional", "LUSC+"),
  c("LUAD+", "LUSC+")
)

state_fill_colors <- c("LUAD+" = "#3C548880", 
                       "Transitional" = "#B0B0B0", 
                       "LUSC+" = "#E64B3580")

p_box <- ggplot(plot_df, aes(x = Tumor_Label, y = TP63, fill = Tumor_Label)) +
  # 纯箱线图主体
  geom_boxplot(width = 0.5, linewidth = 0.5,
               outlier.shape = 16, outlier.size = 0.8, outlier.alpha = 0.4) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.format", size = 7 / .pt) +
  scale_fill_manual(values = state_fill_colors) +
  labs(x = NULL, y = "TP63 Expression") +
  theme_classic() +
  theme(
    text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
    axis.text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
    axis.title = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    panel.grid.major = element_line(color = "#E5E5E5", linetype = "dotted", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# ==========================================
# 6. 导出 7pt 高精度矢量 PDF
# ==========================================
out_pdf <- file.path(base_dir, "CNV_Tumor_Cells_TP63_PureBoxplot_7pt.pdf")

cairo_pdf(out_pdf, width = 3.2, height = 4, family = "Arial")
print(p_box)
dev.off()

cat(paste0("✅ 完美！清爽版 TP63 箱线图已保存至:\n-> ", out_pdf, "\n"))









