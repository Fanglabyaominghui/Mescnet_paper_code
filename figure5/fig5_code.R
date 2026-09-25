# 载入所需的包
library(Seurat)
library(spacexr)
library(Matrix)

# =======================================================
# 1. 准备单细胞参考集 (Reference)
# =======================================================
print("1. 正在读取单细胞参考数据集...")
# 【修改路径】：指向你刚才导出的 13 大类数据文件夹
ref_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/sp/"

genes <- read.csv(paste0(ref_dir, "reference_genes.csv"), header=FALSE)$V1
cells <- read.csv(paste0(ref_dir, "reference_cells.csv"), header=FALSE)$V1

counts <- readMM(paste0(ref_dir, "reference_counts.mtx"))
counts <- as(counts, "CsparseMatrix")
counts <- round(counts)

rownames(counts) <- genes
colnames(counts) <- cells

meta <- read.csv(paste0(ref_dir, "reference_meta.csv"), row.names=1)
# 清理细胞类型名称中的斜杠和空格（针对 Macrophage/Monocyte 等）
meta$RCTD_annotation <- gsub("/", "_", meta$RCTD_annotation)
meta$RCTD_annotation <- gsub(" ", "_", meta$RCTD_annotation)

cell_types <- as.factor(meta$RCTD_annotation)
names(cell_types) <- rownames(meta)

reference <- Reference(counts, cell_types)

# =======================================================
# 2. 读取 HD 细胞切割数据并提取 Query
# =======================================================
print("2. 正在读取 Visium HD 细胞切割数据...")
# 【修改路径】：指向你的 segmented_outputs 目录
sp_path <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-1913-0017-A1/AXB_1913_Cellbin_Result/outs/segmented_outputs/"

# 【核心修改 1】：直接使用 Load10X_Spatial 读取多边形数据
luad_sp <- Load10X_Spatial(
  data.dir = sp_path,
  bin.size = "polygons",
  image.name = "tissue_hires_image.png"
)

# 提取表达矩阵
sp_counts <- luad_sp[["Spatial"]]$counts

# 提取原始多边形坐标
coords_raw <- GetTissueCoordinates(luad_sp)

# 【核心修改 2】：利用 aggregate 计算多边形的质心 (Centroid)
if ("cell" %in% colnames(coords_raw)) {
  # 如果列名是 cell，按 cell 分组求 x,y 的平均值
  coords <- aggregate(cbind(x, y) ~ cell, data = coords_raw, FUN = mean)
  rownames(coords) <- coords$cell
} else if ("barcode" %in% colnames(coords_raw)) {
  # 如果列名是 barcode，按 barcode 分组求平均
  coords <- aggregate(cbind(x, y) ~ barcode, data = coords_raw, FUN = mean)
  rownames(coords) <- coords$barcode
} else {
  coords <- coords_raw
}

# 提取最终的 x, y 列构建坐标矩阵
coords <- coords[, c("x", "y")]
nUMI <- colSums(sp_counts)

# 构建 RCTD 的 Query 对象
query <- SpatialRNA(coords, sp_counts, nUMI)















# =======================================================
# 2. 读取 HD 细胞切割数据并提取 Query
# =======================================================
print("2. 正在读取 Visium HD 细胞切割数据...")

# 【核心修复】：路径往后退一步，只停留在 outs/ 目录
sp_path <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-1913-0017-A1/AXB_1913_Cellbin_Result/outs/"

# 后面不需要动，直接运行
luad_sp <- Load10X_Spatial(
  data.dir = sp_path,
  bin.size = "polygons",
  image.name = "tissue_hires_image.png"
)

# 提取表达矩阵
sp_counts <- luad_sp[["Spatial.Polygons"]]$counts

# 提取原始多边形坐标
coords_raw <- GetTissueCoordinates(luad_sp)

# 利用 aggregate 计算多边形的质心 (Centroid)
if ("cell" %in% colnames(coords_raw)) {
  coords <- aggregate(cbind(x, y) ~ cell, data = coords_raw, FUN = mean)
  rownames(coords) <- coords$cell
} else if ("barcode" %in% colnames(coords_raw)) {
  coords <- aggregate(cbind(x, y) ~ barcode, data = coords_raw, FUN = mean)
  rownames(coords) <- coords$barcode
} else {
  coords <- coords_raw
}

# 提取最终的 x, y 列构建坐标矩阵
coords <- coords[, c("x", "y")]
nUMI <- colSums(sp_counts)

# 构建 RCTD 的 Query 对象
query <- SpatialRNA(coords, sp_counts, nUMI)





# =======================================================
# 3. 创建并运行 RCTD
# =======================================================
cat("3. 正在初始化并运行 RCTD (UMI_min = 0, Doublet 模式)...\n")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)

# UMI_min = 0 确保所有细胞均参与计算
myRCTD <- create.RCTD(query, reference, max_cores = 24, UMI_min = 0)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

print("RCTD 反卷积运行完成！")

# 保存 RCTD 对象 (存到上级目录，避免弄脏原始数据文件夹)
save_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-1913-0017-A1/"
saveRDS(myRCTD, file = paste0(save_dir, 'luad_Cellbin_RCTD.RDS'))

# =======================================================
# 4. 提取保存 RCTD 结果
# =======================================================
print("4. 正在提取并处理 RCTD 结果...")

# 提取 doublet 模式的详细分类结果
results <- myRCTD@results
results_df <- results$results_df

write.csv(results_df, file = paste0(save_dir, "RCTD_results_df.csv"), quote = FALSE)

print("🎉 细胞切割分类结果已保存！随时可以导入 Scanpy 或 Seurat 画图！")



















#data2


# 载入所需的包
library(Seurat)
library(spacexr)
library(Matrix)

# =======================================================
# 1. 准备单细胞参考集 (Reference)
# =======================================================
print("1. 正在读取单细胞参考数据集...")
# 【修改路径】：指向你刚才导出的 13 大类数据文件夹
ref_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/sp/"

genes <- read.csv(paste0(ref_dir, "reference_genes.csv"), header=FALSE)$V1
cells <- read.csv(paste0(ref_dir, "reference_cells.csv"), header=FALSE)$V1

counts <- readMM(paste0(ref_dir, "reference_counts.mtx"))
counts <- as(counts, "CsparseMatrix")
counts <- round(counts)

rownames(counts) <- genes
colnames(counts) <- cells

meta <- read.csv(paste0(ref_dir, "reference_meta.csv"), row.names=1)
# 清理细胞类型名称中的斜杠和空格（针对 Macrophage/Monocyte 等）
meta$RCTD_annotation <- gsub("/", "_", meta$RCTD_annotation)
meta$RCTD_annotation <- gsub(" ", "_", meta$RCTD_annotation)

cell_types <- as.factor(meta$RCTD_annotation)
names(cell_types) <- rownames(meta)

reference <- Reference(counts, cell_types)

# =======================================================
# 2. 读取 HD 细胞切割数据并提取 Query
# =======================================================
print("2. 正在读取 Visium HD 细胞切割数据...")
# 【修改路径】：指向你的 segmented_outputs 目录
sp_path <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-1913-0017-A1/AXB_1913_Cellbin_Result/outs/segmented_outputs/"

# 【核心修改 1】：直接使用 Load10X_Spatial 读取多边形数据
luad_sp <- Load10X_Spatial(
  data.dir = sp_path,
  bin.size = "polygons",
  image.name = "tissue_hires_image.png"
)

# 提取表达矩阵
sp_counts <- luad_sp[["Spatial"]]$counts

# 提取原始多边形坐标
coords_raw <- GetTissueCoordinates(luad_sp)

# 【核心修改 2】：利用 aggregate 计算多边形的质心 (Centroid)
if ("cell" %in% colnames(coords_raw)) {
  # 如果列名是 cell，按 cell 分组求 x,y 的平均值
  coords <- aggregate(cbind(x, y) ~ cell, data = coords_raw, FUN = mean)
  rownames(coords) <- coords$cell
} else if ("barcode" %in% colnames(coords_raw)) {
  # 如果列名是 barcode，按 barcode 分组求平均
  coords <- aggregate(cbind(x, y) ~ barcode, data = coords_raw, FUN = mean)
  rownames(coords) <- coords$barcode
} else {
  coords <- coords_raw
}

# 提取最终的 x, y 列构建坐标矩阵
coords <- coords[, c("x", "y")]
nUMI <- colSums(sp_counts)

# 构建 RCTD 的 Query 对象
query <- SpatialRNA(coords, sp_counts, nUMI)





# =======================================================
# 2. 读取 HD 细胞切割数据并提取 Query
# =======================================================
print("2. 正在读取 Visium HD 细胞切割数据...")

# 【核心修复】：路径往后退一步，只停留在 outs/ 目录
sp_path <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-1913-0017-A1/AXB_1913_Cellbin_Result/outs/"

# 后面不需要动，直接运行
luad_sp <- Load10X_Spatial(
  data.dir = sp_path,
  bin.size = "polygons",
  image.name = "tissue_hires_image.png"
)

# 提取表达矩阵
sp_counts <- luad_sp[["Spatial.Polygons"]]$counts

# 提取原始多边形坐标
coords_raw <- GetTissueCoordinates(luad_sp)

# 利用 aggregate 计算多边形的质心 (Centroid)
if ("cell" %in% colnames(coords_raw)) {
  coords <- aggregate(cbind(x, y) ~ cell, data = coords_raw, FUN = mean)
  rownames(coords) <- coords$cell
} else if ("barcode" %in% colnames(coords_raw)) {
  coords <- aggregate(cbind(x, y) ~ barcode, data = coords_raw, FUN = mean)
  rownames(coords) <- coords$barcode
} else {
  coords <- coords_raw
}

# 提取最终的 x, y 列构建坐标矩阵
coords <- coords[, c("x", "y")]
nUMI <- colSums(sp_counts)

# 构建 RCTD 的 Query 对象
query <- SpatialRNA(coords, sp_counts, nUMI)





# =======================================================
# 3. 创建并运行 RCTD
# =======================================================
cat("3. 正在初始化并运行 RCTD (UMI_min = 0, Doublet 模式)...\n")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)

# UMI_min = 0 确保所有细胞均参与计算
myRCTD <- create.RCTD(query, reference, max_cores = 24, UMI_min = 0)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

print("RCTD 反卷积运行完成！")

# 保存 RCTD 对象 (存到上级目录，避免弄脏原始数据文件夹)
save_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-1913-0017-A1/"
saveRDS(myRCTD, file = paste0(save_dir, 'luad_Cellbin_RCTD.RDS'))

# =======================================================
# 4. 提取保存 RCTD 结果 (包含 first_type 与 百分比矩阵)
# =======================================================
print("4. 正在提取并处理 RCTD 结果...")

# A. 保存详细分类结果 (包含 first_type, spot_class 等)
results_df <- myRCTD@results$results_df
write.csv(results_df, file = paste0(save_dir, "RCTD_results_first_type.csv"), quote = FALSE)

# B. 提取权重矩阵并归一化为 0-1 之间的百分比 (Percentage)
weights <- myRCTD@results$weights
norm_weights <- sweep(weights, 1, rowSums(weights), "/")
write.csv(as.data.frame(norm_weights), file = paste0(save_dir, "RCTD_weights_percentage.csv"), quote = FALSE)

print("🎉 细胞切割分类结果与细胞比例矩阵已全部保存！")


# 载入所需的包
library(Seurat)
library(spacexr)
library(Matrix)

# =======================================================
# 1. 准备单细胞参考集 (Reference)
# =======================================================
print("1. 正在读取单细胞参考数据集...")
# 指向导出的 13 大类数据文件夹
ref_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/sp/"

genes <- read.csv(paste0(ref_dir, "reference_genes.csv"), header=FALSE)$V1
cells <- read.csv(paste0(ref_dir, "reference_cells.csv"), header=FALSE)$V1

counts <- readMM(paste0(ref_dir, "reference_counts.mtx"))
counts <- as(counts, "CsparseMatrix")
counts <- round(counts)

rownames(counts) <- genes
colnames(counts) <- cells

meta <- read.csv(paste0(ref_dir, "reference_meta.csv"), row.names=1)
# 清理细胞类型名称中的斜杠和空格
meta$RCTD_annotation <- gsub("/", "_", meta$RCTD_annotation)
meta$RCTD_annotation <- gsub(" ", "_", meta$RCTD_annotation)

cell_types <- as.factor(meta$RCTD_annotation)
names(cell_types) <- rownames(meta)

reference <- Reference(counts, cell_types)

# =======================================================
# 2. 读取 HD 细胞切割数据并提取 Query
# =======================================================
print("2. 正在读取 Visium HD 细胞切割数据...")

# 路径指向新样本的 outs/ 目录
sp_path <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/AXB_6123_Cellbin_Result/outs/"

# 直接运行，读取多边形数据
luad_sp <- Load10X_Spatial(
  data.dir = sp_path,
  bin.size = "polygons",
  image.name = "tissue_hires_image.png"
)

# 提取表达矩阵 (自动适配 Spatial 或 Spatial.Polygons 的命名)
assay_name <- DefaultAssay(luad_sp)
sp_counts <- luad_sp[[assay_name]]$counts

# 提取原始多边形坐标
coords_raw <- GetTissueCoordinates(luad_sp)

# 利用 aggregate 计算多边形的质心 (Centroid)
if ("cell" %in% colnames(coords_raw)) {
  coords <- aggregate(cbind(x, y) ~ cell, data = coords_raw, FUN = mean)
  rownames(coords) <- coords$cell
} else if ("barcode" %in% colnames(coords_raw)) {
  coords <- aggregate(cbind(x, y) ~ barcode, data = coords_raw, FUN = mean)
  rownames(coords) <- coords$barcode
} else {
  coords <- coords_raw
}

# 提取最终的 x, y 列构建坐标矩阵
coords <- coords[, c("x", "y")]
nUMI <- colSums(sp_counts)

# 构建 RCTD 的 Query 对象
query <- SpatialRNA(coords, sp_counts, nUMI)

# =======================================================
# 3. 创建并运行 RCTD
# =======================================================
cat("3. 正在初始化并运行 RCTD (UMI_min = 0, Doublet 模式)...\n")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)

# UMI_min = 0 确保所有细胞均参与计算
myRCTD <- create.RCTD(query, reference, max_cores = 24, UMI_min = 0)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

print("RCTD 反卷积运行完成！")

# 保存 RCTD 对象到对应的样本分析目录
save_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/"
saveRDS(myRCTD, file = paste0(save_dir, 'AXB_6123_Cellbin_RCTD.RDS'))

# =======================================================
# 4. 提取保存 RCTD 结果 (包含 first_type 与 百分比矩阵)
# =======================================================
print("4. 正在提取并处理 RCTD 结果...")

# A. 保存详细分类结果 (包含 first_type, spot_class 等)
results_df <- myRCTD@results$results_df
write.csv(results_df, file = paste0(save_dir, "RCTD_results_first_type.csv"), quote = FALSE)

# B. 提取权重矩阵并归一化为 0-1 之间的百分比 (Percentage)
weights <- myRCTD@results$weights
norm_weights <- sweep(weights, 1, rowSums(weights), "/")
write.csv(as.data.frame(norm_weights), file = paste0(save_dir, "RCTD_weights_percentage.csv"), quote = FALSE)

print("🎉 细胞切割分类结果与细胞比例矩阵已全部保存！")






#fig1bc
setwd('D:/deskup/dk/sp/')
list.files()




# 1. 设置工作目录并读取数据
setwd('D:/deskup/dk/sp/')
res <- read.csv("result_LUSC_vs_LUAD.edgeR.lrt.csv", stringsAsFactors = FALSE)

# 2. 计算缺失的 FDR (利用 BH 也就是 FDR 方法对 PValue 进行多重假设检验校正)
res$FDR <- p.adjust(res$PValue, method = "BH")

# 3. 科学去重：如果有重复的基因名，保留 PValue 最小（最显著）的那一个
res <- res[!is.na(res$gene) & res$gene != "", ] # 剔除基因名为空的行
res <- res[order(res$PValue), ]                 # 按 PValue 从小到大排序
res <- res[!duplicated(res$gene), ]             # 去除重复项，自动保留最显著的
rownames(res) <- res$gene                       # 现在可以安全地将基因名设为行名了

# 4. 根据 Fang 等人的标准进行严格筛选：绝对 logFC > 2 且 FDR < 0.001
sig_res <- res[!is.na(res$FDR) & !is.na(res$logFC) & 
                 abs(res$logFC) > 2 & res$FDR < 0.001, ]

# 5. 构建肺鳞癌 (LUSC) 和肺腺癌 (LUAD) 的基因特征集
# logFC > 2 代表在 LUSC 中高表达
LUSC_signature <- sig_res$gene[sig_res$logFC > 2]

# logFC < -2 代表在 LUAD 中高表达
LUAD_signature <- sig_res$gene[sig_res$logFC < -2]

# 6. 统计并输出结果
cat("筛选完成！\n")
cat("LUSC (鳞癌) 特征基因数量:", length(LUSC_signature), "\n")
cat("LUAD (腺癌) 特征基因数量:", length(LUAD_signature), "\n")

# 7. 保存结果，方便后续单细胞或空间分析打分
save(LUSC_signature, LUAD_signature, file = "LUAD_LUSC_signatures.Rdata")
write.table(LUSC_signature, file = "LUSC_signature_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(LUAD_signature, file = "LUAD_signature_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


library(Seurat)
library(UCell)
library(BiocParallel)

print("1. 读取 MTX 纯文本矩阵，彻底绕过底层依赖包...")

# 指向截图中确切的 MTX 文件夹路径
mtx_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/AXB_6123_Cellbin_Result/outs/segmented_outputs/filtered_feature_cell_matrix"

# 使用 Read10X 直接读取文件夹
counts <- Read10X(data.dir = mtx_dir)
luad_sp <- CreateSeuratObject(counts = counts)

print("2. 正在执行标准化...")
luad_sp <- NormalizeData(luad_sp)

print("3. 准备基因集并启动 96 线程 UCell 计算...")

signatures <- list(
  ADC_score = adc_genes,
  SCC_score = scc_genes
)

luad_sp <- AddModuleScore_UCell(
  obj = luad_sp,
  features = signatures,
  ncores = 96,
  name = ""
)

print("4. 正在提取分数矩阵并保存为 CSV...")

ucell_matrix <- luad_sp@meta.data[, names(signatures)]

save_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/"
out_csv <- paste0(save_dir, "AXB_6123_UCell_ADC_SCC_scores.csv")
write.csv(ucell_matrix, file = out_csv, quote = FALSE)

print(paste0("运行完毕！分数矩阵已保存至: ", out_csv))





#fig
library(Seurat)
library(SPARK)
library(Matrix)
library(dplyr)
library(hdWGCNA)

# =======================================================
# 0. 路径与环境初始化 (针对 AXB-6123)
# =======================================================
sample_name <- "AXB_6123"
data_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/AXB_6123_Cellbin_Result/outs/"
spark_save_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/SPARK_Results/"
hdwgcna_save_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/hdWGCNA_Results/"

dir.create(spark_save_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(hdwgcna_save_dir, recursive = TRUE, showWarnings = FALSE)

# =======================================================
# 1. 运行 SPARK-X 获取 Top 500 空间高可变基因
# =======================================================
message("\n============ 1. 开始 SPARK-X 分析 ============")
seurat_obj <- Load10X_Spatial(data.dir = data_dir, bin.size = "polygons", image.name = "tissue_hires_image.png")
seurat_obj$orig.ident <- sample_name

sp_counts <- seurat_obj[["Spatial.Polygons"]]$counts
coords <- GetTissueCoordinates(seurat_obj)

# 提取细胞质心坐标
if ("cell" %in% colnames(coords) && nrow(coords) > ncol(sp_counts)) {
  coords <- coords %>% group_by(cell) %>% summarize(x = mean(x), y = mean(y)) %>% as.data.frame()
  rownames(coords) <- coords$cell
} else if ("cell" %in% colnames(coords)) {
  rownames(coords) <- coords$cell
}

locs <- as.matrix(coords[, c("x", "y")])
common_cells <- intersect(colnames(sp_counts), rownames(locs))
sp_counts <- sp_counts[, common_cells]
locs <- locs[common_cells, ]

# 剔除全零基因并运行 SPARK-X
sp_counts_filtered <- sp_counts[rowSums(sp_counts) > 0, ]
sparkx_res <- sparkx(sp_counts_filtered, locs, numCores = 10, option = "mixture")

res <- as.data.frame(sparkx_res$res_mtest)
rownames(res) <- rownames(sp_counts_filtered)

# 剔除线粒体基因并提取 Top 500
mt_idx <- grep("^MT-|^mt-", rownames(res))
if (length(mt_idx) > 0) res <- res[-mt_idx, ]

res_clean <- res[res$combinedPval < 1, ]
res_sorted <- res_clean[order(res_clean$adjustedPval), ]
AXB_6123_svgs <- head(res_sorted, min(500, nrow(res_sorted)))

# 落地保存 SPARK-X 结果
save(AXB_6123_svgs, file = paste0(spark_save_dir, sample_name, "_top500_SVGs.Rdata"))
rm(seurat_obj, sparkx_res, res, res_clean, res_sorted, sp_counts_filtered); gc()

# =======================================================
# 2. 构建空间元细胞 (Metacell Pipeline)
# =======================================================
message("\n============ 2. 开始构建空间元细胞 ============")
obj <- Load10X_Spatial(data.dir = data_dir, bin.size = "polygons")
obj$orig.ident <- sample_name

# 基础降维聚类
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.8, verbose = FALSE)

# 对齐 SPARK-X 提取的高可变基因
valid_genes <- intersect(rownames(AXB_6123_svgs), rownames(obj))

# 物理坐标降维构建
coords <- GetTissueCoordinates(obj)
locs <- as.matrix(coords[, c("x", "y")])
if ("cell" %in% colnames(coords)) {
  rownames(locs) <- coords$cell
} else if (!is.null(rownames(coords))) {
  rownames(locs) <- rownames(coords)
}
colnames(locs) <- c("spcoord_1", "spcoord_2")
obj[["SPATIAL"]] <- CreateDimReducObject(embeddings = locs, key = "spcoord_", assay = "Spatial.Polygons")

# hdWGCNA 初始化与元细胞构建
obj <- SetupForWGCNA(obj, gene_select = "custom", gene_list = valid_genes, wgcna_name = sample_name)
obj <- MetacellsByGroups(
  seurat_obj = obj, group.by = "seurat_clusters", ident.group = "seurat_clusters",
  reduction = "SPATIAL", k = 25, max_shared = 10, assay = "Spatial.Polygons"
)
obj <- NormalizeMetacells(obj)

# 提取表达矩阵
mc_obj <- GetMetacellObject(obj)
valid_clusters <- as.character(na.omit(unique(mc_obj$seurat_clusters)))
obj <- SetDatExpr(obj, group_name = valid_clusters, group.by = "seurat_clusters", assay = "Spatial.Polygons")

# =======================================================
# 3. 落地保存对象，准备对接后续分析
# =======================================================
out_rdata <- paste0(hdwgcna_save_dir, sample_name, "_Cellbin_Metacell.RData")
save(obj, file = out_rdata)
message("🎉 完毕！已成功保存至: ", out_rdata)






#fig4
library(Seurat)
library(UCell)
# 如果有 MEscnet 专门的包也要记得 library

# ==========================================
# 0. 载入之前算好的 hdWGCNA 元细胞数据
# ==========================================
sample_name <- "AXB_6123"
rdata_path <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/hdWGCNA_Results/AXB_6123_Cellbin_Metacell.RData"

# 这一步 load 进来后，你的环境变量里会多出一个叫做 'obj' 的对象
load(rdata_path) 


obj_6123<-obj
obj_6123@misc$MEscnet$datExpr<-obj_6123@misc$AXB_6123$datExpr

obj_6123<-pick_power_scale_free_power(obj_6123)

obj_6123<-ComputeMEscnetModules(obj_6123,min_genes_per_module = 0,resolution = 0.7,number = 2)
modules<-obj_6123@misc$MEscnet$MEscnet_modules$module_ids
module_genes<-obj_6123@misc$MEscnet$MEscnet_modules$gene_lists




seurat_obj<-obj_6123
datExpr <- seurat_obj@misc$MEscnet$datExpr
moduleColors <- seurat_obj@misc$MEscnet$moduleColors
################################
# 2. 计算module eigengene
################################

MEs0 <- moduleEigengenes(
  datExpr,
  colors = moduleColors
)$eigengenes


MEs <- orderMEs(MEs0)



################################
# 3. 计算kME
################################

datKME <- signedKME(
  datExpr,
  MEs,
  outputColumnName="kME_MM."
)


################################
# 4. 获取module名字
################################

modules <- unique(moduleColors)
modules <- modules[modules!="grey"]


modules
length(modules)



################################
# 5. 生成 module-gene-kME list
################################
module_kME_list <- list()
for(mod in modules){
  # 当前module基因
  genes <- names(moduleColors[
    moduleColors==mod
  ])
  
  
  # kME列名
  
  kcol <- paste0(
    "kME_MM.",
    mod
  )
  tmp <- data.frame(
    gene=genes,
    kME=datKME[
      genes,
      kcol
    ],
    stringsAsFactors = FALSE
  )
  # kME排序
  tmp <- tmp[
    order(
      tmp$kME,
      decreasing=TRUE
    ),
  ]
  rownames(tmp)<-NULL
  module_kME_list[[mod]] <- tmp
}


library(dplyr)
library(UCell)
library(Seurat)

cat("⏳ 正在根据 kME > 0.6 过滤模块基因集...\n")

# ==============================================================================
# 1. 提取并过滤每个 Module 中 kME > 0.6 的基因
# ==============================================================================
# 使用 lapply 遍历列表，筛选出符合条件的基因并转化为纯字符向量列表
filtered_signatures <- lapply(module_kME_list, function(df) {
  df %>%
    filter(kME > 0.6) %>%
    pull(gene)
})

# 简单检查一下过滤后的基因集大小
cat("✅ 过滤完成！各模块保留的强连通基因数量如下：\n")
print(lengths(filtered_signatures))

# ==============================================================================
# 2. 使用 UCell 并行计算模块得分 (启用 24 线程)
# ==============================================================================
cat("⏳ 正在启动 UCell 多线程打分引擎 (ncores = 24)...\n")

# 注意：seurat_obj 需要是你当前环境中已经载入的 Seurat 对象
seurat_obj <- UCell::AddModuleScore_UCell(
  obj = seurat_obj,
  features = filtered_signatures,
  ncores = 96,
  name = "_UCell" # UCell 会自动在模块名后缀加上 "_UCell" 以防与原基因名冲突
)

cat("🎉 UCell 评分计算完毕！结果已无缝添加至 seurat_obj@meta.data 中。\n")


# ==========================================
# 1. 提取 UCell 分数并重命名为 M1-M7
# ==========================================
cat("⏳ 1. 正在提取并重命名 Module 分数...\n")

# 获取带有 Module_X_UCell 的列
mod_cols <- grep("Module_[0-9]+_UCell", colnames(seurat_obj@meta.data), value = TRUE)
scores_df <- seurat_obj@meta.data[, mod_cols]

# 将 Module_1_UCell 替换为 M1，以此类推
colnames(scores_df) <- gsub("Module_([0-9]+)_UCell", "M\\1", colnames(scores_df))

# ==========================================
# 2. 内部 Scale 与 AST 分数计算
# ==========================================
cat("⚖️ 2. 正在执行 Scale 标准化与 AST 过渡态计算...\n")

# 对 M1 到 M7 分别进行 Z-score 标准化
scaled_scores <- as.data.frame(scale(scores_df))

# 核心逻辑：定义 AST = Scale(M1) - Scale(M3)
scaled_scores$AST <- scaled_scores$M1 - scaled_scores$M3

# ==========================================
# 3. 导出 CSV
# ==========================================
cat("💾 3. 正在保存所有标准化分数...\n")

save_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/"
out_csv <- paste0(save_dir, "AXB_6123_Scaled_Modules_AST.csv")

# 导出，行名即为 cell_id
write.csv(scaled_scores, file = out_csv, quote = FALSE)
cat("🎉 搞定！Scaled 分数与 AST 分数已保存至:", out_csv, "\n")






#fig5
# ==========================================
# 0. 载入必需包
# ==========================================
library(ComplexHeatmap)
library(circlize)

cat("⏳ [Heatmap] 正在读取数据并合并标签与 UCell 分数...\n")

# 1. 读取 Python 导出的轨迹 CSV (包含 ptime)
local_csv_path <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/AXB_6123_Trajectory_Export_for_Local_2.csv"
pt_data <- read.csv(local_csv_path)

# 2. 从本地读取 CNV 和 Scissor 结果，重新构建 Tumor_Label
pt_data$join_id <- as.character(as.numeric(gsub("[^0-9]", "", pt_data$cell_id)))

cnv_files <- Sys.glob("/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/InferCNV_Results/*Cancer_Cells_Only_With_UCell.csv")
if (length(cnv_files) > 0) {
  cnv_df <- read.csv(cnv_files[1])
  cnv_df$join_id <- as.character(as.numeric(gsub("[^0-9]", "", cnv_df[,1])))
  pt_data <- merge(pt_data, cnv_df[, c("join_id", "malignant_status")], by = "join_id", all.x = TRUE)
}

scissor_df <- read.csv("/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/AXB_6123_Scissor_labels.csv")
scissor_df$join_id <- as.character(as.numeric(gsub("[^0-9]", "", scissor_df[,1])))
colnames(scissor_df)[2] <- "Scissor_Subtype"
pt_data <- merge(pt_data, scissor_df[, c("join_id", "Scissor_Subtype")], by = "join_id", all.x = TRUE)

pt_data$Tumor_Label <- "Transitional"
pt_data$Tumor_Label[pt_data$Scissor_Subtype == "LUAD_like"] <- "LUAD+"
pt_data$Tumor_Label[pt_data$Scissor_Subtype == "LUSC_like"] <- "LUSC+"

meta_data <- seurat_obj@meta.data
meta_data$cell_id <- rownames(meta_data)
pt_data <- merge(pt_data, meta_data, by = "cell_id", all.x = TRUE)
pt_data <- pt_data[order(pt_data$ptime), ]

# ==========================================
# 1. 全部先 Scale，再计算 AST
# ==========================================
sig_cols <- paste0("Module_", 1:7, "_UCell")
pt_data <- pt_data[complete.cases(pt_data[, sig_cols]), ]

for (col in sig_cols) {
  pt_data[[col]] <- scale(pt_data[[col]])[, 1]
}

pt_data$AST <- pt_data$Module_1_UCell - pt_data$Module_3_UCell

module_names <- c(sig_cols, "AST")
display_names <- c(paste0("M", 1:7), "AST")

# ==========================================
# 2. Loess 沿轨迹平滑拟合
# ==========================================
pt_seq <- seq(min(pt_data$ptime, na.rm = TRUE), max(pt_data$ptime, na.rm = TRUE), length.out = 100)

smooth_matrix <- matrix(NA, nrow = length(module_names), ncol = 100)
rownames(smooth_matrix) <- display_names

for (i in seq_along(module_names)) {
  m_name <- module_names[i]
  form <- as.formula(paste0(m_name, " ~ ptime"))
  fit <- loess(form, data = pt_data, span = 0.75)
  smooth_matrix[i, ] <- predict(fit, newdata = data.frame(ptime = pt_seq))
}

scaled_matrix <- t(apply(smooth_matrix, 1, scale))
colnames(scaled_matrix) <- paste0("Bin", 1:100)

custom_row_order <- c("M5", "M3", "M6", "M7", "M1", "M4", "M2", "AST")
scaled_matrix <- scaled_matrix[custom_row_order, ]

# ==========================================
# 3. 计算三层分层山脊图的坐标参数
# ==========================================
cat("🎨 [Heatmap] 正在计算伪时间分层山脊分布...\n")

dens_luad <- density(pt_data$ptime[pt_data$Tumor_Label == "LUAD+"])
dens_trans <- density(pt_data$ptime[pt_data$Tumor_Label == "Transitional"])
dens_lusc <- density(pt_data$ptime[pt_data$Tumor_Label == "LUSC+"])

y_luad <- approx(dens_luad$x, dens_luad$y, xout = pt_seq, rule = 2)$y
y_trans <- approx(dens_trans$x, dens_trans$y, xout = pt_seq, rule = 2)$y
y_lusc <- approx(dens_lusc$x, dens_lusc$y, xout = pt_seq, rule = 2)$y

y_luad <- (y_luad / max(y_luad, na.rm = TRUE)) * 0.9
y_trans <- (y_trans / max(y_trans, na.rm = TRUE)) * 0.35
y_lusc <- (y_lusc / max(y_lusc, na.rm = TRUE)) * 0.9

font_style <- gpar(fontsize = 14, col = "black", fontface = "plain", fontfamily = "sans")

top_anno <- HeatmapAnnotation(
  Density = anno_empty(height = unit(4, "cm"), border = FALSE)
)

bottom_anno <- HeatmapAnnotation(
  axis = anno_empty(height = unit(1, "cm"), border = FALSE)
)

# ==========================================
# 4. 渲染无白缝的平滑渐变热图
# ==========================================
cat("🎨 [Heatmap] 正在渲染无缝隙瀑布流热图...\n")

col_fun <- colorRamp2(c(-2, 0, 2), c("#009ad6", "white", "#d71345"))

ht <- Heatmap(
  scaled_matrix,
  name = "Scaled\nExpr",
  col = col_fun,
  top_annotation = top_anno,
  bottom_annotation = bottom_anno,
  use_raster = TRUE,          # 开启热图主体栅格化
  raster_quality = 5,         # 提高栅格化分辨率
  rect_gp = gpar(col = NA),   # 去除单元格边框线
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 1),
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_title = NULL,
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  row_order = custom_row_order,
  row_title_rot = 0,
  row_title_side = "right",
  row_names_gp = font_style,
  heatmap_legend_param = list(
    title_gp = font_style,
    labels_gp = font_style,
    legend_height = unit(4, "cm"),
    border = "black"
  )
)

# 输出 PDF
cairo_pdf("/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/Heatmap_AST_Trajectory_Seamless_14pt.pdf",
          width = 7, height = 7, family = "sans")
draw(ht)


# ----------------------------------------------------
# 0. 锁定全局字体 (确保与前面 Python 出图完全一致)
# ----------------------------------------------------
font_style <- gpar(fontsize = 14, col = "black", fontface = "plain", fontfamily = "Arial")

# ----------------------------------------------------
# 绘制三层错落的山脊图
# ----------------------------------------------------
decorate_annotation("Density", {
  pushViewport(viewport(xscale = c(0.5, 100.5), yscale = c(0, 3.2)))
  
  # 1. 顶层: LUAD (Base Y = 2) - 深蓝色
  grid.polygon(
    x = c(1:100, 100, 1),
    y = c(2 + y_luad, 2, 2),
    default.units = "native",
    gp = gpar(fill = "#3C548880", col = "#3C5488", lwd = 1.5)
  )
  grid.lines(
    x = c(1, 100),
    y = c(2, 2),
    default.units = "native",
    gp = gpar(col = "#3C5488", lwd = 1.5)
  )
  max_idx_luad <- which.max(y_luad)
  grid.text(
    "LUAD cell",
    x = max_idx_luad,
    y = 2 + y_luad[max_idx_luad] - 0.2,
    default.units = "native",
    just = "center",
    gp = font_style
  )
  
  # 2. 中层: Transitional (Base Y = 1) - 中性灰
  grid.polygon(
    x = c(1:100, 100, 1),
    y = c(1 + y_trans, 1, 1),
    default.units = "native",
    gp = gpar(fill = "#B0B0B080", col = "#B0B0B0", lwd = 1.5)
  )
  grid.lines(
    x = c(1, 100),
    y = c(1, 1),
    default.units = "native",
    gp = gpar(col = "#B0B0B0", lwd = 1.5)
  )
  max_idx_trans <- which.max(y_trans)
  grid.text(
    "Transitional cell",
    x = max_idx_trans,
    y = 1 + y_trans[max_idx_trans] + 0.2,
    default.units = "native",
    just = "center",
    gp = font_style
  )
  
  # 3. 底层: LUSC (Base Y = 0) - 质感红
  grid.polygon(
    x = c(1:100, 100, 1),
    y = c(0 + y_lusc, 0, 0),
    default.units = "native",
    gp = gpar(fill = "#E64B3580", col = "#E64B35", lwd = 1.5)
  )
  grid.lines(
    x = c(1, 100),
    y = c(0, 0),
    default.units = "native",
    gp = gpar(col = "#E64B35", lwd = 1.5)
  )
  max_idx_lusc <- which.max(y_lusc)
  grid.text(
    "LUSC cell",
    x = max_idx_lusc,
    y = 0 + y_lusc[max_idx_lusc] - 0.2,
    default.units = "native",
    just = "center",
    gp = font_style
  )
  
  popViewport()
})
# ----------------------------------------------------
# 底部标注 X 轴
# ----------------------------------------------------
decorate_annotation("axis", {
  grid.text("LUAD", x = 0, y = 0.5, just = "left", gp = font_style)
  grid.text("Pseudotime", x = 0.5, y = 0.5, just = "center", gp = font_style)
  grid.text("LUSC", x = 1, y = 0.5, just = "right", gp = font_style)
})

dev.off()

cat("✅ [Heatmap] 大功告成！无缝隙、三层山脊密度的热图已保存为 Heatmap_AST_Trajectory_Seamless_14pt.pdf\n")



#fig6
# ==========================================
# 0. 载入必需包
# ==========================================
library(ggplot2)
library(dplyr)
library(tidyr)

cat("⏳ [Curve] 正在读取数据、合并 UCell 分数并进行标准化处理...\n")

# 1. 读取 Python 导出的轨迹 CSV (包含 ptime 和 cell_id)
local_csv_path <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/AXB_6123_Trajectory_Export_for_Local_2.csv"
pt_data <- read.csv(local_csv_path)

# 2. 从全量 Seurat 对象中提取 metadata 注入 UCell 分数
meta_data <- seurat_obj@meta.data
meta_data$cell_id <- rownames(meta_data)

# 通过 cell_id 取交集，把 UCell 矩阵和细胞状态带进来
pt_data <- merge(pt_data, meta_data, by = "cell_id", all.x = TRUE)
pt_data <- pt_data[order(pt_data$ptime), ]

# ==========================================
# 1. 核心数学处理：提取 Scale 后的分数并计算 AST
# ==========================================
pt_data$M1_Scaled <- scale(pt_data$Module_1_UCell)[, 1]
pt_data$M3_Scaled <- scale(pt_data$Module_3_UCell)[, 1]

# 计算 AST 分数: M1 (LUSC) - M3 (LUAD)
pt_data$AST_Score <- pt_data$M1_Scaled - pt_data$M3_Scaled

# ==========================================
# 2. 🌟 核心修改：寻找一阶导数最大的"精确转折时间点"
# ==========================================
fit_ast <- loess(AST_Score ~ ptime, data = pt_data, span = 0.75)
seq_pt <- seq(min(pt_data$ptime, na.rm = TRUE), max(pt_data$ptime, na.rm = TRUE), length.out = 2000)
pred_ast <- predict(fit_ast, newdata = data.frame(ptime = seq_pt))

# 计算一阶导数 (变化速率)
deriv_ast <- diff(pred_ast) / diff(seq_pt)
deriv_ast <- c(deriv_ast[1], deriv_ast) # 补齐长度以匹配 seq_pt

# 转折点定义为一阶导数最大值对应的 ptime
turn_ptime <- seq_pt[which.max(deriv_ast)]

cat(paste0("📌 发现最大变化速率转折点：Pseudotime = ", round(turn_ptime, 4), "\n"))

# ==========================================
# 3. 准备绘图颜色与排版参数 (严格锁定 14pt Arial 纯黑)
# ==========================================
state_colors <- c("LUAD+" = "#3C5488", "LUSC+" = "#E64B35", "Transitional" = "#B0B0B0")

font_settings <- theme(
  text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
  axis.text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
  axis.title = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
  legend.text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
  legend.title = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
  legend.position = "right",
  panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
  panel.grid.major = element_line(color = "#E5E5E5", linetype = "dotted"),
  panel.grid.minor = element_blank()
)

# ==========================================
# 4. 使用 ggplot2 渲染动态演化曲线
# ==========================================
cat("🎨 [Curve] 正在渲染 AST 动态演化曲线...\n")

p <- ggplot(pt_data, aes(x = ptime)) +
  # 1. 绘制底层细胞散点
  geom_point(aes(y = AST_Score, color = Tumor_Label), alpha = 0.4, size = 1.5) +
  scale_color_manual(values = state_colors, name = "Cell State (Scatter)") +
  
  # 2. 绘制 M3 (LUAD) 的衰退拟合曲线
  geom_smooth(aes(y = M3_Scaled, linetype = "M3 (LUAD) Trend"), 
              color = "#3C5488", method = "loess", span = 0.75, se = FALSE, linewidth = 1.5) +
  
  # 3. 绘制 M1 (LUSC) 的崛起拟合曲线
  geom_smooth(aes(y = M1_Scaled, linetype = "M1 (LUSC) Trend"), 
              color = "#E64B35", method = "loess", span = 0.75, se = FALSE, linewidth = 1.5) +
  
  # 4. 绘制 AST (M1 - M3) 的综合演化曲线
  geom_smooth(aes(y = AST_Score, linetype = "AST Score (M1 - M3)"), 
              color = "black", method = "loess", span = 0.75, se = FALSE, linewidth = 1.2) +
  
  # 5. 设定线条图例
  scale_linetype_manual(
    name = "Fitted Curves",
    values = c("M3 (LUAD) Trend" = "solid", 
               "M1 (LUSC) Trend" = "solid", 
               "AST Score (M1 - M3)" = "dashed")
  ) +
  
  # 6. 添加横向完美平衡基准线 (Y=0)
  geom_hline(yintercept = 0, color = "black", linetype = "dotted", linewidth = 0.8) +
  
  # 7. 添加"转折时间点"竖向虚线及标注
  geom_vline(xintercept = turn_ptime, color = "black", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = turn_ptime, y = max(pt_data$AST_Score, na.rm = TRUE) * 0.85,
           label = "转折时间点", angle = 90, vjust = -0.8,
           size = 14 / .pt, color = "black", family = "Arial") +
  
  # 8. 标签与主题
  labs(
    title = "Dynamics of Adenosquamous Transition (AST)",
    x = "Pseudotime",
    y = "Scaled Module Score (Z-score)"
  ) +
  theme_classic() +
  font_settings

# ==========================================
# 5. 使用 cairo_pdf 输出高精度 PDF
# ==========================================
cairo_pdf("/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/Curve_AST_Dynamics_Strict_14pt.pdf",
          width = 9, height = 6, family = "Arial")
print(p)
dev.off()

cat("✅ [Curve] 完美！带有十字准星的演化曲线已保存为 Curve_AST_Dynamics_Strict_14pt.pdf\n")

