
#fi7
# ==========================================
# 0. 加载必需包
# ==========================================
library(Seurat)
library(dplyr)
library(UCell)
library(randomForest)
library(tibble)
library(ggplot2)

cat("⏳ 1. 正在读取 Zoom 区域拟时序数据并 Subset 空间对象...\n")

# ==========================================
# 1. 读取拟时序 CSV 并 Subset Seurat 对象
# ==========================================
csv_path <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1/AXB-6123-2431-6976-A1_Zoom_Trajectory_Data_with_ID.csv"
zoom_df <- read.csv(csv_path)

valid_cells <- intersect(zoom_df$cell_id, colnames(obj_6123))
cat(paste("🎯 成功匹配到的 Zoom 区域细胞数：", length(valid_cells), "\n"))

obj_zoom <- subset(obj_6123, cells = valid_cells)

rownames(zoom_df) <- zoom_df$cell_id
obj_zoom <- AddMetaData(obj_zoom, metadata = zoom_df[, c("ptime", "Tumor_Label")])

# ==========================================
# 2. 提取靶向 TF 表达矩阵并构建随机森林
# ==========================================
cat("⏳ 2. 正在提取表达矩阵并训练随机森林模型...\n")

target_tfs <- c("TP63", "BCL6", "SOX2", "NFE2L2", "EHF", "ELF3", "CREB3L1", "SOX9")
valid_tfs <- intersect(target_tfs, rownames(obj_zoom))

expr_matrix <- t(as.matrix(GetAssayData(obj_zoom, assay = "Spatial.Polygons", layer = "data")[valid_tfs, colnames(obj_zoom)]))

rf_data <- data.frame(expr_matrix)
rf_data$Target_AST <- obj_zoom@meta.data$ptime

set.seed(42)
rf_model <- randomForest(
  Target_AST ~ ., 
  data = rf_data, 
  ntree = 500, 
  importance = TRUE
)

# ==========================================
# 3. 计算重要性与皮尔逊相关性
# ==========================================
importance_df <- as.data.frame(importance(rf_model))
importance_df <- tibble::rownames_to_column(importance_df, var = "Gene")
colnames(importance_df)[which(colnames(importance_df) == "%IncMSE")] <- "RF_Importance"

cor_results <- sapply(valid_tfs, function(g) {
  cor(rf_data[[g]], rf_data$Target_AST, method = "pearson")
})

importance_df$Pearson_Cor <- cor_results[importance_df$Gene]

importance_df <- importance_df %>%
  mutate(
    Direction = case_when(
      Pearson_Cor > 0 ~ "Positive Driver (Push to LUSC)",
      Pearson_Cor < 0 ~ "Negative Driver (Anchor in LUAD)",
      TRUE ~ "Neutral"
    )
  )

# 🌟 新增：在控制台直接打印按重要性降序的排行榜
cat("\n🏆 ========================================= 🏆\n")
cat("          转录因子重要性排行榜 (Top Drivers)         \n")
cat("🏆 ========================================= 🏆\n")
ranked_df <- importance_df %>% arrange(desc(RF_Importance))
print(ranked_df[, c("Gene", "RF_Importance", "Pearson_Cor", "Direction")])
cat("\n=============================================\n")

# 按相关性排序用于画图排列
importance_df <- importance_df %>% arrange(Pearson_Cor) 
importance_df$Gene <- factor(importance_df$Gene, levels = importance_df$Gene)

# ==========================================
# 4. 绘制棒棒糖图 (严格锁定 14pt 黑色不加粗)
# ==========================================
cat("🎨 3. 正在绘制 TF 驱动力棒棒糖图...\n")

driver_colors <- c(
  "Positive Driver (Push to LUSC)" = "#E3493B", 
  "Negative Driver (Anchor in LUAD)" = "#3CB4D2"
)

p_lollipop <- ggplot(importance_df, aes(x = Pearson_Cor, y = Gene, color = Direction)) +
  geom_segment(aes(x = 0, xend = Pearson_Cor, y = Gene, yend = Gene), linewidth = 1) +
  geom_point(aes(size = RF_Importance)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_color_manual(values = driver_colors) +
  scale_size_continuous(range = c(3, 10), name = "RF Importance") +
  labs(
    x = "Pearson Correlation (TF Expression vs ptime)", 
    y = NULL,
    title = "AST Transition Drivers"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 14, color = "black", family = "sans", face = "plain"),
    axis.text = element_text(size = 14, color = "black", family = "sans", face = "plain"),
    axis.title = element_text(size = 14, color = "black", family = "sans", face = "plain"),
    plot.title = element_text(size = 14, color = "black", family = "sans", face = "plain", hjust = 0.5),
    legend.text = element_text(size = 14, color = "black", family = "sans", face = "plain"),
    legend.title = element_text(size = 14, color = "black", family = "sans", face = "plain"),
    legend.position = "right",
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black")
  )

out_pdf <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/hdWGCNA_Results/AXB_6123_TF_Lollipop_ptime.pdf"
ggsave(out_pdf, plot = p_lollipop, width = 8, height = 6, device = cairo_pdf)

cat(paste("✅ 大功告成！PDF已保存至：", out_pdf, "\n"))




############################################################
# TP63 Virtual KO Pipeline (Part 1: Data Prep & Virtual KO)
# AXB-6123 CellBin spatial zoom region
# Raw -> Y crop -> RCTD epithelial -> InferCNV tumor -> scTenifoldKnk
############################################################

rm(list = ls())
gc()

# ==========================================
# 0. 载入必需包与路径配置
# ==========================================
library(Seurat)
library(dplyr)
library(stringr)
library(jsonlite)
library(scTenifoldKnk)

cat("⏳ [Step 0] 初始化路径与依赖库...\n")

base_dir <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/seg/AXB-6123-2431-6976-A1"

h5_path <- file.path(base_dir, "AXB_6123_Cellbin_Result", "outs", "segmented_outputs", "filtered_feature_cell_matrix.h5")
geojson_path <- file.path(base_dir, "AXB_6123_Cellbin_Result", "outs", "segmented_outputs", "cell_segmentations.geojson")

# ==========================================
# 1. 加载原始表达矩阵
# ==========================================
cat("⏳ [Step 1] 正在加载原始 H5 表达矩阵...\n")
counts <- Read10X_h5(h5_path)

if(is.list(counts)){
  counts <- counts[["Gene Expression"]]
}

cat("👉 原始全量细胞数:", ncol(counts), "\n")

# ==========================================
# 2. 读取 GeoJSON 空间坐标并重建 Barcode
# ==========================================
cat("⏳ [Step 2] 正在读取 GeoJSON 并提取空间坐标...\n")

geo <- fromJSON(geojson_path, simplifyVector = TRUE)
props <- geo$features$properties
centroid <- props$cell_centroid

if(is.matrix(centroid)){
  X <- centroid[, 1]
  Y <- centroid[, 2]
} else {
  X <- sapply(centroid, function(x) x[[1]])
  Y <- sapply(centroid, function(x) x[[2]])
}

coord_df <- data.frame(
  cell_id = paste0("cellid_", sprintf("%09d", as.numeric(props$cell_id)), "-1"),
  X_coord = as.numeric(X),
  Y_coord = as.numeric(Y),
  stringsAsFactors = FALSE
)

cat("👉 GeoJSON 解析细胞数:", nrow(coord_df), "\n")

# ==========================================
# 3. Y 轴靶区裁剪 (Zoom Region)
# ==========================================
cat("⏳ [Step 3] 正在根据 Y 轴坐标裁剪靶区 (Zoom Region)...\n")

y_max <- max(coord_df$Y_coord, na.rm = TRUE)
y_min <- min(coord_df$Y_coord, na.rm = TRUE)
y_range <- y_max - y_min

y_lower <- y_max - y_range * 0.18
y_upper <- y_max - y_range * 0.01

coord_zoom <- coord_df %>%
  filter(Y_coord > y_lower, Y_coord < y_upper)

cat("👉 靶区裁剪后细胞数:", nrow(coord_zoom), "\n")

# ==========================================
# 4. 矩阵匹配与构建靶区 Seurat 对象
# ==========================================
cat("⏳ [Step 4] 正在匹配表达矩阵并构建 Zoom Seurat 对象...\n")

zoom_cells <- intersect(colnames(counts), coord_zoom$cell_id)

if(length(zoom_cells) == 0){
  stop("❌ 致命错误：细胞 Barcode 未能成功匹配！")
}

counts_zoom <- counts[, zoom_cells]
coord_zoom <- coord_zoom[match(zoom_cells, coord_zoom$cell_id), ]
rownames(coord_zoom) <- coord_zoom$cell_id

obj_zoom <- CreateSeuratObject(counts = counts_zoom, meta.data = coord_zoom)

cat("👉 Zoom Seurat 对象构建完毕，包含细胞数:", ncol(obj_zoom), "\n")

# ==========================================
# 5. RCTD 标签过滤 (仅保留上皮细胞)
# ==========================================
cat("⏳ [Step 5] 正在加载 RCTD 结果并过滤上皮细胞 (Epithelial)...\n")

rctd <- read.csv(file.path(base_dir, "RCTD_results_first_type.csv"), 
                 stringsAsFactors = FALSE, check.names = FALSE)
colnames(rctd)[1] <- "cell_id"
rctd <- rctd[, c("cell_id", "first_type")]

meta <- obj_zoom@meta.data
meta$cell_id <- rownames(meta)
meta <- left_join(meta, rctd, by = "cell_id")

epi_cells <- meta %>%
  filter(first_type == "Epithelial_cell") %>%
  pull(cell_id)

obj_epi <- subset(obj_zoom, cells = epi_cells)

cat("👉 上皮细胞过滤完毕，剩余细胞数:", length(epi_cells), "\n")

# ==========================================
# 6. InferCNV 标签过滤 (仅保留恶性肿瘤细胞)
# ==========================================
cat("⏳ [Step 6] 正在加载 InferCNV 结果并过滤恶性肿瘤细胞 (Cancer/Tumor)...\n")

infer_file <- list.files(
  file.path(base_dir, "InferCNV_Results"),
  pattern = "AXB-6123_InferCNV_malignant_status.csv",
  full.names = TRUE
)[1]

infer <- read.csv(infer_file, stringsAsFactors = FALSE, check.names = FALSE)

tumor_ids <- infer %>%
  filter(malignant_status == "Tumor") %>%
  pull(cell_id)

tumor_cells <- intersect(colnames(obj_epi), tumor_ids)
obj_tumor <- subset(obj_epi, cells = tumor_cells)

cat("👉 恶性上皮肿瘤细胞最终确认为:", length(tumor_cells), "个\n")

# ==========================================
# 7. 预处理与特征基因提取 (HVG + Module Genes)
# ==========================================
cat("⏳ [Step 7] 正在进行标准化并提取高变基因及模块特征基因...\n")

obj_tumor <- NormalizeData(obj_tumor)
obj_tumor <- FindVariableFeatures(obj_tumor, selection.method = "vst", nfeatures = 2000)
hvg_genes <- VariableFeatures(obj_tumor)

# AST 相关模块基因
M1_genes <- c("CLDN1","KRT15","KRT17","TP63","NTS","LMO4","SFN","TNFSF10","FAM43A",
              "CD9","GABRE","FGFR2","BCL6","RGMA","VTCN1","PERP","NFE2L2","LY6E",
              "TMEM44","DSC3","ANKRD62","GSTM3","ADH7","ABCC5","IGF2BP2","EPCAM",
              "KRT19","NTRK2","CES1","TNC","KREMEN1","PTPRZ1","GCLC","ACAP2",
              "DVL3","TFRC","KRT5","RNF7","ATP13A3","PCYT1A","CSTA","EPHB3",
              "NCBP2","FXYD3","PTPRF","MELTF","WNK2","TBL1XR1","GPX2","HSPB1",
              "NDUFB5","EIF4A2","EHF","BDH1","MCCC1","DDR1","DLG1","CASK",
              "TMPRSS4","EHBP1","SOX2","LSG1","CALML3","PITX1","DYNLT2B","BEX3",
              "ZNF639","IFI27","TACSTD2","C12orf75","ABCF3","TNS4","GPC3",
              "NDUFB9","CDK4","EIF2B5","HSP90AB1")

M3_genes <- c("PIGR","SERPINA1","MUC5B","DMBT1","OLFM4","TFF2","IGFBP2","LYZ",
              "FCGBP","AGR2","ATP1B1","TFF3","CA9","SERPINA3","GOLM1","PROM1",
              "CRISP3","MUC6","S100A10","PGC","LCN2","GCNT3","WFDC2","QSOX1",
              "SPINK4","MUC1","CREB3L1","MLPH","CLDN2","SLC44A4","CP","CD24",
              "CEACAM5","CTSE","GP2","LGALS4","SPINK1","MIA","GFUS","TFF1",
              "RASSF7","KRT8","S100A6","MUC5AC","DUOX2","FAM3D","FAM3B","TMC5",
              "BACE2","SLC4A4","ANXA10","ELF3","TMPRSS3","TSPAN8","AQP5","SOX9",
              "MMP1","ERN2","CYSTM1","CEACAM6","NRG1","GMDS","ANXA4")

target_genes <- unique(c(hvg_genes, M1_genes, M3_genes, "TP63"))
target_genes <- intersect(target_genes, rownames(obj_tumor))

cat("👉 参与构建 GRN 网络的总特征基因数:", length(target_genes), "\n")

if (!"TP63" %in% target_genes) {
  stop("❌ 致命错误：目标基因 TP63 不在表达矩阵中！")
}

count_mat <- as.matrix(GetAssayData(obj_tumor, assay = "RNA", layer = "counts")[target_genes, ])

# ==========================================
# 8. 矩阵清洗 (剔除免疫污染基因)
# ==========================================
cat("⏳ [Step 8] 正在扫描并清洗矩阵中的空间溢出污染基因...\n")

all_genes <- rownames(count_mat)
ig_genes <- grep("^IGH|^IGK|^IGL|JCHAIN", all_genes, value = TRUE)

if(length(ig_genes) > 0) {
  cat("🧹 检测到并剔除", length(ig_genes), "个潜在免疫污染基因。\n")
}

clean_genes <- setdiff(all_genes, ig_genes)
mtx_clean <- count_mat[clean_genes, ]

target_gene <- "TP63"

if (!target_gene %in% clean_genes) {
  stop(paste("❌ 致命错误：", target_gene, "在矩阵清洗后意外丢失！"))
}

cat("👉 清洗完成！最终进入虚拟敲除的矩阵维度:", 
    nrow(mtx_clean), "基因 x", ncol(mtx_clean), "细胞\n")

# ==========================================
# 9. 纯净版全局虚拟敲除计算 (scTenifoldKnk)
# ==========================================
cat(paste("\n🚀 [Step 9] 正在启动纯净版", target_gene, "虚拟敲除计算 (关闭 QC)...\n"))

set.seed(2026)
res_knk_clean <- scTenifoldKnk(
  countMatrix = mtx_clean,
  gKO = target_gene,
  qc = FALSE,
  nc_nNet = 10,
  nCores = 64
)

# ==========================================
# 10. 提取与保存虚拟敲除结果
# ==========================================
cat("\n⏳ [Step 10] 提取受影响最显著的 Top 20 基因并保存数据...\n")

diff_reg <- res_knk_clean$diffRegulation
top20_genes <- head(diff_reg[order(diff_reg$p.adj, decreasing = FALSE), ], 20)

print("🏆 Top 20 显著受影响基因列表：")
print(top20_genes[, c("gene", "FC", "p.value", "p.adj")])

save_name <- file.path(base_dir, paste0(target_gene, "_Clean_NoQC_scTenifoldKnk_results.rds"))
saveRDS(res_knk_clean, file = save_name)

csv_name <- file.path(base_dir, paste0(target_gene, "_Clean_Top20_Affected_Genes.csv"))
write.csv(top20_genes, file = csv_name, row.names = FALSE)

cat("\n🎉 [Finished] 虚拟敲除分析顺利结束！对象与 CSV 已保存。\n")


tp63ko<-readRDS('TP63_Clean_NoQC_scTenifoldKnk_results.rds')
str(tp63ko)



# ==========================================
# 0. 载入必需包
# ==========================================
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggrastr)

cat("⏳ [Volcano] 正在处理数据并优化至 2x2 英寸极小画幅...\n")

diff_reg <- tp63ko$diffRegulation

# ==========================================
# 1. 数据处理与分组排序
# ==========================================
plot_data <- diff_reg %>%
  filter(gene != "TP63") %>%
  mutate(
    p.value = ifelse(p.value == 0, min(p.value[p.value > 0], na.rm = TRUE), p.value),
    logFC = log2(FC),
    negLogP = -log10(p.value),
    Group = case_when(
      gene %in% M1_genes ~ "Module 1",
      gene %in% M3_genes ~ "Module 3",
      TRUE ~ "Other"
    )
  )

plot_data$Group <- factor(plot_data$Group, levels = c("Module 1", "Module 3", "Other"))

# 让 Other 沉在底层，M1/M3 浮在顶层
plot_data <- plot_data %>% arrange(desc(Group))

# 只挑最显著的 6 个防重叠
top_genes <- plot_data %>%
  filter(p.value < 0.05 & Group != "Other") %>%
  arrange(desc(logFC)) %>%
  head(6) %>%
  pull(gene)

plot_data <- plot_data %>%
  mutate(label = ifelse(gene %in% top_genes, gene, NA))

# ==========================================
# 2. 极致排版：绝对 7pt Arial, 适配 2x2 英寸
# ==========================================
custom_theme <- theme(
  text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
  axis.text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
  axis.title = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
  legend.text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
  legend.title = element_blank(),
  legend.position = c(0.25, 0.75),
  legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
  legend.key = element_blank(),
  legend.key.size = unit(0.2, "cm"),
  legend.margin = margin(t = 2, r = 3, b = 2, l = 3),
  panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
  panel.grid.major = element_line(color = "#F0F0F0", linetype = "solid", linewidth = 0.3),
  panel.grid.minor = element_blank()
)

# ==========================================
# 3. 绘制并保存火山图
# ==========================================
cat("🎨 [Volcano] 正在渲染...\n")

p_volcano <- ggplot(plot_data, aes(x = logFC, y = negLogP)) +
  rasterise(
    geom_point(aes(color = Group), size = 0.8, alpha = 0.8, shape = 16),
    dpi = 600
  ) +
  scale_color_manual(
    values = c("Module 1" = "#E3493B", "Module 3" = "#3CB4D2", "Other" = "#E5E5E5")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
             color = "grey40", linewidth = 0.3) +
  geom_text_repel(
    aes(label = label),
    size = 7 / .pt,
    color = "black",
    fontface = "plain",
    family = "Arial",
    max.overlaps = 20,
    box.padding = 0.3,
    point.padding = 0.1,
    segment.color = "grey50",
    segment.alpha = 0.8,
    segment.size = 0.2,
    min.segment.length = 0
  ) +
  labs(x = "log(Fold Change)", y = "-log(p-value)") +
  theme_classic() +
  custom_theme +
  guides(color = guide_legend(override.aes = list(size = 2)))

# ==========================================
# 4. 严格输出 2x2 英寸大小 PDF
# ==========================================
out_pdf <- file.path(base_dir, "Volcano_TP63KO_M1_M3_Rasterized_2x2_7pt.pdf")

cairo_pdf(out_pdf, width = 2.5, height = 2.5, family = "Arial")
print(p_volcano)
dev.off()

cat(paste0("✅ 完美！2x2 英寸极小型火山图已保存至:\n-> ", out_pdf, "\n"))





