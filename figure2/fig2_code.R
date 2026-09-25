# ----------------------------------------------------------------------------
# 步骤 0：全局设置（工作目录 + 所需包，仅在此处统一加载一次）
# ----------------------------------------------------------------------------
setwd('E:/deskup/毕业论文/bechmark/raodong/data/')

library(Seurat)
library(Matrix)                 # readMM
library(tidyverse)              # 含 dplyr / tidyr / ggplot2
library(patchwork)
library(SingleCellExperiment)
library(zellkonverter)          # writeH5AD
library(hdWGCNA)                # CreaterMEscnet / ComputeMEscnetModules / pick_power_scale_free_power
library(UCell)                  # ScoreSignatures_UCell
library(jsonlite)               # fromJSON
library(ggsci)                  # pal_npg

# ============================================================================
# 步骤 1：读取原始数据并构建 Seurat 对象
# ============================================================================
message("正在读取数据...")
counts <- readMM(gzfile("GSE178429_PBMCs_stim_scRNAseq_counts.txt.gz"))
genes  <- read.table(gzfile("GSE178429_PBMCs_stim_scRNAseq_geneNames.txt.gz"),
                     header = FALSE, stringsAsFactors = FALSE)
meta   <- read.table(gzfile("GSE178429_PBMCs_stim_scRNAseq_cellMeta.txt.gz"),
                     header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

# 行为基因，列为细胞
if (ncol(genes) == 1) {
  rownames(counts) <- genes$V1
} else {
  rownames(counts) <- genes$V2   # 第二列为 Gene Symbol
}
colnames(counts) <- rownames(meta)

pbmc <- CreateSeuratObject(counts = counts, meta.data = meta, project = "PBMC_Stim")
print(pbmc)

# ============================================================================
# 步骤 2：细胞质控 (QC)
# ============================================================================
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 阈值可根据上方小提琴图微调
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
message("质控过滤后细胞数：", ncol(pbmc))

# ============================================================================
# 步骤 3：标准化、高变基因、PCA、聚类、UMAP
# ============================================================================
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
ElbowPlot(pbmc, ndims = 30)          # 据此确认 PC 数量

pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.6)
pbmc <- RunUMAP(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

# ============================================================================
# 步骤 4：Marker 点图辅助注释 -> 大类注释 -> 保存 pbmc.RDS
# ============================================================================
# 4.1 整理 Marker 基因（去重、与数据集取交集、按大类分组）
marker_genes_raw <- c('CD3E','CD3D','CD3G','TRAC','TRBC1','TRBC2','TRGC1','TRGC2','TRDC','CD4','CD8A','CD8B','PRF1','GZMB','GZMA','GZMK',
                      'CD4','FOXP3','CTLA4','IL2RA',
                      'NCAM1','NCR1','KLRB1','CX3CR1','CXCR1','CXCR2','FCGR3A','KLRC1','KLRC2','KIT','IL2RA',
                      'CD34','RET','THY1','ITGA6','EPCR','GATA2',
                      'ITGAM','CD14','FCGR3A','S100A8','S100A9','VCAN','ITGAX','IRF4','IRF8','IRF7','XCR1','ID2','ZEB2',
                      'TCF4','LILRA4','NRP1','CLEC4A','CLEC10A','CD2','CD1C','AXL','IL3RA','BATF3',
                      'CSF3R','FCGR3B',
                      'TUBB1', 'PF4V1', 'SELP', 'GP1BA', 'GP9', 'PTGS1',
                      'SDC1','CD38','CD19','JCHAIN','MZB1','MS4A1','CD79A','CD79B','IGHG1','IGHG2','IGHG3','IGHG4','IGHD','IGHM',
                      'MS4A2','HDC','TPSAB1','TPSAB2','TPSD',
                      'MKI67')

marker_genes_unique <- unique(marker_genes_raw)
marker_genes_filt <- intersect(marker_genes_unique, rownames(pbmc))

# 提取三大类的经典代表 marker
T_NK_markers        <- intersect(marker_genes_filt, c('CD3D', 'CD3E', 'CD8A', 'CD4', 'NCAM1', 'NKG7', 'GNLY', 'GZMB', 'PRF1'))
B_Plasma_markers    <- intersect(marker_genes_filt, c('CD19', 'MS4A1', 'CD79A', 'CD79B', 'SDC1', 'CD38', 'JCHAIN', 'MZB1'))
Myeloid_Mac_markers <- intersect(marker_genes_filt, c('CD14', 'FCGR3A', 'S100A8', 'S100A9', 'ITGAM', 'ITGAX', 'CD68', 'C1QA'))
grouped_markers <- c(T_NK_markers, B_Plasma_markers, Myeloid_Mac_markers)

# 4.2 注释前：按原始 Cluster 编号查看 Marker 表达
p1 <- DotPlot(pbmc, features = grouped_markers) +
  RotatedAxis() +
  labs(title = "Marker Expression by Original Clusters") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
print(p1)

# 4.3 手动将 Cluster 编号映射为三大类（依据上方 DotPlot 判断，可自行调整）
Idents(pbmc) <- pbmc$seurat_clusters

tnk_clusters <- c("0","1","2","3","4","5","6","8","10","11","15","18","20","23")
b_clusters   <- c("9","13","19","22","24")
mac_clusters <- c("7","12","14","16","17","21")

new.cluster.ids <- rep("Other", length(levels(Idents(pbmc))))
names(new.cluster.ids) <- levels(Idents(pbmc))
new.cluster.ids[tnk_clusters] <- "T_NK_Cell"
new.cluster.ids[b_clusters]   <- "B_Cell"
new.cluster.ids[mac_clusters] <- "Macrophage"

pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$Broad_CellType <- Idents(pbmc)

# 4.4 注释后：按大类查看 Marker 表达
p2 <- DotPlot(pbmc, features = grouped_markers) +
  RotatedAxis() +
  labs(title = "Marker Expression by Broad Cell Types") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
print(p2)

saveRDS(pbmc, file = 'pbmc.RDS')   # 检查点 1

# ============================================================================
# 步骤 5：导出 h5ad（供外部 expiMap / Python 使用）
# ============================================================================
sce <- as.SingleCellExperiment(pbmc, assay = "RNA")

if ("pca" %in% names(pbmc@reductions)) {
  reducedDims(sce)$PCA <- Embeddings(pbmc, "pca")
}
if ("umap" %in% names(pbmc@reductions)) {
  reducedDims(sce)$UMAP <- Embeddings(pbmc, "umap")
}

writeH5AD(sce, file = "pbmc.h5ad", X_name = "counts", compression = "gzip")

# ============================================================================
# 步骤 6：外部 Python（在此暂停，手动执行以下两个算法）
#   a) expiMap：读取 pbmc.h5ad 运行 expiMap，导出"基因-模块"表
#      （表格需含 Gene 与 Module 两列，供步骤 9 使用）
#   b) SPECTRA：运行 SPECTRA，导出 Spectra_Module_Top50_Genes.json
# 步骤 9 会自动读取这两个文件；若实际文件名不同，请修改步骤 9 中的路径。
# ============================================================================

# ============================================================================
# 步骤 7：MEscnet —— Metacell 构建与软阈值
# ============================================================================
# 断点续跑：若单独运行本段（不在同一会话中），请取消下面一行注释
# pbmc <- readRDS("pbmc.RDS")
seurat_obj <- pbmc

# 7.1 用 markers.csv 字典基因 ∩ 高变基因，作为 MEscnet 的输入基因
genesets <- read.csv('markers.csv')
genes <- VariableFeatures(seurat_obj)   # 高变基因已在步骤 3 计算

gene_lists <- strsplit(gsub("\\[|\\]|'", "", genesets$gene_set), ", ")
all_dict_genes <- unique(unlist(gene_lists))
common_genes <- intersect(all_dict_genes, genes)
message("MEscnet 字典基因与高变基因的交集数量: ", length(common_genes))

# 7.2 全部细胞一起构建 Metacell
seurat_obj@meta.data$All_Cells <- "All_Immune"
total_cells <- ncol(seurat_obj)
target_metacells_all <- max(1500, round(total_cells / 3))

seurat_obj <- CreaterMEscnet(
  seurat_obj,
  nfeatures = 1000,
  group.by = "All_Cells",
  reduction = "umap",
  k = 10,
  max_shared = 10,
  target_metacells = target_metacells_all,
  group_name = "All_Immune",
  gene_list = common_genes
)

# 7.3 计算 Scale Free 拓扑软阈值
seurat_obj <- pick_power_scale_free_power(seurat_obj)
message("全局网络软阈值计算完毕！")

# ============================================================================
# 步骤 8：MEscnet —— 参数扫描（resolution x number），保存全部模块
# ============================================================================
# 注意：下游出图只使用 resolution_4_number_4；
#       若为节省时间，可把 resolution_list 精简为 c(4)。
resolution_list <- seq(0, 5, by = 0.1)
number_list <- c(2, 3, 4)

all_MEscnet_modules <- list()

for (res in resolution_list) {
  for (num in number_list) {
    cat("Running: resolution =", res, "number =", num, "\n")

    seurat_tmp <- ComputeMEscnetModules(
      seurat_obj,
      min_genes_per_module = 0,
      resolution = res,
      number = num
    )

    module_genes <- seurat_tmp@misc$MEscnet$MEscnet_modules$gene_lists
    module_ids   <- seurat_tmp@misc$MEscnet$MEscnet_modules$module_ids
    names(module_genes) <- module_ids

    result_name <- paste0("resolution_", res, "_number_", num)
    all_MEscnet_modules[[result_name]] <- module_genes

    rm(seurat_tmp)
    gc()
  }
}

save(all_MEscnet_modules, file = 'mescnet_modules.Rdata')   # 检查点 2

# ============================================================================
# 步骤 9：三大算法 UCell 打分，保存全部得分
# ============================================================================
# 断点续跑：若单独运行本段，请取消下面两行注释
# seurat_obj <- readRDS("pbmc.RDS")
# load("mescnet_modules.Rdata")

cat("⏳ 正在准备表达矩阵与基因列表...\n")
expr_matrix <- GetAssayData(seurat_obj, layer = "data")

# 9.1 expiMap 基因-模块表（外部 Python 输出；如实际文件名不同请修改）
if (!exists("expimap")) {
  expimap <- read.csv("expiMap_module_genes.csv")
}
expimap_list <- split(expimap$Gene, expimap$Module)

# 9.2 SPECTRA 基因列表（外部 Python 输出的 JSON）
spectra_list <- fromJSON("Spectra_Module_Top50_Genes.json")

# 9.3 expiMap 打分
cat("🚀 正在计算 expiMap modules 的 UCell 分数 (ncores = 64)...\n")
expimap_matrix <- ScoreSignatures_UCell(
  matrix = expr_matrix,
  features = expimap_list,
  ncores = 64
)
expimap_scores <- as.list(as.data.frame(expimap_matrix))

# 9.4 SPECTRA 打分
cat("🚀 正在计算 Spectra modules 的 UCell 分数 (ncores = 64)...\n")
spectra_matrix <- ScoreSignatures_UCell(
  matrix = expr_matrix,
  features = spectra_list,
  ncores = 64
)
spectra_scores <- as.list(as.data.frame(spectra_matrix))

# 9.5 MEscnet 打分（遍历全部参数组合；ncores 可按机器核数调整）
cat("🚀 正在计算 MEscnet modules 的 UCell 分数（大列表，请耐心等待）...\n")
mescnet_scores <- list()
for (param in names(all_MEscnet_modules)) {
  cat("   -> 正在处理 MEscnet 参数:", param, "\n")
  module_list <- all_MEscnet_modules[[param]]
  names(module_list) <- paste0("MEscnet_Mod_", seq_along(module_list))

  ucell_matrix <- ScoreSignatures_UCell(
    matrix = expr_matrix,
    features = module_list,
    ncores = 64
  )
  mescnet_scores[[param]] <- as.list(as.data.frame(ucell_matrix))
}

save(expimap_scores, spectra_scores, mescnet_scores,
     file = "Benchmark1_UCell_Scores.Rdata")   # 检查点 3
cat("🎉 全部打分计算完成！\n")

# ============================================================================
# 步骤 10：UMAP 基础图（细胞类型 + 实验条件）
# ============================================================================
# 断点续跑：若单独运行本段，请取消下面一行注释
# seurat_obj <- readRDS("pbmc.RDS")

embed_coords <- Embeddings(seurat_obj, "umap")
df <- data.frame(
  dim_1 = embed_coords[, 1],
  dim_2 = embed_coords[, 2],
  cell_type = seurat_obj$Broad_CellType,
  condition = seurat_obj$Condition
)
df$cell_type <- factor(df$cell_type, levels = c("B_Cell", "T_NK_Cell", "Macrophage"))

# 配色（Nature 风格）
color_celltype <- c(
  "B_Cell"     = "#1e4899ff",
  "T_NK_Cell"  = "#9a3f1fff",
  "Macrophage" = "#1f9a86ff"
)
color_condition <- colorRampPalette(ggsci::pal_npg()(10))(length(unique(df$condition)))

theme_umap <- theme_void() + theme(
  text = element_text(size = 14, color = "black", face = "plain"),
  legend.position = "right",
  legend.title = element_blank(),
  legend.text = element_text(size = 14, color = "black", face = "plain"),
  plot.margin = margin(10, 10, 10, 10)
)

p_celltype <- ggplot(df, aes(x = dim_1, y = dim_2)) +
  geom_point(aes(color = cell_type), size = 0.2, alpha = 0.8) +
  scale_color_manual(values = color_celltype) +
  coord_fixed() +
  theme_umap +
  guides(color = guide_legend(override.aes = list(size = 4)))

p_condition <- ggplot(df, aes(x = dim_1, y = dim_2)) +
  geom_point(aes(color = condition), size = 0.2, alpha = 0.8) +
  scale_color_manual(values = color_condition) +
  coord_fixed() +
  theme_umap +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave("UMAP_CellType_Circled_Nature.pdf", plot = p_celltype, width = 7, height = 4)
ggsave("UMAP_CellType_Circled_Nature.png", plot = p_celltype, width = 7, height = 4, dpi = 600)
ggsave("UMAP_Condition_Uncircled_Nature.pdf", plot = p_condition, width = 7, height = 4)
ggsave("UMAP_Condition_Uncircled_Nature.png", plot = p_condition, width = 7, height = 4, dpi = 600)

# ============================================================================
# 步骤 11：最终面板图 —— 3 算法 (MEscnet/SPECTRA/expiMap) x 3 刺激 (LPS/PMA/IFN)
# ============================================================================
# 断点续跑：若单独运行本段，请取消下面两行注释
# pbmc <- readRDS("pbmc.RDS")
# load("Benchmark1_UCell_Scores.Rdata")

# 11.1 实验设计参数
target_celltypes <- c("B_Cell", "T_NK_Cell", "Macrophage")
my_colors <- c("Unperturbed" = "#22547F", "Perturbed" = "#C05527")

# 11.2 提取 MEscnet 特定参数（下游仅用该组参数）
mescnet_score <- mescnet_scores[["resolution_4_number_4"]]

# 11.3 Min-Max 标准化 (0-1)
min_max_scale <- function(df) {
  df <- as.data.frame(df)
  df_scaled <- as.data.frame(lapply(df, function(x) {
    min_val <- min(x, na.rm = TRUE)
    max_val <- max(x, na.rm = TRUE)
    if (max_val == min_val) return(rep(0, length(x)))
    return((x - min_val) / (max_val - min_val))
  }))
  rownames(df_scaled) <- rownames(df)
  return(df_scaled)
}

message("正在进行 Min-Max 标准化...")
mescnet_scaled <- min_max_scale(mescnet_score)
spectra_scaled <- min_max_scale(spectra_scores)
expimap_scaled <- min_max_scale(expimap_scores)

rownames(mescnet_scaled) <- colnames(pbmc)
rownames(spectra_scaled) <- colnames(pbmc)
rownames(expimap_scaled) <- colnames(pbmc)

# 11.4 Pseudobulk 聚合（按 Donor x Condition x Broad_CellType 取均值）
prepare_pb_input <- function(score_mat, pbmc_obj) {
  score_df <- as.data.frame(score_mat)
  score_df <- score_df[colnames(pbmc_obj), , drop = FALSE]
  score_df$Donor <- pbmc_obj$Donor
  score_df$Condition <- pbmc_obj$Condition
  score_df$Broad_CellType <- as.character(pbmc_obj$Broad_CellType)
  return(score_df)
}

do_aggregate <- function(df) {
  mod_cols <- setdiff(colnames(df), c("Donor", "Condition", "Broad_CellType"))
  df %>%
    group_by(Donor, Condition, Broad_CellType) %>%
    summarise(across(all_of(mod_cols), mean, na.rm = TRUE), .groups = "drop")
}

message("正在进行 Pseudobulk 聚合...")
pb_mescnet <- do_aggregate(prepare_pb_input(mescnet_scaled, pbmc))
pb_spectra <- do_aggregate(prepare_pb_input(spectra_scaled, pbmc))
pb_expimap <- do_aggregate(prepare_pb_input(expimap_scaled, pbmc))
cat("✅ 数据预处理完成！\n")

# 11.5 核心参数映射表（仅限 GolgiPlug 条件）
best_modules <- list(
  MEscnet = list(LPS = "MEscnet_Mod_38_UCell", PMA = "MEscnet_Mod_26_UCell", IFN = "MEscnet_Mod_24_UCell"),
  SPECTRA = list(LPS = "Spectra_Factor_172_UCell", PMA = "Spectra_Factor_183_UCell", IFN = "Spectra_Factor_155_UCell"),
  Expimap = list(LPS = "all_TLR_signaling_UCell", PMA = "T_tcr.activation_UCell", IFN = "all_type.I.ifn.response_UCell")
)

golgi_config <- list(
  LPS = c("ControlGolgiPlug_6h", "LPSGolgiPlug_6h", "LPS perturbation"),
  IFN = c("ControlGolgiPlug_6h", "IFNGolgiPlug_6h", "IFNγ perturbation"),
  PMA = c("ControlGolgiPlug_6h", "PMAGolgiPlug_6h", "TCR perturbation")
)

# 11.6 单算法单行 3 面板图
plot_algorithm_row <- function(algo_name, pb_data, mods) {

  extract_data <- function(stim_type, mod_name) {
    ctrl_cond  <- golgi_config[[stim_type]][1]
    stim_cond  <- golgi_config[[stim_type]][2]
    main_title <- golgi_config[[stim_type]][3]

    df <- pb_data %>% filter(Condition %in% c(ctrl_cond, stim_cond))
    df$Status <- ifelse(df$Condition == ctrl_cond, "Unperturbed", "Perturbed")

    factor_name <- switch(stim_type,
                          "LPS" = "LPS response factor",
                          "PMA" = "TCR activation factor",
                          "IFN" = "IFNγ response factor")

    df$Facet_Group <- paste0(main_title, "\n", factor_name, "\n(", mod_name, ")")
    df$Plot_Score  <- df[[mod_name]]
    return(df)
  }

  lps_data <- extract_data("LPS", mods$LPS)
  pma_data <- extract_data("PMA", mods$PMA)
  ifn_data <- extract_data("IFN", mods$IFN)

  plot_df_final <- bind_rows(lps_data, pma_data, ifn_data)
  plot_df_final$Broad_CellType <- factor(plot_df_final$Broad_CellType, levels = target_celltypes)
  plot_df_final$Status <- factor(plot_df_final$Status, levels = c("Unperturbed", "Perturbed"))
  plot_df_final$Facet_Group <- factor(plot_df_final$Facet_Group, levels = unique(plot_df_final$Facet_Group))

  p <- ggplot(plot_df_final, aes(x = Broad_CellType, y = Plot_Score, fill = Status)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9,
                 position = position_dodge(0.8), color = "black", linewidth = 0.5) +
    geom_point(aes(group = Status),
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               size = 1.2, color = "#000000", alpha = 0.7, show.legend = FALSE) +
    scale_fill_manual(values = my_colors) +
    facet_wrap(~ Facet_Group, ncol = 3, scales = "free_y") +
    theme_classic() +
    theme(
      text = element_text(size = 14, color = "black", face = "plain"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, color = "black", face = "plain"),
      axis.text.y = element_text(size = 14, color = "black", face = "plain"),
      axis.title.y = element_text(size = 14, color = "black", face = "plain"),
      axis.title.x = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = 14, color = "black", face = "plain", lineheight = 1.2),
      legend.position = "none"   # 单行隐藏图例，最后统一收集到底部
    ) +
    labs(y = paste0(algo_name, "\nModule score"))

  return(p)
}

# 11.7 生成 3 行并拼接为最终 PDF
message("正在组装终极面板...")
p_mescnet <- plot_algorithm_row("MEscnet", pb_mescnet, best_modules$MEscnet)
p_spectra <- plot_algorithm_row("SPECTRA", pb_spectra, best_modules$SPECTRA)
p_expimap <- plot_algorithm_row("expiMap", pb_expimap, best_modules$Expimap)

final_plot <- p_mescnet / p_spectra / p_expimap +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black", face = "plain"),
    legend.key.size = unit(0.6, "cm")
  )

pdf_filename <- "Benchmark_GolgiPlug_All_Algorithms_Panel.pdf"
pdf(pdf_filename, width = 15, height = 11)
print(final_plot)
dev.off()

cat("🎉 3x3 终极九宫格图已保存为：", pdf_filename, "\n")











# 1. 加载必要包 -------------------------------------------------
if (!require("jsonlite")) install.packages("jsonlite")
library(jsonlite)

# 2. 检查数据文件是否存在（请确保文件在当前工作目录下）---------
required_files <- c(
  "mescnet_modules.Rdata",
  "Spectra_Module_Top50_Genes.json",
  "EXPIMAP_active_module_top50_genes_UCell.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("❌ 以下文件不存在，请检查路径：", paste(missing_files, collapse = ", "))
}

# 3. 提取 MEscnet 基因 -----------------------------------------
load("mescnet_modules.Rdata")
mescnet_mods <- all_MEscnet_modules[["resolution_4_number_4"]]

mescnet_lps <- mescnet_mods[["38"]]
mescnet_pma <- mescnet_mods[["26"]]
mescnet_ifn <- mescnet_mods[["24"]]

# 4. 提取 SPECTRA 基因 -----------------------------------------
spectra_list <- fromJSON("Spectra_Module_Top50_Genes.json")

spectra_lps <- spectra_list[["Spectra_Factor_172"]]
spectra_pma <- spectra_list[["Spectra_Factor_183"]]
spectra_ifn <- spectra_list[["Spectra_Factor_155"]]

# 5. 提取 expiMap 基因 -----------------------------------------
# 注意：header = FALSE 表示 CSV 无表头，V1=模块名，V2=基因名
expimap_df <- read.csv("EXPIMAP_active_module_top50_genes_UCell.csv", 
                       header = FALSE, 
                       stringsAsFactors = FALSE)

# 使用更新后的正确模块名字进行匹配（已修正为 "-" 连接符）
expimap_lps <- expimap_df$V2[expimap_df$V1 == "all_TLR_signaling"]
expimap_pma <- expimap_df$V2[expimap_df$V1 == "T_tcr-activation"]
expimap_ifn <- expimap_df$V2[expimap_df$V1 == "all_type-I-ifn-response"]

# 6. 检查是否成功匹配到数据 ------------------------------------
if (length(expimap_lps) == 0) warning("⚠️ 未找到 expiMap LPS 模块 'all_TLR_signaling'")
if (length(expimap_pma) == 0) warning("⚠️ 未找到 expiMap PMA 模块 'T_tcr-activation'")
if (length(expimap_ifn) == 0) warning("⚠️ 未找到 expiMap IFN 模块 'all_type-I-ifn-response'")

# 7. 组装最终列表 ---------------------------------------------
final_list <- list(
  "MEscnet_LPS (Mod_38)"        = mescnet_lps,
  "MEscnet_TCR (Mod_26)"        = mescnet_pma,
  "MEscnet_IFN (Mod_24)"        = mescnet_ifn,
  "SPECTRA_LPS (Factor_172)"    = spectra_lps,
  "SPECTRA_TCR (Factor_183)"    = spectra_pma,
  "SPECTRA_IFN (Factor_155)"    = spectra_ifn,
  "expiMap_LPS (TLR_signaling)" = expimap_lps,
  "expiMap_TCR (T_activation)"  = expimap_pma,
  "expiMap_IFN (IFN_response)"  = expimap_ifn
)

# 8. 对齐长度并转换为数据框 ------------------------------------
max_len <- max(sapply(final_list, length))

padded_list <- lapply(final_list, function(x) {
  length(x) <- max_len
  return(x)
})

df_out <- as.data.frame(do.call(rbind, padded_list), stringsAsFactors = FALSE)
df_out <- cbind(Module_Name = rownames(df_out), df_out)
rownames(df_out) <- NULL

# 9. 导出为 CSV ------------------------------------------------
write.table(df_out, 
            file = "Target_9_Modules_Genes.csv", 
            sep = ",", 
            row.names = FALSE, 
            col.names = FALSE, 
            na = "")

cat("✅ 大功告成！CSV 表格已生成：Target_9_Modules_Genes.csv\n")
cat("共导出", nrow(df_out), "个模块，最长模块含", max_len, "个基因。\n")


# ============================================================================
# 合并脚本：EGAD + PPI Benchmark（10 个 hotspot 数据集 x 8 种网络）
# 由 egad_ppi_processing.R + plot.R + sup.R 整理合并
# （删除：重复加载的包、processing 中被 plot.R 取代的旧清洗步骤、
#        sup.R 中与 EGAD/PPI 无关的 UMAP 导出代码、重复的小提琴图循环、
#        被最终版本覆盖的旧图例代码等）
#
# 流程：
#   步骤 1  纯计算函数 (EGAD AUC / PPI Z-score)
#   步骤 2  遍历所有数据集计算 -> benchmark.Rdata（检查点 1）
#   步骤 3  清洗为最终 8 种网络 + 数据集缩写
#   步骤 4  聚合 (median) + 全局排名 + factor 排序
#   步骤 5  主图: EGAD_AUC_rank_final.pdf / PPI_Zscore_rank_final.pdf
#           + 独立图例 Legend_Standalone_EGAD.pdf / PPI.pdf
#           + 保存 benchmark_8networks_cor_matrix_single.Rdata（检查点 2）
#   步骤 6  补充图: 各数据集 EGAD 小提琴图 (EGAD_Violin_<数据集>.pdf)
#
# 断点续跑：每个大步骤开头均保留了“若单独运行本段”的读取检查点。
# ============================================================================

# ----------------------------------------------------------------------------
# 步骤 0：全局设置（工作目录 + 所需包）
# ----------------------------------------------------------------------------
setwd('E:/deskup/毕业论文/bechmark/edga_ppi/data/')

library(tidyverse)      # dplyr / tidyr / ggplot2
library(patchwork)
library(cowplot)        # get_legend / ggdraw
library(grid)           # unit()
library(igraph)
library(org.Hs.eg.db)   # 基因 ID 转换
library(EGAD)           # make_annotations / run_GBA

# ----------------------------------------------------------------------------
# 步骤 1：纯计算函数（不含画图逻辑）
# ----------------------------------------------------------------------------

# 1.1 计算 EGAD AUC（对每个网络运行 GBA，输出 nv_auc / nd_auc）
calc_egad_auc <- function(TOMs, go_data) {
  plot_df <- data.frame()
  goterms <- unique(go_data[, 3])
  
  for (net_name in names(TOMs)) {
    tom <- TOMs[[net_name]]
    gene1 <- rownames(tom)
    
    # 动态识别基因 ID 类型并转换，防止硬报错中断循环
    valid_symbols <- intersect(gene1, keys(org.Hs.eg.db, keytype = "SYMBOL"))
    valid_ensembl <- intersect(gsub("\\..*", "", gene1), keys(org.Hs.eg.db, keytype = "ENSEMBL"))
    valid_entrez  <- intersect(gene1, keys(org.Hs.eg.db, keytype = "ENTREZID"))
    
    if (length(valid_symbols) > 0) {
      gene_symbols <- mapIds(org.Hs.eg.db, keys = gene1, column = "ENTREZID",
                             keytype = "SYMBOL", multiVals = "first")
    } else if (length(valid_ensembl) > 0) {
      gene1_clean <- gsub("\\..*", "", gene1)
      gene_symbols <- mapIds(org.Hs.eg.db, keys = gene1_clean, column = "ENTREZID",
                             keytype = "ENSEMBL", multiVals = "first")
      names(gene_symbols) <- gene1
    } else if (length(valid_entrez) > 0) {
      gene_symbols <- setNames(gene1, gene1)
    } else {
      cat("⚠️ Warning: 网络 [", net_name, "] 的基因名格式无法识别，跳过...\n")
      next
    }
    
    result_df <- data.frame(EntrezID = gene_symbols[gene1], Symbol = gene1)
    gene_table <- na.omit(result_df)
    
    if (nrow(gene_table) < 10) {
      cat("⚠️ Warning: 网络 [", net_name, "] 成功转换的基因过少，跳过...\n")
      next
    }
    
    # 补全下三角
    tom[lower.tri(tom)] <- t(tom)[lower.tri(tom)]
    
    # 筛选对齐基因并转为 Entrez ID
    genelist <- rownames(tom)
    genelist <- genelist[genelist %in% gene_table$Symbol]
    tom <- tom[genelist, genelist]
    
    ix <- match(genelist, gene_table$Symbol)
    gene_ids <- gene_table$EntrezID[ix]
    colnames(tom) <- gene_ids
    rownames(tom) <- gene_ids
    
    annotations <- make_annotations(go_data[, c(2, 3)], gene_ids, goterms)
    
    tryCatch({
      GO_groups_voted <- run_GBA(tom, annotations, max = Inf)
      df <- data.frame(
        nv_auc = as.numeric(GO_groups_voted[[1]][, 1]),  # Neighbor Voting AUROC
        nd_auc = as.numeric(GO_groups_voted[[1]][, 3]),
        network = net_name
      )
      plot_df <- rbind(plot_df, df)
    }, error = function(e) {
      cat("⚠️ Warning: 网络 [", net_name, "] 运行 run_GBA 时出错:", e$message, "\n")
    })
  }
  return(plot_df)
}

# 1.2 计算 PPI Z-score（Top-N 共表达边 与 PPI 边的富集，度保留随机化）
calc_ppi_zscore <- function(TOMs, ppi_edges, top_n = 3000, bootstraps = 100) {
  # 以第一个网络为基准，过滤 PPI 到共有基因
  ref_tom <- TOMs[[1]]
  tom_genes <- rownames(ref_tom)
  
  ppi_filtered <- ppi_edges %>%
    filter(gene1 %in% tom_genes, gene2 %in% tom_genes)
  
  G_ppi <- graph_from_data_frame(ppi_filtered, directed = FALSE)
  edges_ppi <- igraph::as_data_frame(G_ppi, what = "edges") %>%
    mutate(edge_id = paste(pmin(from, to), pmax(from, to), sep = "_"))
  
  results <- data.frame()
  
  for (net in names(TOMs)) {
    tom <- TOMs[[net]]
    tom[lower.tri(tom, diag = TRUE)] <- NA
    
    # 提取 Top N 的共表达边
    edges_coexp <- as.data.frame(as.table(tom)) %>%
      filter(!is.na(Freq)) %>%
      rename(gene1 = Var1, gene2 = Var2, weight = Freq) %>%
      arrange(desc(weight)) %>%
      head(top_n)
    
    G_coexp <- graph_from_data_frame(edges_coexp, directed = FALSE)
    edges_coexp_df <- igraph::as_data_frame(G_coexp, what = "edges") %>%
      mutate(edge_id = paste(pmin(from, to), pmax(from, to), sep = "_"))
    
    obs <- length(intersect(edges_coexp_df$edge_id, edges_ppi$edge_id))
    
    # 度保留随机化
    random <- replicate(bootstraps, {
      G_random <- rewire(G_ppi, keeping_degseq(niter = ecount(G_ppi) * 10))
      edges_random <- igraph::as_data_frame(G_random, what = "edges") %>%
        mutate(edge_id = paste(pmin(from, to), pmax(from, to), sep = "_"))
      length(intersect(edges_coexp_df$edge_id, edges_random$edge_id))
    })
    
    rnd_mean <- mean(random, na.rm = TRUE)
    rnd_sd <- sd(random, na.rm = TRUE)
    z_score <- (obs - rnd_mean) / rnd_sd
    
    results <- rbind(results, data.frame(
      network = net,
      cutoff = top_n,
      obs = obs,
      rnd_mean = rnd_mean,
      rnd_sd = rnd_sd,
      z_score = z_score
    ))
  }
  return(results)
}

# ----------------------------------------------------------------------------
# 步骤 2：加载全局 GO/PPI 数据，遍历所有数据集计算 EGAD + PPI
# ----------------------------------------------------------------------------
cat("⏳ 正在加载全局 GO 和 PPI 数据...\n")
setwd('E:/deskup/Mescnet/data/')
load('go.Rdata')   # 提供 go 数据框

ppi <- read.csv('data2.csv')
ppi_filtered <- ppi %>%
  filter(combined_score > 500) %>%
  select(gene1, gene2) %>%
  distinct()

path <- 'E:/deskup/hotspot/'
files <- list.files(path)

all_benchmark_results <- list()

cat("⏳ 开始遍历计算所有数据集...\n")
for (i in seq_along(files)) {
  dataset_name <- files[i]
  cat("Processing dataset:", dataset_name, "...\n")
  
  path_2 <- file.path(path, dataset_name)
  if (!dir.exists(path_2)) {
    cat("⚠️ 路径不存在，跳过:", path_2, "\n")
    next
  }
  setwd(path_2)
  
  if (!file.exists('TOMs.Rdata')) {
    cat("⚠️ 未找到 TOMs.Rdata，跳过...\n")
    next
  }
  load('TOMs.Rdata')
  
  # 追加额外的 hotspot 网络矩阵
  if (file.exists('output1.csv') & file.exists('output2.csv')) {
    hotspot1 <- read.csv('output1.csv')
    hotspot2 <- read.csv('output2.csv')
    rownames(hotspot1) <- colnames(hotspot1)
    rownames(hotspot2) <- colnames(hotspot2)
    TOMs[['hotspot1']] <- as.matrix(hotspot1)
    TOMs[['hotspot2']] <- as.matrix(hotspot2)
  }
  
  # 1. EGAD 结果
  egad_res <- calc_egad_auc(TOMs, go)
  # 2. PPI 结果
  ppi_res <- calc_ppi_zscore(TOMs, ppi_filtered, top_n = 3000, bootstraps = 100)
  
  all_benchmark_results[[dataset_name]] <- list(egad = egad_res, ppi = ppi_res)
}
cat("✅ 所有数据计算完毕！\n")

# 合并所有数据集
combined_egad <- bind_rows(lapply(names(all_benchmark_results), function(ds) {
  if (nrow(all_benchmark_results[[ds]]$egad) > 0) {
    all_benchmark_results[[ds]]$egad %>% mutate(Dataset = ds)
  }
}))

combined_ppi <- bind_rows(lapply(names(all_benchmark_results), function(ds) {
  if (nrow(all_benchmark_results[[ds]]$ppi) > 0) {
    all_benchmark_results[[ds]]$ppi %>% mutate(Dataset = ds)
  }
}))

setwd('E:/deskup/毕业论文/bechmark/edga_ppi/data/')
save(combined_egad, combined_ppi, all_benchmark_results, file = 'benchmark.Rdata')  # 检查点 1

# ----------------------------------------------------------------------------
# 步骤 3：清洗（最终 8 种网络 + 数据集缩写）
# ----------------------------------------------------------------------------
# 断点续跑：若单独运行本段，请取消下面一行注释
# load('benchmark.Rdata')

dataset_mapping <- c(
  "scp1039" = "BRCA", "scp1265" = "TAA", "scp1289" = "COVID19", "scp1303" = "CM",
  "scp1526" = "ISL", "scp1671" = "HBM", "scp1731" = "HBA", "scp1849" = "ICM",
  "scp1852" = "CHD", "scp1903" = "ADI"
)

network_mapping <- c(
  "cor_matrix_single" = "pearson correlation",
  "pcor"              = "partial correlation",
  "bayes_cor1"        = "bayes_cor1",
  "bayes_cor2"        = "bayes_cor2",
  "bayes_cor3"        = "bayes_cor3",
  "hotspot1"          = "hotspot",
  "hdwgcna"           = "hdWGCNA",
  "random"            = "random network"
)

keep_networks <- c("pearson correlation", "partial correlation", "bayes_cor1",
                   "bayes_cor2", "bayes_cor3", "hotspot", "hdWGCNA", "random network")

egad_clean <- combined_egad %>%
  filter(network %in% names(network_mapping)) %>%
  mutate(
    network = recode(network, !!!network_mapping),
    Dataset = recode(Dataset, !!!dataset_mapping)
  ) %>%
  filter(network %in% keep_networks)

ppi_clean <- combined_ppi %>%
  filter(network %in% names(network_mapping)) %>%
  mutate(
    network = recode(network, !!!network_mapping),
    Dataset = recode(Dataset, !!!dataset_mapping)
  ) %>%
  filter(network %in% keep_networks)

# 校验：最终必须恰好是 8 种网络
if (length(unique(egad_clean$network)) != 8) warning("⚠️ EGAD 最终不是 8 种网络，请检查！")
if (length(unique(ppi_clean$network))  != 8) warning("⚠️ PPI 最终不是 8 种网络，请检查！")

# ----------------------------------------------------------------------------
# 步骤 4：聚合 (median) + 全局排名 + factor 顺序
# ----------------------------------------------------------------------------
# 4.1 EGAD / PPI 按 网络 x 数据集 聚合成中位数
heat_egad <- egad_clean %>%
  group_by(network, Dataset) %>%
  summarise(median_auc = median(nv_auc, na.rm = TRUE), .groups = "drop")

heat_ppi <- ppi_clean %>%
  group_by(network, Dataset) %>%
  summarise(median_z = median(z_score, na.rm = TRUE), .groups = "drop")

# 4.2 全局排名（AUC / Z-score 降序）
rank_egad <- heat_egad %>%
  group_by(network) %>%
  summarise(overall_auc = median(median_auc, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(overall_auc)) %>%
  mutate(Rank = row_number())

rank_ppi <- heat_ppi %>%
  group_by(network) %>%
  summarise(overall_z = median(median_z, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(overall_z)) %>%
  mutate(Rank = row_number())

# 4.3 网络纵轴顺序：Rank 1 放最上面
egad_network_order <- rank_egad$network
ppi_network_order  <- rank_ppi$network

# 4.4 数据集横轴顺序：按 PPI Rank 1 算法的 Z-score 升序
top1_ppi_algo <- ppi_network_order[1]

dataset_order <- heat_ppi %>%
  filter(network == top1_ppi_algo) %>%
  arrange(median_z) %>%
  pull(Dataset) %>%
  as.character()

# 防止某些数据集缺失
all_datasets <- unique(c(as.character(egad_clean$Dataset),
                         as.character(ppi_clean$Dataset)))
dataset_order <- c(dataset_order, setdiff(all_datasets, dataset_order))

# 4.5 写入 factor 顺序
egad_clean$network <- factor(egad_clean$network, levels = rev(egad_network_order))
heat_egad$network <- factor(heat_egad$network, levels = rev(egad_network_order))
rank_egad$network <- factor(rank_egad$network, levels = rev(egad_network_order))
egad_clean$Dataset <- factor(egad_clean$Dataset, levels = dataset_order)
heat_egad$Dataset <- factor(heat_egad$Dataset, levels = dataset_order)

ppi_clean$network <- factor(ppi_clean$network, levels = rev(ppi_network_order))
heat_ppi$network <- factor(heat_ppi$network, levels = rev(ppi_network_order))
rank_ppi$network <- factor(rank_ppi$network, levels = rev(ppi_network_order))
ppi_clean$Dataset <- factor(ppi_clean$Dataset, levels = dataset_order)
heat_ppi$Dataset <- factor(heat_ppi$Dataset, levels = dataset_order)

# ----------------------------------------------------------------------------
# 步骤 5：主图（heatmap + 右侧 boxplot/barplot + 独立图例）
# ----------------------------------------------------------------------------
# 5.1 主题（14号 黑色 不加粗）
nature_theme <- theme_classic() +
  theme(
    text = element_text(size = 14, color = "black", face = "plain"),
    axis.text.x = element_text(size = 14, color = "black", face = "plain",
                               angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, color = "black", face = "plain"),
    axis.title = element_text(size = 14, color = "black", face = "plain"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.text = element_text(size = 14, color = "black", face = "plain"),
    legend.title = element_text(size = 14, color = "black", face = "plain")
  )

heatmap_theme <- nature_theme +
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top"
  )

# 5.2 EGAD heatmap（图例：3 个刻度 + 横向 colorbar）
p_heat_egad <- ggplot(heat_egad, aes(x = Dataset, y = network, fill = median_auc)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(
    colors = c("#0B1D4A", "#265A8C", "#5895A9", "#A2BCA4", "#D5B08A", "#FCE2D4"),
    name = "Median\nAUC",
    breaks = c(0.50, 0.60, 0.65),
    labels = c("0.50", "0.60", "0.65"),
    guide = guide_colorbar(direction = "horizontal", title.position = "left",
                           title.hjust = 0.5, label.position = "bottom",
                           barwidth = unit(6, "cm"), barheight = unit(0.6, "cm"),
                           ticks = TRUE)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed() +
  heatmap_theme

# 5.3 EGAD 右侧 Boxplot + Rank 标签
label_x_egad <- max(egad_clean$nv_auc, na.rm = TRUE) + 0.05
x_min_egad <- min(egad_clean$nv_auc, na.rm = TRUE)

p_box_egad <- ggplot(egad_clean, aes(x = nv_auc, y = network)) +
  geom_boxplot(fill = "#C9B28C", color = "black", linewidth = 0.4,
               width = 0.65, outlier.shape = NA) +
  geom_text(data = rank_egad,
            aes(x = label_x_egad, y = network, label = paste0("#", Rank)),
            inherit.aes = FALSE, size = 5, color = "black") +
  scale_x_continuous(name = "AUC",
                     limits = c(x_min_egad - 0.02, label_x_egad + 0.05),
                     expand = c(0, 0)) +
  nature_theme +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dotted", linewidth = 0.4),
    plot.margin = margin(5, 25, 5, 0)
  )

# 5.4 PPI heatmap（图例：3 个刻度 + 横向 colorbar）
ppi_min <- min(heat_ppi$median_z, na.rm = TRUE)
ppi_max <- max(heat_ppi$median_z, na.rm = TRUE)
ppi_breaks <- c(ppi_min, (ppi_min + ppi_max) / 2, ppi_max)

p_heat_ppi <- ggplot(heat_ppi, aes(x = Dataset, y = network, fill = median_z)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "white", high = "#542788", name = "Median\nZ-score",
    breaks = ppi_breaks,
    labels = round(ppi_breaks, 1),
    guide = guide_colorbar(direction = "horizontal", title.position = "left",
                           title.hjust = 0.5, label.position = "bottom",
                           barwidth = unit(6, "cm"), barheight = unit(0.6, "cm"),
                           ticks = TRUE)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed() +
  heatmap_theme

# 5.5 PPI 右侧 Rank Barplot
min_z_val <- min(rank_ppi$overall_z, na.rm = TRUE)
max_z_val <- max(rank_ppi$overall_z, na.rm = TRUE)
limit_lower <- if (min_z_val < 0) min_z_val * 1.30 - 2 else 0
limit_upper <- if (max_z_val > 0) max_z_val * 1.20 else 1

p_bar_ppi <- ggplot(rank_ppi, aes(x = overall_z, y = network)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  geom_col(fill = "#8E6BBE", width = 0.7, color = "black", linewidth = 0.4) +
  geom_text(aes(label = paste0("#", Rank),
                hjust = ifelse(overall_z >= 0, -0.2, 1.2)),
            size = 5, color = "black") +
  scale_x_continuous(name = "Median Z-score",
                     limits = c(limit_lower, limit_upper),
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  nature_theme +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dotted", linewidth = 0.4),
    plot.margin = margin(5, 25, 5, 0)
  )

# 5.6 单独导出图例（14号 黑色 不加粗）
legend_theme <- theme(
  legend.position = "top",
  legend.title = element_text(size = 14, color = "black", face = "plain", lineheight = 0.9),
  legend.text = element_text(size = 14, color = "black", face = "plain"),
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0)
)

leg_egad <- cowplot::ggdraw(cowplot::get_legend(p_heat_egad + legend_theme))
ggsave("Legend_Standalone_EGAD.pdf", plot = leg_egad, width = 4.5, height = 1.2)

leg_ppi <- cowplot::ggdraw(cowplot::get_legend(p_heat_ppi + legend_theme))
ggsave("Legend_Standalone_PPI.pdf", plot = leg_ppi, width = 4.5, height = 1.2)

# 5.7 拼接最终图（heatmap : 右侧图 = 4 : 1）
p_heat_egad_noleg <- p_heat_egad + theme(legend.position = "none")
p_heat_ppi_noleg  <- p_heat_ppi  + theme(legend.position = "none")

p_final_egad <- p_heat_egad_noleg + p_box_egad + plot_layout(widths = c(4, 1))
ggsave("EGAD_AUC_rank_final.pdf", plot = p_final_egad, width = 10, height = 6.5)

p_final_ppi <- p_heat_ppi_noleg + p_bar_ppi + plot_layout(widths = c(4, 1))
ggsave("PPI_Zscore_rank_final.pdf", plot = p_final_ppi, width = 10, height = 6.5)

# 5.8 保存 benchmark 的 clean 数据
save(egad_clean, ppi_clean, heat_egad, heat_ppi, rank_egad, rank_ppi,
     file = "benchmark_8networks_cor_matrix_single.Rdata")   # 检查点 2

cat("✅ 主图全部完成（Pearson = cor_matrix_single，最终保留 8 种网络）！\n")

# ----------------------------------------------------------------------------
# 步骤 6：补充图 —— 各数据集的 EGAD 小提琴图（8 种网络固定 Nature 配色）
# ----------------------------------------------------------------------------
# 断点续跑：若单独运行本段，请取消下面一行注释
# load('benchmark.Rdata')   # 需包含 all_benchmark_results
# 直接复用步骤 2 已算好的 EGAD 结果，无需重新计算

keep_networks_sup <- c("Bayescor3", "Bayescor2", "Bayescor1", "pCor",
                       "hdWGCNA", "Cor", "Hotspot", "Random")

nature_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                   "#F39B7F", "#8491B4", "#7E6148", "#B09C85")
names(nature_colors) <- keep_networks_sup

network_mapping_sup <- c(
  "cor_matrix_single" = "Cor", "cor" = "Cor",
  "pcor" = "pCor", "partial correlation" = "pCor",
  "bayes_cor1" = "Bayescor1", "bayes_cor2" = "Bayescor2", "bayes_cor3" = "Bayescor3",
  "hotspot1" = "Hotspot", "hotspot" = "Hotspot",
  "hdwgcna" = "hdWGCNA", "hdWGCNA" = "hdWGCNA",
  "random" = "Random", "random network" = "Random"
)

benchmark_data_dir <- 'E:/deskup/毕业论文/bechmark/edga_ppi/data/'

for (ds in names(all_benchmark_results)) {
  cat("\n🚀 正在处理数据集:", ds, "\n")
  
  plot_df <- all_benchmark_results[[ds]]$egad
  if (nrow(plot_df) == 0) {
    cat("⚠️ 当前数据集有效计算结果为空，跳过绘图...\n")
    next
  }
  
  # 清洗 + 按中位数升序排列（Rank 低的放左边）
  plot_df_clean <- plot_df %>%
    filter(network %in% names(network_mapping_sup)) %>%
    mutate(network = recode(network, !!!network_mapping_sup)) %>%
    filter(network %in% keep_networks_sup, !is.na(nv_auc))
  
  auc_summary <- plot_df_clean %>%
    group_by(network) %>%
    summarise(median_auc = median(nv_auc, na.rm = TRUE), .groups = "drop") %>%
    arrange(median_auc)
  
  plot_df_clean$network <- factor(plot_df_clean$network, levels = auc_summary$network)
  
  # Nature 风格小提琴图（Arial 7号 纯黑 不加粗）
  p_egad_vln <- ggplot(plot_df_clean, aes(x = network, y = nv_auc, fill = network)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.5, color = "grey45") +
    geom_violin(adjust = 3, color = "black", linewidth = 0.4, trim = TRUE,
                width = 0.85, alpha = 0.85) +
    geom_boxplot(width = 0.16, fill = "white", color = "black",
                 linewidth = 0.4, outlier.shape = NA) +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.5,
                 fill = "black", color = "black") +
    scale_y_continuous(name = " ", expand = expansion(mult = c(0.02, 0.04))) +
    scale_fill_manual(values = nature_colors) +
    labs(x = NULL, y = " ") +
    theme_classic() +
    theme(
      text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
      axis.text.x = element_text(size = 7, color = "black", face = "plain",
                                 angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 7, color = "black", face = "plain"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 7, color = "black", face = "plain"),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # 输出到各数据集目录（如需统一输出，可改这里的目录）
  setwd(file.path('E:/deskup/hotspot/', ds))
  output_pdf_name <- paste0("EGAD_Violin_", ds, ".pdf")
  ggsave(output_pdf_name, plot = p_egad_vln, width = 5, height = 2.2, device = cairo_pdf)
  cat("✅ 已保存:", output_pdf_name, "\n")
  
  setwd(benchmark_data_dir)
}

cat("\n🎉 所有数据集的 EGAD 小提琴图已全部生成完毕！\n")


############################################################
# MescNet Benchmark
# EGAD + PPI benchmark visualization
#
# Pearson correlation = cor_matrix_single
#
# 最终保留 8 种网络：
# 1. pearson correlation
# 2. partial correlation
# 3. bayes_cor1
# 4. bayes_cor2
# 5. bayes_cor3
# 6. hotspot
# 7. hdWGCNA
# 8. random network
############################################################

rm(list = ls())

library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

# ==========================================================
# 0. 加载数据
# ==========================================================

setwd("E:/deskup/毕业论文/bechmark/edga_ppi/data/")

load("benchmark.Rdata")

cat("\n========================================\n")
cat("原始 EGAD network:\n")
print(unique(combined_egad$network))

cat("\n原始 PPI network:\n")
print(unique(combined_ppi$network))

cat("\n========================================\n")
cat("开始重新清洗数据...\n")


# ==========================================================
# 1. 数据集名称映射
# ==========================================================

dataset_mapping <- c(
  "scp1039" = "BRCA",
  "scp1265" = "TAA",
  "scp1289" = "COVID19",
  "scp1303" = "CM",
  "scp1526" = "ISL",
  "scp1671" = "HBM",
  "scp1731" = "HBA",
  "scp1849" = "ICM",
  "scp1852" = "CHD",
  "scp1903" = "ADI"
)


# ==========================================================
# 2. 网络名称映射
# ==========================================================
#
# 重点：
# cor_matrix_single = 你要的 Pearson correlation
# cor_matrix        = 不要
#
# ==========================================================

network_mapping <- c(
  
  # Pearson
  "cor_matrix_single" = "pearson correlation",
  
  # Partial correlation
  "pcor" = "partial correlation",
  
  # Bayes correlation
  "bayes_cor1" = "bayes_cor1",
  "bayes_cor2" = "bayes_cor2",
  "bayes_cor3" = "bayes_cor3",
  
  # Hotspot
  "hotspot1" = "hotspot",
  
  # hdWGCNA
  "hdwgcna" = "hdWGCNA",
  
  # Random
  "random" = "random network"
)


# ==========================================================
# 3. 最终保留的 8 个算法
# ==========================================================

keep_networks <- c(
  "pearson correlation",
  "partial correlation",
  "bayes_cor1",
  "bayes_cor2",
  "bayes_cor3",
  "hotspot",
  "hdWGCNA",
  "random network"
)


# ==========================================================
# 4. 清洗 EGAD
# ==========================================================

egad_clean <- combined_egad %>%
  
  # 先只保留我们真正需要的原始网络
  filter(
    network %in% names(network_mapping)
  ) %>%
  
  mutate(
    
    # 网络重命名
    network = recode(
      network,
      !!!network_mapping
    ),
    
    # 数据集重命名
    Dataset = recode(
      Dataset,
      !!!dataset_mapping
    )
  ) %>%
  
  # 再保险
  filter(
    network %in% keep_networks
  )


# ==========================================================
# 5. 清洗 PPI
# ==========================================================

ppi_clean <- combined_ppi %>%
  
  filter(
    network %in% names(network_mapping)
  ) %>%
  
  mutate(
    
    network = recode(
      network,
      !!!network_mapping
    ),
    
    Dataset = recode(
      Dataset,
      !!!dataset_mapping
    )
  ) %>%
  
  filter(
    network %in% keep_networks
  )


# ==========================================================
# 6. 检查清洗结果
# ==========================================================

cat("\n========================================\n")
cat("清洗后的 EGAD network:\n")
print(unique(egad_clean$network))

cat("\n清洗后的 PPI network:\n")
print(unique(ppi_clean$network))

cat("\nEGAD 每种网络数据量:\n")
print(table(egad_clean$network))

cat("\nPPI 每种网络数据量:\n")
print(table(ppi_clean$network))

cat("\nEGAD 数据集:\n")
print(unique(egad_clean$Dataset))

cat("\nPPI 数据集:\n")
print(unique(ppi_clean$Dataset))


# 必须确认最终是 8 种算法
if(length(unique(egad_clean$network)) != 8){
  warning("⚠️ EGAD 最终不是 8 种网络，请检查！")
}

if(length(unique(ppi_clean$network)) != 8){
  warning("⚠️ PPI 最终不是 8 种网络，请检查！")
}


# ==========================================================
# 7. 聚合 EGAD
# ==========================================================

cat("\n========================================\n")
cat("正在计算 EGAD median AUC...\n")

heat_egad <- egad_clean %>%
  group_by(
    network,
    Dataset
  ) %>%
  summarise(
    median_auc = median(
      nv_auc,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ==========================================================
# 8. 聚合 PPI
# ==========================================================

cat("正在计算 PPI median Z-score...\n")

heat_ppi <- ppi_clean %>%
  group_by(
    network,
    Dataset
  ) %>%
  summarise(
    median_z = median(
      z_score,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ==========================================================
# 9. EGAD 全局排名
# ==========================================================

rank_egad <- heat_egad %>%
  
  group_by(network) %>%
  
  summarise(
    overall_auc = median(
      median_auc,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  
  arrange(
    desc(overall_auc)
  ) %>%
  
  mutate(
    Rank = row_number()
  )


cat("\n========================================\n")
cat("EGAD global ranking:\n")
print(rank_egad)


# ==========================================================
# 10. PPI 全局排名
# ==========================================================

rank_ppi <- heat_ppi %>%
  
  group_by(network) %>%
  
  summarise(
    overall_z = median(
      median_z,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  
  arrange(
    desc(overall_z)
  ) %>%
  
  mutate(
    Rank = row_number()
  )


cat("\nPPI global ranking:\n")
print(rank_ppi)


# ==========================================================
# 11. 网络纵轴排序
# ==========================================================
#
# Rank 1 放最上面
#
# ==========================================================

egad_network_order <- rank_egad$network

ppi_network_order <- rank_ppi$network


# ==========================================================
# 12. 数据集横轴排序
# ==========================================================
#
# 按照 PPI Rank 1 算法的 Z-score
# 从左到右：低 -> 高
#
# ==========================================================

top1_ppi_algo <- ppi_network_order[1]

cat("\n========================================\n")
cat("PPI Rank #1 algorithm:\n")
print(top1_ppi_algo)


dataset_order <- heat_ppi %>%
  
  filter(
    network == top1_ppi_algo
  ) %>%
  
  arrange(
    median_z
  ) %>%
  
  pull(Dataset) %>%
  
  as.character()


# 防止某些数据集缺失
all_datasets <- unique(
  c(
    as.character(egad_clean$Dataset),
    as.character(ppi_clean$Dataset)
  )
)

dataset_order <- c(
  dataset_order,
  setdiff(
    all_datasets,
    dataset_order
  )
)


cat("\n最终 Dataset 排序:\n")
print(dataset_order)


# ==========================================================
# 13. 设置 factor 顺序
# ==========================================================

# -----------------------------
# EGAD
# -----------------------------

egad_clean$network <- factor(
  egad_clean$network,
  levels = rev(egad_network_order)
)

heat_egad$network <- factor(
  heat_egad$network,
  levels = rev(egad_network_order)
)

rank_egad$network <- factor(
  rank_egad$network,
  levels = rev(egad_network_order)
)

egad_clean$Dataset <- factor(
  egad_clean$Dataset,
  levels = dataset_order
)

heat_egad$Dataset <- factor(
  heat_egad$Dataset,
  levels = dataset_order
)


# -----------------------------
# PPI
# -----------------------------

ppi_clean$network <- factor(
  ppi_clean$network,
  levels = rev(ppi_network_order)
)

heat_ppi$network <- factor(
  heat_ppi$network,
  levels = rev(ppi_network_order)
)

rank_ppi$network <- factor(
  rank_ppi$network,
  levels = rev(ppi_network_order)
)

ppi_clean$Dataset <- factor(
  ppi_clean$Dataset,
  levels = dataset_order
)

heat_ppi$Dataset <- factor(
  heat_ppi$Dataset,
  levels = dataset_order
)


# ==========================================================
# 14. Nature 风格主题
# ==========================================================

nature_theme <- theme_classic() +
  
  theme(
    
    text = element_text(
      size = 14,
      color = "black",
      face = "plain"
    ),
    
    axis.text.x = element_text(
      size = 14,
      color = "black",
      face = "plain",
      angle = 45,
      hjust = 1
    ),
    
    axis.text.y = element_text(
      size = 14,
      color = "black",
      face = "plain"
    ),
    
    axis.title = element_text(
      size = 14,
      color = "black",
      face = "plain"
    ),
    
    axis.line = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    legend.text = element_text(
      size = 14,
      color = "black",
      face = "plain"
    ),
    
    legend.title = element_text(
      size = 14,
      color = "black",
      face = "plain"
    )
  )


heatmap_theme <- nature_theme +
  
  theme(
    
    axis.title = element_blank(),
    
    axis.line = element_blank(),
    
    axis.ticks = element_blank(),
    
    legend.position = "top"
  )


# ==========================================================
# 15. EGAD Heatmap
# ==========================================================

cat("\n🎨 正在绘制 EGAD heatmap...\n")


p_heat_egad <- ggplot(
  heat_egad,
  aes(
    x = Dataset,
    y = network,
    fill = median_auc
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  
  scale_fill_gradientn(
    
    colors = c(
      "#0B1D4A",
      "#265A8C",
      "#5895A9",
      "#A2BCA4",
      "#D5B08A",
      "#FCE2D4"
    ),
    
    name = "Median\nAUC"
  ) +
  
  scale_x_discrete(
    expand = c(0, 0)
  ) +
  
  scale_y_discrete(
    expand = c(0, 0)
  ) +
  
  coord_fixed() +
  
  heatmap_theme


# ==========================================================
# 16. EGAD 右侧 Boxplot
# ==========================================================

label_x_egad <- max(
  egad_clean$nv_auc,
  na.rm = TRUE
) + 0.05


x_min_egad <- min(
  egad_clean$nv_auc,
  na.rm = TRUE
)


p_box_egad <- ggplot(
  egad_clean,
  aes(
    x = nv_auc,
    y = network
  )
) +
  
  geom_boxplot(
    
    fill = "#C9B28C",
    
    color = "black",
    
    linewidth = 0.4,
    
    width = 0.65,
    
    outlier.shape = NA
  ) +
  
  geom_text(
    
    data = rank_egad,
    
    aes(
      x = label_x_egad,
      y = network,
      label = paste0("#", Rank)
    ),
    
    inherit.aes = FALSE,
    
    size = 5,
    
    color = "black"
  ) +
  
  scale_x_continuous(
    
    name = "AUC",
    
    limits = c(
      x_min_egad - 0.02,
      label_x_egad + 0.05
    ),
    
    expand = c(0, 0)
  ) +
  
  nature_theme +
  
  theme(
    
    axis.title.y = element_blank(),
    
    axis.text.y = element_blank(),
    
    axis.ticks.y = element_blank(),
    
    panel.grid.major.y = element_line(
      color = "grey80",
      linetype = "dotted",
      linewidth = 0.4
    ),
    
    plot.margin = margin(
      5,
      25,
      5,
      0
    )
  )


# ==========================================================
# 17. PPI Heatmap
# ==========================================================

cat("🎨 正在绘制 PPI heatmap...\n")


p_heat_ppi <- ggplot(
  heat_ppi,
  aes(
    x = Dataset,
    y = network,
    fill = median_z
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  
  scale_fill_gradient(
    
    low = "white",
    
    high = "#542788",
    
    name = "Median\nZ-score"
  ) +
  
  scale_x_discrete(
    expand = c(0, 0)
  ) +
  
  scale_y_discrete(
    expand = c(0, 0)
  ) +
  
  coord_fixed() +
  
  heatmap_theme


# ==========================================================
# 18. PPI 右侧 Rank barplot
# ==========================================================

min_z_val <- min(
  rank_ppi$overall_z,
  na.rm = TRUE
)

max_z_val <- max(
  rank_ppi$overall_z,
  na.rm = TRUE
)


# 给负值和 rank label 留空间
if(min_z_val < 0){
  
  limit_lower <- min_z_val * 1.30 - 2
  
} else {
  
  limit_lower <- 0
}


if(max_z_val > 0){
  
  limit_upper <- max_z_val * 1.20
  
} else {
  
  limit_upper <- 1
}


p_bar_ppi <- ggplot(
  rank_ppi,
  aes(
    x = overall_z,
    y = network
  )
) +
  
  geom_vline(
    xintercept = 0,
    color = "black",
    linewidth = 0.4
  ) +
  
  geom_col(
    
    fill = "#8E6BBE",
    
    width = 0.7,
    
    color = "black",
    
    linewidth = 0.4
  ) +
  
  geom_text(
    
    aes(
      
      label = paste0("#", Rank),
      
      hjust = ifelse(
        overall_z >= 0,
        -0.2,
        1.2
      )
    ),
    
    size = 5,
    
    color = "black"
  ) +
  
  scale_x_continuous(
    
    name = "Median Z-score",
    
    limits = c(
      limit_lower,
      limit_upper
    ),
    
    expand = c(0, 0)
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  nature_theme +
  
  theme(
    
    axis.title.y = element_blank(),
    
    axis.text.y = element_blank(),
    
    axis.ticks.y = element_blank(),
    
    panel.grid.major.y = element_line(
      color = "grey80",
      linetype = "dotted",
      linewidth = 0.4
    ),
    
    plot.margin = margin(
      5,
      25,
      5,
      0
    )
  )


# ==========================================================
# 19. 单独导出图例
# ==========================================================

cat("\n🎨 正在导出独立图例...\n")


leg_egad <- cowplot::ggdraw(
  cowplot::get_legend(
    p_heat_egad
  )
)


ggsave(
  "Legend_Standalone_EGAD.pdf",
  plot = leg_egad,
  width = 3,
  height = 1
)


leg_ppi <- cowplot::ggdraw(
  cowplot::get_legend(
    p_heat_ppi
  )
)


ggsave(
  "Legend_Standalone_PPI.pdf",
  plot = leg_ppi,
  width = 3,
  height = 1
)


# ==========================================================
# 20. 去掉主图 legend
# ==========================================================

p_heat_egad_noleg <- p_heat_egad +
  theme(
    legend.position = "none"
  )


p_heat_ppi_noleg <- p_heat_ppi +
  theme(
    legend.position = "none"
  )


# ==========================================================
# 21. 拼接 EGAD
# ==========================================================

p_final_egad <-
  
  p_heat_egad_noleg +
  
  p_box_egad +
  
  plot_layout(
    widths = c(4, 1)
  )


ggsave(
  
  "EGAD_AUC_rank_final.pdf",
  
  plot = p_final_egad,
  
  width = 10,
  
  height = 6.5
)


# ==========================================================
# 22. 拼接 PPI
# ==========================================================

p_final_ppi <-
  
  p_heat_ppi_noleg +
  
  p_bar_ppi +
  
  plot_layout(
    widths = c(4, 1)
  )


ggsave(
  
  "PPI_Zscore_rank_final.pdf",
  
  plot = p_final_ppi,
  
  width = 10,
  
  height = 6.5
)


# ==========================================================
# 23. 保存本次真正用于 benchmark 的 clean data
# ==========================================================

save(
  
  egad_clean,
  
  ppi_clean,
  
  heat_egad,
  
  heat_ppi,
  
  rank_egad,
  
  rank_ppi,
  
  file = "benchmark_8networks_cor_matrix_single.Rdata"
)


cat("\n========================================\n")
cat("✅ 全部完成！\n")
cat("✅ Pearson correlation = cor_matrix_single\n")
cat("✅ cor_matrix 已删除\n")
cat("✅ bayes_cor4/5/6 已删除\n")
cat("✅ hotspot2 已删除\n")
cat("✅ 最终保留 8 种算法\n")
cat("========================================\n")






library(ggplot2)
library(cowplot)
library(grid)

# ==========================================================
# 重新覆盖 p_heat_egad 的颜色图例
# ==========================================================

p_heat_egad <- p_heat_egad +
  scale_fill_gradientn(
    colors = c(
      "#0B1D4A",
      "#265A8C",
      "#5895A9",
      "#A2BCA4",
      "#D5B08A",
      "#FCE2D4"
    ),
    
    name = "Median\nAUC",
    
    # 只显示三个数字
    breaks = c(0.50, 0.60, 0.65),
    labels = c("0.50", "0.60", "0.65"),
    
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "left",
      title.hjust = 0.5,
      label.position = "bottom",
      
      barwidth = unit(6, "cm"),
      barheight = unit(0.6, "cm"),
      
      ticks = TRUE
    )
  )


# ==========================================================
# 提取新图例
# ==========================================================

leg_egad <- cowplot::get_legend(
  p_heat_egad +
    theme(
      legend.position = "top",
      
      legend.title = element_text(
        size = 14,
        color = "black",
        face = "plain",
        lineheight = 0.9
      ),
      
      legend.text = element_text(
        size = 14,
        color = "black",
        face = "plain"
      ),
      
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
)

leg_egad_plot <- ggdraw(leg_egad)

ggsave(
  "Legend_Standalone_EGAD.pdf",
  plot = leg_egad_plot,
  width = 4.5,
  height = 1.2
)















scale_fill_gradient(
  low = "white",
  high = "#542788",
  name = "Median\nZ-score"
)
library(ggplot2)
library(cowplot)
library(grid)

# ==========================================================
# 重新覆盖 PPI 图例，只显示 3 个刻度
# ==========================================================

ppi_min <- min(heat_ppi$median_z, na.rm = TRUE)
ppi_max <- max(heat_ppi$median_z, na.rm = TRUE)

ppi_breaks <- c(
  ppi_min,
  (ppi_min + ppi_max) / 2,
  ppi_max
)

p_heat_ppi <- p_heat_ppi +
  scale_fill_gradient(
    low = "white",
    high = "#542788",
    name = "Median\nZ-score",
    
    breaks = ppi_breaks,
    labels = round(ppi_breaks, 1),
    
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "left",
      title.hjust = 0.5,
      label.position = "bottom",
      barwidth = unit(6, "cm"),
      barheight = unit(0.6, "cm"),
      ticks = TRUE
    )
  )

# ==========================================================
# 提取新的 PPI 图例
# ==========================================================

leg_ppi <- cowplot::get_legend(
  p_heat_ppi +
    theme(
      legend.position = "top",
      
      legend.title = element_text(
        size = 14,
        color = "black",
        face = "plain",
        lineheight = 0.9
      ),
      
      legend.text = element_text(
        size = 14,
        color = "black",
        face = "plain"
      ),
      
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
)

leg_ppi_plot <- cowplot::ggdraw(leg_ppi)

ggsave(
  "Legend_Standalone_PPI.pdf",
  plot = leg_ppi_plot,
  width = 4.5,
  height = 1.2
)

#Benchmark1_egad
#运行全局变量和加载函数
# Load EGAD and the data files
library(ggplot2)
library(viridis)
#metacell
#先运行function

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(ppcor)
library(igraph)
library(dplyr)
library(stringr)
library(stringr)
library(grid)  # 处理自定义标签
library(org.Hs.eg.db) 
library(clusterProfiler) 
library(DOSE)
library(RColorBrewer) 
library(reshape2)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(scales)

# 加载必要的库
library(circlize)

library(tidyverse)


library(EGAD)













# ==============================================================================
# 加载必要的包
# ==============================================================================
library(Seurat)
library(tidyverse)
library(WGCNA)
library(igraph)
library(org.Hs.eg.db)
library(EGAD)
library(ggplot2)
library(viridis)

# ==============================================================================
# 模块 1: 纯计算函数 (不包含任何画图逻辑)
# ==============================================================================

# 计算 EGAD AUC 结果
calc_egad_auc <- function(TOMs, go_data) {
  
}





# ==============================================================================
# 模块 2: 加载全局数据并执行循环
# ==============================================================================
cat("⏳ 正在加载全局 GO 和 PPI 数据...\n")
setwd('E:/deskup/Mescnet/data/')
load('go.Rdata') # 提供 go 数据框


path <- 'E:/deskup/hotspot/'
files <- list.files(path)

# 初始化存储所有数据集结果的“大列表”
all_benchmark_results <- list()


i<-1




i<-1+i
dataset_name <- files[i]
path_i<-paste0(path,dataset_name)
setwd(path_i)
load('seurat_obj.Rdata')



library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrastr)
library(ggsci)   
library(cowplot) 

cat("⏳ 正在提取 UMAP 坐标并准备绘图数据...\n")

# ==============================================================================
# 1. 提取 UMAP 坐标和元数据
# ==============================================================================
embed_coords <- Embeddings(seurat_obj, "umap")
df <- data.frame(
  UMAP_1 = embed_coords[, 1],
  UMAP_2 = embed_coords[, 2],
  celltype = seurat_obj@meta.data$cell_type
)

# ==============================================================================
# 2. 动态生成配色 (自动适应不同数据集的细胞群数量)
# ==============================================================================
cell_types <- unique(df$celltype)
num_types <- length(cell_types)


# 使用 Nature 风格的 npg 调色板并进行自动扩展
dynamic_colors <- colorRampPalette(ggsci::pal_npg()(10))(num_types)
names(dynamic_colors) <- cell_types

# ==============================================================================
# 3. 统一全局主题 (主图隐藏图例)
# ==============================================================================
umap_nature_theme <- theme_void() + 
  theme(
    text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
    legend.position = "none", # 主图隐藏图例，保持纯净分布图
    plot.margin = margin(10, 10, 10, 10)
  )

# ==============================================================================
# 4. 绘制主图 (栅格化处理)
# ==============================================================================
cat("🎨 正在绘制栅格化 UMAP...\n")
p_umap <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
  # 将散点光栅化为 600 DPI 像素，Inkscape 打开秒顺滑，外框保留矢量
  rasterise(geom_point(aes(color = celltype), size = 0.2, alpha = 0.8), dpi = 600) +
  scale_color_manual(values = dynamic_colors) +
  coord_fixed() +
  umap_nature_theme

# 保存主图
output_main <- paste0("UMAP_", dataset_name, "_Raster.pdf")
ggsave(output_main, plot = p_umap, width = 5, height = 5, device = cairo_pdf)

# ==============================================================================
# 5. 单独提取并导出图例 (Arial, 7号, 纯黑, 不加粗)
# ==============================================================================
cat("🎨 正在提取独立矢量图例...\n")
p_for_legend <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = celltype), size = 3) +
  scale_color_manual(values = dynamic_colors, name = "Cell Type") +
  theme_void() +
  theme(
    # 精准控制：Arial 7号 纯黑 不加粗
    text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
    legend.text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
    legend.title = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
    legend.key.size = unit(0.4, "cm") # 缩小色块尺寸以适配 7 号小字体
  )

# 剥离图例并保存
leg <- get_legend(p_for_legend)
output_legend <- paste0("Legend_", dataset_name, ".pdf")
ggsave(output_legend, plot = leg, width = 2.5, height = 3, device = cairo_pdf)

cat(paste0("🎉 完成！主图 (", output_main, ") 和独立图例 (", output_legend, ") 均已完美导出！\n"))












library(Seurat)
library(tidyverse)
library(WGCNA)
library(igraph)
library(org.Hs.eg.db)
library(EGAD)
library(ggplot2)

# ==============================================================================
# 0. 环境与路径初始化
# ==============================================================================
cat("⏳ 正在加载全局 GO 数据...\n")
setwd('E:/deskup/Mescnet/data/')
load('go.Rdata') # 提供 go 数据框
goterms <- unique(go[, 3])

path <- 'E:/deskup/hotspot/'
files <- list.files(path)

# ==============================================================================
# 1. 核心修改：为 8 种固定算法指定全球唯一的 Nature 风格配色
# ==============================================================================
keep_networks <- c("Bayescor3", "Bayescor2", "Bayescor1", "pCor", "hdWGCNA", "Cor", "Hotspot", "Random")

nature_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488",  
  "#F39B7F", "#8491B4", "#7E6148", "#B09C85"
)
# 🌟 强行绑定：让每一个算法名字对应一个固定色号
names(nature_colors) <- keep_networks

# 网络名称映射字典
network_mapping <- c(
  "cor_matrix_single" = "Cor",
  "cor"               = "Cor",
  "pcor"              = "pCor",
  "partial correlation" = "pCor",
  "bayes_cor1"        = "Bayescor1",
  "bayes_cor2"        = "Bayescor2",
  "bayes_cor3"        = "Bayescor3",
  "hotspot1"          = "Hotspot",
  "hotspot"           = "Hotspot",
  "hdwgcna"           = "hdWGCNA",
  "hdWGCNA"           = "hdWGCNA",
  "random"            = "Random",
  "random network"    = "Random"
)


# ==============================================================================
# 2. 开始大循环：遍历所有数据集
# ==============================================================================
cat("⏳ 开始循环计算并绘制所有数据集的 EGAD 小提琴图...\n")

for (i in 1:length(files)) {
  dataset_name <- files[i]
  cat("\n--------------------------------------------------\n")
  cat("🚀 正在处理第", i, "/", length(files), "个数据集:", dataset_name, "\n")
  
  path_i <- file.path(path, dataset_name)
  
  # 检查并进入目录
  if(!dir.exists(path_i)) {
    cat("⚠️ 路径不存在，跳过:", path_i, "\n")
    next
  }
  setwd(path_i)
  
  if(!file.exists('TOMs.Rdata')) {
    cat("⚠️ 未找到 TOMs.Rdata，跳过...\n")
    next
  }
  
  load('TOMs.Rdata')
  
  # 追加额外的网络矩阵
  if (file.exists('output1.csv') & file.exists('output2.csv')) {
    hotspot1 <- read.csv('output1.csv')
    hotspot2 <- read.csv('output2.csv')
    rownames(hotspot1) <- colnames(hotspot1)
    rownames(hotspot2) <- colnames(hotspot2)
    TOMs[['hotspot1']] <- as.matrix(hotspot1)
    TOMs[['hotspot2']] <- as.matrix(hotspot2)
  }
  
  # --------------------------------------------------------------------------
  # 2.1 计算当前数据集的 EGAD 结果
  # --------------------------------------------------------------------------
  plot_df <- data.frame()
  
  for (net_name in names(TOMs)) {
    tom <- TOMs[[net_name]]
    gene1 <- rownames(tom)
    
    valid_symbols <- intersect(gene1, keys(org.Hs.eg.db, keytype = "SYMBOL"))
    valid_ensembl <- intersect(gsub("\\..*", "", gene1), keys(org.Hs.eg.db, keytype = "ENSEMBL"))
    valid_entrez <- intersect(gene1, keys(org.Hs.eg.db, keytype = "ENTREZID"))
    
    if (length(valid_symbols) > 0) {
      gene_symbols <- mapIds(org.Hs.eg.db, keys = gene1, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    } else if (length(valid_ensembl) > 0) {
      gene1_clean <- gsub("\\..*", "", gene1)
      gene_symbols <- mapIds(org.Hs.eg.db, keys = gene1_clean, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
      names(gene_symbols) <- gene1
    } else if (length(valid_entrez) > 0) {
      gene_symbols <- setNames(gene1, gene1)
    } else {
      next
    }
    
    result_df <- data.frame(EntrezID = gene_symbols[gene1], Symbol = gene1)
    gene_table <- na.omit(result_df)
    
    if (nrow(gene_table) < 10) next
    
    tom[lower.tri(tom)] <- t(tom)[lower.tri(tom)]
    
    genelist <- rownames(tom)
    genelist <- genelist[(genelist %in% gene_table$Symbol)]
    tom <- tom[genelist, genelist]
    
    ix <- match(genelist, gene_table$Symbol)
    gene_ids <- gene_table$EntrezID[ix]
    colnames(tom) <- gene_ids
    rownames(tom) <- gene_ids
    
    annotations <- make_annotations(go[, c(2, 3)], gene_ids, goterms)
    
    tryCatch({
      GO_groups_voted <- run_GBA(tom, annotations, max = Inf)
      df <- data.frame(
        nv_auc = as.numeric(GO_groups_voted[[1]][, 1]),
        nd_auc = as.numeric(GO_groups_voted[[1]][, 3]),
        network = net_name
      )
      plot_df <- rbind(plot_df, df)
    }, error = function(e) {})
  }
  
  if(nrow(plot_df) == 0) {
    cat("⚠️ 当前数据集有效计算结果为空，跳过绘图...\n")
    next
  }
  
  # --------------------------------------------------------------------------
  # 2.2 数据清洗与升序排序
  # --------------------------------------------------------------------------
  plot_df_clean <- plot_df %>%
    filter(network %in% names(network_mapping)) %>%
    mutate(network = recode(network, !!!network_mapping)) %>%
    filter(network %in% keep_networks, !is.na(nv_auc))
  
  auc_summary <- plot_df_clean %>%
    group_by(network) %>%
    summarise(median_auc = median(nv_auc, na.rm = TRUE), .groups = "drop") %>%
    arrange(median_auc) # 升序排列
  
  plot_df_clean$network <- factor(
    plot_df_clean$network,
    levels = auc_summary$network
  )
  
  # --------------------------------------------------------------------------
  # 2.3 绘制 Nature 风格圆润小提琴图（颜色自动根据绑定的名字匹配）
  # --------------------------------------------------------------------------
  p_egad_vln <- ggplot(plot_df_clean, aes(x = network, y = nv_auc, fill = network)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.5, color = "grey45") +
    geom_violin(adjust = 3, color = "black", linewidth = 0.4, trim = TRUE, width = 0.85, alpha = 0.85) +
    geom_boxplot(width = 0.16, fill = "white", color = "black", linewidth = 0.4, outlier.shape = NA) +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.5, fill = "black", color = "black") +
    scale_y_continuous(name = " ", expand = expansion(mult = c(0.02, 0.04))) +
    scale_fill_manual(values = nature_colors) + # 此处会根据命名好的 vector 自动对号入座
    labs(x = NULL, y = " ") +
    theme_classic() +
    theme(
      text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
      axis.text.x = element_text(size = 7, color = "black", face = "plain", angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 7, color = "black", face = "plain"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 7, color = "black", face = "plain"),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # --------------------------------------------------------------------------
  # 2.4 保存为独立 PDF
  # --------------------------------------------------------------------------
  output_pdf_name <- paste0("EGAD_Violin_", dataset_name, ".pdf")
  ggsave(
    output_pdf_name,
    plot = p_egad_vln,
    width = 5,
    height = 2.2,
    device = cairo_pdf
  )
  
  cat("✅ 已成功保存:", output_pdf_name, "\n")
}

cat("\n🎉 所有的单数据集 EGAD 小提琴图已全部循环生成完毕，且各算法颜色已完美保持全局一致！\n")




















library(Seurat)
library(tidyverse)
library(WGCNA)
library(igraph)
library(org.Hs.eg.db)
library(EGAD)
library(ggplot2)

# ==============================================================================
# 0. 环境与路径初始化
# ==============================================================================
cat("⏳ 正在加载全局 GO 数据...\n")
setwd('E:/deskup/Mescnet/data/')
load('go.Rdata') # 提供 go 数据框
goterms <- unique(go[, 3])

path <- 'E:/deskup/hotspot/'
files <- list.files(path)

# ==============================================================================
# 1. 核心修改：为 8 种固定算法指定全球唯一的 Nature 风格配色
# ==============================================================================
keep_networks <- c("Bayescor3", "Bayescor2", "Bayescor1", "pCor", "hdWGCNA", "Cor", "Hotspot", "Random")

nature_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488",  
  "#F39B7F", "#8491B4", "#7E6148", "#B09C85"
)
# 🌟 强行绑定：让每一个算法名字对应一个固定色号
names(nature_colors) <- keep_networks

# 网络名称映射字典
network_mapping <- c(
  "cor_matrix_single" = "Cor",
  "cor"               = "Cor",
  "pcor"              = "pCor",
  "partial correlation" = "pCor",
  "bayes_cor1"        = "Bayescor1",
  "bayes_cor2"        = "Bayescor2",
  "bayes_cor3"        = "Bayescor3",
  "hotspot1"          = "Hotspot",
  "hotspot"           = "Hotspot",
  "hdwgcna"           = "hdWGCNA",
  "hdWGCNA"           = "hdWGCNA",
  "random"            = "Random",
  "random network"    = "Random"
)



# ==============================================================================
# 2. 开始大循环：遍历所有数据集
# ==============================================================================
cat("⏳ 开始循环计算并绘制所有数据集的 EGAD 小提琴图...\n")

for (i in 1:length(files)) {
  dataset_name <- files[i]
  cat("\n--------------------------------------------------\n")
  cat("🚀 正在处理第", i, "/", length(files), "个数据集:", dataset_name, "\n")
  
  path_i <- file.path(path, dataset_name)
  
  # 检查并进入目录
  if(!dir.exists(path_i)) {
    cat("⚠️ 路径不存在，跳过:", path_i, "\n")
    next
  }
  setwd(path_i)
  
  if(!file.exists('TOMs.Rdata')) {
    cat("⚠️ 未找到 TOMs.Rdata，跳过...\n")
    next
  }
  
  load('TOMs.Rdata')
  
  # 追加额外的网络矩阵
  if (file.exists('output1.csv') & file.exists('output2.csv')) {
    hotspot1 <- read.csv('output1.csv')
    hotspot2 <- read.csv('output2.csv')
    rownames(hotspot1) <- colnames(hotspot1)
    rownames(hotspot2) <- colnames(hotspot2)
    TOMs[['hotspot1']] <- as.matrix(hotspot1)
    TOMs[['hotspot2']] <- as.matrix(hotspot2)
  }
  
  # --------------------------------------------------------------------------
  # 2.1 计算当前数据集的 EGAD 结果
  # --------------------------------------------------------------------------
  plot_df <- data.frame()
  
  for (net_name in names(TOMs)) {
    tom <- TOMs[[net_name]]
    gene1 <- rownames(tom)
    
    valid_symbols <- intersect(gene1, keys(org.Hs.eg.db, keytype = "SYMBOL"))
    valid_ensembl <- intersect(gsub("\\..*", "", gene1), keys(org.Hs.eg.db, keytype = "ENSEMBL"))
    valid_entrez <- intersect(gene1, keys(org.Hs.eg.db, keytype = "ENTREZID"))
    
    if (length(valid_symbols) > 0) {
      gene_symbols <- mapIds(org.Hs.eg.db, keys = gene1, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    } else if (length(valid_ensembl) > 0) {
      gene1_clean <- gsub("\\..*", "", gene1)
      gene_symbols <- mapIds(org.Hs.eg.db, keys = gene1_clean, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
      names(gene_symbols) <- gene1
    } else if (length(valid_entrez) > 0) {
      gene_symbols <- setNames(gene1, gene1)
    } else {
      next
    }
    
    result_df <- data.frame(EntrezID = gene_symbols[gene1], Symbol = gene1)
    gene_table <- na.omit(result_df)
    
    if (nrow(gene_table) < 10) next
    
    tom[lower.tri(tom)] <- t(tom)[lower.tri(tom)]
    
    genelist <- rownames(tom)
    genelist <- genelist[(genelist %in% gene_table$Symbol)]
    tom <- tom[genelist, genelist]
    
    ix <- match(genelist, gene_table$Symbol)
    gene_ids <- gene_table$EntrezID[ix]
    colnames(tom) <- gene_ids
    rownames(tom) <- gene_ids
    
    annotations <- make_annotations(go[, c(2, 3)], gene_ids, goterms)
    
    tryCatch({
      GO_groups_voted <- run_GBA(tom, annotations, max = Inf)
      df <- data.frame(
        nv_auc = as.numeric(GO_groups_voted[[1]][, 1]),
        nd_auc = as.numeric(GO_groups_voted[[1]][, 3]),
        network = net_name
      )
      plot_df <- rbind(plot_df, df)
    }, error = function(e) {})
  }
  
  if(nrow(plot_df) == 0) {
    cat("⚠️ 当前数据集有效计算结果为空，跳过绘图...\n")
    next
  }
  
  # --------------------------------------------------------------------------
  # 2.2 数据清洗与升序排序
  # --------------------------------------------------------------------------
  plot_df_clean <- plot_df %>%
    filter(network %in% names(network_mapping)) %>%
    mutate(network = recode(network, !!!network_mapping)) %>%
    filter(network %in% keep_networks, !is.na(nv_auc))
  
  auc_summary <- plot_df_clean %>%
    group_by(network) %>%
    summarise(median_auc = median(nv_auc, na.rm = TRUE), .groups = "drop") %>%
    arrange(median_auc) # 升序排列
  
  plot_df_clean$network <- factor(
    plot_df_clean$network,
    levels = auc_summary$network
  )
  
  # --------------------------------------------------------------------------
  # 2.3 绘制 Nature 风格圆润小提琴图（颜色自动根据绑定的名字匹配）
  # --------------------------------------------------------------------------
  p_egad_vln <- ggplot(plot_df_clean, aes(x = network, y = nv_auc, fill = network)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.5, color = "grey45") +
    geom_violin(adjust = 3, color = "black", linewidth = 0.4, trim = TRUE, width = 0.85, alpha = 0.85) +
    geom_boxplot(width = 0.16, fill = "white", color = "black", linewidth = 0.4, outlier.shape = NA) +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.5, fill = "black", color = "black") +
    scale_y_continuous(name = " ", expand = expansion(mult = c(0.02, 0.04))) +
    scale_fill_manual(values = nature_colors) + # 此处会根据命名好的 vector 自动对号入座
    labs(x = NULL, y = " ") +
    theme_classic() +
    theme(
      text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
      axis.text.x = element_text(size = 7, color = "black", face = "plain", angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 7, color = "black", face = "plain"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 7, color = "black", face = "plain"),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # --------------------------------------------------------------------------
  # 2.4 保存为独立 PDF
  # --------------------------------------------------------------------------
  output_pdf_name <- paste0("EGAD_Violin_", dataset_name, ".pdf")
  ggsave(
    output_pdf_name,
    plot = p_egad_vln,
    width = 5,
    height = 2.2,
    device = cairo_pdf
  )
  
  cat("✅ 已成功保存:", output_pdf_name, "\n")
}

cat("\n🎉 所有的单数据集 EGAD 小提琴图已全部循环生成完毕，且各算法颜色已完美保持全局一致！\n")