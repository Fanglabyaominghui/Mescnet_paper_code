suppressPackageStartupMessages({
  library(Seurat)
  library(schard)
  library(UCell)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(patchwork)
  library(ggrastr) # 🌟 新增：专门用于栅格化密集的散点图层
})

# ==========================================
# 0. 全局环境与字体设置 (全局强制：14pt 纯黑不加粗)
# ==========================================
available_fonts <- names(pdfFonts())
if ("Arial" %in% available_fonts) {
  family_use <- "Arial"
} else if ("ArialMT" %in% available_fonts) {
  family_use <- "ArialMT"
} else if ("Helvetica" %in% available_fonts) {
  family_use <- "Helvetica"
} else {
  family_use <- "sans"
}
cat("✅ 使用字体:", family_use, "\n")

# 强制三等分断点函数，保证左、中、右绝对对齐，只显示3个数字
exact_3_breaks <- function(x) c(x[1], (x[1] + x[2]) / 2, x[2])

# ==========================================
# 1. 核心单图绘制函数 (栅格化散点 + q70阈值)
# ==========================================
plot_single_module <- function(df, mod_col, mod_name, color_high) {
  
  # 🌟 严格计算 q70 - q95，强力过滤背景底噪
  q_low <- quantile(df[[mod_col]], 0.70, na.rm = TRUE)
  q_high <- quantile(df[[mod_col]], 0.95, na.rm = TRUE)
  
  # 独立横向图例 (Bar图)，放在图底
  my_guide <- guide_colorbar(
    barwidth = unit(4.5, "cm"),
    barheight = unit(0.5, "cm"),
    title.position = "top",
    title.hjust = 0.5,
    label.position = "bottom",
    direction = "horizontal"
  )
  
  # 按表达量排序，保证高表达点在最上层
  df_sorted <- df %>% arrange(.data[[mod_col]])
  
  ggplot(df_sorted, aes(x = Dim_1, y = Dim_2, color = .data[[mod_col]])) +
    # 🌟 核心替换：使用 geom_point_rast 替代 geom_point，并设置 600 DPI 高清栅格化
    geom_point_rast(size = 0.4, stroke = 0, raster.dpi = 600) +
    scale_color_gradient(
      low = "#E8E8E8",
      high = color_high,
      name = mod_name,
      limits = c(q_low, q_high),
      oob = squish,
      breaks = exact_3_breaks,
      labels = function(x) sprintf("%.2f", x),
      guide = my_guide
    ) +
    coord_fixed() +
    theme_void() +
    theme(
      text = element_text(size = 14, color = "black", face = "plain", family = family_use),
      plot.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 14, color = "black", face = "plain", family = family_use),
      legend.title = element_text(size = 14, color = "black", face = "plain", family = family_use),
      legend.margin = margin(t = 5, b = 10),
      plot.margin = margin(10, 10, 10, 10),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# ==========================================
# 2. 批量处理切片文件流水线
# ==========================================
# 强制为 filtered_signatures 命名，与提取列名完美对应
names(filtered_signatures) <- paste0("Module_", 1:4)

# 仅保留 T353 和 T354
files_list <- c("Mouse2_T353.h5ad", "Mouse2_T354.h5ad")

# 颜色设定
colors_list <- c(
  "MOD1" = alpha("#aa363d", 0.9),   # 红色 (Module 1)
  "MOD2" = alpha("#228fbd", 0.9),   # 紫色 (Module 2)
  "MOD3" = alpha("#9b95c9", 0.95),  # 绿色 (Module 3)
  "MOD4" = alpha("#769149", 0.9)    # 蓝色 (Module 4)
)

for (f in files_list) {
  slice_name <- gsub(".h5ad", "", f)
  cat("\n========================================\n")
  cat("🚀 正在处理切片:", slice_name, "...\n")
  
  # A. 读取与计算
  obj <- schard::h5ad2seurat(f)
  cat("   ⏳ 正在计算 UCell 分数...\n")
  
  obj <- UCell::AddModuleScore_UCell(
    obj = obj,
    features = filtered_signatures,
    ncores = 24, 
    name = "_UCell"
  )
  
  # B. 提取数据
  cat("   🎨 正在提取数据并绘制栅格化 2x2 网格分面图...\n")
  plot_df <- data.frame(
    Dim_1 = obj@reductions$Xspatial_@cell.embeddings[, 1],
    Dim_2 = obj@reductions$Xspatial_@cell.embeddings[, 2],
    Module_1 = obj@meta.data$Module_1_UCell,
    Module_2 = obj@meta.data$Module_2_UCell,
    Module_3 = obj@meta.data$Module_3_UCell,
    Module_4 = obj@meta.data$Module_4_UCell
  )
  
  rm(obj)
  gc()
  
  # C. 独立生成 4 张图
  p1 <- plot_single_module(plot_df, "Module_1", "Module 1", colors_list["MOD1"])
  p2 <- plot_single_module(plot_df, "Module_2", "Module 2", colors_list["MOD2"])
  p3 <- plot_single_module(plot_df, "Module_3", "Module 3", colors_list["MOD3"])
  p4 <- plot_single_module(plot_df, "Module_4", "Module 4", colors_list["MOD4"])
  
  # D. Patchwork 完美拼接 (2x2网格)
  final_plot <- wrap_plots(p1, p2, p3, p4, ncol = 2)
  
  # E. 导出
  out_pdf <- paste0(slice_name, "_SeparateModules_Raster_Grid.pdf")
  
  cairo_pdf(out_pdf, height = 8.5, width = 8.5, family = family_use)
  print(final_plot)
  dev.off()
  
  cat("   💾", out_pdf, "已成功生成 (8.5×8.5英寸，q70过滤，散点已栅格化)！\n")
}

cat("\n🎉 双样本栅格化分离绘图完美执行完毕！\n")





#fig4

library(uwot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggrastr)

# ==========================================
# 0. 跨平台字体检测与全局 14pt 设置
# ==========================================
available_fonts <- names(pdfFonts())
family_use <- if ("Arial" %in% available_fonts) "Arial" else "sans"

# 全局强制 14号、纯黑、不加粗
font_style_14pt <- element_text(size = 14, color = "black", face = "plain", family = family_use)

# ==========================================
# 1. 核心降维函数 (无监督版)
# ==========================================
RunMEscnetUMAP <- function(tom, kme_list, n_hubs = 10, n_neighbors = 15, min_dist = 0.1, spread = 1) {
  
  cat("⏳ 正在整合 kME 列表并提取特征...\n")
  
  # A. 转换 list 为统一的 df，并将 Module_2-5 映射为 M1-M4
  gene_info <- bind_rows(lapply(names(kme_list), function(mod_name) {
    df <- kme_list[[mod_name]]
    new_mod <- paste0("M", as.numeric(gsub("Module_", "", mod_name)) - 1)
    df$module <- new_mod
    return(df)
  }))
  
  # B. 过滤出在 TOM 中真实存在的基因
  valid_genes <- intersect(gene_info$gene, rownames(tom))
  gene_info <- gene_info %>% filter(gene %in% valid_genes)
  
  # C. 提取 Hub 基因
  hub_genes <- gene_info %>%
    group_by(module) %>%
    slice_max(order_by = kME, n = n_hubs, with_ties = FALSE) %>%
    pull(gene)
  
  valid_hubs <- intersect(hub_genes, colnames(tom))
  
  # D. 构建特征矩阵
  feature_mat <- tom[valid_genes, valid_hubs, drop = FALSE]
  
  cat("🗺️ 正在执行无监督 UMAP 降维...\n")
  set.seed(42)
  
  hub_umap <- uwot::umap(
    X = feature_mat,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    spread = spread,
    metric = "cosine"
  )
  
  # E. 组装结果 dataframe
  plot_df <- data.frame(
    gene = valid_genes,
    UMAP1 = hub_umap[, 1],
    UMAP2 = hub_umap[, 2]
  ) %>%
    left_join(gene_info, by = "gene")
  
  # 对每个模块内部的 kME 进行缩放，控制点大小
  plot_df <- plot_df %>%
    group_by(module) %>%
    mutate(kME_scaled = (kME - min(kME)) / (max(kME) - min(kME) + 1e-6)) %>%
    ungroup()
  
  return(plot_df)
}

# ==========================================
# 2. 执行降维计算 (读取 bayes_cor2)
# ==========================================
tom_matrix <- seurat_obj@misc$MEscnet$bayes_cor2

umap_df <- RunMEscnetUMAP(
  tom = tom_matrix,
  kme_list = module_kME_list,
  n_hubs = 10,
  n_neighbors = 15,
  min_dist = 0.1,
  spread = 1
)

# ==========================================
# 3. 绘制 14pt 顶刊格式图表 (栅格化散点 + 隐去坐标轴)
# ==========================================
cat("🎨 正在生成纯净版栅格化无监督 UMAP 图...\n")

mod_palette <- c("M1" = "#aa363d", "M2" = "#228fbd", "M3" = "#9b95c9", "M4" = "#769149")

hubs_to_label <- umap_df %>%
  group_by(module) %>%
  slice_max(order_by = kME, n = 3, with_ties = FALSE)

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  # 开启 600 DPI 栅格化
  geom_point_rast(
    aes(color = module, size = kME_scaled * 2),
    alpha = 0.9,
    raster.dpi = 600,
    stroke = 0
  ) +
  geom_text_repel(
    data = hubs_to_label,
    aes(label = gene),
    size = 14 / .pt,
    fontface = "plain",
    family = family_use,
    color = "black",
    bg.color = "white",
    bg.r = 0.15,
    max.overlaps = Inf,
    segment.color = "grey50"
  ) +
  scale_color_manual(values = mod_palette) +
  scale_size_identity() +
  # 使用 theme_void() 彻底隐去坐标轴、刻度和背景网格
  theme_void() +
  theme(
    text = font_style_14pt,
    legend.text = font_style_14pt,
    legend.title = element_blank(),
    legend.position = "right",
    aspect.ratio = 1
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

# ==========================================
# 4. 导出 PDF
# ==========================================
out_pdf <- "Figure_Unsupervised_Module_UMAP_Raster_NoAxis.pdf"

cairo_pdf(out_pdf, width = 4.5, height = 3.5, family = family_use)
print(p_umap)
dev.off()

cat("✅ 搞定！去坐标轴、散点栅格化版的 UMAP 已经准备好被拖进 Inkscape 了！\n")










library(Seurat)
library(corrplot)

# ==========================================
# 0. 跨平台字体检测与全局 14pt 设置
# ==========================================
available_fonts <- names(pdfFonts())
if ("Arial" %in% available_fonts) {
  family_use <- "Arial"
} else if ("ArialMT" %in% available_fonts) {
  family_use <- "ArialMT"
} else if ("Helvetica" %in% available_fonts) {
  family_use <- "Helvetica"
} else {
  family_use <- "sans"
}
cat("✅ 使用字体:", family_use, "\n")

# ==========================================
# 1. 定义列名并提取、缩放 (Scale / Z-score)
# ==========================================
# 提取模块并重命名为 M1-M4
old_mods <- c("Module_1_UCell", "Module_1_UCell.1", "Module_2_UCell", "Module_3_UCell")
new_mods <- c("M1", "M2", "M3", "M4")
mat_modules <- scale(as.matrix(seurat_obj@meta.data[, old_mods]))
colnames(mat_modules) <- new_mods

# 提取疾病
dis_cols <- c("ad_UCell", "adca_UCell", "adhd_UCell", "asd_UCell", "bd_UCell", 
              "DWS_UCell", "dys_UCell", "et_UCell", "hb_UCell", "JS1_UCell", 
              "MDD_UCell", "msa_UCell", "pcd_UCell", "pd_UCell", "scz_UCell", 
              "sps_UCell")
mat_diseases <- scale(as.matrix(seurat_obj@meta.data[, dis_cols]))
colnames(mat_diseases) <- toupper(gsub("_UCell", "", dis_cols))

# 提取细胞类型
cell_cols <- c("Astrocyte", "Bergmann", "Choroid", "Endothelial_mural", "Endothelial_stalk", 
               "Ependymal", "Fibroblast", "Golgi", "Granule", "MLI1", "MLI2", 
               "Macrophage", "Microglia", "ODC", "OPC", "PLI", "Purkinje", "UBC")
mat_cells <- scale(as.matrix(seurat_obj@meta.data[, cell_cols]))

# ==========================================
# 2. 核心计算函数 (相关性 r 与 P 值)
# ==========================================
calc_cor_matrices <- function(mat_row, mat_col) {
  
  cor_mat <- cor(mat_row, mat_col, method = "pearson", use = "pairwise.complete.obs")
  
  p_mat <- matrix(1, nrow = ncol(mat_row), ncol = ncol(mat_col))
  
  for(i in 1:ncol(mat_row)) {
    for(j in 1:ncol(mat_col)) {
      try({
        test <- cor.test(mat_row[, i], mat_col[, j], method = "pearson")
        p_mat[i, j] <- test$p.value
      }, silent = TRUE)
    }
  }
  
  rownames(cor_mat) <- colnames(mat_row)
  colnames(cor_mat) <- colnames(mat_col)
  rownames(p_mat) <- colnames(mat_row)
  colnames(p_mat) <- colnames(mat_col)
  
  # 极端容错
  cor_mat[is.na(cor_mat)] <- 0
  p_mat[is.na(p_mat)] <- 1
  
  return(list(cor = cor_mat, p = p_mat))
}

cat("⏳ 正在计算 5 组皮尔逊相关性矩阵与 P 值...\n")

# A. Module (行) vs Cells (列)
res_mod_cell <- calc_cor_matrices(mat_modules, mat_cells)

# B. Module (行) vs Diseases (列)
res_mod_dis <- calc_cor_matrices(mat_modules, mat_diseases)

# C. Cells (行) vs Diseases (列)
res_cell_dis <- calc_cor_matrices(mat_cells, mat_diseases)

# D. Cells (行) vs Cells (列)
res_cell_cell <- calc_cor_matrices(mat_cells, mat_cells)

# E. Diseases (行) vs Diseases (列)
res_dis_dis <- calc_cor_matrices(mat_diseases, mat_diseases)

# ==========================================
# 3. 画图输出函数 (全局锁定 14pt 黑体不加粗)
# ==========================================
my_palette <- colorRampPalette(c("#4393C3", "#FFFFFF", "#D6604D"))(200)

plot_corr_heatmap <- function(res_list, pdf_name, w, h) {
  
  cairo_pdf(pdf_name, width = w, height = h, family = family_use)
  
  # 强制 14号字体 (ps=14), 不加粗 (font=1)
  par(ps = 14, family = family_use, font = 1, col.axis = "black")
  
  corrplot::corrplot(
    res_list$cor,
    method = "circle",
    order = "original",
    tl.pos = "lt",
    tl.col = "black",
    tl.cex = 1,
    tl.srt = 45,
    cl.pos = "r",
    cl.cex = 1,
    cl.offset = 0.5,
    col = my_palette,
    p.mat = res_list$p,
    sig.level = 0.05,
    insig = "blank",
    mar = c(1, 1, 1, 1)
  )
  
  dev.off()
  cat(" 💾 已保存:", pdf_name, "\n")
}

# ==========================================
# 4. 执行绘图
# ==========================================
cat("🎨 正在生成纯矢量 14pt 相关性热图...\n")

plot_corr_heatmap(res_mod_cell, "Figure_Cor_Heatmap_Modules_vs_Cells.pdf", w = 12, h = 5)
plot_corr_heatmap(res_mod_dis, "Figure_Cor_Heatmap_Modules_vs_Diseases.pdf", w = 11, h = 5)
plot_corr_heatmap(res_cell_dis, "Figure_Cor_Heatmap_Cells_vs_Diseases.pdf", w = 11, h = 10)
plot_corr_heatmap(res_cell_cell, "Figure_Cor_Heatmap_Cells_vs_Cells.pdf", w = 12, h = 12)
plot_corr_heatmap(res_dis_dis, "Figure_Cor_Heatmap_Diseases_vs_Diseases.pdf", w = 12, h = 12)

# ==========================================
# 5. 保存计算结果到 RData
# ==========================================
cat("📦 正在将 5 组相关性计算结果保存至 RData...\n")

save(res_mod_cell, res_mod_dis, res_cell_dis, res_cell_cell, res_dis_dis, 
     file = "Figure5_Correlation_Matrices_Complete.RData")

cat("🎉 大圆满！5 个高清热图及矩阵数据文件已全部生成。\n")









setwd('GSE290806_RAW/')
list.files()

#

# ==========================================
# 2. 合并数据集与基础 QC
# ==========================================
cat("🔄 正在合并所有样本...\n")

seurat_obj <- merge(seurat_list[[1]], y = seurat_list[-1], 
                    add.cell.ids = names(seurat_list))

rm(seurat_list, counts, sobj)
gc()

# 计算线粒体比例并过滤
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# 把 Age 转为 factor 保证画图顺序
seurat_obj$Age <- factor(seurat_obj$Age, levels = c("3", "12", "18", "24"))

# ==========================================
# 3. 标准化、降维与聚类
# ==========================================
cat("⏳ 正在进行标准化、PCA 与 UMAP...\n")

seurat_obj <- NormalizeData(seurat_obj, 
                            normalization.method = "LogNormalize", 
                            scale.factor = 10000)

seurat_obj <- FindVariableFeatures(seurat_obj, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)

seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

seurat_obj <- RunPCA(seurat_obj, 
                     features = VariableFeatures(object = seurat_obj), 
                     verbose = FALSE)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# ==========================================
# 4. 绘制 UMAP 与 Marker 验证图
# ==========================================
cat("🎨 正在绘制 UMAP 概览与 DotPlot...\n")

p_cluster <- DimPlot(seurat_obj, reduction = "umap", 
                     group.by = "seurat_clusters", 
                     label = TRUE, label.size = 14/.pt) +
  ggtitle("UMAP by Clusters") +
  custom_theme

p_age <- DimPlot(seurat_obj, reduction = "umap", 
                 group.by = "Age", label = FALSE) +
  ggtitle("UMAP by Age") +
  custom_theme

p_combined <- p_cluster + p_age + plot_layout(ncol = 2)

cairo_pdf("Figure_UMAP_PreAnnotation_14pt.pdf", 
          width = 12, height = 5.5, family = family_use)
print(p_combined)
dev.off()

# 准备文献 Marker 基因
marker_genes <- c(
  "Aqp4", "Mrc1", "Gdf10", "Ttr", "Pecam1", "Kcnj8", "Col1a1",
  "Lgi2", "Gabra6", "P2ry12", "Lypd6", "Mobp", "Pdgfra", "Gad2",
  "Ppp1r17", "Eomes"
)

genes_use <- marker_genes[marker_genes %in% rownames(seurat_obj)]

p_dot <- DotPlot(seurat_obj, features = genes_use, group.by = "seurat_clusters") +
  theme_classic() +
  theme(
    text = font_style_14pt,
    axis.text.x = element_text(size = 14, color = "black", face = "plain",
                               family = family_use, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, color = "black", face = "plain", family = family_use),
    axis.title = font_style_14pt,
    legend.text = font_style_14pt,
    legend.title = font_style_14pt,
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  labs(x = "Marker Genes", y = "Seurat Clusters")

cairo_pdf("Figure_Annotation_DotPlot_14pt.pdf", width = 12, height = 7, family = family_use)
print(p_dot)
dev.off()

cat("✅ 全部搞定！基础处理完毕\n")

library(Seurat)
library(ggplot2)

# ==========================================
# 0. 跨平台字体检测与全局 14pt 设置
# ==========================================
available_fonts <- names(pdfFonts())
family_use <- if ("Arial" %in% available_fonts) "Arial" else "sans"

font_style_14pt <- element_text(size = 14, color = "black", face = "plain", family = family_use)

# ==========================================
# 1. 完整细胞注释字典
# ==========================================
cluster_annotations <- c(
  "0"  = "Granule",
  "1"  = "Granule",
  "2"  = "Oligo",
  "3"  = "MLI",
  "4"  = "Oligo",
  "5"  = "Granule",
  "6"  = "Golgi",
  "7"  = "Granule",
  "8"  = "Bergmann",
  "9"  = "Astrocyte",
  "10" = "Oligo",
  "11" = "Astrocyte",
  "12" = "PLI",
  "13" = "Microglia",
  "14" = "OPC",
  "15" = "OPC",
  "16" = "Endo",
  "17" = "UBC",
  "18" = "Microglia",
  "19" = "Fibroblast",
  "20" = "Endo. Mural",
  "21" = "Purkinje",
  "22" = "PLI",
  "23" = "Choroid",
  "24" = "Golgi"
)

# 写入 metadata 新列
seurat_obj$cell_type <- unname(cluster_annotations[as.character(seurat_obj$seurat_clusters)])

cat("✅ 细胞注释映射完成！\n")

# ==========================================
# 2. 绘制并导出最终版注释 UMAP
# ==========================================
Idents(seurat_obj) <- "cell_type"

p_umap_final <- DimPlot(seurat_obj, reduction = "umap",
                        label = TRUE, label.size = 14/.pt, repel = TRUE) +
  theme_classic() +
  theme(
    text = font_style_14pt,
    axis.text = font_style_14pt,
    axis.title = font_style_14pt,
    legend.text = font_style_14pt,
    legend.title = element_blank()
  ) +
  ggtitle("UMAP by Cell Type")

cairo_pdf("Figure_UMAP_Annotated_14pt.pdf", width = 10, height = 7, family = family_use)
print(p_umap_final)
dev.off()

cat("🎉 最终注释版 UMAP 图已保存为 Figure_UMAP_Annotated_14pt.pdf！\n")

library(Seurat)
library(UCell)

# ==========================================
# 1. 规范化模块命名 (Module_2-5 -> M1-M4)
# ==========================================
# 这样算出来的列名会自动变成 M1_UCell, M2_UCell 等
names(filtered_signatures) <- paste0("M", 1:4)

# ==========================================
# 2. 96核全速计算 UCell 评分
# ==========================================
cat("🚀 正在开启 96 核全速计算 UCell 模块分数...\n")

seurat_obj <- AddModuleScore_UCell(
  obj = seurat_obj,
  features = filtered_signatures,
  ncores = 96,
  name = "_UCell"  # 后缀，结果列名为 M1_UCell
)

cat("✅ UCell 极速计算完成！\n")



setwd('E:/deskup/毕业论文/mouse/')
load('month_seurat.Rdata')


library(Seurat)
library(dplyr)
library(tidyr)
library(Mfuzz)
library(ggplot2)
library(scales)

# ==========================================
# 1. 提取 Purkinje 细胞子集并锁定 Module 4 基因
# ==========================================
cat("⏳ 正在提取 Purkinje 细胞子集并构建表达矩阵...\n")
sub_purkinje <- subset(seurat_obj, cell_type == "Purkinje")

mod4_genes <- c(
  "Pcp4", "Car8", "Igsf5", "Calb1", "Pcp2", "Pvalb", "Gng13", "Itpr1", 
  "Nsg1", "Ywhah", "Slc1a3", "Ppp1r17", "Homer3", "Dner", "Inpp5a", 
  "Fam107a", "Atp1b1", "Prkcg", "Atp2a2", "Rgs8", "Icmt", "Gad1", 
  "Cck", "Id2", "Gpr37l1", "Gdf10", "Lhx1os", "Thy1", "Slc1a6", 
  "Itm2b", "Gabra1", "Ckb", "Sptbn2"
)

genes_use <- mod4_genes[mod4_genes %in% rownames(sub_purkinje)]

# ==========================================
# 2. 提取 DotPlot 表达矩阵并跑标准 Mfuzz 获取聚类
# ==========================================
dp <- DotPlot(sub_purkinje, features = genes_use, group.by = "Age", assay = "RNA")
dp_df <- dp$data

mat <- dp_df %>%
  select(features.plot, id, avg.exp) %>%
  pivot_wider(names_from = id, values_from = avg.exp) %>%
  as.data.frame()

rownames(mat) <- mat$features.plot
mat$features.plot <- NULL
mat <- as.matrix(mat)
mat <- mat[, c("3", "12", "18", "24"), drop = FALSE]

# 清理标准差为 0 的基因
gene_sd <- apply(mat, 1, sd, na.rm = TRUE)
bad_genes <- names(gene_sd)[gene_sd == 0 | !is.finite(gene_sd)]
if (length(bad_genes) > 0) {
  mat <- mat[!rownames(mat) %in% bad_genes, ]
}

eset <- new("ExpressionSet", exprs = mat)
eset <- standardise(eset)
m <- mestimate(eset)

set.seed(2026)
clust <- mfuzz(eset, c = 4, m = m)

# ==========================================
# 3. 整理用于 ggplot2 渲染的实心渐变数据
# ==========================================
zmat <- exprs(eset)
gene_cluster <- data.frame(Gene = rownames(zmat), Cluster = clust$cluster, stringsAsFactors = FALSE)

gene_cluster$Membership <- sapply(seq_len(nrow(gene_cluster)), function(i) {
  clust$membership[gene_cluster$Gene[i], gene_cluster$Cluster[i]]
})

plot_df <- as.data.frame(zmat) %>%
  mutate(Gene = rownames(.)) %>%
  pivot_longer(cols = -Gene, names_to = "Age", values_to = "Expression") %>%
  left_join(gene_cluster, by = "Gene") %>%
  mutate(
    Age = as.numeric(Age),
    Cluster = factor(Cluster, levels = 1:4, labels = paste0("Cluster ", 1:4))
  ) %>%
  arrange(Membership) # 让高 Membership（核心基因）后画，压在最上层

# 计算每个 Cluster 的加权中心拟合线 (Centroid)
centroid_df <- lapply(1:4, function(k) {
  g_k <- gene_cluster$Gene[gene_cluster$Cluster == k]
  if (length(g_k) == 0) return(NULL)
  m_k <- clust$membership[g_k, k]
  mat_k <- zmat[g_k, , drop = FALSE]
  centroid <- apply(mat_k, 2, function(x) weighted.mean(x, m_k))
  data.frame(Age = as.numeric(colnames(mat_k)), Expression = centroid, Cluster = paste0("Cluster ", k))
}) %>% bind_rows() %>%
  mutate(Cluster = factor(Cluster, levels = paste0("Cluster ", 1:4)))

# ==========================================
# 4. 绘制高颜值、实心填满且拟合线呈红色的趋势图
# ==========================================
p <- ggplot() +
  # 垂直时间虚线
  geom_vline(xintercept = c(12, 18, 24), linetype = "dashed", linewidth = 0.4, color = "grey50") +
  
  # 密集基因轨迹：通过极高透明度叠加，自动形成类似实心热图的区块感
  geom_line(data = plot_df, aes(x = Age, y = Expression, group = Gene, color = Membership), linewidth = 0.6, alpha = 0.7) +
  
  # 中心加权拟合线：用醒目的红色高亮呈现
  geom_line(data = centroid_df, aes(x = Age, y = Expression, group = Cluster), color = "#D73027", linewidth = 1.3) +
  
  # 四个 Cluster 纵向分面
  facet_wrap(~ Cluster, ncol = 1, scales = "fixed") +
  
  # 自定义渐变色：边缘淡黄/亮黄 -> 中间高 membership 呈现炽热红
  scale_color_gradientn(
    colours = c("#FFFF99", "#FFCC00", "#FF6600", "#D73027", "#800026"),
    limits = c(0, 1),
    oob = scales::squish
  ) +
  
  scale_x_continuous(breaks = c(3, 12, 18, 24), limits = c(3, 24), expand = c(0, 0)) +
  labs(x = "Time (Month)", y = "Expression change") +
  
  theme_bw(base_size = 14, base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "#FFFDF0", color = "black", linewidth = 0.6), # 经典米黄色/浅底
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.title = element_text(size = 14, color = "black", family = "Arial"),
    axis.text = element_text(size = 14, color = "black", family = "Arial"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0.04, "lines"),
    legend.position = "none"
  )

# 在每个分面内部加上 Cluster 标签
label_df <- data.frame(Cluster = factor(paste0("Cluster ", 1:4), levels = paste0("Cluster ", 1:4)), Age = 14.5, Expression = 0)
p <- p + geom_text(data = label_df, aes(x = Age, y = Expression, label = Cluster), size = 4, family = "Arial", color = "black", fontface = "bold")

# ==========================================
# 5. 导出矢量 PDF
# ==========================================
ggsave(filename = "Purkinje_Module4_Mfuzz_PerfectFilled.pdf", plot = p, width = 3.5, height = 8.5, device = cairo_pdf)

cat("🎉 完美搞定！完全绕过内置报错、呈现满分实心渐变与红色中心线的 Mfuzz 图已保存为 Purkinje_Module4_Mfuzz_PerfectFilled.pdf\n")




library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

# ==========================================
# 0. 跨平台字体检测与全局 14pt 设置
# ==========================================
available_fonts <- names(pdfFonts())
family_use <- if ("Arial" %in% available_fonts) "Arial" else "sans"

# ==========================================
# 1. 提取 Purkinje 细胞子集
# ==========================================
sub_purkinje <- subset(seurat_obj, cell_type == "Purkinje")

# ==========================================
# 2. 根据 Mfuzz 结果定义基因的 Cluster 排序
# ==========================================
cluster_mapping <- data.frame(
  Gene = c(
    "Slc1a3", "Ppp1r17", "Homer3", "Fam107a", "Atp1b1", "Icmt", "Gad1", 
    "Gpr37l1", "Gdf10", "Thy1",                              # Cluster 1 (10个)
    "Itpr1", "Dner", "Inpp5a", "Prkcg", "Sptbn2",           # Cluster 2 (5个)
    "Calb1", "Atp2a2", "Rgs8", "Lhx1os", "Slc1a6",          # Cluster 3 (5个)
    "Pcp4", "Car8", "Pcp2", "Pvalb", "Gng13", "Nsg1", 
    "Ywhah", "Cck", "Id2", "Itm2b", "Gabra1", "Ckb"         # Cluster 4 (12个)
  ),
  Cluster = c(
    rep("Cluster 1", 10),
    rep("Cluster 2", 5),
    rep("Cluster 3", 5),
    rep("Cluster 4", 12)
  )
)

# 确保只挑选当前子集中存在的基因
genes_ordered <- cluster_mapping$Gene[cluster_mapping$Gene %in% rownames(sub_purkinje)]

# ==========================================
# 3. 提取 DotPlot 数据并应用新顺序
# ==========================================
dp_res <- DotPlot(sub_purkinje, features = genes_ordered, group.by = "Age", assay = "RNA")
plot_df <- dp_res$data

# 强制固定 X 轴基因顺序
plot_df$features.plot <- factor(plot_df$features.plot, levels = genes_ordered)

# Y 轴年龄从上到下是 24 -> 3 个月
plot_df$id <- factor(plot_df$id, levels = rev(c("3", "12", "18", "24")))

# 计算每个 Cluster 边界的位置
cluster_counts <- table(cluster_mapping$Cluster[cluster_mapping$Gene %in% genes_ordered])
v_lines <- cumsum(cluster_counts)[-length(cluster_counts)] + 0.5

# ==========================================
# 4. 绘制气泡图
# ==========================================
p <- ggplot(plot_df, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  # Cluster 垂直分隔线
  geom_vline(xintercept = v_lines, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  # 红蓝配色
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, 0, 1),
    oob = scales::squish,
    name = "Average Expression"
  ) +
  scale_size(
    range = c(0, 6),
    limits = c(0, 100),
    breaks = c(20, 50, 80),
    name = "Percent Expressed"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 14, color = "black", face = "plain", family = family_use),
    axis.text.x = element_text(size = 14, color = "black", face = "plain", 
                               family = family_use, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, color = "black", face = "plain", family = family_use),
    axis.title = element_blank(),
    legend.text = element_text(size = 14, color = "black", face = "plain", family = family_use),
    legend.title = element_text(size = 14, color = "black", face = "plain", family = family_use),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "right"
  )

# ==========================================
# 5. 导出 14pt 矢量 PDF
# ==========================================
cairo_pdf("DotPlot_Purkinje_Module4_MfuzzOrdered_14pt.pdf", 
          width = 13, height = 5.5, family = family_use)
print(p)
dev.off()

cat("🎉 气泡图已保存为 DotPlot_Purkinje_Module4_MfuzzOrdered_14pt.pdf！\n")





