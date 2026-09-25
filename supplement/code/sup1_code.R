#supplement1
# ==========================================
# 0. 环境准备与数据加载
# ==========================================
library(Seurat)
library(ggplot2)
library(dplyr)

setwd('E:/deskup/毕业论文/bechmark/raodong/data/')
message("正在加载 pbmc.RDS ...")
seurat_obj <- readRDS('pbmc.RDS')
library(Seurat)
library(ggplot2)

# 强制注册 Arial 字体
windowsFonts(Arial = windowsFont("Arial"))

# ==========================================
# 1. 规范化细胞群体顺序与提取 Marker
# ==========================================
cell_order <- c("B_Cell", "T_NK_Cell", "Macrophage")
seurat_obj$Broad_CellType <- factor(seurat_obj$Broad_CellType, levels = cell_order)
Idents(seurat_obj) <- "Broad_CellType"

classic_markers <- c(
  # B cells
  'CD19', 'MS4A1', 'CD79A', 'CD79B', 'SDC1', 'CD38', 'JCHAIN', 'MZB1',
  # T/NK cells
  'CD3D', 'CD3E', 'CD8A', 'CD4', 'NCAM1', 'NKG7', 'GNLY', 'GZMB', 'PRF1',
  # Macrophages
  'CD14', 'FCGR3A', 'S100A8', 'S100A9', 'ITGAM', 'ITGAX', 'CD68', 'C1QA'
)

valid_markers <- intersect(classic_markers, rownames(seurat_obj))

# ==========================================
# 2. 绘制 Nature 风格 + 7pt 铁律的 DotPlot
# ==========================================
message("正在渲染 Nature 风格 Marker DotPlot ...")

font_style_7pt <- element_text(size = 7, color = "black", face = "plain", family = "Arial")

p_dot <- DotPlot(seurat_obj, features = valid_markers, dot.scale = 4) +
  scale_color_gradientn(colors = c("#f7f7f7", "#f8c8cd", "#ed8590", "#d71345")) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    text = font_style_7pt,
    axis.text.x = element_text(size = 7, color = "black", angle = 45, 
                               hjust = 1, vjust = 1, face = "plain", family = "Arial"),
    axis.text.y = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
    axis.title = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.text = font_style_7pt,
    legend.title = font_style_7pt,
    legend.key.size = unit(0.3, "cm"),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

# ==========================================
# 3. 导出高质量 PDF
# ==========================================
cairo_pdf("Fig_CellType_Marker_DotPlot_Nature_7pt.pdf", 
          width = 8, height = 3, family = "Arial")
print(p_dot)
dev.off()

message("🎉 Nature 风格 Marker 气泡图 (7pt) 已成功保存！")

# ==========================================
# 0. 加载必备的包
# ==========================================
setwd('E:/deskup/毕业论文/bechmark/raodong/data/')
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

load('Benchmark1_UCell_Scores.Rdata')
pbmc<-seurat_obj

# ==========================================
# 1. 提取目标数据并执行 0-1 Min-Max 缩放
# ==========================================
mescnet_score_raw <- mescnet_scores[["resolution_4_number_4"]]

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

message("正在执行 Min-Max 0-1 归一化...")
mescnet_scaled <- min_max_scale(mescnet_score_raw)
spectra_scaled <- min_max_scale(spectra_scores)
expimap_scaled <- min_max_scale(expimap_scores)

# ==========================================
# 2. 将缩放后的得分注入 Seurat Metadata
# ==========================================
rownames(mescnet_scaled) <- colnames(pbmc)
rownames(spectra_scaled) <- colnames(pbmc)
rownames(expimap_scaled) <- colnames(pbmc)

pbmc$MEscnet_LPS <- mescnet_scaled$MEscnet_Mod_38_UCell
pbmc$MEscnet_PMA <- mescnet_scaled$MEscnet_Mod_26_UCell
pbmc$MEscnet_IFN <- mescnet_scaled$MEscnet_Mod_24_UCell

pbmc$SPECTRA_LPS <- spectra_scaled$Spectra_Factor_172_UCell
pbmc$SPECTRA_PMA <- spectra_scaled$Spectra_Factor_183_UCell
pbmc$SPECTRA_IFN <- spectra_scaled$Spectra_Factor_155_UCell

pbmc$expiMap_LPS <- expimap_scaled$all_TLR_signaling_UCell
pbmc$expiMap_PMA <- expimap_scaled$T_tcr.activation_UCell
pbmc$expiMap_IFN <- expimap_scaled$all_type.I.ifn.response_UCell

# ==========================================
# 3. 定义统一的颜色条（固定范围 0-1）
# ==========================================
# 所有图共享同一个颜色条标准：0 到 1
common_colors <- c("#e0e0e0", "#ffdf91", "#e34a33", "#a80016")

# ==========================================
# 4. 绘制单个 FeaturePlot 并保存为 PDF（带图例）
# ==========================================
plot_and_save <- function(feature_col, title_text, filename) {
  
  p <- FeaturePlot(pbmc, features = feature_col, pt.size = 0.1, order = TRUE) +
    scale_color_gradientn(
      colors = common_colors,
      limits = c(0, 1),           # ✅ 固定范围 0-1
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "Module Score"
    ) +
    labs(title = title_text) +
    theme_void() +
    theme(
      text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
      plot.title = element_text(size = 16, color = "black", face = "bold", family = "Arial", hjust = 0.5),
      legend.text = element_text(size = 12, color = "black", family = "Arial"),
      legend.title = element_text(size = 12, color = "black", family = "Arial"),
      legend.position = "right",
      legend.key.width = unit(0.8, "cm"),
      legend.key.height = unit(1.5, "cm")
    )
  
  # 保存为 PDF
  ggsave(filename, plot = p, width = 7, height = 5, dpi = 600)
  message(paste("✅ 已保存:", filename))
}

# ==========================================
# 5. 批量生成 9 个 PDF
# ==========================================

# 定义 9 个图的信息
plot_info <- list(
  # MEscnet (第一行)
  list(feature = "MEscnet_LPS", title = "MEscnet: LPS", file = "FIG_MEscnet_LPS.pdf"),
  list(feature = "MEscnet_PMA", title = "MEscnet: TCR", file = "FIG_MEscnet_TCR.pdf"),
  list(feature = "MEscnet_IFN", title = "MEscnet: IFNγ", file = "FIG_MEscnet_IFNγ.pdf"),
  
  # SPECTRA (第二行)
  list(feature = "SPECTRA_LPS", title = "SPECTRA: LPS", file = "FIG_SPECTRA_LPS.pdf"),
  list(feature = "SPECTRA_PMA", title = "SPECTRA: TCR", file = "FIG_SPECTRA_TCR.pdf"),
  list(feature = "SPECTRA_IFN", title = "SPECTRA: IFNγ", file = "FIG_SPECTRA_IFNγ.pdf"),
  
  # expiMap (第三行)
  list(feature = "expiMap_LPS", title = "expiMap: LPS", file = "FIG_expiMap_LPS.pdf"),
  list(feature = "expiMap_PMA", title = "expiMap: TCR", file = "FIG_expiMap_TCR.pdf"),
  list(feature = "expiMap_IFN", title = "expiMap: IFNγ", file = "FIG_expiMap_IFNγ.pdf")
)

# 循环生成
for (info in plot_info) {
  plot_and_save(info$feature, info$title, info$file)
}

message("🎉 全部完成！共生成 9 个 PDF 文件，所有图共享相同的颜色条标准 (0-1)")

# ==========================================
# 单独生成颜色条图例 PDF
# ==========================================












library(ggplot2)
library(cowplot) # 专门用来提取图例的神器

setwd('E:/deskup/毕业论文/bechmark/raodong/data/')

# 1. 建立一个极简的虚拟数据集（只要包含 0 和 1 就能撑起完整的 Color Bar）
dummy_data <- data.frame(x = c(1, 2), y = c(1, 2), score = c(0, 1))
common_colors <- c("#e0e0e0", "#ffdf91", "#e34a33", "#a80016")

# 2. 画一个带有完美图例的虚拟图
p_dummy <- ggplot(dummy_data, aes(x = x, y = y, color = score)) +
  geom_point() +
  scale_color_gradientn(
    colors = common_colors,
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    name = "Scaled\nScore",
    # 精细控制颜色条的宽高和标题位置，方便放进 Inkscape
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(1, "cm"),
      barheight = unit(6, "cm"),
      ticks.linewidth = 1,
      frame.colour = "black" # 给颜色条加个纯黑外框更清晰
    )
  ) +
  theme_void() +
  theme(
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, color = "black")
  )

# 3. 剥离图例
legend_only <- get_legend(p_dummy)

# 4. 单独保存为一个极其干净的 PDF
ggsave("Shared_ColorBar_Legend.pdf", plot = legend_only, width = 2, height = 4)

message("🎉 搞定！单独的颜色条已保存为 Shared_ColorBar_Legend.pdf，可以直接拖入 Inkscape。")