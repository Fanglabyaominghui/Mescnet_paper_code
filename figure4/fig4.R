library(Seurat)
setwd('E:/deskup/fig2.0/fig5/data/')
load('T53_anno.Rdata')

library(Seurat)
library(schard)
library(SPARK)
library(tools)

# ==========================================
# 1. 设定目标文件路径
# ==========================================
h5ad_dir <- "E:/deskup/SP/data/"
target_file <- file.path(h5ad_dir, "Mouse2_T353.h5ad")

if (!file.exists(target_file)) {
  stop("找不到文件: ", target_file, "，请检查文件夹里的文件名！")
}

cat("\n========================================\n")
message("开始单独处理切片: Mouse2_T353")
message("正在读取并转换 h5ad: ", target_file)

# ==========================================
# 2. 读取数据与提取坐标
# ==========================================
obj <- schard::h5ad2seurat(target_file)

message("提取 Count 矩阵与空间坐标...")
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
locs <- obj@reductions$Xspatial_@cell.embeddings[, 1:2]

# ==========================================
# 3. 运行 SPARK-X 进行空间高可变基因分析
# ==========================================
message("运行 SPARK-X...")
sparkx_out <- sparkx(counts, locs, numCores = 1, option = "mixture")
res <- sparkx_out$res_mtest

# ==========================================
# 4. 数据清洗：剔除线粒体基因
# ==========================================
message("在全量结果中剔除线粒体基因...")
mt_idx <- grep("^MT-|^mt-", rownames(res))

if (length(mt_idx) > 0) {
  res_clean <- res[-mt_idx, ]
  message("成功拦截并剔除了 ", length(mt_idx), " 个线粒体基因。")
} else {
  res_clean <- res
}

# ==========================================
# 5. 排序并提取 Top 2000 核心基因
# ==========================================
message("提取 Top 2000 高空间可变基因...")
res_sorted <- res_clean[order(res_clean$adjustedPval), ]
top_n <- min(2000, nrow(res_sorted))
t353_svg_genes <- rownames(res_sorted)[1:top_n]

# ==========================================
# 6. 打印结果
# ==========================================
cat("\n🎉 Mouse2_T353 切片处理完毕！\n")
cat("共提取高空间可变基因数量:", length(t353_svg_genes), "\n")
cat("前 20 个核心基因预览:\n")
print(head(t353_svg_genes, 20))

# 1. 按照 adjustedPval 从小到大（显著性从高到低）进行排序
res_sorted <- res_clean[order(res_clean$adjustedPval), ]

# 2. 提取排在最前面的 500 个基因名称
top500_genes <- rownames(res_sorted)[1:500]

# 3. 打印前 20 个看看效果
print(head(top500_genes, 20))






core_svg_genes<-top500_genes
# hdwgcna构建元细胞 (单样本 T53 平铺版)
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(schard)
library(Mescnet)
enableWGCNAThreads(nThreads = 8)
theme_set(theme_cowplot())


if (!file.exists(target_file)) {
  stop("找不到文件: ", target_file, "，请检查路径和文件名！")
}

cat("\n========================================\n")
message("开始单独处理切片: Mouse2_T53")

# ==========================================
# 2. 读取数据并锁定默认 Assay
# ==========================================
message("正在读取 h5ad 文件...")
snhx <- schard::h5ad2seurat(target_file)
DefaultAssay(snhx) <- "RNA"

# ==========================================
# 3. 初始化 hdWGCNA
# ==========================================
message("正在初始化 hdWGCNA...")
snhx <- SetupForWGCNA(
  snhx,
  gene_select = "custom", 
  gene_list = core_svg_genes,  # ⚠️如果你刚才存的是 top500_genes，请换成那个变量名
  wgcna_name = "MEscnet"
)

# ==========================================
# 4. 华大平台专属空间聚合 (构建 Metacells)
# ==========================================
message("正在构建空间元点 (Metacells)...")
snhx <- MetacellsByGroups(
  seurat_obj = snhx,
  group.by = "annotation", 
  ident.group = "annotation",
  reduction = "Xspatial_", 
  k = 25, 
  max_shared = 10
)

# ==========================================
# 5. 标准化空间元点
# ==========================================
message("标准化 Metacells...")
snhx <- NormalizeMetacells(snhx)

# ==========================================
# 6. 获取合法的 annotation (防止 NULL 报错)
# ==========================================
valid_annotations <- as.character(unique(snhx$annotation))
valid_annotations <- valid_annotations[!is.na(valid_annotations)]

# ==========================================
# 7. 提取 datExpr 表达矩阵
# ==========================================
message("正在提取 datExpr 矩阵...")
snhx <- SetDatExpr(
  snhx,
  group_name = valid_annotations, 
  group.by = "annotation",
  assay = "RNA",                  
  slot = "data"                   
)

# 获取最终矩阵
datExpr_T53 <- GetDatExpr(snhx)
cat("🎉 成功提取！当前 T53 矩阵维度为：", dim(datExpr_T53)[1], "个空间元点 x", dim(datExpr_T53)[2], "个基因\n")


base_obj<-snhx
base_obj<-pick_power_scale_free_power(base_obj)
# 算网络模块
base_obj <- ComputeMEscnetModules(base_obj, min_genes_per_module = 15, resolution = 1, number = 3)
modules_ids <- base_obj@misc$MEscnet$MEscnet_modules$module_ids
module_genes <- base_obj@misc$MEscnet$MEscnet_modules$gene_lists


setwd('E:/deskup/毕业论文/mouse/')
save(base_obj,file='base_obj.Rdata')





#figA
library(ggplot2)
library(dplyr)
library(ggrastr)

# 强制注册 Arial 字体
windowsFonts(Arial = windowsFont("Arial"))

cat("🎨 正在为您绘制带栅格化的空间图 (严格 7pt Arial)...\n")

# ==========================================
# 1. 提取坐标和注释
# ==========================================
plot_data <- data.frame(
  X = base_obj@reductions$Xspatial_@cell.embeddings[, 1],
  Y = base_obj@reductions$Xspatial_@cell.embeddings[, 2],
  annotation = base_obj@meta.data$annotation
)

# ==========================================
# 2. 强制图层绘制顺序 (Z-index)
# ==========================================
plot_data$annotation <- factor(plot_data$annotation, 
                               levels = c("molecular layer", "granular layer", "white matter", "purkinje layer"))

plot_data <- plot_data %>% arrange(annotation)

# ==========================================
# 3. HSV 统一色彩体系 (S=82, V=82)
# ==========================================
my_colors_illustration <- c(
  "molecular layer" = "#25D185",  
  "granular layer"  = "#ADD125",  
  "white matter"    = "#D12C25",  
  "purkinje layer"  = "#C525D1"   
)

# ==========================================
# 4. 绘图与输出 (栅格化点图层 + 严格 7pt 字体)
# ==========================================
font_style_7pt <- element_text(size = 7, color = "black", face = "plain", family = "Arial")

p_spatial <- ggplot(plot_data, aes(X, Y, color = annotation)) +
  # 🌟 核心修改：利用 rasterise 将点图层以 600 DPI 栅格化，防止进 Inkscape 卡死
  rasterise(geom_point(size = 0.6, alpha = 1, stroke = 0), dpi = 600) + 
  scale_color_manual(values = my_colors_illustration) +
  coord_fixed() +
  theme_void() +
  theme(
    text = font_style_7pt,
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = font_style_7pt,
    legend.key.size = unit(0.3, "cm"),
    plot.title = element_text(hjust = 0.5, size = 7, color = "black", face = "plain", family = "Arial", margin = margin(b = 5)),
    plot.margin = margin(2, 2, 2, 2),
    panel.background = element_rect(fill = "white", color = NA) 
  ) +
  # 缩小图例点大小以完美匹配 7pt 字体
  guides(color = guide_legend(override.aes = list(size = 2))) +
  labs(title = "Cerebellum Spatial Architecture")

# 🌟 绝杀修复：换用 cairo_pdf 引擎，完美识别系统 Arial 字体
cairo_pdf('figA_Cerebellum_IllustrationStyle_Raster_7pt.pdf', height = 3.5, width = 4.5, family = "Arial")
print(p_spatial)
dev.off()

cat("🎉 搞定！带有 600 DPI 栅格化图层且全部重置为 7pt 字体的高清 PDF 已生成！\n")





#fig2

datExpr <- seurat_obj@misc$MEscnet$datExpr
seurat_obj<-base_obj
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
  ncores = 24,
  name = "_UCell" # UCell 会自动在模块名后缀加上 "_UCell" 以防与原基因名冲突
)

cat("🎉 UCell 评分计算完毕！结果已无缝添加至 seurat_obj@meta.data 中。\n")

# ==============================================================================
# 3. 检查生成的 UCell 列名
# ==============================================================================
# 查看刚刚加进去的评分列，它们通常看起来像 "Module_1_UCell", "Module_2_UCell" 等
ucell_cols <- grep("_UCell$", colnames(seurat_obj@meta.data), value = TRUE)
print(head(ucell_cols))






#fig2
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggnewscale)
library(scales)
library(cowplot)

# ==========================================
# 0. 跨平台字体检测 (7pt Arial)
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

font_style_7pt <- element_text(size = 7, color = "black", face = "plain", family = family_use)

# ==========================================
# 1. 提取坐标与 Module 评分数据
# ==========================================
plot_df <- data.frame(
  Dim_1 = seurat_obj@reductions$Xspatial_@cell.embeddings[, 1],
  Dim_2 = seurat_obj@reductions$Xspatial_@cell.embeddings[, 2],
  MOD1_UCell = seurat_obj@meta.data$Module_2_UCell,
  MOD2_UCell = seurat_obj@meta.data$Module_3_UCell,
  MOD3_UCell = seurat_obj@meta.data$Module_4_UCell,
  MOD4_UCell = seurat_obj@meta.data$Module_5_UCell
)

# ==========================================
# 2. 颜色设定
# ==========================================
col_mod1 <- alpha("#aa363d", 0.9)   # 红色
col_mod2 <- alpha("#228fbd", 0.9)   # 紫色
col_mod3 <- alpha("#9b95c9", 0.95)  # 绿色 (置顶层)
col_mod4 <- alpha("#769149", 0.9)   # 蓝色 (打底层)

# ==========================================
# 3. 阈值设定 (q70 - q95)
# ==========================================
q_low <- 0.70
q_high <- 0.95

q30_mod1 <- quantile(plot_df$MOD1_UCell, q_low, na.rm = TRUE)
q95_mod1 <- quantile(plot_df$MOD1_UCell, q_high, na.rm = TRUE)

q30_mod2 <- quantile(plot_df$MOD2_UCell, q_low, na.rm = TRUE)
q95_mod2 <- quantile(plot_df$MOD2_UCell, q_high, na.rm = TRUE)

q30_mod3 <- quantile(plot_df$MOD3_UCell, q_low, na.rm = TRUE)
q95_mod3 <- quantile(plot_df$MOD3_UCell, q_high, na.rm = TRUE)

q30_mod4 <- quantile(plot_df$MOD4_UCell, q_low, na.rm = TRUE)
q95_mod4 <- quantile(plot_df$MOD4_UCell, q_high, na.rm = TRUE)

# ==========================================
# 4. 图例参数配置 (横向、强制3个数字绝对对齐)
# ==========================================
my_guide <- guide_colorbar(
  barwidth = unit(2.5, "cm"),
  barheight = unit(0.25, "cm"),
  title.position = "right",
  title.vjust = 0.5,
  label.position = "bottom",
  direction = "horizontal"
)

exact_3_breaks <- function(x) c(x[1], (x[1] + x[2]) / 2, x[2])
fmt_mod <- function(x) sprintf("%.2f", x)

# ==========================================
# 5. 绘制完整基础图 (Z轴：4 -> 1 -> 2 -> 3)
# ==========================================
p_base <- ggplot(plot_df, aes(x = Dim_1, y = Dim_2)) +
  # 灰色背景点
  geom_point(color = "#E8E8E8", size = 0.2, stroke = 0) +
  # ---------- 最底层：Layer 1 (MOD4 蓝色) ----------
geom_point(data = plot_df %>% filter(MOD4_UCell > q30_mod4) %>% arrange(MOD4_UCell), 
           aes(color = MOD4_UCell), size = 0.3, stroke = 0) +
  scale_color_gradient(low = "transparent", high = col_mod4, name = "4",
                       limits = c(q30_mod4, q95_mod4), oob = squish, 
                       guide = my_guide, breaks = exact_3_breaks, labels = fmt_mod) +
  new_scale_color() +
  # ---------- 次底层：Layer 2 (MOD1 红色) ----------
geom_point(data = plot_df %>% filter(MOD1_UCell > q30_mod1) %>% arrange(MOD1_UCell), 
           aes(color = MOD1_UCell), size = 0.3, stroke = 0) +
  scale_color_gradient(low = "transparent", high = col_mod1, name = "1",
                       limits = c(q30_mod1, q95_mod1), oob = squish, 
                       guide = my_guide, breaks = exact_3_breaks, labels = fmt_mod) +
  new_scale_color() +
  # ---------- 次顶层：Layer 3 (MOD2 紫色) ----------
geom_point(data = plot_df %>% filter(MOD2_UCell > q30_mod2) %>% arrange(MOD2_UCell), 
           aes(color = MOD2_UCell), size = 0.3, stroke = 0) +
  scale_color_gradient(low = "transparent", high = col_mod2, name = "2",
                       limits = c(q30_mod2, q95_mod2), oob = squish, 
                       guide = my_guide, breaks = exact_3_breaks, labels = fmt_mod) +
  new_scale_color() +
  # ---------- 最顶层：Layer 4 (MOD3 绿色) ----------
geom_point(data = plot_df %>% filter(MOD3_UCell > q30_mod3) %>% arrange(MOD3_UCell), 
           aes(color = MOD3_UCell), size = 0.3, stroke = 0) +
  scale_color_gradient(low = "transparent", high = col_mod3, name = "3",
                       limits = c(q30_mod3, q95_mod3), oob = squish, 
                       guide = my_guide, breaks = exact_3_breaks, labels = fmt_mod) +
  coord_fixed() +
  theme_void() +
  theme(
    text = font_style_7pt,
    plot.title = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = font_style_7pt,
    legend.title = font_style_7pt,
    legend.margin = margin(l = 5),
    legend.spacing.y = unit(0.2, "cm"),
    plot.margin = margin(5, 5, 5, 5),
    panel.background = element_rect(fill = "white", color = NA)
  )

# ==========================================
# 6. 分离与导出
# ==========================================
p_main_only <- p_base + theme(legend.position = "none")
p_legend_only <- get_legend(p_base)

cairo_pdf("Figure_5e_Modules_Overlay_AutoBreaks_MainPlot.pdf", 
          height = 10, width = 10, family = family_use)
print(p_main_only)
dev.off()

cairo_pdf("Figure_5e_Modules_Overlay_AutoBreaks_Legend.pdf", 
          height = 10, width = 2, family = family_use)
print(ggdraw(p_legend_only))
dev.off()

cat("✅ Module 图层叠加图已生成！\n")












library(Seurat)
library(ggplot2)
library(dplyr)
library(ggnewscale)
library(scales)
library(cowplot)

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

# 强制全局 14号、纯黑、不加粗
font_style_14pt <- element_text(size = 14, color = "black", face = "plain", family = family_use)

cat("⏳ 正在进行标准化 (NormalizeData) 并绘制纯矢量基因叠加图...\n")

# ==========================================
# 1. 先显式标准化，确保 data 层数据正确
# ==========================================
seurat_obj <- NormalizeData(seurat_obj, 
                            normalization.method = "LogNormalize", 
                            scale.factor = 10000, 
                            verbose = FALSE)

# 提取坐标与标准化后的基因表达量
genes_expr <- as.data.frame(t(as.matrix(GetAssayData(seurat_obj, layer = "data")[c("Plp1", "Pcp4", "Snap25", "Ptgds"), ])))

plot_df <- data.frame(
  Dim_1 = seurat_obj@reductions$Xspatial_@cell.embeddings[, 1],
  Dim_2 = seurat_obj@reductions$Xspatial_@cell.embeddings[, 2],
  Plp1 = genes_expr$Plp1,
  Pcp4 = genes_expr$Pcp4,
  Snap25 = genes_expr$Snap25,
  Ptgds = genes_expr$Ptgds
)

# ==========================================
# 2. 颜色设定
# ==========================================
col_g1 <- alpha("#aa363d", 0.9)   # 红色 (Plp1)
col_g2 <- alpha("#9b95c9", 0.9)   # 紫色 (Pcp4)
col_g3 <- alpha("#769149", 0.95)  # 绿色 (Snap25)
col_g4 <- alpha("#228fbd", 0.9)   # 蓝色 (Ptgds)

# ==========================================
# 3. 阈值设定 (q5 - q95)
# ==========================================
q_low <- 0.70
q_high <- 0.99

q30_g1 <- quantile(plot_df$Plp1, q_low, na.rm = TRUE)
q95_g1 <- quantile(plot_df$Plp1, q_high, na.rm = TRUE)

q30_g2 <- quantile(plot_df$Pcp4, q_low, na.rm = TRUE)
q95_g2 <- quantile(plot_df$Pcp4, q_high, na.rm = TRUE)

q30_g3 <- quantile(plot_df$Snap25, q_low, na.rm = TRUE)
q95_g3 <- quantile(plot_df$Snap25, q_high, na.rm = TRUE)

q30_g4 <- quantile(plot_df$Ptgds, q_low, na.rm = TRUE)
q95_g4 <- quantile(plot_df$Ptgds, q_high, na.rm = TRUE)

# ==========================================
# 4. 图例参数配置 (横向、14pt 字体、强制3个数字对齐)
# ==========================================
my_guide <- guide_colorbar(
  barwidth = unit(3.5, "cm"),
  barheight = unit(0.35, "cm"),
  title.position = "right",
  title.vjust = 0.5,
  label.position = "bottom",
  direction = "horizontal"
)

exact_3_breaks <- function(x) c(x[1], (x[1] + x[2]) / 2, x[2])
fmt_gene <- function(x) sprintf("%.1f", x)

# ==========================================
# 5. 绘制完整基础图 (纯矢量 geom_point，Z轴: 4 -> 1 -> 2 -> 3)
# ==========================================
p_base <- ggplot(plot_df, aes(x = Dim_1, y = Dim_2)) +
  # 灰色背景点
  geom_point(color = "#E8E8E8", size = 0.2, stroke = 0) +
  # ---------- 最底层：Layer 1 (Ptgds 蓝色) ----------
geom_point(data = plot_df %>% filter(Ptgds > q30_g4) %>% arrange(Ptgds), 
           aes(color = Ptgds), size = 0.3, stroke = 0) +
  scale_color_gradient(low = "transparent", high = col_g4, name = "Ptgds",
                       limits = c(q30_g4, q95_g4), oob = squish, 
                       guide = my_guide, breaks = exact_3_breaks, labels = fmt_gene) +
  new_scale_color() +
  # ---------- 次底层：Layer 2 (Plp1 红色) ----------
geom_point(data = plot_df %>% filter(Plp1 > q30_g1) %>% arrange(Plp1), 
           aes(color = Plp1), size = 0.3, stroke = 0) +
  scale_color_gradient(low = "transparent", high = col_g1, name = "Plp1",
                       limits = c(q30_g1, q95_g1), oob = squish, 
                       guide = my_guide, breaks = exact_3_breaks, labels = fmt_gene) +
  new_scale_color() +
  # ---------- 次顶层：Layer 3 (Pcp4 紫色) ----------
geom_point(data = plot_df %>% filter(Pcp4 > q30_g2) %>% arrange(Pcp4), 
           aes(color = Pcp4), size = 0.3, stroke = 0) +
  scale_color_gradient(low = "transparent", high = col_g2, name = "Pcp4",
                       limits = c(q30_g2, q95_g2), oob = squish, 
                       guide = my_guide, breaks = exact_3_breaks, labels = fmt_gene) +
  new_scale_color() +
  # ---------- 最顶层：Layer 4 (Snap25 绿色) ----------
geom_point(data = plot_df %>% filter(Snap25 > q30_g3) %>% arrange(Snap25), 
           aes(color = Snap25), size = 0.3, stroke = 0) +
  scale_color_gradient(low = "transparent", high = col_g3, name = "Snap25",
                       limits = c(q30_g3, q95_g3), oob = squish, 
                       guide = my_guide, breaks = exact_3_breaks, labels = fmt_gene) +
  coord_fixed() +
  theme_void() +
  theme(
    text = font_style_14pt,
    plot.title = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = font_style_14pt,
    legend.title = font_style_14pt,
    legend.margin = margin(l = 5),
    legend.spacing.y = unit(0.3, "cm"),
    plot.margin = margin(10, 10, 10, 10),
    panel.background = element_rect(fill = "white", color = NA)
  )

# ==========================================
# 6. 分离与导出 (纯矢量 PDF)
# ==========================================
p_main_only <- p_base + theme(legend.position = "none")
p_legend_only <- get_legend(p_base)

cairo_pdf("Figure_5e_Genes_Overlay_MainPlot_Vector.pdf", 
          height = 10, width = 10, family = family_use)
print(p_main_only)
dev.off()

cairo_pdf("Figure_5e_Genes_Overlay_Legend_Vector.pdf", 
          height = 10, width = 2.5, family = family_use)
print(ggdraw(p_legend_only))
dev.off()

cat("✅ 标准化完成，纯矢量基因叠加图已成功导出！\n")



#fig3
suppressPackageStartupMessages({
  library(Seurat)
  library(schard)
  library(UCell)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(tidyr)
})

# ==========================================
# 0. 全局环境与字体设置 (7pt Arial)
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

font_style_7pt <- element_text(size = 7, color = "black", face = "plain", family = family_use)

# ==========================================
# 1. 重算 UCell 流水线 (使用 filtered_signatures)
# ==========================================
cat("🚀 正在为 filtered_signatures 命名并重新计算 UCell...\n")

# 强制命名列表，确保 UCell 输出列名为 Module_1_UCell 等
names(filtered_signatures) <- paste0("Module_", 1:4)
sig_cols <- paste0(names(filtered_signatures), "_UCell")

# 切片列表变更为仅保留 T353, T354
files_list <- c("Mouse2_T353.h5ad", "Mouse2_T354.h5ad")
meta_list <- list()

for (f in files_list) {
  reg_name <- gsub(".h5ad", "", f)
  cat(" 正在处理:", reg_name, "...\n")
  
  # A. 读取
  obj <- schard::h5ad2seurat(f)
  
  # B. 计算 UCell (92核强力输出)
  obj <- UCell::AddModuleScore_UCell(
    obj = obj,
    features = filtered_signatures,
    ncores = 92,
    name = "_UCell"
  )
  
  # C. 补充 region 信息并提取 metadata
  meta <- obj@meta.data
  meta$region <- reg_name
  meta_list[[reg_name]] <- meta
  
  # D. 释放内存
  rm(obj)
  gc()
}

meta_all <- bind_rows(meta_list)
cat("✅ Meta 数据已聚合完毕！\n")

# ==========================================
# 2. KDE 密度分布核心算法
# ==========================================
weighted_var <- function(x, w) {
  w <- w / sum(w)
  mu <- sum(w * x)
  sum(w * (x - mu)^2)
}

scott_bw_1d <- function(x, w) {
  w <- w / sum(w)
  n_eff <- 1 / sum(w^2)
  sd_w <- sqrt(weighted_var(x, w))
  if (!is.finite(sd_w) || sd_w == 0) return(NA_real_)
  sd_w * n_eff^(-1/5)
}

kde_gaussian_1d <- function(x, w, grid, bw) {
  w <- w / sum(w)
  z <- outer(grid, x, "-") / bw
  as.vector((dnorm(z) / bw) %*% w)
}

ucell_density_from_meta <- function(
    meta_df,
    signatures = sig_cols,
    dist_col = "distance",
    area_col = "area",
    keep_area = c("granular", "molecular"),
    group_col = "region",
    bin_width = 0.01,
    grid = NULL,
    n_grid = 300
) {
  df <- meta_df %>%
    filter(.data[[area_col]] %in% keep_area, !is.na(.data[[dist_col]])) %>%
    transmute(
      distance = as.numeric(.data[[dist_col]]),
      area = .data[[area_col]],
      group = as.character(.data[[group_col]]),
      across(all_of(signatures), as.numeric)
    )
  
  if (is.null(grid)) {
    rng <- range(df$distance, na.rm = TRUE)
    grid <- c(rng[1], rng[2])
  }
  
  grid_x <- seq(grid[1], grid[2], length.out = n_grid)
  
  df <- df %>%
    mutate(distance_bin = as.integer(distance / bin_width) * bin_width)
  
  agg <- df %>%
    group_by(group, distance_bin) %>%
    summarise(
      num = n(),
      across(all_of(signatures), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    rename(distance = distance_bin)
  
  out <- agg %>%
    group_by(group) %>%
    group_modify(~{
      dat <- .x
      map_dfr(signatures, function(sig) {
        x <- dat$distance
        w_expr <- dat[[sig]] + 0.001
        bw <- scott_bw_1d(x, w_expr)
        if (!is.finite(bw) || bw <= 0 || length(unique(x)) < 2) {
          return(tibble(distance = grid_x, value = NA_real_, signature = sig))
        }
        dens_expr <- kde_gaussian_1d(x, w_expr, grid_x, bw)
        dens_num <- kde_gaussian_1d(x, dat$num, grid_x, bw)
        tibble(distance = grid_x, value = dens_expr / dens_num, signature = sig)
      })
    }) %>%
    ungroup()
  
  return(out)
}

calculate_spatial_conservation <- function(density_df) {
  density_df %>%
    group_by(signature) %>%
    summarise(
      label = {
        wide_dat <- pivot_wider(
          pick(distance, group, value),
          names_from = group,
          values_from = value
        ) %>%
          drop_na()
        mat <- as.matrix(wide_dat %>% select(-distance))
        if (ncol(mat) >= 2) {
          cor_res <- cor(mat, method = "pearson")
          mean_r <- mean(cor_res[upper.tri(cor_res)])
          sprintf("Mean R = %.2f\nP < 0.001", mean_r)
        } else {
          ""
        }
      },
      .groups = "drop"
    )
}

cat("📊 正在计算平滑密度与切片间相关性...\n")
density_df <- ucell_density_from_meta(meta_df = meta_all)
stats_df <- calculate_spatial_conservation(density_df)

# ==========================================
# 3. Nature 级绘图 (恢复严格 7pt + 完美比例)
# ==========================================
# 重命名 facet 标签为 "Module 1" 等
density_df$sig_label <- factor(
  density_df$signature,
  levels = sig_cols,
  labels = paste0("Module ", 1:4)
)

stats_df$asig_label <- factor(
  stats_df$signature,
  levels = sig_cols,
  labels = paste0("Module ", 1:4)
)

# 完全匹配图中的 5 个等分刻度
exact_5_breaks <- function(x) seq(x[1], x[2], length.out = 5)

plot_nature_density <- function(density_df, stats_df) {
  
  # 清洗样本名称，抹除 "Mouse2_" 前缀
  density_df$group <- gsub("Mouse2_", "", density_df$group)
  
  x_min <- min(density_df$distance, na.rm = TRUE)
  
  ggplot(density_df, aes(x = distance, y = value, color = group)) +
    # 使用半透明 (alpha=0.4) 的特调色值，透出虚线
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
             fill = "#E6EDCF", alpha = 0.4) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "#C8E1DC", alpha = 0.4) +
    # 0 点分割虚线
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", color = "#888888") +
    # 主体折线图 (置于最上层)
    geom_line(linewidth = 0.7, na.rm = TRUE, alpha = 0.9) +
    # 统计学文字标注 (锁死 7pt)
    geom_text(
      data = stats_df,
      aes(x = x_min + 0.05, y = Inf, label = label),
      inherit.aes = FALSE,
      vjust = 1.3,
      hjust = 0,
      size = 7 / .pt,
      color = "black",
      family = family_use
    ) +
    facet_wrap(~ sig_label, ncol = 2, scales = "free_y") +
    theme_classic() +
    theme(
      # 全局锁死 7pt 字体
      text = element_text(size = 7, color = "black", face = "plain", family = family_use),
      axis.text = element_text(size = 7, color = "black", face = "plain", family = family_use),
      axis.title = element_text(
        size = 7, color = "black", face = "plain",
        family = family_use, margin = margin(t = 5, r = 5)
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        size = 7, color = "black", face = "plain",
        family = family_use, hjust = 0, margin = margin(b = 3)
      ),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 7, color = "black", face = "plain", family = family_use),
      legend.key.size = unit(0.4, "cm"),
      legend.margin = margin(b = -5),
      # 轴线和虚线
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      panel.grid.major.y = element_line(
        color = "#A0A0A0", linetype = "dashed", linewidth = 0.4
      ),
      # 强制宽度比高度大 1/4 (aspect.ratio = 0.8)
      aspect.ratio = 0.8,
      panel.spacing = unit(1, "lines")
    ) +
    scale_y_continuous(breaks = exact_5_breaks, labels = function(x) sprintf("%.2f", x)) +
    # 颜色映射字典仅保留 T353 和 T354
    scale_color_manual(values = c("T353" = "#aa363d", "T354" = "#228fbd")) +
    labs(x = "Distance from boundary", y = "Smoothed UCell scores")
}

# ==========================================
# 4. 重新导出矢量 PDF (画幅放大 50%)
# ==========================================
p_final <- plot_nature_density(density_df, stats_df)

cat("💾 正在通过 cairo_pdf 导出 1.5倍 大小的终极矢量图...\n")

# 宽度 9.75，高度 8.25
cairo_pdf("Figure_Density_Filtered_Signatures_Large_7pt.pdf",
          width = 9.75, height = 8.25, family = family_use)
print(p_final)
dev.off()

cat("🎉 T353 和 T354 的双样本对比图已生成！\n")







#fig5
library(Seurat)
library(UCell)
library(stringr)
library(homologene)

# ============================================================
# 0. 科学同源转换函数 (Human -> Mouse homolog mapping)
# ============================================================
get_human_mouse_map <- function(human_genes, seurat_mouse) {
  
  trans_df <- homologene(
    human_genes,
    inTax = 9606,   # Human
    outTax = 10090  # Mouse
  )
  
  trans_df <- trans_df[
    !is.na(trans_df$`9606`) & !is.na(trans_df$`10090`),
  ]
  
  # 使用 rownames(seurat_mouse) 兼容各种默认 Assay
  trans_df <- trans_df[
    trans_df$`10090` %in% rownames(seurat_mouse),
  ]
  
  trans_df <- trans_df[
    !duplicated(trans_df$`9606`),
  ]
  
  return(trans_df)
}

# ==========================================
# 1. 顺延修改 Module 列名 (2345 改为 1234)
# ==========================================
meta_names <- colnames(seurat_obj@meta.data)
meta_names[meta_names == "Module_2_UCell"] <- "Module_1_UCell"
meta_names[meta_names == "Module_3_UCell"] <- "Module_2_UCell"
meta_names[meta_names == "Module_4_UCell"] <- "Module_3_UCell"
meta_names[meta_names == "Module_5_UCell"] <- "Module_4_UCell"
colnames(seurat_obj@meta.data) <- meta_names

# ==========================================
# 2. 批量读取 txt 并进行 Homolog 科学转换
# ==========================================
txt_files <- list.files(pattern = "\\.txt$")
disease_signatures <- list()

for (f in txt_files) {
  disease_name <- gsub("\\.txt$", "", f)
  
  # 逐行读取文件内容
  lines <- readLines(f, warn = FALSE)
  
  # 核心正则：只保留以数字开头、且后面跟着空白字符的行
  valid_lines <- lines[grepl("^[0-9]+\\s+", lines)]
  
  if (length(valid_lines) > 0) {
    # 提取人类基因名
    genes_human <- sapply(strsplit(valid_lines, "\\s+"), function(x) x[2])
    
    # 调用 homologene 进行严格转换与过滤
    map_res <- get_human_mouse_map(genes_human, seurat_obj)
    genes_mouse <- map_res$`10090`
    
    # 只有存在有效同源基因时才加入列表
    if (length(genes_mouse) > 0) {
      disease_signatures[[disease_name]] <- genes_mouse
    } else {
      cat("⚠️ 警告:", disease_name, "中的基因在转换为小鼠并过滤后为空，已跳过。\n")
    }
  }
}

# ==========================================
# 3. 生成科学映射后的疾病与基因 CSV 文件
# ==========================================
cat("📁 正在生成基于 homologene 转换后的 CSV 文件...\n")

# 获取最长基因集的长度
max_len <- max(lengths(disease_signatures))

# 将 list 转换为按行排列的矩阵
csv_data <- do.call(rbind, lapply(names(disease_signatures), function(d) {
  genes <- disease_signatures[[d]]
  c(d, genes, rep(NA, max_len - length(genes)))
}))

# 赋予列名
colnames(csv_data) <- c("Disease", paste0("Gene_", 1:max_len))

# 写入文件
write.csv(csv_data, file = "Disease_Mouse_Homolog_Genes.csv", row.names = FALSE, na = "")

cat("✅ CSV 文件 'Disease_Mouse_Homolog_Genes.csv' 已保存！\n")

# ==========================================
# 4. 运行 UCell 计算疾病评分
# ==========================================
cat("📊 正在计算疾病基因集 UCell 分数...\n")

seurat_obj <- AddModuleScore_UCell(
  obj = seurat_obj,
  features = disease_signatures,
  ncores = 96,
  name = "_UCell"
)

cat("🎉 同源基因严谨转换、CSV 导出及 UCell 评分全部完成！\n")



#fig4
load('Figure5_Correlation_Matrices_Complete.RData')
library(ggplot2)
library(dplyr)

# ==========================================
# 0. 跨平台字体检测与全局 7pt 设置
# ==========================================
available_fonts <- names(pdfFonts())
family_use <- if ("Arial" %in% available_fonts) "Arial" else "sans"

font_style_7pt <- element_text(size = 7, color = "black", face = "plain", family = family_use)

cat("⏳ 正在加载关联矩阵数据...\n")
load("Figure5_Correlation_Matrices_Complete.RData")

# ==========================================
# 1. 提取名称与维度参数
# ==========================================
dis_names <- colnames(res_dis_dis$cor)
cell_names <- colnames(res_cell_cell$cor)
mod_names <- rownames(res_mod_dis$cor)

n_dis <- length(dis_names)   # 16
n_cell <- length(cell_names) # 18

# ==========================================
# 2. 全局绝对坐标与画幅锁定
# ==========================================
# 强制统一 4 个 Module 的绝对 Y 坐标
y_mod_shared <- seq(-4, -16, length.out = length(mod_names))

# 背景底框的统一大小
bg_box_size <- 4

# 强制统一 Y 轴跨度
global_ylim <- c(-19, 1)

# 强制统一 X 轴跨度
left_xlim <- c(-2, 28)
right_xlim <- c(-8, 22)

# 统一的 PDF 画布大小
pdf_w <- 6
pdf_h <- 5.5

# ======================================================================
# Part A: 左半部分 (Diseases)
# ======================================================================
df_dis_box <- data.frame()
for(i in 1:n_dis) {
  for(j in 1:n_dis) {
    if(i > j) {
      df_dis_box <- rbind(df_dis_box, data.frame(
        x = j,
        y = -i,
        cor = res_dis_dis$cor[i, j],
        abs_cor = abs(res_dis_dis$cor[i, j])
      ))
    }
  }
}

df_diag_dis <- data.frame(x = 1:n_dis, y = -(1:n_dis))
df_label_dis <- data.frame(x = 0, y = -(1:n_dis), label = dis_names)

x_mod_right <- n_dis + 5
df_mod_dis <- data.frame(x = x_mod_right, y = y_mod_shared, label = mod_names)

df_link_dis <- data.frame()
for(m in 1:length(mod_names)) {
  for(d in 1:n_dis) {
    if(res_mod_dis$p[m, d] < 0.05) {
      df_link_dis <- rbind(df_link_dis, data.frame(
        x_mod = x_mod_right,
        y_mod = y_mod_shared[m],
        x_diag = d,
        y_diag = -d,
        cor = res_mod_dis$cor[m, d],
        abs_cor = abs(res_mod_dis$cor[m, d])
      ))
    }
  }
}

p_disease <- ggplot() +
  geom_curve(data = df_link_dis, 
             aes(x = x_mod, y = y_mod, xend = x_diag, yend = y_diag, 
                 color = cor, linewidth = abs_cor),
             curvature = 0.15, alpha = 0.5) +
  # 满框底纹
  geom_point(data = df_dis_box, aes(x = x, y = y), 
             shape = 22, size = bg_box_size, 
             color = "grey60", fill = "white", stroke = 0.4) +
  geom_point(data = df_dis_box, aes(x = x, y = y, color = cor, size = abs_cor), 
             shape = 15) +
  geom_point(data = df_diag_dis, aes(x = x, y = y), size = 1, color = "red") +
  geom_point(data = df_mod_dis, aes(x = x, y = y), size = 1.5, color = "red") +
  geom_text(data = df_label_dis, aes(x = x, y = y, label = label), 
            hjust = 1, size = 7/.pt, family = family_use, color = "black") +
  geom_text(data = df_mod_dis, aes(x = x + 0.5, y = y, label = label), 
            hjust = 0, size = 7/.pt, family = family_use, color = "black") +
  scale_color_gradient2(low = "#4393C3", mid = "white", high = "#D6604D", 
                        midpoint = 0, limits = c(-1, 1), name = "Pearson's r") +
  scale_size_continuous(range = c(0.5, bg_box_size - 0.5), name = "|r|") +
  scale_linewidth_continuous(range = c(0.2, 1), guide = "none") +
  coord_fixed(xlim = left_xlim, ylim = global_ylim, clip = "off") +
  theme_void() +
  theme(
    text = font_style_7pt,
    legend.position = "left",
    legend.title = font_style_7pt,
    legend.text = font_style_7pt,
    legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  )

cairo_pdf("Figure5_Mantel_Left_Diseases_Grid.pdf", width = pdf_w, height = pdf_h, family = family_use)
print(p_disease)
dev.off()

# ======================================================================
# Part B: 右半部分 (Cells)
# ======================================================================
df_cell_box <- data.frame()
for(i in 1:n_cell) {
  for(j in 1:n_cell) {
    if(i < j) {
      df_cell_box <- rbind(df_cell_box, data.frame(
        x = j,
        y = -i,
        cor = res_cell_cell$cor[i, j],
        abs_cor = abs(res_cell_cell$cor[i, j])
      ))
    }
  }
}

df_diag_cell <- data.frame(x = 1:n_cell, y = -(1:n_cell))
df_label_cell <- data.frame(x = n_cell + 1, y = -(1:n_cell), label = cell_names)

x_mod_left <- -5
df_mod_cell <- data.frame(x = x_mod_left, y = y_mod_shared, label = mod_names)

df_link_cell <- data.frame()
for(m in 1:length(mod_names)) {
  for(c in 1:n_cell) {
    if(res_mod_cell$p[m, c] < 0.05) {
      df_link_cell <- rbind(df_link_cell, data.frame(
        x_mod = x_mod_left,
        y_mod = y_mod_shared[m],
        x_diag = c,
        y_diag = -c,
        cor = res_mod_cell$cor[m, c],
        abs_cor = abs(res_mod_cell$cor[m, c])
      ))
    }
  }
}

p_cell <- ggplot() +
  geom_curve(data = df_link_cell, 
             aes(x = x_mod, y = y_mod, xend = x_diag, yend = y_diag, 
                 color = cor, linewidth = abs_cor),
             curvature = -0.15, alpha = 0.5) +
  geom_point(data = df_cell_box, aes(x = x, y = y), 
             shape = 22, size = bg_box_size, 
             color = "grey60", fill = "white", stroke = 0.4) +
  geom_point(data = df_cell_box, aes(x = x, y = y, color = cor, size = abs_cor), 
             shape = 15) +
  geom_point(data = df_diag_cell, aes(x = x, y = y), size = 1, color = "red") +
  geom_point(data = df_mod_cell, aes(x = x, y = y), size = 1.5, color = "red") +
  geom_text(data = df_label_cell, aes(x = x, y = y, label = label), 
            hjust = 0, size = 7/.pt, family = family_use, color = "black") +
  geom_text(data = df_mod_cell, aes(x = x - 0.5, y = y, label = label), 
            hjust = 1, size = 7/.pt, family = family_use, color = "black") +
  scale_color_gradient2(low = "#4393C3", mid = "white", high = "#D6604D", 
                        midpoint = 0, limits = c(-1, 1), name = "Pearson's r") +
  scale_size_continuous(range = c(0.5, bg_box_size - 0.5), name = "|r|") +
  scale_linewidth_continuous(range = c(0.2, 1), guide = "none") +
  coord_fixed(xlim = right_xlim, ylim = global_ylim, clip = "off") +
  theme_void() +
  theme(
    text = font_style_7pt,
    legend.position = "right",
    legend.title = font_style_7pt,
    legend.text = font_style_7pt,
    legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  )

cairo_pdf("Figure5_Mantel_Right_Cells_Grid.pdf", width = pdf_w, height = pdf_h, family = family_use)
print(p_cell)
dev.off()

cat("✅ 物理对齐完毕！导出的两个 PDF 长宽尺寸和内部缩放比完全一样了！\n")









library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)

# ==========================================
# 1. 下载 GSE8397
# ==========================================
gse <- getGEO(
  "GSE8397",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)

length(gse)

# 看两个平台
sapply(gse, annotation)

library(GEOquery)
library(Biobase)
library(dplyr)

# ============================================================
# 1. 先使用 GPL96
# ============================================================

eset96 <- gse[[1]]

meta96 <- pData(eset96)

# 看全部样本名称
data.frame(
  GSM   = meta96$geo_accession,
  Title = meta96$title
)

# ============================================================
# 2. 提取 medial substantia nigra
# ============================================================

keep_medial <- grepl(
  "medial",
  meta96$title,
  ignore.case = TRUE
) & grepl(
  "nigra",
  meta96$title,
  ignore.case = TRUE
)

meta_medial <- meta96[keep_medial, , drop = FALSE]

# ============================================================
# 3. 根据 title 定义 PD / Control
# ============================================================

meta_medial$Group <- ifelse(
  grepl(
    "control",
    meta_medial$title,
    ignore.case = TRUE
  ),
  "Control",
  "PD"
)

meta_medial$Group <- factor(
  meta_medial$Group,
  levels = c("Control", "PD")
)

# ============================================================
# 4. 检查最终样本
# ============================================================

sample_info <- data.frame(
  GSM   = meta_medial$geo_accession,
  Title = meta_medial$title,
  Group = meta_medial$Group
)

print(sample_info)

cat("\n========== Sample number ==========\n")
print(table(meta_medial$Group))
















library(GEOquery)
library(Biobase)
library(dplyr)
library(tidyr)

# ==========================================
# 1. 使用 GPL96
# ==========================================
eset96 <- gse[[1]]

expr96 <- exprs(eset96)
meta96 <- pData(eset96)
anno96 <- fData(eset96)

# ==========================================
# 2. 提取 medial substantia nigra
# ==========================================
keep_medial <- grepl(
  "medial",
  meta96$title,
  ignore.case = TRUE
) &
  grepl(
    "nigra",
    meta96$title,
    ignore.case = TRUE
  )

expr_medial <- expr96[, keep_medial, drop = FALSE]
meta_medial <- meta96[keep_medial, , drop = FALSE]

# ==========================================
# 3. PD / Control
# ==========================================
meta_medial$Group <- ifelse(
  grepl(
    "control",
    meta_medial$title,
    ignore.case = TRUE
  ),
  "Control",
  "PD"
)

meta_medial$Group <- factor(
  meta_medial$Group,
  levels = c("Control", "PD")
)

print(table(meta_medial$Group))

# ==========================================
# 4. 检查表达值范围
# ==========================================
summary(as.numeric(expr_medial))

quantile(
  as.numeric(expr_medial),
  probs = c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1),
  na.rm = TRUE
)









library(Seurat)
library(ggplot2)
library(homologene)
library(dplyr)








# ==========================================
# 1. 手动定义 4 个 Mfuzz cluster
# ==========================================

mfuzz_modules_mouse <- list(
  
  Cluster1 = c(
    "Slc1a3", "Ppp1r17", "Homer3", "Fam107a", "Atp1b1",
    "Icmt", "Gad1", "Gpr37l1", "Gdf10", "Thy1"
  ),
  
  Cluster2 = c(
    "Itpr1", "Dner", "Inpp5a", "Prkcg", "Gabra1", "Sptbn2"
  ),
  
  Cluster3 = c(
    "Pcp4", "Calb1", "Atp2a2", "Rgs8", "Lhx1os", "Slc1a6"
  ),
  
  Cluster4 = c(
    "Car8", "Pcp2", "Pvalb", "Gng13", "Nsg1",
    "Ywhah", "Cck", "Id2", "Itm2b", "Ckb"
  )
)


# ==========================================
# 2. mouse -> human homologene
# ==========================================

all_mouse_genes <- unique(
  unlist(mfuzz_modules_mouse)
)

gene_map <- homologene(
  genes = all_mouse_genes,
  inTax = 10090,
  outTax = 9606
)

head(gene_map)



# ==========================================
# 3. 建立 mouse-human 对照表
# ==========================================

map_df <- data.frame(
  Mouse = gene_map$`10090`,
  Human = gene_map$`9606`,
  stringsAsFactors = FALSE
)

map_df$Human <- toupper(map_df$Human)

print(map_df)


# ==========================================
# 4. 每个 Mfuzz cluster 转成人类 gene
# ==========================================

mfuzz_modules_human <- lapply(
  mfuzz_modules_mouse,
  function(x) {
    
    y <- map_df$Human[
      match(x, map_df$Mouse)
    ]
    
    y <- y[
      !is.na(y) &
        y != ""
    ]
    
    unique(y)
  }
)

mfuzz_modules_human




















############################################################
# 从这里继续
# 已经存在：
# expr_medial
# meta_medial
# anno96
# mfuzz_modules_human
#
# 目标：
# GPL96 probe -> human gene
# -> ssGSEA 4个Mfuzz module
# -> Medial substantia nigra
# -> PD vs Control
# -> Boxplot + P value
############################################################

library(Biobase)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GSVA)
library(ggpubr)


# ==========================================================
# 5. GPL96 Probe -> Human Gene Symbol
# ==========================================================

# GEO AnnotGPL 已经有 Gene Symbol
gene_symbol <- as.character(
  anno96$`Gene Symbol`
)

# 有些 probe 对应多个基因，例如：
# DDR1 /// MIR4640
# 这里取第一个 gene symbol
gene_symbol <- sub(
  " ///.*$",
  "",
  gene_symbol
)

# 清理
keep_gene <- !is.na(gene_symbol) &
  gene_symbol != "" &
  gene_symbol != "---"

expr_tmp <- expr_medial[
  keep_gene,
  ,
  drop = FALSE
]

gene_symbol_tmp <- gene_symbol[
  keep_gene
]


# ==========================================================
# 6. 多个 probe 对应一个 gene 时取平均
# ==========================================================

expr_df <- data.frame(
  Gene = gene_symbol_tmp,
  expr_tmp,
  check.names = FALSE
)

expr_gene_df <- expr_df %>%
  group_by(Gene) %>%
  summarise(
    across(
      where(is.numeric),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# 转 matrix
gene_names <- expr_gene_df$Gene

expr_gene_df$Gene <- NULL

expr_gene_mean <- as.matrix(
  expr_gene_df
)

rownames(expr_gene_mean) <- toupper(
  gene_names
)

# 列名改成 GSM
colnames(expr_gene_mean) <- meta_medial$geo_accession


# ==========================================================
# 7. 检查表达矩阵
# ==========================================================

dim(expr_gene_mean)

expr_gene_mean[1:5, 1:5]

summary(
  as.numeric(expr_gene_mean)
)


# ==========================================================
# 8. 建立样本 metadata
# ==========================================================

sample_meta <- data.frame(
  GSM = meta_medial$geo_accession,
  Group = meta_medial$Group,
  Title = meta_medial$title,
  stringsAsFactors = FALSE
)

sample_meta$Group <- factor(
  sample_meta$Group,
  levels = c(
    "Control",
    "PD"
  )
)

rownames(sample_meta) <- sample_meta$GSM

print(
  table(sample_meta$Group)
)

# 应该：
# Control = 8
# PD      = 15


# ==========================================================
# 9. 检查4个Mfuzz module在人类黑质数据中的基因
# ==========================================================

mfuzz_modules_human_detected <- lapply(
  mfuzz_modules_human,
  function(x) {
    
    intersect(
      toupper(x),
      rownames(expr_gene_mean)
    )
  }
)


# 每个 module 检测到多少基因
cat("\n========== Module genes detected ==========\n")

for (i in names(mfuzz_modules_human_detected)) {
  
  cat(
    "\n",
    i,
    ": ",
    length(mfuzz_modules_human_detected[[i]]),
    "/",
    length(mfuzz_modules_human[[i]]),
    "\n",
    sep = ""
  )
  
  print(
    mfuzz_modules_human_detected[[i]]
  )
}


# ==========================================================
# 10. 删除完全没有基因的 module
# 正常情况下应该不会有
# ==========================================================

mfuzz_modules_human_detected <- mfuzz_modules_human_detected[
  lengths(mfuzz_modules_human_detected) > 0
]


# ==========================================================
# 11. 保存 mouse-human-module 对应关系
# ==========================================================

module_gene_table <- bind_rows(
  lapply(
    names(mfuzz_modules_mouse),
    function(clu) {
      
      mouse <- mfuzz_modules_mouse[[clu]]
      
      human <- map_df$Human[
        match(
          mouse,
          map_df$Mouse
        )
      ]
      
      data.frame(
        Cluster = clu,
        Mouse_gene = mouse,
        Human_gene = human,
        stringsAsFactors = FALSE
      )
    }
  )
)

module_gene_table$Detected_GSE8397 <- (
  toupper(module_gene_table$Human_gene) %in%
    rownames(expr_gene_mean)
)

write.csv(
  module_gene_table,
  "GSE8397_Mfuzz4_Mouse_Human_GeneMapping.csv",
  row.names = FALSE
)

print(module_gene_table)


# ==========================================================
# 12. 准备 ssGSEA 输入
# ==========================================================

expr_ssgsea <- expr_gene_mean

storage.mode(expr_ssgsea) <- "numeric"

# 删除异常值
expr_ssgsea <- expr_ssgsea[
  apply(
    expr_ssgsea,
    1,
    function(x) all(is.finite(x))
  ),
  ,
  drop = FALSE
]

dim(expr_ssgsea)


# ==========================================================
# 13. ssGSEA
# 自动兼容新版 / 旧版 GSVA
# ==========================================================

if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
  
  # ----------------------------
  # 新版 GSVA
  # ----------------------------
  
  ssgsea_param <- GSVA::ssgseaParam(
    exprData = expr_ssgsea,
    geneSets = mfuzz_modules_human_detected,
    normalize = TRUE
  )
  
  ssgsea_score <- GSVA::gsva(
    ssgsea_param,
    verbose = FALSE
  )
  
} else {
  
  # ----------------------------
  # 老版 GSVA
  # ----------------------------
  
  ssgsea_score <- GSVA::gsva(
    expr = expr_ssgsea,
    gset.idx.list = mfuzz_modules_human_detected,
    method = "ssgsea",
    kcdf = "Gaussian",
    abs.ranking = TRUE,
    verbose = FALSE
  )
}


# ==========================================================
# 14. 查看ssGSEA结果
# ==========================================================

print(
  ssgsea_score
)

dim(
  ssgsea_score
)

# 行 = Cluster
# 列 = Sample


# ==========================================================
# 15. 转成长格式
# ==========================================================

score_df <- as.data.frame(
  t(ssgsea_score)
)

score_df$GSM <- rownames(
  score_df
)

# 加入 PD / Control
score_df$Group <- sample_meta[
  score_df$GSM,
  "Group"
]

score_long <- score_df %>%
  pivot_longer(
    cols = all_of(
      rownames(ssgsea_score)
    ),
    names_to = "Cluster",
    values_to = "ssGSEA_score"
  )


# 顺序固定
score_long$Group <- factor(
  score_long$Group,
  levels = c(
    "Control",
    "PD"
  )
)

score_long$Cluster <- factor(
  score_long$Cluster,
  levels = c(
    "Cluster1",
    "Cluster2",
    "Cluster3",
    "Cluster4"
  )
)


# ==========================================================
# 16. 看一下数据
# ==========================================================

print(
  head(score_long)
)

print(
  table(
    score_long$Cluster,
    score_long$Group
  )
)


# ==========================================================
# 17. 每个Cluster做 Wilcoxon
# PD vs Control
# ==========================================================

stat_df <- score_long %>%
  group_by(Cluster) %>%
  summarise(
    
    Control_mean = mean(
      ssGSEA_score[
        Group == "Control"
      ],
      na.rm = TRUE
    ),
    
    PD_mean = mean(
      ssGSEA_score[
        Group == "PD"
      ],
      na.rm = TRUE
    ),
    
    Control_median = median(
      ssGSEA_score[
        Group == "Control"
      ],
      na.rm = TRUE
    ),
    
    PD_median = median(
      ssGSEA_score[
        Group == "PD"
      ],
      na.rm = TRUE
    ),
    
    P_value = wilcox.test(
      ssGSEA_score ~ Group,
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  )


# BH校正
stat_df$FDR <- p.adjust(
  stat_df$P_value,
  method = "BH"
)

print(stat_df)


# 保存
write.csv(
  stat_df,
  "GSE8397_MedialSN_Mfuzz4_ssGSEA_PD_vs_Control_statistics.csv",
  row.names = FALSE
)


# ==========================================================
# 18. 保存每个样本的 ssGSEA score
# ==========================================================

write.csv(
  score_long,
  "GSE8397_MedialSN_Mfuzz4_ssGSEA_scores.csv",
  row.names = FALSE
)


# ==========================================================
# 19. Boxplot
# Control 灰色
# PD 蓝色
# 标注 Wilcoxon P值
# ==========================================================

p <- ggplot(
  score_long,
  aes(
    x = Group,
    y = ssGSEA_score,
    fill = Group
  )
) +
  
  geom_boxplot(
    width = 0.55,
    linewidth = 0.6,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    aes(
      color = Group
    ),
    width = 0.12,
    size = 1.8,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  
  stat_compare_means(
    comparisons = list(
      c(
        "Control",
        "PD"
      )
    ),
    method = "wilcox.test",
    label = "p.format",
    size = 4.2
  ) +
  
  facet_wrap(
    ~ Cluster,
    nrow = 1,
    scales = "free_y"
  ) +
  
  scale_fill_manual(
    values = c(
      "Control" = "#BDBDBD",
      "PD" = "#8491B4"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "Control" = "#808080",
      "PD" = "#4DBBD5"
    )
  ) +
  
  labs(
    x = NULL,
    y = "ssGSEA score"
  ) +
  
  theme_classic(
    base_size = 14,
    base_family = "Arial"
  ) +
  
  theme(
    
    text = element_text(
      size = 14,
      color = "black",
      face = "plain"
    ),
    
    axis.text = element_text(
      size = 14,
      color = "black",
      face = "plain"
    ),
    
    axis.title = element_text(
      size = 14,
      color = "black",
      face = "plain"
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 14,
      color = "black",
      face = "plain"
    ),
    
    legend.position = "none",
    
    axis.line = element_line(
      color = "black",
      linewidth = 0.6
    )
  )


# ==========================================================
# 20. 输出 4个横排 boxplot
# ==========================================================

ggsave(
  "GSE8397_MedialSN_Mfuzz4_ssGSEA_PD_vs_Control.pdf",
  p,
  width = 10,
  height = 4,
  device = cairo_pdf
)

print(p)


# ==========================================================
# 21. 如果你想和你发的示意图一样
# 每个 Cluster 单独纵向排列
# ==========================================================

p_vertical <- ggplot(
  score_long,
  aes(
    x = Group,
    y = ssGSEA_score,
    fill = Group
  )
) +
  
  geom_boxplot(
    width = 0.55,
    linewidth = 0.6,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    aes(
      color = Group
    ),
    width = 0.12,
    size = 1.5,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  
  stat_compare_means(
    comparisons = list(
      c(
        "Control",
        "PD"
      )
    ),
    method = "wilcox.test",
    label = "p.format",
    size = 4
  ) +
  
  facet_wrap(
    ~ Cluster,
    ncol = 1,
    scales = "free_y"
  ) +
  
  scale_fill_manual(
    values = c(
      "Control" = "#BDBDBD",
      "PD" = "#8491B4"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "Control" = "#808080",
      "PD" = "#4DBBD5"
    )
  ) +
  
  labs(
    x = NULL,
    y = "ssGSEA score"
  ) +
  
  theme_classic(
    base_size = 14,
    base_family = "Arial"
  ) +
  
  theme(
    
    text = element_text(
      size = 14,
      color = "black",
      face = "plain"
    ),
    
    axis.text = element_text(
      size = 14,
      color = "black",
      face = "plain"
    ),
    
    axis.title = element_text(
      size = 14,
      color = "black",
      face = "plain"
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 14,
      color = "black",
      face = "plain"
    ),
    
    legend.position = "none"
  )


ggsave(
  "GSE8397_MedialSN_Mfuzz4_ssGSEA_PD_vs_Control_vertical.pdf",
  p_vertical,
  width = 3.6,
  height = 10,
  device = cairo_pdf
)

print(p_vertical)


























library(ggplot2)
library(ggpubr)
library(dplyr)

# ==========================================
# 1. Windows 注册 Arial
# ==========================================

windowsFonts(
  Arial = windowsFont("Arial")
)

# 检查
windowsFonts()


# ==========================================
# 2. 保证分组顺序
# ==========================================

score_long$Group <- factor(
  score_long$Group,
  levels = c("Control", "PD")
)

score_long$Cluster <- factor(
  score_long$Cluster,
  levels = c(
    "Cluster1",
    "Cluster2",
    "Cluster3",
    "Cluster4"
  )
)


# ==========================================
# 3. 逐个 Cluster 作图
# ==========================================

for (clu in levels(score_long$Cluster)) {
  
  df_plot <- score_long %>%
    filter(Cluster == clu)
  
  p <- ggplot(
    df_plot,
    aes(
      x = Group,
      y = ssGSEA_score,
      fill = Group
    )
  ) +
    
    geom_boxplot(
      width = 0.55,
      linewidth = 0.45,
      outlier.shape = NA
    ) +
    
    stat_compare_means(
      comparisons = list(
        c("Control", "PD")
      ),
      method = "wilcox.test",
      label = "p.format",
      size = 2.5,
      family = "Arial"
    ) +
    
    scale_fill_manual(
      values = c(
        "Control" = "#C8C2B8",
        "PD" = "#7FA8B8"
      )
    ) +
    
    labs(
      x = NULL,
      y = "ssGSEA score",
      title = paste0(
        "GSE8397 (",
        clu,
        ")"
      )
    ) +
    
    theme_classic(
      base_size = 7,
      base_family = "Arial"
    ) +
    
    theme(
      text = element_text(
        family = "Arial",
        size = 7,
        color = "black",
        face = "plain"
      ),
      
      plot.title = element_text(
        family = "Arial",
        size = 7,
        color = "black",
        face = "plain",
        hjust = 0
      ),
      
      axis.title = element_text(
        family = "Arial",
        size = 7,
        color = "black",
        face = "plain"
      ),
      
      axis.text = element_text(
        family = "Arial",
        size = 7,
        color = "black",
        face = "plain"
      ),
      
      axis.line = element_line(
        linewidth = 0.45,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.4,
        color = "black"
      ),
      
      legend.position = "none"
    )
  
  # ========================================
  # 输出每个 Cluster 单独 PDF
  # ========================================
  
  ggsave(
    filename = paste0(
      "GSE8397_MedialSN_",
      clu,
      "_ssGSEA_PD_vs_Control.pdf"
    ),
    plot = p,
    width = 2.2,
    height = 2.4,
    device = cairo_pdf
  )
  
  print(p)
}


############################################################
############################################################
# GSE48350
# Alzheimer's disease validation
# RAW CEL -> RMA -> Gene expression
# -> Metadata
# -> Mouse-to-human Mfuzz clusters
# -> ssGSEA
# -> Control vs AD
# -> 4 independent PDF figures
############################################################
############################################################

rm(list = ls())

options(stringsAsFactors = FALSE)

############################################################
# 0. 路径
############################################################

base_dir <- "D:/deskup/毕业论文"

setwd(base_dir)

raw_tar <- file.path(
  base_dir,
  "GSE48350_RAW.tar"
)

raw_dir <- file.path(
  base_dir,
  "GSE48350_RAW"
)


############################################################
# 1. 加载包
############################################################

library(affy)
library(Biobase)

library(hgu133plus2.db)
library(AnnotationDbi)

library(GEOquery)

library(homologene)

library(GSVA)

library(tidyr)
library(ggplot2)


############################################################
# 2. 解压 RAW CEL
############################################################

dir.create(
  raw_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

cel_files <- list.files(
  raw_dir,
  pattern = "\\.CEL(\\.gz)?$",
  full.names = TRUE,
  ignore.case = TRUE
)

# 如果还没解压
if (length(cel_files) == 0) {
  
  untar(
    raw_tar,
    exdir = raw_dir
  )
  
  cel_files <- list.files(
    raw_dir,
    pattern = "\\.CEL(\\.gz)?$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
}

cat(
  "CEL file number:",
  length(cel_files),
  "\n"
)

stopifnot(
  length(cel_files) == 253
)


############################################################
# 3. 读取 CEL
############################################################

affy_raw <- ReadAffy(
  filenames = cel_files
)

affy_raw


############################################################
# 4. RMA
############################################################

eset_rma <- rma(
  affy_raw
)

expr_probe <- exprs(
  eset_rma
)

dim(expr_probe)

quantile(
  expr_probe,
  probs = c(
    0,
    0.01,
    0.25,
    0.5,
    0.75,
    0.99,
    1
  ),
  na.rm = TRUE
)

# 注意：
# RMA 本身已经是 log2 scale
# 不要再 log2


############################################################
# 5. Probe -> Gene Symbol
############################################################

probe_ids <- rownames(
  expr_probe
)

anno_df <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_ids,
  keytype = "PROBEID",
  columns = "SYMBOL"
)

anno_df <- anno_df[
  !is.na(anno_df$SYMBOL) &
    anno_df$SYMBOL != "",
  ,
  drop = FALSE
]

# 一个 probe 重复 annotation 时只保留一次
anno_df <- anno_df[
  !duplicated(anno_df$PROBEID),
  ,
  drop = FALSE
]

idx_probe <- match(
  anno_df$PROBEID,
  rownames(expr_probe)
)

expr_anno <- expr_probe[
  idx_probe,
  ,
  drop = FALSE
]

gene_symbol <- toupper(
  anno_df$SYMBOL
)


############################################################
# 6. 多 probe -> gene mean
############################################################

expr_gene_sum <- rowsum(
  expr_anno,
  group = gene_symbol,
  reorder = FALSE
)

probe_number <- table(
  gene_symbol
)

expr_gene_mean <- expr_gene_sum /
  as.numeric(
    probe_number[
      rownames(expr_gene_sum)
    ]
  )

dim(expr_gene_mean)

head(
  rownames(expr_gene_mean)
)


############################################################
# 7. 读取本地 GSE48350 series matrix metadata
############################################################

series_files <- list.files(
  base_dir,
  pattern = "GSE48350_series_matrix",
  full.names = TRUE,
  ignore.case = TRUE
)

print(series_files)

stopifnot(
  length(series_files) >= 1
)

series_file <- series_files[1]

gse_meta <- GEOquery::getGEO(
  filename = series_file,
  GSEMatrix = TRUE
)

if (
  inherits(
    gse_meta,
    "ExpressionSet"
  )
) {
  
  eset_meta <- gse_meta
  
} else {
  
  eset_meta <- gse_meta[[1]]
  
}

meta <- Biobase::pData(
  eset_meta
)

dim(meta)

colnames(meta)


############################################################
# 8. 合并 characteristics
############################################################

char_cols <- grep(
  "^characteristics_ch1",
  colnames(meta),
  value = TRUE
)

meta$all_char <- apply(
  meta[
    ,
    char_cols,
    drop = FALSE
  ],
  1,
  function(x) {
    
    paste(
      x[
        !is.na(x)
      ],
      collapse = " | "
    )
    
  }
)


############################################################
# 9. 建立 Control / AD 分组
############################################################
# GSE48350:
# title 中 indiv = normal control
# title 中 AD    = Alzheimer's disease
############################################################

meta$Group <- NA_character_

meta$Group[
  grepl(
    "_indiv",
    meta$title,
    ignore.case = TRUE
  )
] <- "Control"

meta$Group[
  grepl(
    "_AD",
    meta$title,
    ignore.case = TRUE
  )
] <- "AD"

meta$Group <- factor(
  meta$Group,
  levels = c(
    "Control",
    "AD"
  )
)

cat(
  "\nGroup distribution:\n"
)

print(
  table(
    meta$Group,
    useNA = "ifany"
  )
)

# 应该：
# Control = 173
# AD      = 80


############################################################
# 10. Brain region
############################################################

meta$BrainRegion <- tolower(
  as.character(
    meta$`brain region:ch1`
  )
)

# 统一两种写法
meta$BrainRegion[
  meta$BrainRegion %in% c(
    "postcentral gyrus",
    "post-central gyrus"
  )
] <- "post-central gyrus"

cat(
  "\nBrain region distribution:\n"
)

print(
  table(
    meta$BrainRegion,
    meta$Group,
    useNA = "ifany"
  )
)


############################################################
# 11. 从 CEL 列名提取 GSM
############################################################
# 注意：
# AD CEL:
# GSM1176196_xxxxx.CEL.gz
#
# Control CEL:
# GSM300166.CEL.gz
#
# 因此不能只按 "_" 截断
############################################################

sample_names <- colnames(
  expr_gene_mean
)

gsm <- sub(
  "^(GSM[0-9]+).*$",
  "\\1",
  basename(
    sample_names
  )
)

sample_map <- data.frame(
  CEL = sample_names,
  GSM = gsm,
  stringsAsFactors = FALSE
)

cat(
  "\nUnique GSM:",
  length(
    unique(
      sample_map$GSM
    )
  ),
  "\n"
)

cat(
  "Matched GSM:",
  sum(
    sample_map$GSM %in%
      meta$geo_accession
  ),
  "\n"
)

cat(
  "Unmatched GSM:\n"
)

print(
  setdiff(
    sample_map$GSM,
    meta$geo_accession
  )
)

# 应该：
# Unique GSM  = 253
# Matched GSM = 253
# Unmatched   = character(0)


############################################################
# 12. CEL 与 metadata 完整匹配
############################################################

sample_info <- data.frame(
  CEL = sample_map$CEL,
  GSM = sample_map$GSM,
  stringsAsFactors = FALSE
)

idx <- match(
  sample_info$GSM,
  meta$geo_accession
)

cat(
  "\nMissing metadata:",
  sum(
    is.na(idx)
  ),
  "\n"
)

stopifnot(
  sum(
    is.na(idx)
  ) == 0
)


############################################################
# 加入 metadata
############################################################

sample_info$Group <-
  meta$Group[idx]

sample_info$BrainRegion <-
  meta$BrainRegion[idx]

sample_info$Title <-
  meta$title[idx]

sample_info$Individual <-
  meta$`individual:ch1`[idx]

sample_info$Braak <-
  meta$`braak stage:ch1`[idx]

sample_info$Age <-
  meta$`age (yrs):ch1`[idx]

sample_info$Gender <-
  meta$`gender:ch1`[idx]


############################################################
# 检查
############################################################

cat(
  "\nFinal groups:\n"
)

print(
  table(
    sample_info$Group,
    useNA = "ifany"
  )
)

cat(
  "\nFinal brain regions:\n"
)

print(
  table(
    sample_info$BrainRegion,
    sample_info$Group,
    useNA = "ifany"
  )
)


############################################################
# 13. 最终 expression matrix
############################################################

expr_use <- expr_gene_mean[
  ,
  sample_info$CEL,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(expr_use),
    sample_info$CEL
  )
)

cat(
  "\nExpression matrix dimension:\n"
)

print(
  dim(
    expr_use
  )
)


############################################################
# 14. Mfuzz 4 clusters
############################################################

mfuzz_modules_mouse <- list(
  
  Cluster1 = c(
    "Slc1a3",
    "Ppp1r17",
    "Homer3",
    "Fam107a",
    "Atp1b1",
    "Icmt",
    "Gad1",
    "Gpr37l1",
    "Gdf10",
    "Thy1"
  ),
  
  Cluster2 = c(
    "Itpr1",
    "Dner",
    "Inpp5a",
    "Prkcg",
    "Gabra1",
    "Sptbn2"
  ),
  
  Cluster3 = c(
    "Pcp4",
    "Calb1",
    "Atp2a2",
    "Rgs8",
    "Lhx1os",
    "Slc1a6"
  ),
  
  Cluster4 = c(
    "Car8",
    "Pcp2",
    "Pvalb",
    "Gng13",
    "Nsg1",
    "Ywhah",
    "Cck",
    "Id2",
    "Itm2b",
    "Ckb"
  )
)


############################################################
# 15. Mouse -> Human homolog
############################################################

all_mouse_genes <- unique(
  unlist(
    mfuzz_modules_mouse
  )
)

gene_map <- homologene(
  genes = all_mouse_genes,
  inTax = 10090,
  outTax = 9606
)

map_df <- data.frame(
  Mouse = gene_map$`10090`,
  Human = toupper(
    gene_map$`9606`
  ),
  stringsAsFactors = FALSE
)

mfuzz_modules_human <- lapply(
  mfuzz_modules_mouse,
  function(x) {
    
    y <- map_df$Human[
      match(
        x,
        map_df$Mouse
      )
    ]
    
    unique(
      y[
        !is.na(y) &
          y != ""
      ]
    )
    
  }
)

cat(
  "\nHuman modules:\n"
)

print(
  mfuzz_modules_human
)


############################################################
# 16. 只保留 GSE48350 中检测到的 genes
############################################################

mfuzz_modules_detected <- lapply(
  mfuzz_modules_human,
  function(x) {
    
    intersect(
      x,
      rownames(
        expr_use
      )
    )
    
  }
)

cat(
  "\nDetected genes:\n"
)

print(
  mfuzz_modules_detected
)

cat(
  "\nDetected gene number:\n"
)

print(
  sapply(
    mfuzz_modules_detected,
    length
  )
)


############################################################
# 17. ssGSEA
############################################################

storage.mode(
  expr_use
) <- "numeric"

ssgsea_param <- GSVA::ssgseaParam(
  exprData = expr_use,
  geneSets = mfuzz_modules_detected,
  normalize = TRUE
)

ssgsea_score <- GSVA::gsva(
  ssgsea_param,
  verbose = FALSE
)

cat(
  "\nssGSEA matrix dimension:\n"
)

print(
  dim(
    ssgsea_score
  )
)

# 应该：
# 4 x 253


############################################################
# 18. 整理 score
############################################################

score_df <- as.data.frame(
  t(
    ssgsea_score
  )
)

score_df$CEL <- rownames(
  score_df
)

idx_score <- match(
  score_df$CEL,
  sample_info$CEL
)

stopifnot(
  sum(
    is.na(
      idx_score
    )
  ) == 0
)

score_df$GSM <-
  sample_info$GSM[
    idx_score
  ]

score_df$Group <-
  sample_info$Group[
    idx_score
  ]

score_df$BrainRegion <-
  sample_info$BrainRegion[
    idx_score
  ]

score_df$Age <-
  sample_info$Age[
    idx_score
  ]

score_df$Gender <-
  sample_info$Gender[
    idx_score
  ]

score_df$Braak <-
  sample_info$Braak[
    idx_score
  ]


############################################################
# Wide -> Long
############################################################

score_long <- tidyr::pivot_longer(
  score_df,
  cols = c(
    Cluster1,
    Cluster2,
    Cluster3,
    Cluster4
  ),
  names_to = "Cluster",
  values_to = "ssGSEA_score"
)

score_long$Group <- factor(
  score_long$Group,
  levels = c(
    "Control",
    "AD"
  )
)

cat(
  "\nSamples per cluster:\n"
)

print(
  table(
    score_long$Cluster,
    score_long$Group
  )
)


############################################################
# 19. All brain regions pooled
# Wilcoxon Control vs AD
############################################################

clusters <- c(
  "Cluster1",
  "Cluster2",
  "Cluster3",
  "Cluster4"
)

stat_list <- lapply(
  clusters,
  function(clu) {
    
    dat <- score_long[
      score_long$Cluster == clu,
      ,
      drop = FALSE
    ]
    
    x <- dat$ssGSEA_score[
      dat$Group == "Control"
    ]
    
    y <- dat$ssGSEA_score[
      dat$Group == "AD"
    ]
    
    p <- wilcox.test(
      x,
      y,
      exact = FALSE
    )$p.value
    
    data.frame(
      
      Cluster = clu,
      
      N_Control = length(
        x
      ),
      
      N_AD = length(
        y
      ),
      
      Control_mean = mean(
        x,
        na.rm = TRUE
      ),
      
      AD_mean = mean(
        y,
        na.rm = TRUE
      ),
      
      Control_median = median(
        x,
        na.rm = TRUE
      ),
      
      AD_median = median(
        y,
        na.rm = TRUE
      ),
      
      P_value = p
    )
    
  }
)

stat_df <- do.call(
  rbind,
  stat_list
)

stat_df$FDR <- p.adjust(
  stat_df$P_value,
  method = "BH"
)

stat_df$Direction <- ifelse(
  stat_df$AD_mean <
    stat_df$Control_mean,
  "AD lower",
  "AD higher"
)

cat(
  "\nFinal statistics:\n"
)

print(
  stat_df
)


############################################################
# 20. 保存结果
############################################################

write.csv(
  stat_df,
  file.path(
    base_dir,
    "GSE48350_AllBrainRegions_MfuzzClusters_ssGSEA_statistics.csv"
  ),
  row.names = FALSE
)

write.csv(
  score_long,
  file.path(
    base_dir,
    "GSE48350_AllBrainRegions_MfuzzClusters_ssGSEA_scores.csv"
  ),
  row.names = FALSE
)


############################################################
# 21. Arial
############################################################

if (
  .Platform$OS.type == "windows"
) {
  
  windowsFonts(
    Arial = windowsFont(
      "Arial"
    )
  )
  
}


############################################################
# 22. Colors
############################################################

group_colors <- c(
  "Control" = "#C8C2B8",
  "AD" = "#B56B79"
)


############################################################
# 23. 4 independent PDF
# 无散点
############################################################

for (
  clu in clusters
) {
  
  df_plot <- score_long[
    score_long$Cluster == clu,
    ,
    drop = FALSE
  ]
  
  x <- df_plot$ssGSEA_score[
    df_plot$Group == "Control"
  ]
  
  y <- df_plot$ssGSEA_score[
    df_plot$Group == "AD"
  ]
  
  p_value <- wilcox.test(
    x,
    y,
    exact = FALSE
  )$p.value
  
  
  ##########################################################
  # P value label
  ##########################################################
  
  if (
    p_value < 0.0001
  ) {
    
    p_label <- "p < 0.0001"
    
  } else {
    
    p_label <- paste0(
      "p = ",
      formatC(
        p_value,
        format = "f",
        digits = 4
      )
    )
    
  }
  
  
  ##########################################################
  # Y range
  ##########################################################
  
  ymax <- max(
    df_plot$ssGSEA_score,
    na.rm = TRUE
  )
  
  ymin <- min(
    df_plot$ssGSEA_score,
    na.rm = TRUE
  )
  
  yrange <- ymax - ymin
  
  if (
    !is.finite(
      yrange
    ) ||
    yrange == 0
  ) {
    
    yrange <- 1
    
  }
  
  bracket_y <-
    ymax +
    0.08 * yrange
  
  text_y <-
    ymax +
    0.15 * yrange
  
  
  ##########################################################
  # Plot
  ##########################################################
  
  p <- ggplot(
    df_plot,
    aes(
      x = Group,
      y = ssGSEA_score,
      fill = Group
    )
  ) +
    
    geom_boxplot(
      width = 0.52,
      linewidth = 0.4,
      outlier.shape = NA
    ) +
    
    annotate(
      "segment",
      x = 1,
      xend = 2,
      y = bracket_y,
      yend = bracket_y,
      linewidth = 0.35
    ) +
    
    annotate(
      "segment",
      x = 1,
      xend = 1,
      y = bracket_y -
        0.025 * yrange,
      yend = bracket_y,
      linewidth = 0.35
    ) +
    
    annotate(
      "segment",
      x = 2,
      xend = 2,
      y = bracket_y -
        0.025 * yrange,
      yend = bracket_y,
      linewidth = 0.35
    ) +
    
    annotate(
      "text",
      x = 1.5,
      y = text_y,
      label = p_label,
      family = "Arial",
      size = 7 / ggplot2::.pt
    ) +
    
    scale_fill_manual(
      values = group_colors
    ) +
    
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0.05,
          0.23
        )
      )
    ) +
    
    labs(
      x = NULL,
      y = "ssGSEA score",
      title = paste0(
        "All brain regions pooled (",
        clu,
        ")"
      )
    ) +
    
    theme_classic(
      base_size = 7,
      base_family = "Arial"
    ) +
    
    theme(
      
      text = element_text(
        family = "Arial",
        size = 7,
        color = "black",
        face = "plain"
      ),
      
      plot.title = element_text(
        family = "Arial",
        size = 7,
        color = "black",
        face = "plain",
        hjust = 0
      ),
      
      axis.title = element_text(
        family = "Arial",
        size = 7,
        color = "black",
        face = "plain"
      ),
      
      axis.text = element_text(
        family = "Arial",
        size = 7,
        color = "black",
        face = "plain"
      ),
      
      legend.position = "none"
    )
  
  
  ##########################################################
  # Save PDF
  ##########################################################
  
  ggsave(
    filename = file.path(
      base_dir,
      paste0(
        "GSE48350_",
        clu,
        "_AllBrainRegions_AD_vs_Control.pdf"
      )
    ),
    plot = p,
    width = 2.2,
    height = 2.4,
    device = cairo_pdf
  )
  
  print(
    p
  )
  
}


############################################################
# 24. 单独查看 Cluster3
############################################################

cat(
  "\nCluster3 result:\n"
)

print(
  stat_df[
    stat_df$Cluster ==
      "Cluster3",
    ,
    drop = FALSE
  ]
)


############################################################
# END
############################################################















