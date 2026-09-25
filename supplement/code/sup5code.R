#fig1
setwd('E:/deskup/毕业论文/结直肠癌/data/')
load('figd_suerat.Rdata')
# ==============================================================================
# 加载必要的包
# ==============================================================================
library(Seurat)
library(slingshot)
library(ggplot2)
library(dplyr)
library(ggrastr)
library(scales)
library(patchwork)

# 强制注册 Arial 字体
windowsFonts(Arial = windowsFont("Arial"))

# ==============================================================================
# 0. 环境设置与数据加载
# ==============================================================================
setwd('E:/deskup/毕业论文/结直肠癌/data/')
load('figd_suerat.Rdata')
setwd('E:/deskup/毕业论文/结直肠癌/figs/figs/')

# 定义全局 7pt 主题
font_style_7pt <- element_text(size = 7, color = "black", face = "plain", family = "Arial")

my_theme_7pt <- theme(
  text = font_style_7pt,
  plot.title = element_text(size = 7, color = "black", face = "plain", 
                            family = "Arial", hjust = 0.5),
  legend.text = font_style_7pt,
  legend.title = font_style_7pt,
  legend.key.size = unit(0.3, "cm"),
  plot.margin = margin(2, 2, 2, 2)
)

# ==============================================================================
# 1. Slingshot 轨迹推断与绘图 (proCSC -> CSC -> revCSC)
# ==============================================================================
cat("⏳ 正在运行 Slingshot 轨迹推断并生成轨迹图...\n")

target_linear <- c("proCSC", "CSC", "revCSC")
obj_linear <- subset(seurat_obj, curatedCLUST %in% target_linear)
obj_linear$curatedCLUST <- droplevels(obj_linear$curatedCLUST)

phate_coords_linear <- seurat_obj@reductions[["phate_epi"]]@cell.embeddings
sub_cells <- colnames(obj_linear)

df_linear <- data.frame(
  PHATE_1 = phate_coords_linear[sub_cells, 1],
  PHATE_2 = phate_coords_linear[sub_cells, 2],
  Cluster = as.character(obj_linear$curatedCLUST)
)

set.seed(888)
sling_linear <- slingshot(
  data = as.matrix(df_linear[, 1:2]),
  clusterLabels = df_linear$Cluster,
  start.clus = "proCSC",
  end.clus = "revCSC",
  extend = 'n'
)

curve_data <- as.data.frame(slingCurves(sling_linear)[[1]]$s)
colnames(curve_data) <- c("PHATE_1", "PHATE_2")

new_cluster_colors <- c(
  "CSC" = "#F39B7FFF",
  "proCSC" = "#C77EB5",
  "revCSC" = "#F391A9"
)

p_traj <- ggplot(df_linear, aes(x = PHATE_1, y = PHATE_2)) +
  ggrastr::rasterise(geom_point(aes(color = Cluster), size = 1.5, alpha = 0.8), dpi = 600) +
  geom_path(data = curve_data, aes(x = PHATE_1, y = PHATE_2), 
            color = "black", linewidth = 1, 
            arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
  scale_color_manual(values = new_cluster_colors) +
  theme_void() +
  theme(
    text = font_style_7pt,
    legend.text = font_style_7pt,
    legend.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

cairo_pdf("Slingshot_Trajectory_UpdatedColors_7pt.pdf", width = 5, height = 4, family = "Arial")
print(p_traj)
dev.off()

cat("✔ 轨迹图生成完毕！\n\n")

# ==============================================================================
# 2. 生成 25 个 BayesCor1 模块的 5x5 网格图 (仅保留 Stem 细胞群，宽:高=2:1)
# ==============================================================================
cat("⏳ 正在提取 Stem 细胞 (proCSC, CSC, revCSC) 并生成 5x5 模块打分图...\n")

seurat_clean <- subset(seurat_obj, curatedCLUST %in% target_linear)

phate_coords_clean <- Embeddings(seurat_clean, "phate_epi")

df_plot <- data.frame(
  dim_1 = phate_coords_clean[, 1],
  dim_2 = phate_coords_clean[, 2]
)

module_cols <- paste0("BayesCor1_Module", 1:25, "_UCell")
df_plot <- cbind(df_plot, seurat_clean@meta.data[, module_cols])

global_max <- max(df_plot[, module_cols], na.rm = TRUE)
threshold <- 0.05
my_breaks <- c(0.1, 0.2, 0.3, 0.4, 0.5)

plot_list <- list()

for (i in 1:25) {
  mod_name <- module_cols[i]
  clean_title <- paste0("Module ", i)
  
  df_sub <- df_plot[order(df_plot[[mod_name]]), ]
  
  p <- ggplot(df_sub, aes(x = dim_1, y = dim_2, color = .data[[mod_name]])) +
    rasterise(geom_point(size = 0.28, alpha = 0.9), dpi = 600) +
    scale_color_gradientn(
      colors = c("#d3d7d4", "#ed1941"),
      limits = c(threshold, global_max),
      breaks = my_breaks,
      na.value = "#d3d7d4",
      oob = scales::squish,
      name = "UCell Score"
    ) +
    theme_void() +
    my_theme_7pt +
    ggtitle(clean_title)
  
  plot_list[[i]] <- p
}

combined_5x5 <- wrap_plots(plot_list, ncol = 5) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

cairo_pdf("FeaturePlots_BayesCor1_25Modules_5x5_7pt_StemOnly_2to1.pdf", 
          width = 16, height = 8, family = "Arial")
print(combined_5x5)
dev.off()

cat("✅ 完美搞定！底图已经解除了长宽比锁定，成功拉伸。\n")


#
setwd('E:/deskup/毕业论文/结直肠癌/data/')
load('Human_icms.Rdata')

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggrastr)

# 强制注册 Arial 字体
windowsFonts(Arial = windowsFont("Arial"))

# ==========================================
# 1. 设定分类并计算所有组的差异基因 (DEGs)
# ==========================================
Idents(scRNA) <- "iCMS_msi"
cat("⏳ 正在计算所有 Cluster 的差异表达基因，请稍候...\n")
degs <- FindAllMarkers(scRNA, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)

# ==========================================
# 2. 数据过滤与极性分类
# ==========================================
sig_degs <- degs %>%
  filter(p_val_adj < 0.05) %>%
  mutate(type = ifelse(avg_log2FC > 0, "sigUp", "sigDown"))

# 获取准确的 Cluster 名称及其数量
clusters <- levels(scRNA@meta.data$iCMS_msi)
if (is.null(clusters)) clusters <- unique(as.character(scRNA@meta.data$iCMS_msi))
cluster_mapping <- setNames(1:length(clusters), clusters)

# ==========================================
# 3. 为基因分配 X 轴坐标 (基础坐标 + 随机抖动)
# ==========================================
set.seed(42)
sig_degs$x_base <- cluster_mapping[as.character(sig_degs$cluster)]
# 在每个 Cluster 的中心点左右 0.4 的范围内随机抖动
sig_degs$x_jitter <- sig_degs$x_base + runif(nrow(sig_degs), -0.4, 0.4)

# ==========================================
# 4. 提取每个 Cluster 的 Top 基因用于打标签
# ==========================================
# 默认提取每个 Cluster 上调和下调各自 Top 5 的基因
top_genes <- sig_degs %>%
  group_by(cluster, type) %>%
  top_n(n = 5, wt = abs(avg_log2FC)) %>%
  ungroup()

# ==========================================
# 5. 构建中间彩色标签条块和背景框的数据
# ==========================================
# 动态调整中间条块的高度 (取最大 Log2FC 绝对值的 10% 左右)
y_range <- max(abs(sig_degs$avg_log2FC)) * 0.12 

bar_df <- data.frame(
  cluster = clusters,
  x_center = 1:length(clusters),
  xmin = (1:length(clusters)) - 0.5,
  xmax = (1:length(clusters)) + 0.5,
  ymin = -y_range,
  ymax = y_range
)

# ==========================================
# 6. 生成多组火山图 (Multi-group Volcano)
# ==========================================
# 严格锁定字体为 Arial 7pt
font_style <- element_text(size = 7, color = "black", face = "plain", family = "Arial")

# 设定的细胞类型自定义颜色
custom_colors <- c("iCMS2_MSS" = "#c88400", "iCMS3_MSI-H" = "#00a6ac", "iCMS3_MSS" = "#009ad6")

p_volcano <- ggplot() +
  # A. 添加浅灰色背景柱状阴影区
  geom_rect(data = bar_df, aes(xmin = xmin + 0.05, xmax = xmax - 0.05, ymin = -Inf, ymax = Inf), 
            fill = "grey95", alpha = 0.7) +
  
  # B. 绘制散点 (栅格化处理，600 DPI，彻底杜绝矢量图卡顿，点大小适配缩小至 0.5)
  ggrastr::rasterise(
    geom_point(data = sig_degs, aes(x = x_jitter, y = avg_log2FC, color = type), size = 0.5, alpha = 0.8), 
    dpi = 600
  ) +
  
  # C. 绘制中间的彩色条块区 (linewidth 适当减细适配 7pt 排版)
  geom_rect(data = bar_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster), 
            color = "black", linewidth = 0.3) +
  
  # D. 绘制条块中的 Cluster 文字标签 (Arial 7pt)
  geom_text(data = bar_df, aes(x = x_center, y = 0, label = cluster), 
            family = "Arial", size = 7 / .pt, color = "black") +
  
  # E. 绘制 Top 基因名称的排斥标签 (Arial 7pt)
  geom_text_repel(data = top_genes, aes(x = x_jitter, y = avg_log2FC, label = gene),
                  family = "Arial", size = 7 / .pt, color = "black",
                  box.padding = 0.5, max.overlaps = Inf,
                  segment.color = "black", segment.size = 0.2) +
  
  # F. 颜色与标尺设定
  scale_color_manual(values = c("sigUp" = "#d71345", "sigDown" = "#009ad6")) +
  scale_fill_manual(values = custom_colors) + # 应用指定的模块颜色
  scale_x_continuous(breaks = bar_df$x_center, labels = bar_df$cluster) +
  
  # G. 坐标轴与主题打磨
  theme_classic() +
  labs(x = "Clusters", y = "Average log2FoldChange", color = "type") +
  theme(
    text = font_style,
    axis.text.x = element_blank(), # 隐藏底部原生的 X 轴文字，用中间条块代替
    axis.ticks.x = element_blank(),
    axis.text.y = font_style,
    axis.title = font_style,
    legend.position = "right",
    legend.text = font_style,
    legend.title = font_style,
    legend.key.size = unit(0.3, "cm") # 图例色块缩小适配 7pt 字体
  ) +
  guides(fill = "none") # 隐藏不必要的中间彩色条块的图例

# ==========================================
# 7. 导出为高精度 PDF
# ==========================================
# 画布尺寸压缩，匹配 7pt 字体，防止留白过度
cairo_pdf("MultiGroup_Volcano_Plot_iCMS.pdf", width = 7, height = 5, family = "Arial")
print(p_volcano)
dev.off()

cat("✅ 多组级联火山图生成完毕！\n")



library(dplyr)
library(UCell)

cat("⏳ 1. 正在提取每个 Cluster 的 Top 50 上调基因...\n")

# ==========================================
# 1. 提取 Top 50 差异基因并转换为 List
# ==========================================
top50_genes_list <- degs %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE) %>%
  ungroup() %>%
  # 先提取基因列和分组列，再用 split 转换
  { split(.$gene, .$cluster) }

# 打印检查每个 list 的基因数量，确认是否都取到了 50 个
cat("提取的基因集数量检查：\n")
print(sapply(top50_genes_list, length))

cat("\n⏳ 2. 正在启动 UCell 多线程打分引擎 (ncores = 24)...\n")

# ==========================================
# 2. 运行 UCell 打分
# ==========================================
scRNA <- UCell::AddModuleScore_UCell(
  obj = scRNA,
  features = top50_genes_list,
  ncores = 24,
  name = "_UCell" # UCell 会自动在末尾拼接这个后缀
)

cat("✅ UCell 评分计算完毕！\n")

# ==========================================
# 3. 查看新增的元数据列
# ==========================================
# 新增的列名将是类似 "iCMS2_MSS_UCell", "iCMS3_MSI-H_UCell"
new_ucell_cols <- paste0(names(top50_genes_list), "_UCell")
print("新增的 UCell 打分列如下：")
print(new_ucell_cols)

# 简单看一眼前几行数据
head(scRNA@meta.data[, new_ucell_cols])



library(Seurat)
library(ggplot2)
library(ggrastr)

# 强制注册 Arial 字体
windowsFonts(Arial = windowsFont("Arial"))

# ==========================================
# 1. 统一定义 tSNE 图的主题 (Arial, 14pt, 纯黑, 不加粗)
# ==========================================
my_theme_tsne <- theme(
  text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
  plot.title = element_text(size = 14, color = "black", face = "plain", family = "Arial", hjust = 0.5),
  legend.text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
  legend.title = element_text(size = 14, color = "black", face = "plain", family = "Arial")
)

# 要画的三个 UCell 评分列名 (如果 UCell 自动把 - 改成了 . ，请检查列名)
features_to_plot <- c("iCMS2_MSS_UCell", "iCMS3_MSI-H_UCell", "iCMS3_MSS_UCell")

# ==========================================
# 2. 循环绘制并输出 PDF
# ==========================================
cat("⏳ 正在重新生成高对比度 tSNE 投影图...\n")

for (feat in features_to_plot) {
  
  # 绘制 FeaturePlot
  p <- FeaturePlot(
    scRNA,
    reduction = "tsne",
    features = feat,
    pt.size = 0.5,
    order = TRUE,           # 🌟 核心修复 1：强制把红色的高分细胞画在最上面，防止被灰色掩盖！
    min.cutoff = "q1",      # 🌟 核心修复 2：去掉底部 1% 的极低值干扰
    max.cutoff = "q99",     # 🌟 核心修复 3：去掉顶部 1% 的极高值干扰，让红色充分显现
    raster = TRUE,          # 栅格化防止卡顿
    raster.dpi = c(600, 600)
  ) +
    scale_color_gradientn(
      colors = c("#d3d7d4", "#ed1941"),
      na.value = "#d3d7d4",
      name = "UCell Score"
    ) +
    NoAxes() +
    my_theme_tsne +
    ggtitle(paste0(feat, " (tSNE)"))
  
  # 拼装输出文件名 (替换掉可能引起文件系统问题的特殊符号)
  safe_feat_name <- gsub("-", "_", feat)
  out_pdf <- paste0("tsne_scRNA_", safe_feat_name, ".pdf")
  
  cairo_pdf(out_pdf, width = 7, height = 6, family = "Arial")
  print(p)
  dev.off()
  
  cat("✔ 成功生成:", out_pdf, "\n")
}

cat("✅ 全部 tSNE 降维图绘制完毕！快去看看红点是不是出来了。\n")

















library(Seurat)
library(ggplot2)
library(dplyr)

setwd('E:/deskup/毕业论文/结直肠癌/data/')
load('to_ucell.Rdata')
############################
# 1. 提取 kME > 0.6 gene
############################

module_gene_ucell <- lapply(
  module_kME_list,
  function(x){
    
    x %>%
      filter(kME > 0.6) %>%
      pull(gene)
    
  }
)

# 查看每个module多少gene
sapply(
  module_gene_ucell,
  length
)
############################
# 2. 去除过小module
############################
module_gene_ucell <- module_gene_ucell[
  sapply(module_gene_ucell,length) >= 0
]
length(module_gene_ucell)











# ==========================================
# 0. Load Required Packages
# ==========================================
library(Seurat)
library(UCell)
library(homologene)
module_gene_ucell<-str(filtered_signatures)
print("Step 1: Converting Mouse gene modules to Human gene modules...")

# Initialize an empty list to store the converted human genes
human_modules <- list()


names(module_gene_ucell)<-1:25
# Loop through each module in the mouse gene list
for (mod_name in names(module_gene_ucell)) {
  
  mouse_genes <- module_gene_ucell[[mod_name]]
  
  # Convert Mouse (10090) to Human (9606) using homologene
  trans_df <- homologene(mouse_genes, inTax = 10090, outTax = 9606)
  
  # Extract the unique human genes mapped by the package
  h_genes <- unique(trans_df$`9606`)
  
  # Add to the new list only if mapping was successful for at least one gene
  if (length(h_genes) > 0) {
    human_modules[[mod_name]] <- h_genes
  }
}

print("Conversion complete. Showing the first few human genes of Module_1:")
print(head(human_modules[["Module_1"]]))

print("Step 2: Calculating UCell scores on human scRNA object with 24 cores...")









mouse_modules<-filtered_signatures

# 加载必要的包
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(dplyr)

# 假设你的两个 list 已经分别命名为 human_modules 和 mouse_modules (即 filtered_signatures)

# ==========================================
# 1. 人类基因集富集 (Human GOBP)
# ==========================================
human_enrich <- compareCluster(
  geneClusters = human_modules, 
  fun = "enrichGO", 
  OrgDb = org.Hs.eg.db, 
  keyType = "SYMBOL", 
  ont = "BP", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05
)

# 转化为 dataframe 方便后续操作
human_res <- as.data.frame(human_enrich)

# ==========================================
# 2. 小鼠基因集富集 (Mouse GOBP)
# ==========================================
mouse_enrich <- compareCluster(
  geneClusters = mouse_modules, 
  fun = "enrichGO", 
  OrgDb = org.Mm.eg.db, 
  keyType = "SYMBOL", 
  ont = "BP", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05
)

mouse_res <- as.data.frame(mouse_enrich)

# ==========================================
# 3. 提取双物种保守的 GO Terms (Intersection)
# ==========================================
# 统一列名以方便合并，主要依据 Cluster (模块名) 和 ID (GO Term ID)
conserved_go <- inner_join(
  human_res, 
  mouse_res, 
  by = c("Cluster", "ID"),
  suffix = c("_Human", "_Mouse")
)

# 查看双物种共同显著富集的通路
head(conserved_go[, c("Cluster", "ID", "Description_Human", "p.adjust_Human", "p.adjust_Mouse")])

# 可选：保存结果
write.csv(conserved_go, "Conserved_GOBP_Modules.csv", row.names = FALSE)














# ==========================================
# 2. Run UCell Scoring on Human scRNA Data
# ==========================================
# AddModuleScore_UCell automatically appends the string provided in 'name'
scRNA <- AddModuleScore_UCell(
  obj = scRNA,
  features = human_modules,
  ncores = 24,
  name = "_UCell"
)



library(Seurat)
library(UCell)
library(corrplot)

# 强制注册 Arial 字体
windowsFonts(Arial = windowsFont("Arial"))

# ==========================================
# 1. 暴力清理旧数据，防止任何残留覆盖问题
# ==========================================
cat("⏳ 正在清理旧的 Module 打分...\n")

cols_to_remove <- grep("^Module_.*_UCell$|^M[0-9]+_UCell$", colnames(scRNA@meta.data), value = TRUE)

if(length(cols_to_remove) > 0) {
  scRNA@meta.data[, cols_to_remove] <- NULL
}

# ==========================================
# 2. 从源头重命名基因集并重新计算 UCell
# ==========================================
# 直接把原本编号跳跃的 human_modules 强制重命名为连续的 M1 到 M25
names(human_modules) <- paste0("M", 1:25)

cat("⏳ 正在启动 UCell 多线程引擎重新打分...\n")

scRNA <- UCell::AddModuleScore_UCell(
  obj = scRNA,
  features = human_modules,
  ncores = 24,
  name = "_UCell"          # UCell 会自动加上后缀，所以新列名会变成 M1_UCell 等
)

# ==========================================
# 3. 提取新打分进行 Scale
# ==========================================
new_modules <- paste0("M", 1:25, "_UCell")
icms_cols <- c("iCMS2_MSS_UCell", "iCMS3_MSI-H_UCell", "iCMS3_MSS_UCell")

mat_modules <- scale(as.matrix(scRNA@meta.data[, new_modules]))
mat_icms <- scale(as.matrix(scRNA@meta.data[, icms_cols]))

# ==========================================
# 4. 计算相关性矩阵与 P 值
# ==========================================
cat("⏳ 正在计算 3x25 相关性矩阵...\n")

cor_sub <- cor(mat_icms, mat_modules, method = "pearson", use = "pairwise.complete.obs")

p_sub <- matrix(1, nrow = ncol(mat_icms), ncol = ncol(mat_modules))

for(i in 1:ncol(mat_icms)) {
  for(j in 1:ncol(mat_modules)) {
    try({
      test <- cor.test(mat_icms[, i], mat_modules[, j], method = "pearson")
      p_sub[i, j] <- test$p.value
    }, silent = TRUE)
  }
}

# ==========================================
# 5. 格式化图表标签 (洗掉 _UCell，干干净净)
# ==========================================
rownames(cor_sub) <- c("iCMS2_MSS", "iCMS3_MSI-H", "iCMS3_MSS")
colnames(cor_sub) <- paste0("M", 1:25)
rownames(p_sub) <- rownames(cor_sub)
colnames(p_sub) <- colnames(cor_sub)

# 极端容错：兜底可能的 NA 值
cor_sub[is.na(cor_sub)] <- 0
p_sub[is.na(p_sub)] <- 1

# ==========================================
# 6. 画图输出 (Arial 7pt, 黑体不加粗)
# ==========================================
my_palette <- colorRampPalette(c("#4393C3", "#FFFFFF", "#D6604D"))(200)

cairo_pdf("Figure_Cor_Heatmap_3x25_Circles.pdf", width = 10, height = 2.5, family = "Arial")

par(ps = 7, family = "Arial", font = 1, col.axis = "black")

corrplot::corrplot(
  cor_sub,
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
  p.mat = p_sub,
  sig.level = 0.05,
  insig = "blank",
  mar = c(1, 1, 1, 1)
)

dev.off()

cat("✅ 重构完成！热图已保存。\n")



#fig5
# ==========================================
# 1. 定义增强版 GMT 读取函数
# ==========================================
read_gmt_base <- function(file) {
  lines <- readLines(file, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]  # 去除注释行
  
  gene_list <- lapply(lines, function(x) {
    parts <- strsplit(x, "\t")[[1]]
    if (length(parts) < 3) return(character(0))
    genes <- parts[3:length(parts)]
    genes[genes != ""]
  })
  
  names(gene_list) <- sapply(lines, function(x) {
    parts <- strsplit(x, "\t")[[1]]
    parts[1]
  })
  
  return(gene_list)
}

# ==========================================
# 2. 读取 GMT 文件
# ==========================================
cat("⏳ 正在读取 MSigDB GMT 文件...\n")

ctnnb1_dn_gmt <- read_gmt_base("FEVR_CTNNB1_TARGETS_DN.v2026.1.Hs.gmt")
ctnnb1_up_gmt <- read_gmt_base("FEVR_CTNNB1_TARGETS_UP.v2026.1.Hs.gmt")

# 查看有哪些基因集
cat("DN 文件包含的基因集:\n")
print(names(ctnnb1_dn_gmt))
cat("UP 文件包含的基因集:\n")
print(names(ctnnb1_up_gmt))

# 取第一个基因集（如果是多个，可以按名字取）
genes_dn <- ctnnb1_dn_gmt[[1]]
genes_up <- ctnnb1_up_gmt[[1]]


genes_m2<-human_modules$Module_2

# ==========================================
# 4. 计算交集
# ==========================================
overlap_dn <- intersect(genes_m2, genes_dn)
overlap_up <- intersect(genes_m2, genes_up)

# ==========================================
# 5. 打印统计结果
# ==========================================
cat("\n================ 统计结果 ================\n")
cat(sprintf("M2 模块基因总数: %d\n", length(genes_m2)))
cat(sprintf("CTNNB1_TARGETS_DN 基因总数: %d\n", length(genes_dn)))
cat(sprintf("CTNNB1_TARGETS_UP 基因总数: %d\n", length(genes_up)))
cat("------------------------------------------\n")
cat(sprintf("🔥 M2 ∩ CTNNB1_UP: %d 个基因\n", length(overlap_up)))
if (length(overlap_up) > 0) {
  cat("   ", paste(overlap_up, collapse = ", "), "\n")
}
cat(sprintf("❄️ M2 ∩ CTNNB1_DN: %d 个基因\n", length(overlap_dn)))
if (length(overlap_dn) > 0) {
  cat("   ", paste(overlap_dn, collapse = ", "), "\n")
}
cat("==========================================\n")







# ==============================================================================
# 脚本: 生成 M2 与 CTNNB1_TARGETS_UP/DN 的 Venn 图
# 功能: 绘制两个 Venn 图，分别展示 M2 与 CTNNB1_DN（蓝色）和
#       M2 与 CTNNB1_UP（红色）的基因交集
# 输出: Venn_M2_vs_CTNNB1_DN_Blue.pdf
#       Venn_M2_vs_CTNNB1_UP_Red.pdf
# ==============================================================================

# 1. 加载必要包 ----------------------------------------------------------------
if (!require("ggvenn")) {
  install.packages("ggvenn")
  library(ggvenn)
} else {
  library(ggvenn)
}

if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
} else {
  library(ggplot2)
}

# 2. 字体注册（跨平台兼容）-----------------------------------------------------
if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
} else {
  tryCatch({
    systemfonts::register_font("Arial", plain = "Arial")
  }, error = function(e) {
    cat("⚠️ 未找到 Arial 字体，将使用系统默认无衬线字体\n")
  })
}

# 3. 检查输入数据是否存在 ------------------------------------------------------
if (!exists("genes_m2")) {
  stop("❌ 错误: 变量 'genes_m2' 不存在，请先加载 M2 基因集！")
}
if (!exists("genes_dn")) {
  stop("❌ 错误: 变量 'genes_dn' 不存在，请先加载 CTNNB1_TARGETS_DN 基因集！")
}
if (!exists("genes_up")) {
  stop("❌ 错误: 变量 'genes_up' 不存在，请先加载 CTNNB1_TARGETS_UP 基因集！")
}

cat("✅ 输入数据检查通过\n")
cat("   M2 基因数:", length(genes_m2), "\n")
cat("   CTNNB1_DN 基因数:", length(genes_dn), "\n")
cat("   CTNNB1_UP 基因数:", length(genes_up), "\n\n")

# 4. 统一定义全局字体样式 (14pt, 黑色, 不加粗, Arial) --------------------------
font_style <- element_text(
  size = 14, 
  color = "black", 
  face = "plain", 
  family = "Arial"
)

# ==========================================
# 5. M2 vs CTNNB1_DN (M2 绿色，DN 蓝色)
# ==========================================
list_dn <- list(
  "Module 2" = genes_m2,
  "CTNNB1_TARGETS_DN" = genes_dn
)

p_dn <- ggvenn(
  list_dn,
  fill_color = c("#4DAF4A", "#E41A1C"),  # M2: 绿色, DN: 蓝色
  stroke_size = 0.5,
  stroke_color = "black",
  set_name_size = 14 / .pt,              # 转换为 ggplot 适用的 pt 尺寸
  text_size = 14 / .pt,
  text_color = "black"
) +
  theme_void() +
  theme(text = font_style)

# 导出为 PDF
cairo_pdf("Venn_M2_vs_CTNNB1_DN_Blue.pdf", width = 5, height = 5, family = "Arial")
print(p_dn)
dev.off()

cat("✅ 第一个 Venn 图 (蓝色主题) 已生成。\n")

# ==========================================
# 6. M2 vs CTNNB1_UP (M2 绿色，UP 红色)
# ==========================================
list_up <- list(
  "Module 2" = genes_m2,
  "CTNNB1_TARGETS_UP" = genes_up
)

p_up <- ggvenn(
  list_up,
  fill_color = c("#4DAF4A", "#377EB8"),  # M2: 绿色, UP: 红色
  stroke_size = 0.5,
  stroke_color = "black",
  set_name_size = 14 / .pt,
  text_size = 14 / .pt,
  text_color = "black"
) +
  theme_void() +
  theme(text = font_style)

# 导出为 PDF
cairo_pdf("Venn_M2_vs_CTNNB1_UP_Red.pdf", width = 5, height = 5, family = "Arial")
print(p_up)
dev.off()

cat("✅ 第二个 Venn 图 (红色主题) 已生成。\n")

# ==========================================
# 7. 完成提示
# ==========================================
cat("\n============================================================\n")
cat("🎉 全部 Venn 图生成完毕！\n")
cat("文件列表:\n")
cat("   1. Venn_M2_vs_CTNNB1_DN_Blue.pdf\n")
cat("   2. Venn_M2_vs_CTNNB1_UP_Red.pdf\n")
cat("保存路径:", getwd(), "\n")
cat("============================================================\n")


