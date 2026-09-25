# 强制将 Windows 系统的 Arial 字体注册到 R 的数据库中
windowsFonts(Arial = windowsFont("Arial"))

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrastr)
library(mascarade)
library(Mescnet)
setwd('E:/deskup/Mescnet/data/epi/')

# ==============================================================================
# 0. 数据加载与预处理
# ==============================================================================
seurat_obj <- readRDS('INTepi.rds')

# 剔除 Outlier 噪点细胞
seurat_obj <- subset(seurat_obj, subset = curatedCLUST %in% c("Outlier 1 (Stressed?)", "Outlier 2 (CRC?)"), invert = TRUE)


seurat_obj$cell_type <- seurat_obj$curatedCLUST


#fig1
# ==============================================================================
# 1. 提取 PHATE 坐标
# ==============================================================================
phate_coords <- Embeddings(seurat_clean, "phate_epi")

df <- data.frame(
  dim_1 = phate_coords[, 1],
  dim_2 = phate_coords[, 2],
  cell_type = seurat_clean$cell_type
)

# 保持与之前完全一致的 Cell order
cell_order <- c(
  "CSC", "proCSC", "revCSC",
  "TA 1", "TA 2",
  "Early Enterocyte", "Late Enterocyte",
  "Goblet / DCS", "ER Stress"
)
df$cell_type <- factor(df$cell_type, levels = cell_order)

# 保持与之前完全一致的配色
color_use <- c(
  "CSC"              = "#F39B7FFF",
  "proCSC"           = "#C77EB5",
  "revCSC"           = "#F391A9",
  "TA 1"             = "#3C5488FF",
  "TA 2"             = "#4DBBD5FF",
  "Early Enterocyte" = "#00A087FF",
  "Late Enterocyte"  = "#91D1C2FF",
  "Goblet / DCS"     = "#7E6148FF",
  "ER Stress"        = "#8491B4FF"
)

# ==============================================================================
# 2. 重新规划主群 (用于计算虚线边界)
# ==============================================================================
df <- df %>%
  mutate(Major_Lineage = case_when(
    cell_type %in% c("CSC", "proCSC", "revCSC") ~ "Stem Cells",
    cell_type %in% c("TA 1", "TA 2") ~ "TA Cells",
    cell_type %in% c("Early Enterocyte", "Late Enterocyte") ~ "Enterocytes",
    cell_type %in% c("Goblet / DCS") ~ "Secretory Cells",
    cell_type %in% c("ER Stress") ~ "Stressed Cells",
    TRUE ~ "Other"
  ))

# ==============================================================================
# 3. 生成虚线边界 (低阈值、高平滑度)
# ==============================================================================
cat("✅ 正在生成精准虚线边界...\n")
maskTable <- generateMask(
  dims = df[, 1:2],
  cluster = df$Major_Lineage, 
  minDensity = 0.03,   # 降低阈值扩圈
  smoothSigma = 0.18   # 增加平滑度融圈
)

# ==============================================================================
# 4. ggplot2 终极组合出图 (栅格化 + 去除文字与坐标轴)
# ==============================================================================
cat("🎨 正在绘制终极版带圈 PHATE 图...\n")

p_final <- ggplot(df, aes(x = dim_1, y = dim_2)) +
  
  # A. 底层散点：恢复点大小 0.28 和透明度 0.9，并进行光栅化 (600 DPI)
  rasterise(
    geom_point(aes(color = cell_type), size = 0.28, alpha = 0.9), 
    dpi = 600
  ) +
  scale_color_manual(values = color_use) +
  
  # B. 中层：保留谱系外部边界虚线圈
  geom_path(data = maskTable, aes(group = group),
            linewidth = 0.4, linetype = 2, color = "black") +
  
  # C. 极致断舍离主题 (去除内部标签和手工坐标轴)
  coord_fixed(ratio = 1) + 
  theme_void() + 
  theme(
    text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

# 打印图像
print(p_final)

# ==============================================================================
# 5. 规范化文件导出
# ==============================================================================
setwd('E:/deskup/毕业论文/结直肠癌/figs/')

ggsave("PHATE_Epi_Raster_Clean.pdf", plot = p_final, width = 6.5, height = 5, device = cairo_pdf)
ggsave("PHATE_Epi_Raster_Clean.png", plot = p_final, width = 6.5, height = 5, dpi = 600)

cat("🎉 完美！PHATE 图已生成：保留了虚线圈，去除了多余文字和坐标轴，应用了光栅化排版！\n")



#fig2
setwd('E:/deskup/Mescnet/data/epi/benchmark/')
load('zhongjianshuju/merged_data.Rdata')
seurat_obj@misc$MEscnet$datExpr<-merged_data

library(Mescnet)





library(ggplot2)
library(patchwork)

# 确保系统 Arial 字体可用
windowsFonts(Arial = windowsFont("Arial"))
library(ggplot2)
library(patchwork)

# 确保系统 Arial 字体可用
windowsFonts(Arial = windowsFont("Arial"))

pick_power_scale_free_power <- function (seurat_obj, powerVector = 1:20, RsquaredCut = 0.8, nBreaks = 10) 
{
  # --- 1. 矩阵计算部分保持不变 ---
  pcor_matrix <- pcor(seurat_obj@misc$MEscnet$datExpr)$estimate
  diag(pcor_matrix) <- 0
  colnames(pcor_matrix) <- colnames(seurat_obj@misc$MEscnet$datExpr)
  rownames(pcor_matrix) <- colnames(seurat_obj@misc$MEscnet$datExpr)
  seurat_obj@misc$MEscnet$MEscnet_ppcor <- as.matrix(pcor_matrix)
  
  cor_matrix <- cor(seurat_obj@misc$MEscnet$datExpr, method = "pearson")
  diag(cor_matrix) <- 0
  colnames(cor_matrix) <- colnames(seurat_obj@misc$MEscnet$datExpr)
  rownames(cor_matrix) <- colnames(seurat_obj@misc$MEscnet$datExpr)
  
  bayes_cor1 <- abs(pcor_matrix) * abs(cor_matrix)
  bayes_cor1[is.na(bayes_cor1)] <- 0
  seurat_obj@misc$MEscnet$bayes_cor1 <- as.matrix(bayes_cor1)
  
  bayes_cor2 <- pmax(pcor_matrix, 0) * pmax(cor_matrix, 0)
  bayes_cor2[is.na(bayes_cor2)] <- 0
  seurat_obj@misc$MEscnet$bayes_cor2 <- as.matrix(bayes_cor2)
  
  cor_norm <- (cor_matrix - min(cor_matrix, na.rm = TRUE))/(max(cor_matrix, na.rm = TRUE) - min(cor_matrix, na.rm = TRUE))
  pcor_norm <- (pcor_matrix - min(pcor_matrix, na.rm = TRUE))/(max(pcor_matrix, na.rm = TRUE) - min(pcor_matrix, na.rm = TRUE))
  bayes_cor3 <- pcor_norm * cor_norm
  bayes_cor3[is.na(bayes_cor3)] <- 0
  seurat_obj@misc$MEscnet$bayes_cor3 <- as.matrix(bayes_cor3)
  
  cor_norm_alt <- (1 + cor_matrix)/2
  pcor_norm_alt <- (1 + pcor_matrix)/2
  bayes_cor4 <- pcor_norm_alt * cor_norm_alt
  bayes_cor4[is.na(bayes_cor4)] <- 0
  seurat_obj@misc$MEscnet$bayes_cor4 <- as.matrix(bayes_cor4)
  
  cor_z <- (cor_matrix - mean(cor_matrix, na.rm = TRUE))/sd(cor_matrix, na.rm = TRUE)
  pcor_z <- (pcor_matrix - mean(pcor_matrix, na.rm = TRUE))/sd(pcor_matrix, na.rm = TRUE)
  bayes_cor5 <- pcor_z * cor_z
  bayes_cor5[is.na(bayes_cor5)] <- 0
  seurat_obj@misc$MEscnet$bayes_cor5 <- as.matrix(bayes_cor5)
  
  cor_norm2 <- (1 + pmax(pcor_matrix, 0))/2
  bayes_cor6 <- cor_norm2 * pmax(cor_matrix, 0)
  bayes_cor6[is.na(bayes_cor6)] <- 0
  seurat_obj@misc$MEscnet$bayes_cor6 <- as.matrix(bayes_cor6)
  
  cor_matrices <- list(pcor = pcor_matrix, bayes_cor1 = bayes_cor1, 
                       bayes_cor2 = bayes_cor2, bayes_cor3 = bayes_cor3, 
                       bayes_cor4 = bayes_cor4, bayes_cor5 = bayes_cor5, 
                       bayes_cor6 = bayes_cor6)
  
  result_tables <- list()
  recommended_powers <- list()
  adj_matrices <- list()
  
  # 🌟 统一学术排版主题：全局 Arial, 7pt, 纯黑, 不加粗
  nature_theme <- theme_classic() +
    theme(
      text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
      axis.text = element_text(size = 7, color = "black", face = "plain", family = "Arial"),
      plot.title = element_text(size = 7, color = "black", face = "plain", family = "Arial", hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  for (name in names(cor_matrices)) {
    result <- data.frame(Power = powerVector, Rsquared = NA, Slope = NA, MeanK = NA, MedianK = NA, MaxK = NA)
    for (i in seq_along(powerVector)) {
      power <- powerVector[i]
      adj <- abs(cor_matrices[[name]])^power
      diag(adj) <- 0
      k <- rowSums(adj, na.rm = TRUE)
      hist_k <- hist(k, breaks = nBreaks, plot = FALSE)
      counts <- hist_k$counts
      mids <- hist_k$mids
      valid <- (mids > 0) & (counts > 0)
      if (length(mids[valid]) > 1) {
        log_k <- log10(mids[valid])
        log_p <- log10(counts[valid]/sum(counts[valid]))
        fit <- lm(log_p ~ log_k)
        result$Slope[i] <- coef(fit)[2]
        result$Rsquared[i] <- summary(fit)$r.squared
      }
      result$MeanK[i] <- mean(k, na.rm = TRUE)
      result$MedianK[i] <- median(k, na.rm = TRUE)
      result$MaxK[i] <- max(k, na.rm = TRUE)
    }
    result_tables[[name]] <- result
    idx_r2 <- which(result$Rsquared > RsquaredCut)[1]
    power_selected <- ifelse(is.na(idx_r2), NA, result$Power[idx_r2])
    recommended_powers[[name]] <- power_selected
    
    if (!is.na(power_selected)) {
      adj <- abs(cor_matrices[[name]])^power_selected
      diag(adj) <- 0
      adj_matrices[[name]] <- adj
    }
    
    # --- 2. 升级为 ggplot2 绘图 ---
    colors <- c("#00a6ac", "#45b97c", "#009ad6", "#769149")
    
    # 图 1: R-squared
    p1 <- ggplot(result, aes(x = Power, y = Rsquared)) +
      geom_point(size = 1.5, color = colors[1]) +
      geom_hline(yintercept = RsquaredCut, color = "red", linetype = "dashed") +
      labs(x = "Soft Threshold (power)", y = "Scale Free Topology Fit", title = paste(name, "R-squared vs Power")) +
      nature_theme
    
    if (!is.na(power_selected)) {
      p1 <- p1 + 
        geom_point(data = subset(result, Power == power_selected), color = "red", size = 2.5) +
        annotate("text", x = power_selected, y = result$Rsquared[idx_r2], 
                 label = paste0("Power=", power_selected), vjust = -1.5, color = "red", 
                 size = 7 / .pt, family = "Arial", fontface = "plain") # 文本同步改为 7pt
    }
    
    # 图 2: Mean Connectivity
    p2 <- ggplot(result, aes(x = Power, y = MeanK)) +
      geom_point(size = 1.5, color = colors[2]) +
      labs(x = "Soft Threshold (power)", y = "Mean Connectivity", title = paste(name, "Mean Connectivity")) +
      nature_theme
    if (!is.na(power_selected)) {
      p2 <- p2 + geom_point(data = subset(result, Power == power_selected), color = "red", size = 2.5)
    }
    
    # 图 3: Median Connectivity
    p3 <- ggplot(result, aes(x = Power, y = MedianK)) +
      geom_point(size = 1.5, color = colors[3]) +
      labs(x = "Soft Threshold (power)", y = "Median Connectivity", title = paste(name, "Median Connectivity")) +
      nature_theme
    if (!is.na(power_selected)) {
      p3 <- p3 + geom_point(data = subset(result, Power == power_selected), color = "red", size = 2.5)
    }
    
    # 图 4: Max Connectivity
    p4 <- ggplot(result, aes(x = Power, y = MaxK)) +
      geom_point(size = 1.5, color = colors[4]) +
      labs(x = "Soft Threshold (power)", y = "Max Connectivity", title = paste(name, "Max Connectivity")) +
      nature_theme
    if (!is.na(power_selected)) {
      p4 <- p4 + geom_point(data = subset(result, Power == power_selected), color = "red", size = 2.5)
    }
    
    # 将4张图横向拼接
    combined_plot <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)
    
    # 🌟 缩小画布尺寸适配 7 号字体，宽度 8 英寸，高度 2 英寸
    pdf_filename <- paste0("ScaleFree_1x4_", name, ".pdf")
    ggsave(filename = pdf_filename, plot = combined_plot, width = 8, height = 2, device = cairo_pdf)
    
    cat("✔", name, "Recommended power =", power_selected, "| Plot saved to", pdf_filename, "\n")
  }
  
  seurat_obj@misc$MEscnet$result_tables <- result_tables
  seurat_obj@misc$MEscnet$recommended_powers <- recommended_powers
  seurat_obj@misc$MEscnet$adj_matrices <- adj_matrices
  seurat_obj@misc$MEscnet$tom <- adj_matrices[["pcor"]]
  
  return(seurat_obj)
}
#sup
seurat_obj<-pick_power_scale_free_power(seurat_obj)
seurat_obj<-ComputeMEscnetModules(seurat_obj,min_genes_per_module = 20,resolution = 2,number = 2)
modules<-seurat_obj@misc$MEscnet$MEscnet_modules$module_ids
module_genes<-seurat_obj@misc$MEscnet$MEscnet_modules$gene_lists








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

################################
# 6. 查看
################################

length(module_kME_list)
names(module_kME_list)
head(module_kME_list[[1]])






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




#figc
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(reshape2)

# 强制注册 Arial 字体，防止报错
windowsFonts(Arial = windowsFont("Arial"))

cat("⏳ 正在提取模块计算相关性，并生成光栅化圈圈热图...\n")

# ==============================================================================
# 1. 自动校对列名 (使用正确的 UCell 格式)
# ==============================================================================
all_cols <- colnames(seurat_obj@meta.data)

my_sigs <- c(
  "Module_1_UCell", "Module_2_UCell", "Module_3_UCell", 
  "Module_4_UCell", "Module_5_UCell", "Module_6_UCell", 
  "Module_7_UCell", "Module_8_UCell", "Module_9_UCell", 
  "Module_11_UCell", "Module_12_UCell", "Module_13_UCell", 
  "Module_14_UCell", "Module_15_UCell", "Module_16_UCell", 
  "Module_17_UCell", "Module_18_UCell", "Module_19_UCell", 
  "Module_20_UCell", "Module_21_UCell", "Module_25_UCell", 
  "Module_28_UCell", "Module_29_UCell", "Module_30_UCell", 
  "Module_32_UCell"
)
my_sigs_exist <- my_sigs[my_sigs %in% all_cols]

if(length(my_sigs_exist) < 25) {
  cat("⚠️ 警告：只找到了", length(my_sigs_exist), "个模块，缺失的将被自动跳过。\n")
}

# 作者最终使用的 Signatures 列表
pub_sigs_all <- c(
  "FetalHAN20_UCell","RevscHAN20_UCell","IscHAN20_UCell","YapHAN20_UCell","WntHAN20_UCell",
  "IscGREGORIEFF15_UCell","RepGREGORIEFF15_UCell","YapGREGORIEFF15_UCell","hEphb2MERLOS11_UCell",
  "hLgr5MERLOS11_UCell","hProlifMERLOS11_UCell","IscMERLOS11_UCell","Lgr5MERLOS11_UCell",
  "ProlifMERLOS11_UCell","FetalYUI18_UCell","Ssc1AYYAZ19_UCell","Ssc2AYYAZ19_UCell",
  "Ssc2aAYYAZ19_UCell","Ssc2bAYYAZ19_UCell","Ssc2cAYYAZ19_UCell","StemBUES22_UCell",
  "RsBUES22_UCell","Notch1MOURAO19_UCell","Lgr5MOURAO19_UCell","IscLI17_UCell",
  "StemcorrLI17_UCell","TumourLI17_UCell","WntLI17_UCell","TgfbLI17_UCell",
  "StemtaLI17_UCell","ProgenstemDALERBA11_UCell","ImmatureDALERBA11_UCell",
  "CancerDALERBA11_UCell","ProlifDALERBA11_UCell","StemtaPELKA21_UCell",
  "StemtasecPELKA21_UCell","StemtaprolifPELKA21_UCell","OWNwntreceptors_UCell",
  "OWNsig_stem_UCell","OWNsig_stemO_UCell","OWNsig_stemS_UCell","OWNsig_prolif_UCell",
  "OWNdeclust_stem_UCell","OWNdeclust_stemO_UCell","OWNdeclust_stemS_UCell",
  "LGR5munoz_ALVAREZ22_UCell","FETALmustata_ALVAREZ22_UCell","MEX3Abarriga_ALVAREZ22_UCell",
  "LGR5short_ALVAREZ22_UCell","MEX3A_ALVAREZ22_UCell","LGR5_ALVAREZ22_UCell",
  "YAPsign_ALVAREZ22_UCell","hEpiHrCANELLAS22_UCell","hWntMORRAL20_UCell",
  "hYapWANG18_UCell","hiCMS2_UCell","hiCMS3_UCell","mRegenscGIL22_UCell",
  "hTgfbKEGG_UCell","hKrasGSEA_UCell","hMapkKEGG_UCell","hMapkGO_UCell","hPi3kGO_UCell"
)

sigs_to_remove <- c(
  "ImmatureDALERBA11_UCell", "OWNwntreceptors_UCell", "hKrasGSEA_UCell",
  "hTgfbKEGG_UCell", "hMapkKEGG_UCell", "hMapkGO_UCell", "hPi3kGO_UCell",
  "hiCMS2_UCell", "hiCMS3_UCell"
)
pub_sigs_final <- setdiff(pub_sigs_all, sigs_to_remove)
pub_sigs_exist <- pub_sigs_final[pub_sigs_final %in% all_cols]

# ==============================================================================
# 2. 提取数据并计算排序
# ==============================================================================
# 锁定干细胞和 TA 细胞
meta_sub <- seurat_obj@meta.data %>%
  filter(curatedCLUST %in% c("CSC", "proCSC", "revCSC", "TA 1", "TA 2"))

# 独立聚类获得内部顺序
hc_pub <- hclust(as.dist(1 - cor(scale(as.matrix(meta_sub[, pub_sigs_exist])), method = "pearson")), method = "complete")
order_pub <- pub_sigs_exist[hc_pub$order]

hc_my <- hclust(as.dist(1 - cor(scale(as.matrix(meta_sub[, my_sigs_exist])), method = "pearson")), method = "complete")
order_my <- my_sigs_exist[hc_my$order]

# ==============================================================================
# 3. 计算相关性与 P 值，并将其转换为 ggplot2 长格式数据
# ==============================================================================
final_order <- c(order_my, order_pub)
mat_all <- scale(as.matrix(meta_sub[, final_order]))
cor_all <- cor(mat_all, method = "pearson")
res_all <- corrplot::cor.mtest(mat_all, method = "pearson")

# 使用 reshape2::melt 转化矩阵
cor_df <- melt(cor_all, varnames = c("Var1", "Var2"), value.name = "Cor")
p_df <- melt(res_all$p, varnames = c("Var1", "Var2"), value.name = "p_value")
cor_df$p_value <- p_df$p_value

# 屏蔽不显著的值
cor_df <- cor_df %>% mutate(
  plot_size = ifelse(p_value < 0.05, abs(Cor), NA),
  plot_fill = ifelse(p_value < 0.05, Cor, NA)
)

# 锁定轴的顺序（Var2 倒序排列）
cor_df$Var1 <- factor(cor_df$Var1, levels = final_order)
cor_df$Var2 <- factor(cor_df$Var2, levels = rev(final_order))

# 🌟 核心修改：换成柔和的学术低饱和度配色
# 柔和钢化蓝：#4393C3 -> 纯白：#FFFFFF -> 柔和玫瑰红：#D6604D

my_palette <- colorRampPalette(c("#4393C3", "#FFFFFF", "#D6604D"))(200)

png('Figure_Cor_Heatmap_Final_Circles.png', width = 28, height = 28, units = "in", res = 600)

corrplot::corrplot(
  cor_all,
  method = "circle",
  order = "original",
  tl.pos = "n",
  cl.pos = "r",
  cl.cex = 0.00001,
  cl.offset = 20,
  col = my_palette,
  p.mat = res_all$p,
  sig.level = 0.05,
  insig = "blank",
  mar = c(1, 1, 1, 1)
)

dev.off()

cat("Done!\n")




colnames(seurat_obj@meta.data)





#figd
library(Seurat)
library(slingshot)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggplot2)
library(ggrastr)

# 强制注册 Arial 字体
windowsFonts(Arial = windowsFont("Arial"))
setwd('E:/deskup/毕业论文/结直肠癌/data/')
load('figd_suerat.Rdata')


setwd('E:/deskup/毕业论文/结直肠癌/figs/figs/')
# ==========================================
# 1. 设定前缀并提取数据
# ==========================================
prefix <- "BayesCor1" 
n_modules <- 25 

module_names <- paste0(prefix, "_Module", 1:n_modules, "_UCell")

target_linear <- c("proCSC", "CSC", "revCSC")
obj_linear <- subset(seurat_obj, curatedCLUST %in% target_linear)
obj_linear$curatedCLUST <- droplevels(obj_linear$curatedCLUST)

df_linear <- data.frame(
  PHATE_1 = Embeddings(obj_linear, "phate_epi")[, 1],
  PHATE_2 = Embeddings(obj_linear, "phate_epi")[, 2],
  Cluster = as.character(obj_linear$curatedCLUST)
)

# ==========================================
# 2. 运行 Slingshot 轨迹推断并绘制轨迹图
# ==========================================
set.seed(888)
sling_linear <- slingshot(
  data = as.matrix(df_linear[, 1:2]), 
  clusterLabels = df_linear$Cluster, 
  start.clus = "proCSC", 
  end.clus = "revCSC", 
  extend = 'n' 
)

pt_values <- slingPseudotime(sling_linear)[,1]
valid_cells <- !is.na(pt_values)
pt_values <- pt_values[valid_cells]

pt_order <- order(pt_values)
pt_values <- pt_values[pt_order]
ordered_cells <- names(pt_values)

# --- 绘制并保存 Slingshot 轨迹线 PDF ---
curve_data <- as.data.frame(slingCurves(sling_linear)[[1]]$s)
colnames(curve_data) <- c("PHATE_1", "PHATE_2")

p_traj <- ggplot(df_linear, aes(x = PHATE_1, y = PHATE_2)) +
  # 将点栅格化，防止进 Inkscape 卡顿 (600 DPI)
  ggrastr::rasterise(geom_point(aes(color = Cluster), size = 2, alpha = 0.7), dpi = 600) +
  # 叠加全矢量的黑色轨迹线与箭头
  geom_path(data = curve_data, aes(x = PHATE_1, y = PHATE_2), 
            color = "black", linewidth = 1.2, 
            arrow = arrow(type = "closed", length = unit(0.15, "inches"))) +
  # 使用 theme_void() 彻底扒光坐标系
  theme_void() +
  theme(
    text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
    legend.text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
    legend.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

cairo_pdf(paste0("Slingshot_Trajectory_", prefix, ".pdf"), width = 7, height = 6, family = "Arial")
print(p_traj)
dev.off()
cat("✅ 轨迹图生成完毕！\n")

# ==========================================
# 3. 提取 Modules 并拟合平滑曲线
# ==========================================
pt_data <- data.frame(Pseudotime = pt_values, obj_linear@meta.data[ordered_cells, module_names])

pt_seq <- seq(min(pt_data$Pseudotime), max(pt_data$Pseudotime), length.out = 100)
smooth_matrix <- matrix(NA, nrow = length(module_names), ncol = 100)
rownames(smooth_matrix) <- module_names

for (i in seq_along(module_names)) {
  m_name <- module_names[i]
  form <- as.formula(paste0(m_name, " ~ Pseudotime"))
  fit <- loess(form, data = pt_data, span = 0.75)
  smooth_matrix[i, ] <- predict(fit, newdata = data.frame(Pseudotime = pt_seq))
}

scaled_matrix <- t(apply(smooth_matrix, 1, scale))
colnames(scaled_matrix) <- paste0("Bin", 1:100)

# ==========================================
# ==========================================
# 4. 寻找峰值并实现对角线排序
# ==========================================
peak_idx <- apply(scaled_matrix, 1, which.max)

proCSC_mods <- paste0(prefix, "_Module", c(11, 10, 15, 4, 18, 22, 6, 8, 24, 3, 23), "_UCell")
CSC_mods    <- paste0(prefix, "_Module", c(14, 2, 16), "_UCell")
revCSC_mods <- paste0(prefix, "_Module", c(9, 1, 12, 7, 13, 17, 20, 21, 25, 5, 19), "_UCell")

proCSC_mods_sorted <- proCSC_mods[order(peak_idx[proCSC_mods])]
CSC_mods_sorted    <- CSC_mods[order(peak_idx[CSC_mods])]
revCSC_mods_sorted <- revCSC_mods[order(peak_idx[revCSC_mods])]

final_ordered_mods <- c(proCSC_mods_sorted, CSC_mods_sorted, revCSC_mods_sorted)
scaled_matrix <- scaled_matrix[final_ordered_mods, ]

# 将行名洗回 M1, M2 的格式
rownames(scaled_matrix) <- gsub(".*_Module([0-9]+)_UCell", "M\\1", rownames(scaled_matrix))

split_vector <- c(rep("proCSC Modules", length(proCSC_mods_sorted)),
                  rep("CSC Modules", length(CSC_mods_sorted)),
                  rep("revCSC Modules", length(revCSC_mods_sorted)))
split_vector <- factor(split_vector, levels = c("proCSC Modules", "CSC Modules", "revCSC Modules"))

# ==========================================
# 5. 画热图并输出 PDF (已修复白缝问题)
# ==========================================
col_fun <- colorRamp2(c(-2, 0, 2), c("#009ad6", "white", "#d71345"))
font_style <- gpar(fontsize = 14, col = "black", fontface = "plain", fontfamily = "Arial")

ht <- Heatmap(
  scaled_matrix,
  name = "Scaled\nExpression",
  col = col_fun,
  
  # --- 🌟 核心修复区：消除白色网格缝隙 ---
  use_raster = TRUE,         # 开启热图主体栅格化，颜色完美融合
  raster_quality = 4,        # 提高栅格化分辨率保证高清
  rect_gp = gpar(col = NA),  # 强制去除任何潜在的单元格边框线
  # ---------------------------------------
  
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 1),
  row_gap = unit(2, "mm"),
  
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_title = "Cells ordered along proCSC -> CSC -> revCSC trajectory",
  column_title_side = "bottom",
  
  cluster_rows = FALSE, 
  show_row_dend = FALSE,
  row_split = split_vector,
  row_title_rot = 0,
  row_title_side = "right",
  
  column_title_gp = font_style,
  row_names_gp = font_style,
  row_title_gp = font_style,
  
  heatmap_legend_param = list(
    title_gp = font_style,
    labels_gp = font_style,
    legend_height = unit(4, "cm"), 
    border = "black"
  )
)

output_file <- paste0("Heatmap_Pseudotime_", prefix, ".pdf")
cairo_pdf(output_file, width = 9, height = 8, family = "Arial")
draw(ht)
dev.off()
cat("✅ 热图生成完毕，烦人的白缝已经消失啦！\n")



#save(seurat_obj,file='figd_suerat.Rdata')
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrastr)

# 强制注册 Arial 字体
windowsFonts(Arial = windowsFont("Arial"))

# 1. Module annotations dictionary
module_anno <- c(
  "M1" = "Peptide hormone response",
  "M2" = "Wnt signaling",
  "M3" = "Ribosome biogenesis",
  "M4" = "Cytoplasmic translation",
  "M5" = "TGF-beta response",
  "M6" = "Chromosome segregation",
  "M7" = "Glutathione metabolism",
  "M8" = "Cell junction",
  "M9" = "Cell chemotaxis",
  "M10" = "DNA replication",
  "M11" = "Telomere maintenance",
  "M12" = "Cell-matrix adhesion",
  "M13" = "Multivesicular body assembly",
  "M14" = "Gastrointestinal epithelial maintenance",
  "M15" = "Smoothened signaling pathway",
  "M16" = "Sterol metabolism",
  "M17" = "Purine nucleotide biosynthesis",
  "M18" = "Developmental patterning regulation",
  "M19" = "Hypoxia response",
  "M20" = "Antiviral immunity",
  "M21" = "Lipid raft organization secretion",
  "M22" = "RNA splicing",
  "M23" = "Protein acetylation stemness control",
  "M24" = "Amino acid biosynthesis metabolism",
  "M25" = "Actin cytoskeleton remodeling migration"
)

# 🌟 核心修改：完全按照你截图里的生物学轨迹顺序进行硬编码恢复
final_ordered_mods <- paste0("M", c(11, 15, 4, 18, 22, 6, 8, 3, 23, 10, 24, 
                                    2, 14, 16, 
                                    1, 7, 17, 13, 20, 9, 12, 21, 25, 5, 19))

# 3. Process correlation matrix
cor_mat <- t(cor_sub)
cor_df <- as.data.frame(cor_mat)
cor_df$Module <- rownames(cor_df)

plot_data <- pivot_longer(cor_df, cols = -Module, names_to = "Subtype", values_to = "Correlation")

# 4. Filter thresholds and create variables for legend
plot_data <- plot_data %>%
  filter(abs(Correlation) >= 0.2) %>%
  mutate(
    Direction = ifelse(Correlation > 0, "Positive", "Negative"),
    Strength = ifelse(abs(Correlation) >= 0.4, 
                      "Strong (|r| >= 0.4)", 
                      "Moderate (0.2 <= |r| < 0.4)"),
    line_alpha = ifelse(abs(Correlation) >= 0.4, 0.8, 0.5)
  )

# 5. Build coordinate system
n_mods <- length(final_ordered_mods)

df_nodes_left <- data.frame(
  Module = final_ordered_mods,
  x = 0,
  y = rev(1:n_mods)
)
df_nodes_left$Label <- module_anno[as.character(df_nodes_left$Module)]

subtype_order <- c("iCMS2_MSS", "iCMS3_MSI-H", "iCMS3_MSS")

df_nodes_right <- data.frame(
  Subtype = subtype_order,
  x = 2.1,              # 保持向左平移，缩短连线距离
  y = c(21, 13, 5)
)

df_edges <- plot_data %>%
  inner_join(df_nodes_left, by = "Module") %>%
  dplyr::rename(x_left = x, y_left = y) %>%
  inner_join(df_nodes_right, by = "Subtype") %>%
  dplyr::rename(x_right = x, y_right = y)

# 6. Strict formatting theme (全局 7pt，纯黑，Arial，不加粗)
font_style <- element_text(size = 7, color = "black", face = "plain", family = "Arial")

my_theme_network <- theme_void() +
  theme(
    text = font_style,
    legend.position = "right",
    legend.title = font_style,
    legend.text = font_style,
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.1, "cm")
  )

# 7. Generate right network plot
p_right_network <- ggplot() +
  geom_curve(
    data = df_edges,
    aes(x = 1.1, y = y_left, xend = x_right - 0.15, yend = y_right,
        color = Direction, linetype = Strength, alpha = line_alpha),
    curvature = 0.15,
    linewidth = 0.5
  ) +
  scale_color_manual(
    name = "Direction", 
    values = c("Positive" = "#d71345", "Negative" = "#009ad6")
  ) +
  scale_linetype_manual(
    name = "Correlation",
    values = c("Strong (|r| >= 0.4)" = "solid", 
               "Moderate (0.2 <= |r| < 0.4)" = "dashed")
  ) +
  scale_alpha_identity() +
  geom_label(
    data = df_nodes_left,
    aes(x = x, y = y, label = Label),
    hjust = 0,
    size = 7 / .pt,
    color = "black",
    fill = "white",
    family = "Arial",
    fontface = "plain",
    label.padding = unit(0.15, "lines"),
    label.size = 0.2
  ) +
  geom_label(
    data = df_nodes_right,
    aes(x = x, y = y, label = Subtype),
    size = 7 / .pt,
    color = "black",
    fill = "white",
    family = "Arial",
    fontface = "plain",
    label.padding = unit(0.2, "lines"),
    label.size = 0.3
  ) +
  coord_cartesian(xlim = c(0, 2.7)) +
  my_theme_network +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1, alpha = 1)),
    linetype = guide_legend(override.aes = list(linewidth = 1, alpha = 1))
  )

# 🌟 按要求导出：高度 4.2，宽度 5
cairo_pdf("Right_Network_Modules_to_iCMS_Final_7pt.pdf", width = 5, height = 4.2, family = "Arial")
print(p_right_network)
dev.off()
cat("✅ 最终版网络图生成完毕！节点左移，连线缩短，画布比例优化 (5.4 x 4.5)！\n")

#fige
setwd('E:/deskup/毕业论文/结直肠癌/data/')
load('Human_icms.Rdata')

unique(scRNA@meta.data$iCMS_msi)


library(Seurat)
library(ggplot2)
library(patchwork)
library(ggrastr)
library(scales) # 必须加载 scales 以调用 squish 处理越界值

# 强制注册 Arial 字体
windowsFonts(Arial = windowsFont("Arial"))

# ==========================================
# 统一的主题 (根据你提供的代码设定为 14号, 黑色, Arial, 不加粗)
# ==========================================
font_style <- element_text(size = 14, color = "black", face = "plain", family = "Arial")

my_theme_tsne <- theme(
  text = font_style,
  plot.title = element_text(size = 14, color = "black", face = "plain", family = "Arial", hjust = 0.5),
  legend.text = font_style,
  legend.title = font_style
)

# ==========================================
# 1. t-SNE 降维图 (解决透明/颜色变淡问题)
# ==========================================
p_dim <- DimPlot(
  scRNA, 
  reduction = 'tsne', 
  group.by = "iCMS_msi",
  pt.size = 0.5,
  raster = FALSE # 🌟 核心：关闭 Seurat 自带的劣质栅格化
) + 
  scale_color_manual(values = c("iCMS2_MSS" = "#c88400", "iCMS3_MSI-H" = "#00a6ac", "iCMS3_MSS" = "#009ad6")) +
  NoAxes() + 
  my_theme_tsne

# 🌟 使用 ggrastr 进行专业栅格化，色彩零损失
p_dim <- rasterise(p_dim, layers = "Point", dpi = 600)

cairo_pdf("tsne_plot_no_axes_Fixed.pdf", width = 8, height = 6, family = "Arial")
print(p_dim)
dev.off()

# ==========================================
# 2. UCell 打分图 (解决大量空白消失问题)
# ==========================================
rm(seurat_obj)
gc()



library(Seurat)
library(ggplot2)
library(patchwork)

# 强制注册 Arial 字体
windowsFonts(Arial = windowsFont("Arial"))

# ==========================================
# 统一的主题 (严格 14号, 黑色, Arial, 不加粗)
# ==========================================
font_style <- element_text(size = 14, color = "black", face = "plain", family = "Arial")

my_theme_tsne <- theme(
  text = font_style,
  plot.title = element_text(size = 14, color = "black", face = "plain", family = "Arial", hjust = 0.5),
  legend.text = font_style,
  legend.title = font_style
)

# ==========================================
# 1. t-SNE 降维图 (解决透明发虚问题)
# ==========================================
# 调大 pt.size 至 1.2，抵消栅格化带来的透明/褪色感
p_dim <- DimPlot(
  scRNA,
  reduction = 'tsne',
  group.by = "iCMS_msi",
  pt.size = 1.2,
  raster = TRUE,            # 开启栅格化
  raster.dpi = c(600, 600)
) +
  scale_color_manual(
    values = c("iCMS2_MSS" = "#c88400", 
               "iCMS3_MSI-H" = "#00a6ac", 
               "iCMS3_MSS" = "#009ad6")
  ) +
  NoAxes() +
  my_theme_tsne

# 使用 cairo_pdf 确保 Arial 字体完美写入
cairo_pdf("tsne_plot_no_axes_Fixed.pdf", width = 8, height = 6, family = "Arial")
print(p_dim)
dev.off()


scRNA <- UCell::AddModuleScore_UCell(
  obj = scRNA,
  features = human_modules,
  ncores = 24,
  name = "_UCell" # UCell 会自动在模块名后缀加上 "_UCell" 以防与原基因名冲突
)
# ==========================================
# 2. UCell 打分图 (解决过红问题)
# ==========================================
DefaultAssay(scRNA) <- "UCell_M"

max_score_M14 <- max(GetAssayData(scRNA, assay = "UCell_M")["M14", ])
max_score_M23 <- max(GetAssayData(scRNA, assay = "UCell_M")["M23", ])

global_max <- max(max_score_M14, max_score_M23)
threshold <- 0.05
my_breaks <- c(0.1, 0.2, 0.3, 0.4, 0.5)

p_M14 <- FeaturePlot(
  scRNA,
  reduction = "tsne",
  features = "M14",
  pt.size = 0.5,
  # 彻底去掉 order = TRUE，让红灰自然叠加，反映真实分布密度
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  scale_color_gradientn(
    colors = c("#d3d7d4", "#ed1941"),
    limits = c(threshold, global_max),
    breaks = my_breaks,
    na.value = "#d3d7d4",
    name = "UCell Score"
  ) +
  NoAxes() +
  my_theme_tsne +
  ggtitle("Module 14 (tSNE)")

p_M23 <- FeaturePlot(
  scRNA,
  reduction = "tsne",
  features = "M23",
  pt.size = 0.5,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  scale_color_gradientn(
    colors = c("#d3d7d4", "#ed1941"),
    limits = c(threshold, global_max),
    breaks = my_breaks,
    na.value = "#d3d7d4",
    name = "UCell Score"
  ) +
  NoAxes() +
  my_theme_tsne +
  ggtitle("Module 23 (tSNE)")

cairo_pdf("tsne_scRNA_M14_Filtered_Fixed.pdf", width = 8, height = 6, family = "Arial")
print(p_M14)
dev.off()

cairo_pdf("tsne_scRNA_M23_Filtered_Fixed.pdf", width = 8, height = 6, family = "Arial")
print(p_M23)
dev.off()

cat("✅ 所有 t-SNE 图件生成完毕！\n")




library(Seurat)
library(ggplot2)
library(homologene)
library(dplyr)

windowsFonts(Arial = windowsFont("Arial"))

# ============================================================
# 1. Human Module 2
# ============================================================
M2_human <- human_modules$Module_2

# 正文 5 个核心基因
main_genes_human <- c(
  "AXIN2",
  "EPHB2",
  "MYC",
  "LGR5",
  "RNF43"
)

# 剩余基因
remaining_genes_human <- setdiff(
  M2_human,
  main_genes_human
)


# ============================================================
# 2. Human -> Mouse homolog mapping
# ============================================================
get_human_mouse_map <- function(human_genes, seurat_mouse) {
  
  trans_df <- homologene(
    human_genes,
    inTax = 9606,
    outTax = 10090
  )
  
  trans_df <- trans_df[
    !is.na(trans_df$`9606`) &
      !is.na(trans_df$`10090`),
  ]
  
  trans_df <- trans_df[
    trans_df$`10090` %in% rownames(seurat_mouse[["RNA"]]),
  ]
  
  trans_df <- trans_df[
    !duplicated(trans_df$`9606`),
  ]
  
  trans_df
}


# ============================================================
# 3. 固定细胞类型顺序
# ============================================================
seurat_sub$curatedCLUST <- factor(
  seurat_sub$curatedCLUST,
  levels = c(
    "proCSC",
    "CSC",
    "revCSC"
  )
)

scRNA$iCMS_msi <- factor(
  scRNA$iCMS_msi,
  levels = c(
    "iCMS2_MSS",
    "iCMS3_MSI-H",
    "iCMS3_MSS"
  )
)


# ============================================================
# 4. 通用函数（支持自定义字号与旋转角度）
# ============================================================
make_combined_dotplot <- function(
    human_genes,
    output_file,
    pdf_width = 8,
    pdf_height = 6,
    font_size = 14,
    x_angle = 0,
    x_hjust = 0.5
) {
  
  # ----------------------------------------------------------
  # A. Human -> Mouse
  # ----------------------------------------------------------
  gene_map <- get_human_mouse_map(
    human_genes,
    seurat_sub
  )
  
  mouse_genes <- gene_map$`10090`
  human_labels_mouse <- gene_map$`9606`
  
  # ----------------------------------------------------------
  # B. Human genes 在 scRNA 中存在的部分
  # ----------------------------------------------------------
  human_genes_use <- human_genes[
    human_genes %in% rownames(scRNA[["RNA"]])
  ]
  
  # ----------------------------------------------------------
  # C. Mouse DotPlot data
  # ----------------------------------------------------------
  dp_mouse <- DotPlot(
    seurat_sub,
    features = mouse_genes,
    group.by = "curatedCLUST",
    assay = "RNA"
  )
  
  mouse_df <- dp_mouse$data
  
  mouse_to_human <- setNames(
    human_labels_mouse,
    mouse_genes
  )
  
  mouse_df$Gene <- mouse_to_human[
    as.character(mouse_df$features.plot)
  ]
  
  mouse_df$Group <- as.character(mouse_df$id)
  
  # ----------------------------------------------------------
  # D. Human DotPlot data
  # ----------------------------------------------------------
  dp_human <- DotPlot(
    scRNA,
    features = human_genes_use,
    group.by = "iCMS_msi",
    assay = "RNA"
  )
  
  human_df <- dp_human$data
  
  human_df$Gene <- as.character(
    human_df$features.plot
  )
  
  human_df$Group <- as.character(
    human_df$id
  )
  
  # ----------------------------------------------------------
  # E. 合并
  # ----------------------------------------------------------
  plot_df <- bind_rows(
    mouse_df[, c("Gene", "Group", "pct.exp", "avg.exp.scaled")],
    human_df[, c("Gene", "Group", "pct.exp", "avg.exp.scaled")]
  )
  
  # ----------------------------------------------------------
  # F. 固定 Gene 顺序
  # ----------------------------------------------------------
  gene_order <- human_genes[
    human_genes %in% unique(plot_df$Gene)
  ]
  
  plot_df$Gene <- factor(
    plot_df$Gene,
    levels = gene_order
  )
  
  # ----------------------------------------------------------
  # G. 固定纵轴顺序
  # ----------------------------------------------------------
  group_order <- c(
    "proCSC",
    "CSC",
    "revCSC",
    "iCMS2_MSS",
    "iCMS3_MSI-H",
    "iCMS3_MSS"
  )
  
  plot_df$Group <- factor(
    plot_df$Group,
    levels = rev(group_order)
  )
  
  # ----------------------------------------------------------
  # H. 画图
  # ----------------------------------------------------------
  p <- ggplot(
    plot_df,
    aes(
      x = Gene,
      y = Group
    )
  ) +
    geom_point(
      aes(
        size = pct.exp,
        color = avg.exp.scaled
      )
    ) +
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
    geom_hline(
      yintercept = 3.5,
      linewidth = 0.5,
      color = "black"
    ) +
    theme_classic() +
    theme(
      text = element_text(
        size = font_size,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      axis.text.x = element_text(
        size = font_size,
        color = "black",
        face = "plain",
        family = "Arial",
        angle = x_angle,
        hjust = x_hjust
      ),
      axis.text.y = element_text(
        size = font_size,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      axis.title = element_blank(),
      legend.text = element_text(
        size = font_size,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      legend.title = element_text(
        size = font_size,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      legend.position = "right"
    )
  
  # ----------------------------------------------------------
  # I. 保存
  # ----------------------------------------------------------
  cairo_pdf(
    output_file,
    width = pdf_width,
    height = pdf_height,
    family = "Arial"
  )
  
  print(p)
  dev.off()
  
  cat("✔ 已保存:", output_file, "\n")
}


# ============================================================
# 5. PDF 1: 5 个核心 WNT genes (保持 14pt，平铺)
# ============================================================
make_combined_dotplot(
  human_genes = main_genes_human,
  output_file = "DotPlot_M2_5_WNT_Genes.pdf",
  pdf_width = 8,
  pdf_height = 6,
  font_size = 14,
  x_angle = 0,
  x_hjust = 0.5
)


# ============================================================
# 6. PDF 2: 剩余 Module 2 genes (7pt，宽度拉长至 25，基因名旋转 45 度)
# ============================================================
make_combined_dotplot(
  human_genes = remaining_genes_human,
  output_file = "DotPlot_M2_Remaining_Genes.pdf",
  pdf_width = 12,
  pdf_height = 6,
  font_size = 7,
  x_angle = 45,
  x_hjust = 1
)

cat("✅ 两个合并 DotPlot 全部生成完成！\n")