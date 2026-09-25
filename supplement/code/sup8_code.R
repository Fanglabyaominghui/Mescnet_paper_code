library(CellChat)
library(ggplot2)
library(dplyr)
library(patchwork)


#figb
# 读取已保存的 CellChat 对象
cellchat <- readRDS("CellChat_native_spatial_Ycropped_complete.rds")

# 计算网络中心性
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# 指定通路名称
pathway_show <- "COLLAGEN"

centr_data <- cellchat@netP$centr[[pathway_show]]

# 1. Sender 数据处理：升序排列
# 使得画图时最大值在最后面（也就是 Y 轴的最上面），实现从上到下“多到少”
sender_score <- data.frame(
  Cell = names(centr_data$outdeg),
  Score = centr_data$outdeg
) %>% filter(Score > 0) %>% arrange(Score)

sender_score$Cell <- factor(sender_score$Cell, levels = sender_score$Cell)

# 2. Receiver 数据处理：降序排列 (修改了这里)
# 使得画图时最小值在最后面（也就是 Y 轴的最上面），实现从上到下“少到多”
receiver_score <- data.frame(
  Cell = names(centr_data$indeg),
  Score = centr_data$indeg
) %>% filter(Score > 0) %>% arrange(desc(Score))

receiver_score$Cell <- factor(receiver_score$Cell, levels = receiver_score$Cell)

# 绘图主题设置 (严格14号、黑色、不加粗)
my_theme <- theme_classic() +
  theme(
    text = element_text(size = 14, color = "black", face = "plain", family = "sans"),
    axis.text = element_text(size = 14, color = "black", face = "plain"),
    axis.title = element_text(size = 14, color = "black", face = "plain"),
    plot.title = element_text(size = 14, color = "black", face = "plain", hjust = 0.5),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
    axis.line = element_blank(),
    axis.ticks.y = element_blank()
  )

# 绘制 Sender 图
p_sender <- ggplot(sender_score, aes(x = Score, y = Cell)) +
  geom_segment(aes(x = 0, xend = Score, y = Cell, yend = Cell), color = "#2b8cbe", linewidth = 1) +
  geom_point(color = "#2b8cbe", size = 6) +
  labs(title = paste(pathway_show, "Signaling Sender"), x = "Sender Weighted Score", y = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  my_theme

# 绘制 Receiver 图 
p_receiver <- ggplot(receiver_score, aes(x = Score, y = Cell)) +
  geom_segment(aes(x = 0, xend = Score, y = Cell, yend = Cell), color = "#e6550d", linewidth = 1) +
  geom_point(color = "#e6550d", size = 6) +
  labs(title = paste(pathway_show, "Signaling Receiver"), x = "Receiver Weighted Score", y = NULL) +
  scale_x_reverse(expand = expansion(mult = c(0.1, 0))) +
  scale_y_discrete(position = "right") +
  my_theme

# 拼图并保存
combined_plot <- p_sender + p_receiver + plot_layout(widths = c(1, 1))

pdf("FigE_Signaling_Lollipop_Reversed.pdf", width = 8, height = 5)
print(combined_plot)
dev.off()



#figc
library(CellChat)
library(ggplot2)
library(dplyr)

# 1. 读取 CellChat 对象
cellchat <- readRDS("CellChat_native_spatial_Ycropped_complete.rds")

# 2. 提取 COLLAGEN 通路中 Fibroblast 为发送方的信号
df_col <- subsetCommunication(cellchat, signaling = "COLLAGEN")
df_fibro <- df_col[df_col$source == "Fibroblast", ]

# 3. 按配体-受体对聚合通讯强度
pair_strength <- aggregate(prob ~ interaction_name_2, data = df_fibro, sum)

# 4. 降序排列并计算百分比
pair_strength <- pair_strength %>%
  arrange(desc(prob)) %>%
  mutate(
    Percentage = prob / sum(prob) * 100,
    Label = sprintf("%s\n%.1f%%", interaction_name_2, Percentage)
  )

# 5. 截取前 5 个
top_n <- 5
if (nrow(pair_strength) > top_n) {
  pair_strength_top <- pair_strength[1:top_n, ]
} else {
  pair_strength_top <- pair_strength
}

pair_strength_top$interaction_name_2 <- factor(
  pair_strength_top$interaction_name_2,
  levels = pair_strength_top$interaction_name_2
)

# 6. 绘制环形图
p_donut <- ggplot(pair_strength_top, aes(x = 2, y = Percentage, fill = interaction_name_2)) +
  geom_bar(stat = "identity", color = "white", linewidth = 1) +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    size = 5,
    color = "black",
    family = "Arial",
    fontface = "plain"
  ) +
  theme_void() +
  theme(
    text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
    legend.text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
    legend.title = element_blank(),
    legend.position = "right"
  )

# 7. 导出 PDF
cairo_pdf("FigF_Fibroblast_COLLAGEN_Donut_Arial.pdf", width = 7, height = 5, family = "Arial")
print(p_donut)
dev.off()


























############################################################
# TCGA LUAD Feature Selection: Univariate Cox Only
# Formatting: Arial, 14pt, Black, Plain
############################################################

rm(list = ls())
gc()

# 加载必要包
library(survival)
library(dplyr)
library(ggplot2)

# ==========================================================
# 全局绘图主题 (强制 Arial, 14号, 黑色, 不加粗)
# ==========================================================
my_theme <- theme_classic() +
  theme(
    text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
    axis.text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
    axis.title = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
    plot.title = element_text(size = 14, color = "black", face = "plain", hjust = 0.5, family = "Arial"),
    legend.text = element_text(size = 14, color = "black", face = "plain", family = "Arial"),
    legend.title = element_text(size = 14, color = "black", face = "plain", family = "Arial")
  )

# ==========================================================
# 1. 加载数据
# ==========================================================
setwd('E:/deskup/毕业论文/肺癌/tcga/')
load("TCGA_LUAD_M1_M3_Cox_analysis_workspace.Rdata")

# 修复：截取 TCGA Barcode 第 14-15 位，"01" 代表肿瘤样本
is_tumor <- substr(colnames(expr_log), 14, 15) == "01"
tumor_cols <- colnames(expr_log)[is_tumor]

if(length(tumor_cols) > 0) {
  expr_log <- expr_log[, tumor_cols]
  cat(sprintf("Successfully extracted %d tumor samples.\n", length(tumor_cols)))
} else {
  stop("Error: No tumor samples found. Please check column names.")
}

# ==========================================================
# 2. 提取 M1/M3 基因并合并临床信息
# ==========================================================
M1_genes <- c(
  "CLDN1","KRT15","KRT17","TP63","NTS","LMO4","SFN","TNFSF10","FAM43A",
  "CD9","GABRE","FGFR2","BCL6","RGMA","VTCN1","PERP","NFE2L2","LY6E",
  "TMEM44","DSC3","ANKRD62","GSTM3","ADH7","ABCC5","IGF2BP2","EPCAM",
  "KRT19","NTRK2","CES1","TNC","KREMEN1","PTPRZ1","GCLC","ACAP2",
  "DVL3","TFRC","KRT5","RNF7","ATP13A3","PCYT1A","CSTA","EPHB3",
  "NCBP2","FXYD3","PTPRF","MELTF","WNK2","TBL1XR1","GPX2","HSPB1",
  "NDUFB5","EIF4A2","EHF","BDH1","MCCC1","DDR1","DLG1","CASK",
  "TMPRSS4","EHBP1","SOX2","LSG1","CALML3","PITX1","DYNLT2B","BEX3",
  "ZNF639","IFI27","TACSTD2","C12orf75","ABCF3","TNS4","GPC3",
  "NDUFB9","CDK4","EIF2B5","HSP90AB1"
)

M3_genes <- c(
  "PIGR","SERPINA1","MUC5B","DMBT1","OLFM4","TFF2","IGFBP2","LYZ",
  "FCGBP","AGR2","ATP1B1","TFF3","CA9","SERPINA3","GOLM1","PROM1",
  "CRISP3","MUC6","S100A10","PGC","LCN2","GCNT3","WFDC2","QSOX1",
  "SPINK4","MUC1","CREB3L1","MLPH","CLDN2","SLC44A4","CP","CD24",
  "CEACAM5","CTSE","GP2","LGALS4","SPINK1","MIA","GFUS","TFF1",
  "RASSF7","KRT8","S100A6","MUC5AC","DUOX2","FAM3D","FAM3B","TMC5",
  "BACE2","SLC4A4","ANXA10","ELF3","TMPRSS3","TSPAN8","AQP5","SOX9",
  "MMP1","ERN2","CYSTM1","CEACAM6","NRG1","GMDS","ANXA4"
)

candidate_genes <- intersect(unique(c(M1_genes, M3_genes)), rownames(expr_log))

# 转置并附上病人 ID (TCGA Barcode 前 12 位)
expr_tumor <- as.data.frame(t(expr_log[candidate_genes, ]))
expr_tumor$patient <- substr(rownames(expr_tumor), 1, 12)

# 处理同一患者有多个测序样本的情况 (取均值)
expr_patient <- expr_tumor %>%
  group_by(patient) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

# 合并临床信息
model_data <- inner_join(clinical_patient, expr_patient, by = "patient") %>%
  filter(is.finite(OS), OS > 30, vital %in% c(0, 1))

cat(sprintf("Final merged dataset contains %d patients.\n", nrow(model_data)))

# ==========================================================
# 3. 单因素 Cox 分析
# ==========================================================
cat("\nRunning Univariate Cox Analysis...\n")

cox_result <- data.frame()
for(g in candidate_genes){
  fit <- coxph(as.formula(paste0("Surv(OS, vital) ~ `", g, "`")), data = model_data)
  tmp <- summary(fit)
  cox_result <- rbind(cox_result, data.frame(
    gene = g,
    HR = tmp$coefficients[1,2],
    lower = tmp$conf.int[1,3],
    upper = tmp$conf.int[1,4],
    p = tmp$coefficients[1,5]
  ))
}

# 标记基因来源模块
cox_result$module <- ifelse(cox_result$gene %in% M1_genes, "M1 LUSC", "M3 LUAD")

# 筛选显著基因
cox_plot_data <- cox_result %>%
  filter(p < 0.05) %>%
  arrange(HR)

# 锁定因子顺序，保证从上到下按 HR 排序
cox_plot_data$gene <- factor(cox_plot_data$gene, levels = cox_plot_data$gene)

cat(sprintf("Found %d significant genes (p < 0.05).\n", nrow(cox_plot_data)))

# ==========================================================
# 4. 绘制单因素 Cox 森林图
# ==========================================================
p1_forest <- ggplot(cox_plot_data, aes(x = HR, y = gene)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = HR > 1), height = 0.2, linewidth = 1) +
  geom_point(aes(shape = module, color = HR > 1), size = 4) +
  scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"),
                     labels = c("HR < 1", "HR > 1")) +
  geom_vline(xintercept = 1, linetype = 2, color = "black") +
  labs(title = "Univariate Cox (p < 0.05)",
       x = "Hazard Ratio", y = NULL,
       color = "Hazard Ratio", shape = "Module") +
  my_theme



cairo_pdf("Plot1_Forest_Only.pdf", width = 8, height = 8, family = "Arial")
print(p1_forest)
dev.off()

cat("✅ 森林图已生成完毕，保存在 Plot1_Forest_Only.pdf\n")





# ==========================================================
# 5. 拆分训练集与验证集，准备 LASSO
# ==========================================================
library(glmnet)
library(caret)
library(tidyr)

cat("\nSplitting data into Train and Test sets (2:1)...\n")
set.seed(123)

# 提取 Cox 回归中显著的基因
cox_sig_genes <- as.character(cox_plot_data$gene)

train_index <- createDataPartition(model_data$vital, p = 2/3, list = FALSE)
train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

x_train <- as.matrix(train_data[, cox_sig_genes])
y_train <- Surv(train_data$OS, train_data$vital)

# ==========================================================
# 6. 执行 LASSO 降维与特征提取
# ==========================================================
cat("Running LASSO Regression...\n")

# 拟合 LASSO 模型与交叉验证
fit_lasso <- glmnet(x_train, y_train, family = "cox", alpha = 1)
cv_fit <- cv.glmnet(x_train, y_train, family = "cox", alpha = 1, nfolds = 10)

# 获取最佳 Lambda 值 (误差最小)
best_lambda <- cv_fit$lambda.min

# 提取核心基因及其对应系数
coef_lasso <- coef(cv_fit, s = best_lambda)
lasso_genes <- rownames(coef_lasso)[coef_lasso[, 1] != 0]
weights <- as.numeric(coef_lasso[coef_lasso[, 1] != 0])

cat(sprintf("LASSO selected %d core genes.\n", length(lasso_genes)))

# ==========================================================
# 7. 绘制 LASSO 轨迹图 (Plot 2)
# ==========================================================
beta_data <- as.data.frame(as.matrix(t(fit_lasso$beta)))
beta_data$LogLambda <- log(fit_lasso$lambda)
beta_long <- pivot_longer(beta_data, cols = -LogLambda, names_to = "Gene", values_to = "Coefficient")

p2_lasso_path <- ggplot(beta_long, aes(x = LogLambda, y = Coefficient, color = Gene)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = log(best_lambda), linetype = "dashed", color = "black", linewidth = 1) +
  labs(title = "LASSO Coefficient Path", x = "Log(Lambda)", y = "Coefficients") +
  my_theme + theme(legend.position = "none")

cairo_pdf("Plot2_LASSO_Path.pdf", width = 6, height = 5, family = "Arial")
print(p2_lasso_path)
dev.off()

# ==========================================================
# 8. 绘制 LASSO 交叉验证图 (Plot 3)
# ==========================================================
cv_data <- data.frame(
  LogLambda = log(cv_fit$lambda),
  CVM = cv_fit$cvm, CVUP = cv_fit$cvup, CVLO = cv_fit$cvlo
)

p3_lasso_cv <- ggplot(cv_data, aes(x = LogLambda, y = CVM)) +
  geom_errorbar(aes(ymin = CVLO, ymax = CVUP), width = 0.05, color = "grey50") +
  geom_point(color = "#E64B35", size = 2) +
  geom_vline(xintercept = log(best_lambda), linetype = "dashed", color = "black", linewidth = 1) +
  labs(title = "LASSO Cross-Validation", x = "Log(Lambda)", y = "Partial Likelihood Deviance") +
  my_theme

cairo_pdf("Plot3_LASSO_CV.pdf", width = 6, height = 5, family = "Arial")
print(p3_lasso_cv)
dev.off()

cat("✅ LASSO 轨迹图与交叉验证图已生成完毕。\n")












# ==========================================================
# 8. 绘制训练集和测试集 KM 生存曲线
# 年为单位 + 置信区间 + Number at risk + 0-10年刻度 + 固定阈值防泄露
# ==========================================================

library(survival)
library(survminer)
library(ggplot2)
library(patchwork)

# ==========================================================
# 8.5 [核心修复] 在训练集上计算全局固定的最优截断值
# ==========================================================
train_tmp <- train_data
train_tmp$OS_year <- train_tmp$OS / 365.25

cut_obj <- surv_cutpoint(
  train_tmp,
  time = "OS_year",
  event = "vital",
  variables = "RiskScore"
)

# 提取固定的阈值数值
fixed_cutoff <- cut_obj$cutpoint$cutpoint[1]
cat(sprintf("全局固定最佳 RiskScore 截断值为: %.4f\n", fixed_cutoff))

# ==========================================================
# KM 绘图函数 (引入 cutoff_value 参数)
# ==========================================================
get_km_plot <- function(data, title_name, cutoff_value) {
  
  # ----------------------------------------------------------
  # 1. OS：天 -> 年
  # ----------------------------------------------------------
  d_tmp <- data
  d_tmp$OS_year <- d_tmp$OS / 365.25
  
  # ----------------------------------------------------------
  # 2 & 3. 统一使用固定阈值划分 High / Low Risk
  # ----------------------------------------------------------
  d_tmp$RiskGroup <- ifelse(d_tmp$RiskScore > cutoff_value, "High Risk", "Low Risk")
  d_tmp$RiskGroup <- factor(d_tmp$RiskGroup, levels = c("High Risk", "Low Risk"))
  
  # ----------------------------------------------------------
  # 4. KM 拟合
  # ----------------------------------------------------------
  fit <- survfit(
    Surv(OS_year, vital) ~ RiskGroup,
    data = d_tmp
  )
  
  # ----------------------------------------------------------
  # 5. Log-rank P value
  # ----------------------------------------------------------
  diff <- survdiff(
    Surv(OS_year, vital) ~ RiskGroup,
    data = d_tmp
  )
  
  pValue <- 1 - pchisq(
    diff$chisq,
    df = 1
  )
  
  if (pValue < 0.001) {
    pValue_str <- "p < 0.001"
  } else {
    pValue_str <- paste0(
      "p = ",
      sprintf("%.3f", pValue)
    )
  }
  
  # ----------------------------------------------------------
  # 6. KM 主图 + Risk table
  # ----------------------------------------------------------
  sur <- ggsurvplot(
    
    fit,
    data = d_tmp,
    
    # ==========================
    # 生存曲线
    # ==========================
    conf.int = TRUE,
    conf.int.alpha = 0.20,
    
    pval = pValue_str,
    pval.size = 5,
    
    censor = TRUE,
    censor.size = 2,
    
    size = 1.2,
    
    # ==========================
    # Legend
    # ==========================
    legend.labs = c(
      "High Risk",
      "Low Risk"
    ),
    
    legend.title = title_name,
    
    # ==========================
    # 坐标轴
    # ==========================
    xlab = "Time (years)",
    ylab = "Survival probability",
    
    break.time.by = 1,
    
    # ==========================
    # Risk table
    # ==========================
    risk.table = "absolute",
    
    risk.table.title = "Number at risk",
    
    risk.table.height = 0.28,
    
    risk.table.y.text = FALSE,
    
    risk.table.col = "strata",
    
    # ==========================
    # 颜色
    # ==========================
    palette = c(
      "#E64B35",
      "#4DBBD5"
    ),
    
    # ==========================
    # 主图主题
    # 注意：这里不再写 family="Arial"
    # 避免 grid 把 Arial 错当成颜色
    # ==========================
    ggtheme = theme_classic(
      base_size = 14
    ) +
      theme(
        
        text = element_text(
          size = 14,
          color = "black",
          face = "plain"
        ),
        
        panel.grid.major = element_line(
          color = "gray90",
          linewidth = 0.4
        ),
        
        panel.grid.minor = element_blank(),
        
        axis.title = element_text(
          face = "plain",
          size = 13,
          color = "black"
        ),
        
        axis.text = element_text(
          size = 11,
          color = "black"
        ),
        
        legend.title = element_text(
          face = "plain",
          size = 12,
          color = "black"
        ),
        
        legend.text = element_text(
          size = 11,
          color = "black"
        ),
        
        plot.margin = margin(
          10,
          10,
          10,
          10
        )
      ),
    
    # ==========================
    # Risk table 主题
    # 同样删除 family="Arial"
    # ==========================
    tables.theme = theme_classic(
      base_size = 14
    ) +
      theme(
        
        text = element_text(
          size = 14,
          color = "black",
          face = "plain"
        ),
        
        axis.text.x = element_text(
          size = 11,
          color = "black"
        ),
        
        axis.title.x = element_text(
          size = 13,
          color = "black",
          face = "plain"
        ),
        
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        
        plot.title = element_text(
          size = 14,
          color = "black",
          face = "plain",
          hjust = 0
        ),
        
        plot.margin = margin(
          5,
          10,
          5,
          10
        )
      )
  )
  
  # ==========================================================
  # 7. 主图横坐标：显示 0-10 年
  # ==========================================================
  
  sur$plot <- sur$plot +
    scale_x_continuous(
      breaks = seq(0, 10, by = 1)
    ) +
    coord_cartesian(
      xlim = c(0, 10)
    )
  
  # ==========================================================
  # 8. Risk table 横坐标：0-10 年
  # ==========================================================
  
  sur$table <- sur$table +
    scale_x_continuous(
      breaks = seq(0, 10, by = 1)
    ) +
    coord_cartesian(
      xlim = c(0, 10)
    )
  
  return(sur)
}


# ==========================================================
# 9. Training cohort (传入固定阈值 fixed_cutoff)
# ==========================================================
p_km_train_obj <- get_km_plot(
  train_data,
  "Training Cohort",
  fixed_cutoff
)


# ==========================================================
# 10. Validation cohort (传入固定阈值 fixed_cutoff)
# ==========================================================
p_km_test_obj <- get_km_plot(
  test_data,
  "Validation Cohort",
  fixed_cutoff
)


# ==========================================================
# 11. 左右拼图
# ==========================================================
p_km_combined <- arrange_ggsurvplots(
  
  list(
    p_km_train_obj,
    p_km_test_obj
  ),
  
  print = FALSE,
  ncol = 2,
  nrow = 1
)


# ==========================================================
# 12. 保存 PDF
# ==========================================================
grDevices::cairo_pdf(
  filename = "Plot_KM_Train_Test_ByYear_Combined_Fixed.pdf",
  width = 9,
  height = 6.2,
  family = "Arial"
)

print(p_km_combined)

dev.off()

cat("✅ 修正完毕：两组已强制采用相同的最优截断值，避免了过拟合和人数比例失衡。\n")