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