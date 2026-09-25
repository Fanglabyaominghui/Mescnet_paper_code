#sup
# ============================================================
# Fixed Module Color Library
# ============================================================

module_color_library <- c(
  
  # 1-10
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  
  # 11-20
  "#C77CFF", "#00BFC4", "#F8766D", "#7CAE00", "#619CFF",
  "#C49A00", "#00BA38", "#F564E3", "#00B0F6", "#A58AFF",
  
  # 21-30
  "#E76F51", "#2A9D8F", "#E9C46A", "#457B9D", "#A8DADC",
  "#F4A261", "#6D597A", "#B56576", "#355070", "#81B29A",
  
  # 31-40
  "#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#666666", "#8DD3C7", "#BEBADA",
  
  # 41-50
  "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
  "#BC80BD", "#CCEBC5", "#FFED6F", "#A6CEE3", "#1F78B4",
  
  # 51-60
  "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
  "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
  
  # 61-70
  "#8DA0CB", "#FC8D62", "#66C2A5", "#E78AC3", "#A6D854",
  "#FFD92F", "#E5C494", "#B3B3B3", "#1F968B", "#440154",
  
  # 71-80
  "#31688E", "#35B779", "#FDE725", "#D53E4F", "#3288BD",
  "#ABDDA4", "#FEE08B", "#FDAE61", "#5E4FA2", "#9E0142",
  
  # 81-90
  "#542788", "#998EC3", "#D8DAEB", "#F1A340", "#B35806",
  "#01665E", "#5AB4AC", "#C7EAE5", "#D8B365", "#8C510A",
  
  # 91-100
  "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8", "#D9F0D3",
  "#A6DBA0", "#5AAE61", "#1B7837", "#E08214", "#8073AC"
)

stopifnot(length(module_color_library) >= 100)
#fig1
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
#fig1
seurat_obj<-pick_power_scale_free_power(seurat_obj)







#fig2
library(uwot)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggrepel)

# 强制注册 Arial 字体
windowsFonts(Arial = windowsFont("Arial"))

library(uwot)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(grid)

windowsFonts(
  Arial = windowsFont("Arial")
)


# ============================================================
# 工具函数：从各种字符串中提取原始 module 编号
#
# M78          -> M78
# ME78         -> M78
# kME_M78      -> M78
# kME_MM.M78   -> M78
# 无数字       -> NA
# ============================================================

normalize_module_id <- function(x) {
  
  x <- as.character(x)
  
  has_number <- grepl(
    "[0-9]+",
    x
  )
  
  out <- rep(
    NA_character_,
    length(x)
  )
  
  out[has_number] <- paste0(
    "M",
    sub(
      ".*?([0-9]+)[^0-9]*$",
      "\\1",
      x[has_number]
    )
  )
  
  out
}


# ============================================================
# 主函数
#
# 关键新增参数：
# module_order_raw
#
# 它必须是“原始 module 的正确顺序”
# 例如：
# M78, M5, M20, M3, ...
#
# 函数内部自动转换：
# M78 -> M1
# M5  -> M2
# M20 -> M3
# ...
# ============================================================

plot_supervised_module_network <- function(
    tom,
    datKME,
    moduleColors,
    module_assignment_df,
    module_order_raw,
    n_hubs = 10,
    n_genes_per_mod = 50,
    edge_threshold = 0.25,
    edge_sample = 0.01,
    output_name = "Supervised_Network"
) {
  
  
  # ==========================================================
  # 1. 原始 module -> 连续 module
  # ==========================================================
  
  cat(
    "\n========================================\n"
  )
  
  cat(
    "1. 建立原始 module -> 连续 module 映射\n"
  )
  
  
  module_order_raw <- normalize_module_id(
    module_order_raw
  )
  
  module_order_raw <- module_order_raw[
    !is.na(module_order_raw)
  ]
  
  module_order_raw <- unique(
    module_order_raw
  )
  
  
  n_modules <- length(
    module_order_raw
  )
  
  
  new_module_names <- paste0(
    "M",
    seq_len(n_modules)
  )
  
  
  module_map <- setNames(
    new_module_names,
    module_order_raw
  )
  
  
  mapping_df <- data.frame(
    Original_Module = module_order_raw,
    New_Module = new_module_names,
    stringsAsFactors = FALSE
  )
  
  
  cat(
    "\nModule mapping:\n"
  )
  
  print(
    mapping_df,
    row.names = FALSE
  )
  
  
  cat(
    "\nTotal modules:",
    n_modules,
    "\n"
  )
  
  
  # ==========================================================
  # 2. 检查监督 assignment
  # ==========================================================
  
  module_assignment_df$Module <- paste0(
    "M",
    as.numeric(
      gsub(
        "[^0-9]",
        "",
        module_assignment_df$Module
      )
    )
  )
  
  
  duplicated_assignment <- module_assignment_df %>%
    count(
      Module
    ) %>%
    filter(
      n > 1
    )
  
  
  if (
    nrow(duplicated_assignment) > 0
  ) {
    
    cat(
      "\n❌ module_assignment 中存在重复 module：\n"
    )
    
    print(
      duplicated_assignment
    )
    
    stop(
      "同一个 module 被分配到了多个 Assigned_Cell，请先修正。"
    )
  }
  
  
  expected_modules <- paste0(
    "M",
    seq_len(n_modules)
  )
  
  
  missing_assignment <- setdiff(
    expected_modules,
    module_assignment_df$Module
  )
  
  
  extra_assignment <- setdiff(
    module_assignment_df$Module,
    expected_modules
  )
  
  
  if (
    length(missing_assignment) > 0
  ) {
    
    cat(
      "\n❌ 没有 Assigned_Cell 的 module：\n"
    )
    
    print(
      missing_assignment
    )
    
    stop(
      "module_assignment 不完整。"
    )
  }
  
  
  if (
    length(extra_assignment) > 0
  ) {
    
    cat(
      "\n❌ module_assignment 中存在不存在的 module：\n"
    )
    
    print(
      extra_assignment
    )
    
    stop(
      "module_assignment 包含额外 module。"
    )
  }
  
  
  cat(
    "\n✅ module_assignment 完整：",
    length(unique(module_assignment_df$Module)),
    "个 module\n"
  )
  
  
  # ==========================================================
  # 3. 正确对齐 TOM 和 moduleColors
  #
  # 不再使用之前危险的：
  #
  # names(moduleColors) <-
  #   colnames(tom)[match(...)]
  #
  # ==========================================================
  
  if (
    is.null(names(moduleColors))
  ) {
    
    stop(
      "moduleColors 必须有 gene names。"
    )
  }
  
  
  if (
    is.null(rownames(tom)) ||
    is.null(colnames(tom))
  ) {
    
    stop(
      "TOM 必须具有 gene rownames/colnames。"
    )
  }
  
  
  common_genes <- intersect(
    names(moduleColors),
    rownames(tom)
  )
  
  
  cat(
    "\nTOM genes:",
    nrow(tom),
    "\n"
  )
  
  cat(
    "moduleColors genes:",
    length(moduleColors),
    "\n"
  )
  
  cat(
    "common genes:",
    length(common_genes),
    "\n"
  )
  
  
  if (
    length(common_genes) == 0
  ) {
    
    stop(
      "TOM 与 moduleColors 完全没有共同 gene。"
    )
  }
  
  
  moduleColors_use <- moduleColors[
    common_genes
  ]
  
  
  tom_use <- tom[
    common_genes,
    common_genes,
    drop = FALSE
  ]
  
  
  # ==========================================================
  # 4. moduleColors：
  #
  # 原始 M78/M5/... -> 新 M1/M2/...
  # ==========================================================
  
  raw_gene_modules <- normalize_module_id(
    moduleColors_use
  )
  
  
  gene_modules <- unname(
    module_map[
      raw_gene_modules
    ]
  )
  
  
  names(gene_modules) <- names(
    moduleColors_use
  )
  
  
  # 不属于目标 module 的 gene 删除
  valid_gene_module <- !is.na(
    gene_modules
  )
  
  
  gene_modules <- gene_modules[
    valid_gene_module
  ]
  
  
  cat(
    "\nGenes per NEW module:\n"
  )
  
  print(
    table(
      factor(
        gene_modules,
        levels = expected_modules
      )
    )
  )
  
  
  # ==========================================================
  # 5. datKME：
  #
  # 原始 kME_M78 -> 新 kME_M1
  # 原始 kME_M5  -> 新 kME_M2
  # ...
  # ==========================================================
  
  datKME_use <- datKME
  
  
  raw_kme_module <- normalize_module_id(
    colnames(datKME_use)
  )
  
  
  mapped_kme_module <- unname(
    module_map[
      raw_kme_module
    ]
  )
  
  
  keep_kme <- !is.na(
    mapped_kme_module
  )
  
  
  datKME_use <- datKME_use[
    ,
    keep_kme,
    drop = FALSE
  ]
  
  
  mapped_kme_module <- mapped_kme_module[
    keep_kme
  ]
  
  
  colnames(datKME_use) <- paste0(
    "kME_",
    mapped_kme_module
  )
  
  
  cat(
    "\nConverted datKME columns:\n"
  )
  
  print(
    colnames(datKME_use)
  )
  
  
  # ==========================================================
  # 6. Gene info
  # ==========================================================
  
  cat(
    "\n2. 对齐 kME 与 supervision\n"
  )
  
  
  gene_info <- data.frame(
    gene = names(gene_modules),
    module = unname(gene_modules),
    stringsAsFactors = FALSE
  )
  
  
  gene_info <- gene_info %>%
    filter(
      !is.na(module)
    )
  
  
  # ----------------------------------------------------------
  # 每个 gene 取其所属 NEW module 对应的 kME
  # ----------------------------------------------------------
  
  gene_info$kME <- vapply(
    seq_len(nrow(gene_info)),
    function(i) {
      
      g <- gene_info$gene[i]
      
      m <- gene_info$module[i]
      
      col_name <- paste0(
        "kME_",
        m
      )
      
      
      if (
        g %in% rownames(datKME_use) &&
        col_name %in% colnames(datKME_use)
      ) {
        
        return(
          as.numeric(
            datKME_use[
              g,
              col_name
            ]
          )
        )
        
      } else {
        
        return(
          NA_real_
        )
      }
    },
    numeric(1)
  )
  
  
  cat(
    "\nGenes with valid kME:",
    sum(!is.na(gene_info$kME)),
    "/",
    nrow(gene_info),
    "\n"
  )
  
  
  gene_info <- gene_info %>%
    filter(
      !is.na(kME)
    )
  
  
  # ==========================================================
  # 7. 加 supervision
  # ==========================================================
  
  gene_info <- gene_info %>%
    left_join(
      module_assignment_df[
        ,
        c(
          "Module",
          "Assigned_Cell"
        )
      ],
      by = c(
        "module" = "Module"
      )
    )
  
  
  if (
    any(is.na(gene_info$Assigned_Cell))
  ) {
    
    stop(
      "有 module 没有匹配到 Assigned_Cell。"
    )
  }
  
  
  cat(
    "\nModules retained after kME + supervision:\n"
  )
  
  print(
    table(
      factor(
        gene_info$module,
        levels = expected_modules
      )
    )
  )
  
  
  # ==========================================================
  # 8. 每个 module 选 hub / Top genes
  # ==========================================================
  
  cat(
    "\n3. 选择 hub genes 和 Top genes\n"
  )
  
  
  hub_genes <- gene_info %>%
    group_by(
      module
    ) %>%
    slice_max(
      order_by = kME,
      n = n_hubs,
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    pull(
      gene
    )
  
  
  selected_genes <- gene_info %>%
    group_by(
      module
    ) %>%
    slice_max(
      order_by = kME,
      n = n_genes_per_mod,
      with_ties = FALSE
    ) %>%
    ungroup()
  
  
  valid_selected <- intersect(
    selected_genes$gene,
    rownames(tom_use)
  )
  
  
  valid_hubs <- intersect(
    hub_genes,
    colnames(tom_use)
  )
  
  
  selected_genes <- selected_genes %>%
    filter(
      gene %in% valid_selected
    )
  
  
  cat(
    "\nSelected genes per module:\n"
  )
  
  print(
    table(
      factor(
        selected_genes$module,
        levels = expected_modules
      )
    )
  )
  
  
  retained_modules <- unique(
    selected_genes$module
  )
  
  
  if (
    length(retained_modules) != n_modules
  ) {
    
    missing_after_selection <- setdiff(
      expected_modules,
      retained_modules
    )
    
    cat(
      "\n❌ 以下 module 在选择基因以后消失：\n"
    )
    
    print(
      missing_after_selection
    )
    
    stop(
      "部分 module 没有有效 selected genes。"
    )
  }
  
  
  # ==========================================================
  # 9. Supervised UMAP
  # ==========================================================
  
  feature_mat <- tom_use[
    valid_selected,
    valid_hubs,
    drop = FALSE
  ]
  
  
  y_labels <- selected_genes$Assigned_Cell[
    match(
      valid_selected,
      selected_genes$gene
    )
  ]
  
  
  set.seed(
    42
  )
  
  
  hub_umap <- uwot::umap(
    X = feature_mat,
    y = as.factor(y_labels),
    n_neighbors = 15,
    min_dist = 0.4,
    metric = "cosine",
    spread = 4
  )
  
  
  plot_df <- as.data.frame(
    hub_umap
  )
  
  
  colnames(plot_df) <- c(
    "UMAP1",
    "UMAP2"
  )
  
  
  plot_df$gene <- valid_selected
  
  
  plot_df <- plot_df %>%
    left_join(
      selected_genes,
      by = "gene"
    )
  
  
  # ==========================================================
  # 10. 检查最终是不是所有 module 都存在
  # ==========================================================
  
  cat(
    "\nFINAL modules in plot_df:\n"
  )
  
  
  print(
    table(
      factor(
        plot_df$module,
        levels = expected_modules
      )
    )
  )
  
  
  missing_final <- setdiff(
    expected_modules,
    unique(plot_df$module)
  )
  
  
  if (
    length(missing_final) > 0
  ) {
    
    cat(
      "\n❌ FINAL missing modules:\n"
    )
    
    print(
      missing_final
    )
    
    stop(
      "最终 plot_df 仍有 module 缺失。"
    )
  }
  
  
  # ==========================================================
  # 11. kME point size
  # ==========================================================
  
  kme_min <- min(
    plot_df$kME,
    na.rm = TRUE
  )
  
  
  kme_max <- max(
    plot_df$kME,
    na.rm = TRUE
  )
  
  
  if (
    kme_max == kme_min
  ) {
    
    plot_df$kME_scaled <- 0.5
    
  } else {
    
    plot_df$kME_scaled <- (
      plot_df$kME - kme_min
    ) / (
      kme_max - kme_min
    )
  }
  
  
  plot_df$module <- factor(
    plot_df$module,
    levels = expected_modules
  )
  
  
  # ==========================================================
  # 12. Edges
  #
  # edge_threshold / edge_sample 只影响边
  # 不再影响 module
  # ==========================================================
  
  sub_tom <- tom_use[
    valid_selected,
    valid_selected,
    drop = FALSE
  ]
  
  
  edge_df <- as.data.frame(
    as.table(
      as.matrix(
        sub_tom
      )
    )
  )
  
  
  colnames(edge_df) <- c(
    "from",
    "to",
    "weight"
  )
  
  
  edge_df <- edge_df %>%
    filter(
      as.character(from) < as.character(to),
      weight > edge_threshold
    )
  
  
  if (
    nrow(edge_df) > 0 &&
    edge_sample < 1
  ) {
    
    set.seed(
      42
    )
    
    edge_df <- edge_df %>%
      sample_frac(
        edge_sample
      )
  }
  
  
  gene_to_mod <- setNames(
    as.character(plot_df$module),
    plot_df$gene
  )
  
  
  edge_df$mod1 <- gene_to_mod[
    as.character(edge_df$from)
  ]
  
  
  edge_df$mod2 <- gene_to_mod[
    as.character(edge_df$to)
  ]
  
  
  edge_df$color <- ifelse(
    edge_df$mod1 == edge_df$mod2,
    edge_df$mod1,
    "grey90"
  )
  
  
  edge_df <- edge_df %>%
    left_join(
      dplyr::select(
        plot_df,
        gene,
        x = UMAP1,
        y = UMAP2
      ),
      by = c(
        "from" = "gene"
      )
    ) %>%
    left_join(
      dplyr::select(
        plot_df,
        gene,
        xend = UMAP1,
        yend = UMAP2
      ),
      by = c(
        "to" = "gene"
      )
    )
  
  
  if (
    nrow(edge_df) > 0
  ) {
    
    edge_df <- edge_df %>%
      group_by(
        color
      ) %>%
      mutate(
        alpha_val =
          (
            weight - min(weight)
          ) /
          (
            max(weight) -
              min(weight) +
              1e-6
          )
      ) %>%
      ungroup()
  }
  
  
  # ==========================================================
  # 13. 每个 module Top3 标签
  # ==========================================================
  
  hubs_to_label <- plot_df %>%
    group_by(
      module
    ) %>%
    slice_max(
      order_by = kME,
      n = 3,
      with_ties = FALSE
    ) %>%
    ungroup()
  
  
  # ==========================================================
  # 14. 固定颜色
  # ==========================================================
  
  mod_palette <- setNames(
    module_color_library[
      seq_len(n_modules)
    ],
    expected_modules
  )
  
  
  edge_palette <- c(
    mod_palette,
    "grey90" = "grey90"
  )
  
  
  # ==========================================================
  # 15. Plot
  # ==========================================================
  
  p_final <- ggplot()
  
  
  if (
    nrow(edge_df) > 0
  ) {
    
    p_final <- p_final +
      geom_segment(
        data = edge_df,
        aes(
          x = x,
          y = y,
          xend = xend,
          yend = yend,
          color = color,
          alpha = alpha_val
        ),
        linewidth = 0.15
      )
  }
  
  
  p_final <- p_final +
    
    geom_point(
      data = plot_df,
      aes(
        x = UMAP1,
        y = UMAP2,
        color = module,
        size = kME_scaled
      ),
      alpha = 0.9
    ) +
    
    geom_text_repel(
      data = hubs_to_label,
      aes(
        x = UMAP1,
        y = UMAP2,
        label = gene
      ),
      size = 7 / .pt,
      fontface = "plain",
      family = "Arial",
      color = "black",
      bg.color = "white",
      bg.r = 0.1,
      max.overlaps = Inf,
      segment.color = "grey50",
      segment.size = 0.2
    ) +
    
    scale_color_manual(
      values = edge_palette,
      breaks = expected_modules,
      drop = FALSE
    ) +
    
    scale_alpha_identity() +
    
    scale_size_continuous(
      range = c(
        0.3,
        1.5
      ),
      guide = "none"
    ) +
    
    coord_fixed() +
    
    theme_void() +
    
    theme(
      
      text = element_text(
        size = 7,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      
      legend.text = element_text(
        size = 7,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      
      legend.title = element_blank(),
      
      legend.position = "right",
      
      legend.key.size = unit(
        0.3,
        "cm"
      ),
      
      plot.margin = margin(
        10,
        10,
        10,
        10
      )
    ) +
    
    guides(
      color = guide_legend(
        override.aes = list(
          size = 2,
          alpha = 1
        ),
        ncol = 1
      )
    )
  
  
  ggsave(
    paste0(
      output_name,
      ".pdf"
    ),
    plot = p_final,
    width = 5.5,
    height = 4.5,
    device = cairo_pdf
  )
  
  
  ggsave(
    paste0(
      output_name,
      ".png"
    ),
    plot = p_final,
    width = 5.5,
    height = 4.5,
    dpi = 600
  )
  
  
  cat(
    "\n🎉 完成：",
    output_name,
    "\n"
  )
  
  
  return(
    invisible(
      list(
        plot = p_final,
        mapping = mapping_df,
        plot_df = plot_df,
        gene_info = gene_info
      )
    )
  )
}



#bayescor1
tom<-seurat_obj@misc$MEscnet$bayes_cor1


seurat_obj<-ComputeMEscnetModules(seurat_obj,min_genes_per_module = 20,resolution = 2,number = 2)
modules1<-seurat_obj@misc$MEscnet$MEscnet_modules$module_ids

module_order1_raw <- modules1
module_genes<-seurat_obj@misc$MEscnet$MEscnet_modules$gene_lists
moduleColors <- seurat_obj@misc$MEscnet$moduleColors
MEs0 <- moduleEigengenes(
  datExpr,
  colors = moduleColors
)$eigengenes

MEs <- orderMEs(MEs0)
datKME <- signedKME(
  datExpr,
  MEs,
  outputColumnName="kME_MM."
)



# 1. 根据你的参考图片，手动建立严格的监督映射字典
proCSC_modules <- c("M11", "M15", "M4", "M18", "M22", "M6", "M8", "M3", "M23", "M10", "M24")
CSC_modules    <- c("M2", "M14", "M16")
revCSC_modules <- c("M1", "M7", "M17", "M13", "M20", "M9", "M12", "M21", "M25", "M5", "M19")

# 2. 将字典转化为标准的 dataframe
module_assignment <- data.frame(
  Module = c(proCSC_modules, CSC_modules, revCSC_modules),
  Assigned_Cell = c(
    rep("proCSC", length(proCSC_modules)),
    rep("CSC", length(CSC_modules)),
    rep("revCSC", length(revCSC_modules))
  ),
  stringsAsFactors = FALSE
)


# 打印检查一下，确保 25 个模块一个不少
print(table(module_assignment$Assigned_Cell))
# 假设你之前算的 module_assignment (模块归属关系表) 还在环境中
plot_supervised_module_network(
  tom = tom, 
  datKME = datKME, 
  moduleColors = moduleColors, 
  module_assignment_df = module_assignment, # 包含 Module 和 Assigned_Cell 两列
  output_name = "Fig_Supervised1",
  module_order_raw = module_order1_raw
)


#bayescor2
tom<-seurat_obj@misc$MEscnet$bayes_cor2



resoltuion=2.9
seurat_obj<-ComputeMEscnetModules(seurat_obj,min_genes_per_module = 20,resolution = resoltuion,number = 3)
modules2<-seurat_obj@misc$MEscnet$MEscnet_modules$module_ids


module_order2_raw <- modules2
module2_genes<-seurat_obj@misc$MEscnet$MEscnet_modules$gene_lists
moduleColors <- seurat_obj@misc$MEscnet$moduleColors
MEs0 <- moduleEigengenes(
  datExpr,
  colors = moduleColors
)$eigengenes

MEs <- orderMEs(MEs0)
datKME <- signedKME(
  datExpr,
  MEs,
  outputColumnName="kME_MM."
)



# 1. 根据你的参考图片，手动建立严格的监督映射字典
proCSC_modules <- c("M1", "M22", "M10", "M12", "M17", "M7", "M6", "M16", "M24", "M9", "M19")
CSC_modules    <- c("M2", "M20")
revCSC_modules <- c("M4", "M13", "M5", "M23", "M3", "M21", "M8", "M11", "M15","M18","14")


# 2. 将字典转化为标准的 dataframe
module_assignment <- data.frame(
  Module = c(proCSC_modules, CSC_modules, revCSC_modules),
  Assigned_Cell = c(
    rep("proCSC", length(proCSC_modules)),
    rep("CSC", length(CSC_modules)),
    rep("revCSC", length(revCSC_modules))
  ),
  stringsAsFactors = FALSE
)


# 打印检查一下，确保 25 个模块一个不少
print(table(module_assignment$Assigned_Cell))
# 假设你之前算的 module_assignment (模块归属关系表) 还在环境中
plot_supervised_module_network(
  tom = tom, 
  datKME = datKME, 
  moduleColors = moduleColors, 
  module_order_raw = module_order2_raw,
  module_assignment_df = module_assignment, # 包含 Module 和 Assigned_Cell 两列
  output_name = "Fig_Supervised2"
)






#bayescor3
tom<-seurat_obj@misc$MEscnet$bayes_cor3
resoltuion=2.9
seurat_obj<-ComputeMEscnetModules(seurat_obj,min_genes_per_module = 20,resolution = resoltuion,number = 4)
modules3<-seurat_obj@misc$MEscnet$MEscnet_modules$module_ids
module3_genes<-seurat_obj@misc$MEscnet$MEscnet_modules$gene_lists


module_order3_raw <- modules3

moduleColors <- seurat_obj@misc$MEscnet$moduleColors
MEs0 <- moduleEigengenes(
  datExpr,
  colors = moduleColors
)$eigengenes

MEs <- orderMEs(MEs0)
datKME <- signedKME(
  datExpr,
  MEs,
  outputColumnName="kME_MM."
)



# 1. 根据你的参考图片，手动建立严格的监督映射字典
proCSC_modules <- c("M13", "M6", "M9", "M16", "M14", "M2", "M3", "M15")
CSC_modules    <- c("M12", "M24")
revCSC_modules <- c("M1", "M21", "M17", "M10", "M19", "M22", "M18", "M7", "M20", "M23", "M11","M5","M4","M8")


# 2. 将字典转化为标准的 dataframe
module_assignment <- data.frame(
  Module = c(proCSC_modules, CSC_modules, revCSC_modules),
  Assigned_Cell = c(
    rep("proCSC", length(proCSC_modules)),
    rep("CSC", length(CSC_modules)),
    rep("revCSC", length(revCSC_modules))
  ),
  stringsAsFactors = FALSE
)


# 打印检查一下，确保 25 个模块一个不少
print(table(module_assignment$Assigned_Cell))
# 假设你之前算的 module_assignment (模块归属关系表) 还在环境中
plot_supervised_module_network(
  tom = tom, 
  datKME = datKME, 
  moduleColors = moduleColors, 
  module_order_raw = module_order3_raw,
  module_assignment_df = module_assignment, # 包含 Module 和 Assigned_Cell 两列
  output_name = "Fig_Supervised3"
)














old_sigs <- c(
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

new_sigs <- paste0("bayescor1_module", 1:25)

col_indices <- match(old_sigs, colnames(seurat_obj@meta.data))

if(all(!is.na(col_indices))) {
  colnames(seurat_obj@meta.data)[col_indices] <- new_sigs
  cat("Done!\n")
}

my_sigs <- new_sigs

















resoltuion=2.9
seurat_obj<-ComputeMEscnetModules(seurat_obj,min_genes_per_module = 20,resolution = resoltuion,number = 3)
modules2<-seurat_obj@misc$MEscnet$MEscnet_modules$module_ids
module2_genes<-seurat_obj@misc$MEscnet$MEscnet_modules$gene_lists
str(module2_genes)




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


################################
# 5. 生成 module-gene-kME list
################################
module2_kME_list <- list()

for(mod in modules){
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
  module2_kME_list[[mod]] <- tmp
}









resoltuion=2.9
seurat_obj<-ComputeMEscnetModules(seurat_obj,min_genes_per_module = 20,resolution = resoltuion,number = 4)
modules3<-seurat_obj@misc$MEscnet$MEscnet_modules$module_ids
module3_genes<-seurat_obj@misc$MEscnet$MEscnet_modules$gene_lists
str(module3_genes)



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


################################
# 5. 生成 module-gene-kME list
################################
module3_kME_list <- list()

for(mod in modules){
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
  module3_kME_list[[mod]] <- tmp
}
################################
# 6. 查看






library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)

windowsFonts(
  Arial = windowsFont("Arial")
)

#fig3
plot_OR_heatmap_twoPDF <- function(
    list1,
    list2,
    name1 = "BayesCor1",
    name2 = "BayesCor2",
    heatmap_pdf = "OR_Heatmap_Text.pdf",
    colorbar_pdf = "OR_Heatmap_ColorBars.pdf"
) {
  
  # ==========================================================
  # 1. Basic information
  # ==========================================================
  
  all_universe <- unique(
    c(
      unlist(list1),
      unlist(list2)
    )
  )
  
  N <- length(all_universe)
  
  n1 <- length(list1)
  n2 <- length(list2)
  
  mod1_names <- paste0(
    "M",
    seq_len(n1)
  )
  
  mod2_names <- paste0(
    "M",
    seq_len(n2)
  )
  
  cat(
    "\nNumber of ",
    name1,
    " modules: ",
    n1,
    "\n",
    sep = ""
  )
  
  cat(
    "Number of ",
    name2,
    " modules: ",
    n2,
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 2. OR + Fisher exact test
  # ==========================================================
  
  or_mat <- matrix(
    0,
    nrow = n1,
    ncol = n2,
    dimnames = list(
      mod1_names,
      mod2_names
    )
  )
  
  p_mat <- matrix(
    1,
    nrow = n1,
    ncol = n2,
    dimnames = list(
      mod1_names,
      mod2_names
    )
  )
  
  
  for (i in seq_len(n1)) {
    
    g1 <- unique(
      list1[[i]]
    )
    
    for (j in seq_len(n2)) {
      
      g2 <- unique(
        list2[[j]]
      )
      
      q <- length(
        intersect(
          g1,
          g2
        )
      )
      
      if (q > 0) {
        
        ft <- fisher.test(
          matrix(
            c(
              q,
              length(g1) - q,
              length(g2) - q,
              N - length(g1) - length(g2) + q
            ),
            nrow = 2
          ),
          alternative = "greater"
        )
        
        or_mat[i, j] <- unname(
          ft$estimate
        )
        
        p_mat[i, j] <- ft$p.value
      }
    }
  }
  
  
  # ==========================================================
  # 3. Greedy column matching
  #
  # IMPORTANT:
  # 不再使用 rep(NA, n1)
  #
  # 最终 matched_cols 长度严格等于 n2
  # ==========================================================
  
  matched_cols <- character(0)
  
  avail_cols <- mod2_names
  
  for (i in seq_len(n1)) {
    
    # 所有 column 已经匹配完以后直接停止
    if (length(avail_cols) == 0) {
      break
    }
    
    row_name <- mod1_names[i]
    
    vals <- or_mat[
      row_name,
      avail_cols,
      drop = TRUE
    ]
    
    best_col <- avail_cols[
      which.max(vals)
    ]
    
    matched_cols <- c(
      matched_cols,
      best_col
    )
    
    avail_cols <- setdiff(
      avail_cols,
      best_col
    )
  }
  
  
  # 理论上如果还有列没匹配，追加
  if (length(avail_cols) > 0) {
    
    matched_cols <- c(
      matched_cols,
      avail_cols
    )
  }
  
  
  # 强制安全检查
  stopifnot(
    length(matched_cols) == n2
  )
  
  stopifnot(
    setequal(
      matched_cols,
      mod2_names
    )
  )
  
  
  cat(
    "\nFinal column order:\n"
  )
  
  print(
    matched_cols
  )
  
  
  # ==========================================================
  # 4. Long-format data
  # ==========================================================
  
  df_plot <- reshape2::melt(
    or_mat,
    varnames = c(
      "Mod1",
      "Mod2"
    ),
    value.name = "OddsRatio"
  )
  
  df_plot$Pval <- reshape2::melt(
    p_mat
  )$value
  
  df_plot$OR_capped <- pmin(
    df_plot$OddsRatio,
    50
  )
  
  
  df_plot$Signif <- ifelse(
    df_plot$Pval < 0.001,
    "***",
    
    ifelse(
      df_plot$Pval < 0.01,
      "**",
      
      ifelse(
        df_plot$Pval < 0.05,
        "*",
        ""
      )
    )
  )
  
  
  # ==========================================================
  # 5. Factor ordering
  # ==========================================================
  
  # M1 在最上方
  df_plot$Mod1 <- factor(
    df_plot$Mod1,
    levels = rev(
      mod1_names
    )
  )
  
  # matched columns
  df_plot$Mod2 <- factor(
    df_plot$Mod2,
    levels = matched_cols
  )
  
  
  # ==========================================================
  # 6. PDF 1
  # OR heatmap + M1/M2/... labels
  # ==========================================================
  
  p_text <- ggplot(
    df_plot,
    aes(
      x = Mod2,
      y = Mod1,
      fill = OR_capped
    )
  ) +
    
    geom_tile(
      color = "white",
      linewidth = 0.5
    ) +
    
    geom_text(
      aes(
        label = Signif
      ),
      color = "black",
      size = 14 / .pt,
      family = "Arial",
      vjust = 0.75
    ) +
    
    scale_fill_gradient(
      low = "#FFFFFF",
      high = "#CB181D",
      name = "Odds ratio",
      limits = c(
        0,
        50
      )
    ) +
    
    labs(
      x = paste0(
        name2,
        " modules"
      ),
      y = paste0(
        name1,
        " modules"
      )
    ) +
    
    coord_fixed() +
    
    theme_classic() +
    
    theme(
      
      text = element_text(
        size = 14,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      
      axis.text.x = element_text(
        size = 14,
        color = "black",
        face = "plain",
        family = "Arial",
        angle = 45,
        hjust = 1
      ),
      
      axis.text.y = element_text(
        size = 14,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      
      axis.title.x = element_text(
        size = 14,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      
      axis.title.y = element_text(
        size = 14,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      
      legend.text = element_text(
        size = 14,
        color = "black",
        family = "Arial"
      ),
      
      legend.title = element_text(
        size = 14,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      
      axis.ticks = element_blank(),
      
      panel.grid = element_blank(),
      
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.5
      )
    )
  
  
  ggsave(
    heatmap_pdf,
    plot = p_text,
    width = 12,
    height = 10,
    device = cairo_pdf
  )
  
  
  # ==========================================================
  # 7. Prepare module colors
  # ==========================================================
  
  if (
    n1 > length(module_color_library) ||
    n2 > length(module_color_library)
  ) {
    
    stop(
      "Number of modules exceeds module_color_library."
    )
  }
  
  
  # --------------------------
  # Row module colors
  # --------------------------
  
  row_colors <- setNames(
    module_color_library[
      seq_len(n1)
    ],
    mod1_names
  )
  
  
  # --------------------------
  # Column module colors
  #
  # 注意：
  # M1 永远是颜色1
  # M2 永远是颜色2
  #
  # 只是根据 matched_cols 重排
  # --------------------------
  
  col_colors_all <- setNames(
    module_color_library[
      seq_len(n2)
    ],
    mod2_names
  )
  
  col_colors <- col_colors_all[
    matched_cols
  ]
  
  
  # ==========================================================
  # 8. PDF 2
  #
  # ONLY module color bars
  #
  # 重点：
  #
  # 中间白色区域严格为
  #
  #       n2 columns × n1 rows
  #
  # 左边 = 1 column
  # 底部 = 1 row
  #
  # 因此每一个色块的物理尺寸与 heatmap cell 完全一致
  # ==========================================================
  
  pdf(
    colorbar_pdf,
    width = 12,
    height = 10
  )
  
  
  grid.newpage()
  
  
  # ----------------------------------------------------------
  # Grid layout
  #
  #     left bar | heatmap-sized blank area
  #     ---------|---------------------------
  #              | bottom bar
  #
  # widths:
  # 1 cell + n2 cells
  #
  # heights:
  # n1 cells + 1 cell
  # ----------------------------------------------------------
  
  lay <- grid.layout(
    nrow = 2,
    ncol = 2,
    
    widths = unit(
      c(
        1,
        n2
      ),
      "null"
    ),
    
    heights = unit(
      c(
        n1,
        1
      ),
      "null"
    )
  )
  
  
  pushViewport(
    viewport(
      layout = lay
    )
  )
  
  
  # ==========================================================
  # 8A. Left vertical bar
  # ==========================================================
  
  pushViewport(
    viewport(
      layout.pos.row = 1,
      layout.pos.col = 1,
      
      xscale = c(
        0,
        1
      ),
      
      yscale = c(
        0,
        n1
      )
    )
  )
  
  
  # M1 在最上面
  for (i in seq_len(n1)) {
    
    y_bottom <- n1 - i
    
    grid.rect(
      x = unit(
        0.5,
        "npc"
      ),
      
      y = unit(
        (y_bottom + 0.5) / n1,
        "npc"
      ),
      
      width = unit(
        1,
        "npc"
      ),
      
      height = unit(
        1 / n1,
        "npc"
      ),
      
      gp = gpar(
        fill = row_colors[
          mod1_names[i]
        ],
        col = NA
      )
    )
  }
  
  
  popViewport()
  
  
  # ==========================================================
  # 8B. Bottom horizontal bar
  # ==========================================================
  
  pushViewport(
    viewport(
      layout.pos.row = 2,
      layout.pos.col = 2,
      
      xscale = c(
        0,
        n2
      ),
      
      yscale = c(
        0,
        1
      )
    )
  )
  
  
  for (j in seq_len(n2)) {
    
    module_name <- matched_cols[j]
    
    grid.rect(
      x = unit(
        (j - 0.5) / n2,
        "npc"
      ),
      
      y = unit(
        0.5,
        "npc"
      ),
      
      width = unit(
        1 / n2,
        "npc"
      ),
      
      height = unit(
        1,
        "npc"
      ),
      
      gp = gpar(
        fill = col_colors[
          module_name
        ],
        col = NA
      )
    )
  }
  
  
  popViewport()
  
  
  # ==========================================================
  # 8C. Optional blank central panel
  #
  # 仅画透明/白色，不画边框
  # ==========================================================
  
  pushViewport(
    viewport(
      layout.pos.row = 1,
      layout.pos.col = 2
    )
  )
  
  grid.rect(
    gp = gpar(
      fill = NA,
      col = NA
    )
  )
  
  popViewport()
  
  
  popViewport()
  
  dev.off()
  
  
  # ==========================================================
  # 9. Final checks
  # ==========================================================
  
  cat(
    "\n====================================\n"
  )
  
  cat(
    "PDF 1:",
    heatmap_pdf,
    "\n"
  )
  
  cat(
    "PDF 2:",
    colorbar_pdf,
    "\n"
  )
  
  cat(
    "Rows:",
    n1,
    "\n"
  )
  
  cat(
    "Columns:",
    n2,
    "\n"
  )
  
  cat(
    "Bottom colors:",
    length(col_colors),
    "\n"
  )
  
  cat(
    "Left colors:",
    length(row_colors),
    "\n"
  )
  
  cat(
    "====================================\n"
  )
  
  
  return(
    invisible(
      list(
        plot = p_text,
        OR = or_mat,
        P = p_mat,
        matched_cols = matched_cols,
        row_colors = row_colors,
        col_colors = col_colors
      )
    )
  )
}




res <- plot_OR_heatmap_twoPDF(
  module2_genes,
  module3_genes,
  name1 = "BayesCor1",
  name2 = "BayesCor2",
  heatmap_pdf = "Fig_OR_Text_4.pdf",
  colorbar_pdf = "Fig_OR_ColorBars_4.pdf"
)

#fig4





library(dplyr)
library(UCell)
library(Seurat)
library(corrplot)
library(ggplot2)
library(ggrastr)
library(reshape2)

windowsFonts(Arial = windowsFont("Arial"))
library(dplyr)
library(UCell)
library(Seurat)
library(corrplot)
library(ggplot2)

windowsFonts(Arial = windowsFont("Arial"))

# ==============================================================================
# 1. 正确处理 module：
#    不按 gene number 排序
#    保持 module_kME_list 原始顺序
#    只重新连续编号
# ==============================================================================

process_module_list <- function(mod_list, prefix) {
  
  original_names <- names(mod_list)
  
  # 保持原顺序，只过滤 kME > 0.6
  filtered <- lapply(
    mod_list,
    function(df) {
      df %>%
        filter(kME > 0.6) %>%
        pull(gene)
    }
  )
  
  # 只连续重新编号，不排序
  new_names <- paste0(
    prefix,
    "_Module",
    seq_along(filtered)
  )
  
  names(filtered) <- new_names
  
  # 打印映射检查
  mapping <- data.frame(
    Original = original_names,
    New = new_names,
    Ngenes = lengths(filtered)
  )
  
  cat("\n", prefix, " mapping:\n", sep = "")
  print(mapping)
  
  return(filtered)
}


# ==============================================================================
# 2. 重新生成正确顺序的 signatures
# ==============================================================================

bc1_sigs <- process_module_list(
  module_kME_list,
  "BayesCor1"
)

bc2_sigs <- process_module_list(
  module2_kME_list,
  "BayesCor2"
)

bc3_sigs <- process_module_list(
  module3_kME_list,
  "BayesCor3"
)


# ==============================================================================
# 3. 非常重要：
#    删除之前按错误顺序计算出来的旧 UCell columns
# ==============================================================================

old_bayes_cols <- grep(
  "^BayesCor[123]_Module[0-9]+_UCell$",
  colnames(seurat_obj@meta.data),
  value = TRUE
)

cat("\n旧 BayesCor UCell columns:\n")
print(old_bayes_cols)

if (length(old_bayes_cols) > 0) {
  seurat_obj@meta.data[
    ,
    old_bayes_cols
  ] <- NULL
}

cat(
  "\n旧 BayesCor UCell columns 删除后剩余数量:",
  length(
    grep(
      "^BayesCor[123]_Module[0-9]+_UCell$",
      colnames(seurat_obj@meta.data)
    )
  ),
  "\n"
)


# ==============================================================================
# 4. 重新计算 UCell
# ==============================================================================

all_new_sigs <- c(
  bc1_sigs,
  bc2_sigs,
  bc3_sigs
)

cat(
  "\n重新计算 UCell，总 signatures:",
  length(all_new_sigs),
  "\n"
)

seurat_obj <- UCell::AddModuleScore_UCell(
  obj = seurat_obj,
  features = all_new_sigs,
  ncores = 24,
  name = "_UCell"
)

cat("✅ 新 UCell 评分计算完毕！\n")


# ==============================================================================
# 5. 检查新列是否真的存在
# ==============================================================================

bc1_ucell_cols <- paste0(
  names(bc1_sigs),
  "_UCell"
)

bc2_ucell_cols <- paste0(
  names(bc2_sigs),
  "_UCell"
)

bc3_ucell_cols <- paste0(
  names(bc3_sigs),
  "_UCell"
)

cat("\nBC1 columns:\n")
print(bc1_ucell_cols)

cat("\nBC2 columns:\n")
print(bc2_ucell_cols)

cat("\nBC3 columns:\n")
print(bc3_ucell_cols)

stopifnot(
  all(
    c(
      bc1_ucell_cols,
      bc2_ucell_cols,
      bc3_ucell_cols
    ) %in% colnames(seurat_obj@meta.data)
  )
)


# ==============================================================================
# 6. public signatures
# ==============================================================================

all_cols <- colnames(seurat_obj@meta.data)

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
  "ImmatureDALERBA11_UCell",
  "OWNwntreceptors_UCell",
  "hKrasGSEA_UCell",
  "hTgfbKEGG_UCell",
  "hMapkKEGG_UCell",
  "hMapkGO_UCell",
  "hPi3kGO_UCell",
  "hiCMS2_UCell",
  "hiCMS3_UCell"
)

pub_sigs_final <- setdiff(
  pub_sigs_all,
  sigs_to_remove
)

pub_sigs_exist <- pub_sigs_final[
  pub_sigs_final %in% all_cols
]


# ==============================================================================
# 7. 热图函数
#    这里暂时保持你的原逻辑不动
# ==============================================================================

plot_cor_heatmap <- function(
    seurat_obj,
    my_sigs,
    pub_sigs_exist,
    out_filename
) {
  
  my_sigs_exist <- my_sigs[
    my_sigs %in% colnames(seurat_obj@meta.data)
  ]
  
  meta_sub <- seurat_obj@meta.data %>%
    filter(
      curatedCLUST %in% c(
        "CSC",
        "proCSC",
        "revCSC",
        "TA 1",
        "TA 2"
      )
    )
  
  # public signatures 单独聚类
  hc_pub <- hclust(
    as.dist(
      1 - cor(
        scale(
          as.matrix(
            meta_sub[, pub_sigs_exist, drop = FALSE]
          )
        ),
        method = "pearson"
      )
    ),
    method = "complete"
  )
  
  order_pub <- pub_sigs_exist[
    hc_pub$order
  ]
  
  # BayesCor modules 单独聚类
  hc_my <- hclust(
    as.dist(
      1 - cor(
        scale(
          as.matrix(
            meta_sub[, my_sigs_exist, drop = FALSE]
          )
        ),
        method = "pearson"
      )
    ),
    method = "complete"
  )
  
  order_my <- my_sigs_exist[
    hc_my$order
  ]
  
  cat("\nmy_sigs input:\n")
  print(my_sigs_exist)
  
  cat("\norder_my after hclust:\n")
  print(order_my)
  
  # 固定为：
  # BayesCor modules + public signatures
  final_order <- c(
    order_my,
    order_pub
  )
  
  mat_all <- scale(
    as.matrix(
      meta_sub[
        ,
        final_order,
        drop = FALSE
      ]
    )
  )
  
  cor_all <- cor(
    mat_all,
    method = "pearson"
  )
  
  res_all <- corrplot::cor.mtest(
    mat_all,
    method = "pearson"
  )
  
  cat("\nFinal matrix order check:\n")
  print(
    identical(
      final_order,
      colnames(cor_all)
    )
  )
  
  my_palette <- colorRampPalette(
    c(
      "#4393C3",
      "#FFFFFF",
      "#D6604D"
    )
  )(200)
  
  png(
    out_filename,
    width = 28,
    height = 28,
    units = "in",
    res = 600
  )
  
  par(
    family = "Arial",
    cex = 1.4,
    col.axis = "black",
    font = 1
  )
  
  corrplot::corrplot(
    cor_all,
    method = "circle",
    order = "original",
    tl.pos = "lt",
    tl.col = "black",
    tl.cex = 1,
    tl.srt = 45,
    cl.pos = "r",
    cl.cex = 1,
    cl.offset = 1,
    col = my_palette,
    p.mat = res_all$p,
    sig.level = 0.05,
    insig = "blank",
    mar = c(1, 1, 1, 1)
  )
  
  dev.off()
  
  cat(
    "✔ 成功生成热图:",
    out_filename,
    "\n"
  )
}


# ==============================================================================
# 8. 重新画
# ==============================================================================

plot_cor_heatmap(
  seurat_obj,
  bc1_ucell_cols,
  pub_sigs_exist,
  "Figure_Cor_Heatmap_BayesCor1_NEW.png"
)

plot_cor_heatmap(
  seurat_obj,
  bc2_ucell_cols,
  pub_sigs_exist,
  "Figure_Cor_Heatmap_BayesCor2_NEW.png"
)

plot_cor_heatmap(
  seurat_obj,
  bc3_ucell_cols,
  pub_sigs_exist,
  "Figure_Cor_Heatmap_BayesCor3_NEW.png"
)

#fig4_2

# ============================================================
# 2. Safely obtain BayesCor UCell columns
#
# Important:
# Extract module number specifically from "_ModuleXX_UCell"
# ============================================================

get_bayescor_cols <- function(
    seurat_obj,
    prefix
) {
  
  all_cols <- colnames(
    seurat_obj@meta.data
  )
  
  pattern <- paste0(
    "^",
    prefix,
    "_Module[0-9]+_UCell$"
  )
  
  x <- grep(
    pattern,
    all_cols,
    value = TRUE
  )
  
  # Extract ONLY module number
  module_number <- as.numeric(
    sub(
      paste0(
        "^",
        prefix,
        "_Module([0-9]+)_UCell$"
      ),
      "\\1",
      x
    )
  )
  
  x <- x[
    order(module_number)
  ]
  
  return(x)
}


bc1_cols <- get_bayescor_cols(
  seurat_obj,
  "BayesCor1"
)

bc2_cols <- get_bayescor_cols(
  seurat_obj,
  "BayesCor2"
)

bc3_cols <- get_bayescor_cols(
  seurat_obj,
  "BayesCor3"
)


cat("\nBayesCor1:\n")
print(bc1_cols)

cat("\nBayesCor2:\n")
print(bc2_cols)

cat("\nBayesCor3:\n")
print(bc3_cols)


# ============================================================
# 3. Main Pearson correlation plotting function
# ============================================================

plot_UCell_Circle_CorrectOrder <- function(
    seurat_obj,
    cols1,
    cols2,
    name1 = "BayesCor1",
    name2 = "BayesCor2",
    heatmap_pdf = "Cor_CircleHeatmap_Text.pdf",
    colorbar_pdf = "Cor_CircleHeatmap_ColorBars.pdf"
) {
  
  
  # ==========================================================
  # 3.1 Select cells
  # ==========================================================
  
  meta_sub <- seurat_obj@meta.data %>%
    filter(
      curatedCLUST %in% c(
        "CSC",
        "proCSC",
        "revCSC",
        "TA 1",
        "TA 2"
      )
    )
  
  
  # ==========================================================
  # 3.2 UCell score Z-score
  # ==========================================================
  
  mat1 <- scale(
    as.matrix(
      meta_sub[
        ,
        cols1,
        drop = FALSE
      ]
    )
  )
  
  mat2 <- scale(
    as.matrix(
      meta_sub[
        ,
        cols2,
        drop = FALSE
      ]
    )
  )
  
  
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  
  
  # ----------------------------------------------------------
  # Here M1/M2/... correspond to the corrected continuous
  # module identity.
  # ----------------------------------------------------------
  
  mod1_names <- paste0(
    "M",
    seq_len(n1)
  )
  
  mod2_names <- paste0(
    "M",
    seq_len(n2)
  )
  
  
  colnames(mat1) <- mod1_names
  colnames(mat2) <- mod2_names
  
  
  # ==========================================================
  # 3.3 Pearson correlation + P value
  # ==========================================================
  
  cor_mat <- matrix(
    NA_real_,
    nrow = n1,
    ncol = n2,
    dimnames = list(
      mod1_names,
      mod2_names
    )
  )
  
  p_mat <- matrix(
    NA_real_,
    nrow = n1,
    ncol = n2,
    dimnames = list(
      mod1_names,
      mod2_names
    )
  )
  
  
  for (i in seq_len(n1)) {
    
    for (j in seq_len(n2)) {
      
      test_res <- cor.test(
        mat1[, i],
        mat2[, j],
        method = "pearson"
      )
      
      cor_mat[i, j] <- unname(
        test_res$estimate
      )
      
      p_mat[i, j] <- test_res$p.value
    }
  }
  
  
  # ==========================================================
  # 3.4 Internal hclust
  #
  # This ONLY changes visualization order.
  # It no longer changes module identity.
  # ==========================================================
  
  hc1 <- hclust(
    as.dist(
      1 - cor(
        mat1,
        method = "pearson",
        use = "pairwise.complete.obs"
      )
    ),
    method = "complete"
  )
  
  ordered_mod1 <- mod1_names[
    hc1$order
  ]
  
  
  hc2 <- hclust(
    as.dist(
      1 - cor(
        mat2,
        method = "pearson",
        use = "pairwise.complete.obs"
      )
    ),
    method = "complete"
  )
  
  ordered_mod2 <- mod2_names[
    hc2$order
  ]
  
  
  # ==========================================================
  # 3.5 Keep your confirmed display direction
  #
  # BayesCor1 = normal
  # BayesCor2 / BayesCor3 = reverse
  # ==========================================================
  
  if (name1 != "BayesCor1") {
    
    ordered_mod1 <- rev(
      ordered_mod1
    )
  }
  
  
  if (name2 != "BayesCor1") {
    
    ordered_mod2 <- rev(
      ordered_mod2
    )
  }
  
  
  # ==========================================================
  # 3.6 Long format
  # ==========================================================
  
  df_cor <- reshape2::melt(
    cor_mat,
    varnames = c(
      "Mod1",
      "Mod2"
    ),
    value.name = "Correlation"
  )
  
  
  df_p <- reshape2::melt(
    p_mat,
    varnames = c(
      "Mod1",
      "Mod2"
    ),
    value.name = "Pval"
  )
  
  
  df_plot <- merge(
    df_cor,
    df_p,
    by = c(
      "Mod1",
      "Mod2"
    ),
    sort = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Keep original P < 0.05 display rule
  # ----------------------------------------------------------
  
  df_plot$Plot_Cor <- ifelse(
    df_plot$Pval < 0.05,
    df_plot$Correlation,
    NA_real_
  )
  
  
  # ----------------------------------------------------------
  # X:
  # left -> right = ordered_mod2
  #
  # Y:
  # top -> bottom = ordered_mod1
  # ----------------------------------------------------------
  
  df_plot$Mod1 <- factor(
    df_plot$Mod1,
    levels = rev(
      ordered_mod1
    )
  )
  
  df_plot$Mod2 <- factor(
    df_plot$Mod2,
    levels = ordered_mod2
  )
  
  
  # ==========================================================
  # 4. PDF 1: Text version
  # ==========================================================
  
  p_text <- ggplot(
    df_plot,
    aes(
      x = Mod2,
      y = Mod1
    )
  ) +
    
    geom_tile(
      fill = "white",
      color = "grey90",
      linewidth = 0.2
    ) +
    
    geom_point(
      aes(
        size = abs(Plot_Cor),
        fill = Plot_Cor
      ),
      shape = 21,
      color = "transparent"
    ) +
    
    scale_fill_gradient2(
      low = "#4393C3",
      mid = "#FFFFFF",
      high = "#D6604D",
      midpoint = 0,
      limits = c(
        -1,
        1
      ),
      na.value = "transparent",
      name = "Pearson"
    ) +
    
    scale_size_continuous(
      range = c(
        0,
        10
      ),
      limits = c(
        0,
        1
      ),
      guide = "none"
    ) +
    
    labs(
      x = paste0(
        name2,
        " modules"
      ),
      y = paste0(
        name1,
        " modules"
      )
    ) +
    
    coord_fixed() +
    
    theme_minimal() +
    
    theme(
      
      text = element_text(
        size = 14,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      
      axis.text.x = element_text(
        size = 14,
        color = "black",
        face = "plain",
        family = "Arial",
        angle = 45,
        hjust = 1
      ),
      
      axis.text.y = element_text(
        size = 14,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      
      axis.title = element_text(
        size = 14,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      
      legend.text = element_text(
        size = 14,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      
      legend.title = element_text(
        size = 14,
        color = "black",
        face = "plain",
        family = "Arial"
      ),
      
      panel.grid = element_blank(),
      
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.5
      )
    )
  
  
  ggsave(
    heatmap_pdf,
    plot = p_text,
    width = 12,
    height = 10,
    device = cairo_pdf
  )
  
  
  # ==========================================================
  # 5. Module color mapping
  #
  # M1 ALWAYS color 1
  # M2 ALWAYS color 2
  # ...
  #
  # hclust only changes POSITION.
  # ==========================================================
  
  if (
    n1 > length(module_color_library) ||
    n2 > length(module_color_library)
  ) {
    
    stop(
      "Number of modules exceeds module_color_library."
    )
  }
  
  
  row_colors_all <- setNames(
    module_color_library[
      seq_len(n1)
    ],
    mod1_names
  )
  
  
  col_colors_all <- setNames(
    module_color_library[
      seq_len(n2)
    ],
    mod2_names
  )
  
  
  # ==========================================================
  # 6. PDF 2: color bars only
  # ==========================================================
  
  pdf(
    colorbar_pdf,
    width = 12,
    height = 10
  )
  
  
  grid.newpage()
  
  
  lay <- grid.layout(
    
    nrow = 2,
    ncol = 2,
    
    widths = unit(
      c(
        1,
        n2
      ),
      "null"
    ),
    
    heights = unit(
      c(
        n1,
        1
      ),
      "null"
    )
  )
  
  
  pushViewport(
    viewport(
      layout = lay
    )
  )
  
  
  # ==========================================================
  # 6A. LEFT COLOR BAR
  #
  # Top -> bottom = ordered_mod1
  # ==========================================================
  
  pushViewport(
    viewport(
      layout.pos.row = 1,
      layout.pos.col = 1
    )
  )
  
  
  for (i in seq_len(n1)) {
    
    module_name <- ordered_mod1[i]
    
    y_bottom <- n1 - i
    
    grid.rect(
      
      x = unit(
        0.5,
        "npc"
      ),
      
      y = unit(
        (y_bottom + 0.5) / n1,
        "npc"
      ),
      
      width = unit(
        1,
        "npc"
      ),
      
      height = unit(
        1 / n1,
        "npc"
      ),
      
      gp = gpar(
        fill = row_colors_all[
          module_name
        ],
        col = NA
      )
    )
  }
  
  
  popViewport()
  
  
  # ==========================================================
  # 6B. BOTTOM COLOR BAR
  #
  # Left -> right = ordered_mod2
  # ==========================================================
  
  pushViewport(
    viewport(
      layout.pos.row = 2,
      layout.pos.col = 2
    )
  )
  
  
  for (j in seq_len(n2)) {
    
    module_name <- ordered_mod2[j]
    
    grid.rect(
      
      x = unit(
        (j - 0.5) / n2,
        "npc"
      ),
      
      y = unit(
        0.5,
        "npc"
      ),
      
      width = unit(
        1 / n2,
        "npc"
      ),
      
      height = unit(
        1,
        "npc"
      ),
      
      gp = gpar(
        fill = col_colors_all[
          module_name
        ],
        col = NA
      )
    )
  }
  
  
  popViewport()
  
  
  # ==========================================================
  # 6C. Central blank heatmap-sized area
  # ==========================================================
  
  pushViewport(
    viewport(
      layout.pos.row = 1,
      layout.pos.col = 2
    )
  )
  
  
  grid.rect(
    gp = gpar(
      fill = NA,
      col = NA
    )
  )
  
  
  popViewport()
  
  
  # ==========================================================
  # 6D. Lower-left blank
  # ==========================================================
  
  pushViewport(
    viewport(
      layout.pos.row = 2,
      layout.pos.col = 1
    )
  )
  
  
  grid.rect(
    gp = gpar(
      fill = NA,
      col = NA
    )
  )
  
  
  popViewport()
  
  
  popViewport()
  
  dev.off()
  
  
  # ==========================================================
  # 7. Output checks
  # ==========================================================
  
  cat(
    "\n====================================\n"
  )
  
  cat(
    name1,
    "vs",
    name2,
    "\n"
  )
  
  cat(
    "Rows:",
    n1,
    "\n"
  )
  
  cat(
    "Columns:",
    n2,
    "\n"
  )
  
  cat(
    "\nRow order:\n"
  )
  
  print(
    ordered_mod1
  )
  
  cat(
    "\nColumn order:\n"
  )
  
  print(
    ordered_mod2
  )
  
  cat(
    "\nRow colors:\n"
  )
  
  print(
    row_colors_all[
      ordered_mod1
    ]
  )
  
  cat(
    "\nColumn colors:\n"
  )
  
  print(
    col_colors_all[
      ordered_mod2
    ]
  )
  
  cat(
    "====================================\n"
  )
  
  
  return(
    invisible(
      list(
        plot = p_text,
        cor = cor_mat,
        p = p_mat,
        order1 = ordered_mod1,
        order2 = ordered_mod2,
        row_colors = row_colors_all[
          ordered_mod1
        ],
        col_colors = col_colors_all[
          ordered_mod2
        ]
      )
    )
  )
}


# ============================================================
# 8. BayesCor1 vs BayesCor2
# ============================================================

res_12 <- plot_UCell_Circle_CorrectOrder(
  
  seurat_obj = seurat_obj,
  
  cols1 = bc1_cols,
  cols2 = bc2_cols,
  
  name1 = "BayesCor1",
  name2 = "BayesCor2",
  
  heatmap_pdf =
    "Fig_UCellCircle_Rev_1vs2.pdf",
  
  colorbar_pdf =
    "Fig_UCellCircle_Rev_ColorBars_1vs2.pdf"
)


# ============================================================
# 9. BayesCor1 vs BayesCor3
# ============================================================

res_13 <- plot_UCell_Circle_CorrectOrder(
  
  seurat_obj = seurat_obj,
  
  cols1 = bc1_cols,
  cols2 = bc3_cols,
  
  name1 = "BayesCor1",
  name2 = "BayesCor3",
  
  heatmap_pdf =
    "Fig_UCellCircle_Rev_1vs3.pdf",
  
  colorbar_pdf =
    "Fig_UCellCircle_Rev_ColorBars_1vs3.pdf"
)


# ============================================================
# 10. BayesCor2 vs BayesCor3
# ============================================================

res_23 <- plot_UCell_Circle_CorrectOrder(
  
  seurat_obj = seurat_obj,
  
  cols1 = bc2_cols,
  cols2 = bc3_cols,
  
  name1 = "BayesCor2",
  name2 = "BayesCor3",
  
  heatmap_pdf =
    "Fig_UCellCircle_Rev_2vs3.pdf",
  
  colorbar_pdf =
    "Fig_UCellCircle_Rev_ColorBars_2vs3.pdf"
)
#fig2






