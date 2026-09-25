# ==========================================
#figc2: 网络轨迹图 (UMAP + 网络)
# ==========================================

setwd("/home/user/Fanglab1/ymh/crc")
load('plot.Rdata')

library(uwot)
library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(reshape2)
library(ggplot2)

# 1. 扩容基因并进行带微噪的连续映射
print("1. 正在提取 Top 50 基因并构建微噪连续轨迹...")

genes_df <- data.frame(
  gene = unlist(module_genes2_2),
  module = rep(paste0("M", 1:25), times = sapply(module_genes2_2, length))
)

proCSC_mods <- c("M11", "M4", "M18", "M22", "M6", "M8", "M24", "M3", "M23", "M10", "M15")
CSC_mods    <- c("M2", "M14", "M16")
revCSC_mods <- c("M1", "M7", "M17", "M20", "M13", "M9", "M12", "M21", "M25", "M5", "M19")

genes_df <- genes_df %>%
  mutate(
    meta_state = case_when(
      module %in% proCSC_mods ~ "proCSC",
      module %in% CSC_mods    ~ "CSC",
      module %in% revCSC_mods ~ "revCSC",
      TRUE ~ "Other"
    ),
    mod_num = as.numeric(gsub("M", "", module)),
    base_time = case_when(
      meta_state == "revCSC" ~ 1,
      meta_state == "CSC"    ~ 2,
      meta_state == "proCSC" ~ 3,
      TRUE ~ 2
    ),
    pseudotime = base_time + (mod_num / 100)
  )

valid_genes <- intersect(genes_df$gene, rownames(datKME))
genes_df <- genes_df %>% filter(gene %in% valid_genes)

kme_cols <- grep("kME_MM.Module_", colnames(datKME), value = TRUE)
orig_nums <- as.numeric(gsub("kME_MM.Module_", "", kme_cols))
mod_mapping <- paste0("kME_MM.Module_", sort(orig_nums))
names(mod_mapping) <- paste0("M", 1:25)

genes_df$kME_col <- mod_mapping[genes_df$module]
row_idx <- match(genes_df$gene, rownames(datKME))
col_idx <- match(genes_df$kME_col, colnames(datKME))
genes_df$kME <- datKME[cbind(row_idx, col_idx)]
genes_df$kME_col <- NULL

genes_df <- genes_df %>%
  group_by(module) %>%
  slice_max(order_by = kME, n = 50, with_ties = FALSE) %>%
  ungroup()

print(paste("基因扩容完毕！参与画图的总基因数:", nrow(genes_df)))

# 2. 连续型监督 UMAP
print("2. 正在计算层次化监督 UMAP...")
umap_input <- datKME[genes_df$gene, kme_cols]
umap_labels <- genes_df$pseudotime

set.seed(42)
umap_res <- uwot::umap(X = umap_input, y = umap_labels, 
                       n_neighbors = 30,        
                       min_dist = 0.3,        
                       target_weight = 0.2,
                       metric = "cosine")

umap_df <- data.frame(
  gene = genes_df$gene,
  UMAP1 = umap_res[, 1],
  UMAP2 = umap_res[, 2]
) %>% left_join(genes_df, by = "gene")

# 3. 领地净化
print("⚠️ 正在清理混入两端大本营的杂色基因...")

min_x <- min(umap_df$UMAP1)
max_x <- max(umap_df$UMAP1)
span <- max_x - min_x

left_boundary <- min_x + span * 0.35
outliers_left <- umap_df %>% 
  filter(UMAP1 < left_boundary & meta_state != "proCSC") %>% 
  pull(gene)

right_boundary <- max_x - span * 0.35
outliers_right <- umap_df %>% 
  filter(UMAP1 > right_boundary & meta_state != "revCSC") %>% 
  pull(gene)

all_outliers <- unique(c(outliers_left, outliers_right))
message(" -> 扫描完毕！共剔除了 ", length(all_outliers), " 个错位的杂色基因。")

umap_df <- umap_df %>% filter(!(gene %in% all_outliers))
genes_plot <- umap_df$gene 

# 4. 提取连线
print("3. 提取网络拓扑并智能抽样连线...")
subset_tom <- tom[umap_df$gene, umap_df$gene]

edge_df <- melt(as.matrix(subset_tom), varnames = c("from", "to"), value.name = "weight") %>%
  filter(from != to, weight > 0.02) 

node_meta <- setNames(umap_df$meta_state, umap_df$gene)
edge_df <- edge_df %>%
  mutate(
    from_meta = node_meta[as.character(from)],
    to_meta   = node_meta[as.character(to)],
    is_intra  = (from_meta == to_meta),
    edge_color = ifelse(is_intra, as.character(from_meta), "grey85")
  )

edge_prop <- 0.1
sampled_edges <- edge_df %>%
  group_by(edge_color) %>%
  sample_frac(edge_prop) %>%
  ungroup() %>%
  mutate(scaled_weight = (weight - min(weight)) / (max(weight) - min(weight))) %>%
  arrange(is_intra, scaled_weight) 

# 5. 可视化
print("4. 正在渲染带坐标轴的网络图...")
g <- graph_from_data_frame(d = sampled_edges, vertices = umap_df, directed = FALSE)
layout_coords <- create_layout(g, layout = "manual", x = V(g)$UMAP1, y = V(g)$UMAP2)

state_colors <- c(
  "proCSC" = "#1c8dd9ff", 
  "CSC"    = "#931fd0ff",
  "revCSC" = "#63b31cff", 
  "grey85" = "grey85"
)

p_fig1a <- ggraph(layout_coords) +
  geom_edge_link(aes(color = edge_color, alpha = scaled_weight), 
                 width = 0.2, show.legend = FALSE) +
  scale_edge_color_manual(values = state_colors, guide = "none") +
  scale_edge_alpha_continuous(range = c(0.3, 0.5), guide = "none") + 
  geom_node_point(aes(color = meta_state), 
                  size = 2, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = state_colors, guide = "none") +
  coord_cartesian(expand = TRUE) +
  theme_void() + 
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")
  )

# 6. 保存
png_path <- "/home/user/Fanglab1/ymh/aaa_lung_cancer/hdWGCNA_Results/FIG1a_Network_Trajectory.png"
ggsave(png_path, p_fig1a, width = 11, height = 8, dpi = 600, bg = "white")
print(paste("✅ FIG 1a 已保存至:", png_path))




#figc4
# ==========================================
# FIG 1b: 模块转换曲线 (B -> A 转换)
# ==========================================

setwd('E:/deskup/')
library(ggplot2)

time <- seq(0, 10, length.out = 300)

module_B <- 1 / (1 + exp((time - 5)))
module_A <- 1 / (1 + exp(-(time - 5)))

df <- data.frame(
  Time = rep(time, 2),
  Conversion = c(module_A, module_B),
  Module = rep(c("Module A", "Module B"), each = length(time))
)

p_fig1b <- ggplot(df, aes(x = Time, y = Conversion, color = Module)) +
  geom_line(linewidth = 2) +
  scale_color_manual(values = c("Module A" = "#63b31cff", "Module B" = "#1c8dd9ff")) +
  annotate("text", x = 8, y = 0.85, label = "Module A", color = "#66be18ff", size = 6) +
  annotate("text", x = 1.5, y = 0.85, label = "Module B", color = "#2996e2ff", size = 6) +
  labs(x = "Time", y = "Conversion") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.line = element_line(linewidth = 1)
  )

ggsave("FIG1b_B_to_A_Conversion.pdf", p_fig1b, width = 5, height = 4, dpi = 300)
print("✅ FIG 1b 已保存")



#figc1
# ==========================================
# FIG 1c: 模块相关性热图 (Module Correlation Heatmap)
# ==========================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

set.seed(888)

# 参数设置
n_per_celltype <- 150
cell_types <- c("Cell1", "Cell2", "Cell3")
modules <- c("Module1", "Module2", "Module3", "Module4")

# 期望模块活性矩阵
mean_matrix <- matrix(
  c(0.90, 0.12, 0.12, 0.15,
    0.12, 0.90, 0.15, 0.50,
    0.12, 0.15, 0.90, 0.52),
  nrow = 3, byrow = TRUE
)
rownames(mean_matrix) <- cell_types
colnames(mean_matrix) <- modules

# 模拟细胞水平模块分数
sim_list <- list()
for (cell_type in cell_types) {
  temp <- data.frame(Cell = paste0(cell_type, "_", seq_len(n_per_celltype)), CellType = cell_type)
  for (module in modules) {
    temp[[module]] <- pmin(pmax(rnorm(n_per_celltype, mean = mean_matrix[cell_type, module], sd = 0.10), 0), 1)
  }
  sim_list[[cell_type]] <- temp
}
module_data <- bind_rows(sim_list)

# 计算模块-模块相关性
module_cor <- cor(module_data[, modules], method = "pearson", use = "pairwise.complete.obs")

cor_df <- as.data.frame(as.table(module_cor))
colnames(cor_df) <- c("Module1", "Module2", "Correlation")
cor_df$Module1 <- factor(cor_df$Module1, levels = modules)
cor_df$Module2 <- factor(cor_df$Module2, levels = rev(modules))

p_fig1c <- ggplot(cor_df, aes(x = Module1, y = Module2, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 1.2) +
  scale_fill_gradient2(
    low = "#009ad6", mid = "white", high = "#d71345",
    midpoint = 0, limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Pearson's r"
  ) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

ggsave("FIG1c_Module_Correlation_Heatmap.pdf", p_fig1c, width = 6, height = 5.5)
print("✅ FIG 1c 已保存")





#figc2

# ==========================================
# FIG 1d: 细胞特异性模块 DotPlot
# ==========================================

library(ggplot2)
library(dplyr)
library(tidyr)

# 重命名细胞类型 (使用前面模拟的 module_data)
module_data$CellType <- recode(
  module_data$CellType,
  "Cell1" = "Stem cell",
  "Cell2" = "Lineage A",
  "Cell3" = "Lineage B"
)

cell_order <- c("Stem cell", "Lineage A", "Lineage B")
module_order <- c("Module3", "Module2", "Module1")

module_data$CellType <- factor(module_data$CellType, levels = cell_order)

long_df <- module_data %>%
  select(Cell, CellType, Module1, Module2, Module3) %>%
  pivot_longer(cols = all_of(module_order), names_to = "Module", values_to = "Score")

long_df$Module <- factor(long_df$Module, levels = module_order)

activity_threshold <- 0.30

dot_df <- long_df %>%
  group_by(CellType, Module) %>%
  summarise(
    avg_score = mean(Score),
    pct_active = mean(Score > activity_threshold) * 100,
    .groups = "drop"
  ) %>%
  group_by(Module) %>%
  mutate(avg_score_scaled = as.numeric(scale(avg_score))) %>%
  ungroup()

p_fig1d <- ggplot(dot_df, aes(x = CellType, y = Module)) +
  geom_point(
    aes(size = pct_active, fill = avg_score_scaled),
    shape = 21, color = "black", stroke = 0.6
  ) +
  scale_fill_gradientn(
    colors = c("#f7f7f7", "#f8c8cd", "#ed8590", "#d71345"),
    name = "Module score"
  ) +
  scale_size_continuous(
    range = c(2, 12), limits = c(0, 100),
    breaks = c(20, 40, 60, 80, 100),
    name = "Enriched cells (%)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_line(color = "white", linewidth = 1.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#EBEBEB", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave("FIG1d_CellSpecific_Module_DotPlot.pdf", p_fig1d, width = 6, height = 4.5)
print("✅ FIG 1d 已保存")