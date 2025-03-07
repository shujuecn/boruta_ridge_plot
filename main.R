library(tidyverse)  # 数据处理与可视化
library(Boruta)     # 特征选择
library(ranger)     # 随机森林
library(ggridges)   # 绘制山脊图
library(viridis)    # 色彩调色板
library(scales)     # 透明度设置

set.seed(250307)

simulate_clinical_data <- function(
    num_patients = 5000,
    num_biomarkers = 30,
    effect_increment = 0.02) {
  # 生成严格均衡的临床结局标签
  outcome_labels <- rep(c(1, 0), each = floor(num_patients / 2))
  outcome_labels <- sample(c(outcome_labels, rep(0, num_patients %% 2)))

  # 初始化生物标志物矩阵
  biomarker_data <- matrix(
    nrow = num_patients,
    ncol = num_biomarkers
  )

  # 生成具有渐进相关性的生物标志物
  for (marker_idx in 1:num_biomarkers) {
    # 计算当前特征的影响系数
    effect_size <- marker_idx * effect_increment

    # 生成具有类别差异的正态分布数据
    biomarker_data[, marker_idx] <- rnorm(
      n = num_patients,
      mean = outcome_labels * effect_size, # 正类样本偏移均值
      sd = 1
    )
  }

  # 创建数据框并设置列名
  clinical_data <- as.data.frame(biomarker_data)
  colnames(clinical_data) <- sprintf("biomarker_%02d", 1:num_biomarkers)
  clinical_data$outcome <- factor(outcome_labels)

  return(clinical_data)
}

# 生成验证数据集（效果系数默认0.02）
clinical_dataset <- simulate_clinical_data(num_patients = 500, num_biomarkers = 30)

# 使用 Boruta 算法进行特征选择
boruta_selection <- Boruta(
  outcome ~ .,
  data = clinical_dataset,
  pValue = 0.01,
  mcAdj = TRUE,
  maxRuns = 50,
  doTrace = 3,
  holdHistory = TRUE,
  getImp = getImpRfZ
)

# 提取特征重要性的历史记录，并转换为数据框
imp_history <- as.data.frame(boruta_selection$ImpHistory)

# 区分真实特征与 shadow 特征
shadow_cols <- grep("shadow", colnames(imp_history), value = TRUE)
real_cols <- setdiff(colnames(imp_history), shadow_cols)

# 构造真实特征对应的决策数据框
final_decision <- tibble(
  feature  = real_cols,
  decision = boruta_selection$finalDecision[real_cols]
)

# 为 shadow 特征统一赋值为 "Shadow"
shadow_decision <- tibble(
  feature  = shadow_cols,
  decision = rep("Shadow", length(shadow_cols))
)

# 合并所有特征的决策信息
decision_all <- bind_rows(final_decision, shadow_decision)

# 将 imp_history 数据转换为长格式，并合并决策信息
imp_history_long <- imp_history %>%
  pivot_longer(
    cols = everything(),
    names_to = "feature",
    values_to = "importance"
  ) %>%
  left_join(decision_all, by = "feature")

# 过滤掉非有限值（例如 NA 或 Inf）
imp_history_long_clean <- imp_history_long %>%
  filter(is.finite(importance))

# 根据每个特征的中位数重要性排序，以确定 y 轴因子顺序
feature_order <- imp_history_long_clean %>%
  group_by(feature) %>%
  summarise(median_imp = median(importance, na.rm = TRUE)) %>%
  arrange(median_imp) %>%
  pull(feature)

# 绘制基于 Boruta 的特征重要性山脊图
ggplot(imp_history_long_clean, aes(
  x = importance,
  y = fct_relevel(feature, feature_order),
  fill = decision
)) +
  geom_density_ridges_gradient(
    scale = 3,
    rel_min_height = 0.01,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Confirmed" = alpha("#02BBC1", 0.6),
      "Tentative" = alpha("#FFC107", 0.6),
      "Rejected"  = alpha("#E53935", 0.6),
      "Shadow"    = alpha("#757575", 0.6)
    )
  ) +
  labs(
    title = "Feature Importance Ridge Plot Based on Boruta",
    x = "Importance (Z-Score)",
    y = "Features"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.spacing = unit(0.1, "lines"),
    axis.text.y = element_text(size = 10)
  )

ggsave("boruta_importance_ridge.jpg")
