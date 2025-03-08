library(tidyverse)
library(Boruta)
library(ranger)
library(ggridges)
library(viridis)
library(scales)
library(glue)

set.seed(250307)

simulate_clinical_data = function(num_patients = 5000, num_biomarkers = 30, effect_increment = 0.02) {

  # 生成严格均衡的临床结局标签
  outcome_labels = rep(c(1, 0), each = floor(num_patients / 2))
  outcome_labels = sample(c(outcome_labels, rep(0, num_patients %% 2))) |>
    as.factor()

  # 生成生物标志物数据并正确赋予列名
  biomarker_data = map_dfc(1:num_biomarkers, function(idx) {
    setNames(
      tibble(rnorm(num_patients, mean = as.numeric(outcome_labels) * idx * effect_increment, sd = 1)),
      paste0("biomarker_", sprintf("%02d", idx))
    )
  })

  clinical_data = tibble(outcome = outcome_labels) |>
    bind_cols(biomarker_data)

  return(clinical_data)
}

# 生成验证数据集（效果系数默认0.02）
clinical_dataset = simulate_clinical_data(num_patients = 500, num_biomarkers = 30)

# 使用 Boruta 算法进行特征选择
boruta_selection = Boruta(
  outcome ~ .,
  data = clinical_dataset,
  pValue = 0.01,
  mcAdj = TRUE,
  maxRuns = 1 + 50,  # 最大迭代次数，根据实际任务调整
  doTrace = 3,
  holdHistory = TRUE,
  getImp = getImpRfZ
)

# 特征
features_decision = enframe(
  boruta_selection$finalDecision,
  name = "feature",
  value = "decision"
)

# 影子
shadow_decision = tibble(
  feature = c("shadowMax", "shadowMean", "shadowMin"),
  decision = "Shadow"
)

# 算法对每个特征的决定
all_decision = bind_rows(features_decision, shadow_decision)

# 提取特征重要性的历史记录
imp_history = as_tibble(boruta_selection$ImpHistory) |>
  mutate(iter = row_number(), .before = 1) |>
  pivot_longer(
    cols = -iter,
    names_to = "feature",
    values_to = "importance"
  ) |>
  filter(is.finite(importance)) |>
  left_join(all_decision, by = "feature")

# 根据每个特征的中位数重要性排序
feature_order = imp_history |>
  summarise(
    median_imp = median(importance, na.rm = TRUE),
    .by = feature
  ) |>
  arrange(median_imp) |>
  pull(feature)

# 创建山脊图
ggplot(imp_history, aes(
  x = importance,
  y = fct_relevel(feature, feature_order),
  fill = decision
)) +
  geom_density_ridges_gradient(
    scale = 2.5,
    rel_min_height = 0.01,
    linewidth = 0,  # 取消山脊线
    # color = "black",
  ) +
  stat_summary(
    fun = median,   # 添加中位数点
    geom = "point",
    shape = 21,
    size = 2,
    alpha = 0.8,
  ) +
  scale_fill_manual(
    values = c(
      "Confirmed" = alpha("#02BBC1", 0.6),
      "Tentative" = alpha("#FFC107", 0.6),
      "Rejected"  = alpha("#E53935", 0.6),
      "Shadow"    = alpha("#757575", 0.6)
    )
  ) +
  scale_x_continuous(breaks = breaks_width(2)) +
  labs(
    title = glue("Boruta Feature Importance (pValue: {boruta_selection$pValue}, Iterations: {nrow(boruta_selection$ImpHistory)})"),
    x = "Importance (Z-Score)",
    y = "Features",
    fill = "Decision"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.spacing = unit(0.1, "lines"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

ggsave("./boruta_importance_ridge.jpg", height = 7, width = 9)
