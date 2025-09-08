library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# === Đường dẫn tới thư mục chứa file ===
data_path <- "/Users/melancholy/Desktop/Báo GraphBG/SpatialGlue/run_GraphBGM-multiModals/Simulated_data_csv"

# === Danh sách các file ===
csv_files <- c(
  "sim1992_metrics_results.csv",
  "sim5024_metrics_results.csv",
  "sim10075_metrics_results.csv",
  "sim1992_1000_metrics_results.csv",
  "sim5024_1000_metrics_results.csv",
  "sim10075_1000_metrics_results.csv"
)

# === Tên tương ứng để hiển thị ===
conditions <- c(
  "n=2000, p=100", "n=5000, p=100", "n=10000, p=100",
  "n=2000, p=1000", "n=5000, p=1000", "n=10000, p=1000"
)

# === Đọc và gộp dữ liệu, chỉ giữ dòng 1 và 2 ===
data_all <- bind_rows(lapply(seq_along(csv_files), function(i) {
  file_path <- file.path(data_path, csv_files[i])
  df <- read_csv(file_path) %>% slice(1:2)   # ⬅️ chỉ lấy dòng 1 và 2
  df$Condition <- conditions[i]
  return(df)
}))

# === Chuyển sang long format ===
data_long <- data_all %>%
  pivot_longer(cols = -c(Method, Condition),
               names_to = "Metric", values_to = "Value")

# === Loại bỏ metric ARI ===
data_long <- data_long %>%
  filter(Metric != "ARI")

# === Giữ thứ tự gốc ===
data_long$Method <- factor(data_long$Method, levels = unique(data_long$Method))
data_long$Condition <- factor(data_long$Condition, levels = conditions)

# === Đổi tên hiển thị của metric ===
metric_labels <- c("NMI" = "NMI", "homogeneity" = "HOM", "completeness" = "COM")
data_long$Metric <- factor(data_long$Metric,
                           levels = names(metric_labels),
                           labels = metric_labels)

# === Màu sắc cho các phương pháp ===
method_colors <- c(
  "GraphBG-MM" = "#F8766D",
  "SpatialGlue" = "#7CAE00"
)

# === Vẽ biểu đồ ===
ggplot(data_long, aes(x = Metric, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  facet_wrap(~Condition, ncol = 3, strip.position = "top", scales = "fixed") +
  scale_fill_manual(values = method_colors) +
  labs(x = "", y = "", fill = "Method") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    strip.placement = "inside",
    strip.background = element_blank(),
    strip.text = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing = unit(1, "lines"),
    legend.text = element_text(size = 12)
  )
