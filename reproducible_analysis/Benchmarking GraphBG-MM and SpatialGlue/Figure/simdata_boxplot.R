library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# Đường dẫn tới thư mục chứa file
data_path <- "/Users/melancholy/Desktop/Báo GraphBG/SpatialGlue/run_GraphBGM-multiModals/Simulated_data_csv"

# Danh sách các file
csv_files <- c(
  "sim1992_metrics_results.csv",
  "sim5024_metrics_results.csv",
  "sim10075_metrics_results.csv",
  "sim1992_1000_metrics_results.csv",
  "sim5024_1000_metrics_results.csv",
  "sim10075_1000_metrics_results.csv"
)

# Đọc và gộp dữ liệu (🔹 chỉ lấy 2 dòng đầu tiên của mỗi file)
data_all <- bind_rows(lapply(csv_files, function(file) {
  df <- read_csv(file.path(data_path, file))
  df <- head(df, 2)   # ✅ Chỉ lấy 2 dòng đầu tiên
  return(df)
}))

# Chuyển sang long format
data_long <- data_all %>%
  pivot_longer(cols = -Method, names_to = "Metric", values_to = "Value")

# Đổi thứ tự Metric theo yêu cầu
metric_order <- c("NMI", "homogeneity", "completeness")
metric_labels <- c(
  "NMI" = "NMI",
  "homogeneity" = "HOM",
  "completeness" = "COM"
)
data_long <- data_long %>%
  filter(Metric %in% metric_order) %>%
  mutate(Metric = factor(Metric, levels = metric_order, labels = metric_labels))

# Màu sắc cho các phương pháp
method_colors <- c(
  "GraphBG-MM" = "#F8766D",
  "SpatialGlue" = "#7CAE00"
)

# Vẽ biểu đồ boxplot
ggplot(data_long, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(
    width = 0.6,
    color = "black",
    outlier.shape = 21,
    outlier.size = 2,
    alpha = 0.95,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ Metric, nrow = 1) +
  labs(
    x = "",
    y = NULL,
    fill = "Method",
    title = ""
  ) +
  scale_fill_manual(values = method_colors) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "#d3d3d3", size = 0.4),
    panel.grid.minor = element_line(color = "#d3d3d3", size = 0.2),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 17),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    panel.spacing = unit(0.2, "lines"),
    plot.margin = margin(5, 5, 5, 5),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 12)
  )

# Lưu file ảnh
ggsave("/Users/melancholy/Desktop/simdata_boxplot.png", width = 10, height = 3.5, dpi = 300)
