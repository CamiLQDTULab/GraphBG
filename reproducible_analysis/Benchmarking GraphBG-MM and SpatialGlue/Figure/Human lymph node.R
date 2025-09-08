library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# Đường dẫn tới file CSV duy nhất
csv_file <- "/Users/melancholy/Desktop/Báo GraphBG/SpatialGlue/run_GraphBGM-multiModals/human_metrics_results.csv"

# Đặt tên điều kiện tương ứng (dùng cho tiêu đề)
condition_name <- ""

# Đọc dữ liệu
data_all <- read_csv(csv_file)
data_all$Condition <- condition_name

# Chuyển sang long format
data_long <- data_all %>%
  pivot_longer(cols = -c(Method, Condition),
               names_to = "Metric", values_to = "Value")

# Giữ thứ tự ban đầu
data_long$Method <- factor(data_long$Method, levels = unique(data_long$Method))

# Đổi tên hiển thị của metric
metric_order <- c("NMI", "homogeneity", "completeness")
metric_labels <- c(
  "NMI" = "NMI",
  "homogeneity" = "HOM",
  "completeness" = "COM"
)

# Gán lại nhãn
data_long <- data_long %>%
  filter(Metric %in% metric_order) %>%
  mutate(Metric = factor(Metric, levels = metric_order,
                         labels = metric_labels))

# Màu sắc cho phương pháp
method_colors <- c(
  "GraphBG-MM" = "#F8766D",
  "SpatialGlue" = "#7CAE00"
)


ggplot(data_long, aes(x = Method, y = Value, fill = Method)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.7)) +
  facet_wrap(~ Metric, nrow = 1) +
  scale_fill_manual(values = method_colors) +
  labs(title = condition_name, x = "", y = "", fill = "Method") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_blank(),          # Ẩn chữ dưới trục x
    axis.ticks.x = element_blank(),         # Ẩn vạch chia trục x
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(face = "bold", size = 14),
    #legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.text = element_text(size = 12)
  )
ggsave("/Users/melancholy/Desktop/human_lymph_node_A1_barplot.png", width = 10, height = 3.5, dpi = 300)

