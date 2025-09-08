library(ggplot2)
library(readr)
library(dplyr)

# Đường dẫn tới các file CSV
files <- c(
  "/Users/melancholy/Desktop/Báo GraphBG/SpatialGlue/run_GraphBGM-multiModals/Thymus_csv/thymus_moran_results1.csv",
  "/Users/melancholy/Desktop/Báo GraphBG/SpatialGlue/run_GraphBGM-multiModals/Thymus_csv/thymus_moran_results2.csv",
  "/Users/melancholy/Desktop/Báo GraphBG/SpatialGlue/run_GraphBGM-multiModals/Thymus_csv/thymus_moran_results3.csv",
  "/Users/melancholy/Desktop/Báo GraphBG/SpatialGlue/run_GraphBGM-multiModals/Thymus_csv/thymus_moran_results4.csv",
  "/Users/melancholy/Desktop/Báo GraphBG/SpatialGlue/run_GraphBGM-multiModals/Spleen_csv/spleen_moran_results1.csv",
  "/Users/melancholy/Desktop/Báo GraphBG/SpatialGlue/run_GraphBGM-multiModals/Spleen_csv/spleen_moran_results2.csv",
  "/Users/melancholy/Desktop/Báo GraphBG/SpatialGlue/run_GraphBGM-multiModals/Brain_csv/brain_moran_results_ATAC.csv",
  "/Users/melancholy/Desktop/Báo GraphBG/SpatialGlue/run_GraphBGM-multiModals/Brain_csv/brain_moran_results_H3K4me3.csv"
)

# Điều kiện
conditions <- c("mouse_thymus_1", "mouse_thymus_2", "mouse_thymus_3", "mouse_thymus_4", "mouse_spleen_1", "mouse_spleen_2", "mouse_brain_ATAC","mouse_brain_H3K4me3")

# Đọc và gộp dữ liệu
data_all <- bind_rows(lapply(1:8, function(i) {
  df <- read_csv(files[i])
  df$Condition <- conditions[i]
  return(df)
}))

# Giữ thứ tự cho các yếu tố
data_all$Method <- factor(data_all$Method, levels = unique(data_all$Method))
data_all$Condition <- factor(data_all$Condition, levels = conditions)

# Màu sắc cho các phương pháp
method_colors <- c(
  "GraphBG-MM" = "#F8766D",
  "SpatialGlue" = "#7CAE00"
)

# Vẽ biểu đồ chia thành 2 hàng, mỗi hàng 4 biểu đồ
ggplot(data_all, aes(x = Method, y = Moran_I, fill = Method)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~Condition, nrow = 2, ncol = 4, strip.position = "bottom", scales = "fixed") +
  scale_fill_manual(values = method_colors) +
  labs(x = NULL, y = "Moran's I", fill = "Method") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.title.y = element_text(, size = 14),
    axis.text.y = element_text(size = 12),
    #legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing = unit(0.5, "lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
ggsave("/Users/melancholy/Desktop/thymus_spleen_brain_barplot2.png", width = 10, height = 5, dpi = 300)

