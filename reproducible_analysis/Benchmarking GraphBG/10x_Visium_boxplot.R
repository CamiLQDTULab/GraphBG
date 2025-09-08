library(ggplot2)
library(readr)
library(dplyr)
library(forcats)

# Đọc dữ liệu
df <- read_csv("/Users/melancholy/Desktop/metrics_combined.csv")

# Chuyển GraphBG lên đầu thủ công nếu cần
df$Method <- fct_relevel(df$Method, "GraphBG", after = 0)

com <- ggplot(df, aes(x = Method, y = COM, fill = Method)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 23) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 23),  # Tăng cỡ chữ tên phương pháp
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(
    title = "",
    x = "",
    y = "COM"
  )
ggsave(filename = "/Users/melancholy/Desktop/COM.png", plot = com, width = 10, height = 7, dpi = 300)

hom <- ggplot(df, aes(x = Method, y = HOM, fill = Method)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 23) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 23),  # Tăng cỡ chữ tên phương pháp
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(
    title = "",
    x = "",
    y = "HOM"
  )
ggsave(filename = "/Users/melancholy/Desktop/HOM.png", plot = hom, width = 10, height = 7, dpi = 300)


nmi <- ggplot(df, aes(x = Method, y = NMI, fill = Method)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 23) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 23),  # Tăng cỡ chữ tên phương pháp
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(
    title = "",
    x = "",
    y = "NMI"
  )
ggsave(filename = "/Users/melancholy/Desktop/NMI.png", plot = nmi, width = 10, height = 7, dpi = 300)


install.packages("magick")  # nếu bạn chưa cài
library(magick)
img1 <- image_read("/Users/melancholy/Desktop/NMI.png")
img2 <- image_read("/Users/melancholy/Desktop/HOM.png")
img3 <- image_read("/Users/melancholy/Desktop/COM.png")
combined_img <- image_append(c(img1, img2, img3))
print(combined_img)

# Lưu thành file PNG
image_write(combined_img, path = "/Users/melancholy/Desktop/10x_Visium_NMI_HOM_COM.png", format = "png")

