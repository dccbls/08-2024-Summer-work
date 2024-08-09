#载入依赖库
library(readxl)
library(ggplot2)
library(ggsignif)

# 读取 Excel 文件
df <- read_excel("/Users/yding/NIP_9311_TE_UTRlength.xlsx", sheet = 3)
df_selected <- df[, 1:6] #读取前六列数据

# 重命名列
colnames(df_selected) <- c("NIP_trans", "NIP_TE", "NIP_UTR", "9311_trans", "9311_TE", "9311_UTR")

# 创建长格式数据框
df_utr_long <- data.frame(
  Group = rep(c("NIP", "9311"), each = nrow(df_selected)),
  UTR = c(df_selected$NIP_UTR, df_selected$`9311_UTR`)
)

df_te_long <- data.frame(
  Group2 = rep(c("NIP","9311"), each = nrow(df_selected)),
  TE = c(df_selected$NIP_TE,df_selected$`9311_TE`)
)
"""
#计算离群值
iqr <- IQR(df_utr_long$UTR,na.rm = TRUE)
lower_bound <- quantile(df_utr_long$UTR, 0.25, na.rm = TRUE) - 1.5 * iqr
upper_bound <- quantile(df_utr_long$UTR, 0.75, na.rm = TRUE) + 1.5 * iqr

df_utr_clean <- df_utr_long[df_utr_long$UTR >= lower_bound & df_utr_long$UTR <= upper_bound, ]
"""

# 绘制UTR箱形图
y_pos <- max(df_utr_long$UTR,na.rm = TRUE) * 0.15 #调整显著值高度

ggplot(df_utr_long, aes(x = Group, y = UTR, fill = Group)) +
  geom_boxplot(outlier.shape = NA,width = 0.2) +  # 隐藏离群值，设置箱子宽度
  geom_signif(comparisons = list(c("NIP", "9311")),
              map_signif_level = TRUE,
              textsize = 3,
              y_position = y_pos,
              tip_length = 0) +
  labs(title = "3'UTR Lengths",
       x = "Group",
       y = "Length") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),  # 标题居中
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(linewidth = 1,color = 'black'),
    panel.border = element_rect(color = 'black',linewidth = 1.5)
    ) +
  scale_fill_manual(values = c("skyblue", "salmon")) +
  #coord_cartesian(ylim = c(0,1000))
  coord_cartesian(ylim = c(min(df_utr_long$UTR, na.rm = TRUE), y_pos * 1.5))  # 确保显著性标记在显示范围内

#绘制TE箱形图

ggplot(df_te_long, aes(x = Group2, y = TE, fill = Group2)) +
  geom_boxplot(outlier.shape = NA,width = 0.2) +  # 隐藏离群值，设置箱子宽度
  geom_signif(comparisons = list(c("NIP", "9311")),
              map_signif_level = TRUE,
              textsize = 3,
              y_position = 1, #显著性位置
              tip_length = 0) + #取消显著性竖线
  labs(title = "TE",
       x = "Group",
       y = "values") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),  # 标题居中
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.ticks.length = unit(0.25,"cm"),
    axis.ticks = element_line(linewidth = 1,color = 'black'),
    panel.border = element_rect(color = "black",linewidth = 1.5),#加粗边框
    panel.grid.major = element_blank(),#取消网格线
    panel.grid.minor = element_blank()#取消网格线
  ) +
  scale_fill_manual(values = c("lightgreen", "orange")) +
  coord_cartesian(ylim = c(-1, 3.5 ))  # 设置 y 轴的范围，以确保显著性标记显示在图中