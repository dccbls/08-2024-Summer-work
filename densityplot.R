# 加载必要的库
library(readxl)
library(ggplot2)
library(ggsignif)
library(ggdensity)

# 读取数据
df <- read_excel("/Users/yding/NIP_9311_TE_UTRlength.xlsx", sheet = 3)
df_selected <- df[, c(1:6,14:18)]  # 读取前六列数据

# 重命名列
colnames(df_selected) <- c("NIP_trans", "NIP_aveTE", "NIP_3UTR", "9311_trans", "9311_aveTE",'9311_3UTR','NIP_Trans','NIP_logTE','NIP_logUTR','9311_logTE','9311_logUTR')

# 计算相关系数和 p 值
cor_test <- cor.test(df_selected$`NIP_logUTR`, df_selected$`NIP_logTE`)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

# 格式化 p 值
p_value_formatted <- format_p_value(p_value)

# 创建基础密度图
plot <- ggplot(df_selected, aes(x = `NIP_logUTR`, y = `NIP_logTE`)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour = TRUE) +
  scale_fill_gradient(low = '#fdece3', high = '#893508') +
  labs(x = "NIP UTR(Log)", y = "Translation Efficiency(Log)")

# 添加边框和去掉网格线
plot + annotate("text", x = Inf, y = Inf, 
                label = paste("r =", round(r_value, 2), "\np =", p_value_formatted), 
                hjust = 1.1, vjust = 1.5, 
                size = 5, color = "black", 
                fontface = "bold", 
                lineheight = 0.8) +
  theme_bw() +  # 使用经典主题
  theme(
    panel.border = element_rect(color = "black", linewidth = 1.5),  # 添加黑色边框
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.background = element_blank(),  # 去掉面板背景
    plot.background = element_blank(), # 去掉图背景
    axis.text.x = element_text(size = 14, face = "bold"), 
    # 调整 y 轴刻度字体
    axis.text.y = element_text(size = 14, face = "bold"),
    # 如果需要调整轴标题的字体
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold")
  )


