# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(DESeq2)
library(stringr)
library("pheatmap")
library("RColorBrewer")
library(fgsea)
library(ComplexHeatmap)
library(circlize)
library(rstatix)
library(ggpubr)
library(ggrepel)
library(preprocessCore)
library(clusterProfiler)
library(enrichplot)
library(ggbeeswarm)
library(MASS)      # 负二项回归
library(car)       # 用于Anova() 获取交互作用的P值
library(emmeans)   # 用于事后两两比较
library(broom)     # 整理结果

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#raw data preprocessing
Sample_features <- read.table("../../../TGCT_subtyping/Info_temp.txt",sep = "\t",header = T)
samples <- data.frame(Sample_features)
rownames(samples) <- gsub("-",".",samples$AID)
samples$AID <- rownames(samples) 


##########################################
myTable <- samples
df <- myTable[which(!is.na(myTable$RT) & !is.na(myTable$chr1p36) ),]

df$Histology <- "nonYST"
df$Histology[which(df$YST_1st_his == "YST")] <- "YST"
df$CNV <- "no_loss"
df$CNV[which(df$chr1p36 == "Loss")] <- "loss"

df$Histology <- factor(df$Histology, levels = c("nonYST", "YST"))
df$CNV <- factor(df$CNV, levels = c("no_loss", "loss"))
table(df$Histology )
table(df$CNV)
df$L1 <- df$RT
# 1. 拟合负二项回归（含交互项）
model <- glm.nb(L1 ~ Histology * CNV, data = df)

# 2. 查看模型摘要
summary(model)

# 3. 使用 Anova() 获取各因素及交互项的卡方检验P值（Type II/III Wald检验）
# 推荐使用 Type II（如果是平衡设计）或 Type III（考虑缺失组合）
Anova(model, type = "III")
Anova(model, type = "II")

# 计算 Incidence Rate Ratios (IRR) 并画森林图
library(ggforestplotR)  # 或者 library(forestplot)

# 提取所有系数的 IRR 和 CI（排除截距 Intercept）
res <- data.frame(
  Variable = names(coef(model)),          # 自动提取所有系数名
  IRR = exp(coef(model)),
  CI_low = exp(confint(model)[, 1]),
  CI_high = exp(confint(model)[, 2])
)

# 删除截距行（因为我们只关心变量）
res <- res[!grepl("Intercept", res$Variable), ]

print(res)





# 计算均值和标准误
plot_data <- df %>%
  group_by(Histology, CNV) %>%
  summarise(Mean = mean(L1, na.rm = TRUE), 
            SE = sd(L1, na.rm = TRUE) / sqrt(n())) %>%
  ungroup()

x <- ggplot(plot_data, aes(x = Histology, y = Mean, color = CNV, group = CNV)) +
  geom_line(size = 1.2) +                # 连接线，看交互走向
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.1) +
  geom_point(size = 4) +
  labs(y = "Mean Somatic L1 Insertions", color = "1p36 Status") +
  theme_classic() +
  theme(legend.position = "top")
figure_2<-cbind(ggplotGrob(x),size="last")
ggsave(file="./forest.pdf", plot=figure_2,bg = 'white', width =10, height = 9, units = 'cm', dpi = 600)
####################################



