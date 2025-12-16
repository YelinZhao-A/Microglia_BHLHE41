library(Matrix)
library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)
library(broom)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(DHARMa)
library(MASS)   
library(nlme)  
library(lmtest)
library(car)
library(rstatix)
library(ARTool)

rm(list=ls())

output = "output"
if (!dir.exists(output)) dir.create(output)

# Fig2a ####
## Load data####
Fig_2a <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_2a")
head(Fig_2a)
table(Fig_2a$Genotype, Fig_2a$CNS_Region)

Fig_2a$Genotype <- factor(Fig_2a$Genotype, levels = c("WT", "B41HET", "B41KO"))
Fig_2a$CNS_Region <- factor(Fig_2a$CNS_Region)
Fig_2a$Group_Combined <- interaction(Fig_2a$Genotype, Fig_2a$CNS_Region, sep="_")

## Check data distribution####
by(Fig_2a$Microglia, Fig_2a$Group_Combined, shapiro.test)

## Generalized Linear Model (GLM) ######
model_poisson <- glm(Microglia ~ CNS_Region * Genotype, data = Fig_2a, family = poisson(link = "log"))
anova(model_poisson, test="Chisq")

# Check overdispersion
overdispersion <- sum(residuals(model_poisson, type = "pearson")^2) / model_poisson$df.residual

# posthoc test differences between genotype groups
posthoc <- emmeans(model_poisson, pairwise ~ Genotype | CNS_Region, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig2a_Statistics.txt"))
cat("Fig2a Statistics\n")
cat("\nData type: count\n")
cat("\nGeneralized Linear Model (GLM) with poisson:\n")
print(summary(model_poisson))
print(anova(model_poisson, test="Chisq"))
cat("\nOverdispersion parameter:\n")
print(overdispersion)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (EMMeans with Bonferroni adjustment):\n")
print(posthoc$contrasts)
sink()


############################################################
# Fig2c ####
## Load data####
Fig_2c <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_2c")
head(Fig_2c)
table(Fig_2c$Genotype)

Fig_2c$Genotype <- factor(Fig_2c$Genotype, levels = c("WT", "B41cHET", "B41cKO"))

## Check data distribution####
by(Fig_2c$Bhlhe41_Exp, Fig_2c$Genotype, shapiro.test)

## One-way ANOVA####
model <- aov(Bhlhe41_Exp %>% log2 ~ Genotype , data = Fig_2c)
summary(model)

## Check assumptions ####
# Normality of Residuals
Normality = shapiro.test(residuals(model))
pdf(file = file.path(output, "Diagnosis/Fig2c_Residuals_QQPlot.pdf"))
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()
# Homogeneity of Variance
Homogeneity = leveneTest(Bhlhe41_Exp ~ Genotype , data = Fig_2c)

## Switch to Kruskal-Wallis H test ####
# Although the Shapiro-Wilk test was non-significant ($P > 0.05$), likely due to low statistical power ($N=4/5$ per group), 
#the QQ plot revealed a significant departure from normality, particularly in the tails. 
# Therefore, we proceeded with the more robust non-parametric Kruskal-Wallis H test to compare the groups.
KWH = kruskal.test(Bhlhe41_Exp ~ Genotype, data = Fig_2c)

# posthoc test (Dunn's Test) differences between genotype groups
dunn_test_result <- Fig_2c %>% 
  dunn_test(Bhlhe41_Exp ~ Genotype, p.adjust.method = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig2c_Statistics.txt"))
cat("Fig2c Statistics\n")
cat("\nKruskal-Wallis H test:\n")
cat("\nAlthough the Shapiro-Wilk test was non-significant (P > 0.05), likely due to low statistical power (N=4/5 per group), 
    the QQ plot revealed a significant departure from normality, particularly in the tails.Therefore, we proceeded with the more robust non-parametric Kruskal-Wallis H test to compare the groups.\n")
print(KWH)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (Dunn's Test with Bonferroni adjustment):\n")
print(dunn_test_result)
sink()

############################################################
# Fig2d ####
## Load data####
Fig_2d <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_2d")
head(Fig_2d)
table(Fig_2d$Genotype)

Fig_2d$Genotype <- factor(Fig_2d$Genotype, levels = c("WT",  "B41cKO"))

##gMFI ####
## Check data distribution####
by(Fig_2d$gMFI, Fig_2d$Genotype, shapiro.test)

pdf("output/Diagnosis/Fig2d_QQ_Plot_gMFI.pdf", width = 6, height = 6) # 保存为PDF， 6x6英寸
qqnorm(Fig_2d$gMFI, main = "Q-Q Plot of gMFI")  
qqline(Fig_2d$gMFI,  lwd = 2) #  
dev.off()

# Normality of Residuals
Homogeneity = leveneTest(gMFI ~ Genotype, data = Fig_2d)

## Switch to Wilcoxon ####
#the QQ plot revealed a significant departure from normality, particularly in the tails. 
# Therefore, we proceeded with the more robust non-parametric test to compare the groups.
wilcox = wilcox.test(gMFI ~ Genotype, data = Fig_2d)


# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig2d_Statistics.txt")) #append = TRUE
cat("Fig2d gMFI Statistics\n") 
cat("\nthe QQ plot revealed a significant departure from normality, particularly in the tails. Therefore, we proceeded with the more robust non-parametric test to compare the groups.\n")
print(wilcox)
cat("\n=====================\n")
sink()

## Ratio ####
## Check assumptions ####
pdf("output/Diagnosis/Fig2d_QQ_Plot_Ratio.pdf", width = 6, height = 6) # 保存为PDF， 6x6英寸
qqnorm(Fig_2d$Ratio, main = "Q-Q Plot of gMFI") # 可添加标题
qqline(Fig_2d$Ratio,  lwd = 2) # 可自定义参考线颜色和粗细
dev.off()
# Normality of Residuals
Homogeneity = leveneTest(Ratio ~ Genotype, data = Fig_2d)

## Switch to Wilcoxon ####
#the QQ plot revealed a significant departure from normality, particularly in the tails. 
# Therefore, we proceeded with the more robust non-parametric test to compare the groups.
wilcox = wilcox.test(Ratio ~ Genotype, data = Fig_2d)

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig2d_Statistics.txt"),append = TRUE) #append = TRUE
cat("Fig2d Ratio Statistics\n") 
cat("\nthe QQ plot revealed a significant departure from normality, particularly in the tails. Therefore, we proceeded with the more robust non-parametric test to compare the groups.\n")
print(wilcox)
cat("\n=====================\n")
sink()

############################################################
# Fig2f ####
## Load data####
Fig_2f <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_2f")
head(Fig_2f)
table(Fig_2f$Genotype)

Fig_2f$Genotype <- factor(Fig_2f$Genotype, levels = c("WT",  "B41KO"))

## Check data distribution####
by(Fig_2f$AF488_Pos, Fig_2f$Genotype, shapiro.test)

pdf("output/Diagnosis/Fig2f_QQ_Plot.pdf", width = 6, height = 6) # 保存为PDF， 6x6英寸
qqnorm(Fig_2f$AF488_Pos, main = "Q-Q Plot of AF488_Pos") # 可添加标题
qqline(Fig_2f$AF488_Pos, lwd = 2) # 可自定义参考线颜色和粗细
dev.off()

# Normality of Residuals
Homogeneity = leveneTest(AF488_Pos ~ Genotype, data = Fig_2f)

## t test ####
if (Homogeneity$`Pr(>F)`[1] > 0.05) {
  t_test_result <- t.test(AF488_Pos ~ Genotype, data = Fig_2f, var.equal = TRUE)
} else {
  t_test_result <- t.test(AF488_Pos ~ Genotype, data = Fig_2f, var.equal = FALSE)
}

sample_estimates <- t_test_result$estimate
conf_int_lower <- t_test_result$conf.int[1]
conf_int_upper <- t_test_result$conf.int[2]

result_df <- data.frame(
  Variable = deparse(substitute(AF488_Pos)), 
  Test = ifelse(Homogeneity$`Pr(>F)`[1] > 0.05, "Student's t-test", "Welch's t-test"),
  t_statistic = t_test_result$statistic,
  df = t_test_result$parameter,
  p_value = t_test_result$p.value,
  Mean_Group1 = sample_estimates[1],
  Mean_Group2 = sample_estimates[2],
  CI_Lower = conf_int_lower,
  CI_Upper = conf_int_upper,
  stringsAsFactors = FALSE
)%>% t()

group_names <- names(sample_estimates)
rownames(result_df)[6] <- paste0("Mean_", group_names[1])
rownames(result_df)[7] <- paste0("Mean_", group_names[2])

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig2f_Statistics.txt"))
cat("Fig2f AF488_Pos Statistics\n")
cat("\n--------------------\n")
cat("\nCheck assumptions:\n")
cat("\nNormal distribution:\n")
print(by(Fig_2f$AF488_Pos, Fig_2f$Genotype, shapiro.test))
cat("\nHomogeneity of Variance (Levene's Test):\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nt-test Results:\n")
print(result_df)
cat("\n=====================\n")
sink()

############################################################
# Fig2g ####
## Load data####
Fig_2g <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_2g")
head(Fig_2g)
table(Fig_2g$Genotype)
Fig_2g$Genotype <- factor(Fig_2g$Genotype, levels = c("WT",  "B41cKO"))

## Check data distribution####
by(Fig_2g$AF488_Pos, Fig_2g$Genotype, shapiro.test)

## Check assumptions ####
pdf("output/Diagnosis/Fig2g_QQ_Plot.pdf", width = 6, height = 6)  
qqnorm(Fig_2g$AF488_Pos, main = "Q-Q Plot of AF488_MFI")  
qqline(Fig_2g$AF488_Pos, lwd = 2) #  
dev.off()

# Normality of Residuals
Homogeneity = leveneTest(AF488_Pos ~ Genotype, data = Fig_2g)

## t test ####
if (Homogeneity$`Pr(>F)`[1] > 0.05) {
  # Variances are equal, use Student's t-test
  t_test_result <- t.test(AF488_Pos ~ Genotype, data = Fig_2g, var.equal = TRUE)
} else {
  # Variances are not equal, use Welch's t-test
  t_test_result <- t.test(AF488_Pos ~ Genotype, data = Fig_2g, var.equal = FALSE)
}

sample_estimates <- t_test_result$estimate
conf_int_lower <- t_test_result$conf.int[1]
conf_int_upper <- t_test_result$conf.int[2]

result_df <- data.frame(
  Variable = deparse(substitute(AF488_Pos)), 
  Test = ifelse(Homogeneity$`Pr(>F)`[1] > 0.05, "Student's t-test", "Welch's t-test"),
  t_statistic = t_test_result$statistic,
  df = t_test_result$parameter,
  p_value = t_test_result$p.value,
  Mean_Group1 = sample_estimates[1],
  Mean_Group2 = sample_estimates[2],
  CI_Lower = conf_int_lower,
  CI_Upper = conf_int_upper,
  stringsAsFactors = FALSE
) %>% t()

group_names <- names(sample_estimates)
rownames(result_df)[6] <- paste0("Mean_", group_names[1])
rownames(result_df)[7] <- paste0("Mean_", group_names[2])

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig2g_Statistics.txt"))
cat("Fig2g AF488_Pos Statistics\n")
cat("\n--------------------\n")
cat("\nCheck assumptions:\n")
cat("\nNormal distribution:\n")
print(by(Fig_2g$AF488_Pos, Fig_2g$Genotype, shapiro.test))
cat("\nHomogeneity of Variance (Levene's Test):\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nt-test Results:\n")
print(result_df)
cat("\n=====================\n")
sink()

############################################################
# Fig2h ####
## Load data####
Fig_2h <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_2h")
head(Fig_2h)
table(Fig_2h$Genotype)
Fig_2h$Genotype <- factor(Fig_2h$Genotype, levels = c("WT",  "B41KO"))

## Check data distribution####
by(Fig_2h$AF488_Pos, Fig_2h$Genotype, shapiro.test)

pdf("output/Diagnosis/Fig2h_QQ_Plot.pdf", width = 6, height = 6) # 保存为PDF， 6x6英寸
qqnorm(Fig_2h$AF488_Pos , main = "Q-Q Plot of AF488_Pos")  
qqline(Fig_2h$AF488_Pos , lwd = 2)  
dev.off()

## Switch to Wilcoxon ####
#the QQ plot revealed a significant departure from normality. 
# Therefore, we proceeded with the more robust non-parametric test to compare the groups.
wilcox = wilcox.test(AF488_Pos ~ Genotype, data = Fig_2h)

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig2h_Statistics.txt")) #append = TRUE
cat("Fig2h Statistics\n") 
cat("\nthe QQ plot revealed a significant departure from normality, particularly in the tails. Therefore, we proceeded with the more robust non-parametric test to compare the groups.\n")
print(wilcox)
cat("\n=====================\n")
sink()

############################################################
# Fig2i ####
## Load data####
Fig_2i <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_2i")
head(Fig_2i)
table(Fig_2i$Genotype)
Fig_2i$Genotype <- factor(Fig_2i$Genotype, levels = c("WT",  "B41cKO"))

## Check data distribution####
by(Fig_2i$AF488_Pos, Fig_2i$Genotype, shapiro.test)
pdf("output/Diagnosis/Fig2i_QQ_Plot.pdf", width = 6, height = 6)  
qqnorm(Fig_2i$AF488_Pos, main = "Q-Q Plot of AF488_Pos")
qqline(Fig_2i$AF488_Pos, lwd = 2)
dev.off()

## Switch to Wilcoxon ####
#the QQ plot revealed a significant departure from normality.
# Therefore, we proceeded with non-parametric test to compare the groups.
wilcox = wilcox.test(AF488_Pos ~ Genotype, data = Fig_2i)
# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig2i_Statistics.txt")) #append = TRUE
cat("Fig2i Statistics\n")
cat("\nthe QQ plot revealed a significant departure from normality. Therefore, we proceeded with the more robust non-parametric test to compare the groups.\n")
print(wilcox)
cat("\n=====================\n")
sink()
