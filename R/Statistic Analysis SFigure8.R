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

########################################

# Fig_S8a####
Fig_S8a <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S8ab")
head(Fig_S8a)
table(Fig_S8a$Genotype)
Fig_S8a$Genotype <- factor(Fig_S8a$Genotype, levels = c("WT", "B41HET","B41KO"))

# One-way ANOVA ####
model <- aov(AM_CD22 ~ Genotype, data = Fig_S8a)
anova_summary <- summary(model)
posthoc_aov <- TukeyHSD(model)
# Normality of residuals
shapiro_res <- shapiro.test(residuals(model))
pdf("output/Diagnosis/FigS8a_Residuals_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(model, main = "Q-Q Plot")
dev.off()
# Homogeneity of variances
levene_test <- leveneTest(AM_CD22 ~ Genotype, data = Fig_S8a)

# Kruskal-Wallis test ####
# small sample size
KWH = kruskal.test(AM_CD22 ~ Genotype, data = Fig_S8a)
# posthoc test (Dunn's Test) differences between genotype groups
dunn_test_result <- Fig_S8a %>% 
  dunn_test(AM_CD22 ~ Genotype, p.adjust.method = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS8a_Statistics.txt"))
cat("FigS8a AM_CD22 Statistics\n")
cat("\nKruskal-Wallis H test:\n")
cat("\n small sample size, use non-parametric Kruskal-Wallis H test .\n")
print(KWH)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (Dunn's Test with Bonferroni adjustment):\n")
print(dunn_test_result)
cat("\n===============================================\n")
cat("\n===============================================\n")
cat("\n===Alternative: One-way Anova==========\n")
print(model)
print(anova_summary)
cat("\nNormality of Residuals (Shapiro-Wilk test):\n")
print(shapiro_res)
cat("\nHomogeneity of Variances (Levene's Test):\n")
print(levene_test)
cat("\nOne-way ANOVA Post-hoc Test Results (Tukey HSD):\n")
print(posthoc_aov)
sink()

########################################
# Fig_S8b####
Fig_S8b <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S8ab")
head(Fig_S8b)
table(Fig_S8b$Genotype)
Fig_S8b$Genotype <- factor(Fig_S8b$Genotype, levels = c("WT", "B41HET","B41KO"))

# One-way ANOVA ####
model <- aov(B1_CD22 ~ Genotype, data = Fig_S8b)
anova_summary <- summary(model)
posthoc_aov <- TukeyHSD(model)

# Normality of residuals
shapiro_res <- shapiro.test(residuals(model))
pdf("output/Diagnosis/FigS8b_Residuals_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(model, main = "Q-Q Plot")
dev.off()
# Homogeneity of variances
levene_test <- leveneTest(B1_CD22 ~ Genotype, data = Fig_S8b)
# Kruskal-Wallis test ####
# small sample size
KWH = kruskal.test(B1_CD22 ~ Genotype, data = Fig_S8b)
# posthoc test (Dunn's Test) differences between genotype groups
dunn_test_result <- Fig_S8b %>% 
  dunn_test(B1_CD22 ~ Genotype, p.adjust.method = "bonferroni")
# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS8b_Statistics.txt"))
cat("FigS8b B1_CD22 Statistics\n")
cat("\nKruskal-Wallis H test:\n")
cat("\n small sample size, use non-parametric Kruskal-Wallis H test .\n")
print(KWH)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (Dunn's Test with Bonferroni adjustment):\n")
print(dunn_test_result)
cat("\n===============================================\n")
cat("\n===============================================\n")
cat("\n===Alternative: One-way Anova==========\n")
print(model)
print(anova_summary)
cat("\nNormality of Residuals (Shapiro-Wilk test):\n")
print(shapiro_res)
cat("\nHomogeneity of Variances (Levene's Test):\n")
print(levene_test)
cat("\nOne-way ANOVA Post-hoc Test Results (Tukey HSD):\n")
print(posthoc_aov)
sink()









