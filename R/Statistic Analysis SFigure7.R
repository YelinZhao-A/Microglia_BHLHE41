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

###################################
# Fig_S7b####
Fig_S7b <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S7b")
head(Fig_S7b)
table(Fig_S7b$Genotype)
Fig_S7b$Genotype <- factor(Fig_S7b$Genotype, levels = c("WT", "B41HET","B41KO"))

# One-way ANOVA ####
model <- aov(CD22 ~ Genotype, data = Fig_S7b)
anova_summary <- summary(model)
posthoc_aov <- TukeyHSD(model)

# Normality of residuals
shapiro_res <- shapiro.test(residuals(model))
pdf("output/Diagnosis/FigS7b_Residuals_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(model, main = "Q-Q Plot")
dev.off()
# Homogeneity of variances
levene_test <- leveneTest(CD22 ~ Genotype, data = Fig_S7b)

# Kruskal-Wallis test ####
# small sample size
KWH = kruskal.test(CD22 ~ Genotype, data = Fig_S7b)

# posthoc test (Dunn's Test) differences between genotype groups
dunn_test_result <- Fig_S7b %>% 
  dunn_test(CD22 ~ Genotype, p.adjust.method = "bonferroni")


# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS7b_Statistics.txt"))
cat("FigS7b Aurka Statistics\n")
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
# Fig_S7c####
Fig_S7c <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S7c")
head(Fig_S7c)
table(Fig_S7c$Genotype)
Fig_S7c$Genotype <- factor(Fig_S7c$Genotype, levels = c("WT", "B41cHET","B41cKO"))

# One-way ANOVA ####
model <- aov(CD22 ~ Genotype, data = Fig_S7c)
anova_summary <- summary(model)
posthoc_aov <- TukeyHSD(model)

# Normality of residuals
shapiro_res <- shapiro.test(residuals(model))
pdf("output/Diagnosis/FigS7c_Residuals_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(model, main = "Q-Q Plot")
dev.off()
# Homogeneity of variances
levene_test <- leveneTest(CD22 ~ Genotype, data = Fig_S7c)

# Kruskal-Wallis test ####
# small sample size
KWH = kruskal.test(CD22 ~ Genotype, data = Fig_S7c)
# posthoc test (Dunn's Test) differences between genotype groups
dunn_test_result <- Fig_S7c %>% 
  dunn_test(CD22 ~ Genotype, p.adjust.method = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS7c_Statistics.txt"))
cat("FigS7c CD22c Statistics\n")
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







