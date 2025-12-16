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
# Fig_S1d####
Fig_S1d <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S1de")
head(Fig_S1d)

# B1_dTomato
table(Fig_S1d$Genotype)
Fig_S1d$Genotype <- factor(Fig_S1d$Genotype, levels = c("WT", "B41HET","B41KO"))

#One-way ANOVA ####
model <- aov(B1_dTomato ~ Genotype, data = Fig_S1d)
anova_result <- summary(model)
# Check ANOVA assumptions
# Normality of residuals
shapiro_result <- shapiro.test(residuals(model))
pdf("output/Diagnosis/FigS1d_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()
# Homogeneity of variances
Homogeneity <- car::leveneTest(B1_dTomato ~ Genotype, data = Fig_S1d)

# Kruskal–Wallis test ####
kruskal_result <- kruskal.test(B1_dTomato ~ Genotype, data = Fig_S1d)
# Post-hoc - Dunn’s post-hoc test
posthoc_result <- Fig_S1d %>% 
  dunn_test(B1_dTomato ~ Genotype, p.adjust.method = "bonferroni")

# Report results
sink(file = file.path(output, "FigS1d_Statistics.txt"))
cat("FigS1d B1_dTomato Statistics\n")
cat("\nAssumption Checks:\n")
cat("\nNormality of Residuals:\n")
print(shapiro_result)
cat("\nHomogeneity of Variances:\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nKruskal-Wallis Test Results:\n")
print(kruskal_result)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (Dunn’s post-hoc test with Bonferroni correction):\n")
print(posthoc_result)
cat("\n=====================\n")
sink()

###################################
# Fig_S1e####
Fig_S1e <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S1de")
head(Fig_S1e)

# AM_dTomato
table(Fig_S1e$Genotype)
Fig_S1e$Genotype <- factor(Fig_S1e$Genotype, levels = c("WT", "B41HET","B41KO"))
Fig_S1e$AM_dTomato_Log2 <- log2(Fig_S1e$AM_dTomato) 

#One-way ANOVA ####
model <- aov(AM_dTomato_Log2 ~ Genotype, data = Fig_S1e)
anova_result <- summary(model)
# Check ANOVA assumptions
# Normality of residuals
shapiro_result <- shapiro.test(residuals(model))
pdf("output/Diagnosis/FigS1e_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()
# Homogeneity of variances
Homogeneity <- car::leveneTest(AM_dTomato_Log2 ~ Genotype, data = Fig_S1e)

# Kruskal–Wallis test ####
kruskal_result <- kruskal.test(AM_dTomato_Log2 ~ Genotype, data = Fig_S1e)

# Post-hoc - Dunn’s post-hoc test
posthoc_result <- Fig_S1e %>% 
  dunn_test(AM_dTomato_Log2 ~ Genotype, p.adjust.method = "bonferroni")

# Report results
sink(file = file.path(output, "FigS1e_Statistics.txt"))
cat("FigS1e AM_dTomato Statistics\n")
cat("\nAssumption Checks:\n")
cat("\nNormality of Residuals:\n")
print(shapiro_result)
cat("\nHomogeneity of Variances:\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nKruskal-Wallis Test Results:\n")
print(kruskal_result)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (Dunn’s post-hoc test with Bonferroni correction):\n")
print(posthoc_result)
cat("\n=====================\n")
sink()
