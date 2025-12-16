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

#####################################
# Fig_S6a####
Fig_S6a <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S6ab")
head(Fig_S6a)
table(Fig_S6a$Genotype, Fig_S6a$CNS_Region)

Fig_S6a$Genotype <- factor(Fig_S6a$Genotype, levels = c("WT","AurkacHET","AurkacKO"))
Fig_S6a$CNS_Region <- factor(Fig_S6a$CNS_Region, levels = c("Cortex", "Hippo", "SVZ_CC"))

## S6a_Microglia ####
## Generalized Linear Model (GLM) ######
# For count data
model_poisson <- glm(S6a_Microglia ~ Genotype * CNS_Region, data = Fig_S6a, family = poisson(link = "log"))

# Check overdispersion
overdispersion_poisson <- sum(residuals(model_poisson, type = "pearson")^2) / model_poisson$df.residual

# posthoc test differences between genotype groups
posthoc <- emmeans(model_poisson, pairwise ~ Genotype | CNS_Region, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS6a_Statistics.txt"))
cat("FigS6a Statistics\n")
cat("\nData type: count\n")
cat("\nGeneralized Linear Model (GLM) with poisson. Results:\n")
print(summary(model_poisson))
print(anova(model_poisson, test="Chisq"))
cat("\n--------------------\n")
cat("\nOverdispersion Check:\n")
cat(paste("Overdispersion ratio:", round(overdispersion_poisson, 2), "\n"))
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (EMMeans with Bonferroni adjustment):\n")
print(posthoc$contrasts)
cat("\n=====================\n")
sink()

## S6a_BrdU_Pos ####
## Generalized Linear Model (GLM) ######
# For count data
model_brdU <- glm(S6a_BrdU_Pos ~ Genotype * CNS_Region, data = Fig_S6a, family = poisson(link = "log"))
# Check overdispersion
overdispersion_brdU <- sum(residuals(model_brdU, type = "pearson")^2) / model_brdU$df.residual
# posthoc test differences between genotype groups
posthoc <- emmeans(model_brdU, pairwise ~ Genotype | CNS_Region, adjust = "bonferroni")
# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS6a_Statistics.txt"), append = TRUE)
cat("FigS6a BrdU_Pos Statistics\n")
cat("\nData type: count\n")
cat("\nGeneralized Linear Model (GLM) with poisson. Results:\n")
print(summary(model_brdU))
print(anova(model_brdU, test="Chisq"))
cat("\n--------------------\n")
cat("\nOverdispersion Check:\n")
cat(paste("Overdispersion ratio:", round(overdispersion_brdU, 2), "\n"))
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (EMMeans with Bonferroni adjustment):\n")
print(posthoc$contrasts)
cat("\n=====================\n")
sink()
###################################

# Fig_S6b####
Fig_S6b <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S6ab")
head(Fig_S6b)
table(Fig_S6b$Genotype, Fig_S6b$CNS_Region)
Fig_S6b$Genotype <- factor(Fig_S6b$Genotype, levels = c("WT","AurkacHET","AurkacKO"))
Fig_S6b$CNS_Region <- factor(Fig_S6b$CNS_Region, levels = c("Cortex", "Hippo", "SVZ_CC"))
Fig_S6b$S6b_TUNEL_Pos_log2 <- log2(Fig_S6b$S6b_TUNEL_Pos)

## S6b_TUNEL_Pos ####
## Two-Way ANOVA ######
model <- aov(S6b_TUNEL_Pos_log2 ~ Genotype * CNS_Region, data = Fig_S6b)
anova_result <- summary(model)
# Check ANOVA assumptions
# Normality of residuals
shapiro_result <- shapiro.test(residuals(model))
pdf("output/Diagnosis/FigS6b_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()
# Homogeneity of variances
levene_results <- levene_test(S6b_TUNEL_Pos_log2 ~ Genotype * CNS_Region, data = Fig_S6b)
# Post-hoc test differences between genotype groups
posthoc <- emmeans(model, pairwise ~ Genotype | CNS_Region, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS6b_Statistics.txt"))
cat("FigS6b Statistics\n")
cat("\nData type: continuous (log2 transformed)\n")
cat("\nTwo-Way ANOVA Results:\n")
model
print(anova_result)
cat("\n--------------------\n")
cat("\nAssumption Checks:\n")
cat("\nNormality of Residuals (Shapiro-Wilk Test):\n")
print(shapiro_result)
cat("\nHomogeneity of Variances (Levene's Test):\n")
print(levene_results)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (EMMeans with Bonferroni adjustment):\n")
print(posthoc$contrasts)
cat("\n=====================\n")
sink()

################################
# Fig_S6c####
Fig_S6c <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S6c")
head(Fig_S6c)
table(Fig_S6c$Genotype)
Fig_S6c$Genotype <- factor(Fig_S6c$Genotype, levels = c("WT","AurkacHET","AurkacKO"))

# One-way ANOVA ######
model <- aov(pH3_pos ~ Genotype, data = Fig_S6c)
anova_result <- summary(model)
# Check ANOVA assumptions
# Normality of residuals
shapiro_result <- shapiro.test(residuals(model))
pdf("output/Diagnosis/FigS6c_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()
# Homogeneity of variances
levene_results <- levene_test(pH3_pos ~ Genotype, data = Fig_S6c)
# Post-hoc test differences between genotype groups
posthoc <- emmeans(model, pairwise ~ Genotype, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS6c_Statistics.txt"))
cat("FigS6c Statistics\n")
cat("\nData type: continuous\n")
cat("\nOne-Way ANOVA Results:\n")
print(model)
print(anova_result)
cat("\n--------------------\n")
cat("\nAssumption Checks:\n")
cat("\nNormality of Residuals (Shapiro-Wilk Test):\n")
print(shapiro_result)
cat("\nHomogeneity of Variances (Levene's Test):\n")
print(levene_results)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (EMMeans with Bonferroni adjustment):\n")
print(posthoc$contrasts)
cat("\n=====================\n")
sink()
