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

# Fig1a ####
## Load data####
Fig_1a <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_1a")
head(Fig_1a)

Fig_1a$Genotype <- factor(Fig_1a$Genotype, levels = c("WT", "B41HET", "B41KO"))
Fig_1a$Age <- factor(Fig_1a$Age)
Fig_1a$Group_Combined <- interaction(Fig_1a$Genotype, Fig_1a$Age, sep="_")

## Check data distribution####
by(Fig_1a$dTomato_gMFI %>% log2(), Fig_1a$Group_Combined, shapiro.test)

## Two-way ANOVA####
model <- aov(log2(dTomato_gMFI) ~ Genotype * Age, data = Fig_1a)
summary(model)

## Check assumptions ####
# Normality of Residuals
Normality = shapiro.test(residuals(model)) 
pdf(file = file.path(output, "Diagnosis/Fig1a_Residuals_QQPlot.pdf"))
plot(model, 2)
dev.off()

# Homogeneity of Variance
Homogeneity = leveneTest(log2(dTomato_gMFI) ~ Genotype * Age, data = Fig_1a)

# posthoc test differences between genotype groups
posthoc <- emmeans(model, pairwise ~ Genotype | Age, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig1a_Statistics.txt"))
cat("Fig1a Statistics\n")
cat("\nTwo-way ANOVA Results:\n")
print(summary(model))
cat("\n--------------------\n")
cat("\nCheck assumptions:\n")
cat("\nNormality of Residuals (Shapiro-Wilk test):\n")
print(Normality)
cat("\nHomogeneity of Variance (Levene's Test):\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (EMMeans with Bonferroni adjustment):\n")
print(posthoc$contrasts)
sink()

############################################################
# Fig1b ####
## Load data####
Fig_1b <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_1b")
head(Fig_1b)

Fig_1b$Genotype <- factor(Fig_1b$Genotype, levels = c("WT", "B41HET", "B41KO"))

## One-way ANOVA####
model <- aov(Bhlhe41_Exp~ Genotype , data = Fig_1b)
summary(model)

## Check assumptions ####
# Normality of Residuals
Normality = shapiro.test(residuals(model))
pdf(file = file.path(output, "Diagnosis/Fig1b_Residuals_QQPlot.pdf"))
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()
# Homogeneity of Variance
Homogeneity = leveneTest(Bhlhe41_Exp ~ Genotype , data = Fig_1b)

## Switch to Kruskal-Wallis H test ####
# Although the Shapiro-Wilk test was non-significant ($P > 0.05$), likely due to low statistical power ($N=4/5$ per group), 
#the QQ plot revealed a significant departure from normality, particularly in the tails. 
# Therefore, we proceeded with the more robust non-parametric Kruskal-Wallis H test to compare the groups.
KWH = kruskal.test(Bhlhe41_Exp ~ Genotype, data = Fig_1b)

# posthoc test (Dunn's Test) differences between genotype groups
dunn_test_result <- Fig_1b %>% 
  dunn_test(Bhlhe41_Exp ~ Genotype, p.adjust.method = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig1b_Statistics.txt"))
cat("Fig1b Statistics\n")
cat("\nKruskal-Wallis H test:\n")
cat("\nAlthough the Shapiro-Wilk test was non-significant (P > 0.05), likely due to low statistical power (N=4/5 per group), 
    the QQ plot revealed a significant departure from normality, particularly in the tails.Therefore, we proceeded with the more robust non-parametric Kruskal-Wallis H test to compare the groups.\n")
print(KWH)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (Dunn's Test with Bonferroni adjustment):\n")
print(dunn_test_result)
sink()
