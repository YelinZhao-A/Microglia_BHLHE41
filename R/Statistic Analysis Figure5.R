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

###########################################
# Fig_5b ####
Fig_5b <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_5b")
head(Fig_5b)
table(Fig_5b$Genotype, Fig_5b$Age)

Fig_5b$CD22_log2 <- log2(Fig_5b$CD22)
Fig_5b$CD22_sqr <- sqrt(Fig_5b$CD22)
Fig_5b$Genotype <- factor(Fig_5b$Genotype, levels = c("WT", "B41HET","B41KO"))
Fig_5b$Age <- factor(Fig_5b$Age, levels = c("Young", "Middle_Aged"))

## Two way anova####
model <- aov(CD22_log2 ~ Genotype * Age, data = Fig_5b)
anova_result <- summary(model)

# Normality of Residuals
Normality <- shapiro.test(residuals(model))
pdf(file = file.path(output, "Diagnosis/Fig5b_Residuals_QQPlot.pdf"))
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()

# Check homogeneity of variances
Homogeneity <- car::leveneTest(CD22_log2 ~ Genotype * Age, data = Fig_5b)

## Switch to ART ANOVA ####
# The QQ plot revealed a significant departure from normality. Therefore, we proceeded with Aligned Rank Transform (ART) ANOVA.
model_art <- art(CD22  ~ Genotype * Age, data = Fig_5b)
art_anova_result <- anova(model_art)

# posthoc test differences between genotype groups
emm <- emmeans(artlm(model_art, "Genotype:Age"), ~ Genotype : Age)

# Pairwise comparisons of Genotype within each condition
posthoc_results <- pairs(emm, adjust = "bonferroni")
print(posthoc_results)

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig5b_Statistics.txt"))
cat("Fig5b CD22_log2 Statistics\n")
cat("\nART ANOVA Results:\n")
print(art_anova_result)
cat("\n--------------------\n")
cat("\nAssumption Checks:\n")
cat("\nNormality of Residuals:\n")
print(Normality)
cat("\nHomogeneity of Variances:\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (ART Contrasts with Bonferroni adjustment):\n")
print(posthoc_results)
cat("\n=====================\n")
sink()

############################################################
# Fig_5d ####
Fig_5d <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_5d")
head(Fig_5d)
table(Fig_5d$Genotype)

Fig_5d$Genotype <- factor(Fig_5d$Genotype, levels = c("WT", "cKO"))

## Check data distribution####
by(Fig_5d$CD22, Fig_5d$Genotype, shapiro.test)

pdf("output/Diagnosis/Fig5d_QQ_Plot.pdf", width = 6, height = 6)  
qqnorm(Fig_5d$CD22, main = "Q-Q Plot of AF488_Pos")  
qqline(Fig_5d$CD22, lwd = 2)  
dev.off()

## Wilcoxon rank sum test ####
wilcox_test_result <- wilcox.test(CD22 ~ Genotype, data = Fig_5d, exact = FALSE)

# Report results
sink(file = file.path(output, "Fig5d_Statistics.txt"))
cat("Fig5d CD22 Statistics\n")
cat("\n Non-normal / very small n\n")
cat("\nWilcoxon Rank Sum Test Results:\n")
print(wilcox_test_result)
cat("\n=====================\n")
sink()

############################
# Fig_5f ####
Fig_5f <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_5f")
head(Fig_5f)
table(Fig_5f$Genotype, Fig_5f$Age)
Fig_5f$Genotype <- factor(Fig_5f$Genotype, levels = c("HET", "cKO"))
Fig_5f$Age <- factor(Fig_5f$Age, levels = c("Young", "Middle_Aged"))
Fig_5f$CD22_log2 <- log2(Fig_5f$CD22)
## Two way anova####
model <- aov(CD22_log2 ~ Genotype * Age, data = Fig_5f)
anova_result <- summary(model)
# Normality of Residuals
Normality = shapiro.test(residuals(model))
pdf(file = file.path(output, "Diagnosis/Fig5f_Residuals_QQPlot.pdf"))
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()
# Check homogeneity of variances
Homogeneity <- car::leveneTest(CD22_log2 ~ Genotype * Age, data = Fig_5f)

## Switch to ART ANOVA ####
# Consider 4 samples per group combination. We proceeded with non-parametric test: Aligned Rank Transform (ART) ANOVA.
model_art <- art(CD22  ~ Genotype * Age, data = Fig_5f)
art_anova_result <- anova(model_art)

# posthoc test differences between genotype groups
emm <- emmeans(artlm(model_art, "Genotype:Age"), ~ Genotype : Age)

# Pairwise comparisons of Genotype within each condition
posthoc_results <- pairs(emm, adjust = "bonferroni")
print(posthoc_results)

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig5f_Statistics.txt"))
cat("Fig5f CD22_log2 Statistics\n")
cat("\nART ANOVA Results:\n")
print(art_anova_result)
cat("\n--------------------\n")
cat("\nAssumption Checks:\n")
cat("\nNormality of Residuals:\n")
print(Normality)
cat("\n The Q-Q plot indicate a significant departure from normality.\n")
cat("\nHomogeneity of Variances:\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (ART Contrasts with Bonferroni adjustment):\n")
print(posthoc_results)
cat("\n=====================\n")
sink()

###################################
# Fig_5i ####
Fig_5i <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_5i")
head(Fig_5i)
table(Fig_5i$Genotype)
Fig_5i$Genotype <- factor(Fig_5i$Genotype, levels = c("WT", "cKO"))

## Check data distribution####
by(Fig_5i$CD22, Fig_5i$Genotype, shapiro.test)
pdf("output/Diagnosis/Fig5i_QQ_Plot.pdf", width = 6, height = 6)
qqnorm(Fig_5i$CD22)
qqline(Fig_5i$CD22, lwd = 2)
dev.off()

## Wilcoxon rank sum test ####
wilcox_test_result <- wilcox.test(CD22 ~ Genotype, data = Fig_5i, exact = FALSE)
# Report results
sink(file = file.path(output, "Fig5i_Statistics.txt"))
cat("Fig5i CD22 Statistics\n")
cat("\n Non-normal / very small n\n")
cat("\nWilcoxon Rank Sum Test Results:\n")
print(wilcox_test_result)
cat("\n=====================\n")
sink()
###################################
# Fig_5m#####
Fig_5m <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_5m")
head(Fig_5m)
table(Fig_5m$Treatment)

Fig_5m$Treatment <- factor(Fig_5m$Treatment, levels = c("Anti_Bhlhe41", "IgG"))

# Wilcoxon signed rank test ####
wilcox_test_result_P1 <- wilcox.test(P1 ~ Treatment, data = Fig_5m,  exact = FALSE)
wilcox_test_result_P2 <- wilcox.test(P2 ~ Treatment, data = Fig_5m,  exact = FALSE)
wilcox_test_result_P3 <- wilcox.test(P3 ~ Treatment, data = Fig_5m,exact = FALSE)
wilcox_test_result_NC <- wilcox.test(NC ~ Treatment, data = Fig_5m, exact = FALSE)

# Report results
sink(file = file.path(output, "Fig5m_Statistics.txt"))
cat("Fig5m Bhlhe41 Statistics\n")
cat("\n Wilcoxon Rank Sum Test Results:\n")
cat("\nP1:\n")
print(wilcox_test_result_P1)
cat("\nP2:\n")
print(wilcox_test_result_P2)
cat("\nP3:\n")
print(wilcox_test_result_P3)
cat("\nNC:\n")
print(wilcox_test_result_NC)
cat("\n=====================\n")
sink()
