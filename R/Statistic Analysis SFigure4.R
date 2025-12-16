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
# Fig_S4a####
Fig_S4a <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S4a")
head(Fig_S4a)
table(Fig_S4a$Genotype, Fig_S4a$Age)
Fig_S4a$Genotype <- factor(Fig_S4a$Genotype, levels = c("WT","AurkaHET","AurkacKO"))
Fig_S4a$Age <- factor(Fig_S4a$Age, levels = c('P14', 'P21', 'P60'))
Fig_S4a$dTomato_gMFI_log2 <- log2(Fig_S4a$dTomato_gMFI)

## Two-way ANOVA####
model <- aov(dTomato_gMFI_log2 ~ Genotype * Age, data = Fig_S4a)
anova_result <- summary(model)
# Check ANOVA assumptions
# Normality of residuals
shapiro_result <- shapiro.test(residuals(model))
pdf("output/Diagnosis/FigS4a_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()
# Homogeneity of variances
Homogeneity <- car::leveneTest(dTomato_gMFI_log2 ~ Genotype * Age, data = Fig_S4a)

# Post-hoc test differences between genotype groups
posthoc <- emmeans(model, pairwise ~ Genotype | Age, adjust = "bonferroni")

# Report results
sink(file = file.path(output, "FigS4a_Statistics.txt"))
cat("FigS4a dTomato_gMFI Statistics\n")
cat("\nTwo-way ANOVA Results:\n")
print(model)
print(anova_result)
cat("\n--------------------\n")
cat("\nAssumption Checks:\n")
cat("\nNormality of Residuals:\n")
print(shapiro_result)
cat("\nHomogeneity of Variances:\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (EMMeans with Bonferroni adjustment):\n")
print(posthoc$contrasts)
cat("\n=====================\n")
sink()

#####################################
# Fig_S4b####
Fig_S4b <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S4b")
head(Fig_S4b)
table(Fig_S4b$Genotype)
Fig_S4b$Genotype <- factor(Fig_S4b$Genotype, levels = c("WT","AurkacHET","AurkacKO"))

## Aurka####
## kruskal–Wallis test ####
KWH = kruskal.test(Aurka ~ Genotype, data = Fig_S4b)
# posthoc test (Dunn's Test) differences between genotype groups
posthoc <- Fig_S4b %>% 
  dunn_test(Aurka ~ Genotype, p.adjust.method = "bonferroni")

# Report results
sink(file = file.path(output, "FigS4b_Statistics.txt"))
cat("FigS4b Aurka Statistics\n")
cat("\nKruskal-Wallis H test (bounded data):\n")
print(KWH)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (Dunn's Test with Bonferroni adjustment):\n")
print(posthoc)
cat("\n=====================\n")
sink()

##Bhlhe41####
## kruskal–Wallis test ####
KWH = kruskal.test(Bhlhe41 ~ Genotype, data = Fig_S4b)
# posthoc test (Dunn's Test) differences between genotype groups
posthoc <- Fig_S4b %>% 
  dunn_test(Bhlhe41 ~ Genotype, p.adjust.method = "bonferroni")
# Report results
sink(file = file.path(output, "FigS4b_Statistics.txt"),append=TRUE)
cat("FigS4b Bhlhe41 Statistics\n")
cat("\nKruskal-Wallis H test (bounded data):\n")
print(KWH)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (Dunn's Test with Bonferroni adjustment):\n")
print(posthoc)
cat("\n=====================\n")
sink()

#####################################
# Fig_S4c####
Fig_S4c <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S4c")
head(Fig_S4c)
table(Fig_S4c$Genotype)
Fig_S4c$Genotype <- factor(Fig_S4c$Genotype, levels = c("WT","AurkbcHET","AurkbcKO"))

## dTomato_gMFI####
## kruskal–Wallis test ####
KWH = kruskal.test(dTomato_gMFI ~ Genotype, data = Fig_S4c)
# posthoc test (Dunn's Test) differences between genotype groups
posthoc <- Fig_S4c %>% 
  dunn_test(dTomato_gMFI ~ Genotype, p.adjust.method = "bonferroni")
# Report results
sink(file = file.path(output, "FigS4c_Statistics.txt"))
cat("FigS4c dTomato_gMFI Statistics\n")
cat("\nKruskal-Wallis H test (bounded data):\n")
print(KWH)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (Dunn's Test with Bonferroni adjustment):\n")
print(posthoc)
cat("\n=====================\n")
sink()




