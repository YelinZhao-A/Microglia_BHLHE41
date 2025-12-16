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
#Fig_S5a ####
Fig_S5a <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S5a")
head(Fig_S5a)
table(Fig_S5a$Genotype, Fig_S5a$CNS_Region)
Fig_S5a$Genotype <- factor(Fig_S5a$Genotype, levels = c("WT","AurkacHET","AurkacKO"))
Fig_S5a$CNS_Region <- factor(Fig_S5a$CNS_Region, levels = c('Cortex', 'CC', 'CA1', 'CA3' ,  'DG'))

## Generalized Linear Model (GLM)  ######
# For count data
model_poisson <- glm(Microglia ~ CNS_Region * Genotype, data = Fig_S5a, family = poisson(link = "log"))
anova(model_poisson, test="Chisq")

# Check overdispersion
overdispersion <- sum(residuals(model_poisson, type = "pearson")^2) / model_poisson$df.residual

# posthoc test differences between genotype groups
posthoc <- emmeans(model_poisson, pairwise ~ Genotype | CNS_Region, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS5a_Statistics.txt"))
cat("FigS5a Statistics\n")
cat("\nData type: count\n")
cat("\nGeneralized Linear Model (GLM) with Poisson:\n")
print(summary(model_poisson))
print(anova(model_poisson, test="Chisq"))
cat("\n--------------------\n")
cat("\nOverdispersion Check:\n")
cat(paste("Overdispersion ratio:", round(overdispersion, 2), "\n"))
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (EMMeans with Bonferroni adjustment):\n")
print(posthoc$contrasts)
cat("\n=====================\n")
sink()

#######################################
# Fig_S5b####
Fig_S5b <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S5b")
head(Fig_S5b)
table(Fig_S5b$Genotype, Fig_S5b$CNS_Region)
Fig_S5b$Genotype <- factor(Fig_S5b$Genotype, levels = c("WT","AurkacHET","AurkacKO"))
Fig_S5b$CNS_Region <- factor(Fig_S5b$CNS_Region, levels = c('Cortex', 'CC', 'CA1', 'CA3' ,  'DG'))

## Generalized Linear Model (GLM) with negative binomial ######
# For count data
model_poisson <- glm(Microglia ~ CNS_Region * Genotype, data = Fig_S5b, family = poisson(link = "log"))
anova(model_poisson, test="Chisq")
# Check overdispersion
overdispersion <- sum(residuals(model_poisson, type = "pearson")^2) / model_poisson$df.residual
# posthoc test differences between genotype groups
posthoc <- emmeans(model_poisson, pairwise ~ Genotype | CNS_Region, adjust = "bonferroni")
# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS5b_Statistics.txt"))
cat("FigS5b Statistics\n")
cat("\nData type: count\n")
cat("\nGeneralized Linear Model (GLM) with poisson:\n")
print(summary(model_poisson))
print(anova(model_poisson, test="Chisq"))
cat("\n--------------------\n")
cat("\nOverdispersion Check:\n")
cat(paste("Overdispersion ratio:", round(overdispersion, 2), "\n"))
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (EMMeans with Bonferroni adjustment):\n")
print(posthoc$contrasts)
cat("\n=====================\n")
sink()
#####################################




