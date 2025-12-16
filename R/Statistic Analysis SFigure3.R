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
# Fig_S3a####
Fig_S3a <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S3a")
head(Fig_S3a)
table(Fig_S3a$Genotype, Fig_S3a$CNS_Region)

Fig_S3a$Genotype <- factor(Fig_S3a$Genotype, levels = c("WT","B41HET","B41KO"))
Fig_S3a$CNS_Region <- factor(Fig_S3a$CNS_Region, levels = c('Cortex', 'CC', 'CA1', 'CA3' ,  'DG'))

## Generalized Linear Model (GLM)  ######
model_poisson <- glm(Microglia ~ CNS_Region * Genotype, data = Fig_S3a, family = poisson(link = "log"))
anova(model_poisson, test="Chisq")

# Check overdispersion
overdispersion <- sum(residuals(model_poisson, type = "pearson")^2) / model_poisson$df.residual

# posthoc test differences between genotype groups
posthoc <- emmeans(model_poisson, pairwise ~ Genotype | CNS_Region, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS3a_Statistics.txt"))
cat("FigS3a Statistics\n")
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
# Fig_S3c####
Fig_S3c <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S3c")
head(Fig_S3c)
table(Fig_S3c$Genotype, Fig_S3c$CNS_Region)

Fig_S3c$Genotype <- factor(Fig_S3c$Genotype, levels = c("WT","B41cHET","B41cKO"))
Fig_S3c$CNS_Region <- factor(Fig_S3c$CNS_Region, levels = c('Cortex', 'CC', 'CA1', 'CA3' ,  'DG'))

## Generalized Linear Model (GLM) with negative binomial ######
model_poisson <- glm(Microglia ~ CNS_Region * Genotype, data = Fig_S3c, family = poisson(link = "log"))
anova(model_poisson, test="Chisq")

# Check overdispersion
overdispersion <- sum(residuals(model_poisson, type = "pearson")^2) / model_poisson$df.residual

# posthoc test differences between genotype groups
posthoc <- emmeans(model_poisson, pairwise ~ Genotype | CNS_Region, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS3c_Statistics.txt"))
cat("FigS3c Statistics\n")
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


