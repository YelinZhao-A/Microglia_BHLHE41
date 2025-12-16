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
# Fig_S9a####
Fig_S9a <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S9a")
head(Fig_S9a)
table(Fig_S9a$Genotype)
Fig_S9a$Genotype <- factor(Fig_S9a$Genotype, levels = c("WT", "AurkacHET","AurkacKO"))

# Kruskal-Wallis test ####
# bounded data and small sample size (n = 3 per group)
KWH = kruskal.test(Aurka ~ Genotype, data = Fig_S9a)
# posthoc test (Dunn's Test) differences between genotype groups
dunn_test_result <- Fig_S9a %>% 
  dunn_test(Aurka ~ Genotype, p.adjust.method = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS9a_Statistics.txt"))
cat("FigS9a Aurka Statistics\n")
cat("\nKruskal-Wallis H test:\n")
cat("\n bounded data and small sample size (n = 3 per group), 
    use non-parametric Kruskal-Wallis H test .\n")
print(KWH)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (Dunn's Test with Bonferroni adjustment):\n")
print(dunn_test_result)
cat("\n===============================================\n")
sink()
########################################

# Fig_S9b####
Fig_S9b <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S9b")
head(Fig_S9b)
table(Fig_S9b$Genotype, Fig_S9b$CNS_Region)
Fig_S9b$Genotype <- factor(Fig_S9b$Genotype, levels = c( "AurkacHET","AurkacKO"))
Fig_S9b$CNS_Region <- factor(Fig_S9b$CNS_Region, levels = c("Cortex","CC","CA1","CA3","DG"))

# count data, use generalized linear model
## GLM with Poisson distribution ####
model <- glm(Microglia ~ Genotype  *  CNS_Region, data = Fig_S9b, family = poisson(link = "log"))
anova_summary <- anova(model, test="Chisq")

# Check for overdispersion
overdisp_flag <- FALSE
dispersion <- sum(residuals(model, type="pearson")^2) / model$df.residual
if (dispersion > 1.2) {
  overdisp_flag <- TRUE
  model <- MASS::glm.nb(Microglia ~ Genotype + CNS_Region, data = Fig_S9b)
  anova_summary <- anova(model, test="Chisq")
  print(summary(model))
  dispersion <- sum(residuals(model, type="pearson")^2) / model$df.residual
  print(dispersion)
}

# Post-hoc test (emmeans) differences 
posthoc_emm <- emmeans(model, pairwise ~ Genotype | CNS_Region, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS9b_Statistics.txt"))
cat("FigS9b Microglia Statistics\n")
cat("\n--Count data------\n")
cat("\nGLM with Poisson distribution:\n")
cat("\nCheck for overdispersion:\n")
cat(paste0("Dispersion ratio: ", round(dispersion, 2), "\n"))
if (overdisp_flag) {
  cat("Overdispersion detected → Switched to Negative Binomial GLM\n")
} else {
  cat("No overdispersion detected.\n")
}
print(summary(model))
print(anova_summary)
cat("\n-------------\n")
cat("\nPost-hoc Test Results (emmeans with Bonferroni adjustment):\n")
print(posthoc_emm$contrasts)
sink()


########################################
# Fig_S9c####
Fig_S9c <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S9c")
head(Fig_S9c)
table(Fig_S9b$Genotype, Fig_S9b$CNS_Region)
Fig_S9c$Genotype <- factor(Fig_S9c$Genotype, levels = c( "AurkacHET","AurkacKO"))
Fig_S9c$CNS_Region <- factor(Fig_S9c$CNS_Region, levels = c("Cortex","CC","CA1","CA3","DG"))

# count data, use generalized linear model
## GLM with Poisson distribution ####
model <- glm(Microglia ~ Genotype  *  CNS_Region, data = Fig_S9c, family = poisson(link = "log"))
anova_summary <- anova(model, test="Chisq")

# Check for overdispersion
overdisp_flag <- FALSE
dispersion <- sum(residuals(model, type="pearson")^2) / model$df.residual
if (dispersion > 1.2) {
  overdisp_flag <- TRUE
  print("Overdispersion detected → Switching to Negative Binomial GLM")
  model <- MASS::glm.nb(Microglia ~ Genotype + CNS_Region, data = Fig_S9c)
  anova_summary <- anova(model, test="Chisq")
  print(summary(model))
  dispersion <- sum(residuals(model, type="pearson")^2) / model$df.residual
  print(dispersion)
}

# Post-hoc test (emmeans) differences
posthoc_emm <- emmeans(model, pairwise ~ Genotype | CNS_Region, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS9c_Statistics.txt"))
cat("FigS9c Microglia Statistics\n")
cat("\n--Count data------\n")
cat("\nGLM with Poisson distribution:\n")
cat("\nCheck for overdispersion:\n")
cat(paste0("Dispersion ratio: ", round(dispersion, 2), "\n"))
if (overdisp_flag) {
  cat("Overdispersion detected → Switched to Negative Binomial GLM\n")
} else {
  cat("No overdispersion detected.\n")
}
print(summary(model))
print(anova_summary)
cat("\n-------------\n")
cat("\nPost-hoc Test Results (emmeans with Bonferroni adjustment):\n")
print(posthoc_emm$contrasts)
sink()
########################################





