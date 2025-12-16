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
# Fig_S11a####
Fig_S11a <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S11a")
head(Fig_S11a)
table(Fig_S11a$Genotype,Fig_S11a$CNS_Region)
Fig_S11a$Genotype <- factor(Fig_S11a$Genotype, levels = c("B41Tg", "B41HET","B41KO"))
Fig_S11a$CNS_Region <- factor(Fig_S11a$CNS_Region, levels = c("Cortex","CC","CA1","CA3","DG"))

## GLM with Poisson distribution ####
# count data, use generalized linear model
model <- glm(Microglia ~ Genotype  *  CNS_Region, data = Fig_S11a, family = poisson(link = "log"))
anova_summary <- anova(model, test="Chisq")

# Check for overdispersion
overdisp_flag <- FALSE
dispersion <- sum(residuals(model, type="pearson")^2) / model$df.residual
if (dispersion > 1.2) {
  overdisp_flag <- TRUE
  print("Overdispersion detected → Switching to Negative Binomial GLM")
  model <- MASS::glm.nb(Microglia ~ Genotype + CNS_Region, data = Fig_S11a)
  anova_summary <- anova(model, test="Chisq")
  print(summary(model))
  dispersion <- sum(residuals(model, type="pearson")^2) / model$df.residual
  print(dispersion)
}

# Post-hoc test (emmeans) differences
posthoc_emm <- emmeans(model, pairwise ~ Genotype | CNS_Region, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS11a_Statistics.txt"))
cat("FigS11a Microglia Statistics\n")
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

