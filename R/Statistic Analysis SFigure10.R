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
# Fig_S10b####
Fig_S10b<- readxl::read_excel("Source Data.xlsx", sheet = "Fig_S10b")
head(Fig_S10b)
table(Fig_S10b$Genotype)
Fig_S10b$Genotype <- factor(Fig_S10b$Genotype, levels = c("WT","B41KO", "B41Tg"))

# Kruskal-Wallis test ####
# bounded data and small sample size (n = 3 per group)
KWH = kruskal.test(Bhlhe41 ~ Genotype, data = Fig_S10b)
# posthoc test (Dunn's Test) differences between genotype groups
dunn_test_result <- Fig_S10b %>% 
  dunn_test(Bhlhe41 ~ Genotype, p.adjust.method = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "FigS10b_Statistics.txt"))
cat("FigS10b Bhlhe41 Statistics\n")
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