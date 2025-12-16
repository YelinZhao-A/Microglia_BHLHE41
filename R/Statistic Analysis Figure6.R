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
# Fig_6c####
Fig_6c <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_6c")
head(Fig_6c)
table(Fig_6c$Genotype)
Fig_6c$Genotype <- factor(Fig_6c$Genotype, levels = c("B41HET", "B41KO","B41Tg", "WT"))

#CD22 ####
## One-way ANOVA ####
model <- aov(CD22 ~ Genotype, data = Fig_6c)
anova_result <- summary(model)
# Check ANOVA assumptions
# Normality of residuals
shapiro_result <- shapiro.test(residuals(model))
pdf("output/Diagnosis/Fig6c_CD22_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()
# Homogeneity of variances
Homogeneity <- car::leveneTest(CD22 ~ Genotype, data = Fig_6c)
# Post-hoc 
posthoc <- emmeans(model, pairwise ~ Genotype, adjust = "bonferroni")

# Report results
sink(file = file.path(output, "Fig6c_Statistics.txt"))
cat("Fig6c CD22 Statistics\n")
cat("\nOne-way ANOVA Results:\n")
print(anova_result)
cat("\n--------------------\n")
cat("\nAssumption Checks:\n")
cat("\nNormality of Residuals:\n")
print(shapiro_result)
cat("\nHomogeneity of Variances:\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (with bonferroni correction):\n")
print(posthoc)
cat("\n=====================\n")
sink()

#dTomato ####
## One-way ANOVA ####
model <- aov(dTomato %>% log2 ~ Genotype, data = Fig_6c)
anova_result <- summary(model)
# Check ANOVA assumptions
# Normality of residuals
shapiro_result <- shapiro.test(residuals(model))
pdf("output/Diagnosis/Fig6c_dTomato_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()
# Homogeneity of variances
Homogeneity <- car::leveneTest(dTomato ~ Genotype, data = Fig_6c)
# Post-hoc 
posthoc <- emmeans(model, pairwise ~ Genotype, adjust = "bonferroni")

# Report results
sink(file = file.path(output, "Fig6c_Statistics.txt"),append=TRUE)
cat("Fig6c dTomato Statistics\n")
cat("\nOne-way ANOVA Results - log2(dTomato):\n")
print(anova_result)
cat("\n--------------------\n")
cat("\nAssumption Checks:\n")
cat("\nNormality of Residuals:\n")
print(shapiro_result)
cat("\nHomogeneity of Variances:\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (bonferroni):\n")
print(posthoc)
cat("\n=====================\n")
sink()


###################################
# Fig_6e####
Fig_6e <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_6e")
head(Fig_6e)
table(Fig_6e$Genotype)
Fig_6e$Genotype <- factor(Fig_6e$Genotype, levels = c("B41HET", "B41Tg","KO", "KO_B41Tg"))

#CD22 ####
## One-way ANOVA ####
model <- aov(CD22 %>% log2 ~ Genotype, data = Fig_6e)
anova_result <- summary(model)
# Check ANOVA assumptions
# Normality of residuals
shapiro_result <- shapiro.test(residuals(model))
pdf("output/Diagnosis/Fig6e_CD22_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(model, main = "Q-Q Plot")
dev.off()
# Homogeneity of variances
Homogeneity <- car::leveneTest(CD22 ~ Genotype, data = Fig_6e)
# Post-hoc
posthoc <- emmeans(model, pairwise ~ Genotype, adjust = "bonferroni")

# Report results
sink(file = file.path(output, "Fig6e_Statistics.txt"))
cat("Fig6e CD22 Statistics\n")
cat("\nOne-way ANOVA Results - log2(CD22):\n")
print(anova_result)
cat("\n--------------------\n")
cat("\nAssumption Checks:\n")
cat("\nNormality of Residuals:\n")
print(shapiro_result)
cat("\nHomogeneity of Variances:\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (bonferroni):\n")
print(posthoc)
cat("\n=====================\n")
sink()

# dTomato ####
## One-way ANOVA ####
model <- aov(dTomato %>% log2 ~ Genotype, data = Fig_6e)
anova_result <- summary(model)
# Check ANOVA assumptions
# Normality of residuals
shapiro_result <- shapiro.test(residuals(model))
pdf("output/Diagnosis/Fig6e_dTomato_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(model, main = "Q-Q Plot")
dev.off()
# Homogeneity of variances
Homogeneity <- car::leveneTest(dTomato ~ Genotype, data = Fig_6e)
# Post-hoc
posthoc <- emmeans(model, pairwise ~ Genotype, adjust = "bonferroni")

# Report results
sink(file = file.path(output, "Fig6e_Statistics.txt"),append=TRUE)
cat("Fig6e dTomato Statistics\n")
cat("\nOne-way ANOVA Results - log2(dTomatos):\n")
print(anova_result)
cat("\n--------------------\n")
cat("\nAssumption Checks:\n")
cat("\nNormality of Residuals:\n")
print(shapiro_result)
cat("\nHomogeneity of Variances:\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (bonferroni):\n")
print(posthoc)
cat("\n=====================\n")
sink()

###################################
# Fig_6g####
Fig_6g <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_6g")
head(Fig_6g)
table(Fig_6g$Genotype)
Fig_6g$Genotype <- factor(Fig_6g$Genotype, levels = c("KO", "KO_B41Tg"))

## Check data distribution####
by(Fig_6g$Mye_Dil, Fig_6g$Genotype, shapiro.test)

pdf("output/Diagnosis/Fig6g_QQ_Plot.pdf", width = 6, height = 6)  
qqnorm(Fig_6g$Mye_Dil, main = "Q-Q Plot of AF488_Pos")  
qqline(Fig_6g$Mye_Dil, lwd = 2)  
dev.off()

# Normality of Residuals
Homogeneity = leveneTest(Mye_Dil ~ Genotype, data = Fig_6g)

## t test ####
if (Homogeneity$`Pr(>F)`[1] > 0.05) {
  t_test_result <- t.test(Mye_Dil ~ Genotype, data = Fig_6g, var.equal = TRUE)
} else {
  t_test_result <- t.test(Mye_Dil ~ Genotype, data = Fig_6g, var.equal = FALSE)
}

sample_estimates <- t_test_result$estimate
conf_int_lower <- t_test_result$conf.int[1]
conf_int_upper <- t_test_result$conf.int[2]

result_df <- data.frame(
  Variable = deparse(substitute(Mye_Dil)), 
  Test = ifelse(Homogeneity$`Pr(>F)`[1] > 0.05, "Student's t-test", "Welch's t-test"),
  t_statistic = t_test_result$statistic,
  df = t_test_result$parameter,
  p_value = t_test_result$p.value,
  Mean_Group1 = sample_estimates[1],
  Mean_Group2 = sample_estimates[2],
  CI_Lower = conf_int_lower,
  CI_Upper = conf_int_upper,
  stringsAsFactors = FALSE
)%>% t()

group_names <- names(sample_estimates)
rownames(result_df)[6] <- paste0("Mean_", group_names[1])
rownames(result_df)[7] <- paste0("Mean_", group_names[2])

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig6g_Statistics.txt"))
cat("Fi6g Mye_Dil Statistics\n")
cat("\n--------------------\n")
cat("\nCheck assumptions:\n")
cat("\nNormal distribution:\n")
print(by(Fig_6g$Mye_Dil, Fig_6g$Genotype, shapiro.test))
cat("\nHomogeneity of Variance (Levene's Test):\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nt-test Results:\n")
print(result_df)
cat("\n=====================\n")
sink()
