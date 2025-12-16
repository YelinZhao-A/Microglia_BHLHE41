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
################################################
# Fig4a ####
Fig_4a <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_4a")
head(Fig_4a)
table(Fig_4a$Genotype)
Fig_4a$Genotype <- factor(Fig_4a$Genotype, levels = c("WT","cHET", "cKO"))

## Figure 4a Aurka ####
Fig_4a_Aurka <- Fig_4a %>% select(Genotype, Aurka)

## Kruskal–Wallis test ####
# For bounded data, we proceeded with non-parametric Kruskal-Wallis H test to compare the groups.
KWH = kruskal.test(Aurka ~ Genotype, data = Fig_4a)

# posthoc test (Dunn's Test) differences between genotype groups
dunn_test_result <- Fig_4a %>% 
  dunn_test(Aurka ~ Genotype, p.adjust.method = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig4a_Statistics.txt"))
cat("Fig4a Aurka Statistics\n")
cat("\nKruskal-Wallis H test:\n")
cat("\n Bounded data, use non-parametric Kruskal-Wallis H test .\n")
print(KWH)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (Dunn's Test with Bonferroni adjustment):\n")
print(dunn_test_result)
cat("\n===============================================n")
sink()

## Figure 4a Bhlhe41 ####
Fig_4a_Bhlhe41 <- Fig_4a %>% select(Genotype, Bhlhe41)
## Kruskal–Wallis test ####
# For bounded data, we proceeded with non-parametric Kruskal-Wallis H test to compare the groups.
KWH = kruskal.test(Bhlhe41 ~ Genotype, data = Fig_4a_Bhlhe41)

# posthoc test (Dunn's Test) differences between genotype groups
dunn_test_result <- Fig_4a_Bhlhe41 %>% 
  dunn_test(Bhlhe41 ~ Genotype, p.adjust.method = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig4a_Statistics.txt"),append = T)
cat("Fig4a Bhlhe41 Statistics\n")
cat("\nKruskal-Wallis H test:\n")
cat("\n Bounded data, use non-parametric Kruskal-Wallis H test .\n")
print(KWH)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (Dunn's Test with Bonferroni adjustment):\n")
print(dunn_test_result)
sink()

#################################################################
# Fig4b ####
Fig_4b <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_4b")
head(Fig_4b)
table(Fig_4b$Genotype, Fig_4b$CNS_Region)
Fig_4b$Genotype <- factor(Fig_4b$Genotype, levels = c("WT","cKO"))
Fig_4b$CNS_Region <- factor(Fig_4b$CNS_Region, levels = c('Cortex', 'CC', 'CA1', 'CA3' ,  'DG'))

## Generalized Linear Model (GLM)  ######
model_poisson <- glm(Microglia ~ CNS_Region * Genotype, data = Fig_4b, family = poisson(link = "log"))
anova(model_poisson, test="Chisq")
# Check overdispersion
overdispersion <- sum(residuals(model_poisson, type = "pearson")^2) / model_poisson$df.residual

# posthoc test differences between genotype groups
posthoc <- emmeans(model_poisson, pairwise ~ Genotype | CNS_Region, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig4b_Statistics.txt"))
cat("Fig4b Statistics\n")
cat("\nData type: count\n")
cat("\nGeneralized Linear Model (GLM) Results:\n")
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


#################################################################
# Fig4d ####
Fig_4d <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_4d")
head(Fig_4d)
table(Fig_4d$Genotype)
Fig_4d$Genotype <- factor(Fig_4d$Genotype, levels = c("WT","cKO"))

## Check data distribution####
by(Fig_4d$AF488_Pos, Fig_4d$Genotype, shapiro.test)

pdf("output/Diagnosis/Fig4d_QQ_Plot.pdf", width = 6, height = 6) # 保存为PDF， 6x6英寸
qqnorm(Fig_4d$AF488_Pos, main = "Q-Q Plot of AF488_Pos") # 可添加标题
qqline(Fig_4d$AF488_Pos, lwd = 2) # 可自定义参考线颜色和粗细
dev.off()

# Normality of Residuals
Homogeneity = leveneTest(AF488_Pos ~ Genotype, data = Fig_4d)

## t test ####
if (Homogeneity$`Pr(>F)`[1] > 0.05) {
  t_test_result <- t.test(AF488_Pos ~ Genotype, data = Fig_4d, var.equal = TRUE)
} else {
  t_test_result <- t.test(AF488_Pos ~ Genotype, data = Fig_4d, var.equal = FALSE)
}

sample_estimates <- t_test_result$estimate
conf_int_lower <- t_test_result$conf.int[1]
conf_int_upper <- t_test_result$conf.int[2]

result_df <- data.frame(
  Variable = deparse(substitute(AF488_Pos)), 
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
sink(file = file.path(output, "Fig4d_Statistics.txt"))
cat("Fig4d AF488_Pos Statistics\n")
cat("\n--------------------\n")
cat("\nCheck assumptions:\n")
cat("\nNormal distribution:\n")
print(by(Fig_4d$AF488_Pos, Fig_4d$Genotype, shapiro.test))
cat("\nHomogeneity of Variance (Levene's Test):\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nt-test Results:\n")
print(result_df)
cat("\n=====================\n")
sink()

################################################
# Fig4f ####
Fig4f <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_4f")
head(Fig4f)
table(Fig4f$Genotype, Fig4f$Model_Phase)

Fig4f$Genotype <- factor(Fig4f$Genotype, levels = c("WT","cKO"))
Fig4f$Model_Phase <- factor(Fig4f$Model_Phase, levels = c("Control" , "CPZ", "Recovery"))

## LMM Analysis ####
# test difference between Genotype in each model phase and Axon_Status
model <- lmer(Axons ~ Genotype * Model_Phase * Axon_Status + (1|Mouse_ID), data = Fig4f)
summary(model)
Anova(model, type = 3)

# Check model assumptions:
library(car)
sim_res <- simulateResiduals(model)
test_result <- testUniformity(sim_res)# Test uniformity
print(test_result)  # Check p-value
pdf(file = file.path(output, "Diagnosis/Fig4f_LMM_Residuals_QQPlot.pdf"), width = 8, height = 5)
plot(sim_res)
dev.off()
Homogeneity = leveneTest(residuals(model) ~ Genotype, data = Fig4f)# Test if residuals differ by genotype group -use Levene's test # Homogeneity of Variance
dispersion = testDispersion(sim_res)$p.value 

# If homogeneity of variance is violated, fit a model with variance structure
if (Homogeneity$`Pr(>F)`[1] < 0.05){
  # heteroscedastic → nlme
  print("Fitting LMM with variance structure due to heteroscedasticity: Weighted Linear Mixed-Effects Model)")
  model <- nlme::lme(Axons ~ Genotype * Model_Phase * Axon_Status, 
                     random = ~ 1 | Mouse_ID, 
                     weights = varIdent(form = ~ 1 | Genotype),
                     data = Fig4f)
  Anova(model, type = 3)
} 

# posthoc test differences between genotype groups in each model phase and Axon_Status
emms_results <- emmeans(model, ~ Genotype | Model_Phase * Axon_Status)
posthoc <- pairs(emms_results)

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig4f_Statistics.txt"))
cat("Fig4f Statistics\n")
cat("\nLinear Mixed Model Results:\n")
print(model$call)
print(Anova(model, type = 3))
cat("\n--------------------\n")
cat("\nCheck assumptions:\n")
cat("\nTest uniformity (DHARMa):\n")
print(test_result)
cat("\nTest of overdispersion:\n")
print(dispersion)
cat("\nTest of Homogeneity:\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (EMMeans with Bonferroni adjustment):\n")
print(posthoc)
cat("\n=====================\n")
sink()


####################################################
# Fig_4eghd ####
Fig_4eghd <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_4eghd")
head(Fig_4eghd)
table(Fig_4eghd$Genotype, Fig_4eghd$Model_phase)

Fig_4eghd$Genotype <- factor(Fig_4eghd$Genotype, levels = c("WT","cKO"))
Fig_4eghd$Model_phase <- factor(Fig_4eghd$Model_phase, levels = c("Control" , "CPZ", "Recovery"))

## Figure 4e ####
Fig4e <- Fig_4eghd %>% select(Model_phase, Genotype, `4e_BGII`) 
## Check data distribution####
Fig4e$Group_Combined <- interaction(Fig4e$Model_phase, Fig4e$Genotype, sep="_")
by(Fig4e$`4e_BGII`, Fig4e$Group_Combined, shapiro.test)
## Two way anova ####
model <- aov(`4e_BGII` ~ Model_phase * Genotype, data = Fig4e)
anova_result <- summary(model)
# Normality of Residuals
Normality = shapiro.test(residuals(model))
pdf(file = file.path(output, "Diagnosis/Fig4e_Residuals_QQPlot.pdf"))
plot(model, 2)
dev.off()
# Check homogeneity of variances
Homogeneity <- car::leveneTest(`4e_BGII` ~ Model_phase * Genotype, data = Fig4e)
# posthoc test differences between genotype groups
posthoc <- emmeans(model, pairwise ~ Genotype | Model_phase, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig4e_Statistics.txt"))
cat("Fig4e Statistics\n")
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
cat("\n=====================\n")
sink()

## Figure 4g ####
Fig4g <- Fig_4eghd %>% select(Model_phase, Genotype, `4g_CC1`)
## Check data distribution####
Fig4g$Group_Combined <- interaction(Fig4g$Model_phase, Fig4g$Genotype, sep="_")
by(Fig4g$`4g_CC1`, Fig4g$Group_Combined, shapiro.test)

## Generalized Linear Model (GLM) with Poisson or negative binomial ######
mode <- glm(`4g_CC1` ~ Model_phase * Genotype, data = Fig4g, family = poisson(link = "log"))
anova(mode, test="Chisq")

# Check overdispersion
overdispersion <- sum(residuals(mode, type = "pearson")^2) / mode$df.residual

# posthoc test differences between genotype groups
posthoc <- emmeans(mode, pairwise ~ Genotype | Model_phase, adjust = "bonferroni")
# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig4g_Statistics.txt"))
cat("Fig4g Statistics\n")
cat("\nGeneralized Linear Model (poisson) Results:\n")
print(summary(mode))
print(anova(model, test="Chisq"))
cat("\n--------------------\n")
cat("\nOverdispersion Check:\n")
cat("Overdispersion Ratio:", overdispersion, "\n")
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (EMMeans with Bonferroni adjustment):\n")
print(posthoc$contrasts)
cat("\n=====================\n")
sink()


## Figure 4h_dMbp ####
Fig4h <- Fig_4eghd %>% select(Model_phase, Genotype, `4h_dMbp`)
colnames(Fig4h) <- c("Model_phase", "Genotype", "dMbp")
## Check data distribution####
Fig4h$Group_Combined <- interaction(Fig4h$Model_phase, Fig4h$Genotype, sep="_")
by(Fig4h$`dMbp`, Fig4h$Group_Combined, shapiro.test)
## Two way anova ####
model <- aov(`dMbp` ~ Model_phase * Genotype, data = Fig4h)
anova_result <- summary(model)
# Normality of Residuals
Normality = shapiro.test(residuals(model))
pdf(file = file.path(output, "Diagnosis/Fig4h_Residuals_QQPlot.pdf"))
plot(model, 2)
dev.off()

# Check homogeneity of variances
Homogeneity <- car::leveneTest(`dMbp` ~ Model_phase * Genotype, data = Fig4h)

## Switch to ART ANOVA ####
# The QQ plot revealed a significant departure from normality. Therefore, we proceeded with Aligned Rank Transform (ART) ANOVA.
model_art <- art(`dMbp` ~  Genotype*Model_phase, data = Fig4h)
anova(model_art)

# posthoc test differences between genotype groups
posthoc_results <- art.con(model_art, "Genotype:Model_phase", adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig4h_Statistics.txt"))
cat("Fig4h dMbp Statistics\n")
cat("\nART ANOVA Results:\n")
print(anova(model_art))
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (ART Contrasts with Bonferroni adjustment):\n")
print(posthoc_results)
cat("\n=====================\n")
sink()

## Figure 4h_dMbp_Iba1 ####
Fig_4h_dMbp_Iba1 <- Fig_4eghd %>% select(Genotype, Model_phase, `4h_dMbp_Iba1`)
colnames(Fig_4h_dMbp_Iba1) <- c("Genotype","Model_phase", "dMbp_Iba1")
## Check data distribution####
Fig_4h_dMbp_Iba1$Group_Combined <- interaction(Fig_4h_dMbp_Iba1$Model_phase, Fig_4h_dMbp_Iba1$Genotype, sep="_")
by(Fig_4h_dMbp_Iba1$`dMbp_Iba1`, Fig_4h_dMbp_Iba1$Group_Combined, shapiro.test)

## Two way anova ####
model <- aov(`dMbp_Iba1` ~ Model_phase * Genotype, data = Fig_4h_dMbp_Iba1)
anova_result <- summary(model)
# Normality of Residuals
Normality = shapiro.test(residuals(model))
pdf(file = file.path(output, "Diagnosis/Fig4h_dMbp_Iba1_Residuals_QQPlot.pdf"))
plot(model, 2)
dev.off()
# Check homogeneity of variances
Homogeneity <- car::leveneTest(`dMbp_Iba1` ~ Model_phase * Genotype, data = Fig_4h_dMbp_Iba1)

##switch to ART ANOVA ####
# The QQ plot revealed a significant departure from normality. Therefore, we proceeded with Aligned Rank Transform (ART) ANOVA.
model_art <- art(`dMbp_Iba1` ~  Genotype*Model_phase , data = Fig_4h_dMbp_Iba1)
anova(model_art)

# posthoc test differences between genotype groups
posthoc_results <- art.con(model_art, "Genotype:Model_phase", adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig4h_Statistics.txt"),append=T)
cat("Fig4h dMbp_Iba1 Statistics\n")
cat("\nART ANOVA Results:\n")
print(anova(model_art))
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (ART Contrasts with Bonferroni adjustment):\n")
print(posthoc_results)
cat("\n=====================\n")
sink()
