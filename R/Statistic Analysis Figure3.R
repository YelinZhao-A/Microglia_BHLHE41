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
############################################################
#Fig3bcefg ####
## Load data####
Fig3bcefg <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_3bcefg")
head(Fig3bcefg)
table(Fig3bcefg$Genotype,Fig3bcefg$Model_Phase)
Fig3bcefg$Genotype <- factor(Fig3bcefg$Genotype, levels = c("WT",  "B41cKO"))
Fig3bcefg$Model_Phase <- factor(Fig3bcefg$Model_Phase, levels = c("Control",  "CPZ","Recovery"))
Fig3bcefg$Group_Combined <- interaction(Fig3bcefg$Model_Phase, Fig3bcefg$Genotype, sep="_")

###Fig3b ####
Fig3b = Fig3bcefg %>% select(Genotype, Model_Phase, `3b_BGII`,Group_Combined)
head(Fig3b)

## Two-way ANOVA ####
model <- aov(`3b_BGII` ~ Genotype * Model_Phase, data = Fig3b)
anova_results <- summary(model)
# Check assumptions
# Normality of residuals
shapiro_test <- shapiro.test(residuals(model))
pdf("output/Diagnosis/Fig3b_Residuals_QQ_Plot.pdf", width = 6, height = 6)
qqnorm(residuals(model), main = "Q-Q Plot of Residuals for Fig3b")
qqline(residuals(model), lwd = 2)
dev.off()

# Homogeneity of variances
levene_test <- leveneTest(`3b_BGII` ~ Group_Combined, data = Fig3b)

# Post-hoc test  
posthoc <- emmeans(model, pairwise ~ Genotype | Model_Phase, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig3b_Statistics.txt"),append = F) #append = TRUE
cat("Fig3b Statistics\n")
cat("\nTwo-way ANOVA Results:\n")
print(anova_results)
cat("\n----------------------\n")
cat("\nShapiro-Wilk Test for Normality of Residuals:\n")
print(shapiro_test)
cat("\nLevene's Test for Homogeneity of Variances:\n")
print(levene_test)
cat("\n----------------------\n")
cat("\nPosthoc Results with bonferroni adjustment:\n")
print(posthoc$contrasts)
cat("\n=====================\n")
sink()


###Fig3c ####
Fig3c = Fig3bcefg %>% select(Genotype, Model_Phase, `3c_Mbp`,Group_Combined)

## Two-way ANOVA ####
model <- aov(`3c_Mbp` ~ Genotype * Model_Phase, data = Fig3c)
anova_results <- summary(model)
# Check assumptions
# Normality of residuals
shapiro_test <- shapiro.test(residuals(model))
pdf("output/Diagnosis/Fig3c_Residuals_QQ_Plot.pdf", width = 6, height = 6)
qqnorm(residuals(model), main = "Q-Q Plot of Residuals for Fig3c")
qqline(residuals(model), lwd = 2)
dev.off()
# Homogeneity of variances
levene_test <- leveneTest(`3c_Mbp` ~ Group_Combined, data = Fig3c)
# Post-hoc test
posthoc <- emmeans(model, pairwise ~ Genotype | Model_Phase, adjust = "bonferroni")
# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig3c_Statistics.txt"),append = F) #append = TRUE
cat("Fig3c Statistics\n")
cat("\nTwo-way ANOVA Results:\n")
print(anova_results)
cat("\n----------------------\n")
cat("\nShapiro-Wilk Test for Normality of Residuals:\n")
print(shapiro_test)
cat("\nLevene's Test for Homogeneity of Variances:\n")
print(levene_test)
cat("\n----------------------\n")
cat("\nPosthoc Results:\n")
print(posthoc$contrasts)
cat("\n=====================\n")
sink()

###Fig3e ####
Fig3e = Fig3bcefg %>% select(Genotype, Model_Phase, `3e_Olig2`,Group_Combined)
head(Fig3e)

## Generalized Linear Model (GLM) with negative binomial ######
model_nb <- glm.nb(`3e_Olig2` ~ Model_Phase * Genotype, data = Fig3e)
anova(model_nb, test="Chisq")

# Check overdispersion
overdispersion <- sum(residuals(model_nb, type = "pearson")^2) / model_nb$df.residual

# posthoc test differences between genotype groups
posthoc_nb <- emmeans(model_nb, pairwise ~ Genotype | Model_Phase, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig3e_Statistics.txt"))
cat("Fig3e Statistics\n")
cat("\nData type: count\n")
cat("\nGeneralized Linear Model (GLM) with negative binomial Results:\n")
print(summary(model_nb))
print(anova(model_nb, test="Chisq"))
cat("\nOverdispersion parameter:\n")
print(overdispersion)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (EMMeans with Bonferroni adjustment):\n")
print(posthoc_nb$contrasts)
sink()

###Fig3f ####
Fig3f = Fig3bcefg %>% select(Genotype, Model_Phase, `3f_CC1`,Group_Combined)
head(Fig3f)

## Generalized Linear Model (GLM) with negative binomial ######
model_nb <- glm.nb(`3f_CC1` ~ Model_Phase * Genotype, data = Fig3f)
anova(model_nb, test="Chisq")
# Check overdispersion
overdispersion <- sum(residuals(model_nb, type = "pearson")^2) / model_nb$df.residual
# posthoc test differences between genotype groups
posthoc_nb <- emmeans(model_nb, pairwise ~ Genotype | Model_Phase, adjust = "bonferroni")
# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig3f_Statistics.txt"))
cat("Fig3f Statistics\n")
cat("\nData type: count\n")
cat("\nGeneralized Linear Model (GLM) with negative binomial Results:\n")
print(summary(model_nb))
print(anova(model_nb, test="Chisq"))
cat("\nOverdispersion parameter:\n")
print(overdispersion)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (EMMeans with Bonferroni adjustment):\n")
print(posthoc_nb$contrasts)
sink()

###Fig3g ####
Fig3g = Fig3bcefg %>% select(Genotype, Model_Phase, `3g_dMbp`,`3g_dMbp_Iba1`,Group_Combined)
colnames(Fig3g)[colnames(Fig3g) == "3g_dMbp"] <- "dMbp"
colnames(Fig3g)[colnames(Fig3g) == "3g_dMbp_Iba1"] <- "dMbp_Iba1"
head(Fig3g)
Fig3g$dMbp_log2 = log2(Fig3g$dMbp)

####3g_dMbp#### 
#### Two-way ANOVA ####
model <- aov(dMbp_log2  ~ Genotype * Model_Phase, data = Fig3g)
anova_results <- summary(model)
# Check assumptions
# Normality of residuals
shapiro_test <- shapiro.test(residuals(model))
pdf("output/Diagnosis/Fig3g_dMbp_Residuals_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(model, main = "Q-Q Plot")
dev.off()
# Homogeneity of variances
levene_test <- leveneTest(dMbp_log2 ~ Group_Combined, data = Fig3g)

# Violet Homogeneity assumption
#### ART ANOVA - Switch to Aligned Rank Transform (ART) ####
model_art <- art(dMbp_log2 ~ Genotype * Model_Phase, data = Fig3g)
anova_art <- anova(model_art)

# Post hoc analysis
emm <- emmeans(artlm(model_art, "Genotype:Model_Phase"), ~ Genotype | Model_Phase)
posthoc_result = pairs(emm, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig3g_Statistics.txt"))
cat("Fig3g_dMbp Statistics\n")
cat("\nData type: continuous\n")
cat("\n Aligned Rank Transform (ART) ANOVA:\n")
print(anova_art)
cat("\nAssumption Checks:\n")
cat("\nShapiro-Wilk Normality Test of Residuals:\n")
print(shapiro_test)
cat("\nLevene's Test for Homogeneity of Variances:\n")
print(levene_test)
cat("\n--------------------\n")
cat("posthoc with bonferroni adjustment")
print(posthoc_result)
cat("\n======================================\n")
sink()

####3g_dMbp_Iba1####
#### Two-way ANOVA ####
model <- aov(dMbp_Iba1 %>% log2 ~ Genotype * Model_Phase, data = Fig3g)
anova_results <- summary(model)
# Check assumptions
# Normality of residuals
shapiro_test <- shapiro.test(residuals(model))
pdf("output/Diagnosis/Fig3g_dMbp_Iba1_Residuals_QQ_Plot.pdf", width = 6, height = 6)
qqnorm(residuals(model), main = "Q-Q Plot of Residuals for Fig3g_dMbp_Iba1")
qqline(residuals(model), lwd = 2)
dev.off()
# Homogeneity of variances
levene_test <- leveneTest(dMbp_Iba1 ~ Group_Combined, data = Fig3g)

####ART ANOVA ####
# Violet Homogeneity assumption. Therefore, we proceeded with Aligned Rank Transform (ART) ANOVA.
model_art <- art(dMbp_Iba1 %>% log2~ Genotype * Model_Phase, data = Fig3g)
anova_art <- anova(model_art)

# posthoc test differences between genotype groups
emm <- emmeans(artlm(model_art, "Genotype:Model_Phase"), ~ Genotype | Model_Phase)
posthoc_result = pairs(emm, adjust = "bonferroni")

# Pairwise comparisons of Genotype within each condition
posthoc_results <- pairs(emm, adjust = "bonferroni")
print(posthoc_results)

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig3g_Statistics.txt"),append=TRUE)
cat("Fig3g_dMbp_Iba1 Statistics\n")
cat("\nData type: continuous\n")
cat("\nAligned Rank Transform (ART) ANOVA:\n")
print(anova_art)
cat("\nAssumption Checks:\n")
cat("\nShapiro-Wilk Normality Test of Residuals:\n")
print(shapiro_test)
cat("\nLevene's Test for Homogeneity of Variances:\n")
print(levene_test)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (ART Contrasts with Bonferroni adjustment):\n")
print(posthoc_results)
cat("\n============================\n")
sink()


###############################
# Fig3d ####
Fig3d <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_3d")
head(Fig3d)
table(Fig3d$Genotype,Fig3d$Model_Phase)
Fig3d$Genotype <- factor(Fig3d$Genotype, levels = c("WT", "cKO"))
Fig3d$Model_Phase <- factor(Fig3d$Model_Phase, levels = c("Control" , "CPZ", "Recovery"))

## LMM Analysis ####
# test difference between Genotype in each model phase and Axon_Status
model <- lmer(Axons ~ Genotype * Model_Phase * Axon_Status + (1|Mouse_ID), data = Fig3d)
summary(model)
Anova(model, type = 3)

# Check model assumptions:
library(car)
sim_res <- simulateResiduals(model)
test_result <- testUniformity(sim_res)# Test uniformity
print(test_result)  # Check p-value
pdf(file = file.path(output, "Diagnosis/Fig3d_LMM_Residuals_QQPlot.pdf"), width = 8, height = 5)
plot(sim_res)
dev.off()
Homogeneity = leveneTest(residuals(model) ~ Genotype, data = Fig3d)# Test if residuals differ by genotype group -use Levene's test # Homogeneity of Variance
dispersion = testDispersion(sim_res)$p.value 

# If homogeneity of variance is violated, fit a model with variance structure
if (Homogeneity$`Pr(>F)`[1] < 0.05){
  # heteroscedastic → nlme
  print("Fitting LMM with variance structure due to heteroscedasticity: Weighted Linear Mixed-Effects Model)")
  model <- nlme::lme(Axons ~ Genotype * Model_Phase * Axon_Status, 
                     random = ~ 1 | Mouse_ID, 
                     weights = varIdent(form = ~ 1 | Genotype),
                     data = Fig3d)
  Anova(model, type = 3)
} 

# posthoc test differences between genotype groups in each model phase and Axon_Status
emms_results <- emmeans(model, ~ Genotype | Model_Phase * Axon_Status)
posthoc <- pairs(emms_results)

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig3d_Statistics.txt"))
cat("Fig3d Statistics\n")
cat("\nLinear Mixed Model Results:\n")
cat("\nFitting LMM with variance structure due to heteroscedasticity: Weighted Linear Mixed-Effects Model):\n")
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
