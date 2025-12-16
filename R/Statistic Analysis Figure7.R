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
# Fig_7b####
Fig_7b <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_7b")
head(Fig_7b)
table(Fig_7b$Treatment)
Fig_7b$Treatment <- factor(Fig_7b$Treatment, levels = c("Anti-CD22", "Isotype"))

## Check data distribution####
by(Fig_7b$Mye_Dil, Fig_7b$Treatment, shapiro.test)
pdf("output/Diagnosis/Fig7b_QQ_Plot.pdf", width = 6, height = 6)
qqnorm(Fig_7b$Mye_Dil, main = "Q-Q Plot of Mye_Dil")
qqline(Fig_7b$Mye_Dil, lwd = 2)
dev.off()

# Normality of Residuals
Homogeneity = leveneTest(Mye_Dil ~ Treatment, data = Fig_7b)

## Paired t test ####
if (Homogeneity$`Pr(>F)`[1] > 0.05) {
  t_test_result <- t.test(Mye_Dil ~ Treatment, data = Fig_7b, paired = TRUE, var.equal = TRUE)
} else {
  t_test_result <- t.test(Mye_Dil ~ Treatment, data = Fig_7b, paired = TRUE, var.equal = FALSE)
}

sample_estimates <- t_test_result$estimate
conf_int_lower <- t_test_result$conf.int[1]
conf_int_upper <- t_test_result$conf.int[2]

result_df <- data.frame(
  Variable = deparse(substitute(Mye_Dil)), 
  Test = ifelse(Homogeneity$`Pr(>F)`[1] > 0.05, "Paired Student's t-test", "Paired Welch's t-test"),
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
sink(file = file.path(output, "Fig7b_Statistics.txt"))
cat("Fi7b Mye_Dil Statistics\n")
cat("\n--------Paired t test------------\n")
cat("\nCheck assumptions:\n")
cat("\nNormal distribution:\n")
print(by(Fig_7b$Mye_Dil, Fig_7b$Treatment, shapiro.test))
cat("\nHomogeneity of Variance (Levene's Test):\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nPaired t-test Results:\n")
print(result_df)
cat("\n=====================\n")
sink()
###################################
# Fig_7c####
Fig_7c <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_7c")
head(Fig_7c)
table(Fig_7c$Treatment)
Fig_7c$Treatment <- factor(Fig_7c$Treatment, levels = c("Anti-CD22", "Isotype"))

## Check data distribution####
by(Fig_7c$Mye_Dil, Fig_7c$Treatment, shapiro.test)
pdf("output/Diagnosis/Fig7c_QQ_Plot.pdf", width = 6, height = 6)
qqnorm(Fig_7c$Mye_Dil, main = "Q-Q Plot of Mye_Dil")
qqline(Fig_7c$Mye_Dil, lwd = 2)
dev.off()
# Residuals
Homogeneity = leveneTest(Mye_Dil ~ Treatment, data = Fig_7c)
## Paired t test ####
if (Homogeneity$`Pr(>F)`[1] > 0.05) {
  t_test_result <- t.test(Mye_Dil ~ Treatment, data = Fig_7c, paired = TRUE, var.equal = TRUE)
} else {
  t_test_result <- t.test(Mye_Dil ~ Treatment, data = Fig_7c, paired = TRUE, var.equal = FALSE)
}
sample_estimates <- t_test_result$estimate
conf_int_lower <- t_test_result$conf.int[1]
conf_int_upper <- t_test_result$conf.int[2]
result_df <- data.frame(
  Variable = deparse(substitute(Mye_Dil)), 
  Test = ifelse(Homogeneity$`Pr(>F)`[1] > 0.05, "Paired Student's t-test", "Paired Welch's t-test"),
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
sink(file = file.path(output, "Fig7c_Statistics.txt"))
cat("Fi7c Mye_Dil Statistics\n")
cat("\n--------Paired t test------------\n")
cat("\nCheck assumptions:\n")
cat("\nNormal distribution:\n")
print(by(Fig_7c$Mye_Dil, Fig_7c$Treatment, shapiro.test))
cat("\n--------------------\n")
cat("\nHomogeneity of Variance (Levene's Test):\n")
print(Homogeneity)
cat("\n--------------------\n")
cat("\nPaired t-test Results:\n")
print(result_df)
cat("\n=====================\n")
sink()

###################################
# Fig_7efghi####
Fig_7efghi <- readxl::read_excel("Source Data.xlsx", sheet = "Fig_7efghi")
head(Fig_7efghi)
table(Fig_7efghi$Treatment, Fig_7efghi$Model_Phase)

Fig_7efghi$Treatment <- factor(Fig_7efghi$Treatment, levels = c("Anti_CD22", "Isotype"))
Fig_7efghi$Model_Phase <- factor(Fig_7efghi$Model_Phase, levels = c("CPZ", "Recovery"))
Fig_7efghi$Group_combined <- interaction(Fig_7efghi$Treatment, Fig_7efghi$Model_Phase, sep="_")
colnames(Fig_7efghi) <-c("Model_Phase","Treatment","BGII","Mbp","Olig2","CC1","dMbp","dMbp_Iba1"  ,"Group_combined")

# Fig7e ####
Fig7e <- Fig_7efghi %>% select(Model_Phase, Treatment,BGII,Group_combined,BGII)
head(Fig7e)

## Two-way ANOVA ####
model <- aov(BGII ~ Model_Phase *Treatment, data = Fig7e)
anova_results <- summary(model)
# Check assumptions
# Normality of residuals
shapiro_test <- shapiro.test(residuals(model))
pdf("output/Diagnosis/Fig7e_Residuals_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()

# Homogeneity of variances
levene_test <- leveneTest(`BGII` ~ Group_combined, data = Fig7e)

# Post-hoc test  
posthoc <- emmeans(model, pairwise ~ Treatment*Model_Phase, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig7e_Statistics.txt"),append = F) #append = TRUE
cat("Fig7e Statistics\n")
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

#Fig7f####
Fig7f <- Fig_7efghi %>% select(Model_Phase, Treatment,BGII,Group_combined,Mbp)
head(Fig7f)

## Two-way ANOVA ####
model <- aov(Mbp ~ Model_Phase *Treatment, data = Fig7f)
anova_results <- summary(model)
# Check assumptions
# Normality of residuals
shapiro_test <- shapiro.test(residuals(model))
pdf("output/Diagnosis/Fig7f_Residuals_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()

# Homogeneity of variances
levene_test <- leveneTest(Mbp ~ Group_combined, data = Fig7f)

# Post-hoc test  
posthoc <- emmeans(model, pairwise ~ Treatment*Model_Phase, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig7f_Statistics.txt"),append = F) #append = TRUE
cat("Fig7f Statistics\n")
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

#Fig7g ####
Fig7g <- Fig_7efghi %>% select(Model_Phase, Treatment,BGII,Group_combined,Olig2)
head(Fig7g)

## Generalized Linear Model (GLM) with negative binomial ######
model_nb <- glm.nb(Olig2 ~ Model_Phase * Treatment, data = Fig7g)
anova(model_nb, test="Chisq")

# Check overdispersion
overdispersion <- sum(residuals(model_nb, type = "pearson")^2) / model_nb$df.residual

# posthoc test differences between genotype groups
posthoc_nb <- emmeans(model_nb, pairwise ~ Treatment * Model_Phase, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig7g_Statistics.txt"))
cat("Fig7g Statistics\n")
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

#Fig7h ####
Fig7h <- Fig_7efghi %>% select(Model_Phase, Treatment,BGII,Group_combined,CC1)
head(Fig7h)

## Generalized Linear Model (GLM) poisson ######
model <- glm(CC1 ~ Model_Phase * Treatment, data = Fig7h, family = poisson(link = "log"))
anova(model, test="Chisq")

# Check overdispersion
overdispersion <- sum(residuals(model, type = "pearson")^2) / model$df.residual

# posthoc test differences between genotype groups
posthoc <- emmeans(model, pairwise ~ Treatment* Model_Phase, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig7h_Statistics.txt"))
cat("Fig7h Statistics\n")
cat("\nData type: count\n")
cat("\nGeneralized Linear Model (GLM) poisson:\n")
print(summary(model))
print(anova(model, test="Chisq"))
cat("\nOverdispersion parameter:\n")
print(overdispersion)
cat("\n--------------------\n")
cat("\nPost-hoc Test Results (EMMeans with Bonferroni adjustment):\n")
print(posthoc$contrasts)
sink()

# Fig7i ####
## 7i_dMbp ####
Fig7i <- Fig_7efghi %>% select(Model_Phase, Treatment,Group_combined,dMbp)
head(Fig7i)

## Two-way ANOVA ####
model <- aov(dMbp ~ Model_Phase *Treatment, data = Fig7i)
anova_results <- summary(model)
# Check assumptions
# Normality of residuals
shapiro_test <- shapiro.test(residuals(model))
pdf("output/Diagnosis/Fig7i_dMbp_Residuals_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()

# Homogeneity of variances
levene_test <- leveneTest(dMbp ~ Group_combined, data = Fig7i)

# Post-hoc test  
posthoc <- emmeans(model, pairwise ~ Treatment*Model_Phase, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig7i_Statistics.txt"),append = F) #append = TRUE
cat("Fig7i dMbp Statistics\n")
cat("\nTwo-way ANOVA Results:\n")
model
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

## 7i_dMbp_Iba1 ####
Fig7i <- Fig_7efghi %>% select(Model_Phase, Treatment,Group_combined,dMbp,dMbp_Iba1)
head(Fig7i)

## Two-way ANOVA ####
model <- aov(dMbp_Iba1 ~ Model_Phase *Treatment, data = Fig7i)
anova_results <- summary(model)
# Check assumptions
# Normality of residuals
shapiro_test <- shapiro.test(residuals(model))
pdf("output/Diagnosis/Fig7i_dMbp_Iba1_Residuals_QQ_Plot.pdf", width = 6, height = 6)
car::qqPlot(residuals(model), main = "Q-Q Plot")
dev.off()

# Homogeneity of variances
levene_test <- leveneTest(dMbp ~ Group_combined, data = Fig7i)

# Post-hoc test  
posthoc <- emmeans(model, pairwise ~ Treatment*Model_Phase, adjust = "bonferroni")

# Report results including assumption check results, model summary and post hoc results
sink(file = file.path(output, "Fig7i_Statistics.txt"),append = T)  
cat("Fig7i dMbp_Iba1 Statistics\n")
cat("\nTwo-way ANOVA Results:\n")
model
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
