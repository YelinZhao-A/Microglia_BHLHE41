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

rm(list=ls())

output = "output"
if (!dir.exists(output)) dir.create(output)

# Define function ####
analyze_var <- function(data_df, var_name, figure_name,
                        count_vars = c("Number_of_Branches","Total_Intersections","Max_Intersection")) {
  
  model_type <- NA
  overdisp_flag <- FALSE
  hetero_flag <- FALSE
  diagnostics <- list()
  
  # -------------------------
  # 1. Fit model
  # -------------------------
  if (var_name %in% count_vars) {
    # Count variable  GLMM / GLM / NB
    try_model <- try(glmer(
      as.formula(paste(var_name,"~ Genotype + (1 | Mouse_ID)")),
      data = data_df,
      family = poisson(link="log"),
      control = glmerControl(optimizer="bobyqa")
    ), silent = TRUE)
    
    if (class(try_model) == "try-error" || isSingular(try_model)) {
      model <- glm(as.formula(paste(var_name,"~ Genotype")), data = data_df, family=poisson)
      model_type <- "GLM"
      
      # Overdispersion check
      rp <- residuals(model, type="pearson")
      dispersion <- sum(rp^2)/df.residual(model)
      diagnostics$Overdispersion <- dispersion
      if (dispersion > 1.2) {
        model <- MASS::glm.nb(as.formula(paste(var_name,"~ Genotype")), data = data_df)
        model_type <- "NB-GLM"
        overdisp_flag <- TRUE
      }
      
    } else {
      # GLMM fitted → DHARMa check
      sim <- DHARMa::simulateResiduals(try_model)
      disp_test <- DHARMa::testDispersion(sim)
      ks_test <- DHARMa::testUniformity(sim) 
      diagnostics$DHARMa_testDispersion <- disp_test
      diagnostics$DHARMa_KS_pvalue <- ks_test$p.value
      
      if (disp_test$p.value < 0.05 || ks_test$p.value < 0.05) {
        # Overdispersion detected OR Distribution uniformity violated (KS failed)
        message("Distribution check failed (Dispersion or KS test < 0.05) → Switching to Negative Binomial GLMM")
        model <- glmer.nb(as.formula(paste(var_name,"~ Genotype + (1 | Mouse_ID)")),
                          data = data_df,
                          control=glmerControl(optimizer="bobyqa"))
        model_type <- "NB-GLMM"
        overdisp_flag <- TRUE
      } else {
        model <- try_model
        model_type <- "GLMM_Poisson"
      }
    }
    
  } else {
    # Continuous variable → LMM
    model <- lmer(as.formula(paste(var_name,"~ Genotype + (1 | Mouse_ID)")), data = data_df)
    model_type <- "LMM"
    
    # Diagnostics: Residuals vs Fitted variance check  
    fit_vals <- fitted(model)
    res_vals <- resid(model)
    bp <- bptest(resid(model) ~ fitted(model))
    Homogeneity = leveneTest(residuals(model) ~ Genotype, data = data_df)
    
    # 1. remove NA
    valid_idx <- which(!is.na(fit_vals) & !is.na(res_vals))
    fit_vals <- fit_vals[valid_idx]
    res_vals <- res_vals[valid_idx]
    
    if (length(valid_idx) < 3) {
      diagnostics$Var_by_fit_group <- "Data insufficient (<3 points)"
    } else {
      q <- quantile(fit_vals, probs=c(0,0.33,0.66,1))
      
      subset_var <- function(f, r, q_min, q_max) {
        subset_r <- r[f > q_min & f <= q_max]
        if (length(subset_r) < 2) return(NA_real_)
        return(var(subset_r))
      }
      
      var1 <- subset_var(fit_vals, res_vals, q[1], q[2])
      var2 <- subset_var(fit_vals, res_vals, q[2], q[3])
      var3 <- subset_var(fit_vals, res_vals, q[3], q[4])
      
      all_vars <- c(var1, var2, var3)
      diagnostics$Var_by_fit_group <- all_vars
      
      # Check heteroscedasticity
      if (Homogeneity$`Pr(>F)`[1] < 0.05) {
        # heteroscedastic → nlme
        model <- nlme::lme(as.formula(paste(var_name,"~ Genotype")),
                           random = ~1|Mouse_ID,
                           data = data_df,
                           weights = varIdent(form = ~ 1 | Genotype)) 
        model_type <- "nlme-LMM"
        hetero_flag <- TRUE
      }
    }
  }
  
  # -------------------------
  # Save residual diagnostics plots
  # -------------------------
  diag_prefix <- paste0(output,"/",figure_name,"_",var_name)
  if (model_type %in% c("LMM","nlme-LMM")) {
    png(paste0(diag_prefix,"_QQ.png"),8,5,res=300,units="in")
    qqnorm(resid(model)); qqline(resid(model))
    dev.off()
    
    png(paste0(diag_prefix,"_ResidualFit.png"),8,5,res=300,units="in")
    plot(fitted(model), resid(model)); abline(h=0,col="red")
    dev.off()
  } else {
    sim <- DHARMa::simulateResiduals(model)
    png(paste0(diag_prefix,"_DHARMa.png"),8,5,res=300,units="in")
    plot(sim)
    dev.off()
  }
  
  # -------------------------
  # Fixed effects table
  # -------------------------
  if (model_type == "nlme-LMM") {
    raw <- summary(model)$tTable[-1, , drop=FALSE]
    fixed_df <- data.frame(
      Comparison = rownames(raw),
      Figure = figure_name,
      Variable = var_name,
      ModelType = model_type,
      Estimate = raw[, "Value"],
      Std_Error = raw[, "Std.Error"],
      DF = raw[, "DF"],
      Statistic = raw[, "t-value"],
      p_value = raw[, "p-value"]
    )
  } else {
    raw <- as.data.frame(coef(summary(model)))
    raw <- raw[rownames(raw)!="(Intercept)", , drop=FALSE]
    
    Estimate <- if("Estimate" %in% colnames(raw)) raw$Estimate else raw$Value
    Std_Error <- if(any(grepl("Std", colnames(raw)))) raw[, grep("Std", colnames(raw), value=TRUE)] else NA
    DF <- if("df" %in% colnames(raw)) raw$df else NA
    Statistic <- if("t value" %in% colnames(raw)) raw$`t value` else if("z value" %in% colnames(raw)) raw$`z value` else NA
    p_value <- if("Pr(>|t|)" %in% colnames(raw)) raw$`Pr(>|t|)` else if("Pr(>|z|)" %in% colnames(raw)) raw$`Pr(>|z|)` else NA
    
    fixed_df <- data.frame(
      Comparison = rownames(raw),
      Figure = figure_name,
      Variable = var_name,
      ModelType = model_type,
      Estimate = Estimate,
      Std_Error = Std_Error,
      DF = DF,
      Statistic = Statistic,
      p_value = p_value
    )
  }
  fixed_df$FDR <- p.adjust(fixed_df$p_value, method="fdr")
  
  # -------------------------
  # posthoc
  # -------------------------
  emm <- emmeans(model, ~ Genotype)
  posthoc <- summary(pairs(emm, adjust="bonferroni"))
  stat_col <- intersect(colnames(posthoc), c("t.ratio","z.ratio"))[1]
  
  posthoc_df <- data.frame(
    Figure = figure_name,
    Variable = var_name,
    ModelType = model_type,
    Comparison = posthoc$contrast,
    Estimate = posthoc$estimate,
    Std_Error = posthoc$SE,
    DF = posthoc$df,
    Statistic = posthoc[[stat_col]],
    p_value_tukey = posthoc$p.value
  )
  
  # -------------------------
  # 6. Return diagnostics
  # -------------------------
  diag_df <- data.frame(
    Figure = figure_name,
    Variable = var_name,
    ModelType = model_type,
    Overdispersion = overdisp_flag,
    Heteroscedasticity = hetero_flag,
    Diagnostics = I(list(diagnostics))
  )
  
  return(list(Fixed=fixed_df, Posthoc=posthoc_df, Diagnostics=diag_df))
}

# Load data####
Fig2B <- readxl::read_excel("Sholl-related data.xlsx", sheet = "Fig2B")
Fig4C <- readxl::read_excel("Sholl-related data.xlsx", sheet = "Fig4C")
Fig6A <- readxl::read_excel("Sholl-related data.xlsx", sheet = "Fig6A")
Fig6B <- readxl::read_excel("Sholl-related data.xlsx", sheet = "Fig6B")
FigS3B <- readxl::read_excel("Sholl-related data.xlsx", sheet = "FigS3B")
FigS3D <- readxl::read_excel("Sholl-related data.xlsx", sheet = "FigS3D")
FigS6D <- readxl::read_excel("Sholl-related data.xlsx", sheet = "FigS6D")

variables_to_analyze <- c("Number_of_Branches", "Total_Branches", "Total_Intersections", "Max_Intersection")
count_vars <- c("Number_of_Branches","Total_Intersections","Max_Intersection")

data_list <- list(
  Figure2B = list(data = Fig2B, levels = c("WT", "HE", "KO")),
  Figure4C = list(data = Fig4C, levels = c("WT", "KO")),
  Figure6A = list(data = Fig6A, levels = c("HE", "B41Tg")),
  Figure6B = list(data = Fig6B, levels = c("KO", "KOB41Tg")),
  FigureS3B = list(data = FigS3B, levels = c("WT", "HE", "KO")),
  FigureS3D = list(data = FigS3D, levels = c("WT", "HE", "KO")),
  FigureS6D = list(data = FigS6D, levels = c("WT", "HE", "KO"))
)

## Preprocess data####
processed_data_list <- list()
for (name in names(data_list)) {
  df <- data_list[[name]]$data
  levels <- data_list[[name]]$levels
  
  processed_df <- df %>%
    dplyr::mutate(
      Genotype = factor(Genotype, levels = levels),
      Mouse_ID = factor(Mouse_ID),
      Total_Branches_Log = log1p(Total_Branches)
    )
  
  variables_current <- variables_to_analyze
  variables_current[variables_current == "Total_Branches"] <- "Total_Branches_Log"
  
  processed_data_list[[name]] <- list(data = processed_df, variables = variables_current)
}

# Run analyses####
Results_All_Fixed <- data.frame()
Results_All_Posthoc <- data.frame()
Results_All_Diagnostics <- data.frame()

for (figure_name in names(processed_data_list)) {
  data_info <- processed_data_list[[figure_name]]
  data_df <- data_info$data
  variables_current <- data_info$variables
  
  message("###################################################")
  message("STARTING ANALYSIS FOR: ", figure_name)
  message("###################################################")
  
  for (var in variables_current) {
    res <- analyze_var(data_df, var, figure_name, count_vars)
    
    
    Results_All_Fixed <- bind_rows(Results_All_Fixed, res$Fixed)
    Results_All_Posthoc <- bind_rows(Results_All_Posthoc, res$Posthoc)
    Results_All_Diagnostics <- bind_rows(Results_All_Diagnostics, res$Diagnostics)
    
    write.csv(Results_All_Fixed, paste0(output, "/All_Figures_Fixed_Effects.csv"), row.names=FALSE)
    write.csv(Results_All_Posthoc, paste0(output, "/All_Figures_Posthoc.csv"), row.names=FALSE)
    write.csv(Results_All_Diagnostics, paste0(output, "/Diagnosis/All_Figures_Diagnostics_Summary.csv"), row.names=FALSE)
  }
}


