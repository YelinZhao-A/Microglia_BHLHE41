# Aurka-Bhlhe41 axis prevents premature aging-like microglial dysfunction and promotes remyelination
Aging accelerates CNS remyelination failure and neurodegeneration. Microglia promote remyelination by phagocytosing myelin debris, but this function is impaired by aging-related CD22 upregulation. However, the molecular mechanisms counteracting premature aging-related microglial dysfunction and remyelination impairment remain unclear. Here, we report that Aurka-Bhlhe41 axis prevents premature aging-like microglial dysfunction and promotes remyelination by restraining progressive CD22 upregulation. We identified that microglia-enriched Bhlhe41 was negatively autoregulated and inhibited by Aurka loss. Bhlhe41- or Aurka-deficient young mice exhibited aging-like microglial morphology, phagocytic deficits, progressive CD22 upregulation, and remyelination impairment in cuprizone-induced demyelination model. Conversely, ectopic Bhlhe41 expression induced hypertrophic microglia, and counteracted phagocytic deficits and CD22 upregulation in Aurka-deficient microglia. CD22 blockade restored phagocytic function and remyelination in Bhlhe41-deficient mice. Notably, a conserved pattern of CD22 upregulation was observed in human PCDH9high microglia subsets with BHLHE41 downregulation. These findings offer insights into potential therapeutic strategies to combat aging-related neurodegeneration and CNS functional decline.


This repository contains the code used to generate the figures in the manuscript. The code is written in R and Python. The code is organized into the following directories:

- `R/`: R script and source data used to generate the figures
- `Python/`: Python script used to generate the figures

## Requirements
- R (version 4.2.2)
The following R packages are used in the analyses:

| Package       | Version   | Purpose |
|---------------|-----------|---------|
| ARTool        | 0.11.2   | Aligned Rank Transform ANOVA for factorial non-parametric designs |
| rstatix       | 0.7.2    | Tidy statistical functions and assumption checks |
| car           | 3.1-2    | Companion to Applied Regression, includes Levene’s test |
| carData       | 3.0-5    | Example datasets for `car` package |
| lmtest        | 0.9-40   | Diagnostic testing for linear models |
| nlme          | 3.1-162  | Linear and nonlinear mixed-effects models |
| MASS          | 7.3-60   | GLMs, negative binomial models, and other statistical functions |
| DHARMa        | 0.4.7    | Residual diagnostics for (generalized) mixed models |
| performance   | 0.10.5   | Model performance checks, including R², ICC, and effect sizes |
| emmeans       | 1.11.1   | Estimated marginal means for post hoc comparisons |
| lmerTest      | 3.1-3    | p-values for linear mixed-effects models |
| lme4          | 1.1-33   | Fitting linear and generalized linear mixed-effects models |
| broom         | 1.0.4    | Tidy summaries of model outputs |
| tidyr         | 1.3.0    | Data reshaping and tidying |
| ggplot2       | 3.4.4    | Data visualization |
| readxl        | 1.4.2    | Reading Excel files |
| dplyr         | 1.1.2    | Data manipulation and transformation |
| Matrix        | 1.5-4    | Sparse and dense matrix operations |

- Python (version 3.12)
  Python packages: requests, os, scanpy, numpy, pandas, scipy, seaborn, matplotlib, sklearn, gc

## Analyses Included

1. **Parametric Tests**
   - Student’s t-test, Welch’s t-test, paired t-test
   - ANOVA (one-way, two-way) with interaction terms
   - Post hoc pairwise comparisons using EMMeans with Bonferroni correction  

2. **Non-parametric Tests**
   - Wilcoxon rank-sum and signed-rank tests
   - Kruskal-Wallis H test for non-factorial designs
   - Aligned Rank Transform (ART) ANOVA for factorial non-parametric designs  

3. **Count Data**
   - Poisson and negative binomial GLMs
   - Overdispersion checks using Pearson residuals  

4. **Mixed-Effects Models**
   - Linear mixed-effects models (LMMs) for repeated measurements
   - Random intercepts for individual subjects
   - Residual diagnostics using DHARMa
   - Variance-structured mixed-effects models using nlme for heteroscedasticity  

4. **scRNA analysis**
   - Data pre-processing
   - Clustering/dimension reduction
   - Expression visualization
   - Comparison between clusters

## Usage
1. Clone the repository
```bash
git clone https://github.com/YelinZhao-A/Microglia_BHLHE41
```
3. Run the R script in the `R/` directory to generate the figures and statistic reports
4. Run the Python script in the `Python/` directory to generate the figures for scRNA-seq related figures




