AMR and Salmonella Dynamics Pipeline
================
Enrique J. Delgado Suárez
11 July, 2026

## INTEGRATED PIPELINE

## PROJECT: SALMONELLA IN SURFACE WATERS FROM CENTRAL MEXICO (2019-2023)

------------------------------------------------------------------------

## I. LOAD REQUIRED LIBRARIES, IMPORT AND ADJUST THE EXCEL FILES

``` r
library(vegan)        # Multivariate analyses (Jaccard, PERMANOVA, PERMDISP)
library(ggplot2)      # Produce figures
library(dplyr)        # Structured data manipulation and cleaning
library(tidyr)        # Data reshaping matrices
library(purrr)        # Functional batch processing execution (map_df)
library(ggalluvial)   # Longitudinal MDR genotype-phenotype alluvial tracking
library(vcd)          # Cohen's Kappa with ASE
library(lme4)         # Generalized Linear Mixed-Effects Models (GLMM via glmer)
library(svglite)      # Scaled vector graphics engine (SVG export support)
library(readxl)       # Reading Excel files
```

``` r
AMRg_veganb <- read_excel("AMRg_veganb.xlsx")
data_muestras_cruda <- read_excel("Prevalence.xlsx")
```

### Clean missing fields and explicitly convert character vectors into strict continuous numeric elements to eradicate Excel date-formatting artifacts.

``` r
cleaned_master_dataset <- AMRg_veganb %>%
  mutate(
    # Force continuous variables into numeric layout
    pH          = as.numeric(as.character(pH)),
    Temperature = as.numeric(as.character(Temperature)),
    Turbidity   = as.numeric(as.character(Turbidity)),
    
    # Configure categorical structures into strict factors
    ID          = as.character(ID),
    Year        = as.factor(Year),
    State       = as.factor(State),
    Source      = as.factor(Source),
    Season      = as.factor(Season),
    Site        = as.factor(Site),
    
    # Establish binary multidrug-resistance indices
    MDRg        = as.numeric(as.character(MDRg)),
    MDRf        = as.numeric(as.character(MDRf))
  )
```

------------------------------------------------------------------------

## II. UNIFIED MODEL: GLOBAL SALMONELLA PREVALENCE (GLMM FULL DISCRETE + CONTINUOUS)

### 1. Consolidated risk factor analysis for sample-level prevalence

``` r
base_muestras_anual <- data_muestras_cruda %>%
  mutate(
    Year_Fact   = factor(Year, levels = c("2019", "2020", "2021", "2022", "2023")),
    State       = as.factor(State),
    Source      = as.factor(Source),
    Season      = as.factor(Season),
    Site        = as.factor(Site),      
    pH          = as.numeric(as.character(pH)),
    Temperature = as.numeric(as.character(Temperature)),
    Turbidity   = as.numeric(as.character(Turbidity)),
    Result      = as.numeric(as.character(Result))
  ) %>%
  filter(!is.na(Year_Fact) & !is.na(State) & !is.na(Source) & !is.na(Season) & 
           !is.na(Site) & !is.na(Result) & !is.na(pH) & !is.na(Temperature) & !is.na(Turbidity))

print(paste("Total independent samples entering consolidated GLMM for Prevalence:", nrow(base_muestras_anual)))
```

    ## [1] "Total independent samples entering consolidated GLMM for Prevalence: 899"

``` r
model_annual_prevalence <- glmer(
  Result ~ Year_Fact + State + Source + Season + pH + Temperature + Turbidity + (1 | Site), 
  data = base_muestras_anual, 
  family = binomial(link = "logit")
)

summary_ann_prev <- summary(model_annual_prevalence)$coefficients
ci_ann_prev      <- confint(model_annual_prevalence, method = "Wald")[rownames(summary_ann_prev), ]

table_annual_prevalence <- data.frame(
  Variable   = rownames(summary_ann_prev),
  Odds_Ratio = round(exp(summary_ann_prev[, "Estimate"]), 4),
  CI_Lower   = round(exp(ci_ann_prev[, 1]), 4),
  CI_Upper   = round(exp(ci_ann_prev[, 2]), 4),
  P_Value    = round(summary_ann_prev[, "Pr(>|z|)"], 4)
) %>% mutate(Significance = ifelse(P_Value < 0.05, "SIGNIFICANT 🌟", "Stable"))

print("DEFINITIVE UNIFIED SAMPLE-LEVEL SALMONELLA PREVALENCE DRIVERS")
```

    ## [1] "DEFINITIVE UNIFIED SAMPLE-LEVEL SALMONELLA PREVALENCE DRIVERS"

``` r
print(as.data.frame(table_annual_prevalence))
```

    ##                                  Variable Odds_Ratio CI_Lower  CI_Upper P_Value
    ## (Intercept)                   (Intercept)   345.0358  30.8274 3861.8092  0.0000
    ## Year_Fact2020               Year_Fact2020     0.4078   0.2129    0.7813  0.0068
    ## Year_Fact2021               Year_Fact2021     0.5663   0.3397    0.9439  0.0292
    ## Year_Fact2022               Year_Fact2022     0.2593   0.1666    0.4036  0.0000
    ## Year_Fact2023               Year_Fact2023     0.0822   0.0475    0.1421  0.0000
    ## StateMexico City         StateMexico City     1.2329   0.5208    2.9190  0.6339
    ## StateMorelos                 StateMorelos     1.6301   0.8501    3.1256  0.1413
    ## StateState of Mexico StateState of Mexico     1.7348   0.8179    3.6795  0.1510
    ## StateTlaxcala               StateTlaxcala     3.0015   1.3840    6.5094  0.0054
    ## SourceDam                       SourceDam     0.9551   0.4503    2.0258  0.9048
    ## SourcePond                     SourcePond     0.6606   0.3683    1.1848  0.1643
    ## SourceRiver                   SourceRiver     0.9597   0.5961    1.5450  0.8655
    ## SeasonSpring                 SeasonSpring     6.8759   4.1580   11.3705  0.0000
    ## SeasonSummer                 SeasonSummer     1.2839   0.8504    1.9382  0.2344
    ## SeasonWinter                 SeasonWinter     3.9527   2.4158    6.4673  0.0000
    ## pH                                     pH     0.4638   0.3502    0.6144  0.0000
    ## Temperature                   Temperature     1.0071   0.9668    1.0491  0.7332
    ## Turbidity                       Turbidity     1.0012   0.9985    1.0039  0.3799
    ##                        Significance
    ## (Intercept)          SIGNIFICANT 🌟
    ## Year_Fact2020        SIGNIFICANT 🌟
    ## Year_Fact2021        SIGNIFICANT 🌟
    ## Year_Fact2022        SIGNIFICANT 🌟
    ## Year_Fact2023        SIGNIFICANT 🌟
    ## StateMexico City             Stable
    ## StateMorelos                 Stable
    ## StateState of Mexico         Stable
    ## StateTlaxcala        SIGNIFICANT 🌟
    ## SourceDam                    Stable
    ## SourcePond                   Stable
    ## SourceRiver                  Stable
    ## SeasonSpring         SIGNIFICANT 🌟
    ## SeasonSummer                 Stable
    ## SeasonWinter         SIGNIFICANT 🌟
    ## pH                   SIGNIFICANT 🌟
    ## Temperature                  Stable
    ## Turbidity                    Stable

``` r
write.csv(table_annual_prevalence, "Table_4A_Unified_Global_Prevalence_Drivers.csv", row.names = FALSE)
```

### 2. Export raw results

``` r
print("EXPORTING RAW GLMM CONSOLE SUMMARY FOR REPOSITORIES")
```

    ## [1] "EXPORTING RAW GLMM CONSOLE SUMMARY FOR REPOSITORIES"

``` r
sink("Raw_Output_Table_2A_GLMM_Prevalence.txt")
cat("RAW GLMM CONSOLE SUMMARY: GLOBAL SALMONELLA PREVALENCE\n")
cat("Generated automatically via open-science verification pipelines\n\n")
print(summary(model_annual_prevalence)) # Reports AIC, BIC, random effects, and z values
```

    ## 
    ## Correlation matrix not shown by default, as p = 18 > 12.
    ## Use print(...., correlation=TRUE)  or
    ##     vcov(....)        if you need it

``` r
sink()
print("Raw console summary successfully exported as 'Raw_Output_Table_4A_GLMM_Prevalence.txt'!")
```

    ## [1] "Raw console summary successfully exported as 'Raw_Output_Table_4A_GLMM_Prevalence.txt'!"

------------------------------------------------------------------------

## III. UNIFIED MODEL: MULTIDRUG RESISTANCE (GLMM FULL DISCRETE + CONTINUOUS)

### 1. Consolidated risk factor analysis for isolate-level MDR phenotypes

``` r
base_mdrf_consolidada <- AMRg_veganb %>%
  mutate(
    Year_Fact   = factor(Year, levels = c("2019", "2020", "2021", "2022", "2023")),
    State       = as.factor(State),
    Source      = as.factor(Source),
    Season      = as.factor(Season),
    Site        = as.factor(Site),      
    pH          = as.numeric(as.character(pH)),
    Temperature = as.numeric(as.character(Temperature)),
    Turbidity   = as.numeric(as.character(Turbidity)),
    MDRf        = as.numeric(as.character(MDRf))
  ) %>%
  filter(!is.na(Year_Fact) & !is.na(State) & !is.na(Source) & !is.na(Season) & 
           !is.na(Site) & !is.na(MDRf) & !is.na(pH) & !is.na(Temperature) & !is.na(Turbidity))

print(paste("Total isolates entering consolidated GLMM for MDRf:", nrow(base_mdrf_consolidada)))
```

    ## [1] "Total isolates entering consolidated GLMM for MDRf: 508"

``` r
model_mdrf_unified <- glmer(
  MDRf ~ Year_Fact + State + Source + Season + pH + Temperature + Turbidity + (1 | Site), 
  data = base_mdrf_consolidada, 
  family = binomial(link = "logit")
)

summary_uni_mdrf <- summary(model_mdrf_unified)$coefficients
ci_uni_mdrf      <- confint(model_mdrf_unified, method = "Wald")[rownames(summary_uni_mdrf), ]

table_mdrf_unified_output <- data.frame(
  Variable   = rownames(summary_uni_mdrf),
  Odds_Ratio = round(exp(summary_uni_mdrf[, "Estimate"]), 4),
  CI_Lower   = round(exp(ci_uni_mdrf[, 1]), 4),
  CI_Upper   = round(exp(ci_uni_mdrf[, 2]), 4),
  P_Value    = round(summary_uni_mdrf[, "Pr(>|z|)"], 4)
) %>% mutate(Significance = ifelse(P_Value < 0.05, "SIGNIFICANT 🌟", "Stable"))

print("DEFINITIVE UNIFIED MULTIDRUG RESISTANCE (MDRf) DRIVERS")
```

    ## [1] "DEFINITIVE UNIFIED MULTIDRUG RESISTANCE (MDRf) DRIVERS"

``` r
print(as.data.frame(table_mdrf_unified_output))
```

    ##                                  Variable Odds_Ratio CI_Lower CI_Upper P_Value
    ## (Intercept)                   (Intercept)     1.1403   0.0381  34.1450  0.9397
    ## Year_Fact2020               Year_Fact2020     0.3287   0.1091   0.9900  0.0480
    ## Year_Fact2021               Year_Fact2021     1.0918   0.4822   2.4722  0.8332
    ## Year_Fact2022               Year_Fact2022     0.7361   0.3900   1.3893  0.3444
    ## Year_Fact2023               Year_Fact2023     0.6382   0.2801   1.4541  0.2851
    ## StateMexico City         StateMexico City     0.2993   0.0607   1.4765  0.1385
    ## StateMorelos                 StateMorelos     0.1458   0.0356   0.5974  0.0074
    ## StateState of Mexico StateState of Mexico     0.2662   0.0647   1.0960  0.0668
    ## StateTlaxcala               StateTlaxcala     0.2027   0.0477   0.8611  0.0306
    ## SourceDam                       SourceDam     0.4319   0.0961   1.9411  0.2736
    ## SourcePond                     SourcePond     1.1044   0.4164   2.9289  0.8419
    ## SourceRiver                   SourceRiver     1.8146   0.9787   3.3645  0.0585
    ## SeasonSpring                 SeasonSpring     1.5545   0.7748   3.1187  0.2143
    ## SeasonSummer                 SeasonSummer     0.9163   0.4984   1.6844  0.7783
    ## SeasonWinter                 SeasonWinter     1.8671   0.8142   4.2816  0.1403
    ## pH                                     pH     0.9978   0.6725   1.4805  0.9915
    ## Temperature                   Temperature     0.9977   0.9413   1.0575  0.9394
    ## Turbidity                       Turbidity     1.0028   0.9996   1.0059  0.0825
    ##                        Significance
    ## (Intercept)                  Stable
    ## Year_Fact2020        SIGNIFICANT 🌟
    ## Year_Fact2021                Stable
    ## Year_Fact2022                Stable
    ## Year_Fact2023                Stable
    ## StateMexico City             Stable
    ## StateMorelos         SIGNIFICANT 🌟
    ## StateState of Mexico         Stable
    ## StateTlaxcala        SIGNIFICANT 🌟
    ## SourceDam                    Stable
    ## SourcePond                   Stable
    ## SourceRiver                  Stable
    ## SeasonSpring                 Stable
    ## SeasonSummer                 Stable
    ## SeasonWinter                 Stable
    ## pH                           Stable
    ## Temperature                  Stable
    ## Turbidity                    Stable

``` r
write.csv(table_mdrf_unified_output, "Table_4B_Unified_MDRf_Drivers.csv", row.names = FALSE)
```

### 2. Export results

``` r
print("EXPORTING RAW GLMM CONSOLE SUMMARY FOR REPOSITORIES")
```

    ## [1] "EXPORTING RAW GLMM CONSOLE SUMMARY FOR REPOSITORIES"

``` r
sink("Raw_Output_Table_4B_GLMM_MDRf.txt")
cat("RAW GLMM CONSOLE SUMMARY: GLOBAL MDR PHENOTYPES\n")
cat("Generated automatically via open-science verification pipelines\n\n")
print(summary(model_mdrf_unified)) # Reports AIC, BIC, random effects, and z values
```

    ## 
    ## Correlation matrix not shown by default, as p = 18 > 12.
    ## Use print(summary(model_mdrf_unified), correlation=TRUE)  or
    ##     vcov(summary(model_mdrf_unified))        if you need it

``` r
sink()
print("Raw console summary successfully exported as 'Raw_Output_Table_4B_GLMM_MDRf.txt'!")
```

    ## [1] "Raw console summary successfully exported as 'Raw_Output_Table_4B_GLMM_MDRf.txt'!"

------------------------------------------------------------------------

## IV. DRUG-LEVEL GENOTYPE-PHENOTYPE CONCORDANCE (BATCH INDIVIDUAL KAPPA CALCULATION)

Evaluates each antibiotic against its target genetic group, capturing
ASE vectors inline

``` r
individual_kappa_table <- list(
  c(ant = "AMP", gen = "bla"), c(ant = "AMC", gen = "bla"), c(ant = "CRO", gen = "bla"),
  c(ant = "FEP", gen = "bla"), c(ant = "MEM", gen = "bla"), c(ant = "CIP", gen = "gyrA"),
  c(ant = "AZM", gen = "mph"), c(ant = "TET", gen = "tetA"), c(ant = "CHL", gen = "floR"),
  c(ant = "AMK", gen = "aac"), c(ant = "STR", gen = "aac"), c(ant = "SXT", gen = "sul")
) %>% 
  map_df(function(par) {
    tabla_par  <- table(cleaned_master_dataset[[par["gen"]]], cleaned_master_dataset[[par["ant"]]])
    if(nrow(tabla_par) != 2 || ncol(tabla_par) != 2) return(data.frame(Antibiotic = par["ant"], Associated_Gene = par["gen"], Kappa_or_ASE = NA))
    kappa_obj  <- Kappa(tabla_par)
    # Extracts the unweighted vector directly containing both Coeff and ASE sequentially
    return(data.frame(Antibiotic = par["ant"], Associated_Gene = par["gen"], Kappa_or_ASE = round(as.numeric(kappa_obj$Unweighted), 4)))
  })
```

    ## Warning in data.frame(Antibiotic = par["ant"], Associated_Gene = par["gen"], :
    ## row names were found from a short variable and have been discarded
    ## Warning in data.frame(Antibiotic = par["ant"], Associated_Gene = par["gen"], :
    ## row names were found from a short variable and have been discarded
    ## Warning in data.frame(Antibiotic = par["ant"], Associated_Gene = par["gen"], :
    ## row names were found from a short variable and have been discarded
    ## Warning in data.frame(Antibiotic = par["ant"], Associated_Gene = par["gen"], :
    ## row names were found from a short variable and have been discarded
    ## Warning in data.frame(Antibiotic = par["ant"], Associated_Gene = par["gen"], :
    ## row names were found from a short variable and have been discarded
    ## Warning in data.frame(Antibiotic = par["ant"], Associated_Gene = par["gen"], :
    ## row names were found from a short variable and have been discarded
    ## Warning in data.frame(Antibiotic = par["ant"], Associated_Gene = par["gen"], :
    ## row names were found from a short variable and have been discarded
    ## Warning in data.frame(Antibiotic = par["ant"], Associated_Gene = par["gen"], :
    ## row names were found from a short variable and have been discarded
    ## Warning in data.frame(Antibiotic = par["ant"], Associated_Gene = par["gen"], :
    ## row names were found from a short variable and have been discarded
    ## Warning in data.frame(Antibiotic = par["ant"], Associated_Gene = par["gen"], :
    ## row names were found from a short variable and have been discarded
    ## Warning in data.frame(Antibiotic = par["ant"], Associated_Gene = par["gen"], :
    ## row names were found from a short variable and have been discarded
    ## Warning in data.frame(Antibiotic = par["ant"], Associated_Gene = par["gen"], :
    ## row names were found from a short variable and have been discarded

``` r
print("COHEN'S KAPPA DRUG-LEVEL CONCORDANCE TABLES WITH ASE VALUES ")
```

    ## [1] "COHEN'S KAPPA DRUG-LEVEL CONCORDANCE TABLES WITH ASE VALUES "

``` r
print(as.data.frame(individual_kappa_table))
```

    ##    Antibiotic Associated_Gene Kappa_or_ASE
    ## 1         AMP             bla       0.8755
    ## 2         AMP             bla       0.0279
    ## 3         AMC             bla       0.1082
    ## 4         AMC             bla       0.0420
    ## 5         CRO             bla       0.1515
    ## 6         CRO             bla       0.0472
    ## 7         FEP             bla       0.0142
    ## 8         FEP             bla       0.0181
    ## 9         MEM             bla      -0.0111
    ## 10        MEM             bla       0.0063
    ## 11        CIP            gyrA      -0.0182
    ## 12        CIP            gyrA       0.0364
    ## 13        AZM             mph       0.5395
    ## 14        AZM             mph       0.0615
    ## 15        TET            tetA       0.9216
    ## 16        TET            tetA       0.0176
    ## 17        CHL            floR       0.9216
    ## 18        CHL            floR       0.0199
    ## 19        AMK             aac       0.0738
    ## 20        AMK             aac       0.0263
    ## 21        STR             aac       0.7225
    ## 22        STR             aac       0.0325
    ## 23        SXT             sul       0.6751
    ## 24        SXT             sul       0.0352

``` r
write.csv(individual_kappa_table, "OR5_Kappa_Individual_ATB.csv", row.names = FALSE)
```

------------------------------------------------------------------------

## V. BATCH CONTINUOUS GLM: LONGITUDINAL TRENDS FOR ANTIBIOTICS/MDR PHENOTYPES

### 1. Enforce the year vector as a continuous numerical covariate at the isolate level

``` r
base_continuous_trends <- AMRg_veganb %>%
  mutate(
    Year_Cont = as.numeric(as.character(Year)),
    MDRf      = as.numeric(as.character(MDRf))
  )
```

### 2. Define the array of the 12 phenotypic antibiotics plus MDR phenotypes (MDRf)

``` r
target_variables <- c("AMP", "AMC", "CRO", "FEP", "MEM", "CIP", "AZM", "TET", "CHL", "AMK", "STR", "SXT", "MDRf")
```

### 3. Automated mapping function to run continuous binomial GLM models in batch

``` r
run_continuous_trend_glm <- function(variable_name, data) {
  model_formula <- as.formula(paste(variable_name, "~ Year_Cont"))
  
  model <- tryCatch({
    glm(model_formula, data = data, family = binomial(link = "logit"))
  }, error = function(e) { NULL })
  
  if(is.null(model)) {
    return(data.frame(Variable = variable_name, Beta = NA, Odds_Ratio = NA, P_Value = NA, Interpretation = "Monomorphic"))
  }
  
  summary_matrix <- summary(model)$coefficients
  
  if(!"Year_Cont" %in% rownames(summary_matrix)) {
    return(data.frame(Variable = variable_name, Beta = NA, Odds_Ratio = NA, P_Value = NA, Interpretation = "No variance"))
  }
  
  beta_val <- summary_matrix["Year_Cont", "Estimate"]
  p_val    <- summary_matrix["Year_Cont", "Pr(>|z|)"]
  or_val   <- exp(beta_val)
  
  # Categorize long-term multiannual directional behaviors
  interpretation <- ifelse(p_val < 0.05, 
                           ifelse(or_val > 1, "SIGNIFICANT INCREASE 📈", "SIGNIFICANT DECREASE 📉"), 
                           "Stationary Longitudinal Dynamic (Stable)")
  
  data.frame(
    Variable = variable_name,
    Beta = round(beta_val, 4),
    Odds_Ratio = round(or_val, 4),
    P_Value = round(p_val, 4),
    Interpretation = interpretation
  )
}
```

### 4. Execute the batch calculation matrix across all targets

``` r
table_online_resource_8 <- target_variables %>%
  map_df(~ run_continuous_trend_glm(.x, base_continuous_trends))
```

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

### 5. Print definitive metrics directly to the RStudio console

``` r
print("ONLINE RESOURCE 8: ANTIBIOTICS AND MDRf LONGITUDINAL TRENDS")
```

    ## [1] "ONLINE RESOURCE 8: ANTIBIOTICS AND MDRf LONGITUDINAL TRENDS"

``` r
print(as.data.frame(table_online_resource_8))
```

    ##    Variable     Beta Odds_Ratio P_Value
    ## 1       AMP  -0.0643     0.9377  0.3975
    ## 2       AMC  -0.4610     0.6307  0.0646
    ## 3       CRO  -0.6130     0.5417  0.0091
    ## 4       FEP -17.2015     0.0000  0.9956
    ## 5       MEM  -0.4548     0.6346  0.3146
    ## 6       CIP  -0.1716     0.8423  0.2024
    ## 7       AZM  -0.2846     0.7523  0.0027
    ## 8       TET   0.0837     1.0873  0.1670
    ## 9       CHL   0.0257     1.0261  0.7048
    ## 10      AMK   0.2981     1.3473  0.1696
    ## 11      STR   0.0366     1.0373  0.5520
    ## 12      SXT   0.0914     1.0958  0.1898
    ## 13     MDRf   0.0445     1.0455  0.4844
    ##                              Interpretation
    ## 1  Stationary Longitudinal Dynamic (Stable)
    ## 2  Stationary Longitudinal Dynamic (Stable)
    ## 3                   SIGNIFICANT DECREASE 📉
    ## 4  Stationary Longitudinal Dynamic (Stable)
    ## 5  Stationary Longitudinal Dynamic (Stable)
    ## 6  Stationary Longitudinal Dynamic (Stable)
    ## 7                   SIGNIFICANT DECREASE 📉
    ## 8  Stationary Longitudinal Dynamic (Stable)
    ## 9  Stationary Longitudinal Dynamic (Stable)
    ## 10 Stationary Longitudinal Dynamic (Stable)
    ## 11 Stationary Longitudinal Dynamic (Stable)
    ## 12 Stationary Longitudinal Dynamic (Stable)
    ## 13 Stationary Longitudinal Dynamic (Stable)

### 6. Automated spreadsheet export

``` r
write.csv(table_online_resource_8, "OR8_Longitudinal_Trends_PHEN.csv", row.names = FALSE)
print("Unified continuous trend matrix exported successfully as 'OR8_Longitudinal_Trends_PHEN.csv'!")
```

    ## [1] "Unified continuous trend matrix exported successfully as 'OR8_Longitudinal_Trends_PHEN.csv'!"

------------------------------------------------------------------------

## VI. TARGETED CATEGORICAL GLM: ANNUAL TRENDS FOR AZM AND CRO

### 1. Define the core mapping function for categorical multi-level trend tests

``` r
run_targeted_year_glm <- function(variable, data) {
  # Enforce factor layout establishing the pre-pandemic year (2019) as the baseline group
  data_model <- data %>% 
    mutate(Year = factor(Year, levels = c("2019", "2020", "2021", "2022", "2023")))
  
  model_formula <- as.formula(paste(variable, "~ Year"))
  
  model <- tryCatch({ 
    glm(model_formula, data = data_model, family = binomial(link = "logit")) 
  }, error = function(e) { NULL })
  
  if(is.null(model)) return(NULL)
  
  summary_matrix   <- summary(model)$coefficients
  years_to_extract <- c("Year2020", "Year2021", "Year2022", "Year2023")
  
  map_df(years_to_extract, function(yr) {
    if(!yr %in% rownames(summary_matrix)) return(NULL)
    
    p_val  <- summary_matrix[yr, "Pr(>|z|)"]
    or_val <- exp(summary_matrix[yr, "Estimate"])
    
    # Mathematical classification of annual deviations vs 2019 baseline
    interpretation <- ifelse(p_val < 0.05, 
                             ifelse(or_val > 1, "SIGNIFICANT INCREASE 📈", "SIGNIFICANT DECREASE 📉"), 
                             "Stable (No change vs 2019)")
    
    data.frame(
      Antibiotic = variable,
      Year_Comparison = gsub("Year", "", yr),
      Odds_Ratio = round(or_val, 4),
      P_Value = round(p_val, 4),
      Interpretation = interpretation
    )
  })
}
```

### 2. Execute the batch architecture strictly for targeted drugs of interest

``` r
targeted_drugs <- c("AZM", "CRO")

table_targeted_year_trends <- targeted_drugs %>% 
  map_df(~ run_targeted_year_glm(.x, cleaned_master_dataset))
```

### 3. Print definitive metrics to the RStudio console

``` r
print("TARGETED ANNUAL AMR TRENDS (AZM AND CRO DISCRETE ANALYSIS)")
```

    ## [1] "TARGETED ANNUAL AMR TRENDS (AZM AND CRO DISCRETE ANALYSIS)"

``` r
print(as.data.frame(table_targeted_year_trends))
```

    ##   Antibiotic Year_Comparison Odds_Ratio P_Value             Interpretation
    ## 1        AZM            2020     0.2823  0.0127    SIGNIFICANT DECREASE 📉
    ## 2        AZM            2021     0.7760  0.4909 Stable (No change vs 2019)
    ## 3        AZM            2022     0.3905  0.0164    SIGNIFICANT DECREASE 📉
    ## 4        AZM            2023     0.2558  0.0037    SIGNIFICANT DECREASE 📉
    ## 5        CRO            2020     0.4448  0.3132 Stable (No change vs 2019)
    ## 6        CRO            2021     0.6422  0.5219 Stable (No change vs 2019)
    ## 7        CRO            2022     0.1476  0.0733 Stable (No change vs 2019)
    ## 8        CRO            2023     0.0000  0.9873 Stable (No change vs 2019)

### 4. Automated matrix export to active workspace for GitHub documentation

``` r
write.csv(table_targeted_year_trends, "OR9_Targeted_Year_Trend_CRO_AZM.csv", row.names = FALSE)

print("Targeted annual discrete table successfully exported as 'OR9_Targeted_Year_Trend_CRO_AZM.csv'!")
```

    ## [1] "Targeted annual discrete table successfully exported as 'OR9_Targeted_Year_Trend_CRO_AZM.csv'!"

------------------------------------------------------------------------

## VII. BATCH CONTINUOUS GLM: LONGITUDINAL TRENDS FOR GENES/MDR GENOTYPES

### 1. Enforce the year vector as a continuous numerical covariate at the isolate level

``` r
base_genotypic_trends <- AMRg_veganb %>%
  mutate(
    Year_Cont = as.numeric(as.character(Year)),
    MDRg      = as.numeric(as.character(MDRg))
  )
```

### 2. Define the array of the 14 core resistance genes plus the global genotypic MDRg marker

``` r
target_genes <- c("aac", "qnr", "tetA", "sul", "floR", "fosA", "bla", "mph", "gyrA", "bleO", "lnuF", "sat2", "mcr9.1", "arr3", "MDRg")
```

### 3. Automated mapping function to run continuous binomial GLM models in batch

``` r
run_gene_trend_glm <- function(variable_name, data) {
  model_formula <- as.formula(paste(variable_name, "~ Year_Cont"))
  
  model <- tryCatch({
    glm(model_formula, data = data, family = binomial(link = "logit"))
  }, error = function(e) { NULL })
  
  if(is.null(model)) {
    return(data.frame(Variable = variable_name, Beta = NA, Odds_Ratio = NA, P_Value = NA, Trend = "Monomorphic"))
  }
  
  summary_matrix <- summary(model)$coefficients
  
  if(!"Year_Cont" %in% rownames(summary_matrix)) {
    return(data.frame(Variable = variable_name, Beta = NA, Odds_Ratio = NA, P_Value = NA, Trend = "No variance"))
  }
  
  beta_val <- summary_matrix["Year_Cont", "Estimate"]
  p_val    <- summary_matrix["Year_Cont", "Pr(>|z|)"]
  or_val   <- exp(beta_val)
  
  # Categorize long-term multiannual directional behaviors to match your table style
  interpretation <- ifelse(p_val < 0.05, 
                           ifelse(or_val > 1, "Significant Increase 📈", "Significant Decrease 📉"), 
                           "Stable (Non-significant)")
  
  data.frame(
    Variable = variable_name,
    Beta = round(beta_val, 4),
    Odds_Ratio = round(or_val, 4),
    P_Value = round(p_val, 4),
    Trend = interpretation
  )
}
```

### 4. Execute the batch calculation matrix across all target genes and MDRg

``` r
table_online_resource_10 <- target_genes %>%
  map_df(~ run_gene_trend_glm(.x, base_genotypic_trends))
```

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

### 5. Print definitive metrics directly to the RStudio console

``` r
print("ONLINE RESOURCE 10: GENES AND MDRg LONGITUDINAL TRENDS")
```

    ## [1] "ONLINE RESOURCE 10: GENES AND MDRg LONGITUDINAL TRENDS"

``` r
print(as.data.frame(table_online_resource_10))
```

    ##    Variable     Beta Odds_Ratio P_Value                    Trend
    ## 1       aac   0.0440     1.0450  0.4920 Stable (Non-significant)
    ## 2       qnr   0.0235     1.0238  0.6887 Stable (Non-significant)
    ## 3      tetA   0.0387     1.0395  0.5238 Stable (Non-significant)
    ## 4       sul   0.0856     1.0894  0.1692 Stable (Non-significant)
    ## 5      floR   0.0062     1.0063  0.9273 Stable (Non-significant)
    ## 6      fosA   0.0232     1.0235  0.7559 Stable (Non-significant)
    ## 7       bla  -0.0143     0.9858  0.8527 Stable (Non-significant)
    ## 8       mph  -0.1133     0.8929  0.3004 Stable (Non-significant)
    ## 9      gyrA  -0.0095     0.9906  0.9408 Stable (Non-significant)
    ## 10     bleO  -0.0496     0.9516  0.7385 Stable (Non-significant)
    ## 11     lnuF   0.2630     1.3008  0.1283 Stable (Non-significant)
    ## 12     sat2   0.3571     1.4292  0.4006 Stable (Non-significant)
    ## 13   mcr9.1  -0.4548     0.6346  0.3146 Stable (Non-significant)
    ## 14     arr3 -17.2015     0.0000  0.9956 Stable (Non-significant)
    ## 15     MDRg   0.0693     1.0718  0.2507 Stable (Non-significant)

### 6. Automated spreadsheet export to active workspace for GitHub documentation

``` r
write.csv(table_online_resource_10, "OR10_Longitudinal_Linear_Trends_AMRg.csv", row.names = FALSE)

print("Unified genotypic trend matrix exported successfully as 'OR10_Longitudinal_Linear_Trends_AMRg.csv.csv'!")
```

    ## [1] "Unified genotypic trend matrix exported successfully as 'OR10_Longitudinal_Linear_Trends_AMRg.csv.csv'!"

------------------------------------------------------------------------

## VIII. GENE DISTANCE MATRIX PREPARATION (CONTROL FOR MONOMORPHIC SUSCEPTIBILITY)

### 1. Isolate the pure binary matrix of the 14 core target resistance genes

``` r
pure_gene_matrix <- cleaned_master_dataset %>%
  select(aac, qnr, tetA, sul, floR, fosA, bla, mph, gyrA, bleO, lnuF, sat2, mcr9.1, arr3) %>%
  as.matrix()

rownames(pure_gene_matrix) <- cleaned_master_dataset$ID
```

### 2. Filter out completely susceptible isolates to maintain multivariate stability

``` r
genes_per_isolate <- rowSums(pure_gene_matrix)
multivariate_gene_matrix <- pure_gene_matrix[genes_per_isolate > 0, ]
multivariate_metadata    <- cleaned_master_dataset[genes_per_isolate > 0, ]

print(paste("Isolates retained for multivariate ecological analyses:", nrow(multivariate_gene_matrix)))
```

    ## [1] "Isolates retained for multivariate ecological analyses: 360"

### 3. Compute the distance-based matrix utilizing the Jaccard metric

``` r
jaccard_distance_matrix <- vegdist(multivariate_gene_matrix, method = "jaccard", binary = TRUE)
```

------------------------------------------------------------------------

## IX. MULTIVARIATE AMR SPATIAL-TEMPORAL DYNAMICS ANALYSIS (PERMANOVA & PERMDISP)

### 1. Run PERMANOVA x year - adonis2 function

``` r
print("TEMPORAL BASAL STABILITY (GLOBAL YEAR PERMANOVA)")
```

    ## [1] "TEMPORAL BASAL STABILITY (GLOBAL YEAR PERMANOVA)"

``` r
temporal_permanova <- adonis2(jaccard_distance_matrix ~ Year, data = multivariate_metadata, permutations = 999)
print(temporal_permanova)
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = jaccard_distance_matrix ~ Year, data = multivariate_metadata, permutations = 999)
    ##           Df SumOfSqs      R2      F Pr(>F)  
    ## Model      4    1.863 0.01708 1.5423  0.084 .
    ## Residual 355  107.223 0.98292                
    ## Total    359  109.086 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### 2. Run PERMDISP - assess homocedasticity

``` r
print("TEMPORAL HETEROGENEITY ASSESSMENT (PERMDISP BY YEAR)")
```

    ## [1] "TEMPORAL HETEROGENEITY ASSESSMENT (PERMDISP BY YEAR)"

``` r
temporal_dispersion <- betadisper(jaccard_distance_matrix, multivariate_metadata$Year, bias.adjust = TRUE)
```

### 3. Formal permutation test calculation

``` r
temporal_permutest_results <- permutest(temporal_dispersion, permutations = 999)
```

### 4. Extract ANOVA results

``` r
f_stat  <- round(temporal_permutest_results$tab$F[1], 2)
p_value <- round(temporal_permutest_results$tab$`Pr(>F)`[1], 3)
```

### 5. Format label to two decimals

``` r
dynamic_subtitle <- paste0("PERMDISP Homogeneity Test: F = ", f_stat, " (p = ", p_value, ")")
print(dynamic_subtitle)
```

    ## [1] "PERMDISP Homogeneity Test: F = 1.23 (p = 0.291)"

### 6. Export the table to active workspace for GitHub documentation

``` r
permdisp_table_df <- as.data.frame(temporal_permutest_results$tab)
write.csv(permdisp_table_df, "Fig.6B_PERMDISP_Year_Anova.csv", row.names = TRUE)
```

### 7. Build and export the PERMDISP x Year figure

``` r
dispersion_dataframe <- data.frame(
  Year     = as.factor(temporal_dispersion$group),
  Distance = temporal_dispersion$distances
)

fig_6b_permdisp_boxplot <- ggplot(dispersion_dataframe, aes(x = Year, y = Distance, fill = Year)) +
  geom_boxplot(outlier.shape = 16, outlier.alpha = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 1) + 
  scale_fill_brewer(palette = "Set2") + 
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(face = "bold", margin = margin(r = 10))) +
  labs(x = "Year of Sampling", 
       y = "Genomic Distance to Centroid (Jaccard)",
       title = "Salmonella Resistome Dispersion", 
       subtitle = dynamic_subtitle)

ggsave("Fig6B_Boxplot_PERMDISPxYEAR.svg", plot = fig_6b_permdisp_boxplot, device = "svg", units = "in", width = 8, height = 4.5)
```

### 8. Run PERMANOVA x region - adonis2 function

``` r
print("BIOGEOGRAPHICAL REGIONAL SEGREGATION (STATE PERMANOVA)")
```

    ## [1] "BIOGEOGRAPHICAL REGIONAL SEGREGATION (STATE PERMANOVA)"

``` r
regional_permanova <- adonis2(jaccard_distance_matrix ~ State, data = multivariate_metadata, permutations = 999)
print(regional_permanova)
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = jaccard_distance_matrix ~ State, data = multivariate_metadata, permutations = 999)
    ##           Df SumOfSqs      R2     F Pr(>F)  
    ## Model      4    2.141 0.01963 1.777  0.039 *
    ## Residual 355  106.945 0.98037               
    ## Total    359  109.086 1.00000               
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### 9. Compute multivariate distances to centroids by geographic jurisdiction (State)

``` r
print("SPATIAL HETEROGENEITY ASSESSMENT (PERMDISP BY STATE)")
```

    ## [1] "SPATIAL HETEROGENEITY ASSESSMENT (PERMDISP BY STATE)"

``` r
regional_dispersion <- betadisper(jaccard_distance_matrix, multivariate_metadata$State, bias.adjust = TRUE)
```

### 10. Execute the formal permutation test (open-science gold standard)

``` r
regional_permutest_results <- permutest(regional_dispersion, permutations = 999)
```

### 11. Extract exact metrics from the permuted ANOVA matrix for dynamic labels

``` r
f_stat_reg  <- round(regional_permutest_results$tab$F, 2)[1]
p_value_reg <- round(regional_permutest_results$tab$`Pr(>F)`, 3)[1]
```

### 12. Format the dynamic subtitle vector

``` r
dynamic_subtitle_reg <- paste0("PERMDISP Homogeneity Test: F = ", f_stat_reg, " (p = ", p_value_reg, ")")
print(dynamic_subtitle_reg)
```

    ## [1] "PERMDISP Homogeneity Test: F = 2.76 (p = 0.031)"

### 13. Automated spreadsheet export

``` r
permdisp_table_reg_df <- as.data.frame(regional_permutest_results$tab)
write.csv(permdisp_table_reg_df, "Online_Resource_6D_PERMDISP_Regional_Anova.csv", row.names = TRUE)
```

### 14. Build and export the PERMDISP x region

``` r
dispersion_dataframe_reg <- data.frame(
  State = as.factor(regional_dispersion$group), 
  Distance = regional_dispersion$distances
)

fig_6d_permdisp_boxplot_reg <- ggplot(dispersion_dataframe_reg, aes(x = State, y = Distance, fill = State)) +
  geom_boxplot(outlier.shape = 16, outlier.alpha = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 1) + 
  scale_fill_brewer(palette = "Set1") + # Distinct palette to contrast with the annual plot
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(face = "bold", margin = margin(r = 10))) +
  labs(x = "Geographic Jurisdiction (State)", 
       y = "Genomic Distance to Centroid (Jaccard)",
       title = "Salmonella Resistome Dispersion by Region", 
       subtitle = dynamic_subtitle_reg)

ggsave("Fig_6D_Boxplot_PERMDISPxSTATE.svg", plot = fig_6d_permdisp_boxplot_reg, device = "svg", units = "in", width = 8, height = 4.5)
```

------------------------------------------------------------------------

## X. PERMDISP POST-HOC PAIRWISE COMPARISONS (STATE VS STATE)

### 1. Execute Tukey’s Honest Significant Difference test on the regional dispersion object

``` r
print("EXECUTING MULTIPLE PAIRWISE COMPARISONS FOR REGIONAL DISPERSION")
```

    ## [1] "EXECUTING MULTIPLE PAIRWISE COMPARISONS FOR REGIONAL DISPERSION"

``` r
regional_tukey_results <- TukeyHSD(regional_dispersion)
```

### 2. Extract and convert the matrix into a dataframe using position matrix vectors

``` r
tukey_matrix_df <- as.data.frame(regional_tukey_results$group)
```

### 3. Secure variables by column positions to prevent space-naming syntax errors

``` r
tukey_matrix_df$Contrast       <- rownames(tukey_matrix_df)
tukey_matrix_df$P_Adjusted     <- tukey_matrix_df[, 4] # Force extraction of the 4th column (p adj)
tukey_matrix_df$Interpretation <- ifelse(tukey_matrix_df$P_Adjusted < 0.05, 
                                         "Significantly Different 🌟", 
                                         "Statistically Equivalent (Stable)")
```

### 4. Clean and reorder the final structural layout

``` r
tukey_matrix_df <- tukey_matrix_df %>%
  select(Contrast, diff, lwr, upr, P_Adjusted, Interpretation)
```

### 5. Print definitive metrics directly to the RStudio console for inspection

``` r
print("PERMDISP POST-HOC TUKEY HSD MATRIX")
```

    ## [1] "PERMDISP POST-HOC TUKEY HSD MATRIX"

``` r
print(as.data.frame(tukey_matrix_df))
```

    ##                                                Contrast         diff
    ## Mexico City-Hidalgo                 Mexico City-Hidalgo -0.013090583
    ## Morelos-Hidalgo                         Morelos-Hidalgo -0.009566798
    ## State of Mexico-Hidalgo         State of Mexico-Hidalgo -0.025038557
    ## Tlaxcala-Hidalgo                       Tlaxcala-Hidalgo  0.038836555
    ## Morelos-Mexico City                 Morelos-Mexico City  0.003523785
    ## State of Mexico-Mexico City State of Mexico-Mexico City -0.011947974
    ## Tlaxcala-Mexico City               Tlaxcala-Mexico City  0.051927138
    ## State of Mexico-Morelos         State of Mexico-Morelos -0.015471759
    ## Tlaxcala-Morelos                       Tlaxcala-Morelos  0.048403353
    ## Tlaxcala-State of Mexico       Tlaxcala-State of Mexico  0.063875113
    ##                                     lwr        upr  P_Adjusted
    ## Mexico City-Hidalgo         -0.09845886 0.07227769 0.993446594
    ## Morelos-Hidalgo             -0.08817527 0.06904167 0.997319190
    ## State of Mexico-Hidalgo     -0.08331783 0.03324072 0.764009104
    ## Tlaxcala-Hidalgo            -0.02753943 0.10521254 0.495622374
    ## Morelos-Mexico City         -0.08866146 0.09570903 0.999972651
    ## State of Mexico-Mexico City -0.08754758 0.06365163 0.992643139
    ## Tlaxcala-Mexico City        -0.03007651 0.13393079 0.413066361
    ## State of Mexico-Morelos     -0.08334554 0.05240202 0.971006881
    ## Tlaxcala-Morelos            -0.02653761 0.12334432 0.392249399
    ## Tlaxcala-State of Mexico     0.01064608 0.11710415 0.009659112
    ##                                                Interpretation
    ## Mexico City-Hidalgo         Statistically Equivalent (Stable)
    ## Morelos-Hidalgo             Statistically Equivalent (Stable)
    ## State of Mexico-Hidalgo     Statistically Equivalent (Stable)
    ## Tlaxcala-Hidalgo            Statistically Equivalent (Stable)
    ## Morelos-Mexico City         Statistically Equivalent (Stable)
    ## State of Mexico-Mexico City Statistically Equivalent (Stable)
    ## Tlaxcala-Mexico City        Statistically Equivalent (Stable)
    ## State of Mexico-Morelos     Statistically Equivalent (Stable)
    ## Tlaxcala-Morelos            Statistically Equivalent (Stable)
    ## Tlaxcala-State of Mexico           Significantly Different 🌟

### 6. Automated spreadsheet export

``` r
write.csv(tukey_matrix_df, "OR11_PERMDISP_STATE_TukeyHSD.csv", row.names = FALSE)

print("Pairwise regional dispersion table successfully exported as 'OR11_PERMDISP_STATE_TukeyHSD.csv'!")
```

    ## [1] "Pairwise regional dispersion table successfully exported as 'OR11_PERMDISP_STATE_TukeyHSD.csv'!"

------------------------------------------------------------------------

# XI. NMDS ORDINATION ANALYSIS

### 1. Execute NMDS Ordination

``` r
set.seed(42) # For reproducible ordination pathways
nmds_global_ordination <- metaMDS(jaccard_distance_matrix, k = 2, trymax = 100)
```

    ## Run 0 stress 0.07977747 
    ## Run 1 stress 0.07259967 
    ## ... New best solution
    ## ... Procrustes: rmse 0.03594635  max resid 0.1610725 
    ## Run 2 stress 0.1112001 
    ## Run 3 stress 0.06378234 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01513456  max resid 0.142633 
    ## Run 4 stress 0.07449597 
    ## Run 5 stress 0.06337499 
    ## ... New best solution
    ## ... Procrustes: rmse 0.005416965  max resid 0.08913344 
    ## Run 6 stress 0.06591988 
    ## Run 7 stress 0.07404935 
    ## Run 8 stress 0.07004779 
    ## Run 9 stress 0.06412724 
    ## Run 10 stress 0.06359497 
    ## ... Procrustes: rmse 0.006437655  max resid 0.09083835 
    ## Run 11 stress 0.07319657 
    ## Run 12 stress 0.06412897 
    ## Run 13 stress 0.06437485 
    ## Run 14 stress 0.06674819 
    ## Run 15 stress 0.07263819 
    ## Run 16 stress 0.06337441 
    ## ... New best solution
    ## ... Procrustes: rmse 0.00644196  max resid 0.1093249 
    ## Run 17 stress 0.06284137 
    ## ... New best solution
    ## ... Procrustes: rmse 0.007210336  max resid 0.1043885 
    ## Run 18 stress 0.08144745 
    ## Run 19 stress 0.06581399 
    ## Run 20 stress 0.07211682 
    ## Run 21 stress 0.07866687 
    ## Run 22 stress 0.06328271 
    ## ... Procrustes: rmse 0.009396284  max resid 0.1099874 
    ## Run 23 stress 0.07180293 
    ## Run 24 stress 0.07056861 
    ## Run 25 stress 0.1325611 
    ## Run 26 stress 0.1379119 
    ## Run 27 stress 0.1405985 
    ## Run 28 stress 0.07038649 
    ## Run 29 stress 0.06911575 
    ## Run 30 stress 0.06341154 
    ## Run 31 stress 0.06569768 
    ## Run 32 stress 0.07034931 
    ## Run 33 stress 0.06899443 
    ## Run 34 stress 0.08417428 
    ## Run 35 stress 0.07204029 
    ## Run 36 stress 0.07798496 
    ## Run 37 stress 0.07291592 
    ## Run 38 stress 0.07791013 
    ## Run 39 stress 0.06349121 
    ## Run 40 stress 0.06437708 
    ## Run 41 stress 0.06279925 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001129276  max resid 0.01996707 
    ## Run 42 stress 0.06445615 
    ## Run 43 stress 0.06445458 
    ## Run 44 stress 0.07830258 
    ## Run 45 stress 0.111123 
    ## Run 46 stress 0.08164308 
    ## Run 47 stress 0.06384022 
    ## Run 48 stress 0.0632131 
    ## ... Procrustes: rmse 0.005579238  max resid 0.09918055 
    ## Run 49 stress 0.09037566 
    ## Run 50 stress 0.07693719 
    ## Run 51 stress 0.07410031 
    ## Run 52 stress 0.0737297 
    ## Run 53 stress 0.06622196 
    ## Run 54 stress 0.0754978 
    ## Run 55 stress 0.06492249 
    ## Run 56 stress 0.06487799 
    ## Run 57 stress 0.06488426 
    ## Run 58 stress 0.0637298 
    ## Run 59 stress 0.06490943 
    ## Run 60 stress 0.07215971 
    ## Run 61 stress 0.06854598 
    ## Run 62 stress 0.07050867 
    ## Run 63 stress 0.06555281 
    ## Run 64 stress 0.06873796 
    ## Run 65 stress 0.06293966 
    ## ... Procrustes: rmse 0.006117456  max resid 0.1056324 
    ## Run 66 stress 0.06290655 
    ## ... Procrustes: rmse 0.001302236  max resid 0.01468921 
    ## Run 67 stress 0.06846427 
    ## Run 68 stress 0.06972807 
    ## Run 69 stress 0.0906524 
    ## Run 70 stress 0.1417337 
    ## Run 71 stress 0.06541527 
    ## Run 72 stress 0.06606676 
    ## Run 73 stress 0.06348802 
    ## Run 74 stress 0.06335445 
    ## Run 75 stress 0.06570895 
    ## Run 76 stress 0.0732672 
    ## Run 77 stress 0.06404444 
    ## Run 78 stress 0.07688129 
    ## Run 79 stress 0.07835964 
    ## Run 80 stress 0.07167896 
    ## Run 81 stress 0.07339101 
    ## Run 82 stress 0.06321341 
    ## ... Procrustes: rmse 0.005474897  max resid 0.0996498 
    ## Run 83 stress 0.08533661 
    ## Run 84 stress 0.08225404 
    ## Run 85 stress 0.06351637 
    ## Run 86 stress 0.07980575 
    ## Run 87 stress 0.06962592 
    ## Run 88 stress 0.08774022 
    ## Run 89 stress 0.06279895 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001521687  max resid 0.001076278 
    ## ... Similar to previous best
    ## *** Best solution repeated 1 times

``` r
print(paste("Final NMDS Stress Coefficient Value:", round(nmds_global_ordination$stress, 4)))
```

    ## [1] "Final NMDS Stress Coefficient Value: 0.0628"

``` r
print("CALCULATING COMPOSITIONAL PERMANOVA TESTING FRAMES")
```

    ## [1] "CALCULATING COMPOSITIONAL PERMANOVA TESTING FRAMES"

### 2. Visual optimization and figure export. Build coordinate frame for ggplot layouts

``` r
ordination_coordinates <- as.data.frame(scores(nmds_global_ordination, display = "sites"))
ordination_coordinates$Year  <- multivariate_metadata$Year
ordination_coordinates$State <- multivariate_metadata$State
```

### 5. Dynamic subtitle extraction

``` r
r2_temporal <- round(temporal_permanova$R2[1], 4)
p_temporal  <- round(temporal_permanova$`Pr(>F)`[1], 4)
p_label_temp <- ifelse(p_temporal == 0, "p < 0.001", paste0("p = ", p_temporal))
subtitle_nmds_temp <- paste0("PERMANOVA test: R² = ", r2_temporal, " (", p_label_temp, ")")

r2_regional <- round(regional_permanova$R2[1], 4)
p_regional  <- round(regional_permanova$`Pr(>F)`[1], 4)
p_label_reg <- ifelse(p_regional == 0, "p < 0.001", paste0("p = ", p_regional))
subtitle_nmds_reg <- paste0("PERMANOVA test: R² = ", r2_regional, " (", p_label_reg, ")")
```

### 6. FIGURE 6A: NMDS BY SAMPLING YEAR WITH DYNAMIC LABEL

``` r
fig_6A_temporal_nmds <- ggplot(ordination_coordinates, aes(x = NMDS1, y = NMDS2, color = Year)) +
  geom_point(size = 2.5, alpha = 0.7) +
  stat_ellipse(level = 0.95, linewidth = 1) + 
  scale_color_brewer(palette = "Set2") +       
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), axis.title = element_text(face = "bold")) +
  labs(x = "NMDS1", y = "NMDS2", 
       title = "Resistome temporal structure (2019-2023)",
       subtitle = subtitle_nmds_temp,
       color = "Year")

ggsave("Figure_6A_NMDSxYEAR.svg", plot = fig_6A_temporal_nmds, device = "svg", units = "in", width = 8, height = 4.5)
```

### 7. Figure 6C: NMDS by state with dynamic subtitel Optimized utilizing micro-dispersion jittering, transparency, fixed scales, and contour-only shapes

``` r
fig_6C_spatial_nmds <- ggplot(ordination_coordinates, aes(x = NMDS1, y = NMDS2, color = State)) +
  geom_jitter(aes(shape = Year), size = 3.5, alpha = 0.5, width = 0.04, height = 0.04) +
  stat_ellipse(linewidth = 0.8, level = 0.95) + # Contour lines only (prevents solid filling issues in Inkscape)
  coord_cartesian(xlim = c(-3, 3), ylim = c(-2, 2)) + # Scaled matrix limits
  scale_color_brewer(palette = "Set1") + 
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), legend.box = "vertical", axis.title = element_text(face = "bold")) +
  labs(x = "NMDS1", y = "NMDS2", 
       title = "Resistome spatial partitioning across regions",
       subtitle = subtitle_nmds_reg,
       color = "State")

ggsave("Figure_6C_NMDSxSTATE.svg", plot = fig_6C_spatial_nmds, device = "svg", units = "in", width = 8.5, height = 4.8)
```

------------------------------------------------------------------------

## XII. POPULATION LONGITUDINAL FLOW DIAGRAM & GLOBAL AGREEMENT MDR GEN-PHEN

### 1. Calculate Global Cohen’s Kappa for MDR classifications

``` r
global_confusion_matrix <- table(cleaned_master_dataset$MDRg, cleaned_master_dataset$MDRf)
print("GLOBAL COHEN'S KAPPA AGREEMENT RESULT (MDR CATEGORIES)")
```

    ## [1] "GLOBAL COHEN'S KAPPA AGREEMENT RESULT (MDR CATEGORIES)"

``` r
# Compute and secure the formal agreement matrix object
kappa_calculation_object <- Kappa(global_confusion_matrix)
print(kappa_calculation_object)
```

    ##             value     ASE  z   Pr(>|z|)
    ## Unweighted 0.8085 0.02695 30 9.088e-198
    ## Weighted   0.8085 0.02695 30 9.088e-198

### 2. Extraction of metrics for dynamic subtitle

``` r
kappa_value_extracted <- round(kappa_calculation_object$Unweighted[1], 4)
ase_value_extracted   <- round(kappa_calculation_object$Unweighted[2], 4)
```

### 3. Calculate the Z-value and the rounded P-value

``` r
z_statistic_calc <- kappa_value_extracted / ase_value_extracted
p_value_raw      <- 2 * (1 - pnorm(abs(z_statistic_calc)))
p_value_rounded  <- round(p_value_raw, 4) 
```

### 4. Format p-value labels to 4 digits

``` r
p_label_alluvial <- ifelse(p_value_rounded == 0, "p < 0.0001", paste0("p = ", p_value_rounded))
```

### 5. Build final dynamic subtitle with kappa symbol

``` r
dynamic_alluvial_subtitle <- paste0("Cohen's κ = ", kappa_value_extracted, 
                                    " [ASE = ", ase_value_extracted, "], ", p_label_alluvial)
print(dynamic_alluvial_subtitle)
```

    ## [1] "Cohen's κ = 0.8085 [ASE = 0.0269], p < 0.0001"

### 6. Build alluvial plot and export figure

``` r
alluvial_dataframe <- cleaned_master_dataset %>%
  mutate(
    MDRg_Label = factor(MDRg, levels = c(0, 1), labels = c("Genotype Non-MDR", "Genotype MDR")),
    MDRf_Label = factor(MDRf, levels = c(0, 1), labels = c("Phenotype Non-MDR", "Phenotype MDR"))
  ) %>%
  group_by(Year, MDRg_Label, MDRf_Label) %>%
  tally(name = "Frequency")

fig_7_alluvial_flow <- ggplot(data = alluvial_dataframe, aes(axis1 = Year, axis2 = MDRg_Label, axis3 = MDRf_Label, y = Frequency)) +
  scale_x_discrete(limits = c("Year", "Genotypic MDR", "Phenotypic MDR"), expand = c(.1, .1)) +
  geom_alluvium(aes(fill = Year), alpha = 0.6, width = 1/4) +
  geom_stratum(alpha = 0.8, width = 1/4, fill = "grey95", color = "grey30") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(face = "bold")) +
  labs(title = "Longitudinal flow of genotypic and phenotypic MDR",
       subtitle = dynamic_alluvial_subtitle,
       fill = "Year of sampling")
ggsave("Figure_7_Alluvial_MDR_Flow.svg", plot = fig_7_alluvial_flow, device = "svg", units = "in", width = 8.5, height = 4.8)
```

------------------------------------------------------------------------

## XIII. INDEPENDENT VALIDATION FRAMEWORK: STANDARD LOGISTIC REGRESSIONS (GLM). PARALLEL MATHEMATICAL AUDIT TO BENCHMARK PREVALENCE RISK DRIVERS

### 1. Fit an un-nested univariable logistic regression by omitting the random intercept

``` r
print("RUNNING PARALLEL INTERFERENTIAL BENCHMARKING (STANDARD GLM)")
```

    ## [1] "RUNNING PARALLEL INTERFERENTIAL BENCHMARKING (STANDARD GLM)"

``` r
model_validation_glm <- glm(
  Result ~ Year_Fact + State + Source + Season + pH + Temperature + Turbidity, 
  data = base_muestras_anual, 
  family = binomial(link = "logit")
)
```

### 2. Extract standard errors, z-statistics, and exact Wald confidence intervals

``` r
summary_validation_glm <- summary(model_validation_glm)$coefficients
ci_validation_glm      <- confint.default(model_validation_glm)[rownames(summary_validation_glm), ]
```

### 3. Construct the independent validation benchmarking dataframe

``` r
table_validation_glm <- data.frame(
  Variable   = rownames(summary_validation_glm),
  Odds_Ratio = round(exp(summary_validation_glm[, "Estimate"]), 4),
  CI_Lower   = round(exp(ci_validation_glm[, 1]), 4),
  CI_Upper   = round(exp(ci_validation_glm[, 2]), 4),
  P_Value    = round(summary_validation_glm[, "Pr(>|z|)"], 4)
) %>% 
  mutate(Significance = ifelse(P_Value < 0.05, "SIGNIFICANT 🌟", "Stable"))
```

### 4. Print definitive metrics directly to the RStudio console for inspection

``` r
print("PARALLEL AUDIT MATRIX: STANDARD GLM PREVALENCE COEFFICIENTS")
```

    ## [1] "PARALLEL AUDIT MATRIX: STANDARD GLM PREVALENCE COEFFICIENTS"

``` r
print(as.data.frame(table_validation_glm))
```

    ##                                  Variable Odds_Ratio CI_Lower  CI_Upper P_Value
    ## (Intercept)                   (Intercept)   345.0361  30.8204 3862.6917  0.0000
    ## Year_Fact2020               Year_Fact2020     0.4078   0.2129    0.7813  0.0068
    ## Year_Fact2021               Year_Fact2021     0.5663   0.3397    0.9439  0.0292
    ## Year_Fact2022               Year_Fact2022     0.2593   0.1666    0.4036  0.0000
    ## Year_Fact2023               Year_Fact2023     0.0822   0.0475    0.1421  0.0000
    ## StateMexico City         StateMexico City     1.2329   0.5208    2.9191  0.6339
    ## StateMorelos                 StateMorelos     1.6301   0.8501    3.1256  0.1413
    ## StateState of Mexico StateState of Mexico     1.7348   0.8179    3.6796  0.1510
    ## StateTlaxcala               StateTlaxcala     3.0015   1.3840    6.5094  0.0054
    ## SourceDam                       SourceDam     0.9551   0.4503    2.0258  0.9048
    ## SourcePond                     SourcePond     0.6606   0.3683    1.1848  0.1643
    ## SourceRiver                   SourceRiver     0.9597   0.5961    1.5450  0.8655
    ## SeasonSpring                 SeasonSpring     6.8759   4.1580   11.3705  0.0000
    ## SeasonSummer                 SeasonSummer     1.2839   0.8504    1.9382  0.2344
    ## SeasonWinter                 SeasonWinter     3.9527   2.4158    6.4673  0.0000
    ## pH                                     pH     0.4638   0.3502    0.6144  0.0000
    ## Temperature                   Temperature     1.0071   0.9668    1.0491  0.7332
    ## Turbidity                       Turbidity     1.0012   0.9985    1.0039  0.3799
    ##                        Significance
    ## (Intercept)          SIGNIFICANT 🌟
    ## Year_Fact2020        SIGNIFICANT 🌟
    ## Year_Fact2021        SIGNIFICANT 🌟
    ## Year_Fact2022        SIGNIFICANT 🌟
    ## Year_Fact2023        SIGNIFICANT 🌟
    ## StateMexico City             Stable
    ## StateMorelos                 Stable
    ## StateState of Mexico         Stable
    ## StateTlaxcala        SIGNIFICANT 🌟
    ## SourceDam                    Stable
    ## SourcePond                   Stable
    ## SourceRiver                  Stable
    ## SeasonSpring         SIGNIFICANT 🌟
    ## SeasonSummer                 Stable
    ## SeasonWinter         SIGNIFICANT 🌟
    ## pH                   SIGNIFICANT 🌟
    ## Temperature                  Stable
    ## Turbidity                    Stable

### 5. Automated text export

``` r
sink("Raw_Output_Table_2B_Standard_GLM_Prevalence_Validation.txt")
cat("PARALLEL AUDIT MATRIX: STANDARD GLM PREVALENCE DRIVERS\n")
cat("Executed systematically to validate fixed effects under zero random variance thresholds\n\n")
print(summary(model_validation_glm))
sink()

print("Parallel baseline validation complete. Output secured as 'Raw_Output_Table_2B_Standard_GLM_Prevalence_Validation.txt'!")
```

    ## [1] "Parallel baseline validation complete. Output secured as 'Raw_Output_Table_2B_Standard_GLM_Prevalence_Validation.txt'!"

------------------------------------------------------------------------

## XIV. INDEPENDENT VALIDATION FRAMEWORK: STANDARD LOGISTIC REGRESSIONS (GLM). PARALLEL MATHEMATICAL AUDIT TO BENCHMARK MDR PHENOTYPES DRIVERS

### 1. Enforce strict discrete configurations at the isolate level to align benchmarks

``` r
base_isolates_audit <- AMRg_veganb %>%
  mutate(
    Year_Fact   = factor(Year, levels = c("2019", "2020", "2021", "2022", "2023")),
    State       = as.factor(State),
    Source      = as.factor(Source),
    Season      = as.factor(Season),
    pH          = as.numeric(as.character(pH)),
    Temperature = as.numeric(as.character(Temperature)),
    Turbidity   = as.numeric(as.character(Turbidity)),
    MDRf        = as.numeric(as.character(MDRf))
  ) %>%
  filter(!is.na(Year_Fact) & !is.na(State) & !is.na(Source) & !is.na(Season) & 
           !is.na(MDRf) & !is.na(pH) & !is.na(Temperature) & !is.na(Turbidity))

print(paste("Total independent isolates entering parallel GLM audit for MDRf:", nrow(base_isolates_audit)))
```

    ## [1] "Total independent isolates entering parallel GLM audit for MDRf: 508"

### 2. Fit an un-nested standard logistic regression for Phenotypic MDR using the isolate matrix

``` r
model_validation_mdrf_glm <- glm(
  MDRf ~ Year_Fact + State + Source + Season + pH + Temperature + Turbidity, 
  data = base_isolates_audit, 
  family = binomial(link = "logit")
)
```

### 3. Extract standard errors, z-statistics, and exact Wald confidence intervals

``` r
summary_validation_mdrf_glm <- summary(model_validation_mdrf_glm)$coefficients
ci_validation_mdrf_glm      <- confint.default(model_validation_mdrf_glm)[rownames(summary_validation_mdrf_glm), ]
```

### 4. Construct the independent validation benchmarking dataframe for MDRf

``` r
table_validation_mdrf_glm <- data.frame(
  Variable   = rownames(summary_validation_mdrf_glm),
  Odds_Ratio = round(exp(summary_validation_mdrf_glm[, "Estimate"]), 4),
  CI_Lower   = round(exp(ci_validation_mdrf_glm[, 1]), 4),
  CI_Upper   = round(exp(ci_validation_mdrf_glm[, 2]), 4),
  P_Value    = round(summary_validation_mdrf_glm[, "Pr(>|z|)"], 4)
) %>% 
  mutate(Significance = ifelse(P_Value < 0.05, "SIGNIFICANT 🌟", "Stable"))
```

### 5. Print definitive metrics directly to the RStudio console for inspection

``` r
print("PARALLEL AUDIT MATRIX: STANDARD GLM MDRf COEFFICIENTS")
```

    ## [1] "PARALLEL AUDIT MATRIX: STANDARD GLM MDRf COEFFICIENTS"

``` r
print(as.data.frame(table_validation_mdrf_glm))
```

    ##                                  Variable Odds_Ratio CI_Lower CI_Upper P_Value
    ## (Intercept)                   (Intercept)     1.1403   0.0381  34.1486  0.9397
    ## Year_Fact2020               Year_Fact2020     0.3287   0.1091   0.9900  0.0480
    ## Year_Fact2021               Year_Fact2021     1.0918   0.4822   2.4722  0.8332
    ## Year_Fact2022               Year_Fact2022     0.7361   0.3899   1.3893  0.3444
    ## Year_Fact2023               Year_Fact2023     0.6382   0.2801   1.4541  0.2851
    ## StateMexico City         StateMexico City     0.2993   0.0607   1.4764  0.1385
    ## StateMorelos                 StateMorelos     0.1458   0.0356   0.5973  0.0074
    ## StateState of Mexico StateState of Mexico     0.2662   0.0647   1.0959  0.0668
    ## StateTlaxcala               StateTlaxcala     0.2027   0.0477   0.8610  0.0306
    ## SourceDam                       SourceDam     0.4319   0.0961   1.9409  0.2735
    ## SourcePond                     SourcePond     1.1044   0.4164   2.9289  0.8419
    ## SourceRiver                   SourceRiver     1.8146   0.9787   3.3645  0.0585
    ## SeasonSpring                 SeasonSpring     1.5545   0.7748   3.1187  0.2143
    ## SeasonSummer                 SeasonSummer     0.9163   0.4984   1.6844  0.7783
    ## SeasonWinter                 SeasonWinter     1.8671   0.8142   4.2816  0.1403
    ## pH                                     pH     0.9978   0.6725   1.4806  0.9915
    ## Temperature                   Temperature     0.9977   0.9413   1.0575  0.9394
    ## Turbidity                       Turbidity     1.0028   0.9996   1.0059  0.0825
    ##                        Significance
    ## (Intercept)                  Stable
    ## Year_Fact2020        SIGNIFICANT 🌟
    ## Year_Fact2021                Stable
    ## Year_Fact2022                Stable
    ## Year_Fact2023                Stable
    ## StateMexico City             Stable
    ## StateMorelos         SIGNIFICANT 🌟
    ## StateState of Mexico         Stable
    ## StateTlaxcala        SIGNIFICANT 🌟
    ## SourceDam                    Stable
    ## SourcePond                   Stable
    ## SourceRiver                  Stable
    ## SeasonSpring                 Stable
    ## SeasonSummer                 Stable
    ## SeasonWinter                 Stable
    ## pH                           Stable
    ## Temperature                  Stable
    ## Turbidity                    Stable

### 6. Automated text export

``` r
sink("Raw_Output_Table_2D_Standard_GLM_MDRf_Validation.txt")
cat("PARALLEL AUDIT MATRIX: STANDARD GLM MDRf DRIVERS\n")
cat("Executed systematically to validate isolate fixed effects under zero random variance thresholds\n\n")
print(summary(model_validation_mdrf_glm))
sink()

print("Parallel MDRf baseline validation complete. Output secured as 'Raw_Output_Table_2D_Standard_GLM_MDRf_Validation.txt'!")
```

    ## [1] "Parallel MDRf baseline validation complete. Output secured as 'Raw_Output_Table_2D_Standard_GLM_MDRf_Validation.txt'!"

``` r
print("COMPLETE PIPELINE SUCCESSFULLY EXECUTED")
```

    ## [1] "COMPLETE PIPELINE SUCCESSFULLY EXECUTED"
