AMR and Salmonella Dynamics Pipeline
================
Enrique J. Delgado Suárez
10 July, 2026

# ————————————————————————————————————————

# INTEGRATED PIPELINE

# PROJECT: SALMONELLA IN SURFACE WATERS FROM CENTRAL MEXICO (2019-2023)

# ————————————————————————————————————————

# I. LOAD REQUIRED LIBRARIES, IMPORT AND ADJUST THE EXCEL FILES

library(vegan) \# Multivariate analyses (Jaccard, PERMANOVA, PERMDISP)
library(ggplot2) \# Produce figures library(dplyr) \# Structured data
manipulation and cleaning library(tidyr) \# Data reshaping matrices
library(purrr) \# Functional batch processing execution (map_df)
library(ggalluvial) \# Longitudinal MDR genotype-phenotype alluvial
tracking library(vcd) \# Cohen’s Kappa with ASE library(lme4) \#
Generalized Linear Mixed-Effects Models (GLMM via glmer)
library(svglite) \# Scaled vector graphics engine (SVG export support)
library(readxl) \# Reading Excel files

AMRg_veganb \<- read_excel(“AMRg_veganb.xlsx”)

data_muestras_cruda \<- read_excel(“Prevalence.xlsx”)

# Clean missing fields and explicitly convert character vectors into strict continuous numeric elements to eradicate Excel date-formatting artifacts.

cleaned_master_dataset \<- AMRg_veganb %\>% mutate( \# Force continuous
variables into numeric layout pH = as.numeric(as.character(pH)),
Temperature = as.numeric(as.character(Temperature)), Turbidity =
as.numeric(as.character(Turbidity)),

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

# ————————————————————————————————————————

# II. UNIFIED MODEL: GLOBAL SALMONELLA PREVALENCE (GLMM FULL DISCRETE + CONTINUOUS)

# CONSOLIDATED RISK FACTOR ANALYSIS FOR SAMPLE-LEVEL PREVALENCE

base_muestras_anual \<- data_muestras_cruda %\>% mutate( Year_Fact =
factor(Year, levels = c(“2019”, “2020”, “2021”, “2022”, “2023”)), State
= as.factor(State), Source = as.factor(Source), Season =
as.factor(Season), Site = as.factor(Site),  
pH = as.numeric(as.character(pH)), Temperature =
as.numeric(as.character(Temperature)), Turbidity =
as.numeric(as.character(Turbidity)), Result =
as.numeric(as.character(Result)) ) %\>% filter(!is.na(Year_Fact) &
!is.na(State) & !is.na(Source) & !is.na(Season) & !is.na(Site) &
!is.na(Result) & !is.na(pH) & !is.na(Temperature) & !is.na(Turbidity))

print(paste(“Total independent samples entering consolidated GLMM for
Prevalence:”, nrow(base_muestras_anual)))

model_annual_prevalence \<- glmer( Result ~ Year_Fact + State + Source +
Season + pH + Temperature + Turbidity + (1 \| Site), data =
base_muestras_anual, family = binomial(link = “logit”) )

summary_ann_prev \<- summary(model_annual_prevalence)\$coefficients
ci_ann_prev \<- confint(model_annual_prevalence, method =
“Wald”)\[rownames(summary_ann_prev), \]

table_annual_prevalence \<- data.frame( Variable =
rownames(summary_ann_prev), Odds_Ratio = round(exp(summary_ann_prev\[,
“Estimate”\]), 4), CI_Lower = round(exp(ci_ann_prev\[, 1\]), 4),
CI_Upper = round(exp(ci_ann_prev\[, 2\]), 4), P_Value =
round(summary_ann_prev\[, “Pr(\>\|z\|)”\], 4) ) %\>% mutate(Significance
= ifelse(P_Value \< 0.05, “SIGNIFICANT 🌟”, “Stable”))

print(“DEFINITIVE UNIFIED SAMPLE-LEVEL SALMONELLA PREVALENCE DRIVERS”)
print(as.data.frame(table_annual_prevalence))
write.csv(table_annual_prevalence,
“Table_4A_Unified_Global_Prevalence_Drivers.csv”, row.names = FALSE)

print(“EXPORTING RAW GLMM CONSOLE SUMMARY FOR REPOSITORIES”)
sink(“Raw_Output_Table_2A_GLMM_Prevalence.txt”) cat(“RAW GLMM CONSOLE
SUMMARY: GLOBAL SALMONELLA PREVALENCE”) cat(“Generated automatically via
open-science verification pipelines”)
print(summary(model_annual_prevalence)) \# Reports AIC, BIC, random
effects, and z values sink() print(“Raw console summary successfully
exported as ‘Raw_Output_Table_2A_GLMM_Prevalence.txt’!”) \#
————————————————————————————————————————

# III. UNIFIED MODEL: MULTIDRUG RESISTANCE (GLMM FULL DISCRETE + CONTINUOUS)

# CONSOLIDATED RISK FACTOR ANALYSIS FOR ISOLATE-LEVEL MDR PHENOTYPES

base_mdrf_consolidada \<- AMRg_veganb %\>% mutate( Year_Fact =
factor(Year, levels = c(“2019”, “2020”, “2021”, “2022”, “2023”)), State
= as.factor(State), Source = as.factor(Source), Season =
as.factor(Season), Site = as.factor(Site),  
pH = as.numeric(as.character(pH)), Temperature =
as.numeric(as.character(Temperature)), Turbidity =
as.numeric(as.character(Turbidity)), MDRf =
as.numeric(as.character(MDRf)) ) %\>% filter(!is.na(Year_Fact) &
!is.na(State) & !is.na(Source) & !is.na(Season) & !is.na(Site) &
!is.na(MDRf) & !is.na(pH) & !is.na(Temperature) & !is.na(Turbidity))

print(paste(“Total isolates entering consolidated GLMM for MDRf:”,
nrow(base_mdrf_consolidada)))

model_mdrf_unified \<- glmer( MDRf ~ Year_Fact + State + Source +
Season + pH + Temperature + Turbidity + (1 \| Site), data =
base_mdrf_consolidada, family = binomial(link = “logit”) )

summary_uni_mdrf \<- summary(model_mdrf_unified)\$coefficients
ci_uni_mdrf \<- confint(model_mdrf_unified, method =
“Wald”)\[rownames(summary_uni_mdrf), \]

table_mdrf_unified_output \<- data.frame( Variable =
rownames(summary_uni_mdrf), Odds_Ratio = round(exp(summary_uni_mdrf\[,
“Estimate”\]), 4), CI_Lower = round(exp(ci_uni_mdrf\[, 1\]), 4),
CI_Upper = round(exp(ci_uni_mdrf\[, 2\]), 4), P_Value =
round(summary_uni_mdrf\[, “Pr(\>\|z\|)”\], 4) ) %\>% mutate(Significance
= ifelse(P_Value \< 0.05, “SIGNIFICANT 🌟”, “Stable”))

print(“DEFINITIVE UNIFIED MULTIDRUG RESISTANCE (MDRf) DRIVERS”)
print(as.data.frame(table_mdrf_unified_output))
write.csv(table_mdrf_unified_output,
“Table_4B_Unified_MDRf_Drivers.csv”, row.names = FALSE)

print(“EXPORTING RAW GLMM CONSOLE SUMMARY FOR REPOSITORIES”)
sink(“Raw_Output_Table_2C_GLMM_MDRf.txt”) cat(“RAW GLMM CONSOLE SUMMARY:
GLOBAL MDR PHENOTYPES”) cat(“Generated automatically via open-science
verification pipelines”) print(summary(model_mdrf_unified)) \# Reports
AIC, BIC, random effects, and z values sink() print(“Raw console summary
successfully exported as ‘Raw_Output_Table_2C_GLMM_MDRf.txt’!”)

# ————————————————————————————————————————

# IV. DRUG-LEVEL GENOTYPE-PHENOTYPE CONCORDANCE (BATCH INDIVIDUAL KAPPA CALCULATION)

# Evaluates each antibiotic against its target genetic group, capturing ASE vectors inline

individual_kappa_table \<- list( c(ant = “AMP”, gen = “bla”), c(ant =
“AMC”, gen = “bla”), c(ant = “CRO”, gen = “bla”), c(ant = “FEP”, gen =
“bla”), c(ant = “MEM”, gen = “bla”), c(ant = “CIP”, gen = “gyrA”), c(ant
= “AZM”, gen = “mph”), c(ant = “TET”, gen = “tetA”), c(ant = “CHL”, gen
= “floR”), c(ant = “AMK”, gen = “aac”), c(ant = “STR”, gen = “aac”),
c(ant = “SXT”, gen = “sul”) ) %\>% map_df(function(par) { tabla_par \<-
table(cleaned_master_dataset\[\[par\[“gen”\]\]\],
cleaned_master_dataset\[\[par\[“ant”\]\]\]) if(nrow(tabla_par) != 2 \|\|
ncol(tabla_par) != 2) return(data.frame(Antibiotic = par\[“ant”\],
Associated_Gene = par\[“gen”\], Kappa_or_ASE = NA)) kappa_obj \<-
Kappa(tabla_par) \# Extracts the unweighted vector directly containing
both Coeff and ASE sequentially return(data.frame(Antibiotic =
par\[“ant”\], Associated_Gene = par\[“gen”\], Kappa_or_ASE =
round(as.numeric(kappa_obj\$Unweighted), 4))) })

print(“COHEN’S KAPPA DRUG-LEVEL CONCORDANCE TABLES WITH ASE VALUES”)
print(as.data.frame(individual_kappa_table))
write.csv(individual_kappa_table, “OR5_Kappa_Individual_ATB.csv”,
row.names = FALSE)

# ————————————————————————————————————————

# V. BATCH CONTINUOUS GLM: LONGITUDINAL TRENDS FOR ANTIBIOTICS/MDR PHENOTYPES

# 1. Enforce the year vector as a continuous numerical covariate at the isolate level

base_continuous_trends \<- AMRg_veganb %\>% mutate( Year_Cont =
as.numeric(as.character(Year)), MDRf = as.numeric(as.character(MDRf)) )

# 2. Define the array of the 12 phenotypic antibiotics plus MDR phenotypes (MDRf)

target_variables \<- c(“AMP”, “AMC”, “CRO”, “FEP”, “MEM”, “CIP”, “AZM”,
“TET”, “CHL”, “AMK”, “STR”, “SXT”, “MDRf”)

# 3. Automated mapping function to run continuous binomial GLM models in batch

run_continuous_trend_glm \<- function(variable_name, data) {
model_formula \<- as.formula(paste(variable_name, “~ Year_Cont”))

model \<- tryCatch({ glm(model_formula, data = data, family =
binomial(link = “logit”)) }, error = function(e) { NULL })

if(is.null(model)) { return(data.frame(Variable = variable_name, Beta =
NA, Odds_Ratio = NA, P_Value = NA, Interpretation = “Monomorphic”)) }

summary_matrix \<- summary(model)\$coefficients

if(!“Year_Cont” %in% rownames(summary_matrix)) {
return(data.frame(Variable = variable_name, Beta = NA, Odds_Ratio = NA,
P_Value = NA, Interpretation = “No variance”)) }

beta_val \<- summary_matrix\[“Year_Cont”, “Estimate”\] p_val \<-
summary_matrix\[“Year_Cont”, “Pr(\>\|z\|)”\] or_val \<- exp(beta_val)

\# Categorize long-term multiannual directional behaviors interpretation
\<- ifelse(p_val \< 0.05, ifelse(or_val \> 1, “SIGNIFICANT INCREASE 📈”,
“SIGNIFICANT DECREASE 📉”), “Stationary Longitudinal Dynamic (Stable)”)

data.frame( Variable = variable_name, Beta = round(beta_val, 4),
Odds_Ratio = round(or_val, 4), P_Value = round(p_val, 4), Interpretation
= interpretation ) }

# 4. Execute the batch calculation matrix across all targets

table_online_resource_8 \<- target_variables %\>% map_df(~
run_continuous_trend_glm(.x, base_continuous_trends))

# 5. Print definitive metrics directly to the RStudio console

print(“ONLINE RESOURCE 8: ANTIBIOTICS AND MDRf LONGITUDINAL TRENDS”)
print(as.data.frame(table_online_resource_2a))

# 6. Automated spreadsheet export

write.csv(table_online_resource_8, “OR8_Longitudinal_Trends_PHEN.csv”,
row.names = FALSE)

print(“Unified continuous trend matrix exported successfully as
‘OR8_Longitudinal_Trends_PHEN.csv’!”)

# ————————————————————————————————————————

# VI. TARGETED CATEGORICAL GLM: ANNUAL TRENDS FOR AZM AND CRO

# 1. Define the core mapping function for categorical multi-level trend tests

run_targeted_year_glm \<- function(variable, data) { \# Enforce factor
layout establishing the pre-pandemic year (2019) as the baseline group
data_model \<- data %\>% mutate(Year = factor(Year, levels = c(“2019”,
“2020”, “2021”, “2022”, “2023”)))

model_formula \<- as.formula(paste(variable, “~ Year”))

model \<- tryCatch({ glm(model_formula, data = data_model, family =
binomial(link = “logit”)) }, error = function(e) { NULL })

if(is.null(model)) return(NULL)

summary_matrix \<- summary(model)\$coefficients years_to_extract \<-
c(“Year2020”, “Year2021”, “Year2022”, “Year2023”)

map_df(years_to_extract, function(yr) { if(!yr %in%
rownames(summary_matrix)) return(NULL)

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

}) }

# 2. Execute the batch architecture strictly for targeted drugs of interest

targeted_drugs \<- c(“AZM”, “CRO”)

table_targeted_year_trends \<- targeted_drugs %\>% map_df(~
run_targeted_year_glm(.x, cleaned_master_dataset))

# 3. Print definitive metrics to the RStudio console

print(“TARGETED ANNUAL AMR TRENDS (AZM AND CRO DISCRETE ANALYSIS)”)
print(as.data.frame(table_targeted_year_trends))

# 4. Automated matrix export to active workspace for GitHub documentation

write.csv(table_targeted_year_trends,
“OR9_Targeted_Year_Trend_CRO_AZM.csv”, row.names = FALSE)

print(“Targeted annual discrete table successfully exported as
‘OR9_Targeted_Year_Trend_CRO_AZM.csv’!”)

# ————————————————————————————————————————

# VII. BATCH CONTINUOUS GLM: LONGITUDINAL TRENDS FOR GENES/MDR GENOTYPES

# 1. Enforce the year vector as a continuous numerical covariate at the isolate level

base_genotypic_trends \<- AMRg_veganb %\>% mutate( Year_Cont =
as.numeric(as.character(Year)), MDRg = as.numeric(as.character(MDRg)) )

# 2. Define the array of the 14 core resistance genes plus the global genotypic MDRg marker

target_genes \<- c(“aac”, “qnr”, “tetA”, “sul”, “floR”, “fosA”, “bla”,
“mph”, “gyrA”, “bleO”, “lnuF”, “sat2”, “mcr9.1”, “arr3”, “MDRg”)

# 3. Automated mapping function to run continuous binomial GLM models in batch

run_gene_trend_glm \<- function(variable_name, data) { model_formula \<-
as.formula(paste(variable_name, “~ Year_Cont”))

model \<- tryCatch({ glm(model_formula, data = data, family =
binomial(link = “logit”)) }, error = function(e) { NULL })

if(is.null(model)) { return(data.frame(Variable = variable_name, Beta =
NA, Odds_Ratio = NA, P_Value = NA, Trend = “Monomorphic”)) }

summary_matrix \<- summary(model)\$coefficients

if(!“Year_Cont” %in% rownames(summary_matrix)) {
return(data.frame(Variable = variable_name, Beta = NA, Odds_Ratio = NA,
P_Value = NA, Trend = “No variance”)) }

beta_val \<- summary_matrix\[“Year_Cont”, “Estimate”\] p_val \<-
summary_matrix\[“Year_Cont”, “Pr(\>\|z\|)”\] or_val \<- exp(beta_val)

\# Categorize long-term multiannual directional behaviors to match your
table style interpretation \<- ifelse(p_val \< 0.05, ifelse(or_val \> 1,
“Significant Increase 📈”, “Significant Decrease 📉”), “Stable
(Non-significant)”)

data.frame( Variable = variable_name, Beta = round(beta_val, 4),
Odds_Ratio = round(or_val, 4), P_Value = round(p_val, 4), Trend =
interpretation ) }

# 4. Execute the batch calculation matrix across all target genes and MDRg

table_online_resource_10 \<- target_genes %\>% map_df(~
run_gene_trend_glm(.x, base_genotypic_trends))

# 5. Print definitive metrics directly to the RStudio console

print(“ONLINE RESOURCE 10: GENES AND MDRg LONGITUDINAL TRENDS”)
print(as.data.frame(table_online_resource_10))

# 6. Automated spreadsheet export to active workspace for GitHub documentation

write.csv(table_online_resource_10,
“OR10_Longitudinal_Linear_Trends_AMRg.csv”, row.names = FALSE)

print(“Unified genotypic trend matrix exported successfully as
‘OR10_Longitudinal_Linear_Trends_AMRg.csv.csv’!”)

# ————————————————————————————————————————

# VIII. GENE DISTANCE MATRIX PREPARATION (CONTROL FOR MONOMORPHIC SUSCEPTIBILITY)

# 1. Isolate the pure binary matrix of the 14 core target resistance genes

pure_gene_matrix \<- cleaned_master_dataset %\>% select(aac, qnr, tetA,
sul, floR, fosA, bla, mph, gyrA, bleO, lnuF, sat2, mcr9.1, arr3) %\>%
as.matrix()

rownames(pure_gene_matrix) \<- cleaned_master_dataset\$ID

# 2. Filter out completely susceptible isolates to maintain multivariate stability

genes_per_isolate \<- rowSums(pure_gene_matrix) multivariate_gene_matrix
\<- pure_gene_matrix\[genes_per_isolate \> 0, \] multivariate_metadata
\<- cleaned_master_dataset\[genes_per_isolate \> 0, \]

print(paste(“Isolates retained for multivariate ecological analyses:”,
nrow(multivariate_gene_matrix)))

# 3. Compute the distance-based matrix utilizing the Jaccard metric

jaccard_distance_matrix \<- vegdist(multivariate_gene_matrix, method =
“jaccard”, binary = TRUE)

# ————————————————————————————————————————

# IX. MULTIVARIATE AMR SPATIAL-TEMPORAL DYNAMICS ANALYSIS (PERMANOVA & PERMDISP)

# 1. Run PERMANOVA x year - adonis2 function

print(“TEMPORAL BASAL STABILITY (GLOBAL YEAR PERMANOVA)”)
temporal_permanova \<- adonis2(jaccard_distance_matrix ~ Year, data =
multivariate_metadata, permutations = 999) print(temporal_permanova)

# 2. Run PERMDISP - assess homocedasticity

print(“TEMPORAL HETEROGENEITY ASSESSMENT (PERMDISP BY YEAR)”)
temporal_dispersion \<- betadisper(jaccard_distance_matrix,
multivariate_metadata\$Year, bias.adjust = TRUE)

# 3. Formal permutation test calculation

temporal_permutest_results \<- permutest(temporal_dispersion,
permutations = 999)

# 4. Extract ANOVA results

f_stat \<- round(temporal_permutest_results$tab$F\[1\], 2) p_value \<-
round(temporal_permutest_results$tab$`Pr(>F)`\[1\], 3)

# 5. Format label to two decimals

dynamic_subtitle \<- paste0(“PERMDISP Homogeneity Test: F =”, f_stat, ”
(p = “, p_value,”)“) print(dynamic_subtitle)

# 6. Export the table to active workspace for GitHub documentation

permdisp_table_df \<- as.data.frame(temporal_permutest_results\$tab)
write.csv(permdisp_table_df, “Fig.6B_PERMDISP_Year_Anova.csv”, row.names
= TRUE)

# 7. Build and export the PERMDISP x Year figure

# Fig. 6B: Boxplot of PERMDISP analysis

fig_6b_permdisp_boxplot \<- ggplot(dispersion_dataframe, aes(x = Year, y
= Distance, fill = Year)) + geom_boxplot(outlier.shape = 16,
outlier.alpha = 0.5, width = 0.6) + geom_jitter(width = 0.15, alpha =
0.2, size = 1) + scale_fill_brewer(palette = “Set2”) +
theme_minimal(base_size = 14) + theme(panel.grid.minor =
element_blank(), legend.position = “none”, axis.title.x =
element_text(face = “bold”, margin = margin(t = 10)), axis.title.y =
element_text(face = “bold”, margin = margin(r = 10))) + labs(x = “Year
of Sampling”, y = “Genomic Distance to Centroid (Jaccard)”, title =
“Salmonella Resistome Dispersion”, subtitle = dynamic_subtitle)

ggsave(“Fig6B_Boxplot_PERMDISPxYEAR.svg”, plot =
fig_6b_permdisp_boxplot, device = “svg”, units = “in”, width = 8, height
= 4.5)

# 8. Run PERMANOVA x region - adonis2 function

print(“BIOGEOGRAPHICAL REGIONAL SEGREGATION (STATE PERMANOVA)”)
regional_permanova \<- adonis2(jaccard_distance_matrix ~ State, data =
multivariate_metadata, permutations = 999) print(regional_permanova)

print(“SPATIAL HETEROGENEITY ASSESSMENT (PERMDISP BY STATE)”) \# 9.
Compute multivariate distances to centroids by geographic jurisdiction
(State) regional_dispersion \<- betadisper(jaccard_distance_matrix,
multivariate_metadata\$State, bias.adjust = TRUE)

# 10. Execute the formal permutation test (open-science gold standard)

regional_permutest_results \<- permutest(regional_dispersion,
permutations = 999)

# 11. Extract exact metrics from the permuted ANOVA matrix for dynamic labels

f_stat_reg \<- round(regional_permutest_results$tab$F, 2)\[1\]
p_value_reg \<- round(regional_permutest_results$tab$`Pr(>F)`, 3)\[1\]

# 12. Format the dynamic subtitle vector

dynamic_subtitle_reg \<- paste0(“PERMDISP Homogeneity Test: F =”,
f_stat_reg, ” (p = “, p_value_reg,”)“) print(dynamic_subtitle_reg)

# 13. Automated spreadsheet export

permdisp_table_reg_df \<- as.data.frame(regional_permutest_results\$tab)
write.csv(permdisp_table_reg_df,
“Online_Resource_6D_PERMDISP_Regional_Anova.csv”, row.names = TRUE)

# 14. Build and export the PERMDISP x region

dispersion_dataframe_reg \<- data.frame( State =
as.factor(regional_dispersion$group),  Distance = regional_dispersion$distances
)

fig_6d_permdisp_boxplot_reg \<- ggplot(dispersion_dataframe_reg, aes(x =
State, y = Distance, fill = State)) + geom_boxplot(outlier.shape = 16,
outlier.alpha = 0.5, width = 0.6) + geom_jitter(width = 0.15, alpha =
0.2, size = 1) + scale_fill_brewer(palette = “Set1”) + \# Distinct
palette to contrast with the annual plot theme_minimal(base_size = 14) +
theme(panel.grid.minor = element_blank(), legend.position = “none”,
axis.title.x = element_text(face = “bold”, margin = margin(t = 10)),
axis.title.y = element_text(face = “bold”, margin = margin(r = 10))) +
labs(x = “Geographic Jurisdiction (State)”, y = “Genomic Distance to
Centroid (Jaccard)”, title = “Salmonella Resistome Dispersion by
Region”, subtitle = dynamic_subtitle_reg)

ggsave(“Fig_6D_Boxplot_PERMDISPxSTATE.svg”, plot =
fig_6d_permdisp_boxplot_reg, device = “svg”, units = “in”, width = 8,
height = 4.5)

# ————————————————————————————————————————

# X. PERMDISP POST-HOC PAIRWISE COMPARISONS (STATE VS STATE)

print(“EXECUTING MULTIPLE PAIRWISE COMPARISONS FOR REGIONAL DISPERSION”)

# 1. Execute Tukey’s Honest Significant Difference test on the regional dispersion object

regional_tukey_results \<- TukeyHSD(regional_dispersion)

# 2. Extract and convert the matrix into a dataframe using position matrix vectors

tukey_matrix_df \<- as.data.frame(regional_tukey_results\$group)

# 3. Secure variables by column positions to prevent space-naming syntax errors

tukey_matrix_df$Contrast <- rownames(tukey_matrix_df) tukey_matrix_df$P_Adjusted
\<- tukey_matrix_df\[, 4\] \# Force extraction of the 4th column (p adj)
tukey_matrix_df$Interpretation <- ifelse(tukey_matrix_df$P_Adjusted \<
0.05, “Significantly Different 🌟”, “Statistically Equivalent (Stable)”)

# 4. Clean and reorder the final structural layout

tukey_matrix_df \<- tukey_matrix_df %\>% select(Contrast, diff, lwr,
upr, P_Adjusted, Interpretation)

# 5. Print definitive metrics directly to the RStudio console for inspection

print(“PERMDISP POST-HOC TUKEY HSD MATRIX”)
print(as.data.frame(tukey_matrix_df))

# 6. Automated spreadsheet export

write.csv(tukey_matrix_df, “OR11_PERMDISP_STATE_TukeyHSD.csv”, row.names
= FALSE)

print(“Pairwise regional dispersion table successfully exported as
‘OR11_PERMDISP_STATE_TukeyHSD.csv’!”)

\#————————————————————————————————————————-

# XI. NMDS ORDINATION ANALYSIS

# 1. Execute NMDS Ordination

set.seed(42) \# For reproducible ordination pathways
nmds_global_ordination \<- metaMDS(jaccard_distance_matrix, k = 2,
trymax = 100) print(paste(“Final NMDS Stress Coefficient Value:”,
round(nmds_global_ordination\$stress, 4))) print(“CALCULATING
COMPOSITIONAL PERMANOVA TESTING FRAMES”)

# 2. Visual optimization and figure export

# Build coordinate frame for ggplot layouts

ordination_coordinates \<- as.data.frame(scores(nmds_global_ordination,
display = “sites”))
ordination_coordinates$Year <- multivariate_metadata$Year
ordination_coordinates$State <- multivariate_metadata$State

# 5. Dynamic subtitle extraction

r2_temporal \<-
round(temporal_permanova$R2[1], 4) p_temporal <- round(temporal_permanova$`Pr(>F)`\[1\],
4) p_label_temp \<- ifelse(p_temporal == 0, “p \< 0.001”, paste0(“p =”,
p_temporal)) subtitle_nmds_temp \<- paste0(“PERMANOVA test: R² =”,
r2_temporal, ” (“, p_label_temp,”)“)

r2_regional \<-
round(regional_permanova$R2[1], 4) p_regional <- round(regional_permanova$`Pr(>F)`\[1\],
4) p_label_reg \<- ifelse(p_regional == 0, “p \< 0.001”, paste0(“p =”,
p_regional)) subtitle_nmds_reg \<- paste0(“PERMANOVA test: R² =”,
r2_regional, ” (“, p_label_reg,”)“)

# 6. FIGURE 6A: NMDS BY SAMPLING YEAR WITH DYNAMIC LABEL

fig_6A_temporal_nmds \<- ggplot(ordination_coordinates, aes(x = NMDS1, y
= NMDS2, color = Year)) + geom_point(size = 2.5, alpha = 0.7) +
stat_ellipse(level = 0.95, linewidth = 1) + scale_color_brewer(palette =
“Set2”) +  
theme_minimal(base_size = 14) + theme(panel.grid.minor =
element_blank(), axis.title = element_text(face = “bold”)) + labs(x =
“NMDS1”, y = “NMDS2”, title = “Resistome temporal structure
(2019-2023)”, subtitle = subtitle_nmds_temp, color = “Year”)

ggsave(“Figure_6A_NMDSxYEAR.svg”, plot = fig_6A_temporal_nmds, device =
“svg”, units = “in”, width = 8, height = 4.5)

# 7. FIGURE 6C: NMDS BY STATE WITH DYNAMIC SUBTITLE

# Optimized utilizing micro-dispersion jittering, transparency, fixed scales, and contour-only shapes

fig_6C_spatial_nmds \<- ggplot(ordination_coordinates, aes(x = NMDS1, y
= NMDS2, color = State)) + geom_jitter(aes(shape = Year), size = 3.5,
alpha = 0.5, width = 0.04, height = 0.04) + stat_ellipse(linewidth =
0.8, level = 0.95) + \# Contour lines only (prevents solid filling
issues in Inkscape) coord_cartesian(xlim = c(-3, 3), ylim = c(-2, 2)) +
\# Scaled matrix limits scale_color_brewer(palette = “Set1”) +
theme_minimal(base_size = 14) + theme(panel.grid.minor =
element_blank(), legend.box = “vertical”, axis.title = element_text(face
= “bold”)) + labs(x = “NMDS1”, y = “NMDS2”, title = “Resistome spatial
partitioning across regions”, subtitle = subtitle_nmds_reg, color =
“State”)

ggsave(“Figure_6C_NMDSxSTATE.svg”, plot = fig_6C_spatial_nmds, device =
“svg”, units = “in”, width = 8.5, height = 4.8)

# ————————————————————————————————————————

# XII. POPULATION LONGITUDINAL FLOW DIAGRAM & GLOBAL AGREEMENT MDR GEN-PHEN

# 1. Calculate Global Cohen’s Kappa for MDR classifications

global_confusion_matrix \<-
table(cleaned_master_dataset$MDRg, cleaned_master_dataset$MDRf)
print(“GLOBAL COHEN’S KAPPA AGREEMENT RESULT (MDR CATEGORIES)”)
kappa_calculation_object \<- Kappa(global_confusion_matrix)
print(Kappa(global_confusion_matrix))

# 2. Extraction of metrics for dynamic subtitle

kappa_value_extracted \<-
round(kappa_calculation_object$Unweighted[1], 4) ase_value_extracted <- round(kappa_calculation_object$Unweighted\[2\],
4)

# 3. Calculate the Z-value and the rounded P-value

z_statistic_calc \<- kappa_value_extracted / ase_value_extracted
p_value_raw \<- 2 \* (1 - pnorm(abs(z_statistic_calc))) p_value_rounded
\<- round(p_value_raw, 4)

# 4. Format p-value labels to 4 digits

p_label_alluvial \<- ifelse(p_value_rounded == 0, “p \< 0.0001”,
paste0(“p =”, p_value_rounded))

# 5. Build final dynamic subtitle with kappa symbol

dynamic_alluvial_subtitle \<- paste0(“Cohen’s κ =”,
kappa_value_extracted, ” \[ASE = “, ase_value_extracted,”\], “,
p_label_alluvial)

# 6. Build alluvial plot and export figure

alluvial_dataframe \<- cleaned_master_dataset %\>% mutate( MDRg_Label =
factor(MDRg, levels = c(0, 1), labels = c(“Genotype Non-MDR”, “Genotype
MDR”)), MDRf_Label = factor(MDRf, levels = c(0, 1), labels =
c(“Phenotype Non-MDR”, “Phenotype MDR”)) ) %\>% group_by(Year,
MDRg_Label, MDRf_Label) %\>% tally(name = “Frequency”)

fig_7_alluvial_flow \<- ggplot(data = alluvial_dataframe, aes(axis1 =
Year, axis2 = MDRg_Label, axis3 = MDRf_Label, y = Frequency)) +
scale_x_discrete(limits = c(“Year”, “Genotypic MDR”, “Phenotypic MDR”),
expand = c(.1, .1)) + geom_alluvium(aes(fill = Year), alpha = 0.6, width
= 1/4) + geom_stratum(alpha = 0.8, width = 1/4, fill = “grey95”, color =
“grey30”) + geom_text(stat = “stratum”, aes(label =
after_stat(stratum)), size = 3.5, fontface = “bold”) +
scale_fill_brewer(palette = “Set2”) + theme_minimal(base_size = 14) +
theme(panel.grid.major = element_blank(), panel.grid.minor =
element_blank(), axis.text.y = element_blank(), axis.title.y =
element_blank(), axis.text.x = element_text(face = “bold”)) + labs(title
= “Longitudinal flow of genotypic and phenotypic MDR”, subtitle =
dynamic_alluvial_subtitle, fill = “Year of sampling”)
ggsave(“Figure_7_Alluvial_MDR_Flow.svg”, plot = fig_7_alluvial_flow,
device = “svg”, units = “in”, width = 8.5, height = 4.8)

# ————————————————————————————————————————

# XIII. INDEPENDENT VALIDATION FRAMEWORK: STANDARD LOGISTIC REGRESSIONS (GLM)

# PARALLEL MATHEMATICAL AUDIT TO BENCHMARK PREVALENCE RISK DRIVERS

print(“RUNNING PARALLEL INTERFERENTIAL BENCHMARKING (STANDARD GLM)”)

# 1. Fit an un-nested univariable logistic regression by omitting the random intercept

model_validation_glm \<- glm( Result ~ Year_Fact + State + Source +
Season + pH + Temperature + Turbidity, data = base_muestras_anual,
family = binomial(link = “logit”) )

# 2. Extract standard errors, z-statistics, and exact Wald confidence intervals

summary_validation_glm \<- summary(model_validation_glm)\$coefficients
ci_validation_glm \<-
confint.default(model_validation_glm)\[rownames(summary_validation_glm),
\]

# 3. Construct the independent validation benchmarking dataframe

table_validation_glm \<- data.frame( Variable =
rownames(summary_validation_glm), Odds_Ratio =
round(exp(summary_validation_glm\[, “Estimate”\]), 4), CI_Lower =
round(exp(ci_validation_glm\[, 1\]), 4), CI_Upper =
round(exp(ci_validation_glm\[, 2\]), 4), P_Value =
round(summary_validation_glm\[, “Pr(\>\|z\|)”\], 4) ) %\>%
mutate(Significance = ifelse(P_Value \< 0.05, “SIGNIFICANT 🌟”,
“Stable”))

# 4. Print definitive metrics directly to the RStudio console for inspection

print(“PARALLEL AUDIT MATRIX: STANDARD GLM PREVALENCE COEFFICIENTS”)
print(as.data.frame(table_validation_glm))

# 5. Automated text export

sink(“Raw_Output_Table_2B_Standard_GLM_Prevalence_Validation.txt”)
cat(“PARALLEL AUDIT MATRIX: STANDARD GLM PREVALENCE DRIVERS”)
cat(“Executed systematically to validate fixed effects under zero random
variance thresholds”) print(summary(model_validation_glm)) sink()

print(“Parallel baseline validation complete. Output secured as
‘Raw_Output_Table_2B_Standard_GLM_Prevalence_Validation.txt’!”)

# ————————————————————————————————————————

# XIV. INDEPENDENT VALIDATION FRAMEWORK: STANDARD LOGISTIC REGRESSIONS (GLM)

# PARALLEL MATHEMATICAL AUDIT TO BENCHMARK MDR PHENOTYPES DRIVERS

# 1. Enforce strict discrete configurations at the isolate level to align benchmarks

base_isolates_audit \<- AMRg_veganb %\>% mutate( Year_Fact =
factor(Year, levels = c(“2019”, “2020”, “2021”, “2022”, “2023”)), State
= as.factor(State), Source = as.factor(Source), Season =
as.factor(Season), pH = as.numeric(as.character(pH)), Temperature =
as.numeric(as.character(Temperature)), Turbidity =
as.numeric(as.character(Turbidity)), MDRf =
as.numeric(as.character(MDRf)) ) %\>% filter(!is.na(Year_Fact) &
!is.na(State) & !is.na(Source) & !is.na(Season) & !is.na(MDRf) &
!is.na(pH) & !is.na(Temperature) & !is.na(Turbidity))

print(paste(“Total independent isolates entering parallel GLM audit for
MDRf:”, nrow(base_isolates_audit)))

# 2. Fit an un-nested standard logistic regression for Phenotypic MDR using the isolate matrix

model_validation_mdrf_glm \<- glm( MDRf ~ Year_Fact + State + Source +
Season + pH + Temperature + Turbidity, data = base_isolates_audit,
family = binomial(link = “logit”) )

# 3. Extract standard errors, z-statistics, and exact Wald confidence intervals

summary_validation_mdrf_glm \<-
summary(model_validation_mdrf_glm)\$coefficients ci_validation_mdrf_glm
\<-
confint.default(model_validation_mdrf_glm)\[rownames(summary_validation_mdrf_glm),
\]

# 4. Construct the independent validation benchmarking dataframe for MDRf

table_validation_mdrf_glm \<- data.frame( Variable =
rownames(summary_validation_mdrf_glm), Odds_Ratio =
round(exp(summary_validation_mdrf_glm\[, “Estimate”\]), 4), CI_Lower =
round(exp(ci_validation_mdrf_glm\[, 1\]), 4), CI_Upper =
round(exp(ci_validation_mdrf_glm\[, 2\]), 4), P_Value =
round(summary_validation_mdrf_glm\[, “Pr(\>\|z\|)”\], 4) ) %\>%
mutate(Significance = ifelse(P_Value \< 0.05, “SIGNIFICANT 🌟”,
“Stable”))

# 5. Print definitive metrics directly to the RStudio console for inspection

print(“PARALLEL AUDIT MATRIX: STANDARD GLM MDRf COEFFICIENTS”)
print(as.data.frame(table_validation_mdrf_glm))

# 6. Automated text export

sink(“Raw_Output_Table_2D_Standard_GLM_MDRf_Validation.txt”)
cat(“PARALLEL AUDIT MATRIX: STANDARD GLM MDRf DRIVERS”) cat(“Executed
systematically to validate isolate fixed effects under zero random
variance thresholds”) print(summary(model_validation_mdrf_glm)) sink()

print(“Parallel MDRf baseline validation complete. Output secured as
‘Raw_Output_Table_2D_Standard_GLM_MDRf_Validation.txt’!”)

print(“COMPLETE PIPELINE SUCCESSFULLY EXECUTED”)
