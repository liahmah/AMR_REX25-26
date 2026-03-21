
# Loading... --------------------------------------------------------------

library(tidyverse)
library(ggsignif)
library(dunn.test)
library(MASS)

data_suscep <- read_csv("datasets/microbiology_cultures_cohort.csv")

# looking at demographics
data_demo <- read_csv("datasets/microbiology_cultures_demographics.csv")

merged_demo_suscep <- merge(data_demo, data_suscep, 
                            by = c("pat_enc_csn_id_coded", "anon_id", 
                                   "order_proc_id_coded"))


# looking at leukocyte values
data_leuko_median <- read_csv("datasets/microbiology_cultures_labs.csv") |>
  dplyr::select(anon_id:Period_Day, contains("median")) |>
  mutate(across(where(is.character), ~ na_if(., "Null")))

merged_leuko_suscep <- merge(data_leuko_median, data_suscep, 
                              by = c("pat_enc_csn_id_coded", "anon_id", 
                                     "order_proc_id_coded"))

# Age ---------------------------------------------------------------------

# both categorical because age is binned. using ordinal logistic model
merged_demo_suscep$age <- factor(merged_demo_suscep$age,
                         levels = c("25-34 years","35-44 years","45-54 years",
                                    "55-64 years","65-74 years","75-84 years",
                                    "above 90 years"),
                         ordered = TRUE) |>
  as.numeric()

merged_demo_suscep$susceptibility <- as.factor(merged_demo_suscep$susceptibility)

# running model
age_olm_model <- polr(susceptibility ~ age, data = merged_demo_suscep, 
                      Hess = TRUE)
summary(age_olm_model)

# retrieving p-value results
coef_age_model <- coef(summary(age_olm_model))
age_p <- pnorm(abs(coef_age_model[, "t value"]), lower.tail = FALSE) * 2
age_p


# Gender ------------------------------------------------------------------



# Leukocytes --------------------------------------------------------------

## NEUTROPHILS
# make numeric
merged_leuko_suscep$median_neutrophils <- as.numeric(merged_leuko_suscep$median_neutrophils) |>
# sample size too large for shapiro test, so visualize distribution:
hist(na.omit(merged_leuko_suscep$median_neutrophils)) # heavily skewed, non-parametric

# kruskal-wallis test for non-parametric data, continuous to >2 categorical
kruskal.test(merged_leuko_suscep$susceptibility ~ merged_leuko_suscep$median_neutrophils)
# dunn test for pair-wise comparisons tells us WHICH pairs are significant
dunn.test(merged_leuko_suscep$median_neutrophils, merged_leuko_suscep$susceptibility,
          kw = TRUE)

## not WBC
merged_leuko_suscep$median_wbc <- as.numeric(merged_leuko_suscep$median_wbc)
hist(na.omit(merged_leuko_suscep$median_wbc)) # only one bar

## hemoglobin
merged_leuko_suscep$median_hgb <- as.numeric(merged_leuko_suscep$median_hgb)
hist(na.omit(merged_leuko_suscep$median_hgb))

# ANOVA test for parametric data, continuous to >2 categorical
aov_hgb <- summary(aov(lm(susceptibility ~ median_hgb, data = na.omit(merged_leuko_suscep))))

TukeyHSD(aov(lm(susceptibility ~ median_hgb, data = na.omit(merged_leuko_suscep))))


# kruskal-wallis test for non-parametric data, continuous to >2 categorical
kruskal.test(merged_leuko_suscep$susceptibility ~ merged_leuko_suscep$median_neutrophils)
# dunn test for pair-wise comparisons tells us WHICH pairs are significant
dunn.test(merged_leuko_suscep$median_neutrophils, merged_leuko_suscep$susceptibility,
          kw = TRUE)



