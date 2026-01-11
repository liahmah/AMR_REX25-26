
# libraries and datasets --------------------------------------------------

library(tidyverse)
library(MASS)
library(emmeans)
library(brant)

asthma_dataset <- read_csv("datasets/merged_asthma_susceptibility.csv")

# wrangling ---------------------------------------------------------------

pneumonia <- asthma_dataset |>
  count(antibiotic, organism, susceptibility) |> # counts num obersvations per
  rename('count' = 'n') |>
  filter(organism == "STREPTOCOCCUS PNEUMONIAE")

pneumonia$antibiotic <- as.factor(pneumonia$antibiotic)
# changing data type for levelling

# Stats analysis ----------------------------------------------------------

# set ref level as penicillin, considering of resistance concern and largest sample 
pneumonia$antibiotic <- relevel(pneumonia$antibiotic, ref = "Penicillin")

# ordinal logistic regression accounts for intermediate susceptibility
## formatting for order susceptible < intermediate < resistant
pneumonia$susceptibility <- factor(pneumonia$susceptibility,
  levels = c("Susceptible", "Intermediate", "Resistant"),
  ordered = TRUE)

## modelling (proportional odds model)
olr_model <- polr(susceptibility ~ antibiotic, data = pneumonia, Hess = TRUE)
summary(olr_model) # neg coefficients lean towards susceptibilty than resistance
### antibiotics levofloxacin, linezolid, vancomycin most negative

## converting to odds ratio
exp(coef(olr_model)) # <1 signifies better susceptibility than penicillin
### linezolid, vancomycin, levofloxacin best by 7 degrees (e-08 vs e-01). confirms 

## tukey pairwise to penicillin
emmeans(olr_model, pairwise ~ antibiotic, adjust = "tukey",
        type = "response")
### confirms. estimate ~ 18 for all 3 relative to ~ 0 for all others
### for those 3 only, p.value < 0.0001


## confirming significance
brant(olr_model) # some sparse outcomes, limiting the brant test, but does confirm

## can also confirm by binary collapse, so considering intermediate as resistant
pneumonia$non_susceptible <-  pneumonia$susceptibility %in% 
  c("Intermediate", "Resistant")

pneumonia$antibiotic <- relevel(pneumonia$antibiotic,
  ref = "Penicillin")

bin_model <- glm(non_susceptible ~ antibiotic,
  family = binomial,
  data = pneumonia)
summary(bin_model)

exp(coef(bin_model)) # only <1 intercepts (meaning less resistant) are same 3

emmeans(bin_model, pairwise ~ antibiotic, adjust = "tukey",
        type = "response") #confirms

## in both emmeans and glm, varying sample sizes are accounted for

# Discussion --------------------------------------------------------------

# don't trust linezolid, too small sample size
subset(pneumonia, antibiotic == "Linezolid") # only 61 observations

# checking vancomycin and levofloxacin:
subset(pneumonia, antibiotic == "Vancomycin") # 1973 observations
subset(pneumonia, antibiotic == "Levofloxacin") # 1026 observations

# could small sample size of linezolid be affecting analysis? try without
large_pneumonia <- pneumonia |>
  filter(antibiotic != "Linezolid")

large_pneumonia$antibiotic <- relevel(large_pneumonia$antibiotic, 
                                      ref = "Penicillin")

large_pneumonia$susceptibility <- factor(large_pneumonia$susceptibility,
                                   levels = c("Susceptible", "Intermediate", "Resistant"),
                                   ordered = TRUE)

large_model <- polr(susceptibility ~ antibiotic, data = large_pneumonia, Hess = TRUE)
summary(large_model) # confirms levofloxacin and vancomycin, also slightly meropenem

exp(coef(large_model)) # very small e-08 levofloxacin and vancomycin only

emmeans(large_model, pairwise ~ antibiotic, adjust = "tukey",
        type = "response") # confirms levofloxacin and vancomycin only
  ### p < 0.001 for vancomycin and levofloxacin

# demographic stuff
distinct(asthma_dataset, anon_id) # 32,672 patients
distinct(pneumonia, antibiotic) # 12

# Conclusion --------------------------------------------------------------

## Levofloxacin and Vancomycin have significantly better odds of susceptibility
## against S. pneumoniae than penicillin, as compared to 10 other antibiotics
## (in asthmatic individuals)

