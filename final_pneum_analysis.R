
# Set-up ------------------------------------------------------------------

library(tidyverse)
library(ggsignif)
library(MASS)
library(emmeans)

pneumonia_data <- read.csv("datasets/asthmapneum_dataset.csv")

# to be extra certain that susceptibility is retained, coding "null" or missing
# data as resistant.
extra_resis <- pneumonia_data |>
  mutate(susceptibility = if_else(susceptibility == "Null",
                              "Resistant", susceptibility)) |>
  group_by(antibiotic, susceptibility) |>
  summarise(count = sum(count), .groups = "drop")

# set penicillin as reference value for odds ratios, considering its large sample 
# size, concern regarding its resistance, and its common use to treat pneumonia
extra_resis$antibiotic <- as.factor(extra_resis$antibiotic)

extra_resis$antibiotic <- relevel(extra_resis$antibiotic, ref = "Penicillin")

# Quantifying Relationships -----------------------------------------------

#### rate of susceptibility:
proport_suscep <- extra_resis |>
  group_by(antibiotic) |>
  summarise(
    susceptible = sum(count[susceptibility == "Susceptible"]),
    total_tested = sum(count),
    susceptibility_rate = susceptible / total_tested) |>
  arrange(desc(susceptibility_rate))

#### ordinal logistic regression accounts for intermediate susceptibility
## formatting for order susceptible < intermediate < resistant
extra_resis$susceptibility <- factor(extra_resis$susceptibility,
                                   levels = c("Susceptible", "Intermediate", "Resistant"),
                                   ordered = TRUE)

# modelling (proportional odds model:
# odds of antibiotic being in a higher susceptibility category
olr_model <- polr(susceptibility ~ antibiotic, data = extra_resis, Hess = TRUE)

summary(olr_model) # neg coefficients indicate increasing odds of susceptibility
## antibiotics levofloxacin, linezolid, vancomcin

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

