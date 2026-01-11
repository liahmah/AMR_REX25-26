
# libraries and datasets --------------------------------------------------

library(tidyverse)
library(MASS)
library(emmeans)

asthma_dataset <- read_csv("datasets/merged_asthma_susceptibility.csv")

# wrangling ---------------------------------------------------------------

pneumonia <- asthma_dataset |>
  count(antibiotic, organism, susceptibility) |>
  rename('count' = 'n') |>
  filter(organism == "STREPTOCOCCUS PNEUMONIAE") |>
  pivot_wider(names_from = susceptibility,
              values_from = count)

pneumonia$antibiotic <- as.factor(pneumonia$antibiotic)

# Stats analysis ----------------------------------------------------------

# can strep pneumoniae cause other illnesses than pneumonia?

# set levels as ___, most commonly prescribed antibiotic for pneumonia
# sum_anti_wide$antibiotic <- relevel(sum_anti_wide$antibiotic, ref = "Ciprofloxacin")

# then, run binomial regression
pneum_model <- glm(cbind(Resistant, Susceptible) ~ antibiotic,
                      family = binomial,
                      data = pneumonia)
summary(pneum_model)

# do confidence intervals and coefficients if desired

# confirm via Tukey
tukey_pneum <- glm(cbind(Resistant, Susceptible) ~ antibiotic,
                   family = binomial, 
                   data = pneumonia)

pairs(emmeans(tukey_pneum, ~ antibiotic), adjust = "tukey")

