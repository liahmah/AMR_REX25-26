
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
