
# Set-up ------------------------------------------------------------------

library(tidyverse)
library(ggsignif)
library(MASS)
library(emmeans)

pneumonia_data <- read.csv("datasets/asthmapneum_dataset.csv")

# to be extra certain that susceptibility is retained, coding "null" or missing
# data as resistant.
extra_resis <- pneumonia_data |>
  mutate(susceptibility = recode(susceptibility,
                                 "Null" = "Resistant"))|>
  group_by(susceptibility) |>
  

# Quantifying Relationships -----------------------------------------------

proport_suscep <- pneumonia_data |>
  filter(category != "Null") %>%
  group_by(antibiotic) %>%
  summarise(
    susceptible = sum(count[category == "Susceptible"]),
    total_tested = sum(count),
    susceptibility_rate = susceptible / total_tested
  ) %>%
  arrange(desc(susceptibility_rate))

summary_table
