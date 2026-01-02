
# Setup ---------------------------------------------------------------

library(tidyverse)

merged_asthma_suscep <- read_csv("datasets/merged_asthma_susceptibility.csv")


# Further wrangling -------------------------------------------------------

# has a patient returned multiple times?
duplicated(merged_asthma_suscep$anon_id) # yes, anon_id identifies patient

## filtering, maybe resistance developed?
dupl_anon <- merged_asthma_suscep |>
  group_by(anon_id) |>
  filter(n() > 1)
### but repeated anon_id appears to be given to patients at same utc? but other
### variables differ so unclear... moving on

# converting null to NA
merged_asthma_suscep <- merged_asthma_suscep |>
  mutate(across(c(organism, antibiotic, susceptibility), ~ na_if(., "Null")))

# summarizing totals
summarized <- merged_asthma_suscep |>
  count(susceptibility)

# summarizing number of resistant/susceptible per organism by antibiotic
summarized_organism <- merged_asthma_suscep |>
  count(antibiotic, organism, susceptibility) |>
  rename('count' = 'n')

# or without organism:
summarized_antibiotic <- merged_asthma_suscep |>
  count(antibiotic, susceptibility) |>
  rename('count' = 'n')
  
# but can't compare to non-asthma patients. compare to who?

# could compare by antibiotic class:
max(summarized_antibiotic$count)
# Relationship analysis ---------------------------------------------------


