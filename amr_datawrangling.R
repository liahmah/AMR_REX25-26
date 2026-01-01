
library(tidyverse) # includes readr, ggplot2, dplyr packages
library(arrow)  # for lazy CSV reading without loading into memory


# Loading/wrangling potentially useful datasets ---------------------------


# dataset includes patient data, infecting organism w/ confirmed resistance, 
# antibiotic prescribed, removing columns related to sample processing 
data_resistance <- read_csv("datasets/microbiology_cultures_microbial_resistance.csv") |> 
  select(-order_time_jittered_utc, -resistant_time_to_culturetime)

# dataset includes patient data, age and gender demographics
# smaller dataset, no wrangling needed, may limit sample size
data_demographics <- read_csv("datasets/microbiology_cultures_demographics.csv")

# dataset includes blood composition information, selecting for summary values only
data_leukocytes <- read_csv("datasets/microbiology_cultures_labs.csv") |>
  select(anon_id:median_procalcitonin)

# could be useful for investigating if prior exposure predisposes to resistance
data_priorexposure <- read_csv("datasets/microbiology_cultures_antibiotic_class_exposure.csv")
# another subtype dataset exists if interested in that route

# dataset includes patient data, prior medical conditions, removing columns related
# to sample processing
data_priorinfection <- read_csv("datasets/microbiology_cultures_prior_med.csv") |>
  select(-order_time_jittered_utc, -medication_time_to_culturetime)

# dataset includes susceptibility or resistance of different organisms to specified 
# antibiotic, removing columns relating to sample method
data_susceptibility <- read_csv("datasets/microbiology_cultures_cohort.csv") |>
  select(anon_id:order_proc_id_coded, organism:susceptibility) |>
  mutate(across(c("pat_enc_snc_id_coded", "order_proc_id_coded"), as.double))

## loading comorbidity dataset, too large so filtering for asthma

length(count.fields("datasets/microbiology_cultures_comorbidity.csv"))
# 206 547 141 rows incl header

# Using Arrow for lazy CSV reading - only loads matching rows into memory
# Open the dataset once (doesn't load into memory)
comorbidity_lazy <- open_csv_dataset("datasets/microbiology_cultures_comorbidity.csv")

# Filter for Asthma - uses same dplyr syntax, but only loads matching rows
data_comorbidity_asthma <- comorbidity_lazy |>
  filter(comorbidity_component == "Asthma") |>
  collect()
# ~370k matching observations total


# Exploring Relationships -------------------------------------------------



# merging demographics with confirmed resistance
merged_demo_resis <- merge(data_demographics, data_resistance, 
                           by = c("pat_enc_csn_id_coded", "anon_id", 
                                  "order_proc_id_coded"))

# searching for relationship age-resistance
table_age_resis <- table(merged_demo_resis$age, merged_demo_resis$antibiotic)
# doesn't appear at first glance to be strong relationship

# searching for relationship gender-resistance 
table_age_resis <- table(merged_demo_resis$gender, merged_demo_resis$antibiotic)
print(table_age_resis)

## asthma x antibiotic susceptibility

# merging
merged_asthma_suscep <- merge(data_comorbidity_asthma, data_susceptibility, 
                              by = c("pat_enc_csn_id_coded", "anon_id", 
                                     "order_proc_id_coded"))
