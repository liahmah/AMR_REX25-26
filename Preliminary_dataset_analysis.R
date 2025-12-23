# includes readr, ggplot2, dplyr packages
library(tidyverse)

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
  select(anon_id:order_proc_id_coded, organism:susceptibility)

# merging demographics with confirmed resistance
merged_demo_resis <- merge(data_demographics, data_resistance, 
                           by = c("pat_enc_csn_id_coded", "anon_id", 
                                  "order_proc_id_coded"))

# attempting to download comorbidity data in segments 
dataset_comorbidity <- file("datasets/microbiology_cultures_comorbidity.csv",
                            "r")
length(count.fields("datasets/microbiology_cultures_comorbidity.csv"))
# 206 547 141 rows incl header

comorbidity_s1 <- read.csv(dataset_comorbidity, nrows = 10000000, header = FALSE) |>
  filter(if_any(5, ~ .x == "Asthma")) # filtering for comorbidity = asthma
# no matching observations

  comorbidity_s2 <- read.csv(dataset_comorbidity, nrows = 10000000, skip = 10000000, 
                           header = FALSE) |> # analysing rows 10 000 001 - 20M
  filter(if_any(5, ~ .x == "Asthma")) # filtering for comorbidity = asthma
# 107693 matching observations

comorbidity_s3 <- read.csv(dataset_comorbidity, nrows = 10000000, skip = 20000000, 
                           header = FALSE) |> # analysing rows 20 000 001 - 30M
  filter(if_any(5, ~ .x == "Asthma")) # filtering for comorbidity = asthma
# 72481 matching observations

comorbidity_s4 <- read.csv(dataset_comorbidity, nrows = 10000000, skip = 30000000, 
                           header = FALSE) |> # analysing rows 30 000 001 - 40M
  filter(if_any(5, ~ .x == "Asthma")) # filtering for comorbidity = asthma
# 107657 matching observations

comorbidity_s5 <- read.csv(dataset_comorbidity, nrows = 10000000, skip = 40000000, 
                           header = FALSE) |> # analysing rows 40 000 001 - 50M
  filter(if_any(5, ~ .x == "Asthma")) # filtering for comorbidity = asthma
# 83976 matching observations

comorbidity_s6 <- read.csv(dataset_comorbidity, nrows = 10000000, skip = 50000000, 
                           header = FALSE) |> # analysing rows 50 000 001 - 60M
  filter(if_any(5, ~ .x == "Asthma")) # filtering for comorbidity = asthma
# 0 matching observations

comorbidity_s7 <- read.csv(dataset_comorbidity, nrows = 10000000, skip = 60000000, 
                           header = FALSE) |> # analysing rows 60 000 001 - 70M
  filter(if_any(5, ~ .x == "Asthma")) # filtering for comorbidity = asthma
s#  matching observations

comorbidity_s8 <- read.csv(dataset_comorbidity, nrows = 10000000, skip = 70000000, 
                           header = FALSE) |> # analysing rows 70 000 001 - 80M
  filter(if_any(5, ~ .x == "Asthma")) # filtering for comorbidity = asthma
#  matching observations

# searching for relationship age-resistance
table_age_resis <- table(merged_demo_resis$age, merged_demo_resis$antibiotic)
# doesn't appear at first glance to be strong relationship

# searching for relationship gender-resistance 
table_age_resis <- table(merged_demo_resis$gender, merged_demo_resis$antibiotic)
print(table_age_resis)

#