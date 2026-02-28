
# Loading -----------------------------------------------------------------

library(tidyverse)
library(MASS)
library(emmeans)

data_susceptibility <- read_csv("datasets/microbiology_cultures_cohort.csv") |>
  mutate(across(c("pat_enc_csn_id_coded", "order_proc_id_coded"), as.double))

data_demographics <- read_csv("datasets/microbiology_cultures_demographics.csv")

demo_suscep <- merge(data_demographics, data_susceptibility,
                     by = c("order_proc_id_coded", "pat_enc_csn_id_coded", 
                            "anon_id")) |>
  filter(organism != "Null")

## selecting pneumonia-causing species of interest:

species <- as.data.frame(unique(demo_suscep$organism))
# Filtering for Streptococcus pneumoniae, yes
# Haemophilus influenzae, no
# Klebsiella pneumoniae, yes, 3 subvariants
# Pseudomonas aeruginosa, yes, 3 subvariants
# Escherichia coli, yes, 2 subvariants
# Staphylococcus aureus, yes, 3 subvariants

# excluding unidentified species
# susceptibility coded as true (susceptible) or false (resistant)

data_organisms <- demo_suscep |>
  filter(organism %in% c("STREPTOCOCCUS PNEUMONIAE", "STAPHYLOCOCCUS AUREUS",
                       "STAPH AUREUS(COLONY VARIANT - SMALL COLONY OR OTHER MORPHOTYPE)",
                       "STAPH AUREUS {MRSA}", "KLEBSIELLA PNEUMONIAE", 
                       "KLEBSIELLA PNEUMONIAE SSP. OZAENAE", 
                       "KLEBSIELLA PNEUMONIAE (CARBAPENEM RESISTANT)", 
                       "PSEUDOMONAS AERUGINOSA", "MUCOID PSEUDOMONAS AERUGINOSA",
                       "PSEUDOMONAS AERUGINOSA (NON-MUCOID CF)", "ESCHERICHIA COLI",
                       "ESCHERICHIA COLI (CARBAPENEM RESISTANT)"))

count_suscep_organisms <- data_organisms |>
  count(antibiotic, susceptibility, organism)


# STREPTOCOCCUS PNEUMONIAE ------------------------------------------------

strep <- data_susceptibility |>
  filter(organism == "STREPTOCOCCUS PNEUMONIAE")
# 6486 samples, 764 unique patients (multiple antibiotics tested / patient), 16 antibiotics

  strep_n <- strep |>
    mutate(antibiotic = as.factor(antibiotic))

  strep_logical <- strep_n |>
    filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
    mutate(suscep_logical = susceptibility == "Susceptible")


# STAPHYLOCOCCUS AUREUS ---------------------------------------------------

staph <- data_susceptibility |>
  filter(organism == "STAPHYLOCOCCUS AUREUS")
  # 82519 samples, 23 antibiotics
  
  staph_n <- staph |>
    add_count(antibiotic, susceptibility, name = "count") |>
    mutate(antibiotic = as.factor(antibiotic))
  
  staph_logical <- staph_n |>
    filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
    mutate(suscep_logical = susceptibility == "Susceptible")

# MRSA STAPH AUREUS

staph_mrsa <- data_susceptibility |>
  filter(organism == "STAPH AUREUS {MRSA}") 
  # 16674 samples, 823 patients, 19 antibiotics (1 distinct)
  
  staph_mrsa_n <- staph_mrsa |>
    add_count(antibiotic, susceptibility, name = "count") |>
    mutate(antibiotic = as.factor(antibiotic))

  staph_mrsa_logical <- staph_mrsa_n |>
    filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
    mutate(suscep_logical = susceptibility == "Susceptible")

# variant S. Aureus
staph_var <- data_susceptibility |>
    filter(organism == "STAPH AUREUS(COLONY VARIANT - SMALL COLONY OR OTHER MORPHOTYPE)") 
  # 1870 samples,

  staph_var_n <- staph_var |>
    add_count(antibiotic, susceptibility, name = "count") |>
    mutate(antibiotic = as.factor(antibiotic))
  
  staph_var_logical <- staph_var_n |>
    filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
    mutate(suscep_logical = susceptibility == "Susceptible")


# KLEBSIELLA PNEUMONIAE ---------------------------------------------------
kleb <- data_susceptibility |>
    filter(organism == "KLEBSIELLA PNEUMONIAE") 
  # 151 964 samples,

  kleb_n <- kleb |>
    add_count(antibiotic, susceptibility, name = "count") |>
    mutate(antibiotic = as.factor(antibiotic))
  
  kleb_logical <- kleb_n |>
    filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
    mutate(suscep_logical = susceptibility == "Susceptible")

kleb_ozaenae <- data_susceptibility |>
  filter(organism == "KLEBSIELLA PNEUMONIAE SSP. OZAENAE") 
  # 34 samples,

  kleb_ozaenae_n <- kleb_ozaenae |>
    add_count(antibiotic, susceptibility, name = "count") |>
    mutate(antibiotic = as.factor(antibiotic))
  
  kleb_ozaenae_logical <- kleb_ozaenae_n |>
    filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
    mutate(suscep_logical = susceptibility == "Susceptible")
  
kleb_cpresis <- data_susceptibility |>
  filter(organism == "KLEBSIELLA PNEUMONIAE (CARBAPENEM RESISTANT)") 
  # 373 samples,

  kleb_cpresis_n <- kleb_cpresis |>
    add_count(antibiotic, susceptibility, name = "count") |>
    mutate(antibiotic = as.factor(antibiotic))
  
  kleb_cpresis_logical <- kleb_cpresis_n |>
    filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
    mutate(suscep_logical = susceptibility == "Susceptible")
  
# PSEUDOMONAS AERUGINOSA --------------------------------------------------
pseudoaeru <- data_susceptibility |>
    filter(organism == "PSEUDOMONAS AERUGINOSA") 
  # 48959 samples,
  
  pseudoaeru_n <- pseudoaeru |>
    add_count(antibiotic, susceptibility, name = "count") |>
    mutate(antibiotic = as.factor(antibiotic))
  
  pseudoaeru_logical <- pseudoaeru_n |>
    filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
    mutate(suscep_logical = susceptibility == "Susceptible")

pseudoaeru_muc <- data_susceptibility |>
  filter(organism == "MUCOID PSEUDOMONAS AERUGINOSA") 
  # 33037 samples,
  
  pseudoaeru_muc_n <- pseudoaeru_muc |>
    add_count(antibiotic, susceptibility, name = "count") |>
    mutate(antibiotic = as.factor(antibiotic))
  
  pseudoaeru_muc_logical <- pseudoaeru_muc_n |>
    filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
    mutate(suscep_logical = susceptibility == "Susceptible")

pseudoaeru_nonmuc  <- data_susceptibility |>
    filter(organism == "PSEUDOMONAS AERUGINOSA (NON-MUCOID CF)") 
  # 16327 samples,
  
  pseudoaeru_nonmuc_n <- pseudoaeru_nonmuc |>
    add_count(antibiotic, susceptibility, name = "count") |>
    mutate(antibiotic = as.factor(antibiotic))
  
  pseudoaeru_nonmuc_logical <- pseudoaeru_nonmuc_n |>
    filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
    mutate(suscep_logical = susceptibility == "Susceptible")
  
# ESCHERICHIA COLI --------------------------------------------------------
  
ecoli <- data_susceptibility |>
    filter(organism == "ESCHERICHIA COLI") 
  # 874 897 samples,
  
  ecoli_n <- ecoli |>
    add_count(antibiotic, susceptibility, name = "count") |>
    mutate(antibiotic = as.factor(antibiotic))
  
  ecoli_logical <- ecoli_n |>
    filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
    mutate(suscep_logical = susceptibility == "Susceptible")
  
ecoli_cpresis  <- data_susceptibility |>
    filter(organism == "ESCHERICHIA COLI (CARBAPENEM RESISTANT)") 
  # 329 samples,
  
  ecoli_cpresis_n <- ecoli_cpresis |>
    add_count(antibiotic, susceptibility, name = "count") |>
    mutate(antibiotic = as.factor(antibiotic))
  
  ecoli_cpresis_logical <- ecoli_cpresis_n |>
    filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
    mutate(suscep_logical = susceptibility == "Susceptible")

# odds --------------------------------------------------------------------

## binomial logistic regression
bin_model <- glm(suscep_logical ~ order_time_jittered_utc + antibiotic, 
             data = staph_aur_logical, 
             family = binomial)
summary(bin_model)

## STAPH
staph_bin_model <- glm(suscep_logical ~ antibiotic,
                       data = staph_aur_logical,
                       family = binomial)
emmeans(staph_bin_model, pairwise ~ antibiotic, type = "response")

# extracting results:
staph_emm <- emmeans(staph_bin_model, ~ antibiotic, type = "response")
staph_emm_results <- as.data.frame(staph_emm) 

## MRSA STAPH
staph_mrsa_bin_model <- glm(suscep_logical ~ antibiotic,
                       data = staph_aur_mrsa_logical,
                       family = binomial)
emmeans(staph_mrsa_bin_model, pairwise ~ antibiotic, type = "response")

# extracting results:
staph_mrsa_emm <- emmeans(staph_mrsa_bin_model, ~ antibiotic, type = "response")
staph_mrsa_emm_results <- as.data.frame(staph_mrsa_emm)

## STREP

strep_bin_model <- glm(suscep_logical ~ antibiotic,
                       data = strep_pneum_logical,
                       family = binomial)
emmeans(strep_bin_model, pairwise ~ antibiotic, type = "response")

# extracting results:
strep_emm <- emmeans(strep_bin_model, ~ antibiotic, type = "response")
strep_emm_results <- as.data.frame(strep_emm)  


# significance for Cefuroxime, Erythromycin, 

# visualizations ----------------------------------------------------------
strep_emm_results |>
  ggplot(aes(x = reorder(antibiotic, prob), y = prob)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  coord_flip() +  # flip axes for easier reading
  labs(
    x = "Antibiotic",
    y = "Predicted probability of S. pneumoniae susceptibility",
    title = "Predicted S. pneumoniae Susceptibility by Antibiotic"
  ) +
  theme_minimal()

staph_emm_results |>
  ggplot(aes(x = reorder(antibiotic, prob), y = prob)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  coord_flip() +  # flip axes for easier reading
  labs(
    x = "Antibiotic",
    y = "Predicted probability of S. aureus susceptibility",
    title = "Predicted S. aureus Susceptibility by Antibiotic"
  ) +
  theme_minimal()

staph_mrsa_emm_results |>
  ggplot(aes(x = reorder(antibiotic, prob), y = prob)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  coord_flip() +  # flip axes for easier reading
  labs(
    x = "Antibiotic",
    y = "Predicted probability of MRSA susceptibility",
    title = "Predicted Methicillin-Resistant S. Aureus Susceptibility by Antibiotic"
  ) +
  theme_minimal()







# noodling ----------------------------------------------------------------

table(strep_n$antibiotic, strep_n$susceptibility)
# dropping Cloramphenicol, daptomycin, ertapenem, tetracycline for sample sizes

staph_aur_n |>
  ggplot(aes(x = susceptibility, y = count)) +
  geom_point() +
  facet_wrap(facets = vars(antibiotic)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

staph_aur_mrsa_n |>
  ggplot(aes(x = susceptibility, y = count)) +
  geom_point() +
  facet_wrap(facets = vars(antibiotic)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

demo_suscep |>
  ggplot(aes(x = antibiotic, y = organism, colour = susceptibility)) +
  geom
