# Loading -----------------------------------------------------------------

library(tidyverse)
library(MASS)
library(emmeans)
library(viridis)

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
# wrangling
strep <- data_susceptibility |>
  filter(organism == "STREPTOCOCCUS PNEUMONIAE") |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() 
  # remains 6435 obs, 11 antibiotics, 764 patients, 2008-2024
  
  strep_n <- strep |>
    count(antibiotic, susceptibility, name = "count")

  strep_logical <- strep |>
    mutate(antibiotic = as.factor(antibiotic)) |>
    filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
    mutate(suscep_logical = susceptibility == "Susceptible")

# modelling
strep_bin_model <- glm(suscep_logical ~ antibiotic,
                       data = strep_logical,
                       family = binomial)
  emmeans(strep_bin_model, pairwise ~ antibiotic, type = "response")

# extracting results:
strep_emm <- emmeans(strep_bin_model, ~ antibiotic, type = "response")
strep_emm_results <- as.data.frame(strep_emm) 

# visualizing
strep_emm_results |>
  ggplot(aes(x = reorder(antibiotic, prob), y = prob)) +
  geom_col(fill = "royalblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(x = "Antibiotic",
       y = "Predicted probability of S. pneumoniae susceptibility",
       title = "Predicted S. pneumoniae Susceptibility by Antibiotic") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# STAPHYLOCOCCUS AUREUS ---------------------------------------------------

## Staph
# wrangling
staph <- data_susceptibility |>
  filter(organism == "STAPHYLOCOCCUS AUREUS") |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() 
  # remains 82513 obs, 18 antibiotics, > 4000 patients, 2008 - 2024
 
  staph_n <- staph |>
     count(antibiotic, susceptibility, name = "count")

staph_logical <- staph |>
  mutate(antibiotic = as.factor(antibiotic)) |>
  filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
  mutate(suscep_logical = susceptibility == "Susceptible")

# modelling
staph_bin_model <- glm(suscep_logical ~ antibiotic,
                       data = staph_logical,
                       family = binomial)
emmeans(staph_bin_model, pairwise ~ antibiotic, type = "response")

# extracting results:
staph_emm <- emmeans(staph_bin_model, ~ antibiotic, type = "response")
staph_emm_results <- as.data.frame(staph_emm) 

# visualizing
staph_emm_results |>
  ggplot(aes(x = reorder(antibiotic, prob), y = prob)) +
  geom_col(fill = "royalblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(x = "Antibiotic",
       y = "Predicted probability of S. aureus susceptibility",
       title = "Predicted S. aureus Susceptibility by Antibiotic") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## MRSA 
# wrangling
staph_mrsa <- data_susceptibility |>
  filter(organism == "STAPH AUREUS {MRSA}") |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() 
  # remains 16671 obs, 19 antibiotics, 823 patients, 2015 - 2024
  
  staph_mrsa_n <- staph_mrsa |>
    add_count(antibiotic, susceptibility, name = "count") 

staph_mrsa_logical <- staph_mrsa |>
  filter(susceptibility %in% c("Susceptible", "Resistant")) |>
  mutate(antibiotic = as.factor(antibiotic),
    suscep_logical = if_else(susceptibility == "Susceptible", TRUE, FALSE))

# modelling  
staph_mrsa_bin_model <- glm(suscep_logical ~ antibiotic,
                       data = staph_mrsa_logical,
                       family = binomial)
  emmeans(staph_mrsa_bin_model, pairwise ~ antibiotic, type = "response")

# extracting results:
staph_mrsa_emm <- emmeans(staph_mrsa_bin_model, ~ antibiotic, type = "response")
staph_mrsa_emm_results <- as.data.frame(staph_mrsa_emm)

# visualizing
staph_mrsa_emm_results |>
  ggplot(aes(x = reorder(antibiotic, prob), y = prob)) +
  geom_col(fill = "royalblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(x = "Antibiotic",
       y = "Predicted probability of MRSA susceptibility",
       title = "Predicted Methicillin-Resistant S. Aureus Susceptibility by Antibiotic") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## variant S. Aureus
# wrangling
staph_var <- data_susceptibility |>
  filter(organism == "STAPH AUREUS(COLONY VARIANT - SMALL COLONY OR OTHER MORPHOTYPE)") |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() 
  # remains 1837 obs, 7 antibiotics, 118 patients, 2015 - 2023

  staph_var_n <- staph_var |>
    count(antibiotic, susceptibility, name = "count") 
  
staph_var_logical <- staph_var |>
  mutate(antibiotic = as.factor(antibiotic)) |>
  filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
  mutate(suscep_logical = susceptibility == "Susceptible")

# modelling
staph_var_bin_model <- glm(suscep_logical ~ antibiotic,
                       data = staph_var_logical,
                       family = binomial)
  emmeans(staph_var_bin_model, pairwise ~ antibiotic, type = "response")

# extracting results:
staph_var_emm <- emmeans(staph_var_bin_model, ~ antibiotic, type = "response")
staph_var_emm_results <- as.data.frame(staph_var_emm)

# visualizing
staph_var_emm_results |>
  ggplot(aes(x = reorder(antibiotic, prob), y = prob)) +
  geom_col(fill = "royalblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(x = "Antibiotic",
       y = "Predicted probability of S. aureus variant susceptibility",
       title = "Predicted Susceptibility of S. aureus variant by Antibiotic") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# KLEBSIELLA PNEUMONIAE ---------------------------------------------------

# wrangling
kleb <- data_susceptibility |>
  filter(organism == "KLEBSIELLA PNEUMONIAE") |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() 
  # remains 151 859 obs, 28 antibiotics, > 7000 patients, 2008 - 2024

  kleb_n <- kleb |>
    count(antibiotic, susceptibility, name = "count") 
    
kleb_logical <- kleb |>
  mutate(antibiotic = as.factor(antibiotic)) |>
  filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
  mutate(suscep_logical = susceptibility == "Susceptible")

# modelling
  kleb_bin_model <- glm(suscep_logical ~ antibiotic,
                       data = kleb_logical,
                       family = binomial)
  emmeans(kleb_bin_model, pairwise ~ antibiotic, type = "response")

# extracting results:
  kleb_emm <- emmeans(kleb_bin_model, ~ antibiotic, type = "response")
  kleb_emm_results <- as.data.frame(kleb_emm) 

# visualizing
kleb_emm_results |>
  ggplot(aes(x = reorder(antibiotic, prob), y = prob)) +
  geom_col(fill = "royalblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(x = "Antibiotic",
       y = "Predicted probability of K. influenzae susceptibility",
       title = "Predicted K. influenzae Susceptibility by Antibiotic") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
## SSP Ozaenae

# wrangling
kleb_ozaenae <- data_susceptibility |>
  filter(organism == "KLEBSIELLA PNEUMONIAE SSP. OZAENAE")
  # only 34 obs, too small
  
## CP resistant SSP

# wrangling
kleb_cpresis <- data_susceptibility |>
  filter(organism == "KLEBSIELLA PNEUMONIAE (CARBAPENEM RESISTANT)") |>
  count(antibiotic)
  # 373 obs but no more than 20 obs / antibiotic => too small

# PSEUDOMONAS AERUGINOSA --------------------------------------------------

# wrangling
pseudoaeru <- data_susceptibility |>
  filter(organism == "PSEUDOMONAS AERUGINOSA") |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() 
  # remains 48847 obs, 16 antibiotics, > 3000 patients, 2008 - 2024
 
  pseudoaeru_n <- pseudoaeru |>
    count(antibiotic, susceptibility, name = "count") 
    
pseudoaeru_logical <- pseudoaeru |>
  mutate(antibiotic = as.factor(antibiotic)) |>
  filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
  mutate(suscep_logical = susceptibility == "Susceptible")

# modelling
  pseudoaeru_bin_model <- glm(suscep_logical ~ antibiotic,
                       data = pseudoaeru_logical,
                       family = binomial)
  emmeans(pseudoaeru_bin_model, pairwise ~ antibiotic, type = "response")

# extracting results:
  pseudoaeru_emm <- emmeans(pseudoaeru_bin_model, ~ antibiotic, type = "response")
  pseudoaeru_emm_results <- as.data.frame(pseudoaeru_emm) 

# visualizing
pseudoaeru_emm_results |>
  ggplot(aes(x = reorder(antibiotic, prob), y = prob)) +
  geom_col(fill = "royalblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(x = "Antibiotic",
       y = "Predicted probability of P. aeruginosa susceptibility",
       title = "Predicted P. aeruginosa Susceptibility by Antibiotic") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## mucoid
# wrangling
pseudoaeru_muc <- data_susceptibility |>
  filter(organism == "MUCOID PSEUDOMONAS AERUGINOSA") |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() 
  # remains 32894 obs, 15 antibiotics, 823 patients, 2008 - 2024
 
  pseudoaeru_muc_n <- pseudoaeru_muc |>
    count(antibiotic, susceptibility, name = "count") 
  
pseudoaeru_muc_logical <- pseudoaeru_muc |>
  mutate(antibiotic = as.factor(antibiotic)) |>
  filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
  mutate(suscep_logical = susceptibility == "Susceptible")

# modelling
  pseudoaeru_muc_bin_model <- glm(suscep_logical ~ antibiotic,
                                data = pseudoaeru_muc_logical,
                                family = binomial)
  emmeans(pseudoaeru_muc_bin_model, pairwise ~ antibiotic, type = "response")
  
# extracting results:
  pseudoaeru_muc_emm <- emmeans(pseudoaeru_muc_bin_model, ~ antibiotic, 
                                type = "response")
  pseudoaeru_muc_emm_results <- as.data.frame(pseudoaeru_muc_emm) 

# visualizing
pseudoaeru_muc_emm_results |>
  ggplot(aes(x = reorder(antibiotic, prob), y = prob)) +
  geom_col(fill = "royalblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(x = "Antibiotic",
       y = "Predicted probability of mucoid P. aeruginosa susceptibility",
       title = "Predicted mucoid P. aeruginosa Susceptibility by Antibiotic") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## nonmucoid
# wrangling
pseudoaeru_nonmuc  <- data_susceptibility |>
  filter(organism == "PSEUDOMONAS AERUGINOSA (NON-MUCOID CF)") |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() 
  # remains 16252 obs, 15 antibiotics, 315 patients, 2013 - 2023

  pseudoaeru_nonmuc_n <- pseudoaeru_nonmuc |>
    count(antibiotic, susceptibility, name = "count") 
    
  
pseudoaeru_nonmuc_logical <- pseudoaeru_nonmuc |>
  mutate(antibiotic = as.factor(antibiotic)) |>
  filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
  mutate(suscep_logical = susceptibility == "Susceptible")

# modelling
  pseudoaeru_nonmuc_bin_model <- glm(suscep_logical ~ antibiotic,
                                data = pseudoaeru_nonmuc_logical,
                                family = binomial)
  emmeans(pseudoaeru_nonmuc_bin_model, pairwise ~ antibiotic, type = "response")
  
# extracting results:
  pseudoaeru_nonmuc_emm <- emmeans(pseudoaeru_nonmuc_bin_model, ~ antibiotic, 
                                   type = "response")
  pseudoaeru_nonmuc_emm_results <- as.data.frame(pseudoaeru_nonmuc_emm) 

# visualizing
pseudoaeru_nonmuc_emm_results |>
  ggplot(aes(x = reorder(antibiotic, prob), y = prob)) +
  geom_col(fill = "royalblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(x = "Antibiotic",
       y = "Predicted probability of non-mucoid P. aeruginosa susceptibility",
       title = "Predicted non-mucoid P. aeruginosa Susceptibility by Antibiotic") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
# ESCHERICHIA COLI --------------------------------------------------------

# wrangling
ecoli <- data_susceptibility |>
  filter(organism == "ESCHERICHIA COLI") |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() 
  # remains 874 699 obs, 30 antibiotics, > 39 000 patients, 2007 - 2024
  
  ecoli_n <- ecoli |>
    count(antibiotic, susceptibility, name = "count")

  ecoli_logical <- ecoli |>
    mutate(antibiotic = as.factor(antibiotic)) |>
    filter(susceptibility == "Susceptible" | susceptibility == "Resistant") |>
    mutate(suscep_logical = susceptibility == "Susceptible")

# modelling
  ecoli_bin_model <- glm(suscep_logical ~ antibiotic,
                       data = ecoli_logical,
                       family = binomial)
  emmeans(ecoli_bin_model, pairwise ~ antibiotic, type = "response")

# extracting results:
  ecoli_emm <- emmeans(ecoli_bin_model, ~ antibiotic, type = "response")
  ecoli_emm_results <- as.data.frame(ecoli_emm) 

# visualizing
ecoli_emm_results |>
  ggplot(aes(x = reorder(antibiotic, prob), y = prob)) +
  geom_col(fill = "royalblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(x = "Antibiotic",
       y = "Predicted probability of K. influenzae susceptibility",
       title = "Predicted K. influenzae Susceptibility by Antibiotic") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
## Carbapenem resistant strain
# wrangling
ecoli_cpresis  <- data_susceptibility |>
  filter(organism == "ESCHERICHIA COLI (CARBAPENEM RESISTANT)") |>
  count(antibiotic)
  # only 329 obs, no more than 19 / antibiotic => too small sample size

# noodling ----------------------------------------------------------------

# count / antibiotic facet plot

data_organisms <- data_organisms |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() |>
  group_by(organism) |>
  filter(n() > 400) |>
  ungroup() 

data_organisms |>
  ggplot(aes(x = antibiotic, fill = susceptibility)) +
  geom_bar() +
  scale_y_log10() +
  facet_wrap(~organism) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4))

# tying it all together

all_dfs <- list(ecoli_emm_results, kleb_emm_results, pseudoaeru_emm_results, 
                pseudoaeru_muc_emm_results, pseudoaeru_nonmuc_emm_results,
                staph_emm_results,  staph_mrsa_emm_results, staph_var_emm_results,
                strep_emm_results)

labels <- c("E. coli",  "K. pneumoniae", "P. aeruginosa", 
            "Mucosal P. aeruginosa",  "Non-mucosal P. aeruginosa",
            "S. aureus",  "MRSA", "S. aureus variant",
            "S. pneumoniae")


all_dfs <- mapply(function(df, lab) {
  
  # Convert to data frame if not already
  if(!is.data.frame(df)) df <- as.data.frame(df)
  
  # Add dataset label as character
  df$dataset <- lab
  
  return(df)
}, all_dfs, labels, SIMPLIFY = FALSE)

# --- 3. Standardize columns: fill missing columns with NA ---
# Find all unique column names
all_cols <- unique(unlist(lapply(all_dfs, colnames)))

# Make sure each df has all columns
all_dfs <- lapply(all_dfs, function(df) {
  missing <- setdiff(all_cols, colnames(df))
  for(col in missing) df[[col]] <- NA
  # reorder columns consistently
  df[, all_cols]
})

# --- 4. Combine all data frames ---
all_combined <- do.call(rbind, all_dfs)

# --- 5. Optional: ensure dataset is a factor for faceting ---
all_combined$dataset <- factor(all_combined$dataset, levels = labels)

all_combined <- all_combined |>
  group_by(dataset) |>
  mutate(rank_prob = rank(-prob, ties.method = "first"),
         is_top2 = rank_prob <= 2) |>
  ungroup()

all_combined |>
  ggplot(aes(x = antibiotic, y = prob, fill = is_top2)) +
  geom_col() +
  facet_wrap(~ dataset, scales = "free_x") + # one panel per dataset
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  scale_fill_manual(values = c("TRUE"= "firebrick", "FALSE" = "royalblue")) +
  theme_minimal() +
  labs(y = "Predicted Probability", x = "Antibiotic",
       title = "Probability of antibiotic susceptibility by organism") +
  theme(axis.text.x = element_text(angle = 45, size = 7, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "none")
        
          