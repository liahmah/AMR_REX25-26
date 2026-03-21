## pairwise comparisons
comparisons_signif <- data.frame(
  organism = c("MRSA", "MRSA", "MRSA",
               "Non-mucoid P. aeruginosa",
               "P. aeruginosa", "P. aeruginosa",
               "S. aureus", "S. aureus", "S. aureus", "S. aureus", 
               "S. aureus", "S. aureus",
               "S. aureus morphological variant","S. aureus morphological variant",
               "S. aureus morphological variant", "S. aureus morphological variant",
               "S. aureus morphological variant",
               "S. pneumoniae", "S. pneumoniae", "S. pneumoniae", "S. pneumoniae"),
  left_class = c("Aminoglycoside", "Folate inhibitor", "Glycopeptide",
                 "BL/BLI",
                 "BL/BLI", "BL/BLI",
                 "Cephalosporin", "Glycopeptide", "Glycopeptide", "Glycopeptide",
                 "Carbapenem", "Folate inhibitor",
                 "Aminoglycoside", "Folate inhibitor",
                 "Glycopeptide", "Glycopeptide",
                 "Glycopeptide",
                 "Carbapenem", "Fluoroquinolone", "Glycopeptide", "Cephalosporin"),
  right_class = c("Glycopeptide", "Glycopeptide", "Tetracycline",
                  "Polymyxin",
                  "Folate inhibitor", "Tetracycline",
                  "Glycopeptide", "Lipopeptide", "Nitrofuran", "Oxazolidinone",
                  "Glycopeptide", "Glycopeptide",
                  "Glycopeptide", "Glycopeptide",
                  "Lipopeptide", "Oxazolidinone", 
                  "Tetracycline",
                  "Glycopeptide", "Glycopeptide", "Oxazolidinone", "Glycopeptide"),
  annotation = c("*", "*", "**", 
                 "ns",
                 "*", "ns",
                 "ns", "ns", "ns", "ns",
                 "*", "*",
                 "ns", "ns",
                 "ns", "ns",
                 "ns",
                 "ns", "ns", "ns", "*"),
  y_position = c(1.12, 1.05, 1.15, 
                 1.0,
                 1.0, 1.08,
                 1.05, 1.12, 1.1, 1.08,
                 1.0, 1.03,
                 1.02, 1.07,
                 1.14, 1.01,
                 1.0,
                 1.05, 1.1, 1.05, 1.0))

# individual species level
# STREPTOCOCCUS PNEUMONIAE ------------------------------------------------

# analysis
all_strep <- data_susceptibility |>
  filter(organism %in% c("STREPTOCOCCUS PNEUMONIAE"))
n_distinct(all_strep$anon_id)
n_distinct(all_strep$antibiotic)

final_strep <- all_strep |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() 

n_distinct(final_strep$antibiotic)

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
# analysis
all_staph <- data_susceptibility |>
  filter(organism %in% c("STAPHYLOCOCCUS AUREUS",
                         "STAPH AUREUS(COLONY VARIANT - SMALL COLONY OR OTHER MORPHOTYPE)",
                         "STAPH AUREUS {MRSA}"))
n_distinct(all_staph$anon_id)
n_distinct(all_staph$antibiotic)

final_staph <- all_staph |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() 

n_distinct(final_staph$antibiotic)

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
# analysis
all_kleb <- data_susceptibility |>
  filter(organism %in% c("KLEBSIELLA PNEUMONIAE", 
                         "KLEBSIELLA PNEUMONIAE SSP. OZAENAE", 
                         "KLEBSIELLA PNEUMONIAE (CARBAPENEM RESISTANT)"))
n_distinct(all_kleb$anon_id)
n_distinct(all_kleb$antibiotic)

final_kleb <- all_kleb |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() 

n_distinct(final_kleb$antibiotic)

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
# analysis
all_pseud <- data_susceptibility |>
  filter(organism %in% c("PSEUDOMONAS AERUGINOSA", "MUCOID PSEUDOMONAS AERUGINOSA",
                         "PSEUDOMONAS AERUGINOSA (NON-MUCOID CF)"))
n_distinct(all_pseud$anon_id)
n_distinct(all_pseud$antibiotic)

final_pseud <- all_pseud |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() 

n_distinct(final_pseud$antibiotic)

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

# analysis
all_ecoli <- data_susceptibility |>
  filter(organism %in% c("ESCHERICHIA COLI",
                         "ESCHERICHIA COLI (CARBAPENEM RESISTANT)"))
n_distinct(all_ecoli$anon_id)
n_distinct(all_ecoli$antibiotic)

final_ecoli <- all_ecoli |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() 

n_distinct(final_ecoli$antibiotic)

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
