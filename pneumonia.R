# Loading -----------------------------------------------------------------

library(tidyverse)
library(emmeans)
library(scales)
library(showtext)
library(svglite)
library(ggsignif)
library(brglm2)

data_susceptibility <- read_csv("datasets/microbiology_cultures_cohort.csv") |>
  mutate(across(c("pat_enc_csn_id_coded", "order_proc_id_coded"), as.double))

# dataset components
# n_distinct(data_susceptibility$anon_id) # 283715
# n_distinct(data_susceptibility$antibiotic) # 55 
# n_distinct(data_susceptibility$organism) # 315

## identifying pneumonia-causing species of interest:
# species <- as.data.frame(unique(demo_suscep$organism))
# Filtering for Streptococcus pneumoniae, yes
# Haemophilus influenzae, no
# Klebsiella pneumoniae, yes, 3 subvariants
# Pseudomonas aeruginosa, yes, 3 subvariants
# Escherichia coli, yes, 2 subvariants
# Staphylococcus aureus, yes, 3 subvariants

# excluding unidentified species
# susceptibility coded as true (susceptible) or false (resistant)

data_organisms <- data_susceptibility |>
  filter(organism %in% c("STREPTOCOCCUS PNEUMONIAE", 
                         "STAPHYLOCOCCUS AUREUS",
                         "STAPH AUREUS(COLONY VARIANT - SMALL COLONY OR OTHER MORPHOTYPE)",
                         "STAPH AUREUS {MRSA}", 
                         "KLEBSIELLA PNEUMONIAE", 
                         "KLEBSIELLA PNEUMONIAE SSP. OZAENAE", 
                         "KLEBSIELLA PNEUMONIAE (CARBAPENEM RESISTANT)", 
                         "PSEUDOMONAS AERUGINOSA", 
                         "MUCOID PSEUDOMONAS AERUGINOSA",
                         "PSEUDOMONAS AERUGINOSA (NON-MUCOID CF)", 
                         "ESCHERICHIA COLI",
                         "ESCHERICHIA COLI (CARBAPENEM RESISTANT)"))

# dataset components
n_distinct(data_organisms$anon_id) # 51778
n_distinct(data_organisms$antibiotic) # 50 
n_distinct(data_organisms$organism) # 12

# susceptibility distribution plot ---------------------------------------------

data_wrangled <- data_organisms |>
  group_by(antibiotic) |>
  filter(n() > 100) |>  # removing low-frequency antibiotics
  ungroup() |>
  group_by(organism) |>
  filter(n() > 400) |>
  ungroup() 

# abbreviating antibiotic names for improved readability
abx_lookup <- c("Ciprofloxacin" = "CIP", "Nitrofurantoin" = "NIT",  
                "Ampicillin" = "AMP","Trimethoprim/Sulfamethoxazole" = "SXT",  
                "Gentamicin" = "GEN", "Cefazolin" = "CZO", "Levofloxacin" = "LVX", 
                "Piperacillin/Tazobactam" = "TZP", "Cefuroxime" = "CXM", 
                "Ampicillin/Sulbactam" = "SAM", "Erythromycin" = "ERY",
                "Vancomycin" = "VAN", "Clindamycin" = "CLI", "Oxacillin" = "OXA",
                "Penicillin" = "PEN", "Tetracycline" = "TCY", "Amikacin" = "AMK", 
                "Tobramycin" = "TOB", "Linezolid" = "LNZ", "Moxifloxacin" = "MFX",
                "Ertapenem" = "ETP", "Cefepime" = "FEP", "Ceftriaxone" = "CRO",
                "Ceftazidime" = "CAZ", "Aztreonam" = "ATM", "Meropenem" = "MEM",
                "Imipenem" = "IPM", "Colistin" = "COL", "Piperacillin" = "PIP",
                "Ticarcillin" = "TIC", "Doxycycline" = "DOX", "Cefotaxime" = "CTX",
                "Amoxicillin/Clavulanic Acid" = "AMC", 
                "Cephalexin/Cephalothin" = "LEX/CEP", "Daptomycin" = "DAP",
                "Tigecycline" = "TGC", "Cefotetan" = "CTT", "Cefoxitin" = "FOX",
                "Doripenem" = "DOR", "Fosfomycin" = "FOS", "Ceftaroline" = "CPT",
                "Ceftolozane/Tazobactam" = "CZT", "Ceftazidime/Avibactam" = "CZA",
                "Imipenem/Ebactam" = "IMR")

# grouping by class
class_lookup <- c("CIP" = "Fluoroquinolone", "LVX" = "Fluoroquinolone",
                  "MFX" = "Fluoroquinolone", "GEN" = "Aminoglycoside",
                  "AMK" = "Aminoglycoside", "TOB" = "Aminoglycoside",
                  "IPM" = "Carbapenem", "MEM" = "Carbapenem", "DOR" = "Carbapenem",
                  "ETP" = "Carbapenem", "AMC" = "BL/BLI", "SAM" = "BL/BLI",
                  "TZP" = "BL/BLI", "CZT" = "BL/BLI", "CZA" = "BL/BLI",
                  "AMP" = "Penicillin", "PEN" = "Penicillin", "OXA" = "Penicillin",
                  "PIP" = "Penicillin", "TIC" = "Penicillin", 
                  "CZO" = "Cephalosporin", "CXM" = "Cephalosporin",
                  "FOX" = "Cephalosporin", "CTT" = "Cephalosporin", 
                  "CRO" = "Cephalosporin", "CTX" = "Cephalosporin", 
                  "CAZ" = "Cephalosporin", "FEP" = "Cephalosporin",
                  "CPT" = "Cephalosporin", "ATM" = "Monobactam",
                  "VAN" = "Glycopeptide", "DAP" = "Lipopeptide",
                  "LNZ" = "Oxazolidinone", "TCY" = "Tetracycline",
                  "DOX" = "Tetracycline",  "TGC" = "Tetracycline",
                  "ERY" = "Macrolide", "CLI" = "Lincosamide", "IMR" = "BL/BLI",
                  "SXT" = "Folate inhibitor", "NIT" = "Nitrofuran",
                  "COL" = "Polymyxin", "FOS" = "Phosphonic acid", 
                  "LEX/CEP" = "Cephalosporin")



data_wrangled <- data_wrangled |>
  mutate(
    antibiotic = recode(antibiotic, !!!abx_lookup),
    class = recode(antibiotic, !!!class_lookup),
    genus = recode(organism,
                   `ESCHERICHIA COLI` = "E. coli",
                   `KLEBSIELLA PNEUMONIAE` = "K. pneumoniae",
                   `PSEUDOMONAS AERUGINOSA` = "P. aeruginosa",
                   `MUCOID PSEUDOMONAS AERUGINOSA` = "P. aeruginosa",
                   `PSEUDOMONAS AERUGINOSA (NON-MUCOID CF)` = "P. aeruginosa",
                   `STAPH AUREUS {MRSA}` = "S. aureus",
                   `STAPH AUREUS(COLONY VARIANT - SMALL COLONY OR OTHER MORPHOTYPE)`
                        = "S. aureus", 
                   `STAPHYLOCOCCUS AUREUS` = "S. aureus",
                   `STREPTOCOCCUS PNEUMONIAE` = "S. pneumoniae"),
         organism = recode(organism,
                           `ESCHERICHIA COLI` = "E. coli",
                           `KLEBSIELLA PNEUMONIAE` = "K. pneumoniae",
                           `PSEUDOMONAS AERUGINOSA` = "P. aeruginosa",
                           `MUCOID PSEUDOMONAS AERUGINOSA` = "Mucoid P. aeruginosa",
                           `PSEUDOMONAS AERUGINOSA (NON-MUCOID CF)` = 
                          "Non-mucoid P. aeruginosa",
                          `STAPH AUREUS {MRSA}` = "MRSA",
                          `STAPH AUREUS(COLONY VARIANT - SMALL COLONY OR OTHER MORPHOTYPE)`
                          = "S. aureus morphological variant", 
                          `STAPHYLOCOCCUS AUREUS` = "S. aureus",
                          `STREPTOCOCCUS PNEUMONIAE` = "S. pneumoniae"))

genus_n <- data_wrangled |>
  count(genus) |>
  mutate(label = paste0("italic('", genus, "')~'(n = ", n, ")'"))

data_plot <- data_wrangled |>
  left_join(genus_n, by = "genus")
         

showtext_auto(TRUE)
font_add(family = "TimesNewRoman", regular = "C:/Windows/Fonts/times.ttf",
         italic  = "C:/Windows/Fonts/timesi.ttf")
font_add(family = "TimesExtraBold", regular = "C:/Windows/Fonts/timesbd.ttf")

data_organisms_plot <- data_plot |>
  ggplot(aes(x = antibiotic, fill = susceptibility)) +
  geom_bar() +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  facet_wrap(~ label, labeller = label_parsed, 
             ncol = 2, 
             scales = "free_x") +
  scale_fill_manual(breaks = c("Susceptible", "Intermediate", "Resistant", 
                               "Inconclusive", "Null"),
                    values = c("Susceptible" = "#1f78b4",
                               "Intermediate" = "#FDBF00",
                               "Resistant" = "#E31A1C",
                               "Inconclusive" = "#8C564B",
                               "Null" = "grey70")) +   
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size = 15),
        axis.title.x = element_text(face = "bold", 
                                    size = 25, 
                                    family = "TimesExtraBold"),
        axis.title.y = element_text(face = "bold", 
                                    size = 25, 
                                    family = "TimesExtraBold",
                                    margin = margin(r = 8)),
        strip.text = element_text(size = 32, 
                                  face = "italic", 
                                  family = "TimesNewRoman"),
        legend.title = element_text(family = "TimesExtraBold", 
                                    face = "bold", 
                                    size = 25),
        legend.text = element_text(family = "TimesNewRoman", 
                                   size = 20),
        legend.position = "right",
        legend.key.height = unit(1.5, "lines"),
        legend.spacing.y = unit(0.8, "lines"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey85"),
        axis.ticks.length = unit(3, "pt"),
        panel.spacing = unit(2, "lines")) +
  labs(x = "Antibiotic (Abbreviated)",
       y = "Number of Isolates (log[10] scale)",
       fill = "Susceptibility Outcome")
data_organisms_plot

ggsave("SusceptibilityProfiles.svg",
  plot = data_organisms_plot,
  width = 24,
  height = 18,
  units = "in",
  device = svglite::svglite)


# class probability plot -------------------------------------------------------

all_logical <- data_wrangled |>
  filter(susceptibility %in% c("Susceptible", "Resistant")) |>
  mutate(suscep_logical = susceptibility == "Susceptible") |>
  group_by(organism, class) |>
  summarise(susceptible = sum(suscep_logical),
            resistant = sum(!suscep_logical),
            .groups = "drop")

# modelling
all_bin_model <- glm(cbind(susceptible, resistant) ~ organism * class,
                     data = all_logical,
                     family = binomial,
                     method = "brglmFit")
# some antibiotics are fully susceptible or resistant, so probabilities will be 0 or 1

all_emm <- emmeans(all_bin_model, ~ class | organism, type = "response")

emm_df <- as.data.frame(all_emm)

ranked_classes <- emm_df |>
  group_by(organism) |>
  arrange(dplyr::desc(prob), .by_group = TRUE) |>
  filter(!is.na(prob))

pairwise_results <- pairs(all_emm, adjust = "tukey") |>
  as.data.frame() |>
  filter(!is.na(SE)) |>
  mutate(p.signif = case_when(p.value < 0.001 ~ "***",
                              p.value < 0.01 ~ "**",
                              p.value < 0.05 ~ "*",
                              TRUE ~ "ns"))
         

signif_data <- data.frame(
  organism = c("E. coli", "E. coli", "E. coli", "E. coli", 
               "E. coli", "E. coli", "E. coli", "E. coli", 
               "E. coli", "E. coli",
               "K. pneumoniae", "K. pneumoniae", "K. pneumoniae", "K. pneumoniae", 
               "K. pneumoniae", "K. pneumoniae", "K. pneumoniae", "K. pneumoniae", 
               "K. pneumoniae", "K. pneumoniae", 
               "Mucoid P. aeruginosa",
               "Non-mucoid P. aeruginosa", "Non-mucoid P. aeruginosa", 
               "Non-mucoid P. aeruginosa", "Non-mucoid P. aeruginosa", 
               "Non-mucoid P. aeruginosa", "Non-mucoid P. aeruginosa", 
               "P. aeruginosa", "P. aeruginosa", "P. aeruginosa", "P. aeruginosa", 
               "P. aeruginosa", "P. aeruginosa", "P. aeruginosa", "P. aeruginosa",
               "S. aureus", "S. aureus", "S. aureus", "S. aureus", 
               "S. aureus", "S. aureus", "S. aureus", "S. aureus", 
               "S. aureus morphological variant", "S. aureus morphological variant",
               "S. aureus morphological variant", "S. aureus morphological variant",
               "MRSA", "MRSA", "MRSA", "MRSA", 
               "MRSA", "MRSA", "MRSA", "MRSA", 
               "S. pneumoniae", "S. pneumoniae", "S. pneumoniae", "S. pneumoniae",
               "S. pneumoniae", "S. pneumoniae"),
  class = c("Nitrofuran", "Phosphonic acid", "Aminoglycoside", "BL/BLI",
            "Monobactam", "Cephalosporin","Fluoroquinolone","Folate inhibitor",
            "Tetracycline", "Penicillin",
            "Aminoglycoside", "BL/BLI", "Fluoroquinolone", "Monobactam",
            "Cephalosporin", "Folate inhibitor", "Tetracycline", "Phosphonic acid",
            "Nitrofuran", "Penicillin",
            "Polymyxin",
            "Aminoglycoside", "Carbapenem",
            "Cephalosporin", "Fluoroquinolone",
            "Monobactam", "Penicillin",
            "Cephalosporin", "Aminoglycoside", "Carbapenem", "Fluoroquinolone",
            "Monobactam", "Polymyxin", "Penicillin", "Folate inhibitor",
            "Folate inhibitor", "Aminoglycoside", "Tetracycline", "Carbapenem",
            "Fluoroquinolone", "Lincosamide", "Penicillin", "Macrolide",
            "Penicillin", "Macrolide", "Lincosamide","Fluoroquinolone",
            "Folate inhibitor", "Aminoglycoside", "Tetracycline", "Cephalosporin",
            "Lincosamide", "Fluoroquinolone", "Macrolide", "Penicillin",
            "Cephalosporin", "Folate inhibitor", "Tetracycline", "Penicillin", 
            "Lincosamide", "Macrolide"),
  y = c(1.02, 0.98, 0.97, 0.97, 
        0.94, 0.93, 0.81, 0.74, 
        0.68, 0.57,
        1.02, 1.0, 0.97, 0.97, 
        0.97, 0.9, 0.83, 0.75, 
        0.49, 0.05,
        1.1,
        0.66, 0.72, 
        0.77, 0.63,
        0.77, 0.76,
        0.95, 0.94, 0.88, 0.85, 
        0.81, 0.79, 0.74, 0.41,
        1.02, 1.0, 0.96, 0.91,
        0.73, 0.69, 0.61, 0.59,
        0.44, 0.21, 0.32, 0.41,
        1.0, 0.97, 0.96, 0.68, 
        0.44, 0.24, 0.2, 0.07,
        0.95, 0.84, 0.82, 0.81,
        0.8, 0.74),
  label = c("***", "***", "***", "***",
            "***", "***", "***", "***",
            "***", "***",
            "***", "***", "***", "***",
            "***", "***", "***", "***",
            "***", "***",
            "ns",
            "***", "***",
            "***", "***",
            "***", "***",
            "***", "***", "***", "***", 
            "***", "***", "***", "*",
            "*", "***", "***", "*",
            "***", "***", "***", "***", 
            "***", "***", "***", "***", 
            "*", "*", "**", "***",
            "***", "***", "***", "***", 
            "*", "***", "***", "***",
            "***", "***"))

organism_n <- data_wrangled |>
  count(organism) |>
  mutate(facet = paste0("italic('", organism, "')~'(n = ", n, ")'"))

top_classes <- ranked_classes |>
  dplyr::group_by(organism) |>
  dplyr::slice_max(prob, n = 1, with_ties = FALSE) |>
  dplyr::mutate(is_top = TRUE) |>
  dplyr::select(organism, class, is_top)

plot_df <- ranked_classes |>
  left_join(top_classes, by = c("organism","class")) |>
  mutate(is_top = dplyr::if_else(is.na(is_top), FALSE, TRUE)) |>
  left_join(organism_n, by = "organism") |>
  left_join(signif_data, by = c("organism", "class"))
#
susceptibility_probability_plot <- plot_df |>
  ggplot(aes(x = class, y = prob)) +
  geom_col(aes(color = is_top),
           width = 0.8,
           fill = "white",
           linewidth = 1.2) +
  facet_wrap(~ facet, labeller = label_parsed,
             scales = "free_x") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_minimal(base_size = 15) +
  labs(y = "Predicted Probability", 
       x = "Antibiotic Class") +
  geom_text(aes(x = class, y = y, label = label),
            inherit.aes = FALSE,
            size = 8) +
  theme(axis.text.x = element_text(angle = 45, 
                                   size = 15,
                                   hjust = 1),
        axis.title.x = element_text(face = "bold", 
                                    size = 25, 
                                    family = "TimesExtraBold"),
        axis.title.y = element_text(face = "bold", 
                                    size = 25, 
                                    family = "TimesExtraBold"),
        strip.text = element_text(size = 32, 
                                  face = "italic", 
                                  family = "TimesNewRoman"),
        legend.position = "none",
        axis.ticks.length = unit(3, "pt"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(1.5, "lines")) +
  scale_color_manual(values = c("TRUE" = "#D62728", "FALSE" = "#0072B2"))

susceptibility_probability_plot


ggsave("SusceptibilityProbability.svg",
       plot = susceptibility_probability_plot,
       width = 24,
       height = 18,
       units = "in",
       device = svglite::svglite)
