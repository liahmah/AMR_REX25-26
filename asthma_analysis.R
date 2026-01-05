
# Setup ---------------------------------------------------------------

library(tidyverse)
library(MASS)
library(emmeans)

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

# summarizing 

## totals
summarized <- merged_asthma_suscep |>
  count(susceptibility)

## summarizing number of resistant/susceptible per organism by antibiotic
summarized_organism <- merged_asthma_suscep |>
  count(antibiotic, organism, susceptibility) |>
  rename('count' = 'n')

## or without organism:
summarized_antibiotic <- merged_asthma_suscep |>
  count(antibiotic, susceptibility) |>
  rename('count' = 'n')
  
# but can't compare to non-asthma patients. compare to who - compare within

# Relationship analysis ---------------------------------------------------

# trying poisson regression
poisson_antibiotic <- glm(count ~ susceptibility,
             family = poisson(link = "log"),
             data = summarized_antibiotic)

summary(poisson_antibiotic)

## checking data overdispersion
deviance(poisson_antibiotic) # > 1.5, switching to negative binomial
mean(summarized_antibiotic$count)
var(summarized_antibiotic$count) # var >> mean, hence neg binom more appropriate

# trying negative binomial:
neg_bio_poisson <- glm.nb(count ~ susceptibility, data = summarized_antibiotic)
summary(neg_bio_poisson)

# could also use chi square for antibiotic-specific comparisons
penicillin <- summarized_antibiotic |>
  filter(antibiotic == "Penicillin", 
         susceptibility == "Resistant" | susceptibility == "Susceptible")

chisq_penicillin <- table(penicillin$antibiotic, penicillin$susceptibility)
chisq.test(chisq_penicillin)

chisq <- chisq.test(chisq_penicillin)
chisq$expected # expected value of 1, cannot use chisq

# trying Fisher, maybe for resistant only
resistant <- summarized_antibiotic |>
  filter(susceptibility == "Resistant" | susceptibility == "Susceptible") |>
  filter(antibiotic == "Penicillin" | antibiotic == "Vancomycin")

fisher.test(table(resistant$antibiotic, resistant$susceptibility))

# to consider significant difference between counts: rec. to use binomial test

sum_anti_wide <- summarized_antibiotic |>
  pivot_wider(names_from = susceptibility,
              values_from = count)

sum_anti_wide$antibiotic <- as.factor(sum_anti_wide$antibiotic)

typeof(sum_anti_wide$antibiotic) # confirms factor/integer

levels(sum_anti_wide$antibiotic) 
## reference is amikacin (alphabetical); ref used for odds ratios
### setting reference as vancomycin, since is last-resort antibiotic
sum_anti_wide$antibiotic <- relevel(sum_anti_wide$antibiotic, ref = "Vancomycin")

### running binomial regression
sum_anti_model <- glm(cbind(Resistant, Susceptible) ~ antibiotic,
    family = binomial,
    data = sum_anti_wide)

summary(sum_anti_model) # 47/51 significant differences, 43/51 3 degrees signif
### low coefficient (more signif) indicates lower proportion resistance:suscept

### looking at odds ratio:
exp(coef(sum_anti_model)) # lower coeff/intercept = higher odds susceptibility 
#### lowest coefficients: amox/clavulanic acid 0.02, doxycycline 0.03

exp(confint(sum_anti_model)) # <1 = higher odds suscep, 2.5 - 97.5 doesn't cross 1.0
#### confirms daptomycin as best

## or using Tukey for no ref - use this to compare between any/all antibiotics
pairs(emmeans(sum_anti_model, ~ antibiotic), adjust = "tukey")
