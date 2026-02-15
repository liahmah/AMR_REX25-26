
# libraries and datasets --------------------------------------------------

library(tidyverse)
library(MASS)
library(emmeans)
library(brant)
library(broom)
library(ggsignif)

asthma_dataset <- read_csv("datasets/merged_asthma_susceptibility.csv")

# wrangling ---------------------------------------------------------------

# filtering for pneumonia
pneumonia <- asthma_dataset |>
  filter(organism == "STREPTOCOCCUS PNEUMONIAE")

# looking into count of each antibiotic for sample size
meta_pneum <- pneumonia |>
  count(antibiotic) # 12 antibiotics

# seeing if one should be dropped due to outlying sample size, to avoid skew
ggplot(meta_pneum, aes(x = 1, y = n)) + # using dummy x
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.1)) +
  geom_label(aes(label = antibiotic), fill = "white", alpha = 0.7,
             position = position_jitter(width = 0.05), vjust = 0.5) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "count/antibiotic", title = "Boxplot of count/antibiotic")
  # penicillin at upper bound: use as reference
  # linezolid at lower bound: very small sample size, drop

# adding counts of susceptibility per antibiotic, filtering for no linezolid
pneumonia_filtered <- pneumonia |>
  filter(antibiotic != "Linezolid") |>
  count(antibiotic, organism, susceptibility) |> 
  rename('count' = 'n')
  
# coding "null" as resistant to increase value of significantly susceptible results
pneumonia_count  <- pneumonia_filtered |>
  mutate(susceptibility = if_else(susceptibility == "Null",
                                  "Resistant", susceptibility)) |>
  group_by(antibiotic, susceptibility) |>
  summarise(count = sum(count), .groups = "drop")

# changing to factor type for levelling
pneumonia_levelled <- pneumonia_count
pneumonia_levelled$antibiotic <- as.factor(pneumonia_levelled$antibiotic)

# setting pneumonia as reference value, considering sample size and commonality
pneumonia_levelled$antibiotic <- relevel(pneumonia_levelled$antibiotic, 
                                         ref = "Penicillin")

  # Stats analysis ----------------------------------------------------------

# ordinal logistic regression accounts for intermediate susceptibility
## formatting for order susceptible < intermediate < resistant

pneumonia_ordered <- pneumonia_levelled
pneumonia_ordered$susceptibility <- factor(pneumonia_ordered$susceptibility,
  levels = c("Susceptible", "Intermediate", "Resistant"),
  ordered = TRUE)

# proportional odds model
olr_model <- polr(susceptibility ~ antibiotic, data = pneumonia_ordered, 
                  Hess = TRUE)
summary(olr_model) # neg coefficients lean towards susceptibilty than resistance
### antibiotics levofloxacin, linezolid, vancomycin most negative

### visualizing it raw
pneumonia_ordered |>
ggplot(aes(x = antibiotic, fill = susceptibility)) +
  geom_bar(position = "fill") +
  ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## converting to odds ratio
ratio <- exp(coef(olr_model))
ratio # <1 signifies better susceptibility than penicillin
### meropenem, vancomycin, levofloxacin only <1 

## tukey pairwise 
tukey <- emmeans(olr_model, pairwise ~ antibiotic, adjust = "tukey",
        type = "response")
tukey
### confirms. estimate ~ 18 for all 3 relative to ~ 0 for all others
### for those 3 only, p.value < 0.0001

## visualizing tukey reslts w/ signif

### extracting coefficients (coefficient == "estimate")
tukey_df <- as.data.frame(summary(tukey)) |>
  mutate(lower = estimate - 1.96*SE,  # 95% CI lower
    upper = estimate + 1.96*SE,       # 95% CI upper
    comparison = contrast)            # labels for plotting

### plotting
tukey_df |>
  ggplot(aes(x = reorder(comparison, estimate), y = estimate)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Tukey Pairs",
       y = "Log-odds coefficient (Tukey-adjusted)",
       title = "Tukey-adjusted Coefficients from Proportional Odds Model") 
  geom_signif(comparisons = list(c("Moxiflaxin-Vancomycin",""), c("",""), ...)
              y_position = c(20),
                annotations = c("***"))

## confirming significance
brant(olr_model) # some sparse outcomes, limiting the brant test, but does confirm

## can also confirm by binary collapse, so considering intermediate as resistant
pneumonia_binary <- pneumonia_count
pneumonia_binary$non_susceptible <-  pneumonia_binary$susceptibility %in% 
  c("Intermediate", "Resistant")

bin_model <- glm(non_susceptible ~ antibiotic,
  family = binomial,
  data = pneumonia_binary)
summary(bin_model)

exp(coef(bin_model)) # <0 are more susceptible

emmeans(bin_model, pairwise ~ antibiotic, adjust = "tukey",
        type = "response")

## in both emmeans and glm, varying sample sizes are accounted for

# Discussion --------------------------------------------------------------

# checking vancomycin and levofloxacin sample size:
subset(pneumonia, antibiotic == "Vancomycin") # 1973 observations
subset(pneumonia, antibiotic == "Levofloxacin") # 1026 observations

# demographic stuff
length(unique(pneumonia$anon_id)) # 178 unique patients

# is there a relationship between antibiotic and susceptibility in general?
## chi squared test, for multi-category single predictor and multi-category 
## single response, using counts
chisq <- xtabs(count ~ antibiotic + susceptibility, data = pneumonia_count)
chisq.test(chisq) # p-value very small, shows that there is a significant 
                  # relationship between predictor and response
## to see which antibiotic is MOST associated with susceptibility (highest
## standardized residual)
chisq_res <- chisq.test(chisq)
chisq_res$stdres # vancomycin > moxifloxacin > levofloxacin > ceftriaxone for 
                 # positive relationship to susceptible

# Conclusion --------------------------------------------------------------

## Levofloxacin and Vancomycin have significantly better odds of susceptibility
## against S. pneumoniae than penicillin, as compared to 10 other antibiotics
## (in asthmatic individuals) when using odds ratios, but not looking at binomial

