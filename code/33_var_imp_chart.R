## Purpose of script: plot variable importance as bar chart
## Authors: GEOG 274
## Date: Spring, 2025

rm(list=ls())

## 1. Load libraries and read data
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# 1a. Point to the birds evaluation folder (since our target species is a bird)
eval_dir <- here("results", "birds_results", "evaluations")

# 1b. Read the submodel‐level variable importance CSV
submodel_var_imp <- read_csv(file.path(eval_dir, "all_submodel_variable_importance.csv"))

# 1c. Pivot to long format: one row per (Species, Scenario, Submodel, Variable, Importance)
submodel_var_imp_long <- submodel_var_imp %>%
  pivot_longer(
    cols      = c(slope, aspect, flow_acc, dist_coast, solar,
                  bio1, bio15, bio17, bio18, bio3, bio5, bio6),
    names_to  = "Variable",
    values_to = "Importance"
  )

## 2. Specify the species and scenario of interest
target_species <- "Haliaeetus leucocephalus"
target_scenario <- "ssp370"
species_name <- "Bald Eagle"

## 3. Filter for that species & scenario, then pick top 3 vars per Submodel
vip_top3_submodel <- submodel_var_imp_long %>%
  filter(
    Species  == target_species,
    Scenario == target_scenario
  ) %>%
  group_by(Submodel) %>%
  slice_max(order_by = Importance, n = 3) %>%
  ungroup() %>%
  # Reorder Variable within each Submodel so bars go high → low
  group_by(Submodel) %>%
  mutate(Variable = fct_reorder(Variable, Importance, .desc = TRUE)) %>%
  ungroup()

# (Optional) Inspect:
vip_top3_submodel %>% arrange(Submodel, desc(Importance))

## 4. Plot: one facet per Submodel, showing top 3 variables
ggplot(vip_top3_submodel, aes(x = Variable, y = Importance, fill = Variable)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ Submodel, scales = "free_x", nrow = 1) +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  labs(
    title = paste0(species_name, " (", target_scenario, "): Top 3 Variables by Submodel"),
    x     = "Variable (high → low importance)",
    y     = "Percent Importance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text       = element_text(face = "bold"),
    panel.spacing.x  = unit(0.5, "lines")
  )

