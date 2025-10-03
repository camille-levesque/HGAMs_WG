# Script to plot trends and calculate derivative-based indicators on a 
# multi-species dataset with model "GS" from Pedersen et al. (2019)'s paper.

#-----------------------------------------------------------------------------
# STEP 1: PACKAGES & LIBRARIES
#-----------------------------------------------------------------------------

# install.packages(c("mgcv", "dplyr", "ggplot2", "tidyr", "mvtnorm"))

library(mgcv)    # For fitting Generalized Additive Models (GAMs) and Hierarchical GAMs
library(dplyr)   # For data manipulation 
library(ggplot2) # For data visualization
library(tidyr)   # For reshaping and tidying data
library(mvtnorm) # For working with multivariate normal and t-distributions
library(here)    # For handling file paths relative to the project root
library(gratia)

#-----------------------------------------------------------------------------
# STEP 2: IMPORT AND PRE-PROCESS THE DATA
#-----------------------------------------------------------------------------
# Loading the dataset
data_195 <- read.csv("data/clean/data_195.csv") # This is the cleaned data (in HGAMs_WG/data/clean)

# Aggregating biomass data to get a yearly abundance for each species.
community_ts <- data_195 %>%
  filter(!is.na(YEAR) & !is.na(valid_name) & valid_name != "") %>%
  group_by(YEAR, valid_name) %>%
  summarise(ABUNDANCE = n(), .groups = 'drop') %>%
  rename(year = YEAR, species = valid_name, abundance = ABUNDANCE)

# Selecting the top 30 most frequently observed species (i.e., highest in abundance)
top_species <- community_ts %>%
  group_by(species) %>%
  summarise(total_abundance = sum(abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 30) %>%
  pull(species)

community_ts_subset <- community_ts %>%
  filter(species %in% top_species) %>%
  mutate(species = as.factor(species))

print(head(community_ts_subset))

#-----------------------------------------------------------------------------
# STEP 3: FIT POISSON MODEL GS
#-----------------------------------------------------------------------------
# Here, we are fitting model "GS" from Pedersen et al. (2019) paper.
# As we are working with abundance (count) data, we'll use the Poisson family.

# MODEL: The 'GS' Model (Global smooth + species-specific deviations)
gam_model_GS <- gam(
  abundance ~ s(year, bs = "tp") + # Global relationship
    s(year, by = species, bs = "fs") + # Species-specific smoother: a factor-smoother interaction of year and species
    s(species, bs="re"), # Species as random effects (gives an intercept per species)
  data = community_ts_subset,
  family = poisson(), # Using Poisson family
  method = "REML"
)

#-----------------------------------------------------------------------------
# STEP 4: DERIVATIVES AND INDICATORS FOR THE 'GS' MODEL
#-----------------------------------------------------------------------------
# This section repeats the process, but this time for the GS model.

# Generate posterior simulations from the GS model
set.seed(42)
sim_lp_GS <- predict(gam_model_GS, newdata = predict_data, type = "lpmatrix")
sim_coef_GS <- rmvnorm(n_sim, coef(gam_model_GS), vcov(gam_model_GS, unconditional = TRUE))

# Calculate predicted values and derivatives for the GS model
pred_original_GS <- exp(sim_lp_GS %*% t(sim_coef_GS))
pred_p_eps_GS <- exp(predict(gam_model_GS, newdata = predict_data_p_eps, type = "lpmatrix") %*% t(sim_coef_GS))
pred_m_eps_GS <- exp(predict(gam_model_GS, newdata = predict_data_m_eps, type = "lpmatrix") %*% t(sim_coef_GS))
first_derivative_GS <- (pred_p_eps_GS - pred_m_eps_GS) / (2 * eps)
per_capita_rate_GS <- first_derivative_GS / (pred_original_GS + 1e-9)

# Reshape simulation results for the GS model
sim_deriv_long_GS <- as.data.frame(first_derivative_GS) %>%
  mutate(row = 1:n()) %>%
  pivot_longer(-row, names_to = "sim_id", values_to = "derivative")
sim_per_capita_long_GS <- as.data.frame(per_capita_rate_GS) %>%
  mutate(row = 1:n()) %>%
  pivot_longer(-row, names_to = "sim_id", values_to = "per_capita_rate")

sim_results_GS <- predict_data %>%
  mutate(row = 1:n()) %>%
  left_join(sim_deriv_long_GS, by = "row") %>%
  left_join(sim_per_capita_long_GS, by = c("row", "sim_id"))

# Calculate community indicators for the GS model
community_indicators_GS <- sim_results_GS %>%
  group_by(year, sim_id) %>%
  summarise(
    mean_rate_of_change = mean(derivative, na.rm = TRUE),
    mean_per_capita_rate = mean(per_capita_rate, na.rm = TRUE),
    sd_per_capita_rate = sd(per_capita_rate, na.rm = TRUE),
    .groups = 'drop'
  )

# Final summary of indicators for the GS model
final_indicators_GS <- community_indicators_GS %>%
  group_by(year) %>%
  summarise(
    across(
      .cols = c(mean_rate_of_change, mean_per_capita_rate, sd_per_capita_rate),
      .fns = list(
        median = ~median(.x, na.rm = TRUE),
        lower_ci = ~quantile(.x, 0.025, na.rm = TRUE),
        upper_ci = ~quantile(.x, 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = 'drop'
  ) %>%
  mutate(model_type = "GS Model") # Add a label for plotting

print("--- Indicators from GS Model ---")
print(head(final_indicators_GS))

#-----------------------------------------------------------------------------
# STEP 5: PLOT DERIVATIVE-BASED INDICATORS
#-----------------------------------------------------------------------------

# Plot 1: Mean Rate of Change Comparison
plot1_mean_rof <- ggplot(final_indicators_GS, aes(x = year)) +
  geom_ribbon(aes(ymin = mean_rate_of_change_lower_ci, ymax = mean_rate_of_change_upper_ci), alpha = 0.2) +
  geom_line(aes(y = mean_rate_of_change_median), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    y = "Mean rate of change", x = "Year",
  ) +
  theme_minimal()

# Plot 2: Mean Per-Capita Rate of Change Comparison
plot2_mean_percap_rof <- ggplot(final_indicators_GS, aes(x = year)) +
  geom_ribbon(aes(ymin = mean_per_capita_rate_lower_ci, ymax = mean_per_capita_rate_upper_ci), alpha = 0.2) +
  geom_line(aes(y = mean_per_capita_rate_median), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    y = "Mean per-capita rate of change", x = "Year",
  ) +
  theme_minimal()

# Plot 3: SD of Per-Capita Rates Comparison
plot3_SD_percap_rof <- ggplot(final_indicators_GS, aes(x = year)) +
  geom_ribbon(aes(ymin = sd_per_capita_rate_lower_ci, ymax = sd_per_capita_rate_upper_ci), alpha = 0.2) +
  geom_line(aes(y = sd_per_capita_rate_median), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    y = "SD of per-capita rates of change", x = "Year",
  ) +
  theme_minimal()

# Print the comparison plots
print(plot1_compare)
print(plot2_compare)
print(plot3_compare)

#-----------------------------------------------------------------------------
# STEP 6: PLOT INDIVIDUAL SPECIES TRENDS FROM THE 'GS' MODEL
#-----------------------------------------------------------------------------
# Get predictions from the GS model
preds_GS <- predict(gam_model_GS, newdata = predict_data, type = "response", se.fit = TRUE)
plot_data_GS <- predict_data %>%
  mutate(
    fit = preds_GS$fit,
    se = preds_GS$se.fit,
    lower_ci = fit - 1.96 * se,
    upper_ci = fit + 1.96 * se
  ) %>%
  left_join(community_ts_subset, by = c("year", "species"))

# Plot trends from GS model
species_trends_plot_GS <- ggplot(plot_data_GS, aes(x = year)) +
  geom_point(aes(y = abundance), color = "grey60", alpha = 0.8) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkorange", alpha = 0.3) +
  geom_line(aes(y = fit), color = "darkorange", size = 1) +
  facet_wrap(~ species, scales = "free_y") +
  labs(
    title = "Fitted Abundance Trends for Each Species (GS Model)",
    y = "Predicted Abundance", x = "Year"
  ) +
  theme_minimal()

print(species_trends_plot_GS)