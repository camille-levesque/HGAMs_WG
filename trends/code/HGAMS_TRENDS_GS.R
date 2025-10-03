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
# STEP 3: FIT POISSON MODELS (S vs GS) AND COMPARE AIC
#-----------------------------------------------------------------------------
# Here, we are fitting model "GS" from Pedersen et al. (2019) paper.
# As we are working with abundance (count) data, we'll use the Poisson family.

# MODEL 1: The 'S' Model (Separate smooths for each species)
gam_model_S <- gam(
  abundance ~ species +
    s(year, by = species, bs = "fs") + # Species-specific smoother: a factor-smoother interaction of year and species
    s(species, bs = "re"), # Species as random effects (gives an intercept per species)
  data = community_ts_subset,
  family = poisson(), # Using Poisson family
  method = "REML"
)

# MODEL 2: The 'GS' Model (Global smooth + species-specific deviations)
gam_model_GS <- gam(
  abundance ~ s(year, bs = "tp") + # Global relationship
    s(year, by = species, bs = "fs") + # Species-specific smoother: a factor-smoother interaction of year and species
    s(species, bs="re"), # Species as random effects (gives an intercept per species)
  data = community_ts_subset,
  family = poisson(), # Using Poisson family
  method = "REML"
)

# Compare the models using AIC
aic_S <- AIC(gam_model_S)
aic_GS <- AIC(gam_model_GS)

# Print the results for comparison
print(paste("AIC for Poisson S Model:", round(aic_S, 2)))
print(paste("AIC for Poisson GS Model:", round(aic_GS, 2)))
# Best AIC score: Model GS

#-----------------------------------------------------------------------------
# STEP 4: DERIVATIVES AND INDICATORS FOR THE 'S' MODEL
#-----------------------------------------------------------------------------
# This entire section calculates the three community indicators based ONLY on the S model.

# Create a prediction dataset
predict_data <- community_ts_subset %>%
  select(year, species) %>%
  distinct()

# Define a small number 'eps' for numerical differentiation
eps <- 1e-7
predict_data_p_eps <- predict_data %>% mutate(year = year + eps)
predict_data_m_eps <- predict_data %>% mutate(year = year - eps)

# Generate posterior simulations from the S model
n_sim <- 250
set.seed(42)
sim_lp_S <- predict(gam_model_S, newdata = predict_data, type = "lpmatrix")
sim_coef_S <- rmvnorm(n_sim, coef(gam_model_S), vcov(gam_model_S, unconditional = TRUE))

# Calculate predicted values and derivatives for the S model
pred_original_S <- exp(sim_lp_S %*% t(sim_coef_S))
pred_p_eps_S <- exp(predict(gam_model_S, newdata = predict_data_p_eps, type = "lpmatrix") %*% t(sim_coef_S))
pred_m_eps_S <- exp(predict(gam_model_S, newdata = predict_data_m_eps, type = "lpmatrix") %*% t(sim_coef_S))
first_derivative_S <- (pred_p_eps_S - pred_m_eps_S) / (2 * eps)
per_capita_rate_S <- first_derivative_S / (pred_original_S + 1e-9)

# Reshape simulation results for the S model
sim_deriv_long_S <- as.data.frame(first_derivative_S) %>%
  mutate(row = 1:n()) %>%
  pivot_longer(-row, names_to = "sim_id", values_to = "derivative")
sim_per_capita_long_S <- as.data.frame(per_capita_rate_S) %>%
  mutate(row = 1:n()) %>%
  pivot_longer(-row, names_to = "sim_id", values_to = "per_capita_rate")

sim_results_S <- predict_data %>%
  mutate(row = 1:n()) %>%
  left_join(sim_deriv_long_S, by = "row") %>%
  left_join(sim_per_capita_long_S, by = c("row", "sim_id"))

# Calculate community indicators for the S model
community_indicators_S <- sim_results_S %>%
  group_by(year, sim_id) %>%
  summarise(
    mean_rate_of_change = mean(derivative, na.rm = TRUE),
    mean_per_capita_rate = mean(per_capita_rate, na.rm = TRUE),
    sd_per_capita_rate = sd(per_capita_rate, na.rm = TRUE),
    .groups = 'drop'
  )

# Final summary of indicators for the S model
final_indicators_S <- community_indicators_S %>%
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
  mutate(model_type = "S Model") # Add a label for plotting

print("--- Indicators from S Model ---")
print(head(final_indicators_S))


#-----------------------------------------------------------------------------
# STEP 5: DERIVATIVES AND INDICATORS FOR THE 'GS' MODEL
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
# STEP 6: PLOT AND COMPARE INDICATORS FROM BOTH MODELS
#-----------------------------------------------------------------------------
# Combine the indicator results from both models into one dataframe
combined_indicators <- bind_rows(final_indicators_S, final_indicators_GS)

# Plot 1: Mean Rate of Change Comparison
plot1_compare <- ggplot(combined_indicators, aes(x = year, group = model_type)) +
  geom_ribbon(aes(ymin = mean_rate_of_change_lower_ci, ymax = mean_rate_of_change_upper_ci, fill = model_type), alpha = 0.2) +
  geom_line(aes(y = mean_rate_of_change_median, color = model_type), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Comparison: Mean Rate of Change",
    y = "Mean Rate of Change (Abundance / Year)", x = "Year",
    color = "Model Type", fill = "Model Type"
  ) +
  theme_minimal()

# Plot 2: Mean Per-Capita Rate of Change Comparison
plot2_compare <- ggplot(combined_indicators, aes(x = year, group = model_type)) +
  geom_ribbon(aes(ymin = mean_per_capita_rate_lower_ci, ymax = mean_per_capita_rate_upper_ci, fill = model_type), alpha = 0.2) +
  geom_line(aes(y = mean_per_capita_rate_median, color = model_type), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Comparison: Mean Per-Capita Rate of Change",
    y = "Mean Per-Capita Rate (Year⁻¹)", x = "Year",
    color = "Model Type", fill = "Model Type"
  ) +
  theme_minimal()

# Plot 3: SD of Per-Capita Rates Comparison
plot3_compare <- ggplot(combined_indicators, aes(x = year, group = model_type)) +
  geom_ribbon(aes(ymin = sd_per_capita_rate_lower_ci, ymax = sd_per_capita_rate_upper_ci, fill = model_type), alpha = 0.2) +
  geom_line(aes(y = sd_per_capita_rate_median, color = model_type), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Comparison: SD of Per-Capita Rates of Change",
    y = "SD of Per-Capita Rates (Year⁻¹)", x = "Year",
    color = "Model Type", fill = "Model Type"
  ) +
  theme_minimal()

# Print the comparison plots
print(plot1_compare)
print(plot2_compare)
print(plot3_compare)


#-----------------------------------------------------------------------------
# STEP 7: PLOT INDIVIDUAL SPECIES TRENDS FROM THE 'S' MODEL
#-----------------------------------------------------------------------------
# Get predictions from the S model
preds_S <- predict(gam_model_S, newdata = predict_data, type = "response", se.fit = TRUE)
plot_data_S <- predict_data %>%
  mutate(
    fit = preds_S$fit,
    se = preds_S$se.fit,
    lower_ci = fit - 1.96 * se,
    upper_ci = fit + 1.96 * se
  ) %>%
  left_join(community_ts_subset, by = c("year", "species"))

# Plot trends from S model
species_trends_plot_S <- ggplot(plot_data_S, aes(x = year)) +
  geom_point(aes(y = abundance), color = "grey60", alpha = 0.8) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = fit), color = "steelblue", size = 1) +
  facet_wrap(~ species, scales = "free_y") +
  labs(
    title = "Fitted Abundance Trends for Each Species (S Model)",
    y = "Predicted Abundance", x = "Year"
  ) +
  theme_minimal()

print(species_trends_plot_S)


#-----------------------------------------------------------------------------
# STEP 8: PLOT INDIVIDUAL SPECIES TRENDS FROM THE 'GS' MODEL
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