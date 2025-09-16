#-----------------------------------------------------------------------------
# STEP 1: PACKAGES AND LIBRARIES
#-----------------------------------------------------------------------------

install.packages(c("mgcv", "dplyr", "ggplot2", "tidyr", "mvtnorm"))

library(mgcv)
library(dplyr)
library(ggplot2)
library(tidyr)
library(mvtnorm)

#-----------------------------------------------------------------------------
# STEP 2: LOAD AND PRE-PROCESS THE DATA
#-----------------------------------------------------------------------------
# Loading the dataset
data_195 <- read.csv("data_195.csv")

# The Pedersen paper uses aggregated biomass data. This dataset has individual observations.
# Aggregating biomass data to get a yearly abundance for each species.

# Ensure YEAR and valid_name are not NA
community_ts <- data_195 %>%
  filter(!is.na(YEAR) & !is.na(valid_name) & valid_name != "") %>%
  group_by(YEAR, valid_name) %>%
  summarise(ABUNDANCE = n(), .groups = 'drop') %>%
  rename(year = YEAR, species = valid_name, abundance = ABUNDANCE)

# Selecting the top 30 most frequently observed species to make it manageable,
# similar to how the Pedersen paper selected a subset of species.
top_species <- community_ts %>%
  group_by(species) %>%
  summarise(total_abundance = sum(abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 30) %>%
  pull(species)

community_ts_subset <- community_ts %>%
  filter(species %in% top_species) %>%
# The model needs a complete time series for each species.
# Fill in the missing years with 0 abundance.
  complete(species, year = full_seq(year, 1), fill = list(abundance = 0)) %>%
  mutate(species = as.factor(species))

print(head(community_ts_subset))

#-----------------------------------------------------------------------------
# STEP 3: FITTING THE HIERARCHICAL GENERALIZED ADDITIVE MODEL (HGAM)
#-----------------------------------------------------------------------------
# GAM to model the abundance of each species over time.
# This is a "Hierarchical" GAM because we are fitting a model where the smooths for each species are related (they share smoothing parameters).
# 
# Formula:
#   abundance ~ species + s(year, by = species, bs = "tp")
# - 'abundance': Our response variable.
# - 'species': A fixed effect for the mean abundance of each species.
# - 's(year, by = species)': It fits a separate smooth temporal trend (s(year)) for each level of the 'species' factor.
# - 'bs = "tp"' uses a thin plate regression spline.
#
# Family:
#   The paper Pederson et al., 2020 uses a Tweedie distribution, which is good for data with zeros and continuous positive values (like biomass).
#   Here, abundance data is counts, we could use Poisson or Negative Binomial, but Tweedie is still a robust choice.

gam_model <- gam(
  abundance ~ species + s(year, by = species, bs = "tp"),
  data = community_ts_subset,
  family = tw(a = 1.2, b = 1.8), # Tweedie family, parameters chosen for count-like data
  method = "REML" # Recommended method for estimating smoothness
)

summary(gam_model)

#-----------------------------------------------------------------------------
# STEP 4: DERIVATIVES AND COMMUNITY INDICATORS
#-----------------------------------------------------------------------------
# 1. Predict abundance from the model for a range of years.
# 2. Simulate from the model's posterior to get uncertainty estimates.
# 3. Numerically calculate the first derivative (rate of change) for each simulation.
# 4. Calculate the three community-level indicators.

# Create a prediction dataset
predict_data <- community_ts_subset %>%
  select(year, species) %>%
  distinct()

# Define a small number 'eps' for numerical differentiation
eps <- 1e-7

# Predict at year + eps and year - eps to calculate the central difference
predict_data_p_eps <- predict_data %>% mutate(year = year + eps)
predict_data_m_eps <- predict_data %>% mutate(year = year - eps)

# Generate posterior simulations of the linear predictor (lp)
# n = 250 simulations is a good balance of speed and accuracy
n_sim <- 250
set.seed(42) # for reproducibility
sim_lp <- predict(gam_model, newdata = predict_data, type = "lpmatrix")
sim_coef <- rmvnorm(n_sim, coef(gam_model), vcov(gam_model, unconditional = TRUE))

# Calculate predicted values on the response scale for the original data
# and for the epsilon-shifted data
pred_original <- exp(sim_lp %*% t(sim_coef))
pred_p_eps <- exp(predict(gam_model, newdata = predict_data_p_eps, type = "lpmatrix") %*% t(sim_coef))
pred_m_eps <- exp(predict(gam_model, newdata = predict_data_m_eps, type = "lpmatrix") %*% t(sim_coef))

# Calculate the first derivative (rate of change) for each simulation
# Formula: f'(x) ≈ (f(x + ε) - f(x - ε)) / (2ε)
first_derivative <- (pred_p_eps - pred_m_eps) / (2 * eps)

# Calculate per-capita rate of change (r)
# r = (1/A) * dA/dt. We use the predicted abundance for A.
# To avoid division by zero, we add a small constant to the denominator.
per_capita_rate <- first_derivative / (pred_original + 1e-9)

# Now, organize the results into a tidy data frame
results_df <- predict_data %>%
  mutate(pred_abundance = rowMeans(pred_original),
         pred_derivative = rowMeans(first_derivative),
         pred_per_capita = rowMeans(per_capita_rate))

# Reshape the simulation matrices into long format for aggregation
sim_deriv_long <- as.data.frame(first_derivative) %>%
  mutate(row = 1:n()) %>%
  pivot_longer(-row, names_to = "sim_id", values_to = "derivative")

sim_per_capita_long <- as.data.frame(per_capita_rate) %>%
  mutate(row = 1:n()) %>%
  pivot_longer(-row, names_to = "sim_id", values_to = "per_capita_rate")

# Combine with prediction data info
sim_results <- predict_data %>%
  mutate(row = 1:n()) %>%
  left_join(sim_deriv_long, by = "row") %>%
  left_join(sim_per_capita_long, by = c("row", "sim_id"))

# Calculate the three community indicators for each simulation and year
community_indicators <- sim_results %>%
  group_by(year, sim_id) %>%
  summarise(
    # Indicator 1: Mean rate of change
    mean_rate_of_change = mean(derivative, na.rm = TRUE),
    # Indicator 2: Mean per-capita rate of change
    mean_per_capita_rate = mean(per_capita_rate, na.rm = TRUE),
    # Indicator 3: Standard deviation of per-capita rates of change
    sd_per_capita_rate = sd(per_capita_rate, na.rm = TRUE),
    .groups = 'drop'
  )

# Finally, summarize across all simulations to get the median and credible intervals
final_indicators <- community_indicators %>%
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
  )

print(head(final_indicators))


#-----------------------------------------------------------------------------
# STEP 5: PLOT THE INDICATORS
#-----------------------------------------------------------------------------

# Plot 1: Mean Rate of Change
# Periods where the ribbon doesn't cross zero indicate significant community-wide abundance change.
plot1 <- ggplot(final_indicators, aes(x = year)) +
  geom_ribbon(aes(ymin = mean_rate_of_change_lower_ci, ymax = mean_rate_of_change_upper_ci), fill = "grey80") +
  geom_line(aes(y = mean_rate_of_change_median), color = "black", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Indicator 1: Mean Rate of Change",
    y = "Mean Rate of Change (Abundance / Year)",
    x = "Year"
  ) +
  theme_minimal()

# Plot 2: Mean Per-Capita Rate of Change
# Indicates change in overall community size, but less influenced by dominant species
plot2 <- ggplot(final_indicators, aes(x = year)) +
  geom_ribbon(aes(ymin = mean_per_capita_rate_lower_ci, ymax = mean_per_capita_rate_upper_ci), fill = "grey80") +
  geom_line(aes(y = mean_per_capita_rate_median), color = "black", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Indicator 2: Mean Per-Capita Rate of Change",
    y = "Mean Per-Capita Rate (Year⁻¹)",
    x = "Year"
  ) +
  theme_minimal()

# Plot 3: Standard Deviation of Per-Capita Rates
# Elevated periods indicate high compositional turnover (species replacing each other)
plot3 <- ggplot(final_indicators, aes(x = year)) +
  geom_ribbon(aes(ymin = sd_per_capita_rate_lower_ci, ymax = sd_per_capita_rate_upper_ci), fill = "grey80") +
  geom_line(aes(y = sd_per_capita_rate_median), color = "black", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Indicator 3: SD of Per-Capita Rates of Change",
    y = "SD of Per-Capita Rates (Year⁻¹)",
    x = "Year"
  ) +
  theme_minimal()

# Print the plots
print(plot1)
print(plot2)
print(plot3)

#-----------------------------------------------------------------------------
# STEP 6: PLOT INDIVIDUAL SPECIES TRENDS FROM THE HGAM
#-----------------------------------------------------------------------------
# This section visualizes the fitted smooth trend for each species.

# Get predictions and standard errors on the response scale
model_preds <- predict(gam_model, newdata = predict_data, type = "response", se.fit = TRUE)

# Combine predictions with the data for plotting
plot_data <- predict_data %>%
  mutate(
    fit = model_preds$fit,
    se = model_preds$se.fit,
    lower_ci = fit - 1.96 * se,
    upper_ci = fit + 1.96 * se
  ) %>%
  # Add the original abundance data for context
  left_join(community_ts_subset, by = c("year", "species"))

# Create the plot with facets for each species
species_trends_plot <- ggplot(plot_data, aes(x = year)) +
  # Add the original data points
  geom_point(aes(y = abundance), color = "grey60", alpha = 0.8) +
  # Add the confidence interval ribbon
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "steelblue", alpha = 0.3) +
  # Add the fitted trend line
  geom_line(aes(y = fit), color = "steelblue", size = 1) +
  # Create a separate plot for each species
  # scales = "free_y" allows each plot to have its own y-axis range
  facet_wrap(~ species, scales = "free_y") +
  labs(
    title = "Fitted Abundance Trends for Each Species",
    y = "Predicted Abundance",
    x = "Year"
  ) +
  theme_minimal()

# Print the new plot
print(species_trends_plot)

#-----------------------------------------------------------------------------
# STEP 7: FIT AND PLOT A GLOBAL COMMUNITY MODEL (NOT BY SPECIES)
#-----------------------------------------------------------------------------
# This section fits a single GAM to the total community abundance over time.

# First, create a global (total community) time series
community_ts_global <- community_ts_subset %>%
  group_by(year) %>%
  summarise(total_abundance = sum(abundance), .groups = 'drop')

# Fit a simpler GAM to this global data
gam_model_global <- gam(
  total_abundance ~ s(year, bs = "tp"),
  data = community_ts_global,
  family = tw(a = 1.2, b = 1.8),
  method = "REML"
)

# Get predictions from the global model
global_predict_data <- data.frame(year = unique(community_ts_global$year))
global_preds <- predict(gam_model_global, newdata = global_predict_data, type = "response", se.fit = TRUE)

# Combine predictions with the global data for plotting
global_plot_data <- global_predict_data %>%
  mutate(
    fit = global_preds$fit,
    se = global_preds$se.fit,
    lower_ci = fit - 1.96 * se,
    upper_ci = fit + 1.96 * se
  ) %>%
  left_join(community_ts_global, by = "year")

# Create the plot for the global trend
global_trend_plot <- ggplot(global_plot_data, aes(x = year)) +
  geom_point(aes(y = total_abundance), color = "black", alpha = 0.8) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "grey", alpha = 0.3) +
  geom_line(aes(y = fit), color = "grey", size = 1) +
  labs(
    title = "Fitted Global Abundance Trend (All Species Combined)",
    y = "Total Predicted Abundance",
    x = "Year"
  ) +
  theme_minimal()

# Print the global trend plot
print(global_trend_plot)



