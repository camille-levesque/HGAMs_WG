#---
# Title: "Species Abundance Prediction Models using mvgam"
# Author:
# Date:
# Description: "Hierarchical GAM models for predicting species abundance through time"
# Purpose: "Compare different model approaches and generate predictions for new species"
# Dependencies: "mvgam, marginaleffects, patchwork, cmdstanr"
# Output: "Model objects (.rds) and prediction plots (.jpeg)"
#
# TODO:
#   - Incorporate post-stratification in new species predictions
#   - Adjust mathematical model description for current approach
#   - Fix mod$save_object() method (currently using saveRDS workaround)
#   - Standardize y-axis scales across faceted plots
#   - Test ZMVN vs AR1 trend models systematically
#   - Add model comparison metrics (WAIC, LOO-CV)
#
# FIXME:
#   - mod1 vs mod naming inconsistency
#
#---

# Load previous scripts ---------------------------------------------------
source("prediction/code/01_setup.R") # Import & format data and load libraries

# Model fitting  --------------------------------------------------------

set.seed(2505)


# https://github.com/eco4cast/Statistical-Methods-Seminar-Series/blob/main/clark-dynamic_gams/portal_example.R


data_train <- data_train %>%
  droplevels()
data_test <- data_test %>%
  droplevels()

range(data_test$time)

# Plot series
plot_mvgam_series(
  data = data_train,
  newdata = data_test,
  series = "all"
)

# -------------- Penalized splines forecasts
# Set up a State-Space hierarchical GAM with AR1 dynamics for autocorrelation
# Smooths are penalized splines
mod_forecast_ps <- mvgam(
  data = data_train,
  formula = y ~
    s(series, bs = "re"),

  # Hierarchical smooths of time set up as a
  # State-Space model for sampling efficiency
  trend_formula = ~
    s(time, bs = "tp", k = 6) +
      s(time, trend, bs = "sz", k = 6),
  family = poisson(),

  # AR1 for "residual" autocorrelation
  trend_model = AR(p = 1),
  noncentred = TRUE,
  priors = prior(exponential(2),
    class = sigma
  ),
  backend = "cmdstanr"
)
summary(mod_forecast_ps, include_betas = FALSE)
saveRDS(mod_forecast_ps, "prediction/output/mod_forecast_ps.rds")


# Plot the estimated smooth trends
gratia::draw(mod_forecast, trend_effects = TRUE)
conditional_effects(mod_forecast)

forecast_penalized_splines <- plot_predictions(
  mod_forecast,
  by = c("time", "series", "series"),
  newdata = datagrid(
    time = 1:max(data_test$time),
    series = unique
  ),
  type = "expected"
)

ggsave("prediction/figures/forecast_all_species_penalized_splines.png", forecast_penalized_splines)


# Obviously the splines show high extrapolation uncertainty into the
# test time points, but that is ok as it isn't the focus of this exercise. But if
# we wanted better forecasts, let's use GPs in place of the penalized smooths
# https://ecogambler.netlify.app/blog/autocorrelated-gams/

# Look at some of the AR1 estimates
mcmc_plot(mod_forecast, variable = "ar1", regex = TRUE)


# -------------- Gaussian process forecasts
# Using Gaussian Process in place of penalized smooths to get better forecasts
mod_forecast_GP <- mvgam(
  data = data_train,
  formula = y ~
    s(series, bs = "re"),

  # Hierarchical smooths of time set up as a
  # State-Space model for sampling efficiency
  trend_formula = ~
    gp(time, k = 6) +
      s(time, trend, bs = "sz", k = 6),
  family = poisson(),

  # AR1 for "residual" autocorrelation
  trend_model = AR(p = 1),
  noncentred = TRUE,
  priors = prior(exponential(2),
    class = sigma
  ),
  backend = "cmdstanr"
)
summary(mod_forecast_GP, include_betas = FALSE)
saveRDS(mod_forecast_GP, "prediction/output/mod_forecast_GP.rds")

# Add a vertical line where the train test splits (time = 22)
mod_forecast_GP <- readRDS("prediction/output/mod_forecast_GP.rds")
pred_data <- data.frame(
  time = rep(1:max(data_test$time), each = length(unique(data_train$series))),
  series = rep(unique(data_train$series), times = max(data_test$time))
)

predictions <- predictions(
  mod_forecast_GP,
  newdata = pred_data,
  by = c("time", "series", "series"),
  type = "expected"
)

predictions$forecast <- ifelse(predictions$time > 22, "Forecast", "Fitted")

ggplot(predictions, aes(x = time, y = estimate)) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = series),
    alpha = 0.2
  ) +
  geom_line(
    data = subset(predictions, forecast == "Fitted"),
    aes(color = series), linewidth = 1
  ) +
  geom_line(
    data = subset(predictions, forecast == "Forecast"),
    aes(color = series), linewidth = 1, linetype = "dashed"
  ) +
  geom_vline(xintercept = 22, linetype = "dotted") +
  geom_point(data = data_train, aes(x = time, y = y), alpha = 0.2) +
  geom_point(data = data_test, aes(x = time, y = y), alpha = 0.2) +
  facet_wrap(~series, scales = "free_y") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "italic")
  ) +
  labs(y = "Abundance", x = "Time") +
  theme(axis.title = element_text(size = 14))
ggsave("prediction/figures/forecast_all_species_GP_dotted.png",
  plot = last_plot(),
  width = 10,
  height = 6,
  dpi = 300
)

plot_predictions(
  mod_nick_GP,
  by = c("time"),
  newdata = datagrid(
    time = 1:max(data_test$time),
    series = unique
  ),
  type = "expected"
) +
  geom_vline(xintercept = 22, linetype = "dotted") +
  labs(y = "Global latent trend", x = "Time") +
  theme(axis.title = element_text(size = 14))

ggsave("prediction/figures/forecast_latent_trend_GP.png", plot = last_plot())



# -------------- Gaussian process: predicting new species


# Black-and-white Warbler Post-stratification Model ----
# Create training data excluding Mniotilta varia (Black-and-white Warbler)
data_train_noBAWW <- data_noMniotilta %>%
  droplevels()

# Set up State-Space hierarchical GAM excluding BAWW for post-stratification
set.seed(2505)
data_train_noBAWW %>%
  group_by(series) %>%
  summarise(
    n_obs = n(),
    expected_obs = length(unique(data_train_noBAWW$time)),
    missing_obs = expected_obs - n_obs
  ) %>%
  filter(missing_obs > 0)
unique(data_train_noBAWW$series)


colSums(is.na(data_train_noBAWW))

mod_strat_BAWW <- mvgam(
  data = data_train_noBAWW,
  formula = y ~
    s(series, bs = "re"),

  # Hierarchical smooths of time set up as a
  # State-Space model for sampling efficiency
  trend_formula = ~
    gp(time, k = 3, gr = FALSE),
  family = poisson(),

  # AR1 for "residual" autocorrelation
  trend_model = AR(p = 1),
  noncentred = TRUE,
  priors = prior(exponential(2),
    class = sigma
  ),
  backend = "cmdstanr"
)

# Save the model
saveRDS(mod_strat_BAWW, "prediction/output/mod_strat_BAWW_GP_all_years.rds")
mod_strat_BAWW <- readRDS("prediction/output/mod_strat_BAWW_GP_all_years.rds")
# Model summary
summary(mod_strat_BAWW, include_betas = FALSE)

# Post-stratification for Black-and-white Warbler prediction ----
# predict trends for all species and weight these based
# on their "distance" to the new species
unique_species_noBAWW <- levels(data_train_noBAWW$series)

# Add some weights; here a higher value means that species is "closer" to the
# un-modelled target species. This could for example be the inverse of a phylogenetic
# or functional distance metric

# Initialize all weights to 1
species_weights_BAWW <- rep(1, length(unique_species_noBAWW))
names(species_weights_BAWW) <- unique_species_noBAWW

# Assign higher weight to warblers species
# Weight it 10x higher than other species for strong post-stratification

setophaga_species <- grep("Setophaga", unique_species_noBAWW, value = TRUE)
species_weights_BAWW[setophaga_species] <- 10
species_weights_BAWW["Seiurus aurocapilla"] <- 10


# Generate the prediction grid; here we replicate each species' temporal grid
# a number of times, with the number of replications determined by the weight
# vector above
pred_dat_BAWW <- do.call(
  rbind,
  lapply(seq_along(unique_species_noBAWW), function(sp) {
    sp_name <- unique_species_noBAWW[sp]
    weight <- species_weights_BAWW[sp_name]

    do.call(
      rbind,
      replicate(
        weight,
        data.frame(
          time =
            1:max(data_train_noBAWW$time + 1), # time is indexed starting at 0
          series = sp_name
        ),
        simplify = FALSE
      )
    )
  })
) %>%
  dplyr::mutate(
    series = factor(series, levels = levels(data_train_noBAWW$series))
  )

# Generate post-stratified predictions for Black-and-white Warbler
# Marginalize over "time" to compute weighted average predictions
post_strat_BAWW <- marginaleffects::avg_predictions(
  mod_strat_BAWW,
  newdata = pred_dat_BAWW,
  by = "time",
  type = "expected"
)

# Visualization: Post-stratified BAWW predictions ----
# Plot the post-stratified trend predictions for Black-and-white Warbler
# Compare with actual BAWW data from the original dataset
actual_BAWW_data <- dat %>%
  filter(series == "Mniotilta varia")

plot_BAWW_poststrat <- ggplot(post_strat_BAWW, aes(x = time, y = estimate)) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low),
    colour = NA, fill = "steelblue", alpha = 0.4
  ) +
  geom_line(colour = "steelblue", linewidth = 1.2) +
  geom_point(
    data = actual_BAWW_data,
    aes(x = time, y = y),
    colour = "black", alpha = 0.7, size = 2
  ) +
  theme_classic() +
  labs(
    y = "Abundance (Black-and-white Warbler)",
    x = "Time"
  ) +
  theme(
    plot.title = element_text(size = 12),
    plot.subtitle = element_text(size = 10)
  )

print(plot_BAWW_poststrat)

ggsave("prediction/figures/F_BAWW_PostStratified_GP.jpeg",
  plot = plot_BAWW_poststrat, width = 12, height = 8
)


# Anchored Post-stratified BAWW Model ----
# Use first year of BAWW data to calibrate the intercept of post-stratified predictions

# Get the first year BAWW observation for anchoring
first_year_BAWW <- actual_BAWW_data %>%
  filter(time == 0)

# Find the post-stratified prediction for the same time point
first_year_pred <- post_strat_BAWW %>%
  filter(time == 1)

# Calculate the offset needed to match observed abundance
abundance_offset <- first_year_BAWW$y - first_year_pred$estimate

# Apply offset to all post-stratified predictions
post_strat_BAWW_anchored <- post_strat_BAWW %>%
  mutate(
    estimate_original = estimate,
    conf.low_original = conf.low,
    conf.high_original = conf.high,
    estimate = estimate + abundance_offset,
    conf.low = conf.low + abundance_offset,
    conf.high = conf.high + abundance_offset
  )

mvgam_pred_mean <- post_strat_BAWW_anchored$estimate
mae_mvgam <- mean(abs(data_Mniotilta$y - mvgam_pred_mean))
rmse_mvgam <- sqrt(mean((data_Mniotilta$y - mvgam_pred_mean)^2))


# GAM for BAWW
BAWW_gam <- gam(y ~ s(time), data = data_Mniotilta)


# Calculate correlation between predictions
cor_predictions <- cor(mvgam_pred_mean, gam_pred)
cor_predictions

# Direction of trend agreement
gam_trend_direction <- sign(diff(gam_pred))
mvgam_trend_direction <- sign(diff(mvgam_pred_mean))
trend_agreement <- mean(gam_trend_direction == mvgam_trend_direction)
trend_agreement


# Visualization: Anchored vs Original Post-stratified Predictions ----
plot_BAWW_anchored <- ggplot() +
  # Anchored post-stratified prediction
  geom_ribbon(
    data = post_strat_BAWW_anchored,
    aes(x = time, ymin = conf.low, ymax = conf.high),
    fill = "darkgreen", alpha = 0.4
  ) +
  geom_line(
    data = post_strat_BAWW_anchored,
    aes(x = time, y = estimate),
    color = "darkgreen", size = 1.2, linetype = "dashed"
  ) +

  # Actual BAWW data
  geom_point(
    data = actual_BAWW_data,
    aes(x = time, y = y),
    colour = "black", alpha = 0.2, size = 2.5
  ) +
  theme_classic() +
  labs(
    y =
      expression(paste("Predicted ", italic("Mniotilta varia"), " abundance")),
    x = "Time"
  ) +
  theme(
    axis.title = element_text(size = 12),
    legend.position = "none"
  )


ggsave("prediction/figures/G_BAWW_Anchored_GP_all_years.jpeg",
  plot = plot_BAWW_anchored, width = 6, height = 4
)



# # Model Diagnostics ----
# # Inspect the model summary
# summary(mod)

# ## how_to_cite(mod)

# # Sampling diagnostics (see ?mcmc_plot.mvgam for details on the types
# # of {bayesplot} plots that can be used with {mvgam})
# mcmc_plot(mod, type = "rhat_hist")
# mcmc_plot(mod, type = "trace")

# # Pairs plots are also useful for diagnosing non-identifiabilities in
# # Bayesian models. See ?bayesplot::mcmc_pairs for details. Here a pairs plot
# # of the random effect mean and SD parameters shows no worrisome 'funnel'
# # behaviour that can plague hierarchical models:
# pairs(mod, variable = c("mean(series)", "sd(series)"))

# # # Plot the hierarchical smooth components with the S3 'plot' function
# # # (see ?plot.mvgam for more details)
# # plot(mod, type = 'smooths') # No terms to plot-nothing for plot.mvgam() to do
