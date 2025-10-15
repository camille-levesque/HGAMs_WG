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
#   - Hardcoded k=6 and k=5 values need justification
#
#---

# Load previous scripts ---------------------------------------------------
source("prediction/code/01_setup.R") # Import & format data and load libraries
source("prediction/code/02_data_exploration.R") # Data exploration

# Model fitting  --------------------------------------------------------

set.seed(2505)

# First import of Nicholas' code (starting line 200) to adjust it - we will filter what is relevant later
# https://github.com/eco4cast/Statistical-Methods-Seminar-Series/blob/main/clark-dynamic_gams/portal_example.R

# We only want to predict relative abundance of species through time.
# An initial model will attempt to capture variation
# in species' responses to time series, using a
# Hierarchical GAM (HGAM) with no dynamic component
# (see ?mgcv::s, ?mgcv::te, ?brms::gp and ?mvgam::dynamic for
# more information on the kinds of terms supported in {mvgam}
# formulae). This model takes ~ XXXXXX seconds to fit, after compilation


# Primary Model Fitting ----

mod1 <- mvgam(
  y ~ s(time, bs = "tp", k = 6) + # smooth on the time

    # Hierarchical intercepts capture variation in average
    # relative abundances
    s(series, bs = "re"),
  use_lv = TRUE,

  # Condition on the training data
  data = data_train,
  # Automatically compute forecasts for the test data
  newdata = data_test,

  # Beta observations with independent precisions
  family = poisson(),
  trend_model = ZMVN(), # might need to change to RW or ZMVN
  #
  # use_stan = TRUE,
  chains = 4,
  burnin = 500,
  samples = 2000,
  # run_model = FALSE, # because we only want to inspect the priors
  # prior_simulation = TRUE,
  # cmdstanr is highly recommended over rstan as
  # it is much more up-to-date with the latest
  # features in the Stan universe
  backend = "cmdstanr"
)

# not working?
# mod$save_object("mod1.rds")

# to save mod as .rds if mod$save_object doesn't work
saveRDS(mod1, paste0("prediction/output/mvgam_prediction_mod1.rds"))

# Alternative Model (Katherine's approach) ----
# For now decided to stick directly to Katherine's model so i can get the predictions to show on the graph and then we'll fix whatever might need fixing
mod <- mvgam(
  data = data_train,
  formula = y ~ s(time, bs = "tp", k = 6) +
    s(series, bs = "re"),
  use_lv = TRUE,
  family = "poisson",
  trend_model = "AR1",
  use_stan = TRUE,
  chains = 2,
  burnin = 500,
  samples = 2000
)

saveRDS(mod, paste0("prediction/output/mvgam_prediction_mod.rds"))


# https://mc-stan.org/misc/warnings.html#bulk-ess

# mathematical description of the model:
# The model can be described mathematically as follows: ------ TO BE ADJUSTED

#                   for s in 1:N_species ...
#                  for t in 1:N_timepoints...

#                   ## Observation model ##
#  Relative abundance[s, t] ~ Beta(μ[s, t], φ[s])

#                   ## Linear predictor ##           ------ WE TAKE OUT ? (OU JUSTE PARTIE AVEC MINTEMP)
#            logit(μ[s, t]) = α[s] + f(mintemp)_shared[t] +
#                             f(mintemp)_species[s, t]
#                         f = sum(β_smooth * b)

#                      ## Priors ##
#                         α ~ Normal(μ_population, σ_population)
#              μ_population ~ Normal(0, 1)
#              σ_population ~ Student-t(3, 0, 2.5)[0, ] ------ TO BE ADJUSTED?
#                  β_smooth ~ MVNormal(0, (Ω ∗ λ)^−1)
#                         λ ~ Normal(5, 30)[0, ]
#                         φ ~ Gamma(0.01, 0.01)

# where:  ------------------------------------- TO BE ADJUSTED
# f are the penalized smooth functions
# b are thin plate spline basis functions
# Ω are a set of prior precision matrices for the smooths
# λ are the smoothing penalties that prevent overfitting; note that
#   Normal(5, 30)[0, ] indicates a half-normal prior distribution
# φ are species-level precision parameters


# If you would like to see the underlying Stan code, which is fully
# transparent in its use of prior distributions, use the code()
# function:
code(mod)

# Model Diagnostics ----
# Inspect the model summary
summary(mod)

## how_to_cite(mod)

# Sampling diagnostics (see ?mcmc_plot.mvgam for details on the types
# of {bayesplot} plots that can be used with {mvgam})
mcmc_plot(mod, type = "rhat_hist")
mcmc_plot(mod, type = "trace")

# Pairs plots are also useful for diagnosing non-identifiabilities in
# Bayesian models. See ?bayesplot::mcmc_pairs for details. Here a pairs plot
# of the random effect mean and SD parameters shows no worrisome 'funnel'
# behaviour that can plague hierarchical models:
pairs(mod, variable = c("mean(series)", "sd(series)"))

# # Plot the hierarchical smooth components with the S3 'plot' function
# # (see ?plot.mvgam for more details)
# plot(mod, type = 'smooths') # No terms to plot-nothing for plot.mvgam() to do


# Plot the hierarchical intercepts
plot(mod, type = "re")

# More informative plots can be made using plot_predictions() from
# the {marginaleffects} universe to visualise conditional effects
# on the outcome scale. See ?marginaleffects::plot_predictions
# for details, or visit: https://marginaleffects.com/
# plot_predictions(mod,
#                  condition = c('time', 'series', 'series'),
#                  type = "link",
#                  newdata = data_test,
#                  points = 0.5,
#                  rug = TRUE) +
#   theme(legend.position = 'none') +
#   labs(y = 'Abundance', x = 'Time')

# A first plot to show our training data and the trends for predicted data of the model (with true values in black)

# Visualization: Future Predictions ----
# Example 3: Future states -- Extrapolation

plot_predictions(mod,
  newdata = data_test,
  by = c("time", "series", "series"), # by is for predictive trends (marginal conditions)
  points = 0.5
) + # transparency
  geom_vline(
    xintercept = max(data_train$time),
    linetype = "dashed"
  ) + # adding a line to emphasize the switch from training to testing
  geom_point(data = data_test, aes(x = time, y = y), alpha = .5) + # adding the true values for predicted data
  theme(legend.position = "none") +
  labs(y = "Abundance", x = "Time") +
  xlim(c(0, 30))
ggsave("prediction/figures/A_TrendsFuturePredictions.jpeg", width = 15, height = 10)

# Visualization: Training Trends ----
# A second plot where we see the trend of the training data but only the true points for the predicted ones (so not as interesting)
plot_predictions(mod,
  newdata = data_test,
  condition = c("time", "series", "series"), # conditional predictions which truly means what the trend on training data
  points = 0.5
) +
  geom_vline(
    xintercept = max(data_train$time),
    linetype = "dashed"
  ) +
  geom_point(data = data_test, aes(x = time, y = y), alpha = .5) +
  theme(legend.position = "none") +
  labs(y = "Abundance", x = "Time") +
  xlim(c(0, 30))
ggsave("prediction/figures/B_TrendsModel.jpeg", width = 15, height = 10)

# I'd like to get both the trend fitted for the training data and the predictions (might be overkill though)
# start by fixing the y limits for the two graphs so that we can pretend that it's one graph
max_temp <- plyr::round_any(max(plot_predictions(mod,
  newdata = data_test,
  by = c("time", "series"),
  draw = FALSE
)$conf.high), f = ceiling, accuracy = 10)


# Visualization: Combined Plots ----
# matching two graphs into being 'one' but really it's two
plot_predictions(mod,
  by = c("time", "series"),
  points = 0.5
) +
  theme(legend.position = "none") +
  labs(y = "Abundance", x = "Time") +
  ylim(c(0, max_temp)) +
  plot_predictions(mod,
    newdata = data_test,
    by = c("time", "series"),
    draw = TRUE
  ) +
  theme(legend.position = "none") +
  geom_point(data = data_test, aes(x = time, y = y, color = series), alpha = .5) +
  ylim(c(0, max_temp)) +
  labs(y = "", x = "Time") +
  patchwork::plot_layout(widths = c(.7, .3))
ggsave("prediction/figures/C_TrendsBoth_nofacet.jpeg", width = 10, height = 7)


# problem that we're loosing the facet
# too many species to see everything at once with a facet added but could be an option if 2-3 species
# /!\ the scales are free right now, it would be good to do a fixed one per species
plot_predictions(mod,
  by = c("time", "series"),
  points = 0.5
) +
  theme(legend.position = "none") +
  labs(y = "Abundance", x = "Time") +
  facet_grid(series ~ ., scales = "free_y") +
  plot_predictions(mod,
    newdata = data_test,
    by = c("time", "series"),
    draw = TRUE
  ) +
  theme(legend.position = "none") +
  geom_point(data = data_test, aes(x = time, y = y), alpha = .5) +
  facet_grid(series ~ ., scales = "free_y") +
  labs(y = "", x = "Time") +
  patchwork::plot_layout(widths = c(.7, .3))
ggsave("prediction/figures/D_TrendsBoth_facet.jpeg", width = 10, height = 30)





plot_predictions(mod,
  by = c("time", "series", "series"),
  points = 0.5
) +
  theme(legend.position = "none") +
  labs(y = "Abundance", x = "Time") +
  ylim(c(0, max_temp)) +
  plot_predictions(mod,
    newdata = data_test,
    by = c("time", "series", "series"),
    draw = TRUE
  ) +
  theme(legend.position = "none") +
  geom_point(data = data_test, aes(x = time, y = y, color = series), alpha = .5) +
  ylim(c(0, max_temp)) +
  labs(y = "", x = "Time") +
  patchwork::plot_layout(widths = c(.7, .3))

# New Species Analysis ----
# Example 2: new species
# Predict for Mniotilta varia

# Post-strat: weigh Setophaga coronata highest
## mod without the "new species"
data_noMniotilta$series <- droplevels(data_noMniotilta$series)

mod_noMniotilta <- mvgam(
  data = data_noMniotilta,
  formula = y ~ s(time, bs = "tp", k = 5) +
    s(series, bs = "re"),
  use_lv = TRUE,
  family = "poisson",
  trend_model = "AR1",
  use_stan = TRUE,
  chains = 2,
  burnin = 500,
  samples = 2000
)


saveRDS(mod_noMniotilta, paste0("prediction/output/mvgam_prediction_mod_noMniotilta.rds"))

mod_noMniotilta <- read_rds("prediction/output/mvgam_prediction_mod_noMniotilta.rds")



## plot the predictions for a new species

plot_predictions(mod_noMniotilta,
  newdata = data_Ap,
  by = c("time", "series", "series"), # by is for predictive trends (marginal conditions)
  points = 0.5
) + # transparency
  geom_point(data = data_Ap, aes(x = time, y = y), alpha = .5) + # adding the true values for predicted data
  theme(legend.position = "none") +
  labs(y = "Abundance", x = "Time") +
  xlim(c(0, 30))
ggsave("prediction/figures/E_TrendsNewSpecies.jpeg", width = 15, height = 10)

# Example 1: create smooth for Setophaga pinus, a data-poor group
#### ALI TESTING 
# I have an issue that one of the species has 3 observations which is messing the model so I'm removing it
data_noMniotilta <- data_noMniotilta[data_noMniotilta$series!="Setophaga pinus",]
data_noMniotilta$series <- droplevels(data_noMniotilta$series)

# fitting a model with a smooth on time and at the species level
mod_noMniotilta <- mvgam(
  data = data_noMniotilta,
  formula = y ~ s(time, bs = "tp", k = 5) +
    s(series, bs = "re"),
  use_lv = TRUE,
  family = "poisson",
  trend_model = "AR1",
  use_stan = TRUE,
  chains = 2,
  burnin = 500,
  samples = 2000
)

# creating a new df that will contain the time and our new species (as a factor)
data_Ap <- data.frame(time=unique(data_noMniotilta$time)
                      , series = "Mniotilta")
data_Ap$series <- factor(data_Ap$series)
# extracting the values for the predictions of our new species
idk=data.frame(predict(mod_noMniotilta,
                       newdata = data_Ap, type="response"))
# adding the rest of the info to the df to help plot it
idk$time=unique(data_noMniotilta$time)
idk$series = "Mniotilta"
# we can then plot it as well as all the other species
plot_predictions(mod_noMniotilta,
                 newdata = data_Ap,
                 by = c("time", "series", "series"), # by is for predictive trends (marginal conditions)
                 points = 0.5
)+ # transparency
  geom_point(data = idk, aes(x = time, y = Estimate), alpha = .5) + # adding the true values for predicted data
  theme(legend.position = "none") +
  labs(y = "Abundance", x = "Time") +
  xlim(c(0, 30))

# Now let's try to use that post stratification thing?
# best way I can figure out how to do it is add a new variable that says whether my species are close to the one we'll try to predict (= another level of group that is higher than the species)
data_noMniotilta$similar.species <- "no"
# setting it to 'no' for all species except our yellow warbler
data_noMniotilta$similar.species[data_noMniotilta$series=="Setophaga coronata"] <- "yes"
data_noMniotilta$similar.species <- factor(data_noMniotilta$similar.species)

# fitting a model with a smooth on time, on our new group and on the species
mod_noMniotilta <- mvgam(
  data = data_noMniotilta,
  formula = y ~ s(time, bs = "tp", k = 5) +
    s(similar.species, bs="re") +
    s(series, bs = "re"),
  use_lv = TRUE,
  family = "poisson",
  trend_model = "AR1",
  use_stan = TRUE,
  chains = 2,
  burnin = 500,
  samples = 2000
)

# creating a df with out new species along with the value for the new group
data_Ap <- data.frame(time=unique(data_noMniotilta$time)
                      , similar.species = "yes"
                      , series = "Mniotilta")
data_Ap$series <- factor(data_Ap$series)
# extracting the values for the predictions of our new species
idk=data.frame(predict(mod_noMniotilta,
        newdata = data_Ap, type="response"))
# adding the rest of the info to the df to help plot it
idk$time=unique(data_noMniotilta$time)
idk$series = "Mniotilta"
idk$similar.species = "yes"

# and now we're plotting all of our species + the new one
# we'll have in different colors the species we used for our extra smoother
# not sure what the difference is between the blue line and the black dots?
# maybe blue line = smoother for the yellow warbler and black lines = predictions based on global smooth + yellow warbler
plot_predictions(mod_noMniotilta,
                 newdata = data_Ap,
                 by = c("time", "similar.species", "series"), # by is for predictive trends (marginal conditions)
                 points = 0.5
)+ # transparency
  geom_point(data = idk, aes(x = time, y = Estimate), alpha = .5) + # adding the true values for predicted data
  theme(legend.position = "none") +
  labs(y = "Abundance", x = "Time") +
  xlim(c(0, 30))

#  /!\ the two results look very similar BUT pay attention that the abundance range changes from 0-125 to 0-150



#### Nick's additions ####
data_train <- data_train %>%
  droplevels()
data_test <- data_test %>%
  droplevels()

# Plot series
plot_mvgam_series(
  data = data_train,
  newdata = data_test,
  series = "all"
)

# Set up a State-Space hierarchical GAM with AR1 dynamics for autocorrelation
mod_nick <- mvgam(
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
                 class = sigma),
  backend = "cmdstanr"
)
summary(mod_nick, include_betas = FALSE)

# Plot the estimated smooth trends
gratia::draw(mod_nick, trend_effects = TRUE)
conditional_effects(mod_nick)
plot_predictions(
  mod_nick,
  by = c("time", "series", "series"),
  newdata = datagrid(time = 1:max(data_test$time),
                     series = unique),
  type = "expected"
)

# Obviously the splines show high extrapolation uncertainty into the 
# test time points, but that is ok as it isn't the focus of this exercise. But if
# we wanted better forecasts, I'd use GPs in place of the penalized smooths
# https://ecogambler.netlify.app/blog/autocorrelated-gams/

# Look at some of the AR1 estimates
mcmc_plot(mod_nick, variable = "ar1", regex = TRUE)

# Post-stratification: predict trends for all species and weight these based
# on their "distance" to the new species
unique_species <- levels(data_train$series)

# Add some weights; here a higher value means that species is "closer" to the
# un-modelled target species. This could for example be the inverse of a phylogenetic
# or functional distance metric
species_weights <- pmax(
  1, 
  rnbinom(
    nlevels(data_train$series), 
    mu = 5, 
    size = 1
  )
)

# Generate the prediction grid; here we replicate each species' temporal grid
# a number of times, with the number of replications determined by the weight
# vector above
pred_dat <- do.call(
  rbind, 
  lapply(seq_along(unique_species), function(sp){
    do.call(
      rbind, 
      replicate(
        species_weights[sp],
        data.frame(time = 1:max(data_test$time),
                   series = unique_species[sp]),
        simplify = FALSE
      )
    )
  })
) %>%
  dplyr::mutate(series = factor(series, levels = levels(data_train$series)))

# Marginalize over "time" to compute the weighted average predictions, accounting
# for full uncertainty in the species-level trends but IGNORING the AR1 process
post_strat_trend <- marginaleffects::avg_predictions(
  mod_nick,
  newdata = pred_dat,
  by = "time",
  type = "expected"
)

# Plot the post-stratified trend expectations
ggplot(post_strat_trend,
       aes(x = time, y = estimate)) +
  geom_ribbon(aes(ymax = conf.high,
                  ymin = conf.low),
              colour = NA,
              fill = "darkred",
              alpha = 0.4) +
  geom_line(colour = "darkred") +
  theme_classic() +
  labs(y = "Post-stratified trend prediction",
       x = "Time")

