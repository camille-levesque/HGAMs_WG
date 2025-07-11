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
mod1 <- mvgam(y ~ s(time, bs = "tp", k = 6) + # smooth on the time
               
               # Hierarchical intercepts capture variation in average
               # relative abundances
               s(series, bs = 're'),
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
             chains = 2,
             burnin = 500,
             samples = 2000,
             # run_model = FALSE, # because we only want to inspect the priors
             # prior_simulation = TRUE,
             # cmdstanr is highly recommended over rstan as 
             # it is much more up-to-date with the latest 
             # features in the Stan universe
             backend = 'cmdstanr')

# not working?
# mod$save_object("mod1.rds")

# to save mod as .rds if mod$save_object doesn't work
saveRDS(mod1, paste0("prediction/output/mvgam_prediction_mod1.rds")) 


# For now decided to stick directly to Katherine's model so i can get the predictions to show on the graph and then we'll fix whatever might need fixing
mod <- mvgam(data = data_train,
            formula =  y ~ s(time, bs = "tp", k = 6) + 
              s(series, bs = "re"),
            use_lv = TRUE,
            family = "poisson",
            trend_model = 'AR1',
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

# Inspect the model summary
summary(mod)

## how_to_cite(mod)

# Sampling diagnostics (see ?mcmc_plot.mvgam for details on the types
# of {bayesplot} plots that can be used with {mvgam})
mcmc_plot(mod, type = 'rhat_hist')
mcmc_plot(mod, type = 'trace')

# Pairs plots are also useful for diagnosing non-identifiabilities in
# Bayesian models. See ?bayesplot::mcmc_pairs for details. Here a pairs plot
# of the random effect mean and SD parameters shows no worrisome 'funnel' 
# behaviour that can plague hierarchical models:
pairs(mod, variable = c('mean(series)', 'sd(series)'))

# # Plot the hierarchical smooth components with the S3 'plot' function
# # (see ?plot.mvgam for more details)
# plot(mod, type = 'smooths') # No terms to plot-nothing for plot.mvgam() to do


# Plot the hierarchical intercepts
plot(mod, type = 're')

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
plot_predictions(mod, 
                 newdata = data_test,
                 by = c('time', 'series', 'series'), # by is for predictive trends (marginal conditions)
                 points = 0.5) + # transparency
  geom_vline(xintercept = max(data_train$time),
             linetype = 'dashed')+ # adding a line to emphasize the switch from training to testing
  geom_point(data=data_test, aes(x=time, y=y), alpha=.5)+ # adding the true values for predicted data
  theme(legend.position = 'none') +
  labs(y = 'Abundance', x = 'Time') +
  xlim(c(0,30))
ggsave("prediction/figures/A_TrendsFuturePredictions.jpeg", width = 15, height=10)

# A second plot where we see the trend of the training data but only the true points for the predicted ones (so not as interesting)
  plot_predictions(mod, 
                   newdata = data_test,
                   condition = c('time', 'series', 'series'), # conditional predictions which truly means what the trend on training data
                   points = 0.5) +
  geom_vline(xintercept = max(data_train$time),
             linetype = 'dashed')+
  geom_point(data=data_test, aes(x=time, y=y), alpha=.5)+
  theme(legend.position = 'none') +
  labs(y = 'Abundance', x = 'Time')+
    xlim(c(0,30))
  ggsave("prediction/figures/B_TrendsModel.jpeg", width = 15, height=10)
  
# I'd like to get both the trend fitted for the training data and the predictions (might be overkill though)
  # start by fixing the y limits for the two graphs so that we can pretend that it's one graph
  max_temp <- plyr::round_any(max(plot_predictions(mod, 
                                                   newdata = data_test,
                                                   by = c('time', 'series'),
                                                   draw = FALSE)$conf.high), f=ceiling, accuracy=10)
  # matching two graphs into being 'one' but really it's two
  plot_predictions(mod, 
                   by = c('time', 'series'),
                   points = 0.5) +
    theme(legend.position = 'none') +
    labs(y = 'Abundance', x = 'Time') + 
    ylim(c(0,max_temp)) +
    plot_predictions(mod, 
                     newdata = data_test,
                     by = c('time', 'series'),
                     draw = TRUE)+
    theme(legend.position = 'none') +
    geom_point(data=data_test, aes(x=time, y=y, color=series), alpha=.5) +
    ylim(c(0,max_temp)) +
    labs(y = '', x = 'Time') + 
    patchwork::plot_layout(widths=c(.7,.3))
  ggsave("prediction/figures/C_TrendsBoth_nofacet.jpeg", width = 10, height=7)
  
  
  # problem that we're loosing the facet
  # too many species to see everything at once with a facet added but could be an option if 2-3 species
  # /!\ the scales are free right now, it would be good to do a fixed one per species
  plot_predictions(mod, 
                   by = c('time', 'series'),
                   points = 0.5) +
    theme(legend.position = 'none') +
    labs(y = 'Abundance', x = 'Time') + 
    facet_grid(series~., scales="free_y") +
    plot_predictions(mod, 
                     newdata = data_test,
                     by = c('time', 'series'),
                     draw = TRUE)+
    theme(legend.position = 'none') +
    geom_point(data=data_test, aes(x=time, y=y), alpha=.5) +
    facet_grid(series~., scales="free_y") +
    labs(y = '', x = 'Time')+ 
    patchwork::plot_layout(widths=c(.7,.3))
  ggsave("prediction/figures/D_TrendsBoth_facet.jpeg", width = 10, height=30)
  




  plot_predictions(mod, 
                   by = c('time', 'series', 'series'),
                   points = 0.5) +
    theme(legend.position = 'none') +
    labs(y = 'Abundance', x = 'Time') + 
    ylim(c(0,max_temp)) +
    plot_predictions(mod, 
                     newdata = data_test,
                     by = c('time', 'series', 'series'),
                     draw = TRUE)+
    theme(legend.position = 'none') +
    geom_point(data=data_test, aes(x=time, y=y, color=series), alpha=.5) +
    ylim(c(0,max_temp)) +
    labs(y = '', x = 'Time') + 
    patchwork::plot_layout(widths=c(.7,.3))


# Plot for new species predicted
  
## mod without the "new species"
data_noAp$series <- droplevels(data_noAp$series)
  
mod_noAp <- mvgam(data = data_noAp,
                    formula =  y ~ s(time, bs = "tp", k = 5) + 
                      s(series, bs = "re"),
                    use_lv = TRUE,
                    family = "poisson",
                    trend_model = 'AR1',
                    use_stan = TRUE,
                    chains = 2, 
                    burnin = 500,
                    samples = 2000)
  
  
saveRDS(mod_noAp, paste0("prediction/output/mvgam_prediction_mod_noAp.rds")) 
  
mod_noAp <- read_rds("prediction/output/mvgam_prediction_mod_noAp.rds")

## plot the predictions for Agelaius phoeniceus
plot_predictions(mod_noAp, 
                   newdata = data_Ap,
                   by = c('time', 'series', 'series'), # by is for predictive trends (marginal conditions)
                   points = 0.5) + # transparency
    geom_point(data=data_Ap, aes(x=time, y=y), alpha=.5)+ # adding the true values for predicted data
    theme(legend.position = 'none') +
    labs(y = 'Abundance', x = 'Time') +
    xlim(c(0,30))
ggsave("prediction/figures/E_TrendsNewSpecies.jpeg", width = 15, height=10)
  
  
  
  
  
  
  
  
  
  
  
  
