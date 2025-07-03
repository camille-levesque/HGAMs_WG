# Load previous scripts ---------------------------------------------------
source("prediction/code/01_setup.R") # Import & format data and load libraries
source("prediction/code/02_data_exploration.R") # Data exploration
# Model fitting  --------------------------------------------------------

set.seed(2505)

# First import of Nicholas' code (starting line 200) to adjust it - we will filter what is relevant later 

# We only want to predict relative abundance of species through time.
# An initial model will attempt to capture variation 
# in species' responses to time series, using a 
# Hierarchical GAM (HGAM) with no dynamic component 
# (see ?mgcv::s, ?mgcv::te, ?brms::gp and ?mvgam::dynamic for 
# more information on the kinds of terms supported in {mvgam}
# formulae). This model takes ~ XXXXXX seconds to fit, after compilation
mod <- mvgam(y ~ s(time, bs = "tp", k = 6) + # smooth on the time
               
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
mod$save_object("mod.rds")

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

# now we are actually predicting but we don't have the points for the true values yet (still a first win!)
plot_predictions(mod, 
                 newdata = data_test,
                 by = c('time', 'series', 'series'),
                 points = 0.5) +
  geom_vline(xintercept = max(data_train$time),
             linetype = 'dashed')+
  geom_point(data=data_test, aes(x=time, y=y), alpha=.5)+
    theme(legend.position = 'none') +
    labs(y = 'Abundance', x = 'Time')


for (i in 1:length(unique(data_test$series))) {
  plot(mod, type ="forecast", series =i)
  
}





















