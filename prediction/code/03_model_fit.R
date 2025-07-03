# Load previous scripts ---------------------------------------------------
source("prediction/code/01_setup.R") # Import & format data and load libraries
source("prediction/code/02_data_exploration.R") # Data exploration
# Model fitting  --------------------------------------------------------

# First import of Nicholas' code (starting line 200) to adjust it - we will filter what is relevant later 

# We only want to predict relative abundance of species through time.
# An initial model will attempt to capture variation 
# in species' responses to time series, using a 
# Hierarchical GAM (HGAM) with no dynamic component 
# (see ?mgcv::s, ?mgcv::te, ?brms::gp and ?mvgam::dynamic for 
# more information on the kinds of terms supported in {mvgam}
# formulae). This model takes ~ XXXXXX seconds to fit, after compilation
mod <- mvgam(rel_abun ~ 
               
               # Hierarchical intercepts capture variation in average
               # relative abundances
               s(series, bs = 're'),
               
               # A shared smooth of minimum temperature         # NO NEED FOR US
               ## s(mintemp, k = 8) +                            
               
               # Deviation smooths of minimum temperature,      # NO NEED FOR US
               # allowing each species' response to mintemp to vary
               # from the shared smooth
               ## s(mintemp, series, bs = 'sz', k = 8) - 1,        
             
             # Condition on the training data
             data = data_train,
             
             # Automatically compute forecasts for the test data
             newdata = data_test,
             
             # Beta observations with independent precisions
             family = betar(),
             
             # cmdstanr is highly recommended over rstan as 
             # it is much more up-to-date with the latest 
             # features in the Stan universe
             backend = 'cmdstanr')

# Warning messages:
#   1: There were 4 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
# https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 2: Examine the pairs() plot to diagnose sampling problems
# 
# 3: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#bulk-ess 
# 4: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#tail-ess 

mod$save_object("mod.rds")


# Chains finished between 67.2 and 77.9 seconds
# ERROR CODE
# Error in !is.null(csv_contents$metadata$save_warmup) && csv_contents$metadata$save_warmup : 
#   invalid 'y' type in 'x && y'
# SOLUTION SHOULD BE TO UPDATE CMDSTANR

# # we recommend running this in a fresh R session or restarting your current session
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# After updating cmdstanR, Rtools was  an issue, fixed by downloading RTools 4.4

# The model now runs in 153.548

# Warning message:
#   Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
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
plot_predictions(mod, 
                 condition = c('time', 'rel_abun', 'series'),
                 points = 0.5,
                 rug = TRUE) +
  theme(legend.position = 'none') +
  labs(y = 'Relative abundance', x = 'Time')























