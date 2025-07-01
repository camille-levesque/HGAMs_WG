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


# Chains finished between 67.2 and 77.9 seconds
# ERROR CODE
# Error in !is.null(csv_contents$metadata$save_warmup) && csv_contents$metadata$save_warmup : 
#   invalid 'y' type in 'x && y'
# SOLUTION SHOULD BE TO UPDATE CMDSTANR





# mathematical description of the model