# Load packages
library(mgcv)            # Fit GAMs
library(tidyverse)       # Tidy and flexible data manipulation
library(ggplot2)         # Flexible plotting
library(marginaleffects) # Prediction-based inference


# Set up plotting environment
theme_set(
  theme_bw()
)
myhist = function(...){
  geom_histogram(
    col = 'white',
    fill = '#B97C7C', 
    ...
  )
}
hist_theme = function(){
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )
}

# Load the 'aphids' data from the ecostats package. In each of two fields 
# (one of oats, one of wheat) there were eight plots, four with plastic netting 
# to exclude birds, and four without. Aphid abundance was counted on seven 
# different occasions over the first 38 weeks following netting. The hypothesis
# was that aphid numbers would DECREASE when birds were excluded, because an important 
# food source to tree sparrows is aphid predators, hoverflies and ladybird beetles,
# so presence of birds may be limit the effectiveness of a biological control of aphids.
load(
  url('https://github.com/dwarton/ecostats/raw/main/data/aphids.RData')
)

# Bind the two datasets (experimental observations of aphid abundances
# over time in oat and wheat crop plots, under two experimental treatments)
aphid_dat <- aphids$oat %>%
  mutate(crop = 'oat',
         series = paste0('oat_plot_', Plot)) %>%
  bind_rows(
    aphids$wheat %>%
      mutate(
        crop = 'wheat',
        series = paste0('wheat_plot_', Plot)
      )
  )
# View the data structure
glimpse(aphid_dat)

# Wrangle data to improve variable names and create a
# time_since_treat variable
aphid_dat %>%
  mutate(
    series = as.factor(series),
    crop = as.factor(crop),
    birds_excluded = as.factor(
      ifelse(
        Treatment == 'excluded', 'yes', 'no')
    ),
    time_since_treat = Time
  ) %>%
  janitor::clean_names() %>%
  select(-plot, -logcount, -treatment) %>%
  arrange(series, time) -> aphid_ts

# Plot the data
aphid_ts %>%
  ggplot(., aes(x = time_since_treat,
                y = counts,
                col = birds_excluded)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ series) +
  labs(x = 'Time since treatment',
       y = 'Counts of aphids')

# A hierarchical GAM
mod0 <- gam(
  # Observation formula
  formula = counts ~
    
    # Hierarchical (random) intercepts per individual plot
    s(series, bs = 're') +
    
    # Parametric interaction between crop type and bird exclusion,
    # used to centre the hierarchical smooth effects below
    crop * birds_excluded +
    
    # 'Average' nonlinear effect of time since treatment
    s(time_since_treat,
      k = 7,
      bs = 'bs') +
    
    # 'Deviation' nonlinear effects of time since treatment, where
    # every level of the interaction between crop type and bird exclusion
    # has a different smooth function that can deviate from the 'Average' 
    # function above. We use the 'sz' basis to ensure the deviations sum 
    # to zero, which is good for inference but requires that we centre the
    # effects using the parametric interaction effect above
    s(time_since_treat,
      birds_excluded, 
      crop,
      bs = 'sz', 
      k = 7),
  
  # The Aphid data in 'long' format
  data = aphid_ts,
  
  # A Negative Binomial observation family to capture excess
  # overdispersion
  family = nb()
)
summary(mod0)

# Inspect residuals
appraise(mod0, method = 'simulate')

# Inspect fit against the observed data
plot_predictions(
  mod0,
  by = c('time_since_treat',
         'birds_excluded',
         'series'),
  points = 0.5
)

# Average predictions for each treatment
plot_predictions(
  mod0,
  by = c('time_since_treat',
         'birds_excluded',
         'crop')
)

# Smoother predictions on a fine grid
plot_predictions(
  mod0,
  by = c('time_since_treat',
         'birds_excluded'),
  newdata = datagrid(model = mod0,
                     time_since_treat = 3:40,
                     crop = 'oat',
                     birds_excluded = unique)
) +
  labs(title = 'Oat crops') +
plot_predictions(
  mod0,
  by = c('time_since_treat',
         'birds_excluded'),
  newdata = datagrid(model = mod0,
                     time_since_treat = 30:70,
                     crop = 'wheat',
                     birds_excluded = unique)
) +
  labs(title = 'Wheat crops')

# It seems there is evidence to support the hypothesis that
# aphid numbers would DECREASE when birds were excluded (i.e. counts would be
# lower when birds_excluded == 'yes'); But we can more directly assess this
# using comparisons()
plot_comparisons(
  mod0,
  variables = 'birds_excluded',
  by = c('time_since_treat',
         'birds_excluded'),
  newdata = datagrid(model = mod0,
                     time_since_treat = 30:70,
                     crop = 'oat',
                     birds_excluded = unique)
) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(title = 'Oat crops',
       y = 'Difference (birds included - birds excluded') +
  theme(legend.position = 'none') +
  plot_comparisons(
    mod0,
    variables = 'birds_excluded',
    by = c('time_since_treat',
           'birds_excluded'),
    newdata = datagrid(model = mod0,
                       time_since_treat = 30:70,
                       crop = 'wheat',
                       birds_excluded = unique)
  ) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(title = 'Wheat crops',
       y = '') +
  theme(legend.position = 'none')
