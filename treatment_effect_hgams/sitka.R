# Load packages
library(tidyverse); theme_set(theme_classic())
library(mgcv)
library(marginaleffects)
library(MASS)

# Inspect the Sitka data, which contains growth curves for 
# 79 Sitka spruce trees in 1989, 54 of which were grown in 
# ozone-enriched chambers and 25 were controls. The size was 
# measured eight times in 1989, at roughly monthly intervals
data(Sitka)
?Sitka

# Convert the treatment and tree indicators to unordered factors
dat <- Sitka %>%
  mutate(treat = factor(treat),
         tree = factor(tree)) %>%
  janitor::clean_names()
glimpse(dat)

# Plot the curves
ggplot(
  dat,
  aes(x = time,
      y = size,
      group = tree)
) +
  geom_line() +
  facet_wrap(~ treat) +
  labs(y = 'Size (log(height * diameter))',
       x = 'Days since January 1, 1988')

# Fit a hierarchical GAM that allows the growth curves
# to be partially pooled by treatment; use bam()  
# for faster fitting
mod <- bam(
  size ~ 
    s(time, k = 5) +
    s(time, treat, bs = 'fs', k = 5) +
    s(time, tree, bs = 'fs', k = 5),
  data = dat,
  method = 'fREML',
  family = gaussian(link = 'log')
)
summary(mod)

# Draw the component smooths
gratia::draw(mod)

# Plot predicted curves for each tree
plot_predictions(
  mod, by = c('time', 'tree', 'treat'),
  points = 0.5
) +
  theme(legend.position = 'none')

# Plot conditional curves for each treatment
plot_predictions(
  mod, condition = c('time', 'treat'),
)

# Plot marginal effects (first derivatives)
# for each treatment
plot_slopes(
  mod, 
  condition = c('time', 'treat'),
  variables = 'time'
)

# Are these slopes 'significantly' different?
# Use a sequence of times to compute slopes from each group and 
# test whether these slopes differ statistically
hypotheses(
  slopes(
    mod, 
    variables = 'time', 
    by = c('time', 'treat'),
    newdata = datagrid(
      time = seq(150, 250, length.out = 10),
      treat = unique
    )
  ),
  hypothesis = difference ~ sequential | time
)
