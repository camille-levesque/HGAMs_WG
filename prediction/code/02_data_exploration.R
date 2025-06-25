# Load previous scripts ---------------------------------------------------
# source("code/01-set_up.R")
# Data subsets ------------------------------------------------------------

##########################################
# Katherine Hebert's data exploration

# summarise data contents
dls = d_crop |>
  group_by(valid_name) |>
  group_split()

# check for basic population trends
mls = lapply(dls, function(x){lm(ABUNDANCE ~ YEAR, data = x)})

# extract slopes
coefs = mls |> lapply(coef) |> 
  bind_rows() |> 
  mutate("species" = unique(d_crop$valid_name))

# plot to see the slopes
ggplot(data = coefs) +
  geom_point(aes(y = species, x = YEAR))
#ggsave("variance/figures/linear_slopes.png")

##########################################

# Camille and Kim's data exploration (using Nick Clark's portal example)

# Look at the variables' format
dplyr::glimpse(d)
  # We check the formats of the variables ; we do not need to change them.

# NAs
max(d$YEAR)
image(is.na(t(d %>%
                dplyr::arrange(dplyr::desc(YEAR)))), axes = F,
      col = c('grey80', 'darkred'))
axis(3, at = seq(0,1, len = NCOL(d)), 
     labels = colnames(d))
  # No NAs

# Plot all of the time series together
plot_mvgam_series(data = d_crop, y = "rel_abun", series = "all")

# Plot some more in-depth features for individual series
plot_mvgam_series(data = d_crop, y = "rel_abun", series = 1)
plot_mvgam_series(data = d_crop, y = "rel_abun", series = 2)
plot_mvgam_series(data = d_crop, y = "rel_abun", series = 3)
plot_mvgam_series(data = d_crop, y = "rel_abun", series = 4)













