# Load previous scripts ---------------------------------------------------
source("prediction/code/01_setup.R") # Source previous script, where we import & format data and load libraries
# Data exploration --------------------------------------------------------

##########################################

# Camille and Kim's data exploration (using Nick Clark's portal example)

# Look at the variables' format
dplyr::glimpse(dat)
  # We check the formats of the variables ; we do not need to change them.

# NAs
max(d_crop$YEAR)
image(is.na(t(dat %>%
                dplyr::arrange(dplyr::desc(time)))), axes = F,
      col = c('grey80', 'darkred'))
axis(3, at = seq(0,1, len = NCOL(dat)), 
     labels = colnames(dat))
  # No NAs

# Plot all of the time series together
plot_mvgam_series(data = dat, y = "y", series = "all")

# Plot some more in-depth features for individual series
## 4 first species
plot_mvgam_series(data = dat, y = "y", series = 1) 
plot_mvgam_series(data = dat, y = "y", series = 2)
plot_mvgam_series(data = dat, y = "y", series = 3)
plot_mvgam_series(data = dat, y = "y", series = 4)










