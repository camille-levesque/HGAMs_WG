# load packages

library(mvgam) # Fit, interrogate and forecast DGAMs
library(dplyr) # Tidy and flexible data manipulation
library(ggplot2) # Flexible plotting
library(gratia) # Graceful ggplot-based graphics for GAMs
library(marginaleffects) # Compute interpretable model predictions
library(tidyverse)
library(mgcv)
library(cmdstanr)
cmdstan_path <- Sys.getenv("CMDSTAN")
cmdstanr::set_cmdstan_path(cmdstan_path)

# load data

## clean data
d <- read.csv("data/clean/data_195.csv") |>
  select(-c(SAMPLE_DESC, DAY, MONTH, BIOMAS))
# Data was already cleaned by Katherine Hebert (KH) : data cleaning script in the data_cleaning folder
# We import it from the "clean" folder (data/clean)
# Data in the long format

## clean cropped data
d_crop <- readRDS("variance/example_species_df.rds")
# Data was cropped by KH in her script (variance/01_fit-model.R)
# Data regrouped by species per year (for one localization 44.55;-74.4833)

# Add a rare/undersampled species for Example 1
d_crop_rare <- readRDS("prediction/example_rare_species.rds")

d_crop_merged <- rbind(d_crop, d_crop_rare)
# Rename columns and adjust abundance
dat <- d_crop_merged %>%
  rename(
    y = ABUNDANCE,
    series = valid_name, # for the mvgam requirement
    lat = LATITUDE,
    long = LONGITUDE,
    time = YEAR # for the mvgam requirement
  )

dat <- dat %>%
  group_by(time) %>%
  mutate(total_abun = sum(y)) %>%
  mutate(rel_abun = y / total_abun) %>%
  select(-total_abun)

# mvgam format requirements

unique(dat$series)

## Adjust columns for mvgam requirements
dat$y <- as.vector(dat$y)
dat <- filter(dat, series != "Vireo olivaceus")
dat$series <- as.factor(dat$series)
dat$time <- as.integer(dat$time) - min(dat$time)

# subsetting the data in training and testing folds
data_train <- dat[which(d_crop$YEAR <= 1999), ]
data_test <- dat[which(d_crop$YEAR > 2000), ]

# subsetting the data with and without a species for out-of-sample forecasting
# Replacing red-winged black-bird with a species that has a similar ecology to the rest
data_noMniotilta <- filter(dat, series != "Mniotilta varia")
data_Mniotilta <- filter(dat, series == "Mniotilta varia")
