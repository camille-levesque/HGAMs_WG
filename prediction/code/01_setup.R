# load packages

library(mvgam)           # Fit, interrogate and forecast DGAMs
library(dplyr)           # Tidy and flexible data manipulation
library(ggplot2)         # Flexible plotting
library(gratia)          # Graceful ggplot-based graphics for GAMs
library(marginaleffects) # Compute interpretable model predictions
library(tidyverse)
library(mgcv)



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

# Rename columns
d_crop <- d_crop %>%
  rename(
    abun = ABUNDANCE,
    series = valid_name,
    lat = LATITUDE,
    long = LONGITUDE,
    year = YEAR
    ) 

d_crop <- d_crop %>%
  group_by(year) %>%
  mutate(total_abun = sum(abun)) %>%
  mutate(rel_abun = abun/total_abun) %>%
  select(-total_abun)
  