# load packages

library(mvgam)           # Fit, interrogate and forecast DGAMs
library(dplyr)           # Tidy and flexible data manipulation
library(ggplot2)         # Flexible plotting
library(gratia)          # Graceful ggplot-based graphics for GAMs
library(marginaleffects) # Compute interpretable model predictions
library(tidyverse)
library(mgcv)


## devtools: you need to install the package via the GitHub repository

## install.packages("bbsBayes2", repos = c(bbsbayes = "https://bbsbayes.r-universe.dev",
##                                         CRAN = getOption("repos")))
library(bbsBayes2) # v1.1.2.1
library(cmdstanr) # v0.7.1

cmdstanr::install_cmdstan()

cmdstanr::check_cmdstan_toolchain(fix = TRUE)

fetch_bbs_data()



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

# Rename columns and adjust abundance 
d_crop <- d_crop %>%
  rename(
    abun = ABUNDANCE,
    series = valid_name, # for the mvgam requirement
    lat = LATITUDE,
    long = LONGITUDE,
    time = YEAR # for the mvgam requirement
    ) 

d_crop <- d_crop %>%
  group_by(time) %>%
  mutate(total_abun = sum(abun)) %>%
  mutate(rel_abun = abun/total_abun) %>%
  select(-total_abun)

# Conditions

## Adjust columns for mvgam requirements
d_crop %>%
  dplyr::mutate(series = as.factor(series)) -> d_crop
dplyr::glimpse(d_crop) # we now only have integers, factors or numeric columns
levels(d_crop$series) # 28 species

## Verify if there is 0 or 1 values for the relative abundance 
min(d_crop$rel_abun)
max(d_crop$rel_abun)




  