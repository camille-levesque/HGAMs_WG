# load packages

library(mvgam)           # Fit, interrogate and forecast DGAMs
library(dplyr)           # Tidy and flexible data manipulation
library(ggplot2)         # Flexible plotting
library(gratia)          # Graceful ggplot-based graphics for GAMs
library(marginaleffects) # Compute interpretable model predictions
library(tidybayes)       # Tidy manipulation / plots of posterior draws
library(tidyverse)



# load data
data <- read.csv("data/clean/data_195.csv")
    # Data was already cleaned by Katherine Hebert : data cleaning script in the data_cleaning folder
    # We import it from the "clean" folder (data/clean)

# adjust data

  # changer type colonne (as.factor, etc.)
  # formatter
  # creer nouvelles colonnes, etc.
