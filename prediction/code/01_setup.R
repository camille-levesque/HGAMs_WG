# load packages

library(mvgam)           # Fit, interrogate and forecast DGAMs
library(dplyr)           # Tidy and flexible data manipulation
library(ggplot2)         # Flexible plotting
library(gratia)          # Graceful ggplot-based graphics for GAMs
library(marginaleffects) # Compute interpretable model predictions
library(tidybayes)       # Tidy manipulation / plots of posterior draws
library(tidyverse)


## devtools: you need to install the package via the GitHub repository

##install.packages("bbsBayes2", repos = c(bbsbayes = "https://bbsbayes.r-universe.dev",
##                                        CRAN = getOption("repos")))
library(bbsBayes2) # v1.1.2.1
library(cmdstanr) # v0.7.1

cmdstanr::install_cmdstan()

cmdstanr::check_cmdstan_toolchain(fix = TRUE)

fetch_bbs_data()




# load data

# adjust data

  # changer type colonne (as.factor, etc.)
  # formatter
  # creer nouvelles colonnes, etc.

#vignettes/articles/bbsBayes2.Rmd