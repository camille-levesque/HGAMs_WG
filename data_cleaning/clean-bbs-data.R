# Script to clean the BBS data downloaded from BioTIME

# load libraries
library(ggplot2)
library(dplyr)

# load
df = read.csv("data/raw/raw_data_195.csv")

# filter to years in common across many species
df_years = df |>
  group_by(valid_name) |>
  summarise("min_year" = min(YEAR),
            "max_year" = max(YEAR))
# look for the most common minimum and maximum years
df_years$min_year |> table()
df_years$max_year |> table()

# 309 species to keep
sp_to_keep = df_years |> filter(min_year == 1978, max_year == 2007)

# Cut to species that have data between 1978 and 2007
df = df |> filter(valid_name %in% sp_to_keep$valid_name)

# save this dataset
write.csv(df, "data/clean/data_195.csv", row.names = FALSE)