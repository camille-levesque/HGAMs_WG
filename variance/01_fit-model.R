# Script to model population trends from the BBS dataset and to output the
# covariance matrix of species associations

# load libraries

library(mvgam)
library(mgcv)
library(dplyr)
library(ggplot2)

theme_set(theme_minimal())

# load data

d = read.csv("data/clean/data_195.csv") |>
  select(-c(SAMPLE_DESC, DAY, MONTH, BIOMAS))

sites = paste0(d$LONGITUDE,"_", d$LATITUDE)
table(sites)

# crop to a site
d_crop = d[sites %in% "-74.4833_44.55",]

## crop to birds 
species = table(d_crop$valid_name)
species = species[which(species == 30)]
ex = names(species)
#saveRDS(ex, "variance/example_species_list.rds")
ex <- readRDS("variance/example_species_list.rds")

d_crop = dplyr::filter(d_crop, valid_name %in% ex)
d_crop = dplyr::filter(d_crop, valid_name != "Vireo olivaceus")
saveRDS(d_crop, "variance/example_species_df.rds")


## data exploration ------------------------------------------------------------

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
ggsave("variance/figures/linear_slopes.png")


## build a model for each species ----------------------------------------------

# check for basic population trends
gamls = lapply(dls, function(x){gam(ABUNDANCE ~ s(YEAR, k = 4), data = x)})
lapply(gamls, plot)

## build an mvgam model for all species together ------------------------------

# prepare the data ----

# format into long
dat = d_crop |>
  select(valid_name, ABUNDANCE, YEAR) |>
  rename("time" = "YEAR",
         "series" = "valid_name",
         "y" = "ABUNDANCE")
dat$y = as.vector(dat$y)
dat <- filter(dat, series != "Vireo olivaceus")
dat$series <- as.factor(dat$series)
dat$time <- as.integer(dat$time)-min(dat$time)


data_train = dat[which(d_crop$YEAR <= 1999),]
data_test = dat[which(d_crop$YEAR > 2000),]

ggplot(data = data_train) +
  geom_smooth(aes(x = time, y = y, 
                  col = series, fill = series),
              alpha = .1)
ggsave("variance/figures/individual_gams.png")

npops = length(unique(dat$series))

# prepare the priors ----
knots = 4

# mvgam_prior <- mvgam(data = data_train,
#                      formula = y ~ 
#                        # global smoother for all pops over time
#                        s(time, bs = "tp", k = knots) + 
#                        # random intercept per group
#                        s(series, bs = 're', k = npops),
#                      family = "poisson",
#                      trend_model = 'GP',
#                      chains = 3,
#                      use_stan = TRUE,
#                      prior_simulation = TRUE)
# 
# # record the priors
# test_priors <- get_mvgam_priors(y ~ 
#                                   # global smoother for all pops over time
#                                   s(time, bs = "tp", k = knots) + 
#                                   # random intercept per group
#                                   s(series, bs = 're', k = npops),
#                                 family = "poisson",
#                                 data = data_train,
#                                 trend_model = 'GP',
#                                 use_stan = TRUE)
# write.csv(test_priors, "variance/outputs/test_priors.csv")

# # look at the priors
# plot(mvgam_prior, type = 'smooths')
# 
# png("figures/trend_priors.png", width = 1000, height = 1300, type = "cairo")
# par(mfrow = c(6,5))
# for(i in 1:ncol(Year_Geom_Means_all)){
#   plot(mvgam_prior, type = 'trend', series = i)
# } 
# dev.off()
# 
# png("figures/re_priors.png", width = 1000, height = 700, type = "cairo")
# plot(mvgam_prior, type = 're')
# dev.off()

# check number of knots ----
knots = 4

hgam = mgcv::gam(y ~ s(time, bs = "tp", k = 6) + 
                         s(series, bs = 're', k = npops), 
                       family = "poisson",
                       data = data_train)
mgcv::gam.check(hgam)
plot(hgam)

# train the model on data ----
m1 <- mvgam(data = data_train,
              formula =  y ~ s(time, bs = "tp", k = 6) + 
                s(series, bs = "re"),
              use_lv = TRUE,
              family = "poisson",
              trend_model = 'AR1',
              use_stan = TRUE,
              chains = 2, 
              burnin = 500,
              samples = 2000
)
saveRDS(m1, paste0("variance/outputs/mvgam_variance.rds")) 
# model didn't converge very well

plot(m1)


# species associations ----

sp_corr = mvgam::lv_correlations(m1)

# clean up the species names 

colnames(sp_corr$mean_correlations) = gsub("_", " ", colnames(sp_corr$mean_correlations)) |> 
  stringr::str_to_sentence()
rownames(sp_corr$mean_correlations) = gsub("_", " ", rownames(sp_corr$mean_correlations)) |> 
  stringr::str_to_sentence()


# Plot as heatmap --------------------------------------------------------------

png(height=1800, width=1800, file="variance/figures/species_associations.png", type = "cairo")
corrplot::corrplot(sp_corr$mean_correlations, 
                   type = "lower",
                   method = "color", 
                   tl.cex = 2.5, cl.cex = 3, tl.col = "black", font = 3)
dev.off()
