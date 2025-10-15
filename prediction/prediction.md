# Species Abundance Prediction Models using mvgam

## Overview

This document describes the Hierarchical Generalized Additive Models (HGAMs) developed for predicting songbird species abundance through time using North American Breeding Bird Survey (BBS) data. The models are implemented using the `mvgam` package and focus on temporal dynamics and species-level variation.

## Data

### Source
- **Dataset**: North American Breeding Bird Survey data (BioTIME study #195)
- **Temporal coverage**: 30 years (1978-2007)
- **Location**: Single site (latitude: 44.55, longitude: -74.4833)
- **Original records**: 699,142 across 383 species
- **Analysis subset**: Filtered to species with complete time series from 1978-2007

### Species Analyzed
The models focus on several songbird species, including:
- *Agelaius phoeniceus* (Red-winged Blackbird) - most frequently occurring species
- *Setophaga coronata* (Yellow-rumped Warbler)
- *Mniotilta varia* (Black-and-white Warbler) - used for new species prediction examples
- *Setophaga pinus* (Pine Warbler) - data-poor species example
- Additional warbler and songbird species

### Data Processing
- Data split into training (≤1999) and testing (>2000) periods
- Abundance converted to relative abundance within each year
- Time variable centered at zero for model fitting
- Species treated as factor levels for hierarchical modeling

## Models

### Primary Model (mod1)
**Model Type**: Hierarchical GAM with Zero-Mean Vector Normal (ZMVN) trend model

**Formula**:
```r
y ~ s(time, bs = "tp", k = 6) + s(series, bs = "re")
```

**Key Features**:
- **Response**: Poisson-distributed abundance counts
- **Temporal smooth**: Thin plate spline with 6 basis functions
- **Hierarchical structure**: Random intercepts for species (`s(series, bs = "re")`)
- **Latent variables**: `use_lv = TRUE` for dimension reduction
- **Trend model**: ZMVN for capturing temporal dynamics
- **Estimation**: Bayesian inference via Stan (4 chains, 500 burnin, 2000 samples)

### Alternative Model (mod)
**Model Type**: Hierarchical GAM with AR(1) trend model

**Formula**:
```r
y ~ s(time, bs = "tp", k = 6) + s(series, bs = "re")
```

**Key Features**:
- **Response**: Poisson distribution
- **Temporal dynamics**: AR(1) autocorrelation structure
- **Otherwise identical** to primary model
- **Purpose**: Comparison of trend model approaches (AR1 vs ZMVN)

### Advanced Hierarchical Model (mod_nick)
**Model Type**: State-Space Hierarchical GAM with AR(1) dynamics

**Formula**:
```r
# Observation model
y ~ s(series, bs = "re")

# Trend model
trend_formula = ~ s(time, bs = "tp", k = 6) + s(time, trend, bs = "sz", k = 6)
```

**Key Features**:
- **State-Space formulation**: Separates observation and process models
- **Hierarchical temporal smooths**: Both global and series-specifi
- c time effects
- **AR(1) residual correlation**: `trend_model = AR(p = 1)`
- **Non-centered parameterization**: Improved sampling efficiency
- **Priors**: Exponential(2) priors on sigma parameters

## New Species Prediction

### Approach 1: Direct Prediction
**Model**: `mod_noMniotilta` - trained excluding *Mniotilta varia*
- Uses global temporal smooth and species random effects
- Predicts new species based on population-level temporal patterns

### Approach 2: Post-Stratification
**Enhanced Model**: Includes similarity grouping variable
```r
y ~ s(time, bs = "tp", k = 5) + s(similar.species, bs="re") + s(series, bs = "re")
```

**Features**:
- **Similarity grouping**: Species categorized as "similar" or "not similar" to target species
- **Weighted predictions**: *Setophaga coronata* weighted highest for *Mniotilta varia* prediction
- **Hierarchical borrowing**: Information shared across ecologically similar species

### Approach 3: Phylogenetic/Functional Weighting (mod_nick)
**Post-stratification method** using species-specific weights:
- Weights based on ecological/phylogenetic distance to target species
- Species replicated in prediction grid based on similarity weights
- Marginalized predictions account for full uncertainty while ignoring AR(1) process

### Approach 4: State-Space Post-Stratification for BWAW (mod_strat_BAWW)
**Model Type**: State-Space Hierarchical GAM specifically for Black-and-white Warbler prediction

**Formula**:
```r
# Observation model
y ~ s(series, bs = "re")

# Trend model
trend_formula = ~ s(time, bs = "tp", k = 6) + s(time, trend, bs = "sz", k = 6)
```

**Training Data**: Excludes *Mniotilta varia* (Black-and-white Warbler) from `data_train_noBWAW`

**Post-stratification Strategy**:
- **Primary weight**: *Setophaga coronata* (Yellow-rumped Warbler) weighted **10x higher**
- **Secondary weights**: Other *Setophaga* species weighted **3x higher** 
- **Baseline weight**: All other species weighted **1x**
- **Method**: Species-specific replication in prediction grid based on ecological similarity

**Key Features**:
- Same State-Space structure as `mod_nick` (AR1 dynamics, hierarchical temporal smooths)
- Targeted post-stratification for ecologically similar warbler species
- Transparent weighting scheme based on taxonomic/ecological similarity

### Approach 5: Anchored Post-Stratification
**Method**: Calibrates post-stratified predictions using first-year observed data

**Implementation**:
```r
# Calculate offset from first observation
abundance_offset = first_year_BWAW$y - first_year_pred$estimate

# Apply to entire prediction timeline
post_strat_BWAW_anchored = post_strat_BWAW + abundance_offset
```

**Rationale**:
- **Temporal dynamics**: Leverages similarity-weighted trends from related species
- **Abundance calibration**: Anchors absolute abundance scale to observed data
- **Realistic predictions**: Combines model-based temporal patterns with species-specific abundance levels

**Advantages**:
- Reduces bias in absolute abundance predictions
- Maintains uncertainty structure from post-stratification
- Provides clear validation framework against remaining observations
- Intuitive interpretation: "What if BWAW follows similar temporal trends but at its observed abundance level?"

## Model Outputs

### Diagnostic Plots
- **Trace plots**: MCMC convergence assessment
- **R-hat histograms**: Parameter convergence diagnostics
- **Pairs plots**: Hierarchical parameter relationships

### Prediction Visualizations
1. **Future Predictions** (`A_TrendsFuturePredictions.jpeg`): Extrapolation beyond training period
2. **Training Trends** (`B_TrendsModel.jpeg`): Model fit to training data
3. **Combined Plots** (`C_TrendsBoth_nofacet.jpeg`, `D_TrendsBoth_facet.jpeg`): Training and prediction periods together
4. **New Species Predictions** (`E_TrendsNewSpecies.jpeg`): Out-of-sample species predictions
5. **Post-stratified BWAW** (`F_BWAW_PostStratified.jpeg`): State-space post-stratified predictions for Black-and-white Warbler
6. **Anchored vs Original** (`G_BWAW_Anchored_vs_Original.jpeg`): Comparison of anchored and original post-stratified predictions

## Model Comparison Notes

### Pending Analyses
- **Information criteria**: WAIC and LOO-CV comparisons between models
- **Trend model evaluation**: Systematic comparison of ZMVN vs AR(1) approaches
- **Basis function optimization**: Justification for k=6 and k=5 choices
- **Scale standardization**: Fixed y-axis scales across species in faceted plots

### Current Limitations
- Some species with very few observations (*Setophaga pinus* with 3 records) cause model fitting issues
- Post-stratification implementation could be enhanced with formal phylogenetic/functional distance metrics
- Mathematical model description needs updating to reflect current Poisson (not Beta) specification

## Technical Details

### Software Requirements
- **Primary package**: `mvgam` for multivariate GAM fitting
- **Backend**: `cmdstanr` for Stan model compilation and sampling
- **Visualization**: `marginaleffects` and `ggplot2` for prediction plots
- **Model diagnostics**: `bayesplot` integration for MCMC diagnostics

### Computational Specifications
- **Chains**: 2-4 MCMC chains depending on model
- **Burn-in**: 500 iterations
- **Sampling**: 2000 iterations post burn-in
- **Convergence**: Monitored via R-hat statistics and effective sample size

## Model Outputs and Files

### Saved Models
- `mvgam_prediction_mod1.rds`: Primary HGAM with ZMVN trend model
- `mvgam_prediction_mod.rds`: Alternative HGAM with AR(1) trend model
- `mod_strat_BAWW.rds`: State-Space HGAM for BWAW post-stratification

### Data Outputs
- `post_strat_BWAW_anchored.csv`: Anchored post-stratified predictions for BWAW
- `BWAW_prediction_comparison.csv`: Comparison of anchored vs original prediction accuracy

### Performance Metrics
The anchored post-stratification approach provides:
- **Mean Absolute Error** comparison between original and anchored predictions
- **Root Mean Square Error** assessment for both methods
- **Residual analysis** across the full time series
- **Validation framework** using first-year anchoring against remaining observations