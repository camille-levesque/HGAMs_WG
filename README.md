# Working Group: Hierarchical GAMs in ecology 📝

## Data

For the coded examples, we will be using the North American Breeding Bird Survey data. The raw data is stored in `data/raw/raw_data_195.csv` along with a metadata file and a text file containing information about how to cite the data. We then cleaned the data a little bit as described below. If you are running analyses, __please use the cleaned dataset in `data/clean/data_195.csv`__.

### Data description

The raw data summary from BioTIME describes the dataset's characteristics:

Realm: Terrestrial
Climate: Temperate
Biome: Temperate grasslands, savannas and shrublands
Central latitude: 40.809241
Central longitude: -96.187269
Duration: 30 years, from 1978 to 2007
699142 records
383 distinct species
Across the time series, _Agelaius phoeniceus_ is the most frequently occurring species.

### Data source

We downloaded the data on June 17, 2025 from BioTIME ([link](https://biotime.st-andrews.ac.uk/selectStudy.php?study=195)), and the data should be cited as noted in `data/raw/citation.txt`:

> Pardieck, K. L., Ziolkowski Jr, D. J., Hudson, M. A. R., & Campbell, K. (2015). North American breeding bird survey dataset 1966–2014, version 2014.0. US Geological Survey, Patuxent Wildlife Research Center.

The dataset is from BioTIME, which should be cited as:
> Dornelas, M., Antao, L. H., Moyes, F., Bates, A. E., Magurran, A. E., Adam, D., ... & Murphy, G. (2018). BioTIME: A database of biodiversity time series for the Anthropocene. Global Ecology and Biogeography, 27(7), 760-786.

### Data cleaning

The script `data_cleaning/clean-bbs-data.R` shows how the raw data was cleaned to make `data/clean/data_195.csv`. 

For now, cleaning involves a single step:
1. We filtered the dataset to keep species' whose time series begain in 1978 and ended in 2007.

However, if you clean the dataset more as you analyse it (it is still fresh!), please update the script `data_cleaning/clean-bbs-data.R` and the clean dataset `data/clean/data_195.csv`. Be careful of course, but don't worry too much about overwriting something important - we are using GitHub, which allows us to flip back to older versions if we need to.
