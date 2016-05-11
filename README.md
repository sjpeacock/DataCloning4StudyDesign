# Study design and parameter estimability for spatial and temporal ecological models

Submitted May 11, 2016 to Methods in Ecology and Evolution

Stephanie J Peacock (stephanie.j.peacock at gmail.com), Martin Krkosek, Mark Lewis, and Subhash Lele

## Abstract
- The statistical tools available to ecologists are becoming increasingly sophisticated, allowing more complex, mechanistic models to be fit to ecological data.
- Such models have the potential to provide new insights into the processes underlying ecological patterns, but the inferences made are limited by the information in the data.
- Statistical non-estimability of model parameters due to insufficient information in the data is a problem too-often ignored by ecologists employing complex models.
- Here, we show how a new statistical computing method called data cloning can be used to inform study design by assessing the estimability of parameters under different spatial and temporal scales of sampling. 
- A case study of parasite transmission from farmed to wild salmon highlights that assessing the estimability of ecologically relevant parameters should be a key step when designing studies in which fitting complex mechanistic models is the end goal.

*Keywords:* modelling, spatial or time-series, statistics


## Contents of repository

This repository contains three folders with the following files:

1. Code
  - `figures.R`: loads `.RData` workspace containing results from the Workspaces subfolder and creates the figures for the main text of the paper
  - `model.R`: the model function for JAGS, called by the `SeaLice_dclone_xxx.R` files
  - `sim_model.R`: version of `model.R` in `R` syntax that takes parameters and gives expected number of lice per fish
  - `SeaLice_dclone_lessSpread.R`: simulates data for the less-spread data scenario and fits the model using data cloning
  - `SeaLice_dclone_moreSpread.R`: simulates data for the more-spread data scenario and fits the model using data cloning
  - `SeaLice_dclone_original.R`: fits the model to the original dataset using data cloning

2. Data
  - `Leps.txt` contains the sea louse data
  - `Summary.txt` contains the site data, including distance along migration route
  
3. Supplement: contains code, workspaces and figures for the supplemental material looking at 4 different prior assumptions.
