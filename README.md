# DataCloning4StudyDesign
Data and code to accompy the manuscript "Data cloning can guide study design to ensure parameter estimability in complex ecological models".  

This repository contains three folders:
1) Code
  - figures.R: loads .RData workspace containing results from the Workspaces subfolder and creates the figures for the main text of the paper.
  - model.R: the model function for JAGS, called by the SeaLice_dclone_xxx.R files
  - SeaLice_dclone_lessSpread.R: simualtes data for the less-spread data scenario and fits the model using data cloning.
  - SeaLice_dclone_moreSpread.R: simualtes data for the more-spread data scenario and fits the model using data cloning.
  - SeaLice_dclone_original.R: fits the model to the original dataset using data cloning.

2) Data
  - Leps.txt contains the sea louse data
  - Summary.txt contains the site data, including distance along migration route
  
