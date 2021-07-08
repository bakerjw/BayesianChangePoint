# Bayesian Change Point Calculation

Detect change points in a spatial Poisson processes: Application to induced seismicity in Oklahoma

The model is described in the following paper:

Gupta, A., and Baker, J. W. (2017). “Estimating Spatially Varying Event Rates with a Change Point using Bayesian Statistics: Application to Induced Seismicity.” Structural Safety, 65, 1–11. https://doi.org/10.1016/j.strusafe.2016.11.002

# File descriptions

changePointPredictionOklahoma.m - main file that runs change point analysis in Oklahoma

complete.csv - Earthquake catalog obtained from Oklahoma Geological Survey

CatalogFunctions - files to read and process earthquake catalog

ChangePoint - files to calculate Bayes Factor, and posterior distributions of time of change, and event rates

Declustering - files to implement Reasenberg declustering.

MappingData - files that contain information about state boundaries used to process latlong data

