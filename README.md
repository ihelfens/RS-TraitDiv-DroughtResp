# Satellite-based trait-based diversity and drought response: Linear Modeling

This repository contains the script for performing linear modeling on remotely-sensed trait-based functional diversity data. This script performs linear and mixed-effects model analyses. We indicated the chosen model, tested the influence of the order of the contributions and exported the data for further analyses and creating figures.

## Files

- `R_resp_Fdiv_pub.R`: The R script for linear modeling.
- `R_resp_Fdiv_reg.csv`: The CSV file containing the data processed from Sentinel-2 data on trait-based functional diversity and drought response data from Sentinel-2.

## Usage

1. Clone the repository.
2. Ensure you have R installed on your system.
3. Open R and set the working directory to the cloned repository.
4. Run the script using `source("R_resp_Fdiv_pub.R")`.

The script will read the data from `R_resp_Fdiv_reg.csv`, perform the modeling.

## Contributions
Contains modified Copernicus Sentinel data [2017-2020].
