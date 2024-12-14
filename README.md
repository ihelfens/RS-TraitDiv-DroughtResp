# Satellite-based trait-based diversity and drought response: Linear Modeling

This repository contains an R script for performing linear modeling on remotely-sensed trait-based functional diversity data. This script performs various linear and mixed-effects model analyses. We indicated the chosen model, tested the influence of the order of the contributions and exported the data for further analyses in MATLAB. Please find more information on the data acquisition and the analysis in the corresponding study.

 Helfenstein, I. S., Sturm, J. T., Schmid, B., Damm, A., Schuman, M. C., & Morsdorf, F. (2024). Satellite observations reveal positive relationship between trait-based diversity and drought response in temperate forests. EcoEvoRxiv [Preprint]. https://doi.org/10.32942/X24619

## Files

- `R_resp_Fdiv_def_rec_pub.R`: The R script for linear modeling.
- `data.csv`: The CSV file containing the data processed from Sentinel-2 data on trait-based functional diversity and drought response data fro Sentinel-2.

## Usage

1. Clone the repository.
2. Ensure you have R installed on your system.
3. Open R and set the working directory to the cloned repository.
4. Run the script using `source("R_resp_Fdiv_def_rec_pub.R")`.

The script will read the data from `data.csv`, perform the modeling, and save the results to `data_table_rst_eve.txt`.

## Contributions

Contains modified Copernicus Sentinel data [2017-2020].