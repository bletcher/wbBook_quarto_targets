# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(getWBData)

# Set target options:
tar_option_set(
  packages = c("tibble",     
               "tidyverse",
               "lubridate",
               "knitr",
               "targets",
               "tarchetypes",
               "getWBData",
               "getPrepareWBData"
               ), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multiprocess")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

#Connect to the WB database on AWS
#reconnect()

# Run the R scripts in the /R folder with your custom functions:
lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)

list(
  #target_globalVariables
  getEnvData_target,
  getElectroData_target,
  dataCMR_target,
  dataWanding_target,
  dataAntenna_target,
  modelYOY_target
  #tar_quarto(report, need to fill in
)
