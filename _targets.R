# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(getWBData)
library(future)
library(future.callr)

# for parallel processing, https://books.ropensci.org/targets/hpc.html#future-locally
plan(callr)

# Set target options:
tar_option_set(
  packages = c(
     "tibble",     
     "tidyr",
     "lubridate",
     "knitr",
     "targets",
     "tarchetypes",
     "getWBData",
     "getPrepareWBData",
     "daymetr",
     "nimble",
     "nimbleEcology",
     "tidyverse"
  ), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
  #error= "null"
)

# tar_make_clustermq() configuration (okay to leave alone):
#options(clustermq.scheduler = "multiprocess")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

#Connect to the WB database on AWS
#reconnect()

# Run the R scripts in the /R folder with your custom functions:
# recursive = FALSE does not go into subfolders (e.g. `hold`)
lapply(list.files("R", pattern = "*.R", full.names = TRUE, recursive = FALSE), source)

list(
  getEnvData_target,
  modelFDC_target,  
  getElectroData_target,
  dataCMR_WB_2002_2014_target,
  dataCMR_WBbkt_2002_2014_target,
  dataCMR_WBbnt_2002_2014_target,
  dataCMR_OB_2002_2014_target,
  dataWanding_target,
  dataAntenna_target,
  dataAll_target,
  modelYOY_target,
  modelFlow_target,
  modelConditionFactor_target,
  modelGrowthInMass_target,

  modelXGBoost_target,
  modelCMR_NN_OB_target
  

  # 
  # # retrieves the mcmc files
  # # run the files 'by hand' in ./models/cmrFlowWB/modelCMR_ttt_ft_cohort_WB_makeFile.R
  # modelCMR_ttt_WB_target
  
  # turn off for now
  # modelCMR_tt_ft_OB_flow_target,
  # modelCMR_tt_ft_OB_flowByRiver_target,
  # 
  # modelCMR_tt_ft_cohort_OB_flow_target,
  # modelCMR_tt_ft_cohort_OB_flowByRiver_target
  
  # old models
  # modelCMR_tt_OB_target
  
  # modelCMR_tt_OB_flow_target,
  # modelCMR_tt_OB_flowByRiver_target,
  #run_ttt_models_target
  #modelCMR_ttt_ft_cohort_WB_flow_target,
  #modelCMR_ttt_ft_cohort_WB_flowByRiver_target
  
  
  
  #tar_quarto(book) # not exactly sure what this does, except create correct tar_visualize() result
)