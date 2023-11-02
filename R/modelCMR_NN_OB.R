
tar_option_set(packages = c("tidyverse", "nimble", "nimbleEcology", "MCMCvis", "targets"))

# library(getWBData)
# library(lubridate)
# library(kableExtra)
# library(GGally)
# library(nimble)
# library(nimbleEcology)
# library(MCMCvis)
# library(reticulate)
# library(tidyr)
# library(tidyverse)
# library(targets)

load('./models/cmrNN_OB/ss/runsOut/ss_NN_OB_mostRecent.RData')
out_ss_NN_OB <- d

# To get the mMCMCSummaryRMNA function which I adapted to deal with NAs
source('./models/cmrNN_OB/ss/modelCMR_ss_NN_OB_functionsToSource.R')

#Input data
#tar_load(eh_OB_2002_2014_target)
#eh <- eh_OB_2002_2014_target

modelCMR_NN_OB_target <-
  tar_plan(
    summary_ss_target = MCMCsummary(object = out_ss_NN_OB$mcmc, params = c("phiInt", "pInt"), round = 3) %>%
        rownames_to_column(var = "var") |> 
        mutate(mean_01 = exp(mean)/(1+exp(mean))),

    summary_ss_z0_target = MCMCSummaryRMNA(object = out_ss_NN_OB$mcmc, params = c("z"), round = 3) %>%
      rownames_to_column(var = "var") |>
      mutate(
        ind0 = str_match(var, "([0-9]+), ([0-9]+)")[,2],
        t0 =   str_match(var, "([0-9]+), ([0-9]+)")[,3],
        ind = as.numeric(ind0),
        t = as.numeric(t0)
      ) |>
      dplyr::select(-c(t0, ind0)),
    
    # 'phi' has each individual, 'phiInt' does not
    summary_ss_phi0_target = MCMCSummaryRMNA(object = out_ss_NN_OB$mcmc, params = c("phi"), round = 3) %>%
      rownames_to_column(var = "var") |>
      mutate(
        ind0 = str_match(var, "([0-9]+), ([0-9]+)")[,2],
        t0 =   str_match(var, "([0-9]+), ([0-9]+)")[,3],
        ind = as.numeric(ind0),
        t = as.numeric(t0)
      ) |>
      dplyr::select(-c(t0, ind0)),
    
    summary_ss_p0_target = MCMCSummaryRMNA(object = out_ss_NN_OB$mcmc, params = c("p"), round = 3) %>%
      rownames_to_column(var = "var") |>
      mutate(
        ind0 = str_match(var, "([0-9]+), ([0-9]+)")[,2],
        t0 =   str_match(var, "([0-9]+), ([0-9]+)")[,3],
        ind = as.numeric(ind0),
        t = as.numeric(t0)
      ) |>
      dplyr::select(-c(t0, ind0)),
    
    ehLong_target =
      eh_OB_2002_2014_target$eh |>
      as.data.frame() |>
      rownames_to_column("ind0") |>
      pivot_longer(starts_with("ais")) |>
      mutate(
        t = as.numeric(str_match(name, "([0-9]+)")[,1]),
        enc = value,
        ind = as.numeric(ind0)
      ) |>
      dplyr::select(-c(name, value, ind0)),
    
    firstLong_target = eh_OB_2002_2014_target$first |>
      as_tibble() |>
      rownames_to_column("ind") |>
      mutate(ind = as.numeric(ind)) |>
      rename(first = value),
    
    lastLong_target = eh_OB_2002_2014_target$last |>
      as_tibble() |>
      rownames_to_column("ind") |>
      mutate(ind = as.numeric(ind)) |>
      rename(last = value),
    
    cohortLong_target = eh_OB_2002_2014_target$cohort |>
      as_tibble() |>
      rownames_to_column("ind") |>
      mutate(ind = as.numeric(ind)),
    
    #### z
    summary_ss_z_target = summary_ss_z0_target |>
      left_join(ehLong_target) |>
      mutate(
        meanM1 = mean - 1,
        pSurv = ifelse(meanM1 == 0, 1, 1 - meanM1),
        enc8 = ifelse(enc == 1, 0.8, 0)
      ) |>
      left_join(firstLong_target) |>
      left_join(lastLong_target) |>
      left_join(cohortLong_target),
    
    ##### p
    summary_ss_p_target = summary_ss_p0_target |>
      left_join(ehLong_target) |>
      mutate(
        meanM1 = mean - 1,
        pSurv = ifelse(meanM1 == 0, 1, 1 - meanM1),
        enc8 = ifelse(enc == 1, 0.2, 0)
      ) |>
      left_join(firstLong_target) |>
      left_join(lastLong_target) |>
      left_join(cohortLong_target),
    
    #### phi
    summary_ss_phi_target = summary_ss_phi0_target |>
      left_join(ehLong_target) |>
      mutate(
        meanM1 = mean - 1,
        pSurv = ifelse(meanM1 == 0, 1, 1 - meanM1),
        enc8 = ifelse(enc == 1, 0.2, 0)
      ) |>
      left_join(firstLong_target) |>
      left_join(lastLong_target) |>
      left_join(cohortLong_target)
  )

