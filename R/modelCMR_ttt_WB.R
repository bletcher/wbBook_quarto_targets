tar_option_set(packages = c("tidyverse", "nimble", "nimbleEcology", "MCMCvis", "targets"))

load('./models/cmrFlow4rivers/runsOut/fromWorkbench/rstudio-export/mcmc_phiT_pT_psiT_DHMM_dirch_mostRecent.RData',
     tmpEnv <- new.env())

modelCMR_ttt_WB_target <-
  tar_plan(
    toSave_ttt_WB_target = tmpEnv$toSave
  )
    
    
######################################
#### Functions
#####################################

