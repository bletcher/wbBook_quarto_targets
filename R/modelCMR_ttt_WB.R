tar_option_set(packages = c("tidyverse", "nimble", "nimbleEcology", "MCMCvis", "targets"))

# load('./models/cmrFlowWB/runsOut/fromWorkbench/rstudio-export/mcmc_phiT_pT_psiT_DHMM_dirch_mostRecent.RData',
#      tmpEnv <- new.env())
# 
# modelCMR_ttt_WB_target <-
#   tar_plan(
#     toSave_ttt_WB_target = tmpEnv$toSave
#   )
    

#library(targets)
create_output <- function(file) {
  data <- read.csv(file)
  output <- head(data)
  write.csv(output, "output.csv")
  "output.csv"
}

modelCMR_ttt_WB_target <- 
  list(
    tar_target(name = mcmcInMod1, command = "./models/cmrFlowWB/runsOut/mod1/ttt_WB_mcmc_mod1.csv", format = "file"), 
    tar_target(name = mcmcOutMod1, command = create_output(mcmcInMod1), format = "file"),
    tar_target(name = mcmcInMod2, command = "./models/cmrFlowWB/runsOut/mod2/ttt_WB_mcmc_mod2.csv", format = "file"),
    tar_target(name = mcmcOutMod2, command = create_output(mcmcInMod2), format = "file")
  )