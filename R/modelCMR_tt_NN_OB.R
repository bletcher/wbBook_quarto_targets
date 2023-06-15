# not using this for targets. Just loading model output in .qmd file


tar_option_set(packages = c("tidyverse", "nimble", "nimbleEcology", "MCMCvis", "targets"))

# load('./models/cmrFlowWB/runsOut/fromWorkbench/rstudio-export/mcmc_phiT_pT_psiT_DHMM_dirch_mostRecent.RData',
#      tmpEnv <- new.env())
# 
# modelCMR_ttt_WB_target <-
#   tar_plan(
#     toSave_ttt_WB_target = tmpEnv$toSave
#   )
    

#library(targets)
# https://books.ropensci.org/targets/data.html#external-files
create_output <- function(file) {
  #data <- read.csv(file)
  load(file) # loads 'd'
  output <-  d#head(d$mcmc$chain1[,1:10]) #head(data)
  save(output, file = "output.RData")
  "output.RData"
}

modelCMR_tt_NN_OB_target <- 
  list(
    tar_target(name = mcmcIn_NN_OB, command = "./models/cmrNN_OB/tt/runsOut/tt_NN_OB_mostRecent.RData", format = "file"), 
    tar_target(name = mcmcOut_NN_OB, command = create_output(mcmcIn_NN_OB), format = "file")
  )