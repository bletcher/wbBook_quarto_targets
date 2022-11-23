tar_option_set(packages = c("tidyverse", "nimble", "nimbleEcology", "MCMCvis", "targets"))

#load('./models/cmrFlowOB/dataOut/eh_2002200320042005200620072008200920102011201220132014_wb obear.RData')

#eh_OB_target = tar_read(eh_OB_2002_2014_target)

modelCMR_tt_OB_target <-
  tar_plan(
  
    tt_inputData_OB =
      list(
        y = eh_OB_2002_2014_target$eh,
        first = eh_OB_2002_2014_target$first, 
        last = eh_OB_2002_2014_target$last, 
        zInits = tt_initialValues_OB(ncol(eh_OB_2002_2014_target$eh), 
                                          eh_OB_2002_2014_target$eh) 
       # nStates = length(unique(eh_OB_2002_2014_target$data$sizeState)),
       # nRivers = length(unique(eh_OB_2002_2014_target$data$sizeState)) # for now
      ), 
    
    tt_runData_OB = list(
      # Updateable model-specific variables 
      nIter = 5000, 
      nBurnin = 2000, 
      nChains = 2,
      thinRate = 5
    ),  
    
    tt_myConstants_OB = list(
      N = nrow(tt_inputData_OB$y),
      T = ncol(tt_inputData_OB$y),
      first = tt_inputData_OB$first,
      last = tt_inputData_OB$last
      
      #nRivers = tt_inputData_OB$nRivers,
      #length = tt_inputData_OB$last - tt_inputData_OB$first + 1,
      
      #alphaR1 = tt_alpha_OB$alphaR1,
      #alphaR2 = tt_alpha_OB$alphaR2,
      #alphaR3 = tt_alpha_OB$alphaR3,
      
      #deltaProps = tt_inputData_OB$deltaProps,
      #nStates = tt_inputData_OB$nStates
    ),
    
    tt_myData_OB = list(
      y = tt_inputData_OB$y + 1
    ),
      
    #model-specific variables
    # tt model
    tt_modelCode_OB = nimbleCode({
      delta[1] <- 1                    # Pr(alive t = 1) = 1
      delta[2] <- 0                    # Pr(dead t = 1) = 0
      for (t in 1:(T-1)){ # loop over time
        phi[t] ~ dunif(0, 1)           # prior survival
        gamma[1,1,t] <- phi[t]         # Pr(alive t -> alive t+1)
        gamma[1,2,t] <- 1 - phi[t]     # Pr(alive t -> dead t+1)
        gamma[2,1,t] <- 0              # Pr(dead t -> alive t+1)
        gamma[2,2,t] <- 1              # Pr(dead t -> dead t+1)
        p[t] ~ dunif(0, 1)             # prior detection
        omega[1,1,t] <- 1 - p[t]       # Pr(alive t -> non-detected t)
        omega[1,2,t] <- p[t]           # Pr(alive t -> detected t)
        omega[2,1,t] <- 1              # Pr(dead t -> non-detected t)
        omega[2,2,t] <- 0              # Pr(dead t -> detected t)
      }
      # likelihood
      for (i in 1:N){
        z[i,first[i]] ~ dcat(delta[1:2])
        for (j in (first[i]+1):last[i]){
          z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
          y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1])
        }
      }
    }),


    tt_Rmodel_OB = nimbleModel(
      code = tt_modelCode_OB,
      constants = tt_myConstants_OB,
      data = tt_myData_OB,
      inits = tt_initialValues_OB(tt_myConstants_OB$T, tt_inputData_OB$y),
      calculate = FALSE
    ),

    tt_parametersToSave_OB = c("phi", "p"),

    tt_conf_OB = configureMCMC(
      tt_Rmodel_OB,
      monitors = tt_parametersToSave_OB
    ),

    tt_Rmcmc_OB = buildMCMC(tt_conf_OB, useConjugacy = FALSE),
    tt_Cmodel_OB = compileNimble(tt_Rmodel_OB),
    tt_Cmcmc_OB = compileNimble(tt_Rmcmc_OB, project = tt_Rmodel_OB),

    tt_model_OB = runMCMC(
      tt_Cmcmc_OB,
      niter = tt_runData_OB$nIter,
      nburnin = tt_runData_OB$nBurnin,
      thin = tt_runData_OB$thinRate,
      nchains = tt_runData_OB$nChains
    ),

    tt_modelOut_OB =
      list(
        mcmc = tt_model_OB, # "Error : invalid nimbleFunction argument"
        name = "phiT_pT_OB",
        modelCode = tt_modelCode_OB,
        myConstants = tt_myConstants_OB,
        runData = tt_runData_OB
      ),

    tt_save_OB = saveModelOut_tt_OB(tt_modelOut_OB)
  )

###################################################
# Add flow
###################################################

modelCMR_tt_OB_flow_target <-
  tar_plan(
    
    tt_inputData_OB_flow =
      list(
        y = eh_OB_2002_2014_target$eh,
        first = eh_OB_2002_2014_target$first, 
        last = eh_OB_2002_2014_target$last, 
        zInits = tt_initialValues_OB(ncol(eh_OB_2002_2014_target$eh), 
                                     eh_OB_2002_2014_target$eh) 
        # nStates = length(unique(eh_OB_2002_2014_target$data$sizeState)),
        # nRivers = length(unique(eh_OB_2002_2014_target$data$sizeState)) # for now
      ), 
    
    tt_runData_OB_flow = list(
      # Updateable model-specific variables 
      nIter = 5000, 
      nBurnin = 2000, 
      nChains = 2,
      thinRate = 5
    ),  
    
    tt_myConstants_OB_flow = list(
      N = nrow(tt_inputData_OB_flow$y),
      T = ncol(tt_inputData_OB_flow$y),
      first = tt_inputData_OB_flow$first,
      last = tt_inputData_OB_flow$last
      
      #nRivers = tt_inputData_OB_flow$nRivers,
      #length = tt_inputData_OB_flow$last - tt_inputData_OB_flow$first + 1,
      
      #alphaR1 = tt_alpha_OB_flow$alphaR1,
      #alphaR2 = tt_alpha_OB_flow$alphaR2,
      #alphaR3 = tt_alpha_OB_flow$alphaR3,
      
      #deltaProps = tt_inputData_OB_flow$deltaProps,
      #nStates = tt_inputData_OB_flow$nStates
    ),
    
    tt_myData_OB_flow = list(
      y = tt_inputData_OB_flow$y + 1
    ),
    
    #model-specific variables
    # tt model
    ## need to add flow effects using flow
    ## then need to add flow effects using flowByRiver in another set of targets
    tt_modelCode_OB_flow = nimbleCode({
      delta[1] <- 1                    # Pr(alive t = 1) = 1
      delta[2] <- 0                    # Pr(dead t = 1) = 0
      for (t in 1:(T-1)){ # loop over time
        phi[t] ~ dunif(0, 1)           # prior survival
        gamma[1,1,t] <- phi[t]         # Pr(alive t -> alive t+1)
        gamma[1,2,t] <- 1 - phi[t]     # Pr(alive t -> dead t+1)
        gamma[2,1,t] <- 0              # Pr(dead t -> alive t+1)
        gamma[2,2,t] <- 1              # Pr(dead t -> dead t+1)
        p[t] ~ dunif(0, 1)             # prior detection
        omega[1,1,t] <- 1 - p[t]       # Pr(alive t -> non-detected t)
        omega[1,2,t] <- p[t]           # Pr(alive t -> detected t)
        omega[2,1,t] <- 1              # Pr(dead t -> non-detected t)
        omega[2,2,t] <- 0              # Pr(dead t -> detected t)
      }
      # likelihood
      for (i in 1:N){
        z[i,first[i]] ~ dcat(delta[1:2])
        for (j in (first[i]+1):last[i]){
          z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
          y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1])
        }
      }
    }),
    
    
    tt_Rmodel_OB_flow = nimbleModel(
      code = tt_modelCode_OB_flow,
      constants = tt_myConstants_OB_flow,
      data = tt_myData_OB_flow,
      inits = tt_initialValues_OB_flow(tt_myConstants_OB_flow$T, tt_inputData_OB_flow$y),
      calculate = FALSE
    ),
    
    tt_parametersToSave_OB_flow = c("phi", "p"),
    
    tt_conf_OB_flow = configureMCMC(
      tt_Rmodel_OB_flow,
      monitors = tt_parametersToSave_OB_flow
    ),
    
    tt_Rmcmc_OB_flow = buildMCMC(tt_conf_OB_flow, useConjugacy = FALSE),
    tt_Cmodel_OB_flow = compileNimble(tt_Rmodel_OB_flow),
    tt_Cmcmc_OB_flow = compileNimble(tt_Rmcmc_OB_flow, project = tt_Rmodel_OB_flow),
    
    tt_model_OB_flow = runMCMC(
      tt_Cmcmc_OB_flow,
      niter = tt_runData_OB_flow$nIter,
      nburnin = tt_runData_OB_flow$nBurnin,
      thin = tt_runData_OB_flow$thinRate,
      nchains = tt_runData_OB_flow$nChains
    ),
    
    tt_modelOut_OB_flow =
      list(
        mcmc = tt_model_OB_flow, # "Error : invalid nimbleFunction argument"
        name = "phiT_pT_OB_flow",
        modelCode = tt_modelCode_OB_flow,
        myConstants = tt_myConstants_OB_flow,
        runData = tt_runData_OB_flow
      )#,
    
#    tt_save_OB_flow = saveModelOut_tt_OB_flow(tt_modelOut_OB_flow)
  )



    
    
######################################
#### Functions
#####################################
tt_initialValues_OB = function(t, y) {
  list(phi = runif(t - 1, 0, 1),
       p = runif(t - 1, 0, 1),
       z = getInits_tt_OB(y)
      )
  }

getInits_tt_OB <- function(d) {
  d <- d + 1
  d[d == 2] <- 1
  return(d)
}

saveModelOut_tt_OB <- function(d) {
  save(d, file = paste0('./models/cmrFlowOB/runsOut/tt_OB_', substr(Sys.time(),1,13), '.RData'))
}
