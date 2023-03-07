tar_option_set(packages = c("tidyverse", "nimble", "nimbleEcology", "MCMCvis", "targets"))

# This is a 'fork' of modelCMR_tt_ft_cohort_OB, following Daniel's March 2023 model update to the WB (4 river) model

# Run this to load data for testing
#eh_WB_2002_2014_target=tar_read(eh_wB_2002_2014_target)

# Model naming
# mod1 = ttt_ft_cohort_WB
# mod2 = ttt_fByRiverT_cohort_WB

#Source Daniel Turek's dDHMMo2 function, specifically made for this model to reduce model size and allow running on a laptop
source('./r/dDHMMo2.R')

#Change value to force 'restart' of target chains
dummyIn <- 12

#############################################################################
# Add flow
###################################################

modelCMR_ttt_ft_cohort_WB_flow_target <-
  tar_plan(

    mod1_dataIn_target = list(
      d = eh_WBbkt_2002_2014_target
    ),
    
    mod1_input_target =
      list(
        y = mod1_dataIn_target$d$eh * mod1_dataIn_target$d$riverN,
        first = mod1_dataIn_target$d$first, 
        last = mod1_dataIn_target$d$last, 
        zInits = ttt_initialValues_WB(ncol(mod1_dataIn_target$d$eh), 
                                      mod1_dataIn_target$d$eh),
        zInitsNA = ifelse(is.na(mod1_dataIn_target$d$flow), NA, 1),
        ##########################
        flow = mod1_dataIn_target$d$flow,
        ##########################
        temp = mod1_dataIn_target$d$temperature,
        nCohorts = length(unique(mod1_dataIn_target$d$cohorts$cohort)),
        cohort = mod1_dataIn_target$d$cohorts$cohort,
        nSeasons = length(unique(mod1_dataIn_target$d$data$season)),
        seasonArray = c(3,4,1,2,3,4,1,2,3,4,1,2),
        isYOY = mod1_dataIn_target$d$isYOY,
        # nStates = length(unique(mod1_dataIn_target$d$data$sizeState)),
        nRivers = length(unique(mod1_dataIn_target$d$data$riverN)),
        dummy = 0
      ), 
    
    mod1_runData_target = list(
      # Updateable model-specific variables 
      nIter = 12500, 
      nBurnin = 2500, 
      nChains = 2,
      thinRate = 5
    ),  
    
    mod1_myConstants0_target = list(
      N = nrow(mod1_input_target$y),
      T = ncol(mod1_input_target$y),
      first = mod1_input_target$first,
      last = mod1_input_target$last,
      cohort = mod1_input_target$cohort - min(mod1_input_target$cohort) + 1,
      nCohorts = mod1_input_target$nCohorts,
      season = mod1_input_target$seasonArray,

      flow = mod1_input_target$flow,

      temp = mod1_input_target$temp,
      length = mod1_input_target$last - mod1_input_target$first + 1,
      isYOY = mod1_input_target$isYOY,
      indToKeep = which(mod1_input_target$first < 12),
      y1 = mod1_input_target$y[,1] #y[,1],
    ),

    mod1_myConstants_target = list(
      N = length(mod1_myConstants0_target$indToKeep),
      T = ncol(mod1_input_target$y),
      
      first = mod1_myConstants0_target$first[mod1_myConstants0_target$indToKeep],
      last = mod1_myConstants0_target$last[mod1_myConstants0_target$indToKeep],
      cohort = mod1_myConstants0_target$cohort[mod1_myConstants0_target$indToKeep],
      flow = mod1_myConstants0_target$flow[mod1_myConstants0_target$indToKeep,],
      temp = mod1_myConstants0_target$temp[mod1_myConstants0_target$indToKeep,],
      length = mod1_myConstants0_target$length[mod1_myConstants0_target$indToKeep],
      isYOY = mod1_myConstants0_target$isYOY[mod1_myConstants0_target$indToKeep,],
      
      nCohorts = mod1_myConstants0_target$nCohorts,
      nSeasons = mod1_myConstants0_target$nSeasons,
      season =   mod1_myConstants0_target$season,
      
      alpha = c(1,1,1,1),
      #Priors for psi where more likely to stay than move
      alphaR1 = c(0.7, 0.1, 0.1, 0.1),
      alphaR2 = c(0.1, 0.7, 0.1, 0.1),
      alphaR3 = c(0.1, 0.1, 0.7, 0.1),
      alphaR4 = c(0.1, 0.1, 0.1, 0.7),
      
      alphaR = list(
        c(0.7, 0.1, 0.1, 0.1),
        c(0.1, 0.7, 0.1, 0.1),
        c(0.1, 0.1, 0.7, 0.1),
        c(0.1, 0.1, 0.1, 0.7)
      ),
      
      deltaProps = table(mod1_myConstants0_target$y1[mod1_myConstants0_target$y1>0]) / 
                  length(mod1_myConstants0_target$y1[mod1_myConstants0_target$y1>0])
    ),

    mod1_uniques_target = list(
      flow = unique(as.numeric(mod1_myConstants_target$flow)),
      temp = unique(as.numeric(mod1_myConstants_target$temp)),
      cohort = unique(as.numeric(mod1_myConstants_target$cohort)),
      season = unique(as.numeric(mod1_myConstants_target$season)),
      time = as.numeric(1:ncol(mod1_input_target$y)),
      
      cohortArray = matrix(mod1_myConstants_target$cohort, mod1_myConstants_target$N, mod1_myConstants_target$T),
      seasonArray = matrix(mod1_myConstants_target$season, mod1_myConstants_target$N, mod1_myConstants_target$T, byrow = TRUE)
      
    ),
    
    mod1_indexing_target = list(
       ##
     flowIndex = apply(mod1_myConstants_target$flow, c(1,2), function(x) which(x == mod1_uniques_target$flow)),
     tempIndex = apply(mod1_myConstants_target$temp, c(1,2), function(x) which(x == mod1_uniques_target$temp)),

     cohortIndex = apply(mod1_uniques_target$cohortArray, c(1,2), function(x) which(x == mod1_uniques_target$cohort)),
     seasonIndex = apply(mod1_uniques_target$seasonArray, c(1,2), function(x) which(x == mod1_uniques_target$season)),
     
     timeIndex = apply(
       matrix(mod1_uniques_target$time, mod1_myConstants_target$N, mod1_myConstants_target$T, byrow = TRUE), # timeArray
       c(1,2), function(x) which(x == mod1_uniques_target$time))
    ),
    
    mod1_indCombined_target = list(
      
     indCombined = paste0(
       as.numeric(mod1_indexing_target$flowIndex), '-',
       as.numeric(mod1_indexing_target$tempIndex), '-',
       as.numeric(mod1_indexing_target$cohortIndex), '-',
       as.numeric(mod1_indexing_target$seasonIndex), '-',
       as.numeric(mod1_indexing_target$timeIndex))
    ),
    
    mod1_uniqueIndCombined_target = list(
      uniqueIndCombined = unique(mod1_indCombined_target$indCombined)
    ), 

     ###length(uniqueIndCombined)    ## 2385 unique combindations of {flow,temp}
     ###                             ## 4425 unique combindations of {flow,temp,cohort}
     ###                             ## 4464 unique combindations of {flow,temp,cohort,season}
     ###                             ## 4541 unique combindationsof {flow,temp,cohort,season,time}
     ##
    
    mod1_uniqueIndSplitList_target = list(
      uniqueIndSplitList = strsplit(mod1_uniqueIndCombined_target$uniqueIndCombined, '-')
    ),
    
    mod1_kLooupTable_target = list(
     nK = length(mod1_uniqueIndCombined_target$uniqueIndCombined),

     flowK   = mod1_uniques_target$flow[as.numeric(sapply(mod1_uniqueIndSplitList_target$uniqueIndSplitList, `[`, 1))],
     tempK   = mod1_uniques_target$temp[as.numeric(sapply(mod1_uniqueIndSplitList_target$uniqueIndSplitList, `[`, 2))],
     cohortK = mod1_uniques_target$cohort[as.numeric(sapply(mod1_uniqueIndSplitList_target$uniqueIndSplitList, `[`, 3))],
     seasonK = mod1_uniques_target$season[as.numeric(sapply(mod1_uniqueIndSplitList_target$uniqueIndSplitList, `[`, 4))],
     timeK   = as.numeric(sapply(mod1_uniqueIndSplitList_target$uniqueIndSplitList, `[`, 5)),
     ##
     ## index this by kLookupTable[i,t], to get the value of k corresponding to {i,t}
     kLookupTable = array(match(mod1_indCombined_target$indCombined, mod1_uniqueIndCombined_target$uniqueIndCombined), 
                              c(mod1_myConstants_target$N, mod1_myConstants_target$T))
     
    ), 
    
    # not sure how to combine targets without getting a cyclic graph. Just making another target here
    mod1_constantsIn_target = list(
      N = mod1_myConstants_target$N,
      T = mod1_myConstants_target$T,
      
      first = mod1_myConstants_target$first,
      last = mod1_myConstants_target$last,
      cohort = mod1_myConstants_target$cohort,
      #flow = mod1_myConstants_target$flow,
      #temp = mod1_myConstants_target$temp,
      length = mod1_myConstants_target$length,
      isYOY = mod1_myConstants_target$isYOY,
      
      nRivers = mod1_input_target$nRivers,
      nCohorts = mod1_myConstants_target$nCohorts,
      #nSeasons = mod1_myConstants_target$nSeasons,
      #season =   mod1_myConstants_target$season,
      alpha = mod1_myConstants_target$alpha,
      #Priors for psi where more likely to stay than move
      alphaR1 = mod1_myConstants_target$alphaR1,
      alphaR2 = mod1_myConstants_target$alphaR2,
      alphaR3 = mod1_myConstants_target$alphaR3,
      alphaR4 = mod1_myConstants_target$alphaR4,
      
      alphaR = mod1_myConstants_target$alphaR,
      
      deltaProps = mod1_myConstants_target$deltaProps,
      
      nK = mod1_kLooupTable_target$nK,
      flowK = mod1_kLooupTable_target$flowK,
      tempK = mod1_kLooupTable_target$tempK,
      cohortK = mod1_kLooupTable_target$cohortK,
      seasonK = mod1_kLooupTable_target$seasonK,
      timeK = mod1_kLooupTable_target$timeK,
      kLookupTable = mod1_kLooupTable_target$kLookupTable,
      flow = NULL,
      temp = NULL,
      season = NULL,
      dummy = dummyIn
    ),
    
    #####################################
    # 
    mod1_myData_target = list(
      y = mod1_input_target$y[mod1_myConstants0_target$indToKeep,] + 1
    ),
    # 
    mod1_modelCode_target = nimbleCode({
      dummy  ~ dnorm(0,1)
      
      # Initial distribution among rivers
      delta[1] <- deltaProps[1]                  # Pr(alive t = 1 and in river 1) = 0.4
      delta[2] <- deltaProps[2]
      delta[3] <- deltaProps[3]
      delta[4] <- deltaProps[4]
      delta[5] <- 0                    # Pr(dead t = 1) = 0
      ##
      for (t in 1:(T-1)){ # loop over time
        for (c in 1:nCohorts){
          psi[1,1:nRivers,t,c] ~ ddirch(alphaR1[1:nRivers])
          psi[2,1:nRivers,t,c] ~ ddirch(alphaR2[1:nRivers])
          psi[3,1:nRivers,t,c] ~ ddirch(alphaR3[1:nRivers])
          psi[4,1:nRivers,t,c] ~ ddirch(alphaR4[1:nRivers])
        }
      }
      ##
      for (r in 1:nRivers){
        for (s in 1:4){
          for (c in 1:nCohorts){
            betaInt[r,s,c] ~ dnorm(0,1)
          }
          betaFlow[1,r,s] ~ dnorm(0,1)
          betaFlow[2,r,s] ~ dnorm(0,1)
          betaFlow[3,r,s] ~ dnorm(0,1)
        }
      }
      for (r in 1:nRivers){
        for (t in 1:(T-1)){ # loop over time
          for (c in 1:nCohorts){
            betaP[r,t,c] ~ dnorm(0,1)
          }
        }
      }
      ##
      ## NEW:
      for(k in 1:nK) {
        for(r in 1:nRivers) {
          # need to add isYOY
          logit(phi[r,k]) <- 
            betaInt [  r,seasonK[k],cohortK[k]] +
            betaFlow[1,r,seasonK[k]           ] * flowK[k] +
            betaFlow[2,r,seasonK[k]           ] * tempK[k] +
            betaFlow[3,r,seasonK[k]           ] * tempK[k] * flowK[k]
        }
      }
      ##
      ## NEW:
      for(k in 1:nK) {
        gamma[1,1,k] <- phi[1,k] * psi[1,1,timeK[k],cohortK[k]]
        gamma[1,2,k] <- phi[1,k] * psi[1,2,timeK[k],cohortK[k]]
        gamma[1,3,k] <- phi[1,k] * psi[1,3,timeK[k],cohortK[k]]
        gamma[1,4,k] <- phi[1,k] * psi[1,4,timeK[k],cohortK[k]]
        gamma[1,5,k] <- 1 - phi[1,k]
        gamma[2,1,k] <- phi[2,k] * psi[2,1,timeK[k],cohortK[k]]
        gamma[2,2,k] <- phi[2,k] * psi[2,2,timeK[k],cohortK[k]]
        gamma[2,3,k] <- phi[2,k] * psi[2,3,timeK[k],cohortK[k]]
        gamma[2,4,k] <- phi[2,k] * psi[2,4,timeK[k],cohortK[k]]
        gamma[2,5,k] <- 1 - phi[2,k]
        gamma[3,1,k] <- phi[3,k] * psi[3,1,timeK[k],cohortK[k]]
        gamma[3,2,k] <- phi[3,k] * psi[3,2,timeK[k],cohortK[k]]
        gamma[3,3,k] <- phi[3,k] * psi[3,3,timeK[k],cohortK[k]]
        gamma[3,4,k] <- phi[3,k] * psi[3,4,timeK[k],cohortK[k]]
        gamma[3,5,k] <- 1 - phi[3,k]
        gamma[4,1,k] <- phi[4,k] * psi[4,1,timeK[k],cohortK[k]]  ## 0.2
        gamma[4,2,k] <- phi[4,k] * psi[4,2,timeK[k],cohortK[k]]  ## 0.2
        gamma[4,3,k] <- phi[4,k] * psi[4,3,timeK[k],cohortK[k]]  ## 0.2
        gamma[4,4,k] <- phi[4,k] * psi[4,4,timeK[k],cohortK[k]]  ## 0.2
        gamma[4,5,k] <- 1 - phi[4,k]                             ## 0.2
        gamma[5,1,k] <- 0
        gamma[5,2,k] <- 0
        gamma[5,3,k] <- 0
        gamma[5,4,k] <- 0
        gamma[5,5,k] <- 1
        ###### gamma for the last occasion
        ####for (a in 1:(nRivers+1)){
        ####    for (b in 1:nRivers){
        ####        gamma[a,b,T,k] <- 0
        ####    }
        ####    gamma[a,5,T,k] <- 1
        ####}
      }
      ##
      ##
      ## NEW: omega1 for first obs
      omega1[1,1] <- 0          # Pr(alive A t -> non-detected t)
      omega1[1,2] <- 1          # Pr(alive A t -> detected A t)
      omega1[1,3] <- 0          # Pr(alive A t -> detected B t)
      omega1[1,4] <- 0          # Pr(alive A t -> detected C t)
      omega1[1,5] <- 0          # Pr(alive A t -> detected D t)
      omega1[2,1] <- 0          # Pr(alive B t -> non-detected t)
      omega1[2,2] <- 0          # Pr(alive B t -> detected A t)
      omega1[2,3] <- 1          # Pr(alive B t -> detected B t)
      omega1[2,4] <- 0          # Pr(alive B t -> detected C t)
      omega1[2,5] <- 0          # Pr(alive B t -> detected C t)
      omega1[3,1] <- 0          # Pr(alive C t -> non-detected t)
      omega1[3,2] <- 0          # Pr(alive C t -> detected A t)
      omega1[3,3] <- 0          # Pr(alive C t -> detected B t)
      omega1[3,4] <- 1          # Pr(alive C t -> detected C t)
      omega1[3,5] <- 0          # Pr(alive C t -> detected C t)
      omega1[4,1] <- 0          # Pr(dead t -> non-detected t)
      omega1[4,2] <- 0          # Pr(dead t -> detected A t)
      omega1[4,3] <- 0          # Pr(dead t -> detected B t)
      omega1[4,4] <- 0          # Pr(dead t -> detected C t)
      omega1[4,5] <- 1          # Pr(dead t -> detected C t)
      omega1[5,1] <- 1          # Pr(dead t -> non-detected t)
      omega1[5,2] <- 0          # Pr(dead t -> detected A t)
      omega1[5,3] <- 0          # Pr(dead t -> detected B t)
      omega1[5,4] <- 0          # Pr(dead t -> detected C t)
      omega1[5,5] <- 0          # Pr(dead t -> detected D t)
      ##
      ##
      ## NEW: omega[] elements are assigned for each *cohort*, rather than each *individual*
      ##      same for pA, pB, pC, pD.
      for (t in 2:T){ # loop over time
        for (c in 1:nCohorts){
          logit(pA[t,c]) <- betaP[1,t-1,c]
          logit(pB[t,c]) <- betaP[2,t-1,c]
          logit(pC[t,c]) <- betaP[3,t-1,c] 
          logit(pD[t,c]) <- betaP[4,t-1,c]
          ##
          # probabilities of y(t) given z(t)
          # omega[z, y, t, i]
          ##
          # z=1 = alive in River 1, z=2 = alive in River 2...z=5 = dead
          # y=1 = unobserved, y=2 = observed in River 1, y=3 = observed in River 2, etc
          ##
          omega[1,1,t,c] <- 1 - pA[t,c]     # Pr(alive A t -> non-detected t)
          omega[1,2,t,c] <- pA[t,c]         # Pr(alive A t -> detected A t)
          omega[1,3,t,c] <- 0               # Pr(alive A t -> detected B t)
          omega[1,4,t,c] <- 0               # Pr(alive A t -> detected C t)
          omega[1,5,t,c] <- 0               # Pr(alive A t -> detected D t)
          omega[2,1,t,c] <- 1 - pB[t,c]     # Pr(alive B t -> non-detected t)
          omega[2,2,t,c] <- 0               # Pr(alive B t -> detected A t)
          omega[2,3,t,c] <- pB[t,c]         # Pr(alive B t -> detected B t)
          omega[2,4,t,c] <- 0               # Pr(alive B t -> detected C t)
          omega[2,5,t,c] <- 0               # Pr(alive B t -> detected C t)
          omega[3,1,t,c] <- 1 - pC[t,c]     # Pr(alive C t -> non-detected t)
          omega[3,2,t,c] <- 0               # Pr(alive C t -> detected A t)
          omega[3,3,t,c] <- 0               # Pr(alive C t -> detected B t)
          omega[3,4,t,c] <- pC[t,c]         # Pr(alive C t -> detected C t)
          omega[3,5,t,c] <- 0               # Pr(alive C t -> detected C t)
          omega[4,1,t,c] <- 1 - pD[t,c]     # Pr(alive D t -> non-detected t))
          omega[4,2,t,c] <- 0               # Pr(dead D t -> detected A t)
          omega[4,3,t,c] <- 0               # Pr(dead D t -> detected B t)
          omega[4,4,t,c] <- 0               # Pr(dead D t -> detected C t)
          omega[4,5,t,c] <- pD[t,c]         # Pr(alive D t -> detected D t)
          omega[5,1,t,c] <- 1               # Pr(dead t -> non-detected t)
          omega[5,2,t,c] <- 0               # Pr(dead t -> detected A t)
          omega[5,3,t,c] <- 0               # Pr(dead t -> detected B t)
          omega[5,4,t,c] <- 0               # Pr(dead t -> detected C t)
          omega[5,5,t,c] <- 0               # Pr(dead t -> detected D t)
        }
      }
      ##
      ##
      for (i in 1:N) {
        y[i,first[i]:last[i]] ~ dDHMMo2(init = delta[1:5],
                                        f = first[i],
                                        probTrans1  = gamma[1:5, 1:5, kLookupTable[i,1]],
                                        probTrans2  = gamma[1:5, 1:5, kLookupTable[i,2]],
                                        probTrans3  = gamma[1:5, 1:5, kLookupTable[i,3]],
                                        probTrans4  = gamma[1:5, 1:5, kLookupTable[i,4]],
                                        probTrans5  = gamma[1:5, 1:5, kLookupTable[i,5]],
                                        probTrans6  = gamma[1:5, 1:5, kLookupTable[i,6]],
                                        probTrans7  = gamma[1:5, 1:5, kLookupTable[i,7]],
                                        probTrans8  = gamma[1:5, 1:5, kLookupTable[i,8]],
                                        probTrans9  = gamma[1:5, 1:5, kLookupTable[i,9]],
                                        probTrans10 = gamma[1:5, 1:5, kLookupTable[i,10]],
                                        probTrans11 = gamma[1:5, 1:5, kLookupTable[i,11]],
                                        probTrans12 = gamma[1:5, 1:5, kLookupTable[i,12]],
                                        probObs1 = omega1[1:5, 1:5],
                                        probObs = omega[1:5, 1:5, first[i]:last[i], cohort[i]],
                                        len = length[i],
                                        checkRowSums = 1)
      }
    }),
    
    mod1_Rmodel_target = nimbleModel(
      code = mod1_modelCode_target,
      constants = mod1_constantsIn_target,
      data = mod1_myData_target,
      inits = mod_initialValues(mod1_constantsIn_target$nRivers, mod1_constantsIn_target$nCohorts,
                                mod1_constantsIn_target$nK, mod1_constantsIn_target$T,
                                mod1_constantsIn_target, mod1_constantsIn_target$alphaR),
      calculate = FALSE
    ),
    # 
    mod1_parametersToSave_target = c(
      "phi", "betaP", "psi"
    ),

    mod1_conf_target = configureMCMC(
      tar_read(mod1_Rmodel_target),
      monitors = mod1_parametersToSave_target,
      useConjugacy = FALSE
    )
    # 
    # mod1_Rmcmc_target = buildMCMC(mod1_conf_target, useConjugacy = FALSE),
    # mod1_Cmodel_target = compileNimble(mod1_Rmodel_target),
    # mod1_Cmcmc_target = compileNimble(mod1_Rmcmc_target, project = mod1_Rmodel_target),
    # 
    # mod1_model_target = runMCMC(
    #   mod1_Cmcmc_target,
    #   niter = mod1_runData_target$nIter,
    #   nburnin = mod1_runData_target$nBurnin,
    #   thin = mod1_runData_target$thinRate,
    #   nchains = mod1_runData_target$nChains
    # ),
    # 
    # mod1_modelOut_target =
    #   list(
    #     mcmc = mod1_model_target, # "Error : invalid nimbleFunction argument"
    #     name = "mod1",
    #     modelCode = mod1_modelCode_target,
    #     myConstants = mod1_constantsIn_target,
    #     runData = mod1_runData_target
    #   ),
    # 
    # mod1_save_target = saveModelOut_mod1(mod1_modelOut_target)
   )
# 
# 
# # #############################################################################
# # # Switch to flow by river
# # # Flow estimated for each river independently
# # ###################################################
modelCMR_ttt_ft_cohort_WB_flowByRiver_target <-
  tar_plan(

    mod2_input_target =
      list(
        y = mod1_dataIn_target$d$eh * mod1_dataIn_target$d$riverN,
        first = mod1_dataIn_target$d$first,
        last = mod1_dataIn_target$d$last,
        zInits = ttt_initialValues_WB(ncol(mod1_dataIn_target$d$eh),
                                      mod1_dataIn_target$d$eh),
        zInitsNA = ifelse(is.na(mod1_dataIn_target$d$flow), NA, 1),
        ##########################
        flow = mod1_dataIn_target$d$flowByRiver,
        ##########################
        temp = mod1_dataIn_target$d$temperature,
        nCohorts = length(unique(mod1_dataIn_target$d$cohorts$cohort)),
        cohort = mod1_dataIn_target$d$cohorts$cohort,
        nSeasons = length(unique(mod1_dataIn_target$d$data$season)),
        seasonArray = c(3,4,1,2,3,4,1,2,3,4,1,2),
        isYOY = mod1_dataIn_target$d$isYOY,
        # nStates = length(unique(mod1_dataIn_target$d$data$sizeState)),
        nRivers = length(unique(mod1_dataIn_target$d$data$riverN)),
        dummy = 0
      ),

    mod2_runData_target = list(
      # Updateable model-specific variables
      nIter = 12500,
      nBurnin = 2500,
      nChains = 2,
      thinRate = 5
    ),

    mod2_myConstants0_target = list(
      N = nrow(mod2_input_target$y),
      T = ncol(mod2_input_target$y),
      first = mod2_input_target$first,
      last = mod2_input_target$last,
      cohort = mod2_input_target$cohort - min(mod2_input_target$cohort) + 1,
      nCohorts = mod2_input_target$nCohorts,
      season = mod2_input_target$seasonArray,

      flow = mod2_input_target$flow,

      temp = mod2_input_target$temp,
      length = mod2_input_target$last - mod2_input_target$first + 1,
      isYOY = mod2_input_target$isYOY,
      indToKeep = which(mod2_input_target$first < 12),
      y1 = mod2_input_target$y[,1] #y[,1],
    ),

    mod2_myConstants_target = list(
      N = length(mod2_myConstants0_target$indToKeep),
      T = ncol(mod2_input_target$y),

      first = mod2_myConstants0_target$first[mod2_myConstants0_target$indToKeep],
      last = mod2_myConstants0_target$last[mod2_myConstants0_target$indToKeep],
      cohort = mod2_myConstants0_target$cohort[mod2_myConstants0_target$indToKeep],
      flow = mod2_myConstants0_target$flow[mod2_myConstants0_target$indToKeep,],
      temp = mod2_myConstants0_target$temp[mod2_myConstants0_target$indToKeep,],
      length = mod2_myConstants0_target$length[mod2_myConstants0_target$indToKeep],
      isYOY = mod2_myConstants0_target$isYOY[mod2_myConstants0_target$indToKeep,],

      nCohorts = mod2_myConstants0_target$nCohorts,
      nSeasons = mod2_myConstants0_target$nSeasons,
      season =   mod2_myConstants0_target$season,
      alpha = c(1,1,1,1),
      #Priors for psi where more likely to stay than move
      alphaR1 = c(0.7, 0.1, 0.1, 0.1),
      alphaR2 = c(0.1, 0.7, 0.1, 0.1),
      alphaR3 = c(0.1, 0.1, 0.7, 0.1),
      alphaR4 = c(0.1, 0.1, 0.1, 0.7),

      alphaR = list(
        c(0.7, 0.1, 0.1, 0.1),
        c(0.1, 0.7, 0.1, 0.1),
        c(0.1, 0.1, 0.7, 0.1),
        c(0.1, 0.1, 0.1, 0.7)
      ),
      
      deltaProps = table(mod2_myConstants0_target$y1[mod2_myConstants0_target$y1>0]) /
        length(mod2_myConstants0_target$y1[mod2_myConstants0_target$y1>0])
    ),

    mod2_uniques_target = list(
      flow = unique(as.numeric(mod2_myConstants_target$flow)),
      temp = unique(as.numeric(mod2_myConstants_target$temp)),
      cohort = unique(as.numeric(mod2_myConstants_target$cohort)),
      season = unique(as.numeric(mod2_myConstants_target$season)),
      time = as.numeric(1:ncol(mod2_input_target$y)),

      cohortArray = matrix(mod2_myConstants_target$cohort, mod2_myConstants_target$N, mod2_myConstants_target$T),
      seasonArray = matrix(mod2_myConstants_target$season, mod2_myConstants_target$N, mod2_myConstants_target$T, byrow = TRUE)

    ),

    mod2_indexing_target = list(
      ##
      flowIndex = apply(mod2_myConstants_target$flow, c(1,2), function(x) which(x == mod2_uniques_target$flow)),
      tempIndex = apply(mod2_myConstants_target$temp, c(1,2), function(x) which(x == mod2_uniques_target$temp)),

      cohortIndex = apply(mod2_uniques_target$cohortArray, c(1,2), function(x) which(x == mod2_uniques_target$cohort)),
      seasonIndex = apply(mod2_uniques_target$seasonArray, c(1,2), function(x) which(x == mod2_uniques_target$season)),

      timeIndex = apply(
        matrix(mod2_uniques_target$time, mod2_myConstants_target$N, mod2_myConstants_target$T, byrow = TRUE), # timeArray
        c(1,2), function(x) which(x == mod2_uniques_target$time))
    ),

    mod2_indCombined_target = list(

      indCombined = paste0(
        as.numeric(mod2_indexing_target$flowIndex), '-',
        as.numeric(mod2_indexing_target$tempIndex), '-',
        as.numeric(mod2_indexing_target$cohortIndex), '-',
        as.numeric(mod2_indexing_target$seasonIndex), '-',
        as.numeric(mod2_indexing_target$timeIndex))
    ),

    mod2_uniqueIndCombined_target = list(
      uniqueIndCombined = unique(mod2_indCombined_target$indCombined)
    ),

    ###length(uniqueIndCombined)    ## 2385 unique combindations of {flow,temp}
    ###                             ## 4425 unique combindations of {flow,temp,cohort}
    ###                             ## 4464 unique combindations of {flow,temp,cohort,season}
    ###                             ## 4541 unique combindationsof {flow,temp,cohort,season,time}
    ##

    mod2_uniqueIndSplitList_target = list(
      uniqueIndSplitList = strsplit(mod2_uniqueIndCombined_target$uniqueIndCombined, '-')
    ),

    mod2_kLooupTable_target = list(
      nK = length(mod2_uniqueIndCombined_target$uniqueIndCombined),

      flowK   = mod2_uniques_target$flow[as.numeric(sapply(mod2_uniqueIndSplitList_target$uniqueIndSplitList, `[`, 1))],
      tempK   = mod2_uniques_target$temp[as.numeric(sapply(mod2_uniqueIndSplitList_target$uniqueIndSplitList, `[`, 2))],
      cohortK = mod2_uniques_target$cohort[as.numeric(sapply(mod2_uniqueIndSplitList_target$uniqueIndSplitList, `[`, 3))],
      seasonK = mod2_uniques_target$season[as.numeric(sapply(mod2_uniqueIndSplitList_target$uniqueIndSplitList, `[`, 4))],
      timeK   = as.numeric(sapply(mod2_uniqueIndSplitList_target$uniqueIndSplitList, `[`, 5)),
      ##
      ## index this by kLookupTable[i,t], to get the value of k corresponding to {i,t}
      kLookupTable = array(match(mod2_indCombined_target$indCombined, mod2_uniqueIndCombined_target$uniqueIndCombined),
                           c(mod2_myConstants_target$N, mod2_myConstants_target$T))

    ),

    # not sure how to combine targets without getting a cyclic graph. Just making another target here
    mod2_constantsIn_target = list(
      N = mod2_myConstants_target$N,
      T = mod2_myConstants_target$T,

      first = mod2_myConstants_target$first,
      last = mod2_myConstants_target$last,
      cohort = mod2_myConstants_target$cohort,
      #flow = mod2_myConstants_target$flow,
      #temp = mod2_myConstants_target$temp,
      length = mod2_myConstants_target$length,
      isYOY = mod2_myConstants_target$isYOY,

      nRivers = mod2_input_target$nRivers,
      nCohorts = mod2_myConstants_target$nCohorts,
      #nSeasons = mod2_myConstants_target$nSeasons,
      #season =   mod2_myConstants_target$season,
      #alpha = mod2_myConstants_target$alpha,
      #Priors for psi where more likely to stay than move
      alphaR1 = mod2_myConstants_target$alphaR1,
      alphaR2 = mod2_myConstants_target$alphaR2,
      alphaR3 = mod2_myConstants_target$alphaR3,
      alphaR4 = mod2_myConstants_target$alphaR4,

      alphaR = mod2_myConstants_target$alphaR,
      
      deltaProps = mod2_myConstants_target$deltaProps,

      nK = mod2_kLooupTable_target$nK,
      flowK = mod2_kLooupTable_target$flowK,
      tempK = mod2_kLooupTable_target$tempK,
      cohortK = mod2_kLooupTable_target$cohortK,
      seasonK = mod2_kLooupTable_target$seasonK,
      timeK = mod2_kLooupTable_target$timeK,
      kLookupTable = mod2_kLooupTable_target$kLookupTable,
      flow = NULL,
      temp = NULL,
      season = NULL,
      dummy = dummyIn
    ),

    #####################################
    #
    mod2_myData_target = list(
      y = mod2_input_target$y[mod2_myConstants0_target$indToKeep,] + 1
    ),
    #
    mod2_modelCode_target = nimbleCode({
      dummy  ~ dnorm(0,1)

      # Initial distribution among rivers
      delta[1] <- deltaProps[1]                  # Pr(alive t = 1 and in river 1) = 0.4
      delta[2] <- deltaProps[2]
      delta[3] <- deltaProps[3]
      delta[4] <- deltaProps[4]
      delta[5] <- 0                    # Pr(dead t = 1) = 0
      ##
      for (t in 1:(T-1)){ # loop over time
        for (c in 1:nCohorts){
          psi[1,1:nRivers,t,c] ~ ddirch(alphaR1[1:nRivers])
          psi[2,1:nRivers,t,c] ~ ddirch(alphaR2[1:nRivers])
          psi[3,1:nRivers,t,c] ~ ddirch(alphaR3[1:nRivers])
          psi[4,1:nRivers,t,c] ~ ddirch(alphaR4[1:nRivers])
        }
      }
      ##
      for (r in 1:nRivers){
        for (s in 1:4){
          for (c in 1:nCohorts){
            betaInt[r,s,c] ~ dnorm(0,1)
          }
          betaFlow[1,r,s] ~ dnorm(0,1)
          betaFlow[2,r,s] ~ dnorm(0,1)
          betaFlow[3,r,s] ~ dnorm(0,1)
        }
      }
      for (r in 1:nRivers){
        for (t in 1:(T-1)){ # loop over time
          for (c in 1:nCohorts){
            betaP[r,t,c] ~ dnorm(0,1)
          }
        }
      }
      ##
      ## NEW:
      for(k in 1:nK) {
        for(r in 1:nRivers) {
          # need to add isYOY
          logit(phi[r,k]) <-
            betaInt [  r,seasonK[k],cohortK[k]] +
            betaFlow[1,r,seasonK[k]           ] * flowK[k] +
            betaFlow[2,r,seasonK[k]           ] * tempK[k] +
            betaFlow[3,r,seasonK[k]           ] * tempK[k] * flowK[k]
        }
      }
      ##
      ## NEW:
      for(k in 1:nK) {
        gamma[1,1,k] <- phi[1,k] * psi[1,1,timeK[k],cohortK[k]]
        gamma[1,2,k] <- phi[1,k] * psi[1,2,timeK[k],cohortK[k]]
        gamma[1,3,k] <- phi[1,k] * psi[1,3,timeK[k],cohortK[k]]
        gamma[1,4,k] <- phi[1,k] * psi[1,4,timeK[k],cohortK[k]]
        gamma[1,5,k] <- 1 - phi[1,k]
        gamma[2,1,k] <- phi[2,k] * psi[2,1,timeK[k],cohortK[k]]
        gamma[2,2,k] <- phi[2,k] * psi[2,2,timeK[k],cohortK[k]]
        gamma[2,3,k] <- phi[2,k] * psi[2,3,timeK[k],cohortK[k]]
        gamma[2,4,k] <- phi[2,k] * psi[2,4,timeK[k],cohortK[k]]
        gamma[2,5,k] <- 1 - phi[2,k]
        gamma[3,1,k] <- phi[3,k] * psi[3,1,timeK[k],cohortK[k]]
        gamma[3,2,k] <- phi[3,k] * psi[3,2,timeK[k],cohortK[k]]
        gamma[3,3,k] <- phi[3,k] * psi[3,3,timeK[k],cohortK[k]]
        gamma[3,4,k] <- phi[3,k] * psi[3,4,timeK[k],cohortK[k]]
        gamma[3,5,k] <- 1 - phi[3,k]
        gamma[4,1,k] <- phi[4,k] * psi[4,1,timeK[k],cohortK[k]]  ## 0.2
        gamma[4,2,k] <- phi[4,k] * psi[4,2,timeK[k],cohortK[k]]  ## 0.2
        gamma[4,3,k] <- phi[4,k] * psi[4,3,timeK[k],cohortK[k]]  ## 0.2
        gamma[4,4,k] <- phi[4,k] * psi[4,4,timeK[k],cohortK[k]]  ## 0.2
        gamma[4,5,k] <- 1 - phi[4,k]                             ## 0.2
        gamma[5,1,k] <- 0
        gamma[5,2,k] <- 0
        gamma[5,3,k] <- 0
        gamma[5,4,k] <- 0
        gamma[5,5,k] <- 1
        ###### gamma for the last occasion
        ####for (a in 1:(nRivers+1)){
        ####    for (b in 1:nRivers){
        ####        gamma[a,b,T,k] <- 0
        ####    }
        ####    gamma[a,5,T,k] <- 1
        ####}
      }
      ##
      ##
      ## NEW: omega1 for first obs
      omega1[1,1] <- 0          # Pr(alive A t -> non-detected t)
      omega1[1,2] <- 1          # Pr(alive A t -> detected A t)
      omega1[1,3] <- 0          # Pr(alive A t -> detected B t)
      omega1[1,4] <- 0          # Pr(alive A t -> detected C t)
      omega1[1,5] <- 0          # Pr(alive A t -> detected D t)
      omega1[2,1] <- 0          # Pr(alive B t -> non-detected t)
      omega1[2,2] <- 0          # Pr(alive B t -> detected A t)
      omega1[2,3] <- 1          # Pr(alive B t -> detected B t)
      omega1[2,4] <- 0          # Pr(alive B t -> detected C t)
      omega1[2,5] <- 0          # Pr(alive B t -> detected C t)
      omega1[3,1] <- 0          # Pr(alive C t -> non-detected t)
      omega1[3,2] <- 0          # Pr(alive C t -> detected A t)
      omega1[3,3] <- 0          # Pr(alive C t -> detected B t)
      omega1[3,4] <- 1          # Pr(alive C t -> detected C t)
      omega1[3,5] <- 0          # Pr(alive C t -> detected C t)
      omega1[4,1] <- 0          # Pr(dead t -> non-detected t)
      omega1[4,2] <- 0          # Pr(dead t -> detected A t)
      omega1[4,3] <- 0          # Pr(dead t -> detected B t)
      omega1[4,4] <- 0          # Pr(dead t -> detected C t)
      omega1[4,5] <- 1          # Pr(dead t -> detected C t)
      omega1[5,1] <- 1          # Pr(dead t -> non-detected t)
      omega1[5,2] <- 0          # Pr(dead t -> detected A t)
      omega1[5,3] <- 0          # Pr(dead t -> detected B t)
      omega1[5,4] <- 0          # Pr(dead t -> detected C t)
      omega1[5,5] <- 0          # Pr(dead t -> detected D t)
      ##
      ##
      ## NEW: omega[] elements are assigned for each *cohort*, rather than each *individual*
      ##      same for pA, pB, pC, pD.
      for (t in 2:T){ # loop over time
        for (c in 1:nCohorts){
          logit(pA[t,c]) <- betaP[1,t-1,c]
          logit(pB[t,c]) <- betaP[2,t-1,c]
          logit(pC[t,c]) <- betaP[3,t-1,c]
          logit(pD[t,c]) <- betaP[4,t-1,c]
          ##
          # probabilities of y(t) given z(t)
          # omega[z, y, t, i]
          ##
          # z=1 = alive in River 1, z=2 = alive in River 2...z=5 = dead
          # y=1 = unobserved, y=2 = observed in River 1, y=3 = observed in River 2, etc
          ##
          omega[1,1,t,c] <- 1 - pA[t,c]     # Pr(alive A t -> non-detected t)
          omega[1,2,t,c] <- pA[t,c]         # Pr(alive A t -> detected A t)
          omega[1,3,t,c] <- 0               # Pr(alive A t -> detected B t)
          omega[1,4,t,c] <- 0               # Pr(alive A t -> detected C t)
          omega[1,5,t,c] <- 0               # Pr(alive A t -> detected D t)
          omega[2,1,t,c] <- 1 - pB[t,c]     # Pr(alive B t -> non-detected t)
          omega[2,2,t,c] <- 0               # Pr(alive B t -> detected A t)
          omega[2,3,t,c] <- pB[t,c]         # Pr(alive B t -> detected B t)
          omega[2,4,t,c] <- 0               # Pr(alive B t -> detected C t)
          omega[2,5,t,c] <- 0               # Pr(alive B t -> detected C t)
          omega[3,1,t,c] <- 1 - pC[t,c]     # Pr(alive C t -> non-detected t)
          omega[3,2,t,c] <- 0               # Pr(alive C t -> detected A t)
          omega[3,3,t,c] <- 0               # Pr(alive C t -> detected B t)
          omega[3,4,t,c] <- pC[t,c]         # Pr(alive C t -> detected C t)
          omega[3,5,t,c] <- 0               # Pr(alive C t -> detected C t)
          omega[4,1,t,c] <- 1 - pD[t,c]     # Pr(alive D t -> non-detected t))
          omega[4,2,t,c] <- 0               # Pr(dead D t -> detected A t)
          omega[4,3,t,c] <- 0               # Pr(dead D t -> detected B t)
          omega[4,4,t,c] <- 0               # Pr(dead D t -> detected C t)
          omega[4,5,t,c] <- pD[t,c]         # Pr(alive D t -> detected D t)
          omega[5,1,t,c] <- 1               # Pr(dead t -> non-detected t)
          omega[5,2,t,c] <- 0               # Pr(dead t -> detected A t)
          omega[5,3,t,c] <- 0               # Pr(dead t -> detected B t)
          omega[5,4,t,c] <- 0               # Pr(dead t -> detected C t)
          omega[5,5,t,c] <- 0               # Pr(dead t -> detected D t)
        }
      }
      ##
      ##
      for (i in 1:N) {
        y[i,first[i]:last[i]] ~ dDHMMo2(init = delta[1:5],
                                        f = first[i],
                                        probTrans1  = gamma[1:5, 1:5, kLookupTable[i,1]],
                                        probTrans2  = gamma[1:5, 1:5, kLookupTable[i,2]],
                                        probTrans3  = gamma[1:5, 1:5, kLookupTable[i,3]],
                                        probTrans4  = gamma[1:5, 1:5, kLookupTable[i,4]],
                                        probTrans5  = gamma[1:5, 1:5, kLookupTable[i,5]],
                                        probTrans6  = gamma[1:5, 1:5, kLookupTable[i,6]],
                                        probTrans7  = gamma[1:5, 1:5, kLookupTable[i,7]],
                                        probTrans8  = gamma[1:5, 1:5, kLookupTable[i,8]],
                                        probTrans9  = gamma[1:5, 1:5, kLookupTable[i,9]],
                                        probTrans10 = gamma[1:5, 1:5, kLookupTable[i,10]],
                                        probTrans11 = gamma[1:5, 1:5, kLookupTable[i,11]],
                                        probTrans12 = gamma[1:5, 1:5, kLookupTable[i,12]],
                                        probObs1 = omega1[1:5, 1:5],
                                        probObs = omega[1:5, 1:5, first[i]:last[i], cohort[i]],
                                        len = length[i],
                                        checkRowSums = 1)
      }
    }),

    mod2_Rmodel_target = nimbleModel(
      code = mod2_modelCode_target,
      constants = mod2_constantsIn_target,
      data = mod2_myData_target,
      inits = mod_initialValues(mod2_constantsIn_target$nRivers, mod2_constantsIn_target$nCohorts,
                                mod2_constantsIn_target$nK, mod2_constantsIn_target$T,
                                mod2_constantsIn_target, mod2_constantsIn_target$alphaR),
      calculate = FALSE
    ),
    # 
    # 
    mod2_parametersToSave_target = c(
      "phi", "betaP", "psi"
    ),

    mod2_conf_target = configureMCMC(
      tar_read(mod2_Rmodel_target),
      monitors = mod2_parametersToSave_target,
      useConjugacy = FALSE
    )
    # 
    # mod2_Rmcmc_target = buildMCMC(mod2_conf_target, useConjugacy = FALSE),
    # mod2_Cmodel_target = compileNimble(mod2_Rmodel_target),
    # mod2_Cmcmc_target = compileNimble(mod2_Rmcmc_target, project = mod2_Rmodel_target),
    # 
    # mod2_model_target = runMCMC(
    #   mod2_Cmcmc_target,
    #   niter = mod2_runData_target$nIter,
    #   nburnin = mod2_runData_target$nBurnin,
    #   thin = mod2_runData_target$thinRate,
    #   nchains = mod2_runData_target$nChains
    # ),
    # 
    # mod2_modelOut_target =
    #   list(
    #     mcmc = mod2_model_target, # "Error : invalid nimbleFunction argument"
    #     name = "mod2",
    #     modelCode = mod2_modelCode_target,
    #     myConstants = mod2_constantsIn_target,
    #     runData = mod2_runData_target
    #   ),
    # 
    # mod2_save_target = saveModelOut_mod2(mod2_modelOut_target)
  )

#     
######################################
#### Functions
#####################################
ttt_initialValues_WB = function(t, y) {
  list(phi = runif(t - 1, 0, 1),
       p = runif(t - 1, 0, 1),
       z = getInits_ttt_ft_cohort_WB(y)
      )
}



getInits_ttt_ft_cohort_WB <- function(d) {
  d <- d + 1
  d[d == 2] <- 1
  return(d)
}

saveModelOut_tt_ft_cohort_WB <- function(d) {
  save(d, file = paste0('./models/cmrFlowOB/runsOut/tt_ft_cohort_WB_', substr(Sys.time(),1,13), '.RData'))
}
saveModelOut_mod1 <- function(d) {
  save(d, file = paste0('./models/cmrFlowWB/runsOut/mod1', substr(Sys.time(),1,13), '.RData'))
}

saveModelOut_mod2 <- function(d) {
  save(d, file = paste0('./models/cmrFlowOB/runsOut/mod2', substr(Sys.time(),1,13), '.RData'))
}

# fill in entries for dirichlet priors where a[r,,t,c] sums to 1 for any r,t,c combo
getDirchPriors <- function(nRivers,myConstants, nCohorts, alpha){
  ##a = array(rep(0, nRivers * nRivers * (myConstants$T - 1) * nCohorts) , c(nRivers, nRivers, (myConstants$T - 1), nCohorts ))
  a = array(rep(0, nRivers * nRivers * myConstants$T * nCohorts) , c(nRivers, nRivers, myConstants$T, nCohorts ))
  for(r in 1:nRivers){
    ##for(t in 1:(myConstants$T - 1)){
    for(t in 1:myConstants$T){
      for (c in 1:nCohorts){
        dirch <- rdirch(1, alpha)
        for (r2 in 1:nRivers){
          a[r,r2,t,c] <- dirch[r2]
        }
      }
    }
  }
  return(a)
}

getDirchPriorsR <- function(nRivers, myConstants, nCohorts, alphaR){
  ##a = array(rep(0, nRivers * nRivers * (myConstants$T - 1) * nCohorts) , c(nRivers, nRivers, (myConstants$T - 1), nCohorts ))
  a = array(rep(0, nRivers * nRivers * myConstants$T * nCohorts) , c(nRivers, nRivers, myConstants$T, nCohorts))
  for(r in 1:nRivers){
    ##for(t in 1:(myConstants$T - 1)){
    for(t in 1:myConstants$T){
      for (c in 1:nCohorts){
        dirch <- rdirch(1, alphaR[[r]])
        for (r2 in 1:nRivers){
          a[r,r2,t,c] <- dirch[r2]
        }
      }
    }
  }
  return(a)
}

mod_initialValues <- function(r, c, nK, t, constants, a){
  list(
    #betaPhiRiver = array(runif(nRivers, 0, 1), c(nRivers)),
    #betaPhiRiverCohort = array(runif(nRivers * nCohorts, 0, 1), c(nRivers, nCohorts)),
    ##phi = array(rnorm(nRivers * (myConstants$T - 1) * nCohorts , 0, 1), c(nRivers, (myConstants$T - 1), nCohorts)),   ## ORIGINAL CODE
    
    phi = array(rnorm(r * nK, 0, 1), c(r, nK)),
    
    betaFlow = array(rnorm(3 * r * 4 , 0, 1), c(3, r, 4)),
    betaInt = array(rnorm(r * 4 * c , 0, 1), c(r, 4, c)),
    
    #betaPRiver = array(runif(nRivers, 0, 1), c(nRivers)),
    #betaPRiverCohort = array(runif(nRivers * nCohorts, 0, 1), c(nRivers, nCohorts)),
    betaP = array(rnorm(r * (t - 1) * c , 0, 1), c(r, (t - 1), c)),        
    psi = getDirchPriorsR(r, constants, c, a)
  )
}


# 
# inits = mod_initialValues(tmp$nRivers, tmp$nCohorts,
#                           tmp$nK, tmp$T,
#                           tmp, alphaR)
# 
# getDirchPriorsR(tmp$nRivers,tmp, tmp$nCohorts, c(1,1,1,1))
# 
# alphaR <- list()

