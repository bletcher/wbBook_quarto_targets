tar_option_set(packages = c("tidyverse", "nimble", "nimbleEcology", "MCMCvis", "targets"))

# This is a 'fork' of modelCMR_tt_ft_cohort_OB, following Daniel's March 2023 model update to the WB (4 river) model

# Run this to load data for testing
#eh_WB_2002_2014_target=tar_read(eh_wB_2002_2014_target)

# in
# constants0
# indexing
# constants

# Model naming
# mod1 = ttt_ft_cohort_WB
# mod2 = ttt_fByRiverT_cohort_WB

#############################################################################
# Add flow
###################################################

modelCMR_ttt_ft_cohort_WB_flow_target <-
  tar_plan(
    
    mod1_input_target =
      list(
        y = eh_WB_2002_2014_target$eh * eh_WB_2002_2014_target$riverN,
        first = eh_WB_2002_2014_target$first, 
        last = eh_WB_2002_2014_target$last, 
        zInits = ttt_initialValues_WB(ncol(eh_WB_2002_2014_target$eh), 
                                      eh_WB_2002_2014_target$eh),
        zInitsNA = ifelse(is.na(eh_WB_2002_2014_target$flow), NA, 1),
        
        flow = eh_WB_2002_2014_target$flow,
        
        temp = eh_WB_2002_2014_target$temperature,
        nCohorts = length(unique(eh_WB_2002_2014_target$cohorts$cohort)),
        cohort = eh_WB_2002_2014_target$cohorts$cohort,
        nSeasons = length(unique(eh_WB_2002_2014_target$data$season)),
        seasonArray = c(3,4,1,2,3,4,1,2,3,4,1,2),
        isYOY = eh_WB_2002_2014_target$isYOY,
        # nStates = length(unique(eh_WB_2002_2014_target$data$sizeState)),
        nRivers = length(unique(eh_WB_2002_2014_target$data$riverN))
      ), 
    
    mod1_runData_target = list(
      # Updateable model-specific variables 
      nIter = 30, 
      nBurnin = 0, 
      nChains = 1,
      thinRate = 0
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
      #alpha = mod1_myConstants_target$alpha,
      #Priors for psi where more likely to stay than move
      alphaR1 = mod1_myConstants_target$alphaR1,
      alphaR2 = mod1_myConstants_target$alphaR2,
      alphaR3 = mod1_myConstants_target$alphaR3,
      alphaR4 = mod1_myConstants_target$alphaR4,
      
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
      season = NULL
    ),
    
    #####################################
    # 
    mod1_myData_target = list(
      y = mod1_input_target$y[mod1_myConstants0_target$indToKeep,] + 1
    ),
    # 
    mod1_modelCode_target = nimbleCode({
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
    
    mod1_Rmodel = nimbleModel(
      code = mod1_modelCode_target,
      constants = mod1_constantsIn_target,
      data = mod1_myData_target,
      inits = mod1_initialValues(mod1_constantsIn_target$nRivers, mod1_constantsIn_target$nCohorts,
                                                 mod1_constantsIn_target$nK, mod1_constantsIn_target$T,
                                                 mod1_constantsIn_target, mod1_constantsIn_target$alphaR),
      calculate = FALSE
    )
  )
#     
#     tt_parametersToSave_ft_cohort_WB_flow = c(
#       "betaIntTop", "betaFlowTop","betaPTop",  
#       "betaIntYOY", "betaFlowYOY","betaPYOY",
#       "betaIntYOYSeason", "betaPYOYSeason",
#       "betaInt", "betaFlow", "betaP",
#                                     
#       "betaIntYOYOut", "betaPYOYOut",
#       "betaIntYOYSeasonOut", "betaPYOYSeasonOut",
#       "betaIntOut", "betaPOut"
#     ),
#     
#     tt_conf_ft_cohort_WB_flow = configureMCMC(
#       tt_Rmodel_ft_cohort_WB_flow,
#       monitors = tt_parametersToSave_ft_cohort_WB_flow
#     ),
#     
#     tt_Rmcmc_ft_cohort_WB_flow = buildMCMC(tt_conf_ft_cohort_WB_flow, useConjugacy = FALSE),
#     tt_Cmodel_ft_cohort_WB_flow = compileNimble(tt_Rmodel_ft_cohort_WB_flow),
#     tt_Cmcmc_ft_cohort_WB_flow = compileNimble(tt_Rmcmc_ft_cohort_WB_flow, project = tt_Rmodel_ft_cohort_WB_flow),
#     
#     tt_model_ft_cohort_WB_flow = runMCMC(
#       tt_Cmcmc_ft_cohort_WB_flow,
#       niter = tt_runData_ft_cohort_WB_flow$nIter,
#       nburnin = tt_runData_ft_cohort_WB_flow$nBurnin,
#       thin = tt_runData_ft_cohort_WB_flow$thinRate,
#       nchains = tt_runData_ft_cohort_WB_flow$nChains
#     ),
#     
#     tt_modelOut_ft_cohort_WB_flow =
#       list(
#         mcmc = tt_model_ft_cohort_WB_flow, # "Error : invalid nimbleFunction argument"
#         name = "phiT_pT_ft_cohort_WB_flow",
#         modelCode = tt_modelCode_ft_cohort_WB_flow,
#         myConstants = tt_myConstants_ft_cohort_WB_flow,
#         runData = tt_runData_ft_cohort_WB_flow
#       ),
#     
#     tt_save_ft_cohort_WB_flow = saveModelOut_tt_ft_cohort_WB_flow(tt_modelOut_ft_cohort_WB_flow)
#   )
# 
# 
# #############################################################################
# # Add flow by river
# # Flow estimated for each river independently
# ###################################################
# modelCMR_tt_ft_cohort_WB_flowByRiver_target <-
#   tar_plan(
#     
#     tt_inputData_ft_cohort_WB_flowByRiver =
#       list(
#         y = eh_WB_2002_2014_target$eh,
#         first = eh_WB_2002_2014_target$first, 
#         last = eh_WB_2002_2014_target$last, 
#         zInits = tt_initialValues_OB(ncol(eh_WB_2002_2014_target$eh), 
#                                      eh_WB_2002_2014_target$eh),
#         zInitsNA = ifelse(is.na(eh_WB_2002_2014_target$flow), NA, 1),
#         #############
#         flow = eh_WB_2002_2014_target$flowByRiver,
#         ##################
#         temp = eh_WB_2002_2014_target$temperature,
#         nCohorts = length(unique(eh_WB_2002_2014_target$cohorts$cohort)),
#         cohort = eh_WB_2002_2014_target$cohorts$cohort,
#         nSeasons = length(unique(eh_WB_2002_2014_target$data$season)),
#         seasonArray = c(3,4,1,2,3,4,1,2,3,4,1,2),
#         isYOY = eh_WB_2002_2014_target$isYOY
#         # nStates = length(unique(eh_WB_2002_2014_target$data$sizeState)),
#         # nRivers = length(unique(eh_WB_2002_2014_target$data$sizeState)) # for now
#       ), 
#     
#     tt_runData_ft_cohort_WB_flowByRiver = list(
#       # Updateable model-specific variables 
#       nIter = 30000, 
#       nBurnin = 20000, 
#       nChains = 2,
#       thinRate = 5
#     ),  
#     
#     tt_myConstants_ft_cohort_WB_flowByRiver0 = list(
#       N = nrow(tt_inputData_ft_cohort_WB_flowByRiver$y),
#       T = ncol(tt_inputData_ft_cohort_WB_flowByRiver$y),
#       first = tt_inputData_ft_cohort_WB_flowByRiver$first,
#       last = tt_inputData_ft_cohort_WB_flowByRiver$last,
#       cohort = tt_inputData_ft_cohort_WB_flowByRiver$cohort - min(tt_inputData_ft_cohort_WB_flowByRiver$cohort) + 1,
#       nCohorts = tt_inputData_ft_cohort_WB_flowByRiver$nCohorts,
#       season = tt_inputData_ft_cohort_WB_flowByRiver$seasonArray,
#       flow = tt_inputData_ft_cohort_WB_flowByRiver$flow,
#       temp = tt_inputData_ft_cohort_WB_flowByRiver$temp,
#       length = tt_inputData_ft_cohort_WB_flowByRiver$last - tt_inputData_ft_cohort_WB_flowByRiver$first + 1,
#       isYOY = tt_inputData_ft_cohort_WB_flowByRiver$isYOY,
#       indToKeep = which(tt_inputData_ft_cohort_WB_flowByRiver$first < 12)
#     ),
#     
#     tt_myConstants_ft_cohort_WB_flowByRiver = list(
#       N = length(tt_myConstants_ft_cohort_WB_flowByRiver0$indToKeep),
#       T = ncol(tt_inputData_ft_cohort_WB_flowByRiver$y),
#       first = tt_myConstants_ft_cohort_WB_flowByRiver0$first[tt_myConstants_ft_cohort_WB_flowByRiver0$indToKeep],
#       last = tt_myConstants_ft_cohort_WB_flowByRiver0$last[tt_myConstants_ft_cohort_WB_flowByRiver0$indToKeep],
#       cohort = tt_myConstants_ft_cohort_WB_flowByRiver0$cohort[tt_myConstants_ft_cohort_WB_flowByRiver0$indToKeep], 
#       season = tt_myConstants_ft_cohort_WB_flowByRiver0$season,
#       flow = tt_myConstants_ft_cohort_WB_flowByRiver0$flow[tt_myConstants_ft_cohort_WB_flowByRiver0$indToKeep,],
#       temp = tt_myConstants_ft_cohort_WB_flowByRiver0$temp[tt_myConstants_ft_cohort_WB_flowByRiver0$indToKeep,],
#       length = tt_myConstants_ft_cohort_WB_flowByRiver0$length[tt_myConstants_ft_cohort_WB_flowByRiver0$indToKeep],
#       nCohorts = tt_inputData_ft_cohort_WB_flowByRiver$nCohorts,
#       nSeasons = tt_inputData_ft_cohort_WB_flowByRiver$nSeasons,
#       isYOY = tt_inputData_ft_cohort_WB_flowByRiver$isYOY[tt_myConstants_ft_cohort_WB_flowByRiver0$indToKeep,]
#     ),
#     
#     tt_myData_ft_cohort_WB_flowByRiver = list(
#       y = tt_inputData_ft_cohort_WB_flowByRiver$y + 1
#     ),
#     
#     tt_modelCode_ft_cohort_WB_flowByRiver = nimbleCode({
#       dummy  ~ dnorm(0,1)
#       
#       # from https://bletcher.github.io/westBrook-book/models.html#model-phit_pt_cohort_flowcohorthierdhmm
#       
#       delta[1] <- 1                    # Pr(alive t = 1) = 1
#       delta[2] <- 0                    # Pr(dead t = 1) = 0
#       
# 
#       for (i in 1:N){
#         for (t in 1:(T-1)){ # loop over time
#           logit(phi[t,i]) <- 
#             betaInt[   isYOY[i,t],season[t],cohort[i]] +
#          #   betaPhi[   isYOY[i,t],season[t],cohort[i]] + 
#             betaFlow[1,isYOY[i,t],season[t]] * flow[i,t] +
#             betaFlow[2,isYOY[i,t],season[t]] * temp[i,t] +
#             betaFlow[3,isYOY[i,t],season[t]] * temp[i,t] * flow[i,t]
#           # prior survival
#           ##
#           gamma[1,1,t,i] <- phi[t,i]         # Pr(alive t -> alive t+1)
#           gamma[1,2,t,i] <- 1 - phi[t,i]     # Pr(alive t -> dead t+1)
#           gamma[2,1,t,i] <- 0              # Pr(dead t -> alive t+1)
#           gamma[2,2,t,i] <- 1              # Pr(dead t -> dead t+1)
#         }
#         
#         gamma[1,1,T,i] <- 0
#         gamma[1,2,T,i] <- 1
#         gamma[2,1,T,i] <- 0
#         gamma[2,2,T,i] <- 1
#         
#         omega[1,1,first[i],i] <- 0       # Pr(alive t -> non-detected t)
#         omega[1,2,first[i],i] <- 1           # Pr(alive t -> detected t)
#         omega[2,1,first[i],i] <- 1              # Pr(dead t -> non-detected t)
#         omega[2,2,first[i],i] <- 0              # Pr(dead t -> detected t)
#         
#         for(t in (first[i]+1):last[i]) {
#           logit(p[t,i]) <- betaP[isYOY[i,t],season[t-1],cohort[i]]             # prior detection
#           omega[1,1,t,i] <- 1 - p[t,i]       # Pr(alive t -> non-detected t)
#           omega[1,2,t,i] <- p[t,i]           # Pr(alive t -> detected t)
#           omega[2,1,t,i] <- 1              # Pr(dead t -> non-detected t)
#           omega[2,2,t,i] <- 0              # Pr(dead t -> detected t)
#         }
#       }
#       
#       ## 
#       # for (y in 1:2) {
#       #   for (s in 1:nSeasons){
#       #     for (c in 1:nCohorts){
#       #       betaInt[y,s,c] ~ dnorm(betaIntYOY[y],1)
#       #       betaPhi[y,s,c] ~ dnorm(betaPhiYOY[y],1)
#       #       betaP[y,s,c] ~ dnorm(betaPYOY[y],1)
#       #     }
#       #     
#       #     betaFlow[1,y,s] ~ dnorm(betaFlowYOY[1,y],1)
#       #     betaFlow[2,y,s] ~ dnorm(betaFlowYOY[2,y],1)
#       #     betaFlow[3,y,s] ~ dnorm(betaFlowYOY[3,y],1)
#       #     
#       #   }
#       # }
#       
#       for (y in 1:2) {
#         for (s in 1:nSeasons){
#           for (c in 1:nCohorts){
#             betaInt[y,s,c] ~ dnorm(betaIntYOYSeason[y,s],1)
#             #betaPhi[y,s,c] ~ dnorm(betaPhiYOYSeason[y,s],1)
#             betaP[y,s,c] ~ dnorm(betaPYOYSeason[y,s],1)
#           }
#           
#           betaIntYOYSeason[y,s] ~ dnorm(betaIntYOY[y],1)
#           #betaPhiYOYSeason[y,s] ~ dnorm(betaPhiYOY[y],1)
#           betaPYOYSeason[y,s] ~ dnorm(betaPYOY[y],1)
#           
#           betaFlow[1,y,s] ~ dnorm(betaFlowYOY[1,y],1)
#           betaFlow[2,y,s] ~ dnorm(betaFlowYOY[2,y],1)
#           betaFlow[3,y,s] ~ dnorm(betaFlowYOY[3,y],1)
#           
#         }
#       }
#       
#       for (y in 1:2) {
#         
#         betaIntYOY[y] ~ dnorm(betaIntTop,1)
#         #betaPhiYOY[y] ~ dnorm(betaPhiTop,1)
#         betaFlowYOY[1,y] ~ dnorm(betaFlowTop[1],1)
#         betaFlowYOY[2,y] ~ dnorm(betaFlowTop[2],1)
#         betaFlowYOY[3,y] ~ dnorm(betaFlowTop[3],1)
#         
#         betaPYOY[y] ~ dnorm(0,1)
#       }
#       
#       betaIntTop ~ dnorm(0,1)
#       #betaPhiTop ~ dnorm(0,1)
#       betaFlowTop[1] ~ dnorm(0,1)
#       betaFlowTop[2] ~ dnorm(0,1)
#       betaFlowTop[3] ~ dnorm(0,1)
#       
#       betaPTop ~ dnorm(0,1)
#       
#       ##    
#       # back-transform for examining output
#       ##
#       for (y in 1:2) {
#         betaIntYOYOut[y] <- 1/(1 + exp(-betaIntYOY[y]))
#         betaPYOYOut[y] <- 1/(1 + exp(-betaPYOY[y]))
#         
#         for (s in 1:nSeasons) {
#           betaIntYOYSeasonOut[y,s] <- 1/(1 + exp(-betaIntYOYSeason[y,s]))
#           betaPYOYSeasonOut[y,s] <- 1/(1 + exp(-betaPYOYSeason[y,s]))
#           
#           for (c in 1:nCohorts) {
#             betaIntOut[y,s,c] <- 1/(1 + exp(-betaInt[y,s,c]))
#             betaPOut[y,s,c] <- 1/(1 + exp(-betaP[y,s,c]))
#           }
#         }
#       }
#       ##    
#       # likelihood
#       for (i in 1:N){
#         y[i,first[i]:last[i]] ~ dDHMMo(init = delta[1:2],
#                                        probTrans = gamma[1:2, 1:2, first[i]:last[i], i],
#                                        probObs =   omega[1:2, 1:2, first[i]:last[i], i],
#                                        len = length[i],
#                                        checkRowSums = 1)
#       }
#     }),
#     
#     tt_Rmodel_ft_cohort_WB_flowByRiver = nimbleModel(
#       code = tt_modelCode_ft_cohort_WB_flowByRiver,
#       constants = tt_myConstants_ft_cohort_WB_flowByRiver,
#       data = tt_myData_ft_cohort_WB_flowByRiver,
#       inits = tt_initialValues_ft_cohort_WB_flow(tt_myConstants_ft_cohort_WB_flowByRiver$T, tt_myConstants_ft_cohort_WB_flowByRiver$nCohorts, 
#                                           tt_inputData_ft_cohort_WB_flowByRiver$y, tt_inputData_ft_cohort_WB_flowByRiver$zInitsNA),
#       calculate = FALSE
#     ),
#     
#     tt_parametersToSave_ft_cohort_WB_flowByRiver = c(
#       "betaIntTop", "betaFlowTop","betaPTop",  
#       "betaIntYOY", "betaFlowYOY","betaPYOY",
#       "betaIntYOYSeason", "betaPYOYSeason",
#       "betaInt", "betaFlow", "betaP",
#       
#       "betaIntYOYOut", "betaPYOYOut",
#       "betaIntYOYSeasonOut", "betaPYOYSeasonOut",
#       "betaIntOut", "betaPOut"
#     ),
#     
#     tt_conf_ft_cohort_WB_flowByRiver = configureMCMC(
#       tt_Rmodel_ft_cohort_WB_flowByRiver,
#       monitors = tt_parametersToSave_ft_cohort_WB_flowByRiver
#     ),
#     
#     tt_Rmcmc_ft_cohort_WB_flowByRiver = buildMCMC(tt_conf_ft_cohort_WB_flowByRiver, useConjugacy = FALSE),
#     tt_Cmodel_ft_cohort_WB_flowByRiver = compileNimble(tt_Rmodel_ft_cohort_WB_flowByRiver),
#     tt_Cmcmc_ft_cohort_WB_flowByRiver = compileNimble(tt_Rmcmc_ft_cohort_WB_flowByRiver, project = tt_Rmodel_ft_cohort_WB_flowByRiver),
#     
#     tt_model_ft_cohort_WB_flowByRiver = runMCMC(
#       tt_Cmcmc_ft_cohort_WB_flowByRiver,
#       niter = tt_runData_ft_cohort_WB_flowByRiver$nIter,
#       nburnin = tt_runData_ft_cohort_WB_flowByRiver$nBurnin,
#       thin = tt_runData_ft_cohort_WB_flowByRiver$thinRate,
#       nchains = tt_runData_ft_cohort_WB_flowByRiver$nChains
#     ),
#     
#     tt_modelOut_ft_cohort_WB_flowByRiver =
#       list(
#         mcmc = tt_model_ft_cohort_WB_flowByRiver, # "Error : invalid nimbleFunction argument"
#         name = "phiT_pT_ft_cohort_WB_flowByRiver",
#         modelCode = tt_modelCode_ft_cohort_WB_flowByRiver,
#         myConstants = tt_myConstants_ft_cohort_WB_flowByRiver,
#         runData = tt_runData_ft_cohort_WB_flowByRiver
#       ),
#     
#     tt_save_ft_cohort_WB_flowByRiver = saveModelOut_tt_ft_cohort_WB_flowByRiver(tt_modelOut_ft_cohort_WB_flowByRiver)
#   )
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
saveModelOut_tt_ft_cohort_WB_flow <- function(d) {
  save(d, file = paste0('./models/cmrFlowOB/runsOut/tt_ft_cohort_WB_flow', substr(Sys.time(),1,13), '.RData'))
}

saveModelOut_tt_ft_cohort_WB_flowByRiver <- function(d) {
  save(d, file = paste0('./models/cmrFlowOB/runsOut/tt_ft_cohort_WB_flowByRiver', substr(Sys.time(),1,13), '.RData'))
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

getDirchPriorsR <- function(nRivers,myConstants, nCohorts, alphaR){
  ##a = array(rep(0, nRivers * nRivers * (myConstants$T - 1) * nCohorts) , c(nRivers, nRivers, (myConstants$T - 1), nCohorts ))
  a = array(rep(0, nRivers * nRivers * myConstants$T * nCohorts) , c(nRivers, nRivers, myConstants$T, nCohorts ))
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

mod1_initialValues <- function(r, c, nK, t, constants, a){
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



