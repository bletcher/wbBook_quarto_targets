tar_option_set(packages = c("tidyverse", "nimble", "nimbleEcology", "MCMCvis", "targets"))

#load('./models/cmrFlowOB/dataOut/eh_2002200320042005200620072008200920102011201220132014_wb obear.RData')

#eh_OB_target = tar_read(eh_OB_2002_2014_target)
#eh_OB_2002_2014_target=tar_read(eh_OB_2002_2014_target)


#############################################################################
# Add flow
###################################################

modelCMR_tt_ft_cohort_OB_flow_target <-
  tar_plan(
    
    tt_inputData_ft_cohort_OB_flow =
      list(
        y = eh_OB_2002_2014_target$eh,
        first = eh_OB_2002_2014_target$first, 
        last = eh_OB_2002_2014_target$last, 
        zInits = tt_initialValues_OB(ncol(eh_OB_2002_2014_target$eh), 
                                             eh_OB_2002_2014_target$eh),
        zInitsNA = ifelse(is.na(eh_OB_2002_2014_target$flow), NA, 1),
        flow = eh_OB_2002_2014_target$flow,
        temp = eh_OB_2002_2014_target$temperature,
        nCohorts = length(unique(eh_OB_2002_2014_target$cohorts$cohort)),
        cohort = eh_OB_2002_2014_target$cohorts$cohort,
        nSeasons = length(unique(eh_OB_2002_2014_target$data$season)),
        seasonArray = c(3,4,1,2,3,4,1,2,3,4,1,2),
        isYOY = eh_OB_2002_2014_target$isYOY
        # nStates = length(unique(eh_OB_2002_2014_target$data$sizeState)),
        # nRivers = length(unique(eh_OB_2002_2014_target$data$sizeState)) # for now
      ), 
    
    tt_runData_ft_cohort_OB_flow = list(
      # Updateable model-specific variables 
      nIter = 30000, 
      nBurnin = 20000, 
      nChains = 2,
      thinRate = 5
    ),  
    
    tt_myConstants_ft_cohort_OB_flow0 = list(
      N = nrow(tt_inputData_ft_cohort_OB_flow$y),
      T = ncol(tt_inputData_ft_cohort_OB_flow$y),
      first = tt_inputData_ft_cohort_OB_flow$first,
      last = tt_inputData_ft_cohort_OB_flow$last,
      cohort = tt_inputData_ft_cohort_OB_flow$cohort - min(tt_inputData_ft_cohort_OB_flow$cohort) + 1,
      nCohorts = tt_inputData_ft_cohort_OB_flow$nCohorts,
      season = tt_inputData_ft_cohort_OB_flow$seasonArray,
      flow = tt_inputData_ft_cohort_OB_flow$flow,
      temp = tt_inputData_ft_cohort_OB_flow$temp,
      length = tt_inputData_ft_cohort_OB_flow$last - tt_inputData_ft_cohort_OB_flow$first + 1,
      isYOY = tt_inputData_ft_cohort_OB_flow$isYOY,
      indToKeep = which(tt_inputData_ft_cohort_OB_flow$first < 12)
    ),
    
    tt_myConstants_ft_cohort_OB_flow = list(
      N = length(tt_myConstants_ft_cohort_OB_flow0$indToKeep),
      T = ncol(tt_inputData_ft_cohort_OB_flow$y),
      first = tt_myConstants_ft_cohort_OB_flow0$first[tt_myConstants_ft_cohort_OB_flow0$indToKeep],
      last = tt_myConstants_ft_cohort_OB_flow0$last[tt_myConstants_ft_cohort_OB_flow0$indToKeep],
      cohort = tt_myConstants_ft_cohort_OB_flow0$cohort[tt_myConstants_ft_cohort_OB_flow0$indToKeep], 
      season = tt_myConstants_ft_cohort_OB_flow0$season,
      flow = tt_myConstants_ft_cohort_OB_flow0$flow[tt_myConstants_ft_cohort_OB_flow0$indToKeep,],
      temp = tt_myConstants_ft_cohort_OB_flow0$temp[tt_myConstants_ft_cohort_OB_flow0$indToKeep,],
      length = tt_myConstants_ft_cohort_OB_flow0$length[tt_myConstants_ft_cohort_OB_flow0$indToKeep],
      nCohorts = tt_inputData_ft_cohort_OB_flow$nCohorts,
      nSeasons = tt_inputData_ft_cohort_OB_flow$nSeasons,
      isYOY = tt_inputData_ft_cohort_OB_flow$isYOY[tt_myConstants_ft_cohort_OB_flow0$indToKeep,]
    ),
    
    tt_myData_ft_cohort_OB_flow = list(
      y = tt_inputData_ft_cohort_OB_flow$y + 1
    ),
    
    tt_modelCode_ft_cohort_OB_flow = nimbleCode({
      #dummy  ~ dnorm(0,1)
      
      # from https://bletcher.github.io/westBrook-book/models.html#model-phit_pt_cohort_flowcohorthierdhmm
      
      delta[1] <- 1                    # Pr(alive t = 1) = 1
      delta[2] <- 0                    # Pr(dead t = 1) = 0
      
      #log flow data?

      for (i in 1:N){
        for (t in 1:(T-1)){ # loop over time
          logit(phi[t,i]) <- 
            betaInt[   isYOY[i,t],season[t],cohort[i]] +
        #    betaPhi[   isYOY[i,t],season[t],cohort[i]] + 
            betaFlow[1,isYOY[i,t],season[t]] * flow[i,t] +
            betaFlow[2,isYOY[i,t],season[t]] * temp[i,t] +
            betaFlow[3,isYOY[i,t],season[t]] * temp[i,t] * flow[i,t]
          # prior survival
          ##
          gamma[1,1,t,i] <- phi[t,i]         # Pr(alive t -> alive t+1)
          gamma[1,2,t,i] <- 1 - phi[t,i]     # Pr(alive t -> dead t+1)
          gamma[2,1,t,i] <- 0              # Pr(dead t -> alive t+1)
          gamma[2,2,t,i] <- 1              # Pr(dead t -> dead t+1)
        }

        gamma[1,1,T,i] <- 0
        gamma[1,2,T,i] <- 1
        gamma[2,1,T,i] <- 0
        gamma[2,2,T,i] <- 1

        omega[1,1,first[i],i] <- 0       # Pr(alive t -> non-detected t)
        omega[1,2,first[i],i] <- 1           # Pr(alive t -> detected t)
        omega[2,1,first[i],i] <- 1              # Pr(dead t -> non-detected t)
        omega[2,2,first[i],i] <- 0              # Pr(dead t -> detected t)

        for(t in (first[i]+1):last[i]) {
          logit(p[t,i]) <- betaP[isYOY[i,t],season[t-1],cohort[i]]             # prior detection # CHECK t, t-1 here
          omega[1,1,t,i] <- 1 - p[t,i]       # Pr(alive t -> non-detected t)
          omega[1,2,t,i] <- p[t,i]           # Pr(alive t -> detected t)
          omega[2,1,t,i] <- 1              # Pr(dead t -> non-detected t)
          omega[2,2,t,i] <- 0              # Pr(dead t -> detected t)
        }
      }
      ## 
    
      for (y in 1:2) {
        for (s in 1:nSeasons){
          for (c in 1:nCohorts){
            betaInt[y,s,c] ~ dnorm(betaIntYOYSeason[y,s],1)
            betaPhi[y,s,c] ~ dnorm(betaPhiYOYSeason[y,s],1)
            betaP[y,s,c] ~ dnorm(betaPYOYSeason[y,s],1)
          }
          
          betaIntYOYSeason[y,s] ~ dnorm(betaIntYOY[y],1)
          betaPhiYOYSeason[y,s] ~ dnorm(betaPhiYOY[y],1)
          betaPYOYSeason[y,s] ~ dnorm(betaPYOY[y],1)
          
          betaFlow[1,y,s] ~ dnorm(betaFlowYOY[1,y],1)
          betaFlow[2,y,s] ~ dnorm(betaFlowYOY[2,y],1)
          betaFlow[3,y,s] ~ dnorm(betaFlowYOY[3,y],1)
          
        }
      }
          
      for (y in 1:2) {
        
        betaIntYOY[y] ~ dnorm(betaIntTop,1)
        betaPhiYOY[y] ~ dnorm(betaPhiTop,1)
        betaFlowYOY[1,y] ~ dnorm(betaFlowTop[1],1)
        betaFlowYOY[2,y] ~ dnorm(betaFlowTop[2],1)
        betaFlowYOY[3,y] ~ dnorm(betaFlowTop[3],1)
        
        betaPYOY[y] ~ dnorm(0,1)
      }
      
      betaIntTop ~ dnorm(0,1)
      betaPhiTop ~ dnorm(0,1)
      betaFlowTop[1] ~ dnorm(0,1)
      betaFlowTop[2] ~ dnorm(0,1)
      betaFlowTop[3] ~ dnorm(0,1)
      
      betaPTop ~ dnorm(0,1)
            
      ##    
      # back-transform for examining output
      ##
      for (y in 1:2) {
        for (s in 1:nSeasons){
         # for (c in 1:nCohorts){
            betaIntYOYSeasonOut[y,s] <- 1/(1 + exp(-betaIntYOYSeason[y,s]))
            betaPhiYOYSeasonOut[y,s] <- 1/(1 + exp(-betaPhiYOYSeason[y,s]))
          #  betaFlowOut[1,y,s,c] <- 1/(1 + exp(-betaFlow[1,y,s,c]))
          #  betaFlowOut[2,y,s,c] <- 1/(1 + exp(-betaFlow[2,y,s,c]))
          #  betaFlowOut[3,y,s,c] <- 1/(1 + exp(-betaFlow[3,y,s,c]))

            betaPYOYSeasonOut[y,s] <- 1/(1 + exp(-betaPYOYSeason[y,s]))
        #  }
        }
      }
      ##    
      # likelihood
      for (i in 1:N){
        y[i,first[i]:last[i]] ~ dDHMMo(init = delta[1:2],
                                       probTrans = gamma[1:2, 1:2, first[i]:last[i], i],
                                       probObs =   omega[1:2, 1:2, first[i]:last[i], i],
                                       len = length[i],
                                       checkRowSums = 1)
      }
    }),
    
    tt_Rmodel_ft_cohort_OB_flow = nimbleModel(
      code = tt_modelCode_ft_cohort_OB_flow,
      constants = tt_myConstants_ft_cohort_OB_flow,
      data = tt_myData_ft_cohort_OB_flow,
      inits = tt_initialValues_ft_cohort_OB_flow(tt_myConstants_ft_cohort_OB_flow$T, tt_myConstants_ft_cohort_OB_flow$nCohorts, 
                                                 tt_inputData_ft_cohort_OB_flow$y, tt_inputData_ft_cohort_OB_flow$zInitsNA),
      calculate = FALSE
    ),
    
    tt_parametersToSave_ft_cohort_OB_flow = c("betaIntTop", "betaPhiTop","betaFlowTop","betaPTop",  
                                    "betaIntYOY", "betaPhiYOY","betaFlowYOY","betaPYOY",
                                    "betaIntYOYSeason", "betaPhiYOYSeason","betaPYOYSeason",
                                    "betaIntYOYSeasonOut", "betaPhiYOYSeasonOut","betaPYOYSeasonOut",
                                    "betaInt", "betaPhi","betaFlow","betaP"
                                   # "betaIntOut", "betaPhiOut","betaFlowOut","betaPOut"
                                    ),
    
    tt_conf_ft_cohort_OB_flow = configureMCMC(
      tt_Rmodel_ft_cohort_OB_flow,
      monitors = tt_parametersToSave_ft_cohort_OB_flow
    ),
    
    tt_Rmcmc_ft_cohort_OB_flow = buildMCMC(tt_conf_ft_cohort_OB_flow, useConjugacy = FALSE),
    tt_Cmodel_ft_cohort_OB_flow = compileNimble(tt_Rmodel_ft_cohort_OB_flow),
    tt_Cmcmc_ft_cohort_OB_flow = compileNimble(tt_Rmcmc_ft_cohort_OB_flow, project = tt_Rmodel_ft_cohort_OB_flow),
    
    tt_model_ft_cohort_OB_flow = runMCMC(
      tt_Cmcmc_ft_cohort_OB_flow,
      niter = tt_runData_ft_cohort_OB_flow$nIter,
      nburnin = tt_runData_ft_cohort_OB_flow$nBurnin,
      thin = tt_runData_ft_cohort_OB_flow$thinRate,
      nchains = tt_runData_ft_cohort_OB_flow$nChains
    ),
    
    tt_modelOut_ft_cohort_OB_flow =
      list(
        mcmc = tt_model_ft_cohort_OB_flow, # "Error : invalid nimbleFunction argument"
        name = "phiT_pT_ft_cohort_OB_flow",
        modelCode = tt_modelCode_ft_cohort_OB_flow,
        myConstants = tt_myConstants_ft_cohort_OB_flow,
        runData = tt_runData_ft_cohort_OB_flow
      ),
    
    tt_save_ft_cohort_OB_flow = saveModelOut_tt_ft_cohort_OB_flow(tt_modelOut_ft_cohort_OB_flow)
  )


#############################################################################
# Add flow by river
# Flow estimated for each river independently
###################################################
modelCMR_tt_ft_cohort_OB_flowByRiver_target <-
  tar_plan(
    
    tt_inputData_ft_cohort_OB_flowByRiver =
      list(
        y = eh_OB_2002_2014_target$eh,
        first = eh_OB_2002_2014_target$first, 
        last = eh_OB_2002_2014_target$last, 
        zInits = tt_initialValues_OB(ncol(eh_OB_2002_2014_target$eh), 
                                     eh_OB_2002_2014_target$eh),
        zInitsNA = ifelse(is.na(eh_OB_2002_2014_target$flow), NA, 1),
        #############
        flow = eh_OB_2002_2014_target$flowByRiver,
        ##################
        temp = eh_OB_2002_2014_target$temperature,
        nCohorts = length(unique(eh_OB_2002_2014_target$cohorts$cohort)),
        cohort = eh_OB_2002_2014_target$cohorts$cohort,
        nSeasons = length(unique(eh_OB_2002_2014_target$data$season)),
        seasonArray = c(3,4,1,2,3,4,1,2,3,4,1,2),
        isYOY = eh_OB_2002_2014_target$isYOY
        # nStates = length(unique(eh_OB_2002_2014_target$data$sizeState)),
        # nRivers = length(unique(eh_OB_2002_2014_target$data$sizeState)) # for now
      ), 
    
    tt_runData_ft_cohort_OB_flowByRiver = list(
      # Updateable model-specific variables 
      nIter = 30000, 
      nBurnin = 20000, 
      nChains = 2,
      thinRate = 5
    ),  
    
    tt_myConstants_ft_cohort_OB_flowByRiver0 = list(
      N = nrow(tt_inputData_ft_cohort_OB_flowByRiver$y),
      T = ncol(tt_inputData_ft_cohort_OB_flowByRiver$y),
      first = tt_inputData_ft_cohort_OB_flowByRiver$first,
      last = tt_inputData_ft_cohort_OB_flowByRiver$last,
      cohort = tt_inputData_ft_cohort_OB_flowByRiver$cohort - min(tt_inputData_ft_cohort_OB_flowByRiver$cohort) + 1,
      nCohorts = tt_inputData_ft_cohort_OB_flowByRiver$nCohorts,
      season = tt_inputData_ft_cohort_OB_flowByRiver$seasonArray,
      flow = tt_inputData_ft_cohort_OB_flowByRiver$flow,
      temp = tt_inputData_ft_cohort_OB_flowByRiver$temp,
      length = tt_inputData_ft_cohort_OB_flowByRiver$last - tt_inputData_ft_cohort_OB_flowByRiver$first + 1,
      isYOY = tt_inputData_ft_cohort_OB_flowByRiver$isYOY,
      indToKeep = which(tt_inputData_ft_cohort_OB_flowByRiver$first < 12)
    ),
    
    tt_myConstants_ft_cohort_OB_flowByRiver = list(
      N = length(tt_myConstants_ft_cohort_OB_flowByRiver0$indToKeep),
      T = ncol(tt_inputData_ft_cohort_OB_flowByRiver$y),
      first = tt_myConstants_ft_cohort_OB_flowByRiver0$first[tt_myConstants_ft_cohort_OB_flowByRiver0$indToKeep],
      last = tt_myConstants_ft_cohort_OB_flowByRiver0$last[tt_myConstants_ft_cohort_OB_flowByRiver0$indToKeep],
      cohort = tt_myConstants_ft_cohort_OB_flowByRiver0$cohort[tt_myConstants_ft_cohort_OB_flowByRiver0$indToKeep], 
      season = tt_myConstants_ft_cohort_OB_flowByRiver0$season,
      flow = tt_myConstants_ft_cohort_OB_flowByRiver0$flow[tt_myConstants_ft_cohort_OB_flowByRiver0$indToKeep,],
      temp = tt_myConstants_ft_cohort_OB_flowByRiver0$temp[tt_myConstants_ft_cohort_OB_flowByRiver0$indToKeep,],
      length = tt_myConstants_ft_cohort_OB_flowByRiver0$length[tt_myConstants_ft_cohort_OB_flowByRiver0$indToKeep],
      nCohorts = tt_inputData_ft_cohort_OB_flowByRiver$nCohorts,
      nSeasons = tt_inputData_ft_cohort_OB_flowByRiver$nSeasons,
      isYOY = tt_inputData_ft_cohort_OB_flowByRiver$isYOY[tt_myConstants_ft_cohort_OB_flowByRiver0$indToKeep,]
    ),
    
    tt_myData_ft_cohort_OB_flowByRiver = list(
      y = tt_inputData_ft_cohort_OB_flowByRiver$y + 1
    ),
    
    tt_modelCode_ft_cohort_OB_flowByRiver = nimbleCode({
      #dummy  ~ dnorm(0,1)
      
      # from https://bletcher.github.io/westBrook-book/models.html#model-phit_pt_cohort_flowcohorthierdhmm
      
      delta[1] <- 1                    # Pr(alive t = 1) = 1
      delta[2] <- 0                    # Pr(dead t = 1) = 0
      

      for (i in 1:N){
        for (t in 1:(T-1)){ # loop over time
          logit(phi[t,i]) <- 
            betaInt[   isYOY[i,t],season[t],cohort[i]] +
         #   betaPhi[   isYOY[i,t],season[t],cohort[i]] + 
            betaFlow[1,isYOY[i,t],season[t]] * flow[i,t] +
            betaFlow[2,isYOY[i,t],season[t]] * temp[i,t] +
            betaFlow[3,isYOY[i,t],season[t]] * temp[i,t] * flow[i,t]
          # prior survival
          ##
          gamma[1,1,t,i] <- phi[t,i]         # Pr(alive t -> alive t+1)
          gamma[1,2,t,i] <- 1 - phi[t,i]     # Pr(alive t -> dead t+1)
          gamma[2,1,t,i] <- 0              # Pr(dead t -> alive t+1)
          gamma[2,2,t,i] <- 1              # Pr(dead t -> dead t+1)
        }
        
        gamma[1,1,T,i] <- 0
        gamma[1,2,T,i] <- 1
        gamma[2,1,T,i] <- 0
        gamma[2,2,T,i] <- 1
        
        omega[1,1,first[i],i] <- 0       # Pr(alive t -> non-detected t)
        omega[1,2,first[i],i] <- 1           # Pr(alive t -> detected t)
        omega[2,1,first[i],i] <- 1              # Pr(dead t -> non-detected t)
        omega[2,2,first[i],i] <- 0              # Pr(dead t -> detected t)
        
        for(t in (first[i]+1):last[i]) {
          logit(p[t,i]) <- betaP[isYOY[i,t],season[t-1],cohort[i]]             # prior detection
          omega[1,1,t,i] <- 1 - p[t,i]       # Pr(alive t -> non-detected t)
          omega[1,2,t,i] <- p[t,i]           # Pr(alive t -> detected t)
          omega[2,1,t,i] <- 1              # Pr(dead t -> non-detected t)
          omega[2,2,t,i] <- 0              # Pr(dead t -> detected t)
        }
      }
      
      ## 
      # for (y in 1:2) {
      #   for (s in 1:nSeasons){
      #     for (c in 1:nCohorts){
      #       betaInt[y,s,c] ~ dnorm(betaIntYOY[y],1)
      #       betaPhi[y,s,c] ~ dnorm(betaPhiYOY[y],1)
      #       betaP[y,s,c] ~ dnorm(betaPYOY[y],1)
      #     }
      #     
      #     betaFlow[1,y,s] ~ dnorm(betaFlowYOY[1,y],1)
      #     betaFlow[2,y,s] ~ dnorm(betaFlowYOY[2,y],1)
      #     betaFlow[3,y,s] ~ dnorm(betaFlowYOY[3,y],1)
      #     
      #   }
      # }
      
      for (y in 1:2) {
        for (s in 1:nSeasons){
          for (c in 1:nCohorts){
            betaInt[y,s,c] ~ dnorm(betaIntYOYSeason[y,s],1)
            betaPhi[y,s,c] ~ dnorm(betaPhiYOYSeason[y,s],1)
            betaP[y,s,c] ~ dnorm(betaPYOYSeason[y,s],1)
          }
          
          betaIntYOYSeason[y,s] ~ dnorm(betaIntYOY[y],1)
          betaPhiYOYSeason[y,s] ~ dnorm(betaPhiYOY[y],1)
          betaPYOYSeason[y,s] ~ dnorm(betaPYOY[y],1)
          
          betaFlow[1,y,s] ~ dnorm(betaFlowYOY[1,y],1)
          betaFlow[2,y,s] ~ dnorm(betaFlowYOY[2,y],1)
          betaFlow[3,y,s] ~ dnorm(betaFlowYOY[3,y],1)
          
        }
      }
      
      for (y in 1:2) {
        
        betaIntYOY[y] ~ dnorm(betaIntTop,1)
        betaPhiYOY[y] ~ dnorm(betaPhiTop,1)
        betaFlowYOY[1,y] ~ dnorm(betaFlowTop[1],1)
        betaFlowYOY[2,y] ~ dnorm(betaFlowTop[2],1)
        betaFlowYOY[3,y] ~ dnorm(betaFlowTop[3],1)
        
        betaPYOY[y] ~ dnorm(0,1)
      }
      
      betaIntTop ~ dnorm(0,1)
      betaPhiTop ~ dnorm(0,1)
      betaFlowTop[1] ~ dnorm(0,1)
      betaFlowTop[2] ~ dnorm(0,1)
      betaFlowTop[3] ~ dnorm(0,1)
      
      betaPTop ~ dnorm(0,1)
      
      ##    
      # back-transform for examining output
      ##
      for (y in 1:2) {
        for (s in 1:nSeasons){
          # for (c in 1:nCohorts){
          betaIntYOYSeasonOut[y,s] <- 1/(1 + exp(-betaIntYOYSeason[y,s]))
          betaPhiYOYSeasonOut[y,s] <- 1/(1 + exp(-betaPhiYOYSeason[y,s]))
          #  betaFlowOut[1,y,s,c] <- 1/(1 + exp(-betaFlow[1,y,s,c]))
          #  betaFlowOut[2,y,s,c] <- 1/(1 + exp(-betaFlow[2,y,s,c]))
          #  betaFlowOut[3,y,s,c] <- 1/(1 + exp(-betaFlow[3,y,s,c]))
          
          betaPYOYSeasonOut[y,s] <- 1/(1 + exp(-betaPYOYSeason[y,s]))
          #  }
        }
      }
      ##    
      # likelihood
      for (i in 1:N){
        y[i,first[i]:last[i]] ~ dDHMMo(init = delta[1:2],
                                       probTrans = gamma[1:2, 1:2, first[i]:last[i], i],
                                       probObs =   omega[1:2, 1:2, first[i]:last[i], i],
                                       len = length[i],
                                       checkRowSums = 1)
      }
    }),
    
    tt_Rmodel_ft_cohort_OB_flowByRiver = nimbleModel(
      code = tt_modelCode_ft_cohort_OB_flowByRiver,
      constants = tt_myConstants_ft_cohort_OB_flowByRiver,
      data = tt_myData_ft_cohort_OB_flowByRiver,
      inits = tt_initialValues_ft_cohort_OB_flow(tt_myConstants_ft_cohort_OB_flowByRiver$T, tt_myConstants_ft_cohort_OB_flowByRiver$nCohorts, 
                                          tt_inputData_ft_cohort_OB_flowByRiver$y, tt_inputData_ft_cohort_OB_flowByRiver$zInitsNA),
      calculate = FALSE
    ),
    
    tt_parametersToSave_ft_cohort_OB_flowByRiver = c("betaIntTop", "betaPhiTop","betaFlowTop","betaPTop",  
                                       "betaIntYOY", "betaPhiYOY","betaFlowYOY","betaPYOY",
                                       "betaIntYOYSeason", "betaPhiYOYSeason","betaPYOYSeason",
                                       "betaIntYOYSeasonOut", "betaPhiYOYSeasonOut","betaPYOYSeasonOut",
                                       "betaInt", "betaPhi","betaFlow","betaP"
                                       #"betaIntOut", "betaPhiOut","betaFlowOut","betaPOut"
    ),
    
    tt_conf_ft_cohort_OB_flowByRiver = configureMCMC(
      tt_Rmodel_ft_cohort_OB_flowByRiver,
      monitors = tt_parametersToSave_ft_cohort_OB_flowByRiver
    ),
    
    tt_Rmcmc_ft_cohort_OB_flowByRiver = buildMCMC(tt_conf_ft_cohort_OB_flowByRiver, useConjugacy = FALSE),
    tt_Cmodel_ft_cohort_OB_flowByRiver = compileNimble(tt_Rmodel_ft_cohort_OB_flowByRiver),
    tt_Cmcmc_ft_cohort_OB_flowByRiver = compileNimble(tt_Rmcmc_ft_cohort_OB_flowByRiver, project = tt_Rmodel_ft_cohort_OB_flowByRiver),
    
    tt_model_ft_cohort_OB_flowByRiver = runMCMC(
      tt_Cmcmc_ft_cohort_OB_flowByRiver,
      niter = tt_runData_ft_cohort_OB_flowByRiver$nIter,
      nburnin = tt_runData_ft_cohort_OB_flowByRiver$nBurnin,
      thin = tt_runData_ft_cohort_OB_flowByRiver$thinRate,
      nchains = tt_runData_ft_cohort_OB_flowByRiver$nChains
    ),
    
    tt_modelOut_ft_cohort_OB_flowByRiver =
      list(
        mcmc = tt_model_ft_cohort_OB_flowByRiver, # "Error : invalid nimbleFunction argument"
        name = "phiT_pT_ft_cohort_OB_flowByRiver",
        modelCode = tt_modelCode_ft_cohort_OB_flowByRiver,
        myConstants = tt_myConstants_ft_cohort_OB_flowByRiver,
        runData = tt_runData_ft_cohort_OB_flowByRiver
      ),
    
    tt_save_ft_cohort_OB_flowByRiver = saveModelOut_tt_ft_cohort_OB_flowByRiver(tt_modelOut_ft_cohort_OB_flowByRiver)
  )


    
######################################
#### Functions
#####################################
tt_initialValues_OB = function(t, y) {
  list(phi = runif(t - 1, 0, 1),
       p = runif(t - 1, 0, 1),
       z = getInits_tt_ft_cohort_OB(y)
      )
}


tt_initialValues_ft_cohort_OB_flow <- function(t, c, y, z) list(
  betaIntTop = array(rnorm(1, 0, 1), 1),
  betaPhiTop = array(rnorm(1, 0, 1), 1),
  betaFlowTop = array(rnorm(3, 0, 1), 3),
  betaPTop = array(rnorm(1, 0, 1), 1),
  
  betaIntYOY = array(rnorm(2, 0, 1), 2),
  betaPhiYOY = array(rnorm(2, 0, 1), 2),
  betaFlowYOY = array(rnorm(3, 0, 1), c(3,2)),
  betaPYOY = array(rnorm(2, 0, 1), 2),
  
  # betaIntYOYCohort = array(rnorm(2*c, 0, 1), c(2,c)),
  # betaPhiYOYCohort = array(rnorm(2*c, 0, 1), c(2,c)),
  # betaFlowYOYCohort = array(rnorm(2*2*c, 0, 1), c(2,2,c)),
  # betaPYOYCohort = array(rnorm(2*c, 0, 1), c(2,c)),
  
  betaIntYOYSeason = array(rnorm(2*4, 0, 1), c(2,4)),
  betaPhiYOYSeason = array(rnorm(2*4, 0, 1), c(2,4)),
  # betaFlowYOYSeason = array(rnorm(3*2*4, 0, 1), c(3,2,4)),
  betaPYOYSeason = array(rnorm(2*4, 0, 1), c(2,4)),
  
  betaInt = array(rnorm(2*4*c, 0, 1), c(2,4,c)),
  betaPhi = array(rnorm(2*4*c, 0, 1), c(2,4,c)),
  betaFlow = array(rnorm(3*2*4, 0, 1), c(3,2,4)),
  betaP = array(rnorm(2*4*c, 0, 1), c(2,4,c))
)


getInits_tt_ft_cohort_OB <- function(d) {
  d <- d + 1
  d[d == 2] <- 1
  return(d)
}

saveModelOut_tt_ft_cohort_OB <- function(d) {
  save(d, file = paste0('./models/cmrFlowOB/runsOut/tt_ft_cohort_OB_', substr(Sys.time(),1,13), '.RData'))
}
saveModelOut_tt_ft_cohort_OB_flow <- function(d) {
  save(d, file = paste0('./models/cmrFlowOB/runsOut/tt_ft_cohort_OB_flow', substr(Sys.time(),1,13), '.RData'))
}

saveModelOut_tt_ft_cohort_OB_flowByRiver <- function(d) {
  save(d, file = paste0('./models/cmrFlowOB/runsOut/tt_ft_cohort_OB_flowByRiver', substr(Sys.time(),1,13), '.RData'))
}




