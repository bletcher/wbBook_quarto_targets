tar_option_set(packages = c("tidyverse", "nimble", "nimbleEcology", "MCMCvis", "targets"))

#load('./models/cmrFlowOB/dataOut/eh_2002200320042005200620072008200920102011201220132014_wb obear.RData')

#eh_OB_target = tar_read(eh_OB_2002_2014_target)
#eh_OB_2002_2014_target=tar_read(eh_OB_2002_2014_target)

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

#############################################################################
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
                                          eh_OB_2002_2014_target$eh),
        zInitsNA = ifelse(is.na(eh_OB_2002_2014_target$flow), NA, 1),
        flow = eh_OB_2002_2014_target$flow,
        nCohorts = length(unique(eh_OB_2002_2014_target$cohorts$cohort)),
        cohort = eh_OB_2002_2014_target$cohorts$cohort,
        nSeasons = length(unique(eh_OB_2002_2014_target$data$season)),
        seasonArray = c(3,4,1,2,3,4,1,2,3,4,1,2)
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
    
    tt_myConstants_OB_flow0 = list(
      N = nrow(tt_inputData_OB_flow$y),
      T = ncol(tt_inputData_OB_flow$y),
      first = tt_inputData_OB_flow$first,
      last = tt_inputData_OB_flow$last,
      cohort = tt_inputData_OB_flow$cohort - min(tt_inputData_OB_flow$cohort) + 1,
      nCohorts = tt_inputData_OB_flow$nCohorts,
      season = tt_inputData_OB_flow$seasonArray,
      flow = tt_inputData_OB_flow$flow,
      length = tt_inputData_OB_flow$last - tt_inputData_OB_flow$first + 1,
      indToKeep = which(tt_inputData_OB_flow$first < 12)
    ),
    
    tt_myConstants_OB_flow = list(
      N = length(tt_myConstants_OB_flow0$indToKeep),
      T = ncol(tt_inputData_OB_flow$y),
      first = tt_myConstants_OB_flow0$first[tt_myConstants_OB_flow0$indToKeep],
      last = tt_myConstants_OB_flow0$last[tt_myConstants_OB_flow0$indToKeep],
      cohort = tt_myConstants_OB_flow0$cohort[tt_myConstants_OB_flow0$indToKeep], 
      season = tt_myConstants_OB_flow0$season,
      flow = tt_myConstants_OB_flow0$flow[tt_myConstants_OB_flow0$indToKeep,],
      length = tt_myConstants_OB_flow0$length[tt_myConstants_OB_flow0$indToKeep],
      nCohorts = tt_inputData_OB_flow$nCohorts,
      nSeasons = tt_inputData_OB_flow$nSeasons
    ),
    
    tt_myData_OB_flow = list(
      y = tt_inputData_OB_flow$y + 1
    ),
    
    tt_modelCode_OB_flow = nimbleCode({
      # from https://bletcher.github.io/westBrook-book/models.html#model-phit_pt_cohort_flowcohorthierdhmm
      delta[1] <- 1                    # Pr(alive t = 1) = 1
      delta[2] <- 0                    # Pr(dead t = 1) = 0
      
      for (i in 1:N){
        for (t in 1:(T-1)){ # loop over time
          logit(phi[t,i]) <- 
            betaInt +
            betaPhi[t,cohort[i]] + 
            betaFlow[1,season[t],cohort[i]] * flow[i,t] +
            betaFlow[2,season[t],cohort[i]] * flow[i,t] * flow[i,t]
          # prior survival
          ##
          gamma[1,1,t,i] <- phi[t,i]         # Pr(alive t -> alive t+1)
          gamma[1,2,t,i] <- 1 - phi[t,i]     # Pr(alive t -> dead t+1)
          gamma[2,1,t,i] <- 0              # Pr(dead t -> alive t+1)
          gamma[2,2,t,i] <- 1              # Pr(dead t -> dead t+1)
          ##            
          ## DT changes:
          ## definition of omega is moved below, to make it
          ## correctly condition on the first (positive) observation
          ##logit(p[t,i]) <- betaP[t,cohort[i]]             # prior detection
          ##omega[1,1,t,i] <- 1 - p[t,i]       # Pr(alive t -> non-detected t)
          ##omega[1,2,t,i] <- p[t,i]           # Pr(alive t -> detected t)
          ##omega[2,1,t,i] <- 1              # Pr(dead t -> non-detected t)
          ##omega[2,2,t,i] <- 0              # Pr(dead t -> detected t)
        }
        ## DT changes:
        ## need to pad the gamma matrix with an extra t=T row, to ensure it's
        ## always a matrix.  This values are never actually used (except maybe for internal checking of row sums = 1),
        ## but defining them is necessary.
        gamma[1,1,T,i] <- 0
        gamma[1,2,T,i] <- 1
        gamma[2,1,T,i] <- 0
        gamma[2,2,T,i] <- 1
        ## DT changes:
        ## time period t = first[i]: guaranteed detection:
        omega[1,1,first[i],i] <- 0       # Pr(alive t -> non-detected t)
        omega[1,2,first[i],i] <- 1           # Pr(alive t -> detected t)
        omega[2,1,first[i],i] <- 1              # Pr(dead t -> non-detected t)
        omega[2,2,first[i],i] <- 0              # Pr(dead t -> detected t)
        ## DT changes:
        ## time t > first[i]:
        for(t in (first[i]+1):last[i]) {
          logit(p[t,i]) <- betaP[t-1,cohort[i]]             # prior detection
          omega[1,1,t,i] <- 1 - p[t,i]       # Pr(alive t -> non-detected t)
          omega[1,2,t,i] <- p[t,i]           # Pr(alive t -> detected t)
          omega[2,1,t,i] <- 1              # Pr(dead t -> non-detected t)
          omega[2,2,t,i] <- 0              # Pr(dead t -> detected t)
        }
      }
      ##    
      betaInt ~ dnorm(0,1)
      betaFlowTop[1] ~ dnorm(0,1)
      betaFlowTop[2] ~ dnorm(0,1)
      ##    
      for (c in 1:nCohorts){
        # mean values
        betaPhiCohort[c] ~ dnorm(0,1)
        betaPCohort[c] ~ dnorm(0,1)
        betaFlowCohort[1,c] ~ dnorm(betaFlowTop[1],1)
        betaFlowCohort[2,c] ~ dnorm(betaFlowTop[2],1)
        for (t in 1:(T-1)){ 
          betaPhi[t,c] ~ dnorm(betaPhiCohort[c],1)
          betaP[t,c] ~ dnorm(betaPCohort[c],1)
        }
      }
      ##    
      # back-transform for examining output
      for (c in 1:nCohorts){
        betaPhiCohortOut[c] <- 1/(1 + exp(-betaPhiCohort[c]))
        betaPCohortOut[c] <- 1/(1 + exp(-betaPCohort[c]))
        for (t in 1:(T-1)){ 
          betaPhiOut[t,c] <- 1/(1 + exp(-betaPhi[t,c]))
          betaPOut[t,c] <- 1/(1 + exp(-betaP[t,c])) 
        }
      }
      ##    
      for (s in 1:nSeasons){
        for (c in 1:nCohorts){
          betaFlow[1,s,c] ~ dnorm(betaFlowCohort[1,c],1)
          betaFlow[2,s,c] ~ dnorm(betaFlowCohort[2,c],1)
        }   
      }
      ##    
      # likelihood
      for (i in 1:N){
        y[i,first[i]:last[i]] ~ dDHMMo(init = delta[1:2],
                                       probTrans = gamma[1:2, 1:2, first[i]:last[i], i],
                                       probObs = omega[1:2, 1:2, first[i]:last[i], i],
                                       len = length[i],
                                       checkRowSums = 1)
      }
    }),
    
    tt_Rmodel_OB_flow = nimbleModel(
      code = tt_modelCode_OB_flow,
      constants = tt_myConstants_OB_flow,
      data = tt_myData_OB_flow,
      inits = tt_initialValues_OB_flow(tt_myConstants_OB_flow$T, tt_myConstants_OB_flow$nCohorts, 
                                       tt_inputData_OB_flow$y, tt_inputData_OB_flow$zInitsNA),
      calculate = FALSE
    ),
    
    tt_parametersToSave_OB_flow = c("betaInt", 
                                    "betaPhi", "betaP", "betaPhiCohort", "betaPCohort",
                                    "betaPhiOut", "betaPOut", "betaPhiCohortOut", "betaPCohortOut", 
                                    "betaFlow",
                                    "betaFlowCohort", "betaFlowTop"),
    
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
      ),
    
    tt_save_OB_flow = saveModelOut_tt_OB_flow(tt_modelOut_OB_flow)
  )


#############################################################################
# Add flow by river
# Flow estimated for each river independently
###################################################

modelCMR_tt_OB_flowByRiver_target <-
  tar_plan(
    
    tt_inputData_OB_flowByRiver =
      list(
        y = eh_OB_2002_2014_target$eh,
        first = eh_OB_2002_2014_target$first, 
        last = eh_OB_2002_2014_target$last, 
        zInits = tt_initialValues_OB(ncol(eh_OB_2002_2014_target$eh), 
                                     eh_OB_2002_2014_target$eh),
        zInitsNA = ifelse(is.na(eh_OB_2002_2014_target$flow), NA, 1),
        ###### Change input flow data compared to OB_flow model above.
        flow = eh_OB_2002_2014_target$flowByRiver,
        ######
        nCohorts = length(unique(eh_OB_2002_2014_target$cohorts$cohort)),
        cohort = eh_OB_2002_2014_target$cohorts$cohort,
        nSeasons = length(unique(eh_OB_2002_2014_target$data$season)),
        seasonArray = c(3,4,1,2,3,4,1,2,3,4,1,2)
        # nStates = length(unique(eh_OB_2002_2014_target$data$sizeState)),
        # nRivers = length(unique(eh_OB_2002_2014_target$data$sizeState)) # for now
      ), 
    
    tt_runData_OB_flowByRiver = list(
      # Updateable model-specific variables 
      nIter = 5000, 
      nBurnin = 2000, 
      nChains = 2,
      thinRate = 5
    ),  
    
    tt_myConstants_OB_flowByRiver0 = list(
      N = nrow(tt_inputData_OB_flowByRiver$y),
      T = ncol(tt_inputData_OB_flowByRiver$y),
      first = tt_inputData_OB_flowByRiver$first,
      last = tt_inputData_OB_flowByRiver$last,
      cohort = tt_inputData_OB_flowByRiver$cohort - min(tt_inputData_OB_flowByRiver$cohort) + 1,
      nCohorts = tt_inputData_OB_flowByRiver$nCohorts,
      season = tt_inputData_OB_flowByRiver$seasonArray,
      flow = tt_inputData_OB_flowByRiver$flow,
      length = tt_inputData_OB_flowByRiver$last - tt_inputData_OB_flowByRiver$first + 1,
      indToKeep = which(tt_inputData_OB_flowByRiver$first < 12)
    ),
    
    tt_myConstants_OB_flowByRiver = list(
      N = length(tt_myConstants_OB_flowByRiver0$indToKeep),
      T = ncol(tt_inputData_OB_flowByRiver$y),
      first = tt_myConstants_OB_flowByRiver0$first[tt_myConstants_OB_flowByRiver0$indToKeep],
      last = tt_myConstants_OB_flowByRiver0$last[tt_myConstants_OB_flowByRiver0$indToKeep],
      cohort = tt_myConstants_OB_flowByRiver0$cohort[tt_myConstants_OB_flowByRiver0$indToKeep], 
      season = tt_myConstants_OB_flowByRiver0$season,
      flow = tt_myConstants_OB_flowByRiver0$flow[tt_myConstants_OB_flowByRiver0$indToKeep,],
      length = tt_myConstants_OB_flowByRiver0$length[tt_myConstants_OB_flowByRiver0$indToKeep],
      nCohorts = tt_inputData_OB_flowByRiver$nCohorts,
      nSeasons = tt_inputData_OB_flowByRiver$nSeasons
    ),
    
    tt_myData_OB_flowByRiver = list(
      y = tt_inputData_OB_flowByRiver$y + 1
    ),
    
    tt_modelCode_OB_flowByRiver = nimbleCode({
      # from https://bletcher.github.io/westBrook-book/models.html#model-phit_pt_cohort_flowcohorthierdhmm
      delta[1] <- 1                    # Pr(alive t = 1) = 1
      delta[2] <- 0                    # Pr(dead t = 1) = 0
      
      for (i in 1:N){
        for (t in 1:(T-1)){ # loop over time
          logit(phi[t,i]) <- 
            betaInt +
            betaPhi[t,cohort[i]] + 
            betaFlow[1,season[t],cohort[i]] * flow[i,t] +
            betaFlow[2,season[t],cohort[i]] * flow[i,t] * flow[i,t]
          # prior survival
          ##
          gamma[1,1,t,i] <- phi[t,i]         # Pr(alive t -> alive t+1)
          gamma[1,2,t,i] <- 1 - phi[t,i]     # Pr(alive t -> dead t+1)
          gamma[2,1,t,i] <- 0              # Pr(dead t -> alive t+1)
          gamma[2,2,t,i] <- 1              # Pr(dead t -> dead t+1)
          ##            
          ## DT changes:
          ## definition of omega is moved below, to make it
          ## correctly condition on the first (positive) observation
          ##logit(p[t,i]) <- betaP[t,cohort[i]]             # prior detection
          ##omega[1,1,t,i] <- 1 - p[t,i]       # Pr(alive t -> non-detected t)
          ##omega[1,2,t,i] <- p[t,i]           # Pr(alive t -> detected t)
          ##omega[2,1,t,i] <- 1              # Pr(dead t -> non-detected t)
          ##omega[2,2,t,i] <- 0              # Pr(dead t -> detected t)
        }
        ## DT changes:
        ## need to pad the gamma matrix with an extra t=T row, to ensure it's
        ## always a matrix.  This values are never actually used (except maybe for internal checking of row sums = 1),
        ## but defining them is necessary.
        gamma[1,1,T,i] <- 0
        gamma[1,2,T,i] <- 1
        gamma[2,1,T,i] <- 0
        gamma[2,2,T,i] <- 1
        ## DT changes:
        ## time period t = first[i]: guaranteed detection:
        omega[1,1,first[i],i] <- 0       # Pr(alive t -> non-detected t)
        omega[1,2,first[i],i] <- 1           # Pr(alive t -> detected t)
        omega[2,1,first[i],i] <- 1              # Pr(dead t -> non-detected t)
        omega[2,2,first[i],i] <- 0              # Pr(dead t -> detected t)
        ## DT changes:
        ## time t > first[i]:
        for(t in (first[i]+1):last[i]) {
          logit(p[t,i]) <- betaP[t-1,cohort[i]]             # prior detection
          omega[1,1,t,i] <- 1 - p[t,i]       # Pr(alive t -> non-detected t)
          omega[1,2,t,i] <- p[t,i]           # Pr(alive t -> detected t)
          omega[2,1,t,i] <- 1              # Pr(dead t -> non-detected t)
          omega[2,2,t,i] <- 0              # Pr(dead t -> detected t)
        }
      }
      ##    
      betaInt ~ dnorm(0,1)
      betaFlowTop[1] ~ dnorm(0,1)
      betaFlowTop[2] ~ dnorm(0,1)
      ##    
      for (c in 1:nCohorts){
        # mean values
        betaPhiCohort[c] ~ dnorm(0,1)
        betaPCohort[c] ~ dnorm(0,1)
        betaFlowCohort[1,c] ~ dnorm(betaFlowTop[1],1)
        betaFlowCohort[2,c] ~ dnorm(betaFlowTop[2],1)
        for (t in 1:(T-1)){ 
          betaPhi[t,c] ~ dnorm(betaPhiCohort[c],1)
          betaP[t,c] ~ dnorm(betaPCohort[c],1)
        }
      }
      ##    
      # back-transform for examining output
      for (c in 1:nCohorts){
        betaPhiCohortOut[c] <- 1/(1 + exp(-betaPhiCohort[c]))
        betaPCohortOut[c] <- 1/(1 + exp(-betaPCohort[c]))
        for (t in 1:(T-1)){ 
          betaPhiOut[t,c] <- 1/(1 + exp(-betaPhi[t,c]))
          betaPOut[t,c] <- 1/(1 + exp(-betaP[t,c])) 
        }
      }
      ##    
      for (s in 1:nSeasons){
        for (c in 1:nCohorts){
          betaFlow[1,s,c] ~ dnorm(betaFlowCohort[1,c],1)
          betaFlow[2,s,c] ~ dnorm(betaFlowCohort[2,c],1)
        }   
      }
      ##    
      # likelihood
      for (i in 1:N){
        y[i,first[i]:last[i]] ~ dDHMMo(init = delta[1:2],
                                       probTrans = gamma[1:2, 1:2, first[i]:last[i], i],
                                       probObs = omega[1:2, 1:2, first[i]:last[i], i],
                                       len = length[i],
                                       checkRowSums = 1)
      }
    }),
    
    tt_Rmodel_OB_flowByRiver = nimbleModel(
      code = tt_modelCode_OB_flowByRiver,
      constants = tt_myConstants_OB_flowByRiver,
      data = tt_myData_OB_flowByRiver,
      inits = tt_initialValues_OB_flow(tt_myConstants_OB_flowByRiver$T, tt_myConstants_OB_flowByRiver$nCohorts, 
                                       tt_inputData_OB_flowByRiver$y, tt_inputData_OB_flowByRiver$zInitsNA),
      calculate = FALSE
    ),
    
    tt_parametersToSave_OB_flowByRiver = c("betaInt", 
                                    "betaPhi", "betaP", "betaPhiCohort", "betaPCohort",
                                    "betaPhiOut", "betaPOut", "betaPhiCohortOut", "betaPCohortOut", 
                                    "betaFlow",
                                    "betaFlowCohort", "betaFlowTop"),
    
    tt_conf_OB_flowByRiver = configureMCMC(
      tt_Rmodel_OB_flowByRiver,
      monitors = tt_parametersToSave_OB_flowByRiver
    ),
    
    tt_Rmcmc_OB_flowByRiver = buildMCMC(tt_conf_OB_flowByRiver, useConjugacy = FALSE),
    tt_Cmodel_OB_flowByRiver = compileNimble(tt_Rmodel_OB_flowByRiver),
    tt_Cmcmc_OB_flowByRiver = compileNimble(tt_Rmcmc_OB_flowByRiver, project = tt_Rmodel_OB_flowByRiver),
    
    tt_model_OB_flowByRiver = runMCMC(
      tt_Cmcmc_OB_flowByRiver,
      niter = tt_runData_OB_flowByRiver$nIter,
      nburnin = tt_runData_OB_flowByRiver$nBurnin,
      thin = tt_runData_OB_flowByRiver$thinRate,
      nchains = tt_runData_OB_flowByRiver$nChains
    ),
    
    tt_modelOut_OB_flowByRiver =
      list(
        mcmc = tt_model_OB_flowByRiver, # "Error : invalid nimbleFunction argument"
        name = "phiT_pT_OB_flowByRiver",
        modelCode = tt_modelCode_OB_flowByRiver,
        myConstants = tt_myConstants_OB_flowByRiver,
        runData = tt_runData_OB_flowByRiver
      ),
    
    tt_save_OB_flowByRiver = saveModelOut_tt_OB_flowByRiver(tt_modelOut_OB_flowByRiver)
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


tt_initialValues_OB_flow  <- function(t, c, y, z) list(
  betaInt = rnorm(1, 0, 1),
  ## DT change:
  ## don't give phi and p initial values;
  ## they're deterministic nodes, so they'll be calculated
  ## in terms of other variables.  Also, when I made changes to these
  ## in the code, these (unnecessary) initial values were the wrong sizes
  ##phi = array(runif((t - 1) * myConstants$N, 0, 1),c((t - 1), myConstants$N)),
  ##p =   array(runif((t - 1) * myConstants$N, 0, 1),c((t - 1), myConstants$N)),
  #z = getInits_tt_OB(z),
  betaPhi = array(runif((t - 1) * c, 0, 1), c((t - 1), c)),
  betaP =   array(runif((t - 1) * c, 0, 1), c((t - 1), c)),
  betaPhiCohort = array(runif(c, 0, 1),c(c)),
  betaPCohort =   array(runif(c, 0, 1),c(c)),
  betaFlow = array(rnorm(2 * 4 * c, 0, 1), c(2, 4, c)),
  betaFlowCohort = array(rnorm(2 * c, 0, 1), c(2, c)),
  betaFlowTop = rnorm(2, 0, 1)
)


getInits_tt_OB <- function(d) {
  d <- d + 1
  d[d == 2] <- 1
  return(d)
}

saveModelOut_tt_OB <- function(d) {
  save(d, file = paste0('./models/cmrFlowOB/runsOut/tt_OB_', substr(Sys.time(),1,13), '.RData'))
}
saveModelOut_tt_OB_flow <- function(d) {
  save(d, file = paste0('./models/cmrFlowOB/runsOut/tt_OB_flow', substr(Sys.time(),1,13), '.RData'))
}
saveModelOut_tt_OB_flowByRiver <- function(d) {
  save(d, file = paste0('./models/cmrFlowOB/runsOut/tt_OB_flowByRiver', substr(Sys.time(),1,13), '.RData'))
}




