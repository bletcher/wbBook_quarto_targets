tar_option_set(packages = c("tidyverse", "nimble", "nimbleEcology", "MCMCvis", "targets"))

# This is a 'fork' of modelCMR_tt_ft_cohort_OB, following Daniel's March 2023 model update to the WB (4 river) model

# Run this to load data for testing
#eh_WB_2002_2014_target=tar_read(eh_wB_2002_2014_target)

# Model naming
# mod1 = ttt_ft_cohort_WB
# mod2 = ttt_fByRiverT_cohort_WB

#Source Daniel Turek's dDHMMo2 function, specifically made for this model to reduce model size and allow running on a laptop
source('./r/dDHMMo2.R')

run_ttt_models_target <- tar_plan(
  ttt_flow = run_ttt_model(),
  ttt_flowBrRiver = run_ttt_model(flowByRiver = TRUE)
  )


######################################
#### Functions
#####################################


run_ttt_model <- function(eh = tar_read(eh_WBbkt_2002_2014_target), flowByRiver = FALSE) {

  y <- eh$eh * eh$riverN
  (nCohorts <- nrow(unique(eh$cohorts)))
  (nSeasons <- nrow(unique(eh$seasons)))
  (nRivers <- length(unique(eh$data$riverN)))# rivers 1:4
  seasonArray <- c(3,4,1,2,3,4,1,2,3,4,1,2)
  ##
  ##
  first <- eh$first #apply(y, 1, function(x) min(which(x !=0)))
  last <- eh$last
  cohort = ((eh$cohorts) - min(eh$cohorts) + 1)$cohort #can't be a data frame or tibble
  ##
  zinits <- y + 1 # non-detection -> alive
  zinits[zinits == 2] <- 1 # dead -> alive
  zInitsNA <- ifelse(is.na(eh$flow), NA, 1)
  ##
  # Proportion of fish in each river on the first observation
  y1 <- y[,1]
  deltaProps <- table(y1[y1>0]) / length(y1[y1>0])
  ##
  # alpha for the dirichlet prior
  alpha <- c(1,1,1,1)

  ##
  # Priors for psi where more likely to stay than move
  alphaR <- list()
  alphaR[[1]] <- alphaR1 <- c(0.7, 0.1, 0.1, 0.1)
  alphaR[[2]] <- alphaR2 <- c(0.1, 0.7, 0.1, 0.1)
  alphaR[[3]] <- alphaR3 <- c(0.1, 0.1, 0.7, 0.1)
  alphaR[[4]] <- alphaR4 <- c(0.1, 0.1, 0.1, 0.7)
  ##

  ## model code using DHMMo2 distribution
  hmm.phiT_pT_psiT_DHMM_dirch <- nimbleCode({
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
  })
  
  myConstants0 <- list(N = nrow(y),
                       T = ncol(y),
                       first = first,
                       last = last,
                       cohort = cohort,
                       nCohorts = nCohorts,
                       nRivers = nRivers,
                       season = seasonArray,
                       
                       flow = if(flowByRiver) {
                         eh$flowByRiver
                       } else {
                         eh$flow
                       },
                       
                       temp = eh$temperature,
                       length = last - first + 1,
                       alphaR1 = alphaR[[1]],
                       alphaR2 = alphaR[[2]],
                       alphaR3 = alphaR[[3]],
                       alphaR4 = alphaR[[4]],
                       deltaProps = deltaProps
  )
  
  ## DT changes:
  myData0 <- list(###yCJS = eh$eh, #y,    ## data for CJS distribution
    y = y + 1
  )   ## data for DHMM distribution
  
  ## if you change this FALSE to TRUE
  ## this makes the dataset smaller - only 200 observations,
  ## for quicker testing
  makeDataSet20Ind <- FALSE
  if(makeDataSet20Ind) {
    newN <- 20
    oldN <- dim(y)[1]
    set.seed(0)
    indToKeep <- sample(1:oldN, size = newN, replace = FALSE)
  }
  
  ## this removes fish that were only observed on the very last observation
  if(!makeDataSet20Ind) {
    indToKeep <- which(first < ncol(y))
    newN <- length(indToKeep)
  }
  
  myConstants <- list(
    N = newN,
    T = myConstants0$T,
    first = myConstants0$first[indToKeep],
    last = myConstants0$last[indToKeep],
    nRivers = myConstants0$nRivers,
    cohort = myConstants0$cohort[indToKeep],
    nCohorts = myConstants0$nCohorts,
    season = myConstants0$season,
    flow = myConstants0$flow[indToKeep,],
    temp = myConstants0$temp[indToKeep,],
    length = myConstants0$length[indToKeep],
    alphaR1 = alphaR1,
    alphaR2 = alphaR2,
    alphaR3 = alphaR3,
    alphaR4 = alphaR4,
    deltaProps = deltaProps
  )
  
  
  myData <- list(
    ##yCJS = myData0$yCJS[indToKeep,],
    y = myData0$y[indToKeep,]
  )
 
  ## NEW:
  ## code blocks below, which serve to remove {i,t} indexing from phi and gamma
  ## 
  uniqueFlow <- unique(as.numeric(myConstants$flow))
  uniqueTemp <- unique(as.numeric(myConstants$temp))
  uniqueCohort <- unique(as.numeric(myConstants$cohort))
  uniqueSeason <- unique(as.numeric(myConstants$season))
  uniqueTime <- as.numeric(1:myConstants$T)
  ##
  flowIndex <- apply(myConstants$flow, c(1,2), function(x) which(x == uniqueFlow))
  tempIndex <- apply(myConstants$temp, c(1,2), function(x) which(x == uniqueTemp))
  cohortArray <- matrix(myConstants$cohort, myConstants$N, myConstants$T)
  cohortIndex <- apply(cohortArray, c(1,2), function(x) which(x == uniqueCohort))
  seasonArray <- matrix(myConstants$season, myConstants$N, myConstants$T, byrow = TRUE)
  seasonIndex <- apply(seasonArray, c(1,2), function(x) which(x == uniqueSeason))
  timeArray <- matrix(uniqueTime, myConstants$N, myConstants$T, byrow = TRUE)
  timeIndex <- apply(timeArray, c(1,2), function(x) which(x == uniqueTime))
  ##
  indCombined <- paste0(
    as.numeric(flowIndex), '-',
    as.numeric(tempIndex), '-',
    as.numeric(cohortIndex), '-',
    as.numeric(seasonIndex), '-',
    as.numeric(timeIndex))
  ##
  uniqueIndCombined <- unique(indCombined)
  ###length(uniqueIndCombined)    ## 2385 unique combindations of {flow,temp}
  ###                             ## 4425 unique combindations of {flow,temp,cohort}
  ###                             ## 4464 unique combindations of {flow,temp,cohort,season}
  ###                             ## 4541 unique combindationsof {flow,temp,cohort,season,time}
  ##
  nK <- length(uniqueIndCombined)
  ##
  uniqueIndSplitList <- strsplit(uniqueIndCombined, '-')
  ##
  flowK   <- uniqueFlow[as.numeric(sapply(uniqueIndSplitList, `[`, 1))]
  tempK   <- uniqueTemp[as.numeric(sapply(uniqueIndSplitList, `[`, 2))]
  cohortK <- uniqueCohort[as.numeric(sapply(uniqueIndSplitList, `[`, 3))]
  seasonK <- uniqueSeason[as.numeric(sapply(uniqueIndSplitList, `[`, 4))]
  timeK   <- as.numeric(sapply(uniqueIndSplitList, `[`, 5))
  ##
  ## index this by kLookupTable[i,t], to get the value of k corresponding to {i,t}
  kLookupTable <- array(match(indCombined, uniqueIndCombined), c(myConstants$N, myConstants$T))
  ##
  ## check that kLookupTable works correctly:
  if(FALSE) {
    for(i in 1:myConstants$N) {
      for(t in 1:myConstants$T) {
        k <- kLookupTable[i,t]
        if(!all(c(as.numeric(myConstants$flow[i,t]) == flowK[k],
                  as.numeric(myConstants$temp[i,t]) == tempK[k],
                  myConstants$cohort[i] == cohortK[k],
                  myConstants$season[t] == seasonK[k],
                  t == timeK[k]))) stop()
      }
    }
  }
  ##
  ##
  myConstants$nK <- nK
  myConstants$flowK <- flowK
  myConstants$tempK <- tempK
  myConstants$cohortK <- cohortK
  myConstants$seasonK <- seasonK
  myConstants$timeK <- timeK
  myConstants$kLookupTable <- kLookupTable
  myConstants$flow <- NULL
  myConstants$temp <- NULL
  myConstants$season <- NULL
  ##
  
  set.seed(0)
  inits <- initialValues(myConstants$nRivers, myConstants$nCohorts,
                myConstants$nK, myConstants$T,
                myConstants, myConstants$alphaR)
  
  constants <- myConstants
  data <- myData

  (start = Sys.time())
  
  ## you'll get warnings that the data 'yCJS' is not used, and the 'z' initial
  ## values are not in the model.  Those don't cause any problems,
  ## and let us use the same myData and initialValue() for both models.
  system.time(
    Rmodel <- nimbleModel(
      code = hmm.phiT_pT_psiT_DHMM_dirch,
      constants = myConstants,
      data = myData,              
      inits = inits,
      calculate = FALSE
    )
  )
  
  #Rmodel$calculate()   ## 10051.23 (for 20 ind)  ## -97602.24 (for all ind)
  
  parametersToSave <- c("phi", #"betaPhiRiver", "betaPhiRiverCohort", 
                        "betaP", #  "betaPRiver",   "betaPRiverCohort",
                        ##"phiOut",# "betaPhiRiverOut", "betaPhiRiverCohortOut", 
                        #"betaPOut", #  "betaPRiverOut",   "betaPRiverCohortOut",
                        "psi"
  )
  
  nIter <- 12500
  nBurnin <- 2500
  nChains <- 2
  thinRate <- 5
  
  #rm(conf, Rmcmc, Cmodel, Cmcmc) # so old versions don't run if there is an error in an earlier step
  print('conf')
  system.time(
    conf <- configureMCMC(
      Rmodel,
      monitors = parametersToSave,
      useConjugacy = FALSE
    )
  )
  
  print('Rmcmc')
  system.time(Rmcmc <- buildMCMC(conf, useConjugacy = FALSE))
  
  print('Cmodel')
  system.time(Cmodel <- compileNimble(Rmodel))
  
  print('Cmcmc')
  system.time(Cmcmc <- compileNimble(Rmcmc, project = Rmodel))
  
  print('runMCMC')
  system.time(
    mcmc <- runMCMC(
      Cmcmc, 
      niter = nIter, 
      nburnin = nBurnin, 
      thin = thinRate, 
      nchains = nChains
    )
  )
  
  end <- Sys.time()
  (elapsed <- end - start)
  
  toSave <- list(
    mcmc = mcmc, 
    elapsed = elapsed,
    name = paste0('ttt_flow', flowByRiver, sep = "_"),
    flowByRiver = flowByRiver,
    myConstants = myConstants, 
    nIter = nIter, 
    nBurnin = nBurnin,
    thinRate = thinRate, 
    nSeasons = nSeasons, 
    nCohorts = nCohorts,
    nChains = nChains
  )

  return(toSave)
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

#     

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

initialValues <- function(r, c, nK, t, constants, a){
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

