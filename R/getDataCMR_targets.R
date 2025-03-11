tar_option_set(packages = c("tidyverse", "lubridate", "getWBData", "data.table", 
                            "validate", "rlang", "stringi", "writexl", "hrbrthemes", "viridis", "targets"))

# maximum ageInSamples for both createCmrData and getEH
maxAgeInSamples <- 12
maxAgeInSamples2 <- 12 + 8

# https://stackoverflow.com/questions/69583424/using-tidy-eval-for-multiple-arbitrary-filter-conditions
# Assumes LHS is the name of a variable and OP is
# the name of a function
cols <- list("cohort")
ops <-  list("%in%")
vals <- list(2002:2014)

dataCMR_WB_2002_2014_target <-
  tar_plan(
    cdWB_CMR0_target = 
      createCoreData(
        sampleType = "electrofishing", #"stationaryAntenna","portableAntenna"),
        whichDrainage = "west",
        columnsToAdd =
          c("sampleNumber",
            "river",
            "riverMeter",
            "survey",
            "pass",
            'observedLength',
            'observedWeight')
      ) %>%
      addTagProperties(
        columnsToAdd =
          c("cohort",
            "species",
            "dateEmigrated",
            "sex",
            "species")
      ) %>%
      dplyr::filter(!is.na(tag),
                    area %in% c("trib","inside","below","above"),
                    !is.na(sampleNumber)
      ) %>%
      createCmrData(maxAgeInSamples = maxAgeInSamples + 1, # +1 so we get env data for the last interval
                    inside = F,
                    censorDead = F,
                    censorEmigrated = F) %>% # may want to change censorEmigrated = T to = F
      # sample 83 is the last tagging sample

      filter(sampleNumber <= 83) %>%
      addSampleProperties() %>%
      arrange(tag, ageInSamples) |> # so we get correct lagDetectionDate in addEnvironmental3
      addEnvironmental3(envDataWBIn = envDataWB_target, envDataWB_fdcThreshIn = envDataWB_fdcThresh_target) %>%
      # these commented functions below do not work for CMR data - they separate out shock and non-shock samples
      #addEnvironmentalDaily() %>%
      #addEnvironmentalInterval() %>%
      addKnownZ2() %>%
      addFirstLast() %>%
      fillRiver() %>%
      addRiverTagged() %>%
      scaleEnvData() %>%
      addIsYOY() %>%
      addRiverN() %>%
      mutate(
        riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                                     labels = c("West Brook","WB Jimmy","WB Mitchell","WB OBear"), ordered = T),
        sampleInterval = as.numeric(lagDetectionDate - detectionDate)
      ),
    
    #@ with longer maxAgeinSamples
    cdWB_CMR0_2_target = 
      createCoreData(
        sampleType = "electrofishing", #"stationaryAntenna","portableAntenna"),
        whichDrainage = "west",
        columnsToAdd =
          c("sampleNumber",
            "river",
            "riverMeter",
            "survey",
            "pass",
            'observedLength',
            'observedWeight')
      ) %>%
      addTagProperties(
        columnsToAdd =
          c("cohort",
            "species",
            "dateEmigrated",
            "sex",
            "species")
      ) %>%
      dplyr::filter(!is.na(tag),
                    area %in% c("trib","inside","below","above"),
                    !is.na(sampleNumber)
      ) %>%
      createCmrData(maxAgeInSamples = maxAgeInSamples2 + 1, # +1 so we get env data for the last interval
                    inside = F,
                    censorDead = F,
                    censorEmigrated = F) %>% # may want to change censorEmigrated = T to = F
      # sample 83 is the last tagging sample
      
      filter(sampleNumber <= 83) %>%
      addSampleProperties() %>%
      arrange(tag, ageInSamples) |> # so we get correct lagDetectionDate in addEnvironmental3
      addEnvironmental3(envDataWBIn = envDataWB_target, envDataWB_fdcThreshIn = envDataWB_fdcThresh_target) %>%
      # these commented functions below do not work for CMR data - they separate out shock and non-shock samples
      #addEnvironmentalDaily() %>%
      #addEnvironmentalInterval() %>%
      addKnownZ2() %>%
      addFirstLast() %>%
      fillRiver() %>%
      addRiverTagged() %>%
      scaleEnvData() %>%
      addIsYOY() %>%
      addRiverN() %>%
      addRiverAbbrev() |>
      mutate(
        riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                              labels = c("West Brook","WB Jimmy","WB Mitchell","WB OBear"), ordered = T),
        sampleInterval = as.numeric(lagDetectionDate - detectionDate)
      ),

   eh_WB_2002_2014_target = getEH_AIS(cdWB_CMR0_target, cols, ops, vals, maxAgeInSamples),#, maxIndexByCohort = 100)
   eh_WB_1997_2014_target = getEH_AIS(cdWB_CMR0_target, cols, ops, 1997:2014, maxAgeInSamples)
   )  

###############################################################
## WB - BKT fish
# read down through the cols, ops, vals variables for filter conditions
cols_wBbkt <- list("cohort",  "species")
ops_wBbkt <-  list("%in%",    "==")
vals_wBbkt <- list(2002:2014, "bkt")
######################################  


dataCMR_WBbkt_2002_2014_target <-
  tar_plan(
    eh_WBbkt_2002_2014_target = getEH_AIS(cdWB_CMR0_target, cols_wBbkt, ops_wBbkt, vals_wBbkt, maxAgeInSamples)
  )

###############################################################
###############################################################
## Loop over river and cohorts - BKT fish
# read down through the cols, ops, vals variables for filter conditions
cols_bkt <- list("cohort",  "species", "riverAbbrev")
ops_bkt <-  list("%in%",    "==",      "%in%")
vals_bkt <- list(2002:2015, "bkt",     c("WB", "OL", "OS", "IS"))

######################################
# Loop over cohorts and rivers
dataCMR_bkt_target <- tar_map(
  values = tidyr::crossing(
    riverAbbrev = c("WB", "OL", "OS", "IS"), 
    cohort = 2002:2015
  ),
  tar_target(
    name = eh_bkt,#paste0("dataCMR_bkt_river", riverN, "_cohort", cohort),
    getEH_AIS(
      cdWB_CMR0_2_target, 
      cols_bkt, 
      ops_bkt, 
      list(cohort, "bkt", riverAbbrev),
      maxAgeInSamples2
    )
  )
)

################################################################
## WB - BNT fish
# read down through the cols, ops, vals variables for filter conditions
cols_wBbnt <- list("cohort",  "species")
ops_wBbnt <-  list("%in%",    "==")
vals_wBbnt <- list(2002:2014, "bnt")
######################################  


dataCMR_WBbnt_2002_2014_target <-
  tar_plan(
    eh_WBbnt_2002_2014_target = getEH_AIS(cdWB_CMR0_target, cols_wBbnt, ops_wBbnt, vals_wBbnt, maxAgeInSamples)#, maxIndexByCohort = 100)
  )

###############################################################
###############################################################
## Loop over river and cohorts - BNT fish
# read down through the cols, ops, vals variables for filter conditions
cols_bnt <- list("cohort",  "species", "riverAbbrev")
ops_bnt <-  list("%in%",    "==",      "%in%")
vals_bnt <- list(2002:2015, "bnt",     c("WB", "OL", "OS"))

######################################
# Loop over cohorts and rivers
dataCMR_bnt_target <- tar_map(
  values = tidyr::crossing(
    riverAbbrev = c("WB", "OL", "OS"), 
    cohort = 2002:2015
  ),
  tar_target(
    name = eh_bnt,#paste0("dataCMR_bkt_river", riverN, "_cohort", cohort),
    getEH_AIS(
      cdWB_CMR0_2_target, 
      cols_bnt, 
      ops_bnt, 
      list(cohort, "bnt", riverAbbrev),
      maxAgeInSamples2
    )
  )
)

################################################################
## WB - ATS fish
# read down through the cols, ops, vals variables for filter conditions
cols_wBats <- list("cohort",  "species")
ops_wBats <-  list("%in%",    "==")
vals_wBats <- list(1997:2004, "ats")
######################################  

#tmp = tar_read(cdWB_CMR0_target) |> filter(species == "ats")
#table(tmp$cohort)

dataCMR_WBats_1997_2004_target <-
  tar_plan(
    eh_WBats_1997_2004_target = getEH_AIS(cdWB_CMR0_target, cols_wBats, ops_wBats, vals_wBats, maxAgeInSamples)#, maxIndexByCohort = 100)
  )

###############################################################
###############################################################
## Loop over river and cohorts - ATS fish
# read down through the cols, ops, vals variables for filter conditions
cols_ats <- list("cohort",  "species", "riverAbbrev")
ops_ats <-  list("%in%",    "==",      "%in%")
vals_ats <- list(1997:2004, "ats",     c("WB", "OL", "OS"))

######################################
# Loop over cohorts and rivers
dataCMR_ats_target <- tar_map(
  values = tidyr::crossing(
    riverAbbrev = c("WB", "OL", "OS"), 
    cohort = 1997:2004
  ),
  tar_target(
    name = eh_ats,#paste0("dataCMR_bkt_river", riverN, "_cohort", cohort),
    getEH_AIS(
      cdWB_CMR0_2_target, 
      cols_ats, 
      ops_ats, 
      list(cohort, "ats", riverAbbrev),
      maxAgeInSamples2
    )
  )
)

#################################################################################
## OB fish

# #########################################
# # all cohorts from O'Bear 2002:2014

# read down through the cols, ops, vals variables for filter conditions
cols_OB <- list("cohort",  "riverTagged")
ops_OB <-  list("%in%",    "==")
vals_OB <- list(2002:2014, "wb obear")
######################################  

dataCMR_OB_2002_2014_target <-
  tar_plan(
    eh_OB_2002_2014_target = getEH_AIS_Filter_first(cdWB_CMR0_target, cols_OB, ops_OB, vals_OB, maxAgeInSamples),#, maxIndexByCohort = 100)
    eh_OB_2002_2014_AISall_target = getEH_AIS_Filter_first(cdWB_CMR0_2_target, cols_OB, ops_OB, vals_OB, maxAgeInSamples2)
  )  


# #########################################
# # specific cohorts from O'Bear 2002:2014

# read down through the cols, ops, vals variables for filter conditions
cols_OB <- list("cohort",  "riverTagged")
ops_OB <-  list("%in%",    "==")
#vals_OB <- list(2010, "wb obear")
######################################  

dataCMR_OB_singleCohorts_target <-
  tar_plan(
    eh_OB_2002_target = getEH_AIS(cdWB_CMR0_target, cols_OB, ops_OB, list(2002, "wb obear"), maxAgeInSamples), # don't need getEH_AIS_Filter_first() because cohorts are independent and limited to observed AIS values
    eh_OB_2003_target = getEH_AIS(cdWB_CMR0_target, cols_OB, ops_OB, list(2003, "wb obear"), maxAgeInSamples),
    eh_OB_2004_target = getEH_AIS(cdWB_CMR0_target, cols_OB, ops_OB, list(2004, "wb obear"), maxAgeInSamples),
    eh_OB_2005_target = getEH_AIS(cdWB_CMR0_target, cols_OB, ops_OB, list(2005, "wb obear"), maxAgeInSamples),
    eh_OB_2006_target = getEH_AIS(cdWB_CMR0_target, cols_OB, ops_OB, list(2006, "wb obear"), maxAgeInSamples),
    eh_OB_2007_target = getEH_AIS(cdWB_CMR0_target, cols_OB, ops_OB, list(2007, "wb obear"), maxAgeInSamples),
    eh_OB_2008_target = getEH_AIS(cdWB_CMR0_target, cols_OB, ops_OB, list(2008, "wb obear"), maxAgeInSamples),
    eh_OB_2009_target = getEH_AIS(cdWB_CMR0_target, cols_OB, ops_OB, list(2009, "wb obear"), maxAgeInSamples),
    eh_OB_2010_target = getEH_AIS(cdWB_CMR0_target, cols_OB, ops_OB, list(2010, "wb obear"), maxAgeInSamples),
    eh_OB_2011_target = getEH_AIS(cdWB_CMR0_target, cols_OB, ops_OB, list(2011, "wb obear"), maxAgeInSamples),
    eh_OB_2012_target = getEH_AIS(cdWB_CMR0_target, cols_OB, ops_OB, list(2012, "wb obear"), maxAgeInSamples),
    eh_OB_2013_target = getEH_AIS(cdWB_CMR0_target, cols_OB, ops_OB, list(2013, "wb obear"), maxAgeInSamples),
    eh_OB_2014_target = getEH_AIS(cdWB_CMR0_target, cols_OB, ops_OB, list(2014, "wb obear"), maxAgeInSamples)
  ) 



#####################################
## dataCMR functions 
#####################################

#addEnvironmental3() is in 'generalFunctions.R'

getKnown <- function(x) {
  firstObs <- min(which(x == 1))
  lastObs <- max(which(x == 1))
  known <- rep(0, length(x))
  known[firstObs:lastObs] <- 1
  if (lastObs != length(known)) {
    known[(lastObs + 1):length(known)] <- NA
  }
  return(known)
}

addKnownZ2 <- function(d) {
  d %>% 
    group_by(tag) %>%
    arrange(sampleNumber) %>%
    mutate(knownZ = getKnown(enc)) %>%
    ungroup() %>%
    arrange(tag, sampleNumber)
}

addFirstLast <- function(d){
  firstLast <- d %>% 
    group_by(tag) %>%
    filter(knownZ == 1) %>%
    summarize(firstObserved = min(sampleNumber, na.rm = TRUE),
              lastObserved = max(sampleNumber, na.rm = TRUE)) %>%
    ungroup()
  
  left_join(d, firstLast) %>%
    mutate(isFirstObserved = sampleNumber == firstObserved,
           isLastObserved = sampleNumber == lastObserved)
} 

addFirstLastAIS <- function(d){
  firstLast <- d %>% 
    group_by(tag) %>%
    filter(knownZ == 1) %>%
    summarize(firstObservedAIS = min(ageInSamples, na.rm = TRUE),
              lastObservedAIS = max(ageInSamples, na.rm = TRUE)) %>%
    ungroup()
  
  left_join(d, firstLast) %>%
    mutate(isFirstObservedAIS = ageInSamples == firstObservedAIS,
           isLastObservedAIS = ageInSamples == lastObservedAIS)
}

getMaxAISByCohort <- function(d) {
  maxAIS <- d |>
    group_by(cohort, ageInSamples) |>
    summarize(
      count = n()
    ) |>
    group_by(cohort) |>
    summarize(maxAISByCohort = max(ageInSamples)) |>
    ungroup()
  
  d <- left_join(d, maxAIS)
}

removeFishLastPossibleAIS <- function(d) {
  d |> 
    group_by(tag) |> 
    filter(firstObservedAIS < maxAISByCohort - 1) |> # minus one because last is defined as minus one in getEH_AIs
    ungroup()
}

fillRiver <- function (data, location = T){
  fillLocation <- function(location) {
    known <- which(!is.na(location))
    unknown <- which(is.na(location))
    nKnown <- length(unique(location[known]))
    if (nKnown == 1) {
      location[unknown] <- location[known[1]]
    }
    else {
      for (i in unknown) {
        location[i] <- location[known[max(which(i > known))]]
      }
    }
    return(location)
  }
  if (location) {
    data <- data %>% 
      group_by(tag) %>% 
      mutate(river = fillLocation(river)) %>% 
      ungroup()
  }
  
  return(data)
}

scale_this <- function(x){
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

scaleEnvData <- function(d){
  tmp <- d %>%
    group_by(river, season) %>% 
    summarize(meanMeanFlow = mean(meanFlow, na.rm = TRUE),
              sdMeanFlow = sd(meanFlow, na.rm = TRUE),
              meanMeanFlowByRiver = mean(meanFlowByRiver, na.rm = TRUE),
              sdMeanFlowByRiver = sd(meanFlowByRiver, na.rm = TRUE),
              
              meanMeanFlowByArea_flowExt = mean(meanFlowByArea_flowExt, na.rm = TRUE),
              sdMeanFlowByArea_flowExt = sd(meanFlowByArea_flowExt, na.rm = TRUE),
              meanMeanFlowByArea_ByRiver = mean(meanFlowByArea_ByRiver, na.rm = TRUE),
              sdMeanFlowByArea_ByRiver = sd(meanFlowByArea_ByRiver, na.rm = TRUE),
              
              meanMeanTemperature = mean(meanTemperature, na.rm = TRUE),
              sdMeanTemperature = sd(meanTemperature, na.rm = TRUE)
    ) %>%
    ungroup()
  
  out <- left_join(d, tmp) %>%
    mutate(
      meanFlowScaled = (meanFlow - meanMeanFlow) / sdMeanFlow,
      meanFlowByRiverScaled = (meanFlowByRiver - meanMeanFlowByRiver) / sdMeanFlowByRiver,
      meanFlowByArea_flowExtScaled = (meanFlowByArea_flowExt - meanMeanFlowByArea_flowExt) / sdMeanFlowByArea_flowExt,
      meanFlowByArea_ByRiverScaled = (meanFlowByArea_ByRiver - meanMeanFlowByArea_ByRiver) / sdMeanFlowByArea_ByRiver,
      meanTemperatureScaled = (meanTemperature - meanMeanTemperature) / sdMeanTemperature
    )  
}

getNeverCaptured <- function(d, maxOccasionValue){
  d %>%
    #filter(ageInSamples > 0 & ageInSamples <= maxOccasionValue) %>%
    filter(ageInSamples %in% 1:maxOccasionValue) %>%
    group_by(tag) %>%
    summarize(sumEnc = sum(enc, na.rm = TRUE)) %>%
    filter(sumEnc == 0) %>%
    dplyr::select(tag)
}

addRiverTagged <- function(d){
  d1 <- d %>% 
    filter(isFirstObserved) %>%
    mutate(riverTagged = river) %>%
    dplyr::select(tag, riverTagged)
  
  return(left_join(d, d1))
  
}

addIsYOY <- function(d){
  d %>%
    mutate(isYOY = ifelse(ageInSamples <= 3, 1, 2))
  
}

addRiverN <- function(d){
  level_key <- c("west brook" = 1, "wb jimmy" = 2, "wb mitchell" = 3, "wb obear" = 4)
  d %>% 
    mutate(riverN = recode(river, !!!level_key))
}

addRiverAbbrev <- function(d){
  level_key <- c("west brook" = "WB", "wb jimmy" = "OL", "wb mitchell" = "OS", "wb obear" = "IS")
  d %>% 
    mutate(riverAbbrev = recode(river, !!!level_key))
}

`%notin%` <- Negate(`%in%`)

addLaggedSection <- function(d) {
  d0 <- d %>% 
    group_by(tag) %>% 
    arrange(ageInSamples) %>%
    mutate(sectionN = as.numeric(section)) %>%
    mutate(lagSectionN = lead(sectionN)) %>% 
    ungroup()
  
  d <- left_join(d, d0)
}

addLaggedSectionRiverN <- function(d) {
  d0 <- d %>% 
    group_by(tag) %>% 
    arrange(ageInSamples) %>%
    mutate(sectionRiverN = as.numeric(sectionRiverN)) %>%
    mutate(lagSectionRiverN = lead(sectionRiverN)) %>% 
    ungroup()
  
  d <- left_join(d, d0)
}

#####################################
## EH functions 
#####################################
# https://stackoverflow.com/questions/69583424/using-tidy-eval-for-multiple-arbitrary-filter-conditions
# Assumes LHS is the name of a variable and OP is
# the name of a function
op_call <- function(op, lhs, rhs) {
  call(op, sym(lhs), rhs)
}

ehFilter <- function(data, cols, ops, vals) {
  exprs <- purrr::pmap(list(ops, cols, vals), op_call)
  data %>% dplyr::filter(!!!exprs)
}


# var is the variable to put in the encounter history (e.g. 'enc' or 'temp')
# occasionVar is now fixed to ageInSamples
# maxOccasionValue is the maximum value for occasion columns, in units of 'occasionVar'
getEHDataWide_AIS <- function(d, cols, ops, vals, var, maxOccasionValue, valuesFill = 0){   
  d %>%
    ehFilter(cols, ops, vals) %>% 
    #filter(ageInSamples > 0, ageInSamples <= maxOccasionValue) %>%
    arrange(ageInSamples) %>% #need this to get correct order of columns. 
    pivot_wider(
      id_cols = tag,
      names_from = ageInSamples,
      names_prefix = "ais_",
      values_from = eval(substitute(var)),
      values_fill = valuesFill
    )
}

getEH_AIS <- function(dIn, cols, ops, vals, maxOccasionValue, maxIndexByCohort = 1E10){
  
  d <- dIn %>% 
    #filter(ageInSamples > 0, ageInSamples <= maxOccasionValue)
    filter(ageInSamples %in% 1:maxOccasionValue) 
   # addFirstLast() |>
  #  filter(ageInSamples != maxOccasionValue, !isFirstObserved) # remove fish seen for the first time on the last occasion

  # Fish with no observed occasions
  neverCaptured <- getNeverCaptured(d, maxOccasionValue)
  d <- d %>%
    filter(tag %notin% neverCaptured$tag)
  
  # limit data to first 'maxIndexByCohort' individuals for each cohort
  d <- d %>%
    group_by(cohort) %>%
    mutate(indexByCohort = rleid(tag)) %>%
    filter(indexByCohort <= maxIndexByCohort) %>%
    ungroup()
  
  # for sectionRiverN for Xioawei
  d <- d %>%
    mutate(
      sectionTMP = ifelse(riverN == 3 & section < 1, 1, section),
      sectionRiverN = case_when(
        riverN == 1 & area == "below" ~ 0,
        riverN == 1 & area == "inside" ~ as.numeric(sectionTMP),
        riverN == 1 & area == "above" ~ 48,
        riverN == 2 ~ as.numeric(sectionTMP) + 48,
        riverN == 3 ~ as.numeric(sectionTMP) + 63,
        riverN == 4 ~ as.numeric(sectionTMP) + 79,
      )
    )
  #table(dIn2$riverN,dIn2$sectionRiverN, dIn2$area)
  
  encWide <- getEHDataWide_AIS(d, cols, ops, vals, "enc", maxOccasionValue, valuesFill = 0)
  eh <- as.matrix(encWide %>% dplyr::select(-tag), nrow = nrow(encWide), ncol = ncol(encWide) - 1)
  
  flowFill <- 0
  flowWide <- getEHDataWide_AIS(d, cols, ops, vals, "meanFlowScaled", maxOccasionValue, valuesFill = flowFill)
  flowMatrix <- as.matrix(flowWide %>% dplyr::select(-tag), nrow = nrow(flowWide), ncol = ncol(flowWide) - 1)
  flowMatrix <- ifelse(is.finite(flowMatrix), flowMatrix, flowFill)
  
  flowByRiverFill <- 0
  flowByRiverWide <- getEHDataWide_AIS(d, cols, ops, vals, "meanFlowByRiverScaled", maxOccasionValue, valuesFill = flowByRiverFill)
  flowByRiverMatrix <- as.matrix(flowByRiverWide %>% dplyr::select(-tag), nrow = nrow(flowByRiverWide), ncol = ncol(flowByRiverWide) - 1)
  flowByRiverMatrix <- ifelse(is.finite(flowByRiverMatrix), flowByRiverMatrix, flowByRiverFill)
  
  flowByArea_flowExtFill <- 0
  flowByArea_flowExtWide <- getEHDataWide_AIS(d, cols, ops, vals, "meanFlowByArea_flowExtScaled", maxOccasionValue, valuesFill = flowByArea_flowExtFill)
  flowByArea_flowExtMatrix <- as.matrix(flowByArea_flowExtWide %>% dplyr::select(-tag), nrow = nrow(flowByArea_flowExtWide), ncol = ncol(flowByArea_flowExtWide) - 1)
  flowByArea_flowExtMatrix <- ifelse(is.finite(flowByArea_flowExtMatrix), flowByArea_flowExtMatrix, flowByArea_flowExtFill)
  
  flowByArea_ByRiverFill <- 0
  flowByArea_ByRiverWide <- getEHDataWide_AIS(d, cols, ops, vals, "meanFlowByArea_ByRiverScaled", maxOccasionValue, valuesFill = flowByArea_ByRiverFill)
  flowByArea_ByRiverMatrix <- as.matrix(flowByArea_ByRiverWide %>% dplyr::select(-tag), nrow = nrow(flowByArea_ByRiverWide), ncol = ncol(flowByArea_ByRiverWide) - 1)
  flowByArea_ByRiverMatrix <- ifelse(is.finite(flowByArea_ByRiverMatrix), flowByArea_ByRiverMatrix, flowByArea_ByRiverFill)
  
  temperatureFill <- 0
  temperatureWide <- getEHDataWide_AIS(d, cols, ops, vals, "meanTemperatureScaled", maxOccasionValue, valuesFill = temperatureFill)
  temperatureMatrix <- as.matrix(temperatureWide %>% dplyr::select(-tag), nrow = nrow(temperatureWide), ncol = ncol(temperatureWide) - 1)
  temperatureMatrix <- ifelse(is.finite(temperatureMatrix), temperatureMatrix, temperatureFill)
  
  propBelowLoFlowThreshByRiverFill <- 0
  propBelowLoFlowThreshByRiverWide <- getEHDataWide_AIS(d, cols, ops, vals, "propBelowLoFlowThreshByRiver", maxOccasionValue, valuesFill = propBelowLoFlowThreshByRiverFill)
  propBelowLoFlowThreshByRiverMatrix <- as.matrix(propBelowLoFlowThreshByRiverWide %>% dplyr::select(-tag), nrow = nrow(propBelowLoFlowThreshByRiverWide), ncol = ncol(propBelowLoFlowThreshByRiverWide) - 1)
  propBelowLoFlowThreshByRiverMatrix <- ifelse(is.finite(propBelowLoFlowThreshByRiverMatrix), propBelowLoFlowThreshByRiverMatrix, propBelowLoFlowThreshByRiverFill)
  
  propAboveHiFlowThreshByRiverFill <- 0
  propAboveHiFlowThreshByRiverWide <- getEHDataWide_AIS(d, cols, ops, vals, "propAboveHiFlowThreshByRiver", maxOccasionValue, valuesFill = propAboveHiFlowThreshByRiverFill)
  propAboveHiFlowThreshByRiverMatrix <- as.matrix(propAboveHiFlowThreshByRiverWide %>% dplyr::select(-tag), nrow = nrow(propAboveHiFlowThreshByRiverWide), ncol = ncol(propAboveHiFlowThreshByRiverWide) - 1)
  propAboveHiFlowThreshByRiverMatrix <- ifelse(is.finite(propAboveHiFlowThreshByRiverMatrix), propAboveHiFlowThreshByRiverMatrix, propAboveHiFlowThreshByRiverFill)
  
  propBelowLoFlowThreshByArea_flowExtFill <- 0
  propBelowLoFlowThreshByArea_flowExtWide <- getEHDataWide_AIS(d, cols, ops, vals, "propBelowLoFlowThreshByArea_flowExt", maxOccasionValue, valuesFill = propBelowLoFlowThreshByArea_flowExtFill)
  propBelowLoFlowThreshByArea_flowExtMatrix <- as.matrix(propBelowLoFlowThreshByArea_flowExtWide %>% dplyr::select(-tag), nrow = nrow(propBelowLoFlowThreshByArea_flowExtWide), ncol = ncol(propBelowLoFlowThreshByArea_flowExtWide) - 1)
  propBelowLoFlowThreshByArea_flowExtMatrix <- ifelse(is.finite(propBelowLoFlowThreshByArea_flowExtMatrix), propBelowLoFlowThreshByArea_flowExtMatrix, propBelowLoFlowThreshByArea_flowExtFill)
  
  propAboveHiFlowThreshByArea_flowExtFill <- 0
  propAboveHiFlowThreshByArea_flowExtWide <- getEHDataWide_AIS(d, cols, ops, vals, "propAboveHiFlowThreshByArea_flowExt", maxOccasionValue, valuesFill = propAboveHiFlowThreshByArea_flowExtFill)
  propAboveHiFlowThreshByArea_flowExtMatrix <- as.matrix(propAboveHiFlowThreshByArea_flowExtWide %>% dplyr::select(-tag), nrow = nrow(propAboveHiFlowThreshByArea_flowExtWide), ncol = ncol(propAboveHiFlowThreshByArea_flowExtWide) - 1)
  propAboveHiFlowThreshByArea_flowExtMatrix <- ifelse(is.finite(propAboveHiFlowThreshByArea_flowExtMatrix), propAboveHiFlowThreshByArea_flowExtMatrix, propAboveHiFlowThreshByArea_flowExtFill)
  
  riverWide <- getEHDataWide_AIS(d, cols, ops, vals, "river", maxOccasionValue, valuesFill = NA)
  riverMatrix <- as.matrix(riverWide %>% dplyr::select(-tag), nrow = nrow(riverWide), ncol = ncol(riverWide) - 1)
  
  riverNWide <- getEHDataWide_AIS(d, cols, ops, vals, "riverN", maxOccasionValue, valuesFill = 0)
  riverNMatrix <- as.matrix(riverNWide %>% dplyr::select(-tag), nrow = nrow(riverNWide), ncol = ncol(riverNWide) - 1)
  
  sectionRiverNWide <- getEHDataWide_AIS(d, cols, ops, vals, "sectionRiverN", maxOccasionValue, valuesFill = NA)
  sectionRiverNMatrix <- as.matrix(sectionRiverNWide %>% dplyr::select(-tag), nrow = nrow(sectionRiverNWide), ncol = ncol(sectionRiverNWide) - 1)
  
  isYOYWide <- getEHDataWide_AIS(d, cols, ops, vals, "isYOY", maxOccasionValue, valuesFill = 2)
  isYOYMatrix <- as.matrix(isYOYWide %>% dplyr::select(-tag), nrow = nrow(isYOYWide), ncol = ncol(riverWide) - 1)
  
  lengthWide <- getEHDataWide_AIS(d, cols, ops, vals, "observedLength", maxOccasionValue, valuesFill = NA)
  lengthMatrix <- as.matrix(lengthWide %>% dplyr::select(-tag), nrow = nrow(lengthWide), ncol = ncol(lengthWide) - 1)
  
  weightWide <- getEHDataWide_AIS(d, cols, ops, vals, "observedWeight", maxOccasionValue, valuesFill = NA)
  weightMatrix <- as.matrix(weightWide %>% dplyr::select(-tag), nrow = nrow(weightWide), ncol = ncol(weightWide) - 1)
  
  sampleIntervalWide <- getEHDataWide_AIS(d, cols, ops, vals, "sampleInterval", maxOccasionValue, valuesFill = NA)
  sampleIntervalMatrix <- as.matrix(sampleIntervalWide %>% dplyr::select(-tag), nrow = nrow(sampleIntervalWide), ncol = ncol(sampleIntervalWide) - 1)
  
  seasonWide <- getEHDataWide_AIS(d, cols, ops, vals, "season", maxOccasionValue, valuesFill = NA)
  seasonMatrix <- as.matrix(seasonWide %>% dplyr::select(-tag), nrow = nrow(seasonWide), ncol = ncol(seasonWide) - 1)
  
  yearWide <- getEHDataWide_AIS(d, cols, ops, vals, "year", maxOccasionValue, valuesFill = NA)
  yearMatrix <- as.matrix(yearWide %>% dplyr::select(-tag), nrow = nrow(yearWide), ncol = ncol(yearWide) - 1)
  
  tags <- encWide %>% dplyr::select(tag)
  
  data <- d %>%
    ehFilter(cols, ops, vals) %>% 
    #filter(ageInSamples > 0, ageInSamples <= maxOccasionValue) %>%
    filter(ageInSamples %in% 1:maxOccasionValue) %>%
    arrange(tag, ageInSamples)
  
  cohorts <- tags %>% left_join(data %>% dplyr::select(tag, cohort) %>% unique()) %>% dplyr::select(cohort)
  seasons <- tags %>% left_join(data %>% dplyr::select(tag, season) %>% unique()) %>% dplyr::select(season)
  species <- tags %>% left_join(data %>% dplyr::select(tag, species) %>% unique()) %>% dplyr::select(species)
  
  first <- apply(eh, 1, function(x) min(which(x != 0)))
  last <- apply(riverMatrix, 1, function(x) max(which(!is.na(x)))) # riverMatrix fills in river values between frist and last
  last <- ifelse(last == maxOccasionValue, last, last - 1) #commented out because this leads to first=last and backwards indexing. we have last[i]-1 in the model
  
  meanWeight_AIS <- d |> 
    group_by(ageInSamples) |> 
    summarize(meanWeight = mean(observedWeight, na.rm = TRUE))
  
  return(list(eh = eh, flow = flowMatrix, flowByRiver = flowByRiverMatrix, flowByArea_flowExt = flowByArea_flowExtMatrix, 
              flowByArea_ByRiver = flowByArea_ByRiverMatrix, 
              propBelowLoFlowThreshByRiver = propBelowLoFlowThreshByRiverMatrix, propAboveHiFlowThreshByRiver = propAboveHiFlowThreshByRiverMatrix,
              propBelowLoFlowThreshByArea_flowExt = propBelowLoFlowThreshByArea_flowExtMatrix, propAboveHiFlowThreshByArea_flowExt = propAboveHiFlowThreshByArea_flowExtMatrix,
              temperature = temperatureMatrix, river = riverMatrix, section = sectionRiverNMatrix,
              riverN = riverNMatrix, isYOY = isYOYMatrix, length = lengthMatrix, weight = weightMatrix,
              sampleInterval = sampleIntervalMatrix, season = seasonMatrix, year = yearMatrix,
              tags = tags, cohorts = cohorts, seasons = seasons, species = species,
              first = first, last = last, meanWeight_AIS = meanWeight_AIS, data = data))
}

# filter out fish that we caught for the first time on the last or 2nd to last occasion
getEH_AIS_Filter_first <- function(dIn, cols, ops, vals, maxOccasionValue, maxIndexByCohort = 1E10){
  
  d <- dIn %>% 
    #filter(ageInSamples > 0, ageInSamples <= maxOccasionValue)
    filter(ageInSamples %in% 1:maxOccasionValue) |>
    ehFilter(cols, ops, vals) |>
    getMaxAISByCohort() |>
    addFirstLastAIS() |>
    removeFishLastPossibleAIS()
  

#> table(d2$cohort,d2$ageInSamples)

#        1   2   3   4   5   6   7   8   9  10  11  12
# 2002  44  44 139 193 204 212 218 218 219 219 220 220
# 2003 111 145 185 239 256 259 261 262 263 263 263 263
# 2004  61  75  85 152 165 169 172 173 174 174 174 174
# 2005  29  43  50  59  66  68  68  68  68  68  68  68
# 2006  60  90  99 119 121 121 122 124 125 125 125 125
# 2007  38  64  71  88  89  93  93  93  93  93  94  94
# 2008  40  61  67  72  75  80  81  82  82  82  82  82
# 2009 193 281 332 458 488 499 501 503 503 504 504 505
# 2010   8  18  28  39  44  47  47  48  48  48  48  48
# 2011  22  31  35  35  36  37  37  38  38  38  38  38
# 2012 266 391 440 500 520 525 525 526 527 527 527 527
# 2013  54  77  79  86  89  89  89  89  89   0   0   0
# 2014  91 134 148 177 177   0   0   0   0   0   0   0

  
      #filter(ageInSamples > 0, ageInSamples <= maxOccasionValue) %>%

  #addFirstLast() 
  #  filter(ageInSamples != maxOccasionValue, !isFirstObserved) # remove fish seen for the first time on the last occasion
  
  # Fish with no observed occasions
  neverCaptured <- getNeverCaptured(d, maxOccasionValue)
  d <- d %>%
    filter(tag %notin% neverCaptured$tag)
  
  # limit data to first 'maxIndexByCohort' individuals for each cohort
  d <- d %>%
    group_by(cohort) %>%
    mutate(indexByCohort = rleid(tag)) %>%
    filter(indexByCohort <= maxIndexByCohort) %>%
    ungroup()
  
  # for sectionRiverN for Xioawei
  d <- d %>%
    mutate(
      sectionTMP = ifelse(riverN == 3 & section < 1, 1, section),
      sectionRiverN = case_when(
        riverN == 1 & area == "below" ~ 0,
        riverN == 1 & area == "inside" ~ as.numeric(sectionTMP),
        riverN == 1 & area == "above" ~ 48,
        riverN == 2 ~ as.numeric(sectionTMP) + 48,
        riverN == 3 ~ as.numeric(sectionTMP) + 63,
        riverN == 4 ~ as.numeric(sectionTMP) + 79,
      )
    )
  #table(dIn2$riverN,dIn2$sectionRiverN, dIn2$area)
  
  encWide <- getEHDataWide_AIS(d, cols, ops, vals, "enc", maxOccasionValue, valuesFill = 0)
  eh <- as.matrix(encWide %>% dplyr::select(-tag), nrow = nrow(encWide), ncol = ncol(encWide) - 1)
  
  flowFill <- 0
  flowWide <- getEHDataWide_AIS(d, cols, ops, vals, "meanFlowScaled", maxOccasionValue, valuesFill = flowFill)
  flowMatrix <- as.matrix(flowWide %>% dplyr::select(-tag), nrow = nrow(flowWide), ncol = ncol(flowWide) - 1)
  flowMatrix <- ifelse(is.finite(flowMatrix), flowMatrix, flowFill)
  
  flowByRiverFill <- 0
  flowByRiverWide <- getEHDataWide_AIS(d, cols, ops, vals, "meanFlowByRiverScaled", maxOccasionValue, valuesFill = flowByRiverFill)
  flowByRiverMatrix <- as.matrix(flowByRiverWide %>% dplyr::select(-tag), nrow = nrow(flowByRiverWide), ncol = ncol(flowByRiverWide) - 1)
  flowByRiverMatrix <- ifelse(is.finite(flowByRiverMatrix), flowByRiverMatrix, flowByRiverFill)
  
  flowByArea_flowExtFill <- 0
  flowByArea_flowExtWide <- getEHDataWide_AIS(d, cols, ops, vals, "meanFlowByArea_flowExtScaled", maxOccasionValue, valuesFill = flowByArea_flowExtFill)
  flowByArea_flowExtMatrix <- as.matrix(flowByArea_flowExtWide %>% dplyr::select(-tag), nrow = nrow(flowByArea_flowExtWide), ncol = ncol(flowByArea_flowExtWide) - 1)
  flowByArea_flowExtMatrix <- ifelse(is.finite(flowByArea_flowExtMatrix), flowByArea_flowExtMatrix, flowByArea_flowExtFill)
  
  flowByArea_ByRiverFill <- 0
  flowByArea_ByRiverWide <- getEHDataWide_AIS(d, cols, ops, vals, "meanFlowByArea_ByRiverScaled", maxOccasionValue, valuesFill = flowByArea_ByRiverFill)
  flowByArea_ByRiverMatrix <- as.matrix(flowByArea_ByRiverWide %>% dplyr::select(-tag), nrow = nrow(flowByArea_ByRiverWide), ncol = ncol(flowByArea_ByRiverWide) - 1)
  flowByArea_ByRiverMatrix <- ifelse(is.finite(flowByArea_ByRiverMatrix), flowByArea_ByRiverMatrix, flowByArea_ByRiverFill)
  
  temperatureFill <- 0
  temperatureWide <- getEHDataWide_AIS(d, cols, ops, vals, "meanTemperatureScaled", maxOccasionValue, valuesFill = temperatureFill)
  temperatureMatrix <- as.matrix(temperatureWide %>% dplyr::select(-tag), nrow = nrow(temperatureWide), ncol = ncol(temperatureWide) - 1)
  temperatureMatrix <- ifelse(is.finite(temperatureMatrix), temperatureMatrix, temperatureFill)
  
  propBelowLoFlowThreshByRiverFill <- 0
  propBelowLoFlowThreshByRiverWide <- getEHDataWide_AIS(d, cols, ops, vals, "propBelowLoFlowThreshByRiver", maxOccasionValue, valuesFill = propBelowLoFlowThreshByRiverFill)
  propBelowLoFlowThreshByRiverMatrix <- as.matrix(propBelowLoFlowThreshByRiverWide %>% dplyr::select(-tag), nrow = nrow(propBelowLoFlowThreshByRiverWide), ncol = ncol(propBelowLoFlowThreshByRiverWide) - 1)
  propBelowLoFlowThreshByRiverMatrix <- ifelse(is.finite(propBelowLoFlowThreshByRiverMatrix), propBelowLoFlowThreshByRiverMatrix, propBelowLoFlowThreshByRiverFill)
  
  propAboveHiFlowThreshByRiverFill <- 0
  propAboveHiFlowThreshByRiverWide <- getEHDataWide_AIS(d, cols, ops, vals, "propAboveHiFlowThreshByRiver", maxOccasionValue, valuesFill = propAboveHiFlowThreshByRiverFill)
  propAboveHiFlowThreshByRiverMatrix <- as.matrix(propAboveHiFlowThreshByRiverWide %>% dplyr::select(-tag), nrow = nrow(propAboveHiFlowThreshByRiverWide), ncol = ncol(propAboveHiFlowThreshByRiverWide) - 1)
  propAboveHiFlowThreshByRiverMatrix <- ifelse(is.finite(propAboveHiFlowThreshByRiverMatrix), propAboveHiFlowThreshByRiverMatrix, propAboveHiFlowThreshByRiverFill)
  
  propBelowLoFlowThreshByArea_flowExtFill <- 0
  propBelowLoFlowThreshByArea_flowExtWide <- getEHDataWide_AIS(d, cols, ops, vals, "propBelowLoFlowThreshByArea_flowExt", maxOccasionValue, valuesFill = propBelowLoFlowThreshByArea_flowExtFill)
  propBelowLoFlowThreshByArea_flowExtMatrix <- as.matrix(propBelowLoFlowThreshByArea_flowExtWide %>% dplyr::select(-tag), nrow = nrow(propBelowLoFlowThreshByArea_flowExtWide), ncol = ncol(propBelowLoFlowThreshByArea_flowExtWide) - 1)
  propBelowLoFlowThreshByArea_flowExtMatrix <- ifelse(is.finite(propBelowLoFlowThreshByArea_flowExtMatrix), propBelowLoFlowThreshByArea_flowExtMatrix, propBelowLoFlowThreshByArea_flowExtFill)
  
  propAboveHiFlowThreshByArea_flowExtFill <- 0
  propAboveHiFlowThreshByArea_flowExtWide <- getEHDataWide_AIS(d, cols, ops, vals, "propAboveHiFlowThreshByArea_flowExt", maxOccasionValue, valuesFill = propAboveHiFlowThreshByArea_flowExtFill)
  propAboveHiFlowThreshByArea_flowExtMatrix <- as.matrix(propAboveHiFlowThreshByArea_flowExtWide %>% dplyr::select(-tag), nrow = nrow(propAboveHiFlowThreshByArea_flowExtWide), ncol = ncol(propAboveHiFlowThreshByArea_flowExtWide) - 1)
  propAboveHiFlowThreshByArea_flowExtMatrix <- ifelse(is.finite(propAboveHiFlowThreshByArea_flowExtMatrix), propAboveHiFlowThreshByArea_flowExtMatrix, propAboveHiFlowThreshByArea_flowExtFill)
  
  riverWide <- getEHDataWide_AIS(d, cols, ops, vals, "river", maxOccasionValue, valuesFill = NA)
  riverMatrix <- as.matrix(riverWide %>% dplyr::select(-tag), nrow = nrow(riverWide), ncol = ncol(riverWide) - 1)
  
  riverNWide <- getEHDataWide_AIS(d, cols, ops, vals, "riverN", maxOccasionValue, valuesFill = 0)
  riverNMatrix <- as.matrix(riverNWide %>% dplyr::select(-tag), nrow = nrow(riverNWide), ncol = ncol(riverNWide) - 1)
  
  sectionRiverNWide <- getEHDataWide_AIS(d, cols, ops, vals, "sectionRiverN", maxOccasionValue, valuesFill = NA)
  sectionRiverNMatrix <- as.matrix(sectionRiverNWide %>% dplyr::select(-tag), nrow = nrow(sectionRiverNWide), ncol = ncol(sectionRiverNWide) - 1)
  
  isYOYWide <- getEHDataWide_AIS(d, cols, ops, vals, "isYOY", maxOccasionValue, valuesFill = 2)
  isYOYMatrix <- as.matrix(isYOYWide %>% dplyr::select(-tag), nrow = nrow(isYOYWide), ncol = ncol(riverWide) - 1)
  
  lengthWide <- getEHDataWide_AIS(d, cols, ops, vals, "observedLength", maxOccasionValue, valuesFill = NA)
  lengthMatrix <- as.matrix(lengthWide %>% dplyr::select(-tag), nrow = nrow(lengthWide), ncol = ncol(lengthWide) - 1)
  
  weightWide <- getEHDataWide_AIS(d, cols, ops, vals, "observedWeight", maxOccasionValue, valuesFill = NA)
  weightMatrix <- as.matrix(weightWide %>% dplyr::select(-tag), nrow = nrow(weightWide), ncol = ncol(weightWide) - 1)
  
  sampleIntervalWide <- getEHDataWide_AIS(d, cols, ops, vals, "sampleInterval", maxOccasionValue, valuesFill = NA)
  sampleIntervalMatrix <- as.matrix(sampleIntervalWide %>% dplyr::select(-tag), nrow = nrow(sampleIntervalWide), ncol = ncol(sampleIntervalWide) - 1)
  
  seasonWide <- getEHDataWide_AIS(d, cols, ops, vals, "season", maxOccasionValue, valuesFill = NA)
  seasonMatrix <- as.matrix(seasonWide %>% dplyr::select(-tag), nrow = nrow(seasonWide), ncol = ncol(seasonWide) - 1)
  
  yearWide <- getEHDataWide_AIS(d, cols, ops, vals, "year", maxOccasionValue, valuesFill = NA)
  yearMatrix <- as.matrix(yearWide %>% dplyr::select(-tag), nrow = nrow(yearWide), ncol = ncol(yearWide) - 1)
  
  tags <- encWide %>% dplyr::select(tag)
  
  data <- d %>%
    ehFilter(cols, ops, vals) %>% 
    #filter(ageInSamples > 0, ageInSamples <= maxOccasionValue) %>%
    filter(ageInSamples %in% 1:maxOccasionValue) %>%
    arrange(tag, ageInSamples)
  
  cohorts <- tags %>% left_join(data %>% dplyr::select(tag, cohort) %>% unique()) %>% dplyr::select(cohort)
  seasons <- tags %>% left_join(data %>% dplyr::select(tag, season) %>% unique()) %>% dplyr::select(season)
  species <- tags %>% left_join(data %>% dplyr::select(tag, species) %>% unique()) %>% dplyr::select(species)
  
  first <- apply(eh, 1, function(x) min(which(x != 0)))
  last <- apply(riverMatrix, 1, function(x) max(which(!is.na(x))))
  last <- ifelse(last == maxOccasionValue, last, last - 1) #commented out because this leads to first=last and backwards indexing. we have last[i]-1 in the model
  
  meanWeight_AIS <- d |> 
    group_by(ageInSamples) |> 
    summarize(meanWeight = mean(observedWeight, na.rm = TRUE))
  
  return(list(eh = eh, flow = flowMatrix, flowByRiver = flowByRiverMatrix, flowByArea_flowExt = flowByArea_flowExtMatrix, 
              flowByArea_ByRiver = flowByArea_ByRiverMatrix, 
              propBelowLoFlowThreshByRiver = propBelowLoFlowThreshByRiverMatrix, propAboveHiFlowThreshByRiver = propAboveHiFlowThreshByRiverMatrix,
              propBelowLoFlowThreshByArea_flowExt = propBelowLoFlowThreshByArea_flowExtMatrix, propAboveHiFlowThreshByArea_flowExt = propAboveHiFlowThreshByArea_flowExtMatrix,
              temperature = temperatureMatrix, river = riverMatrix, section = sectionRiverNMatrix,
              riverN = riverNMatrix, isYOY = isYOYMatrix, length = lengthMatrix, weight = weightMatrix,
              sampleInterval = sampleIntervalMatrix, season = seasonMatrix, year = yearMatrix,
              tags = tags, cohorts = cohorts, seasons = seasons, species = species,
              first = first, last = last, meanWeight_AIS = meanWeight_AIS, data = data))
}
