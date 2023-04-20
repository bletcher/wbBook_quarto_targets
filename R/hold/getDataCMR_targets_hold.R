tar_option_set(packages = c("tidyverse", "lubridate", "getWBData", "data.table", 
                            "validate", "rlang", "stringi", "writexl", "hrbrthemes", "viridis", "targets"))

# maximum ageInSamples for both createCmrData and getEH
maxAgeInSamples <- 12

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
      addEnvironmental() %>%
      # these commented functions below do not work for CMR data - they separate out shock and non-shock samples
      #addEnvironmentalDaily() %>%
      #addEnvironmentalInterval() %>%
      addKnownZ2() %>%
      addFirstLast() %>%
      fillRiver() %>%
      addRiverTagged() %>%
      scaleEnvData() %>%
      addIsYOY() %>%
      addRiverN()%>%
      mutate(riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                                   labels = c("West Brook","WB Jimmy","WB Mitchell","WB OBear"), ordered = T)
      ),

   eh_WB_2002_2014_target = getEH_AIS(cdWB_CMR0_target, cols, ops, vals, maxAgeInSamples)#, maxIndexByCohort = 100)
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
    eh_WBbkt_2002_2014_target = getEH_AIS(tar_read(cdWB_CMR0_target), cols_wBbkt, ops_wBbkt, vals_wBbkt, maxAgeInSamples)#, maxIndexByCohort = 100)
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
    eh_WBbnt_2002_2014_target = getEH_AIS(tar_read(cdWB_CMR0_target), cols_wBbnt, ops_wBbnt, vals_wBbnt, maxAgeInSamples)#, maxIndexByCohort = 100)
  )



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
    eh_OB_2002_2014_target = getEH_AIS(tar_read(cdWB_CMR0_target), cols_OB, ops_OB, vals_OB, maxAgeInSamples)#, maxIndexByCohort = 100)
  )  



#####################################
## dataCMR functions 
#####################################

# This function is a copy of addEnvironmetal3() in getDataElectro_targets
# Add any updates there
addEnvironmental <- function(coreData, sampleFlow = F, funName = "mean") {
  reconnect()
  
  func <- get(funName)
  whichDrainage <- "west"
  if (all(!unique(coreData$river) %in% c("west brook", "wb jimmy", 
                                         "wb mitchell", "wb obear"))) {
    whichDrainage <- "stanley"
  }
  if (whichDrainage == "west") {
    envData <- 
      tar_read(envDataWB_target) |> 
      dplyr::filter(date <= max(coreData$detectionDate), 
                    date >= min(coreData$detectionDate)) |> 
      data.frame()
    
    envDataWB_fdcThresh <- tar_read(envDataWB_fdcThresh_target)
      
      # tbl(conDplyr, "data_daily_temperature") %>% 
      # collect(n = Inf) %>% 
      # full_join(tbl(conDplyr, "data_flow_extension") %>% 
      #             collect(n = Inf), by = c("river", "date")) %>% 
      # dplyr::select(-source) %>% 
      # dplyr::filter(date <= max(coreData$detectionDate), 
      #               date >= min(coreData$detectionDate)) %>% 
      # rename(temperature = daily_mean_temp, flow = qPredicted) %>% 
      # data.frame() %>%
      # mutate(dateDate = date(date)) %>%
      # left_join(tar_read(flowByRiver_target))
  }
  else {
    envData <- tbl(conDplyr, "stanley_environmental") %>% 
      filter(section == 11) %>% 
      dplyr::select(datetime, temperature, depth) %>% 
      collect(n = Inf) %>% 
      rename(flow = depth, date = datetime) %>% 
      data.frame()
    warning("Depth was inserted into flow column because that is what is available in Stanley")
  }
  
  coreData <- coreData %>% 
    group_by(tag) %>% 
    arrange(ageInSamples) %>%
    mutate(lagDetectionDate = lead(detectionDate)) %>% 
    ungroup()
  
  if (whichDrainage == "west") {
    getIntervalMean <- function(start, end, r, e, fun = func) {
      d <- envData$date
      if (e == "Temperature") {
        envCol <- "temperature"
        if (is.na(r)) 
          meanVar <- fun(envData[d >= start & d <= end, envCol], na.rm = T)
        if (!is.na(r)) 
          meanVar <- fun(envData[d >= start & d <= end & envData$river == r, envCol], na.rm = T)
      }
      # will need to make this river-specific
      if (e == "Flow") {
        envCol <- "flow"
        meanVar <- fun(envData[d >= start & d <= end, envCol], na.rm = T)
      }
      if (e == "FlowByRiver") {
        envCol <- "flowByRiver"
        if (is.na(r)) 
          meanVar <- fun(envData[d >= start & d <= end, envCol], na.rm = T)
        if (!is.na(r)) 
          meanVar <- fun(envData[d >= start & d <= end & envData$river == r, envCol], na.rm = T)
      }
      if (e == "FlowByArea_flowExt") {
        envCol <- "flowByArea_flowExt"
        if (is.na(r)) 
          meanVar <- fun(envData[d >= start & d <= end, envCol], na.rm = T)
        if (!is.na(r)) 
          meanVar <- fun(envData[d >= start & d <= end & envData$river == r, envCol], na.rm = T)
      }
      if (e == "FlowByArea_ByRiver") {
        envCol <- "flowByArea_ByRiver"
        if (is.na(r)) 
          meanVar <- fun(envData[d >= start & d <= end, envCol], na.rm = T)
        if (!is.na(r)) 
          meanVar <- fun(envData[d >= start & d <= end & envData$river == r, envCol], na.rm = T)
      }
      return(meanVar)
    }
    #######Function for flow duration curve data
    getIntervalMeanFDC <- function(start, end, r, e, fun = func) {
      d <- envDataWB_fdcThresh$date
      if(is.na(r)) {
        dat <- envDataWB_fdcThresh[d >= start & d <= end, e]
        prop <- sum(dat, na.rm = TRUE) / sum(!is.na(dat))
      } else {
        dat <- envDataWB_fdcThresh[d >= start & d <= end & envDataWB_fdcThresh$river == r, e]
        prop <- sum(dat, na.rm = TRUE) / sum(!is.na(dat))
      }
      return(prop)
    }
    
    coreDataUniqueDates <- coreData %>% 
      dplyr::select(river, detectionDate, lagDetectionDate) %>% 
      unique() %>% 
      group_by(river, detectionDate, lagDetectionDate) %>% 
      mutate(meanTemperature = getIntervalMean(detectionDate, lagDetectionDate, river, "Temperature"), 
             meanFlow =        getIntervalMean(detectionDate, lagDetectionDate, river, "Flow"),
             meanFlowByRiver = getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByRiver"),
             meanFlowByArea_flowExt  = getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByArea_flowExt"),
             meanFlowByArea_ByRiver  = getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByArea_ByRiver"),
             sdFlow =          getIntervalMean(detectionDate, lagDetectionDate, river, "Flow", get("sd")),
             sdFlowByRiver =   getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByRiver", get("sd")),
             sdFlowByArea_flowExt  =   getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByArea_flowExt", get("sd")),
             sdFlowByArea_ByRiver  =   getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByArea_ByRiver", get("sd")),
             
             meanBelowLowFlowThreshFlowByRiver = getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "belowLowFlowThreshFlowByRiver"),
             meanAboveHighFlowThreshFlowByRiver = getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "aboveHighFlowThreshFlowByRiver"),
             meanBelowLowFlowThreshFlowByArea_flowExt = getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "belowLowFlowThreshFlowByArea_flowExt"),
             meanAboveHighFlowThreshFlowByArea_flowExt = getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "aboveHighFlowThreshFlowByArea_flowExt")
      ) %>% 
      ungroup()
    
    coreData <- left_join(coreData, coreDataUniqueDates, 
                          by = c("detectionDate", "river", "lagDetectionDate"))
  }
  else {
    getIntervalMean <- function(start, end, e, fun = func) {
      d <- envData$date
      meanEnv <- fun(envData[d >= start & d <= end, tolower(e)], 
                     na.rm = T)
      return(meanEnv)
    }
    
    coreDataUniqueDates <- coreData %>% 
      dplyr::select(detectionDate, lagDetectionDate) %>% 
      unique() %>% 
      group_by(detectionDate, lagDetectionDate) %>% 
      mutate(meanTemperature = getIntervalMean(detectionDate, lagDetectionDate, "Temperature"), 
             meanFlow = getIntervalMean(detectionDate, lagDetectionDate, "Flow")) %>% 
      ungroup()
    
    coreData <- left_join(coreData, coreDataUniqueDates, 
                          by = c("detectionDate", "lagDetectionDate"))
  }
  
  if (sampleFlow) {
    coreData <- coreData %>% 
      mutate(date = as.Date(detectionDate)) %>% 
      filter(enc == 1) %>% dplyr::select(sampleName, date) %>% 
      group_by(sampleName, date) %>% summarize(n = n()) %>% 
      ungroup() %>% 
      left_join(envData %>% 
                  filter(!is.na(flow)) %>% 
                  mutate(date = as.Date(date)) %>% 
                  dplyr::select(date, flow) %>% 
                  rename(flowForP = flow) %>% 
                  unique(), by = c("date")) %>% 
      group_by(sampleName) %>% summarize(flowForP = sum(flowForP * n)/(sum(n))) %>% 
      ungroup() %>% 
      right_join(coreData, by = "sampleName")
  }
  # names(coreData)[which(names(coreData) == "meanTemperature")] <- paste0(funName, "Temperature")
  # names(coreData)[which(names(coreData) == "meanFlow")] <- paste0(funName,  "Flow")
  # names(coreData)[which(names(coreData) == "meanFlowByRiver")] <- paste0(funName,  "FlowByRiver")
  # names(coreData)[which(names(coreData) == "meanFlowByArea_flowExt")] <- paste0(funName,  "FlowByArea_flowExt")
  # names(coreData)[which(names(coreData) == "meanFlowByArea_ByRiver")] <- paste0(funName,  "FlowByArea_ByRiver")
  # names(coreData)[which(names(coreData) == "sdFlow")] <- paste0("sd",  "Flow")
  # names(coreData)[which(names(coreData) == "sdFlowByRiver")] <- paste0("sd",  "FlowByRiver")
  # names(coreData)[which(names(coreData) == "sdFlowByArea_flowExt")] <- paste0("sd",  "FlowByArea_flowExt")
  # names(coreData)[which(names(coreData) == "sdFlowByArea_ByRiver")] <- paste0("sd",  "FlowByArea_ByRiver")
  return(coreData)
}

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
  
  meanBelowLowFlowThreshFlowByRiverFill <- 0
  meanBelowLowFlowThreshFlowByRiverWide <- getEHDataWide_AIS(d, cols, ops, vals, "meanBelowLowFlowThreshFlowByRiver", maxOccasionValue, valuesFill = meanBelowLowFlowThreshFlowByRiverFill)
  meanBelowLowFlowThreshFlowByRiverMatrix <- as.matrix(meanBelowLowFlowThreshFlowByRiverWide %>% dplyr::select(-tag), nrow = nrow(meanBelowLowFlowThreshFlowByRiverWide), ncol = ncol(meanBelowLowFlowThreshFlowByRiverWide) - 1)
  meanBelowLowFlowThreshFlowByRiverMatrix <- ifelse(is.finite(meanBelowLowFlowThreshFlowByRiverMatrix), meanBelowLowFlowThreshFlowByRiverMatrix, meanBelowLowFlowThreshFlowByRiverFill)
  
  meanAboveHighFlowThreshFlowByRiverFill <- 0
  meanAboveHighFlowThreshFlowByRiverWide <- getEHDataWide_AIS(d, cols, ops, vals, "meanAboveHighFlowThreshFlowByRiver", maxOccasionValue, valuesFill = meanAboveHighFlowThreshFlowByRiverFill)
  meanAboveHighFlowThreshFlowByRiverMatrix <- as.matrix(meanAboveHighFlowThreshFlowByRiverWide %>% dplyr::select(-tag), nrow = nrow(meanAboveHighFlowThreshFlowByRiverWide), ncol = ncol(meanAboveHighFlowThreshFlowByRiverWide) - 1)
  meanAboveHighFlowThreshFlowByRiverMatrix <- ifelse(is.finite(meanAboveHighFlowThreshFlowByRiverMatrix), meanAboveHighFlowThreshFlowByRiverMatrix, meanAboveHighFlowThreshFlowByRiverFill)
  
  meanBelowLowFlowThreshFlowByArea_flowExtFill <- 0
  meanBelowLowFlowThreshFlowByArea_flowExtWide <- getEHDataWide_AIS(d, cols, ops, vals, "meanBelowLowFlowThreshFlowByArea_flowExt", maxOccasionValue, valuesFill = meanBelowLowFlowThreshFlowByArea_flowExtFill)
  meanBelowLowFlowThreshFlowByArea_flowExtMatrix <- as.matrix(meanBelowLowFlowThreshFlowByArea_flowExtWide %>% dplyr::select(-tag), nrow = nrow(meanBelowLowFlowThreshFlowByArea_flowExtWide), ncol = ncol(meanBelowLowFlowThreshFlowByArea_flowExtWide) - 1)
  meanBelowLowFlowThreshFlowByArea_flowExtMatrix <- ifelse(is.finite(meanBelowLowFlowThreshFlowByArea_flowExtMatrix), meanBelowLowFlowThreshFlowByArea_flowExtMatrix, meanBelowLowFlowThreshFlowByArea_flowExtFill)
  
  meanAboveHighFlowThreshFlowByArea_flowExtFill <- 0
  meanAboveHighFlowThreshFlowByArea_flowExtWide <- getEHDataWide_AIS(d, cols, ops, vals, "meanAboveHighFlowThreshFlowByArea_flowExt", maxOccasionValue, valuesFill = meanAboveHighFlowThreshFlowByArea_flowExtFill)
  meanAboveHighFlowThreshFlowByArea_flowExtMatrix <- as.matrix(meanAboveHighFlowThreshFlowByArea_flowExtWide %>% dplyr::select(-tag), nrow = nrow(meanAboveHighFlowThreshFlowByArea_flowExtWide), ncol = ncol(meanAboveHighFlowThreshFlowByArea_flowExtWide) - 1)
  meanAboveHighFlowThreshFlowByArea_flowExtMatrix <- ifelse(is.finite(meanAboveHighFlowThreshFlowByArea_flowExtMatrix), meanAboveHighFlowThreshFlowByArea_flowExtMatrix, meanAboveHighFlowThreshFlowByArea_flowExtFill)
  
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
  last <- ifelse(last == maxOccasionValue, last, last - 1)
  
  return(list(eh = eh, flow = flowMatrix, flowByRiver = flowByRiverMatrix, flowByArea_flowExt = flowByArea_flowExtMatrix, 
              flowByArea_ByRiver = flowByArea_ByRiverMatrix, 
              meanBelowLowFlowThreshFlowByRiver = meanBelowLowFlowThreshFlowByRiverMatrix, meanAboveHighFlowThreshFlowByRiver = meanAboveHighFlowThreshFlowByRiverMatrix,
              meanBelowLowFlowThreshFlowByArea_flowExt = meanBelowLowFlowThreshFlowByArea_flowExtMatrix, meanAboveHighFlowThreshFlowByArea_flowExt = meanAboveHighFlowThreshFlowByArea_flowExtMatrix,
              temperature = temperatureMatrix, river = riverMatrix, section = sectionRiverNMatrix,
              riverN = riverNMatrix, isYOY = isYOYMatrix, length = lengthMatrix, tags = tags, cohorts = cohorts, seasons = seasons, species = species,
              first = first, last = last, data = data))
}
