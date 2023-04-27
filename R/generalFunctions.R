# Define the facet labeller function
labelsIntYOY <- c(
  "1" = "Growth year 1",
  "2" = "Growth year 2"
)

labelsIntSeason <- c(
  "1" = "Spring",
  "2" = "Summer",
  "3" = "Autumn",
  "4" = "Winter"
)

labelsIntRiver <- c(
  "1" = "West brook",
  "2" = "Open Large",
  "3" = "Open small",
  "4" = "Isolated small"
)

global_labellerInt_WB <- labeller(
  isYOY = labelsIntYOY,
  season = labelsIntSeason,
  river = labelsIntRiver
  #.default = label_both
)

labelsRiver <- c(
  "west brook" = "West brook",
  "wb jimmy" = "Open Large",
  "wb mitchell" = "Open small",
  "wb obear" = "Isolated small"
)

global_labellerRiverSeasonInt_WB <- labeller(
  isYOY = labelsIntYOY,
  season = labelsIntSeason,
  river = labelsRiver
  #.default = label_both
)


global_labellerIntYOYSeason <- labeller(
  isYOY = labelsIntYOY
  # season = labelsIntSeason
  #.default = label_both
)

addGG <- function(d) {
  
  dOut <- 
    d |> 
    mutate(
      riverGG = factor(
        river,
        levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
        labels = c("West Brook","Open Large","Open Small","Isolated Small"),
        ordered = T
      ),
      seasonGG = factor(
        season, 
        labels = c("Spring","Summer","Autumn","Winter"), 
        ordered = T
      ),
      speciesGG = factor(
        species, 
        levels = c('bkt','bnt','ats'), 
        labels = c("Brook trout", "Brown trout", "Atlantic salmon"), 
        ordered = T
      )
    )
  return(dOut)
}

addGG_noSpecies <- function(d) {
  
  dOut <- 
    d |> 
    mutate(
      riverGG = factor(
        river,
        levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
        labels = c("West Brook","Open Large","Open Small","Isolated Small"),
        ordered = T
      ),
      seasonGG = factor(
        season, 
        labels = c("Spring","Summer","Autumn","Winter"), 
        ordered = T
      )
    )
  return(dOut)
}

# This function was in both getDataElectro.R and getDataCMR_targets.R
# Now here only so no need to update func in both R files.
addEnvironmental3 <- function(coreData, sampleFlow = F, funName = "mean") {
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
      #envDataIn |> 
      dplyr::filter(date <= max(coreData$detectionDate), 
                    date >= min(coreData$detectionDate)) |> 
      data.frame()
    
    # Flow duration curve env data
    envDataWB_fdcThresh <- tar_read(envDataWB_fdcThresh_target)

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
    #arrange(ageInSamples) %>% # this step is in this function in getDataCMR-targets
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
    getIntervalMeanFDC <- function(start, end, r, e) {
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
  
    coreDataUniqueDates <- coreData |>
      dplyr::select(river, detectionDate, lagDetectionDate)|>
      unique()|>
      group_by(river, detectionDate, lagDetectionDate)|>
      mutate(
        meanTemperature = getIntervalMean(detectionDate, lagDetectionDate, river, "Temperature") ,
        meanFlow =        getIntervalMean(detectionDate, lagDetectionDate, river, "Flow"),
        meanFlowByRiver = getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByRiver"),
        meanFlowByArea_flowExt = getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByArea_flowExt"),
        meanFlowByArea_ByRiver = getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByArea_ByRiver"),
        
        sdFlow =        getIntervalMean(detectionDate, lagDetectionDate, river, "Flow", get("sd")),
        sdFlowByRiver = getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByRiver", get("sd")),
        sdFlowByArea_flowExt = getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByArea_flowExt", get("sd")),
        sdFlowByArea_ByRiver = getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByArea_ByRiver", get("sd")),
        
        propBelowLoFlowThreshByRiver =        getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "belowLoFlowThreshByRiver"),
        propAboveHiFlowThreshByRiver =        getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "aboveHiFlowThreshByRiver"),
        propBelowLoFlowThreshByArea_flowExt = getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "belowLoFlowThreshByArea_flowExt"),
        propAboveHiFlowThreshByArea_flowExt = getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "aboveHiFlowThreshByArea_flowExt")
      )|>
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
  return(coreData)
}

# # Version from getDataCMR_targets.R
# addEnvironmental <- function(coreData, sampleFlow = F, funName = "mean") {
#   reconnect()
#   
#   func <- get(funName)
#   whichDrainage <- "west"
#   if (all(!unique(coreData$river) %in% c("west brook", "wb jimmy", 
#                                          "wb mitchell", "wb obear"))) {
#     whichDrainage <- "stanley"
#   }
#   if (whichDrainage == "west") {
#     envData <- 
#       tar_read(envDataWB_target) |> 
#       dplyr::filter(date <= max(coreData$detectionDate), 
#                     date >= min(coreData$detectionDate)) |> 
#       data.frame()
#     
#     envDataWB_fdcThresh <- tar_read(envDataWB_fdcThresh_target)
#     
#     # tbl(conDplyr, "data_daily_temperature") %>% 
#     # collect(n = Inf) %>% 
#     # full_join(tbl(conDplyr, "data_flow_extension") %>% 
#     #             collect(n = Inf), by = c("river", "date")) %>% 
#     # dplyr::select(-source) %>% 
#     # dplyr::filter(date <= max(coreData$detectionDate), 
#     #               date >= min(coreData$detectionDate)) %>% 
#     # rename(temperature = daily_mean_temp, flow = qPredicted) %>% 
#     # data.frame() %>%
#     # mutate(dateDate = date(date)) %>%
#     # left_join(tar_read(flowByRiver_target))
#   }
#   else {
#     envData <- tbl(conDplyr, "stanley_environmental") %>% 
#       filter(section == 11) %>% 
#       dplyr::select(datetime, temperature, depth) %>% 
#       collect(n = Inf) %>% 
#       rename(flow = depth, date = datetime) %>% 
#       data.frame()
#     warning("Depth was inserted into flow column because that is what is available in Stanley")
#   }
#   
#   coreData <- coreData %>% 
#     group_by(tag) %>% 
#     arrange(ageInSamples) %>%
#     mutate(lagDetectionDate = lead(detectionDate)) %>% 
#     ungroup()
#   
#   if (whichDrainage == "west") {
#     getIntervalMean <- function(start, end, r, e, fun = func) {
#       d <- envData$date
#       if (e == "Temperature") {
#         envCol <- "temperature"
#         if (is.na(r)) 
#           meanVar <- fun(envData[d >= start & d <= end, envCol], na.rm = T)
#         if (!is.na(r)) 
#           meanVar <- fun(envData[d >= start & d <= end & envData$river == r, envCol], na.rm = T)
#       }
#       # will need to make this river-specific
#       if (e == "Flow") {
#         envCol <- "flow"
#         meanVar <- fun(envData[d >= start & d <= end, envCol], na.rm = T)
#       }
#       if (e == "FlowByRiver") {
#         envCol <- "flowByRiver"
#         if (is.na(r)) 
#           meanVar <- fun(envData[d >= start & d <= end, envCol], na.rm = T)
#         if (!is.na(r)) 
#           meanVar <- fun(envData[d >= start & d <= end & envData$river == r, envCol], na.rm = T)
#       }
#       if (e == "FlowByArea_flowExt") {
#         envCol <- "flowByArea_flowExt"
#         if (is.na(r)) 
#           meanVar <- fun(envData[d >= start & d <= end, envCol], na.rm = T)
#         if (!is.na(r)) 
#           meanVar <- fun(envData[d >= start & d <= end & envData$river == r, envCol], na.rm = T)
#       }
#       if (e == "FlowByArea_ByRiver") {
#         envCol <- "flowByArea_ByRiver"
#         if (is.na(r)) 
#           meanVar <- fun(envData[d >= start & d <= end, envCol], na.rm = T)
#         if (!is.na(r)) 
#           meanVar <- fun(envData[d >= start & d <= end & envData$river == r, envCol], na.rm = T)
#       }
#       return(meanVar)
#     }
#     #######Function for flow duration curve data
#     getIntervalMeanFDC <- function(start, end, r, e, fun = func) {
#       d <- envDataWB_fdcThresh$date
#       if(is.na(r)) {
#         dat <- envDataWB_fdcThresh[d >= start & d <= end, e]
#         prop <- sum(dat, na.rm = TRUE) / sum(!is.na(dat))
#       } else {
#         dat <- envDataWB_fdcThresh[d >= start & d <= end & envDataWB_fdcThresh$river == r, e]
#         prop <- sum(dat, na.rm = TRUE) / sum(!is.na(dat))
#       }
#       return(prop)
#     }
#     
#     coreDataUniqueDates <- coreData %>% 
#       dplyr::select(river, detectionDate, lagDetectionDate) %>% 
#       unique() %>% 
#       group_by(river, detectionDate, lagDetectionDate) %>% 
#       mutate(meanTemperature = getIntervalMean(detectionDate, lagDetectionDate, river, "Temperature"), 
#              meanFlow =        getIntervalMean(detectionDate, lagDetectionDate, river, "Flow"),
#              meanFlowByRiver = getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByRiver"),
#              meanFlowByArea_flowExt  = getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByArea_flowExt"),
#              meanFlowByArea_ByRiver  = getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByArea_ByRiver"),
#              sdFlow =          getIntervalMean(detectionDate, lagDetectionDate, river, "Flow", get("sd")),
#              sdFlowByRiver =   getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByRiver", get("sd")),
#              sdFlowByArea_flowExt  =   getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByArea_flowExt", get("sd")),
#              sdFlowByArea_ByRiver  =   getIntervalMean(detectionDate, lagDetectionDate, river, "FlowByArea_ByRiver", get("sd")),
#              
#              # meanBelowLowFlowThreshFlowByRiver = getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "belowLowFlowThreshFlowByRiver"),
#              # meanAboveHighFlowThreshFlowByRiver = getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "aboveHighFlowThreshFlowByRiver"),
#              # meanBelowLowFlowThreshFlowByArea_flowExt = getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "belowLowFlowThreshFlowByArea_flowExt"),
#              # meanAboveHighFlowThreshFlowByArea_flowExt = getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "aboveHighFlowThreshFlowByArea_flowExt")
#              
#              propBelowLoFlowThreshByRiver =        getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "belowLoFlowThreshByRiver"),
#              propAboveHiFlowThreshByRiver =        getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "aboveHiFlowThreshByRiver"),
#              propBelowLoFlowThreshByArea_flowExt = getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "belowLoFlowThreshByArea_flowExt"),
#              propAboveHiFlowThreshByArea_flowExt = getIntervalMeanFDC(detectionDate, lagDetectionDate, river, "aboveHiFlowThreshByArea_flowExt")
#       ) %>% 
#       ungroup()
#     
#     coreData <- left_join(coreData, coreDataUniqueDates, 
#                           by = c("detectionDate", "river", "lagDetectionDate"))
#   }
#   else {
#     getIntervalMean <- function(start, end, e, fun = func) {
#       d <- envData$date
#       meanEnv <- fun(envData[d >= start & d <= end, tolower(e)], 
#                      na.rm = T)
#       return(meanEnv)
#     }
#     
#     coreDataUniqueDates <- coreData %>% 
#       dplyr::select(detectionDate, lagDetectionDate) %>% 
#       unique() %>% 
#       group_by(detectionDate, lagDetectionDate) %>% 
#       mutate(meanTemperature = getIntervalMean(detectionDate, lagDetectionDate, "Temperature"), 
#              meanFlow = getIntervalMean(detectionDate, lagDetectionDate, "Flow")) %>% 
#       ungroup()
#     
#     coreData <- left_join(coreData, coreDataUniqueDates, 
#                           by = c("detectionDate", "lagDetectionDate"))
#   }
#   
#   if (sampleFlow) {
#     coreData <- coreData %>% 
#       mutate(date = as.Date(detectionDate)) %>% 
#       filter(enc == 1) %>% dplyr::select(sampleName, date) %>% 
#       group_by(sampleName, date) %>% summarize(n = n()) %>% 
#       ungroup() %>% 
#       left_join(envData %>% 
#                   filter(!is.na(flow)) %>% 
#                   mutate(date = as.Date(date)) %>% 
#                   dplyr::select(date, flow) %>% 
#                   rename(flowForP = flow) %>% 
#                   unique(), by = c("date")) %>% 
#       group_by(sampleName) %>% summarize(flowForP = sum(flowForP * n)/(sum(n))) %>% 
#       ungroup() %>% 
#       right_join(coreData, by = "sampleName")
#   }
#   # names(coreData)[which(names(coreData) == "meanTemperature")] <- paste0(funName, "Temperature")
#   # names(coreData)[which(names(coreData) == "meanFlow")] <- paste0(funName,  "Flow")
#   # names(coreData)[which(names(coreData) == "meanFlowByRiver")] <- paste0(funName,  "FlowByRiver")
#   # names(coreData)[which(names(coreData) == "meanFlowByArea_flowExt")] <- paste0(funName,  "FlowByArea_flowExt")
#   # names(coreData)[which(names(coreData) == "meanFlowByArea_ByRiver")] <- paste0(funName,  "FlowByArea_ByRiver")
#   # names(coreData)[which(names(coreData) == "sdFlow")] <- paste0("sd",  "Flow")
#   # names(coreData)[which(names(coreData) == "sdFlowByRiver")] <- paste0("sd",  "FlowByRiver")
#   # names(coreData)[which(names(coreData) == "sdFlowByArea_flowExt")] <- paste0("sd",  "FlowByArea_flowExt")
#   # names(coreData)[which(names(coreData) == "sdFlowByArea_ByRiver")] <- paste0("sd",  "FlowByArea_ByRiver")
#   return(coreData)
# }
# 
# 
