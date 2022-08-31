tar_option_set(packages = c("tidyverse", "lubridate", "getWBData", "getPrepareWBData"))

# default values for createCoreData()
# function (sampleType = "electrofishing", baseColumns = T, 
#    columnsToAdd = NULL, includeUntagged = F, whichDrainage = "west") 

drainage <- 'west'
reconnect() 

getElectroData_target <-
  tar_plan(
    cdWB_electro0_target = createCoreData(
      sampleType = "electrofishing",  #"stationaryAntenna","portableAntenna"
      columnsToAdd = c("sampleNumber",
                       "river",
                       "survey",
                       "pass",
                       "observedLength",
                       "observedWeight",
                       "comments"),
      includeUntagged = TRUE,
      whichDrainage = "west"
    ) %>%
      addTagProperties(
        columnsToAdd = c("cohort",
                         "species",
                         "dateEmigrated",
                         "sex",
                         "species"
        )
      ) %>%
      dplyr::filter(species %in% c( "bkt","bnt","ats"),
                    area %in% c("trib","inside","below","above"),
                    !is.na(sampleNumber)) %>%
      addSampleProperties() %>%
      addEnvironmental2() %>%
      mutate(riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                                   labels = c("West Brook","WB Jimmy","WB Mitchell","WB OBear"), ordered = T)
      ),
    
    # functions in getPrepareWBData library
    cdWB_electro_target = cdWB_electro0_target %>%
      cleanData(drainage) %>%
      mergeSites(drainage) %>%
      addNPasses(drainage) %>%
      mutate(drainage = drainage)
    
  )  


#####################################
## getData functions 
#####################################
# with targets, need 'reconnect()' within each function that accesses the DB
# use this until we add reconnect() to the addEnvironmental package
addEnvironmental2 <- function (coreData, sampleFlow = F, funName = "mean") {
  reconnect()
  
  func <- get(funName)
  whichDrainage <- "west"
  if (all(!unique(coreData$river) %in% c("west brook", "wb jimmy", 
                                         "wb mitchell", "wb obear"))) {
    whichDrainage <- "stanley"
  }
  if (whichDrainage == "west") {
    envData <- tbl(conDplyr, "data_daily_temperature") %>% 
      collect(n = Inf) %>% full_join(tbl(conDplyr, "data_flow_extension") %>% 
                                       collect(n = Inf), by = c("river", "date")) %>% select(-source) %>% 
      dplyr::filter(date <= max(coreData$detectionDate), 
                    date >= min(coreData$detectionDate)) %>% rename(temperature = daily_mean_temp, 
                                                                    flow = qPredicted) %>% data.frame()
  }
  else {
    envData <- tbl(conDplyr, "stanley_environmental") %>% 
      filter(section == 11) %>% select(datetime, temperature, 
                                       depth) %>% collect(n = Inf) %>% rename(flow = depth, 
                                                                              date = datetime) %>% data.frame()
    warning("Depth was inserted into flow column because that is what is available in Stanley")
  }
  coreData <- coreData %>% group_by(tag) %>% mutate(lagDetectionDate = lead(detectionDate)) %>% 
    ungroup()
  if (whichDrainage == "west") {
    getIntervalMean <- function(start, end, r, e, fun = func) {
      d <- envData$date
      if (e == "Temperature") {
        envCol <- "temperature"
        if (is.na(r)) 
          meanTemp <- fun(envData[d >= start & d <= end, 
                                  envCol], na.rm = T)
        if (!is.na(r)) 
          meanTemp <- fun(envData[d >= start & d <= end & 
                                    envData$river == r, envCol], na.rm = T)
      }
      if (e == "Flow") {
        envCol <- "flow"
        meanTemp <- fun(envData[d >= start & d <= end, 
                                envCol], na.rm = T)
      }
      return(meanTemp)
    }
    coreDataUniqueDates <- coreData %>% select(river, detectionDate, 
                                               lagDetectionDate) %>% unique() %>% group_by(river, 
                                                                                           detectionDate, lagDetectionDate) %>% mutate(meanTemperature = getIntervalMean(detectionDate, 
                                                                                                                                                                         lagDetectionDate, river, "Temperature"), meanFlow = getIntervalMean(detectionDate, 
                                                                                                                                                                                                                                             lagDetectionDate, river, "Flow")) %>% ungroup()
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
    coreDataUniqueDates <- coreData %>% select(detectionDate, 
                                               lagDetectionDate) %>% unique() %>% group_by(detectionDate, 
                                                                                           lagDetectionDate) %>% mutate(meanTemperature = getIntervalMean(detectionDate, 
                                                                                                                                                          lagDetectionDate, "Temperature"), meanFlow = getIntervalMean(detectionDate, 
                                                                                                                                                                                                                       lagDetectionDate, "Flow")) %>% ungroup()
    coreData <- left_join(coreData, coreDataUniqueDates, 
                          by = c("detectionDate", "lagDetectionDate"))
  }
  if (sampleFlow) {
    coreData <- coreData %>% mutate(date = as.Date(detectionDate)) %>% 
      filter(enc == 1) %>% select(sampleName, date) %>% 
      group_by(sampleName, date) %>% summarize(n = n()) %>% 
      ungroup() %>% left_join(envData %>% filter(!is.na(flow)) %>% 
                                mutate(date = as.Date(date)) %>% select(date, flow) %>% 
                                rename(flowForP = flow) %>% unique(), by = c("date")) %>% 
      group_by(sampleName) %>% summarize(flowForP = sum(flowForP * 
                                                          n)/(sum(n))) %>% ungroup() %>% right_join(coreData, 
                                                                                                    by = "sampleName")
  }
  names(coreData)[which(names(coreData) == "meanTemperature")] <- paste0(funName, 
                                                                         "Temperature")
  names(coreData)[which(names(coreData) == "meanFlow")] <- paste0(funName, 
                                                                  "Flow")
  return(coreData)
}
