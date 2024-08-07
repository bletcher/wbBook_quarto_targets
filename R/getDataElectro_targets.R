tar_option_set(packages = c("tidyverse", "lubridate", "getWBData", "getPrepareWBData"))

# default values for createCoreData()
# function (sampleType = "electrofishing", baseColumns = T, 
#    columnsToAdd = NULL, includeUntagged = F, whichDrainage = "west") 

drainage <- 'west' # should move this to a target
reconnect() 

getElectroData_target <-
  tar_plan(
    cdWB_electro0_target = createCoreData(
      sampleType = "electrofishing",  #"stationaryAntenna","portableAntenna"
      columnsToAdd = c(
        "sampleNumber",
        "river",
        "survey",
        "pass",
        "observedLength",
        "observedWeight",
        "comments"
      ),
      includeUntagged = TRUE,
      whichDrainage = "west"
    ) %>%
    addTagProperties(
      columnsToAdd = 
        c(
          "cohort",
          "species",
          "dateEmigrated",
          "sex"
        )
    ) %>%
    dplyr::filter(
      species %in% c( "bkt","bnt","ats"),
      area %in% c("trib","inside","below","above"),
      !is.na(sampleNumber)
    ) %>%
    addSampleProperties() %>%
    addEnvironmental3(envDataWBIn = envDataWB_target, envDataWB_fdcThreshIn = envDataWB_fdcThresh_target) %>%
    mutate(
      riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell', "wb obear"),
                                   labels = c("West Brook", "WB Jimmy", "WB Mitchell", "WB OBear"), ordered = T),
      # # variables to allow merging with wanding data
      readerId = NA,
      sectionN = as.numeric(section),
      aliveOrDead = "alive",
      instance = NA,
      quarter = NA,
      leftOrRight = NA,
      habitat = NA,
      cover = NA,
      justification = NA,
      sectionWQuarter = NA,
      j = NA, # not sure what this is - julian date
      
      date = date(detectionDate),
      dummy = 1
    ),
    
    # functions in getPrepareWBData library
    cdWB_electro_target = cdWB_electro0_target %>%
      cleanData2(drainage) %>%
      mergeSites(drainage) %>%
      addNPasses(drainage) %>%
      mutate(drainage = drainage) %>%
      addSizeIndGrowthWeight2() %>% # adapted below so
      addCF() %>%
      addCountPerIndividual() |> 
      addSeasonFT_zScore() |> 
      addSeasonRiverFT_zScore(),
      # 
    
    medianWinterSampleDate_target = cdWB_electro_target |> 
      filter(season == 4) |>
      group_by(year) |> 
      summarize(start0 = median(date)) |> 
      mutate(start = setYday(start0, medianWinterYday_target$medianYday)) |> 
      ungroup() ,
    
    medianWinterYday_target = cdWB_electro_target |> 
      filter(season == 4) |> 
      summarize(medianYday = median(yday(date))),
    
    missingWinterYears_target = 
      tibble(
        river = rep(c('west brook', 'wb jimmy', 'wb mitchell', 'wb obear'), each = 4),
        year = rep(c(2005, 2007, 2012, 2014), 4),
        season = 4
      ) |> 
      left_join(medianWinterSampleDate_target),
    
    medianDates_target = 
      medDate <- cdWB_electro_target |> 
        group_by(river, year, season) |> 
        summarize(start = median(date)) |> 
        ungroup() |> 
        add_row(missingWinterYears_target |> dplyr::select(-start0)) |> 
        arrange(river, start) |> 
        group_by(river) |> 
        mutate(end = lead(start)) |> 
        filter(!is.na(end))
  )  


#####################################
## getData functions 
#####################################
# with targets, need 'reconnect()' within each function that accesses the DB
# use this until we add reconnect() to the addEnvironmental package

#addEnvironmental3() is in 'generalFunctions.R'


# needed to set yday in mutate()
setYday <- function(d, dateYday){
  yday(d) = dateYday
  return(d)
}


addSizeIndGrowthWeight2 <- function(dIn) {
  
  d <- dIn |> filter(sampleNumberDiff == 1)
  
  lmMod <- function(dd){
    lm(log(grWeight) ~ log(observedWeight), data = dd)
  }
  
  weight_grMod <- d |> 
    filter(grWeight > 0) |> # to avoid log() problems
    group_by(species, season, river) |> 
    nest() |> 
    mutate(model = map(data, lmMod))
  
  weight_grMod_GT <- weight_grMod  |>  
    mutate(glanced = map(model, broom::glance),
           tidied = map(model, broom::tidy)) 
  
  weight_grMod_GT_slope <- weight_grMod_GT  |>  
    unnest(tidied) |> 
    filter(term == "log(observedWeight)") |> 
    dplyr::select(species, river, season, estimate) |> 
    rename(wGR_Slope = estimate) |> 
    mutate(wGR_Slope = ifelse(species == 'ats' & river == 'wb jimmy' & season == 2, -0.572, wGR_Slope), # set values = WB values - limited data
           wGR_Slope = ifelse(species == 'ats' & river == 'wb jimmy' & season == 4, -0.139, wGR_Slope))
  
  dOut <- dIn |> 
    left_join(weight_grMod_GT_slope) |> 
    mutate(grWeightS = (lagObservedWeight ^ wGR_Slope -
                          observedWeight ^ wGR_Slope) /
             (as.numeric(lagDetectionDate - detectionDate) * wGR_Slope)
    )
  
  return(dOut)
  
}

addCF <- function(d){
  cfModel = lm(log10(observedWeight) ~ log10(observedLength) * factor(species), data = d)
  
  cf_coefs = tibble(
    species = c("ats", "bkt", "bnt"),
    lw_intercept = c(
      coef(cfModel)[["(Intercept)"]],
      coef(cfModel)[["(Intercept)"]] + coef(cfModel)[["factor(species)bkt"]],
      coef(cfModel)[["(Intercept)"]] + coef(cfModel)[["factor(species)bnt"]]
    ),
    lw_slope = c(
      coef(cfModel)[["log10(observedLength)"]],
      coef(cfModel)[["log10(observedLength)"]] + coef(cfModel)[["log10(observedLength):factor(species)bkt"]],
      coef(cfModel)[["log10(observedLength)"]] + coef(cfModel)[["log10(observedLength):factor(species)bnt"]]
    )
  )

 dOut <- d |> 
   left_join(cf_coefs) |> 
   mutate(
     cf = (10000 * observedWeight) / observedLength ^ 3,
     relCF = (10000 * observedWeight) / observedLength ^ lw_slope
   )
 
 return(dOut)
}


# This function does not include flowByRiver (addEnvironmental3 does)
# addEnvironmental2 <- function (coreData, sampleFlow = F, funName = "mean") {
#   reconnect()
#   
#   func <- get(funName)
#   whichDrainage <- "west"
#   if (all(!unique(coreData$river) %in% c("west brook", "wb jimmy", 
#                                          "wb mitchell", "wb obear"))) {
#     whichDrainage <- "stanley"
#   }
#   if (whichDrainage == "west") {
#     envData <- tbl(conDplyr, "data_daily_temperature") %>% 
#       collect(n = Inf) %>% full_join(tbl(conDplyr, "data_flow_extension") %>% 
#                                        collect(n = Inf), by = c("river", "date")) %>% dplyr::select(-source) %>% 
#       dplyr::filter(date <= max(coreData$detectionDate), 
#                     date >= min(coreData$detectionDate)) %>% rename(temperature = daily_mean_temp, 
#                                                                     flow = qPredicted) %>% data.frame()
#   }
#   else {
#     envData <- tbl(conDplyr, "stanley_environmental") %>% 
#       filter(section == 11) %>% dplyr::select(datetime, temperature, 
#                                        depth) %>% collect(n = Inf) %>% rename(flow = depth, 
#                                                                               date = datetime) %>% data.frame()
#     warning("Depth was inserted into flow column because that is what is available in Stanley")
#   }
#   coreData <- coreData %>% group_by(tag) %>% mutate(lagDetectionDate = lead(detectionDate)) %>% 
#     ungroup()
#   if (whichDrainage == "west") {
#     getIntervalMean <- function(start, end, r, e, fun = func) {
#       d <- envData$date
#       if (e == "Temperature") {
#         envCol <- "temperature"
#         if (is.na(r)) 
#           meanTemp <- fun(envData[d >= start & d <= end, 
#                                   envCol], na.rm = T)
#         if (!is.na(r)) 
#           meanTemp <- fun(envData[d >= start & d <= end & 
#                                     envData$river == r, envCol], na.rm = T)
#       }
#       if (e == "Flow") {
#         envCol <- "flow"
#         meanTemp <- fun(envData[d >= start & d <= end, 
#                                 envCol], na.rm = T)
#       }
#       return(meanTemp)
#     }
#     coreDataUniqueDates <- coreData %>% dplyr::select(river, detectionDate, 
#                                                lagDetectionDate) %>% unique() %>% group_by(river, 
#                                                                                            detectionDate, lagDetectionDate) %>% mutate(meanTemperature = getIntervalMean(detectionDate, 
#                                                                                                                                                                          lagDetectionDate, river, "Temperature"), meanFlow = getIntervalMean(detectionDate, 
#                                                                                                                                                                                                                                              lagDetectionDate, river, "Flow")) %>% ungroup()
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
#     coreDataUniqueDates <- coreData %>% dplyr::select(detectionDate, 
#                                                lagDetectionDate) %>% unique() %>% group_by(detectionDate, 
#                                                                                            lagDetectionDate) %>% mutate(meanTemperature = getIntervalMean(detectionDate, 
#                                                                                                                                                           lagDetectionDate, "Temperature"), meanFlow = getIntervalMean(detectionDate, 
#                                                                                                                                                                                                                        lagDetectionDate, "Flow")) %>% ungroup()
#     coreData <- left_join(coreData, coreDataUniqueDates, 
#                           by = c("detectionDate", "lagDetectionDate"))
#   }
#   if (sampleFlow) {
#     coreData <- coreData %>% mutate(date = as.Date(detectionDate)) %>% 
#       filter(enc == 1) %>% dplyr::select(sampleName, date) %>% 
#       group_by(sampleName, date) %>% summarize(n = n()) %>% 
#       ungroup() %>% left_join(envData %>% filter(!is.na(flow)) %>% 
#                                 mutate(date = as.Date(date)) %>% dplyr::select(date, flow) %>% 
#                                 rename(flowForP = flow) %>% unique(), by = c("date")) %>% 
#       group_by(sampleName) %>% summarize(flowForP = sum(flowForP * 
#                                                           n)/(sum(n))) %>% ungroup() %>% right_join(coreData, 
#                                                                                                     by = "sampleName")
#   }
#   names(coreData)[which(names(coreData) == "meanTemperature")] <- paste0(funName, 
#                                                                          "Temperature")
#   names(coreData)[which(names(coreData) == "meanFlow")] <- paste0(funName, 
#                                                                   "Flow")
#   return(coreData)
# }

# cleanData function is from getPrepareWBData package. Editing here temporarily
cleanData2 <- function (d, drainageIn) {
  d$sectionOriginal <- d$section
  d$section <- as.numeric(d$section)
  if (drainageIn == "west") {
    maxSectionNum <- 47
    d$riverOrdered <- factor(d$river, levels = c("west brook", "wb jimmy", "wb mitchell", "wb obear"), 
                             labels = c("west brook", "wb jimmy", "wb mitchell", "wb obear"), 
                             ordered = T)
    minYear = min(d$year)
  }
  else if (drainageIn == "stanley") {
    maxSectionNum <- 50
    d$riverOrdered <- factor(d$river, levels = c("mainstem", "west", "east"), 
                             labels = c("mainstem", "west", "east"), 
                             ordered = T)
    minYear = min(d$year)
  }
  d$inside <- ifelse(d$section %in% 1:maxSectionNum | d$survey == "stationaryAntenna", 
                     T, F)
  d$year <- year(d$detectionDate)
  d$yday <- yday(d$detectionDate)
  d$ageInSamples <- (d$year - d$cohort) * 4 + (d$season - 2)
  d$isYOY <- ifelse(d$ageInSamples <= 3, TRUE, FALSE)
  dUntagged <- d %>% 
    filter(is.na(tag)) %>% 
    mutate(minSample = min(sampleNumber), 
           maxSample = max(sampleNumber), 
           minYear = minYear, 
           moveDir = 0, 
           sampleInterval = 0)
  
  d <- d %>% 
    filter(!is.na(tag)) %>% 
    group_by(tag) %>% 
    mutate(lagSection = lead(section), 
           distMoved = section - lagSection, 
           observedWeight = ifelse(observedWeight <= -9999, NA, observedWeight), 
           lagObservedWeight = lead(observedWeight), 
           lagObservedLength = lead(observedLength), 
           lagSampleNumber = lead(sampleNumber), 
           sampleNumberDiff = lagSampleNumber - sampleNumber,
           grWeight = (log(lagObservedWeight) - log(observedWeight))/as.numeric((lagDetectionDate - detectionDate)),
           #grWeight = ifelse(grWeight0 <= 0, 0.0001, grWeight0),
           grLength = (lagObservedLength - observedLength)/as.numeric((lagDetectionDate - detectionDate)), 
           minSample = min(sampleNumber), 
           maxSample = max(sampleNumber), 
           minYear = minYear) %>% 
    ungroup()
  
  d$moveDir <- ifelse(d$section == d$lagSection, 
                      0, 
                      ifelse(d$section > d$lagSection, 
                             1, 
                             -1)
                      )
  d$sampleInterval <- as.numeric(d$lagDetectionDate - d$detectionDate)
  d <- bind_rows(d, dUntagged)
  return(d)
}

addSizeIndGrowthWeight <- function(dIn) {
  
  d <- dIn |> filter(sampleNumberDiff == 1)
  
  lmMod <- function(dd){
    lm(log(grWeight) ~ log(observedWeight), data = dd)
  }
  
  weight_grMod <- d |> 
    filter(grWeight > 0) |> # to avoid log() problems
    group_by(species, season, river) |> 
    nest() |> 
    mutate(model = map(data, lmMod))
  
  weight_grMod_GT <- weight_grMod  |>  
    mutate(glanced = map(model, broom::glance),
           tidied = map(model, broom::tidy)) 
  
  weight_grMod_GT_slope <- weight_grMod_GT  |>  
    unnest(tidied) |> 
    filter(term == "log(observedWeight)") |> 
    dplyr::select(species, river, season, estimate) |> 
    rename(wGR_Slope = estimate) |> 
    mutate(wGR_Slope = ifelse(species == 'ats' & river == 'wb jimmy' & season == 2, -0.572, wGR_Slope), # set values = WB values - limited data
           wGR_Slope = ifelse(species == 'ats' & river == 'wb jimmy' & season == 4, -0.139, wGR_Slope))
  
  dOut <- d |> 
    left_join(weight_grMod_GT_slope) |> 
    mutate(grWeightS = (lagObservedWeight ^ wGR_Slope -
                          observedWeight ^ wGR_Slope) /
             (as.numeric(lagDetectionDate - detectionDate) * wGR_Slope)
    )
  
  return(dOut)
  
}

addCountPerIndividual <- function(dIn){
  d <- dIn |> 
    ungroup() |> 
    group_by(tag) |> 
    summarize(nPerInd = n())
  
  left_join(dIn, d)
}

addSeasonFT_zScore <- function(d){
  seasonTF <- d |> 
    group_by(season) |> 
    summarize(meanT = mean(meanTemperature, na.rm = TRUE),
              sdT = sd(meanTemperature, na.rm = TRUE),
              meanF = mean(meanFlow, na.rm = TRUE),
              sdF = sd(meanFlow, na.rm = TRUE),
              meanFByRiver = mean(meanFlowByRiver, na.rm = TRUE),
              sdFByRiver = sd(meanFlowByRiver, na.rm = TRUE),
              meanFByArea_flowExt = mean(meanFlowByArea_flowExt, na.rm = TRUE),
              sdFByArea_flowExt = sd(meanFlowByArea_flowExt, na.rm = TRUE),
              meanFByArea_ByRiver = mean(meanFlowByArea_ByRiver, na.rm = TRUE),
              sdFByArea_ByRiver = sd(meanFlowByArea_ByRiver, na.rm = TRUE)) 
  
  d <- left_join(d, seasonTF) |> 
    mutate(
      meanTemperatureScaledBySeason = (meanTemperature - meanT)/sdT, 
      meanFlowScaledBySeason = (meanFlow - meanF)/sdF,
      meanFlowByRiverScaledBySeason = (meanFlowByRiver - meanFByRiver)/sdFByRiver,
      meanFlowByArea_flowExtScaledBySeason = (meanFlowByArea_flowExt - meanFByArea_flowExt)/sdFByArea_flowExt,
      meanFlowByArea_ByRiverScaledBySeason = (meanFlowByArea_ByRiver - meanFByArea_ByRiver)/sdFByArea_ByRiver
    ) |> 
    dplyr::select(-c(meanT, sdT, meanF, sdF, meanFByRiver, sdFByRiver, meanFByArea_flowExt, meanFByArea_ByRiver, sdFByArea_flowExt, sdFByArea_ByRiver))
  
  return(d)
}

addSeasonRiverFT_zScore <- function(d){
  seasonRiverTF <- d |> 
    group_by(season, river) |> 
    summarize(meanT = mean(meanTemperature, na.rm = TRUE),
              sdT = sd(meanTemperature, na.rm = TRUE),
              meanF = mean(meanFlow, na.rm = TRUE),
              sdF = sd(meanFlow, na.rm = TRUE),
              meanFByRiver = mean(meanFlowByRiver, na.rm = TRUE),
              sdFByRiver = sd(meanFlowByRiver, na.rm = TRUE),
              meanFByArea_flowExt = mean(meanFlowByArea_flowExt, na.rm = TRUE),
              sdFByArea_flowExt = sd(meanFlowByArea_flowExt, na.rm = TRUE),
              meanFByArea_ByRiver = mean(meanFlowByArea_ByRiver, na.rm = TRUE),
              sdFByArea_ByRiver = sd(meanFlowByArea_ByRiver, na.rm = TRUE)) 
  
  d <- left_join(d, seasonRiverTF) |> 
    mutate(
      meanTemperatureScaledBySeasonRiver = (meanTemperature - meanT)/sdT, 
      meanFlowScaledBySeasonRiver = (meanFlow - meanF)/sdF,
      meanFlowByRiverScaledBySeasonRiver = (meanFlowByRiver - meanFByRiver)/sdFByRiver,
      meanFlowByArea_flowExtScaledBySeasonRiver = (meanFlowByArea_flowExt - meanFByArea_flowExt)/sdFByArea_flowExt,
      meanFlowByArea_ByRiverScaledBySeasonRiver = (meanFlowByArea_ByRiver - meanFByArea_ByRiver)/sdFByArea_ByRiver
    ) |> 
    dplyr::select(-c(meanT, sdT, meanF, sdF, meanFByRiver, sdFByRiver, meanFByArea_flowExt, meanFByArea_ByRiver, sdFByArea_flowExt, sdFByArea_ByRiver))
  
  return(d)
}


