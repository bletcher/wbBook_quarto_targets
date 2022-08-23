tar_option_set(packages = c("tidyverse", "lubridate", "getWBData", "targets"))

selectedVariables <- c("tag", "species", "river", "detectionDate", "sampleNumber", "n", 
                       "proportionSampled", "observedLength", "observedWeight", "area", 
                       "section", "season", "isYOY")

cdWB_electro <- tar_read(cdWB_electro_target)
envDataWB <- tar_read(envDataWB_target)

spawn_month <- "11" # spawning
spawn_day <- "15"
emerge_month <- "03" # emergence
emerge_day <- "01"

minYear <- 2000
maxYear <- 2015

modelYOY_target <-
  tar_plan(
    firstObs_noTag_target = cdWB_electro %>%
      filter(is.na(tag), ageInSamples == 1) %>%
      mutate(n = 1) %>%
      dplyr:: select(all_of(selectedVariables)),
    
    firstObs_tag_target = cdWB_electro %>%
        group_by(tag) %>%
        mutate(isFirstObs = detectionDate == min(detectionDate),
               n = n()) %>%
        filter(isFirstObs, ageInSamples == 1) %>%
        dplyr::select(all_of(selectedVariables)) %>%
        ungroup(),

    firstObs0_target = add_row(firstObs_tag_target, firstObs_noTag_target) %>%
      mutate(date = as_date(detectionDate),
             yday = yday(date),
             year = year(date)
             ),
    
    # Then merge results with firstObs0 to create firstObs.
    firstObsDates_target = firstObs0_target %>% 
      distinct(date = date(detectionDate), river),
    
    firstObs_Env_target = firstObsDates_target %>%
      rowwise() %>%
      mutate(
        year = year(date),
        spawnDate = ymd(paste0(year,spawn_month,spawn_day)) - years(1),
        emergeDate = ymd(paste0(year,emerge_month,emerge_day)),
        oneMonthDate = date - days(as.integer(1 * 30.5)), #months(1), 'months gives error when prev month has 30 days and current has 31
        threeMonthDate = date - days(as.integer(3 * 30.5)),
        fiveMonthDate = date - days(as.integer(5 * 30.5)),
        spawn_emerge = list(getEnvMeans(envDataWB, river, spawnDate, emergeDate)),
        emerge_detect = list(getEnvMeans(envDataWB, river, emergeDate, date)),
        spawn_detect = list(getEnvMeans(envDataWB, river, spawnDate, date)),
        oneMonth = list(getEnvMeans(envDataWB, river, oneMonthDate, date)),
        threeMonth = list(getEnvMeans(envDataWB, river, threeMonthDate, date)),
        fiveMonth = list(getEnvMeans(envDataWB, river, fiveMonthDate, date))
      ),
    
    # merge env data into firstObs0
    firstObs_target = firstObs0_target %>%
      left_join(firstObs_Env_target),
    
    # this scales across all individuals - I think this is ok
    firstObsUnnested_target0 = firstObs_target %>% 
      unnest(cols = c(spawn_emerge, emerge_detect, spawn_detect, oneMonth, threeMonth, fiveMonth), names_sep = "_") %>%
      mutate(
        emerge_detect_sumTScaled = getScaled(emerge_detect_sumT),
        emerge_detect_sumFScaled = getScaled(emerge_detect_sumF),
        oneMonth_sumTScaled = getScaled(oneMonth_sumT),
        oneMonth_sumFScaled = getScaled(oneMonth_sumF),
        threeMonth_sumTScaled = getScaled(threeMonth_sumT),
        threeMonth_sumFScaled = getScaled(threeMonth_sumF),
        fiveMonth_sumTScaled = getScaled(fiveMonth_sumT),
        fiveMonth_sumFScaled = getScaled(fiveMonth_sumF),
        ydayScaled = getScaled(yday)
      ),
    
    countsRSY_target = firstObs_target %>%
      filter(year %in% minYear:maxYear) %>%
      group_by(river, species, year) %>%
      summarize(
        count = n(),
        meanPropSampled = mean(proportionSampled, na.rm = TRUE)
      ) %>%
      mutate(countAdj = count / meanPropSampled) %>%
      ungroup() %>%
      group_by(river, species) %>%
      mutate(meanCountRS = mean(count, na.rm = TRUE),
             sdCountRS = sd(count, na.rm = TRUE),
             countRS_Scaled = (count - meanCountRS) / sdCountRS) %>%
      ungroup(), 
    
    countsRY_target = firstObs_target %>%
      filter(year %in% minYear:maxYear) %>%
      group_by(river, year) %>%
      summarize(
        count = n(),
        meanPropSampled = mean(proportionSampled, na.rm = TRUE)
      ) %>%
      mutate(countAdj = count / meanPropSampled) %>%
      ungroup() %>%
      group_by(river) %>%
      mutate(meanCountR = mean(count, na.rm = TRUE),
             sdCountR = sd(count, na.rm = TRUE),
             countR_Scaled = (count - meanCountR) / sdCountR) %>%
      ungroup(),
    
    countsMetaY_target = firstObs_target %>%
      filter(river != "wb obear", year %in% minYear:maxYear) %>%
      group_by(year) %>%
      summarize(
        count = n(),
        meanPropSampled = mean(proportionSampled, na.rm = TRUE)
      ) %>%
      mutate(countAdj = count / meanPropSampled) %>%
      ungroup() %>%
      mutate(meanCount = mean(count, na.rm = TRUE),
             sdCount = sd(count, na.rm = TRUE),
             count_Scaled = (count - meanCount) / sdCount),
    # missing data for tribs in 2000, 2001 - may skew scaled count a bit low - should fix
    
    firstObsUnnested_target = firstObsUnnested_target0 %>%
      left_join(countsMetaY_target %>% dplyr::select(year, count_Scaled)),
    
    firstObsUnnestedWB_target = firstObsUnnested_target %>% filter(river == "west brook"),
    
    envDataWBFlow_target = envDataWB %>%
      mutate(year = year(dateDate),
             yday = yday(dateDate),
             flowNoNAs = ifelse(is.na(flow), 0, flow),
             tempNoNAs = ifelse(is.na(temperature), 0, temperature)) %>%
      filter(year %in% minYear:maxYear, 
             yday > 100, yday < 300,
             river == "west brook") %>%
      group_by(year) %>%
      mutate(cumulFlow = cumsum(flowNoNAs),
             cumulFlow01 = cumulFlow / max(cumulFlow),
             cumulTemp = cumsum(tempNoNAs)) %>%
      ungroup()
    
  )




#####################################
#### Functions
#####################################

# move to getPrepareWBData
getEnvMeans <- function(d, riverIn, start, end) { 
  out <- d %>% 
    filter(river == riverIn, dateDate >= start, dateDate <= end) %>%
    summarize(
      sumT = sum(temperature, na.rm = TRUE),
      meanT = mean(temperature, na.rm = TRUE),
      sdT = sd(temperature, na.rm = TRUE), 
      cvT = sdT/meanT,
      
      sumF = sum(flow, na.rm = TRUE),
      meanF = mean(flow, na.rm = TRUE),
      sdF = sd(flow, na.rm = TRUE),
      cvF = sdF/meanF,
      n = n()
    )
  
  #message(paste(river, start, end,tag))
  return(out)
}

getScaled <- function(d){
  (d - mean(d, na.rm = TRUE)) / sd(d, na.rm = TRUE)
}

getMeansData <- function(d, t, f) {
  means <- d %>% 
    group_by(species, year) %>% 
    summarize(meanLength = mean(observedLength, na.rm = TRUE), 
              meanEmerge_detect_sumTScaled = mean(emerge_detect_sumTScaled, na.rm = TRUE),
              meanEmerge_detect_sumFScaled = mean(emerge_detect_sumTScaled, na.rm = TRUE),
              meanTTime_sumTScaled = mean(get(t), na.rm = TRUE),
              meanFTime_sumFScaled = mean(get(f), na.rm = TRUE),
              meanYdayScaled = mean(ydayScaled, na.rm = TRUE),
              meanCount_Scaled = mean(count_Scaled, na.rm = TRUE)
    )
  return(means)
}

getMeansDataByRiver <- function(d, t, f) {
  means <- d %>% 
    group_by(species, year, river) %>% 
    summarize(meanLength = mean(observedLength, na.rm = TRUE), 
              meanEmerge_detect_sumTScaled = mean(emerge_detect_sumTScaled, na.rm = TRUE),
              meanEmerge_detect_sumFScaled = mean(emerge_detect_sumTScaled, na.rm = TRUE),
              meanTTime_sumTScaled = mean(get(t), na.rm = TRUE),
              meanFTime_sumFScaled = mean(get(f), na.rm = TRUE),
              meanYdayScaled = mean(ydayScaled, na.rm = TRUE),
              meanCount_Scaled = mean(count_Scaled, na.rm = TRUE)
    )
  return(means)
}

plotMeans <- function(means){
  out <- list()
  out[[1]] <- ggplot(means, aes(meanTTime_sumTScaled, meanLength, color = species)) +
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE)
  
  out[[2]] <- ggplot(means, aes(meanFTime_sumFScaled, meanLength, color = species)) +
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE)
  
  out[[3]] <- ggplot(means, aes(meanTTime_sumTScaled, meanFTime_sumFScaled, color = species)) +
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE)
  return(out)
}

runMeanModels <- function(means) {
  modLMMeans1 <- lm(meanLength ~ (factor(species) + meanFTime_sumFScaled + meanTTime_sumTScaled + meanYdayScaled + meanCount_Scaled), data = means)
  modLMMeans2 <- lm(meanLength ~ (factor(species) + meanFTime_sumFScaled + meanTTime_sumTScaled + meanYdayScaled + meanCount_Scaled)^2, data = means)
  modLMMeans3 <- lm(meanLength ~ (factor(species) + meanFTime_sumFScaled + meanTTime_sumTScaled + meanYdayScaled + meanCount_Scaled)^3, data = means)
  return(list(modLMMeans1, modLMMeans2, modLMMeans3))
}
