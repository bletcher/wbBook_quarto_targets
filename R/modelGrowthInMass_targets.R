tar_option_set(packages = c("tidyverse", "lubridate"))

modelGrowthInMass_target <-
  tar_plan(
    cd1_target = cdWB_electro_target |> 
      filter(sampleNumberDiff == 1, # consecutive captures only
             tag %notin% c('1bf20ff490', '1bf20ebe4e', '1c2c582218')) |> 
      mutate(
        negGrowth = grWeight < 0,
        month = month(date),
      )|> 
      addGG() |> # in generalFunctions.R
      left_join(
        firstLast_target
      ) |> 
      mutate(
        samplesBeforeLast = lastObserved - sampleNumber
      ),
    
    # cd1Wide_observedWeight_target = cd1_target |> 
    #   filter(
    #     #!(tag=="1bf18c1e3e" & sampleNumber == 40),
    #     ageInSamples > 0
    #   ) |> 
    #   dplyr::select(speciesGG, tag, cohort, firstObserved, lastObserved, ageInSamples, observedWeight) |> 
    #   pivot_wider(
    #     names_from = ageInSamples, 
    #     values_from = observedWeight, 
    #     names_sort = TRUE, 
    #     names_prefix = "AIS_"
    #   ), 
    
    cd1Wide_observedWeight_target = cd1_target |> 
      filter(
        #!(tag=="1bf18c1e3e" & sampleNumber == 40),# bad AIS because caught in Jan
          ageInSamples > 0
        ) |> 
      mutate(age = year - cohort,
             ageSeason = paste(age, season, sep = "_")) |> 
      dplyr::select(speciesGG, tag, cohort, firstObserved, lastObserved, ageSeason, observedWeight) |> 
      pivot_wider(
        names_from = ageSeason, 
        values_from = observedWeight, 
        names_sort = TRUE 
      ),
    
    cd1Wide_grWeight_target = cd1_target |> 
      filter(!(tag=="1bf18c1e3e" & sampleNumber == 40),# bad AIS because caught in Jan
             ageInSamples > 0) |> 
      mutate(age = year - cohort,
             ageSeason = paste(age, season, sep = "_")) |> 
      dplyr::select(speciesGG, tag, cohort, firstObserved, lastObserved, ageSeason, grWeight) |> 
      pivot_wider(
        names_from = ageSeason, 
        values_from = grWeight, 
        names_sort = TRUE 
      ),
    
    
    firstLast_target = cdWB_CMR0_target |> 
      dplyr::select(tag, firstObserved, lastObserved) |> 
      distinct(),
    
    negGr_beforeLast_target = cd1_target |> 
      group_by(speciesGG, riverGG, seasonGG, samplesBeforeLast, negGrowth) |> 
      summarize(
        meanGR = mean(grWeight, na.rm = TRUE),
        n = n()
      ),
    
    negGr_beforeLastMean_target = cd1_target |> 
      filter(!is.na(negGrowth)) |> 
      group_by(speciesGG, riverGG, seasonGG, negGrowth) |> 
      summarize(
        meanGR = mean(grWeight, na.rm = TRUE),
        sdGR = sd(grWeight, na.rm = TRUE),
        meanSamplesBeforeLast = mean(samplesBeforeLast, na.rm = TRUE),
        sdSamplesBeforeLast = sd(samplesBeforeLast, na.rm = TRUE),
        n = n()
      ),
    
    negGr_beforeLastMean_year_target = cd1_target |> 
      filter(!is.na(negGrowth)) |> 
      group_by(speciesGG, riverGG, seasonGG, year, negGrowth) |> 
      summarize(
        meanGR = mean(grWeight, na.rm = TRUE),
        sdGR = sd(grWeight, na.rm = TRUE),
        meanSamplesBeforeLast = mean(samplesBeforeLast, na.rm = TRUE),
        sdSamplesBeforeLast = sd(samplesBeforeLast, na.rm = TRUE),
        n = n()
      ),
    
    dGAM_target = cd1_target |> 
      filter(!is.na(grWeight)) |> 
      mutate(age = year - cohort,
             ageF = factor(age),
             cohortF = factor(cohort)) |> 
      dplyr::select(riverGG, seasonGG, speciesGG, observedWeight, cohort, age, cohortF, ageF, 
                    #meanTemperatureScaledBySeason, meanFlowScaledBySeason,  meanFlowByRiverScaledBySeason,
                    #meanTemperatureScaledBySeasonRiver, meanFlowScaledBySeasonRiver, meanFlowByRiverScaledBySeasonRiver,
                    starts_with("mean"),
                    starts_with("propA"),
                    starts_with("propB"),
                    grWeight) |> 
      rename(tempS = meanTemperatureScaledBySeason, 
             flowS = meanFlowScaledBySeason, 
             flowByRiverS = meanFlowByRiverScaledBySeason,
             flowByAreaS = meanFlowByArea_flowExtScaledBySeason, #using just ..._flowExt for now (not adding _ByRiver)
             tempSR = meanTemperatureScaledBySeasonRiver, 
             flowSR = meanFlowScaledBySeasonRiver, 
             flowByRiverSR = meanFlowByRiverScaledBySeasonRiver,
             flowByAreaSR = meanFlowByArea_flowExtScaledBySeasonRiver
             ) |> 
      filter(tempSR > -4),
    
    ###############################################
    envIn_target = envDataWB_target |> 
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
        ),
    
    propNegSRS_target = cd1_target |> 
      group_by(speciesGG, riverGG, seasonGG) |> 
      summarize(numNeg = sum(negGrowth, na.rm = TRUE),
                n = n()
      ) |> 
      mutate(
        numPos = n - numNeg,
        propPos = numPos/n,
        propNeg = numNeg/n
      ),
    
    propNegSRsN_target = cd1_target |> 
      group_by(speciesGG, riverGG, sampleNumber, seasonGG, year) |> 
      summarize(numNeg = sum(negGrowth, na.rm = TRUE),
                n = n(),
                mT = mean(meanTemperature, na.rm = TRUE),
                mF = mean(meanFlowByRiver, na.rm = TRUE),
                mF_log10 = mean(log10(meanFlowByRiver), na.rm = TRUE)
      ) |> 
      mutate(numPos = n - numNeg,
             propPos = numPos/n,
             propNeg = numNeg/n) |> 
      ungroup() |> 
      left_join(indCountsBySpp_target) |> 
      left_join(indCounts_target),
    
   #### 
   meanNegSRsN_target = cd1_target |> 
      group_by(speciesGG, riverGG, seasonGG, sampleNumber, year, negGrowth) |> 
      summarize(meanNegPos = mean(grWeight, na.rm = TRUE),
                n = n()) |> 
      ungroup(),  
   
   meanNegSRsNWide_target = meanNegSRsN_target |> 
      filter(speciesGG == "Brook trout", !is.na(negGrowth)) |> 
      dplyr::select(-n) |> 
      pivot_wider(names_from = negGrowth, values_from = meanNegPos),
   ### 
 
    
    
    envIn_propNeg_target = envIn_target |> 
      left_join(
        propNegSRsN_target |> 
          dplyr::select(riverGG, year, seasonGG, speciesGG, propNeg)
      ) |> 
      mutate(
        yearSeason = paste(year, seasonGG, sep = "_")
      ),
    
    # Merge propNeg into daily environmental data   
    # So we can style the graph by propNeg growth for each year/season
    indCountsBySpp_target = cd1_target |> 
      group_by(riverGG, seasonGG, speciesGG, year) |> 
      summarize(nIndBySpp = n()),
    
    indCounts_target = cd1_target |> 
      group_by(riverGG, seasonGG, year) |> 
      summarize(nInd = n()),
    
    propNegLabels_target = envIn_propNeg_target |> 
      group_by(riverGG, seasonGG, speciesGG, year) |> 
      summarize(propNeg = mean(propNeg),
                minTemp = min(temperature),
                minFlow = min(flowByRiver),
                maxTemp = max(temperature),
                maxFlow = max(flowByRiver)) |> 
      left_join(indCountsBySpp_target) |> 
      left_join(indCounts_target)
    
    
    )
