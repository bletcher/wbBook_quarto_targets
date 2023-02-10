tar_option_set(packages = c("tidyverse", "lubridate"))

modelGrowthInMass_target <-
  tar_plan(
    cd1_target = cdWB_electro_target |> 
      filter(sampleNumberDiff == 1,
             tag %notin% c('1bf20ff490', '1bf20ebe4e', '1c2c582218')) |> 
      mutate(
        negGrowth = grWeight < 0,
        month = month(date),
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
      ),
    
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
                minFlow = min(flowByRiverm3s),
                maxTemp = max(temperature),
                maxFlow = max(flowByRiverm3s)) |> 
      left_join(indCountsBySpp_target) |> 
      left_join(indCounts_target)
    
    
    )