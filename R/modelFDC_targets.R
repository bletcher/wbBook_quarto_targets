tar_option_set(packages = c("tidyverse", "lubridate", "hydroTSM"))

modelFDC_target <-
  tar_plan(
    flowVars_target = c("flowByArea_ByRiver", "flowByArea_flowExt", "flowByRiver"),
    fdcStats_target =
      envDataWB_target |> 
        filter(!is.na(flowByRiverWB_WithTribs)) |> 
        dplyr::select(river, date, all_of(flowVars_target)) |> 
        pivot_longer(cols = flowVars_target) |> 
        rename(flow = value, variable = name) |> 
        group_by(river, variable) |>
        reframe(
          getFDC(flow)
        ),
    fdcStatsY_target =
      envDataWB_target |> 
        filter(!is.na(flowByRiverWB_WithTribs)) |> 
        dplyr::select(river, date, year, all_of(flowVars_target)) |> 
        pivot_longer(cols = flowVars_target) |> 
        rename(flow = value, variable = name) |> 
        group_by(river, year, variable) |>
        reframe(
          getFDC(flow)
        ),
    fdcStatsS_target =
      envDataWB_target |>
        filter(!is.na(flowByRiverWB_WithTribs)) |> 
        dplyr::select(river, date, season, all_of(flowVars_target)) |> 
        pivot_longer(cols = flowVars_target) |> 
        rename(flow = value, variable = name) |> 
        group_by(river, season, variable) |>
        reframe(
          getFDC(flow)
        ),
    fdcStatsYS_target =
      envDataWB_target |>
        filter(!is.na(flowByRiverWB_WithTribs), !is.na(season)) |> 
        dplyr::select(river, date, year, season, all_of(flowVars_target)) |> 
        pivot_longer(cols = flowVars_target) |> 
        rename(flow = value, variable = name) |> 
        group_by(river, year, season, variable) |>
        reframe(
          getFDC_YS(flow)
        ),
    
    fdcThresholds_target = c(0.1, 0.9),
    
    # do this in getEnvData so we don't need to redo in getDataElectro
    # no, use envDataWB_fdcThresh_target  in getDataElectro
    # envDataWB_fdcThresh_target =
    #   envDataWB_target |> 
    #     addFDC_stats(fdcStatsS_target, fdcThresholds_target),
    
    envDataWB_fdcThresh_target =
      envDataWB_target |> 
      left_join(fdcThreshWide_target) |> 
      mutate(
        belowLoFlowThreshByRiver = flowByRiver < flowByRiver_0.9, #will need to make these names not hardcoded
        aboveHiFlowThreshByRiver = flowByRiver > flowByRiver_0.1,
        belowLoFlowThreshByArea_flowExt = flowByArea_flowExt < flowByArea_flowExt_0.9,
        aboveHiFlowThreshByArea_flowExt = flowByArea_flowExt > flowByArea_flowExt_0.1
#addBack        belowLoFlowThreshByArea_ByArea = flowByArea_ByArea < flowByArea_ByArea_0.9,
#        aboveHiFlowThreshByArea_ByArea = flowByArea_ByArea > flowByArea_ByArea_0.1
      ),
    
    fdcThreshWide_target =   
      fdcStatsS_target |>
      filter(stat %in% c(fdcThresholds_target[1], fdcThresholds_target[2]), !is.na(season)) |> 
      pivot_wider(
        names_from = c(variable, stat),
        values_from = perc
      ), 
    
    # use this only for seasonal comparisons. Is not fish specific.
    propFDC_aboveBelow_target = 
      envDataWB_fdcThresh_target |> 
        filter(!is.na(season)) |> 
        group_by(river, year, season) |> 
        summarize(
          propBelowLoFlowThreshByRiver = sum(belowLoFlowThreshByRiver, na.rm = TRUE) / sum(!is.na(belowLoFlowThreshByRiver)),
          propAboveHiFlowThreshRiver = sum(aboveHiFlowThreshByRiver, na.rm = TRUE) / sum(!is.na(aboveHiFlowThreshByRiver)),
          propBelowLoFlowThreshByArea_flowExt = sum(belowLoFlowThreshByArea_flowExt, na.rm = TRUE) / sum(!is.na(belowLoFlowThreshByArea_flowExt)),
          propAboveHiFlowThreshByArea_flowExt = sum(aboveHiFlowThreshByArea_flowExt, na.rm = TRUE) / sum(!is.na(aboveHiFlowThreshByArea_flowExt))
          
# addBack         propBelowLoFlowThreshByArea_ByArea = sum(belowLoFlowThreshByArea_ByArea, na.rm = TRUE) / sum(!is.na(belowLoFlowThreshByArea_ByArea)),
#          propAboveHiFlowThreshByArea_ByArea = sum(aboveHiFlowThreshByArea_ByArea, na.rm = TRUE) / sum(!is.na(aboveHiFlowThreshByArea_ByArea))
        ) |> 
      ungroup() |> 
      addGG_noSpecies()  # in generalFunctions.R
  )

#### Functions ###########
getFDC <- function(x, stat = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.1, 0.2, 0.3, 0.5, 
                                0.7, 0.8, 0.9, 0.975, 0.99, 0.995, 0.999, 0.9995, 0.9999)
                    ) {
  y <- fdc(x, plot=FALSE)
  f <- approxfun(y, x, method="linear", ties=mean)
  #perc <- f(stat)
  return(
    tibble(
      perc = f(stat),
      stat = stat
    )
  )
}

getFDC_YS <- function(x, stat = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.1, 0.2, 0.3, 0.5, 
                                  0.7, 0.8, 0.9, 0.975, 0.99, 0.995, 0.999, 0.9995, 0.9999)){
  
  if(length(x) > 3) {
    y <- fdc(x, plot=FALSE)
    f <- approxfun(y, x, method="linear", ties=mean)
    #perc <- f(stat)
    
    out <- 
      tibble(
        perc = f(stat),
        stat = stat
      )
  } else {
    out <- 
      tibble(
        perc = NA,
        stat = NA
      )
  }
  return(out)
}

addFDC_stats <- function(d, fdcIn, fdcThresh) {
  
  fdcThrshWide <- 
    fdcIn |>
    filter(stat %in% c(fdcThresh[1], fdcThresh[2]), !is.na(season)) |> 
    pivot_wider(
      names_from = c(variable, stat),
      values_from = perc
    ) 
  
  dOut <- d |> 
    left_join(fdcThrshWide) |> 
    mutate(
      belowLoFlowThreshByRiver = flowByRiver < flowByRiver_0.9, #will need to make these names not hardcoded
      aboveHiFlowThreshByRiver = flowByRiver > flowByRiver_0.1,
      belowLoFlowThreshByArea_flowExt = flowByArea_flowExt < flowByArea_flowExt_0.9,
      aboveHiFlowThreshByArea_flowExt = flowByArea_flowExt > flowByArea_flowExt_0.1
# addBack     belowLoFlowThreshByArea_ByArea = flowByArea_ByArea < flowByArea_ByArea_0.9,
#      aboveHiFlowThreshByArea_ByArea = flowByArea_ByArea > flowByArea_ByArea_0.1
    )
  
  return(dOut)
}


