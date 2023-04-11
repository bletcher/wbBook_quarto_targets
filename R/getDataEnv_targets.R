tar_option_set(packages = c("tidyverse", "lubridate", "getWBData", "daymetr", "fuzzyjoin"))

library(daymetr) # not sure why this is needed here, but without it get can't find download_daymet error



getEnvData_target <-
  tar_plan(
    WB_daymet_target = getDaymet(42.43896889699634, -72.67994313694251, 1997, 2021),
    flowByRiver_target = getFlowByRiver(), # data from Jenn's modeling
    #flowByArea_target = getFlowByArea()
    
    envDataWB_target = 
      tbl(conDplyr, "data_daily_temperature") %>% 
      collect(n = Inf) %>% 
      full_join(tbl(conDplyr, "data_flow_extension") %>% 
                  collect(n = Inf), by = c("river", "date")) %>% 
      dplyr::select(-source) %>% 
      rename(temperature = daily_mean_temp, flow = qPredicted) %>%
      mutate(dateDate = as_date(date),
             yday = yday(dateDate),
             year = year(as_date(dateDate))) %>%
      mutate(riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                                   labels = c("West Brook","WB Jimmy","WB Mitchell","WB OBear"), ordered = T)
      ) %>%
      left_join(WB_daymet_target, by = c('year', 'yday')) %>%
      left_join(flowByRiver_target) |> 
      addSeason(tar_read(medianDates_target)) |> 
      addFlowToTribs() |> 
      addFlowByRiverToTribs() |> 
      getFlowByArea()
  ) 


#### Functions

getDaymet <- function(lat, lon, start_year, end_year) {
  WB_daymet_list <- download_daymet(site = "WestBrook", lat = lat, lon = lon, 
                                    start = start_year, end = end_year)
  out <- WB_daymet_list$data %>%
    rename(dayLength = dayl..s., precip = prcp..mm.day., solarRadiation = srad..W.m.2.,
           swe = swe..kg.m.2., airTempMax = tmax..deg.c., airTempMin = tmin..deg.c.,
           vaporPressure = vp..Pa.) %>%
    mutate(airTempMedian = airTempMin + (airTempMax - airTempMin)/2)
  
  return(out)
}

# loads 'WBFlow'
# for now, fixed dataframe uploaded to ../dataIn. Will want to automate uploading...
getFlowByRiver <- function(){
  load('./dataIn/wbFlow/WBFlow_fromJenn.RData')
  wbFlow2 <- WBFlow %>%
    select(date, ETmm, Pmm, OB_Flow_cfs, WL_Flow_cfs, MB_Flow_cfs, JB_Flow_cfs) %>%
    pivot_longer(
      cols = ends_with("_Flow_cfs"),
      names_to = c("river", "flow", "cfs"),
      names_sep = "_",
      values_to = "flowByRiver_cfs") %>%
    mutate(river = recode(river, "WL" = 'west brook' ,"JB" = 'wb jimmy',"MB" = 'wb mitchell',"OB" = "wb obear"),
           dateDate = as_date(date),
           flowByRiver = flowByRiver_cfs * 0.028316847) %>%
    mutate(
      riverOrdered = factor(river, levels = c('west brook','wb jimmy','wb mitchell',"wb obear"),
                                   labels = c("West Brook","WB Jimmy","WB Mitchell","WB OBear"), ordered = T)
    ) %>%
    select(!c("flow", "cfs", "date")) 
  return(wbFlow2)
}

addFlowToTribs <- function(dIn) {
  flowByDate <- dIn |> 
    filter(river == "west brook") |> 
    rename(flowWithTribs = flow) |> 
    dplyr::select(date, flowWithTribs)
  
  dOut <- left_join(dIn, flowByDate)
  return(dOut)
}

# adds WB flow estimates to tribs for scaling in getFlowByARea()
addFlowByRiverToTribs <- function(dIn) {
  flowByDate <- dIn |>
    filter(river == "west brook") |>
    rename(flowByRiverWB_WithTribs = flowByRiver) |>
    dplyr::select(date, flowByRiverWB_WithTribs)

  dOut <- left_join(dIn, flowByDate)
  return(dOut)
}

getFlowByArea <- function(dIn) {
  dOut <- dIn |> 
    mutate(
      propRiverArea = 
        case_when(
          river == "wb jimmy" ~ 2.460 / 21.756,
          river == "wb mitchell" ~ 1.036 / 21.756,
          river == "wb obear" ~ 1.295 / 21.756,
          river == "west brook" ~ 21.756 / 21.756
        ),
      flowByArea_flowExt = flowWithTribs * propRiverArea,
      flowByArea_ByRiver = flowByRiverWB_WithTribs * propRiverArea,
      dummy = 1
    )
  return(dOut)
}

getFlowByArea_hold <- function(dIn) {
  dOut <- dIn |> 
    mutate(
      propRiverArea = 1,
      flowByArea = flowWithTribs * propRiverArea
    )
}

addSeason <- function(d = d, medDate = tar_read(medianDates_target)){

  medDateSeason <- fuzzy_left_join(
    d |> dplyr::select(river, date), medDate,
    by = c(
      "river" = "river",
      #"year" = "year",
      "date" = "start",
      "date" = "end"
    ),
    match_fun = list(`==`, `>=`, `<`)
  ) %>%
    dplyr::select(river=river.x, date, season, start, end)
  
  return(d |> left_join(medDateSeason))
}
######################################################
# tmp = envDataWB
# flowByRiver = getFlowByRiver()
# tmp2=tmp%>%left_join(flowByRiver)
# ggplot(tmp2%>%filter(year==2002), aes(dateDate,flowByRiverm3s)) +
#   geom_point() +
#   geom_point(aes(dateDate, flow, color = "blue")) +
#   facet_wrap(~river)
