tar_option_set(packages = c("tidyverse", "lubridate", "getWBData", "daymetr"))

library(daymetr) # not sure why this is needed here, but without it get can't find download_daymet error

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
      values_to = "flowByRiver") %>%
    mutate(river = recode(river, "WL" = 'west brook' ,"JB" = 'wb jimmy',"MB" = 'wb mitchell',"OB" = "wb obear"),
           dateDate = as_date(date),
           flowByRiverm3s = flowByRiver * 0.028316847) %>%
    mutate(
      riverOrdered = factor(river, levels = c('west brook','wb jimmy','wb mitchell',"wb obear"),
                                   labels = c("West Brook","WB Jimmy","WB Mitchell","WB OBear"), ordered = T)
    ) %>%
    select(!c("flow", "cfs", "date")) 
  return(wbFlow2)
}

getEnvData_target <-
  tar_plan(
    WB_daymet_target = getDaymet(42.43896889699634, -72.67994313694251, 1997, 2021),
    flowByRiver_target = getFlowByRiver(), # data from Jenn's modeling
    
    envDataWB_target = 
      tbl(conDplyr, "data_daily_temperature") %>% 
      collect(n = Inf) %>% 
      full_join(tbl(conDplyr, "data_flow_extension") %>% 
                  collect(n = Inf), by = c("river", "date")) %>% 
      select(-source) %>% 
      rename(temperature = daily_mean_temp, flow = qPredicted) %>%
      mutate(dateDate = as_date(date),
             yday = yday(dateDate),
             year = year(dateDate)) %>%
      mutate(riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                                   labels = c("West Brook","WB Jimmy","WB Mitchell","WB OBear"), ordered = T)
      ) %>%
      left_join(WB_daymet_target, by = c('year', 'yday')) %>%
      left_join(flowByRiver_target)
  )  

######################################################
# tmp = envDataWB
# flowByRiver = getFlowByRiver()
# tmp2=tmp%>%left_join(flowByRiver)
# ggplot(tmp2%>%filter(year==2002), aes(dateDate,flowByRiverm3s)) +
#   geom_point() +
#   geom_point(aes(dateDate, flow, color = "blue")) +
#   facet_wrap(~river)
