tar_option_set(packages = c("tidyverse", "lubridate", "getWBData"))


getEnvData_target <-
  tar_plan(
    envDataWB_target = tbl(conDplyr, "data_daily_temperature") %>% 
      collect(n = Inf) %>% 
      full_join(tbl(conDplyr, "data_flow_extension") %>% 
                  collect(n = Inf), by = c("river", "date")) %>% 
      select(-source) %>% 
      rename(temperature = daily_mean_temp, flow = qPredicted) %>%
      mutate(dateDate = as_date(date),
             yday = yday(dateDate))%>%
      mutate(riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                                   labels = c("West Brook","WB Jimmy","WB Mitchell","WB OBear"), ordered = T)
      )
  )  


#####################################
## getData functions 
#####################################