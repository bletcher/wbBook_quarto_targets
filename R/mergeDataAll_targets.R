tar_option_set(packages = c("tidyverse", "lubridate", "getWBData"))

dataAll_target <-
  tar_plan(
      cdWB_all = cdWB_electro_target %>% 
      add_row(cdWB_wanding0_target %>% select(-riverOrdered)) %>% 
      add_row(cdWB_antenna0_target %>% select(-riverOrdered)) %>%
      mutate(
        riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                              labels = c("West Brook","WB Jimmy","WB Mitchell","WB OBear"), ordered = T
        )
      )
    
  )
