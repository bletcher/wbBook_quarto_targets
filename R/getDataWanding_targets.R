tar_option_set(packages = c("tidyverse", "lubridate", "getWBData"))

dataWanding_target <-
  tar_plan(
    cdWB_wanding0_target = 
      createCoreData(
        sampleType = "portableAntenna",
        columnsToAdd = c("tag", 
                         "detectionDate", 
                         "river", 
                         "area", 
                         "section", 
                         "survey", 
                         "cohort",
                         "sampleName", 
                         "readerId", 
                         "aliveOrDead", 
                         "instance", 
                         "pass", 
                         "quarter", 
                         "leftOrRight", 
                         "habitat", 
                         "cover", 
                         "justification", 
                         "comment")
      ) %>% 
      addTagProperties() %>%
      dplyr::filter(species %in% c( "bkt","bnt","ats" )) %>%
      mutate(riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                            labels = c("West Brook","WB Jimmy","WB Mitchell","WB OBear"), ordered = T)
      )
  )
