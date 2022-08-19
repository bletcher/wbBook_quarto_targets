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
      dplyr::filter(species %in% c( "bkt","bnt","ats" ))
  )
