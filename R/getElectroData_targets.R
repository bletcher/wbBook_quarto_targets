tar_option_set(packages = c("tidyverse", "lubridate", "getWBData", "getPrepareWBData"))

# default values for createCoreData()
# function (sampleType = "electrofishing", baseColumns = T, 
#    columnsToAdd = NULL, includeUntagged = F, whichDrainage = "west") 

drainage <- 'west'

getElectroData_target <-
  tar_plan(
    cdWB_electro0_target = createCoreData(
      sampleType = "electrofishing",  #"stationaryAntenna","portableAntenna"
      columnsToAdd = c("sampleNumber",
                       "river",
                       "survey",
                       "pass",
                       "observedLength",
                       "observedWeight",
                       "comments"),
      includeUntagged = TRUE,
      whichDrainage = "west"
    ) %>%
      addTagProperties(
        columnsToAdd = c("cohort",
                         "species",
                         "dateEmigrated",
                         "sex",
                         "species"
        )
      ) %>%
      dplyr::filter(species %in% c( "bkt","bnt","ats"),
                    area %in% c("trib","inside","below","above"),
                    !is.na(sampleNumber)) %>%
      addSampleProperties() %>%
      addEnvironmental(),
    
    # functions in getPrepareWBData library
    cdWB_electro_target = cdWB_electro0_target %>%
      cleanData(drainage) %>%
      mergeSites(drainage) %>%
      addNPasses(drainage) %>%
      mutate(drainage = drainage)
    
  )  


#####################################
## getData functions 
#####################################