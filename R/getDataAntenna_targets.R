tar_option_set(packages = c("tidyverse", "lubridate", "getWBData"))


tmp=data.frame(tbl(conDplyr,"raw_tags_antenna"))

dataAntenna_target <-
  tar_plan(
    cdWB_antenna0_target = createCoreData(
        sampleType=c("stationaryAntenna"), 
        whichDrainage = "west",
        columnsToAdd=c(
          "river",
          "riverMeter",
          "survey",
      #"section",
          "readerID",
         # "antenna_name",
          "comments"
        )
      ) %>%  
      filter(!is.na(tag)) %>% # for now
      addTagProperties(
        columnsToAdd = c(
          "cohort",
          "species",
          "dateEmigrated",
          "sex",
          "species"
        )
      ) %>%
      mutate(riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                                         labels = c("West Brook", "WB Jimmy", "WB Mitchell", "WB OBear"), ordered = T)
      ) %>%
      updateAntennaData(sites_target), # not updating with sites_target for now
    
    sites_target = data.frame(tbl(conDplyr,"data_sites")) %>%
      filter(is.na(quarter) & !is.na(quarter_length) & drainage == 'west') %>% 
      select(-quarter) %>%
      #mutate(section = as.numeric(section)) %>%
      rename(riverMeter = river_meter),
    
    antennaDeployment_target = data.frame(tbl(conDplyr,"antenna_deployment")) %>%
      rename(riverMeter = river_meter) %>%
      mutate(
        riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                               labels = c("West Brook", "WB Jimmy", "WB Mitchell", "WB OBear"), ordered = T
        )),
    
    antennaLog_target = data.frame(tbl(conDplyr,"data_antenna_log")) %>%
      rename(riverMeter = river_meter) %>%
      mutate(
        riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                              labels = c("West Brook", "WB Jimmy", "WB Mitchell", "WB OBear"), ordered = T
        )),
    
    antennaStatus_target = data.frame(tbl(conDplyr,"data_antenna_status")) %>%
      rename(riverMeter = river_meter) %>%
      mutate(
        riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                              labels = c("West Brook", "WB Jimmy", "WB Mitchell", "WB OBear"), ordered = T
        ))
  )

#################################
### Functions
#################################

updateAntennaData <- function(d0, s) {

  d <- d0 #left_join(d0, s) # need to make this work with section...
  
  # some formatting fixes
  #d$sectionOriginal <- as.character(d$section)
  #d$section <- as.numeric(d$section)
  #d$inside <- ifelse(d$section %in% 1:47 | d$survey == "stationaryAntenna", T, F)

  d$year <- year(d$detectionDate)
  d$yday <- yday(d$detectionDate)
  d$date <- date(d$detectionDate)

  # d <- d %>%
  #   group_by(tag) %>%
  #   # arrange(tag,sampleNumber) %>%
  #   mutate(lagSection = lead(section),
  #          distMoved = section - lagSection
  #          #minSample = min(sampleNumber),
  #          #maxSample = max(sampleNumber)
  #          ) %>%
  #   ungroup()
  # 
  # d$moveDir <- ifelse(d$section == d$lagSection, 0, ifelse(d$section > d$lagSection, 1, -1))

  d$drainage <- "west"

  #d$sizeForGraph <- ifelse( is.na(d$observedLength), 60, d$observedLength )
  
  return(d)
}

