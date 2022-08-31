tar_option_set(packages = c("tidyverse", "lubridate", "getWBData"))

dataAntenna_target <-
  tar_plan(
    cdWB_antenna0_target = 
      cdWB_antenna0 <- createCoreData(
        sampleType=c("stationaryAntenna"), 
        whichDrainage = "west",
        columnsToAdd=c(
          "river",
          "riverMeter",
          "survey",
          #"section",
          "readerID",
          "comment"
        )
      ) %>%  
      filter(!is.na(tag)) %>% # for now
      addTagProperties(columnsToAdd = c(
        "cohort",
        "species",
        "dateEmigrated",
        "sex",
        "species")
      ) %>%
      mutate(riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                                   labels = c("West Brook","WB Jimmy","WB Mitchell","WB OBear"), ordered = T)
      )
  )


# May need to merge in sites data as below (from originalwestBrook-book project:

# Sites table

# # get sites table
# reconnect()
# sitesIn <- data.frame(tbl(conDplyr,"data_sites") )
# sites <- sitesIn %>% filter(is.na(quarter) & !is.na(quarter_length) & drainage == 'west') %>% select(-quarter)
# sites$section <- as.numeric(sites$section)

# 
# Merge antenna data into tagging data and wanding data

### Prepare data

# 
# # some formatting fixes
# #cdWB_antenna0$sectionOriginal <- cdWB_antenna0$section
# #cdWB_antenna0$section <- as.numeric( cdWB_antenna0$section )
# #cdWB_antenna0$inside <- ifelse( cdWB_antenna0$section %in% 1:47 | cdWB_antenna0$survey == "stationaryAntenna", T, F ) 
# 
# cdWB_antenna0$year <- year(cdWB_antenna0$detectionDate)
# cdWB_antenna0$yday <- yday(cdWB_antenna0$detectionDate)
# 
# cdWB_antenna0$riverOrdered <- factor(cdWB_antenna0$river,levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),labels=c("west brook","wb jimmy","wb mitchell","wb obear"), ordered=T)
# 
# cdWB_antenna0 <- cdWB_antenna0 %>%
#   group_by(tag) %>%
#   # arrange(tag,sampleNumber) %>%
#   mutate(lagSection = lead(section),
#          distMoved = section - lagSection,
#          minSample = min(sampleNumber),
#          maxSample = max(sampleNumber)) %>%
#   ungroup()
# 
# cdWB_antenna0$moveDir <- ifelse( cdWB_antenna0$section == cdWB_antenna0$lagSection, 0, ifelse( cdWB_antenna0$section > cdWB_antenna0$lagSection, 1,-1 ) )
# 
# cdWB_antenna0$drainage <- "west"
# 
# cdWB_antenna0$sizeForGraph <- ifelse( is.na(cdWB_antenna0$observedLength), 60, cdWB_antenna0$observedLength )




