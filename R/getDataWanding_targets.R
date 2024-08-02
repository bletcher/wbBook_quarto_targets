tar_option_set(packages = c("tidyverse", "lubridate", "getWBData"))

dataWanding_target <-
  tar_plan(
    cdWB_wanding0_target = 
      createCoreData(
        sampleType = "portableAntenna",
        columnsToAdd = c(
          "tag", 
          "detectionDate", 
          "river", 
          "area", 
          "section", 
          "riverMeter", #added 5/2/2024
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
          "comments"
        )
      ) %>% 
      addTagProperties(
        columnsToAdd = c(
          "cohort",
          "species",
          "dateEmigrated",
          "sex",
          "species"
        )
      ) %>%
      dplyr::filter(species %in% c( "bkt","bnt","ats" )) %>%
      mutate(
        riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),
                            labels = c("West Brook","WB Jimmy","WB Mitchell","WB OBear"), ordered = T
                            ),
        section = as.numeric(section)
      ) %>%
      updateWandingData(),
    
    cdWB_wanding_target = cdWB_wanding0_target %>%
      group_by(tag) %>%
      mutate(sectionWQuarterLagged = lead(sectionWQuarter),
             detectionDateLagged = lead(detectionDate),
             moveDist = (sectionWQuarter - sectionWQuarterLagged) * 20,
             moveTime = as.numeric((difftime(detectionDateLagged, detectionDate, units="days"))) / 7,
             moveRate = moveDist/moveTime) %>%
      #filter(moveTime >= 1 & moveTime < 365) %>%
      mutate(month = as.numeric(month( detectionDate))) |> 
      ungroup(),
    
    
    #cd2
    cdWB_wandingTribs_target = cdWB_wanding_target %>%
      # group_by(tag) %>%
      # mutate(sectionWQuarterLagged = lead(sectionWQuarter),
      #         detectionDateLagged = lead(detectionDate),
      #         moveDist = (sectionWQuarter - sectionWQuarterLagged) * 20,
      #         moveTime = as.numeric((difftime(detectionDateLagged, detectionDate, units="days"))) / 7,
      #         moveRate = moveDist/moveTime) %>%
      filter(moveTime >= 1 & moveTime < 365),
      # mutate(month = as.numeric(month( detectionDate))),
    
    #cd3
    cdWB_wandingTribs3_target = cdWB_wandingTribs_target %>%
      filter(aliveOrDead != 'dead' & 
               year %in% c(2009,2010) & 
               river %in% c('wb mitchell','wb jimmy') & 
               month %in% 9:11) %>%
      mutate(yoy = ifelse(year == cohort, 1, 0)),
    
    #cd4
    #cdWB_wandingTribs4_target = cdWB_wandingTribs3_target %>%
    #  filter(moveTime < 50 * 7) #%>%
      #mutate(interval = cut_interval(moveTime, 5))

    # cd2wb 
    cdWB_wandingWB2_target = cdWB_wanding_target, #%>%
      # group_by(tag) %>%
      # mutate(sectionWQuarterLagged = lead(sectionWQuarter),
      #         detectionDateLagged = lead(detectionDate),
      #         moveDist = (sectionWQuarter-sectionWQuarterLagged) * 20,
      #         moveTime = as.numeric((difftime(detectionDateLagged, detectionDate, units="days"))) / 7,
      #         moveRate = moveDist/moveTime) %>%
      # #      filter( moveTime >= 1 & moveTime < 365 ) %>%
      # mutate(month = as.numeric(month(detectionDate))),
    
    cdWB_wandingWB3_target = cdWB_wandingWB2_target %>%
      filter(aliveOrDead != 'dead' & 
               year %in% c(2002,2003) & 
               river == 'west brook' & 
               species != 'bnt' & sectionWQuarter > 30 & 
               sectionWQuarter < 35 & 
               j > 50)
  )

#########################
### Functions
#########################

updateWandingData <- function(cd) {
  
  cd %>%
    mutate(
      sectionN = as.numeric(section),
      year = year(detectionDate),
      j = yday(detectionDate), #day of year

      quarter = ifelse(is.na(quarter), 1, quarter), # so we don't lose the data for sectionWQuarter when quarter is NA
      sectionWQuarter = sectionN + as.numeric(quarter) / 4 - 0.25
    )
}

#Code below from GitHub/wandingDataWB/wandingDataWBAnalysis.R

# 
# Thinking about using the wanding data to evaluate whether we can use the seasonal shocking data to talk about patterns of movement. If there is a lot of movement in the wanding data, may be difficult to justify using the 4 seasonal samples to describe movemenmts. If not much movement, may be ok.
# 

# ```
# 
