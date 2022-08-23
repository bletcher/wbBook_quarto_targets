tar_option_set(packages = c("tidyverse", "tidyr", "lubridate", "getWBData", "targets", "tarchetypes"))

modelFlow_target <-
  tar_plan(
    dataFlow_target = read.csv("./dataIn/wbFlow/EcoDrought_Continuous_MA.csv"),
    
    dFlow_target = dataFlow_target %>%
      filter(Site_Name %in% c("Jimmy Brook", "Mitchell Brook", "Obear Brook Lower", "West Brook 0")) %>%
      mutate(date = mdy_hm(DateTime_EST),
             site = recode(Site_Name, "Jimmy Brook" = "OL", "Mitchell Brook" = "OS", "Obear Brook Lower" = "IS", "West Brook 0" = "WB"),
             dischargeLog = log(Discharge_Hobo_cfs + 0.01)),
    
   
    # hard-coded for now
    dFlowWide_target = dFlow_target %>% 
      pivot_wider(id_cols = date, 
                  values_from = dischargeLog, 
                  names_from = site
      ) %>%
      mutate(
        OLScaled = scaleCol(OL),
        ISScaled = scaleCol(IS),
        OSScaled = scaleCol(OS),
        WBScaled = scaleCol(WB),
        yday = yday(date),
        year = year(date)
      )
  )

scaleCol <- function(d){
  return (d - mean(d, na.rm = TRUE)) / sd(d, na.rm = TRUE)
}