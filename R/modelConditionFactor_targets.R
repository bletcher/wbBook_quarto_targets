tar_option_set(packages = c("tidyverse", "lubridate"))
# model code in getDataElectro_targets.R in addCF()

modelConditionFactor_target <-
  tar_plan(
    relCF_byYear_target = cd1_target |> 
      group_by(speciesGG, riverGG, sampleNumber, seasonGG, year) |> 
      summarize(
        mRCF = mean(relCF, na.rm = TRUE),
        sdRCFF = sd(meanFlowByRiver, na.rm = TRUE),
        n = n()
      )
  )