#This target is for reading the data from the original Access database for the non-WB (Catamaran CATAMARAN   Sawmill   SAWMILL    Shorey    SHOREY ) rivers. 
#Stanley is in its own Access DB.

tar_option_set(packages = c("tidyverse", "lubridate", "getWBData", "getPrepareWBData", "stringr"))

#saved `lengths and weights by fish_noWB` as xlsx file from access database `C:\Users\bletcher\OneDrive - DOI\PITTAGMAIN\LWbF table linked`
#Then saved lengths and weights by fish_noWB.csv from the xlsx file to read in below

find_negative_9999 <- function(df) {
  df %>%
    # Map across all columns
    purrr::map_lgl(~ any(. == -9999, na.rm = TRUE)) %>%
    # Get names of columns that returned TRUE
    which() %>%
    names()
}

replace_negative_9999 <- function(df) {
  df %>%
    mutate(across(everything(), ~ifelse(. == -9999, NA, .)))
}

getAllDAta_noWB_target <-
  tar_plan(
    cdNoWB_all_target = read.csv("C:/Users/bletcher/OneDrive - DOI/PITTAGMAIN/lengths and weights by fish_noWB.csv") |>
      filter(River != "west brook") |>
      dplyr::select(
        id = ID,
        source = Source,
        sampleType = Sample.type,
        drainage = Drainage,
        species = Species,
        river = River,
        survey = Survey,
        tag = Tag.number,
        #cohort = cohort,
        firstSeen = FirstSeen,
        lastSeen = LastSeen,
        sampleNumber = Sample.Num,
        section = Section,
        area = Area,
        sampleName = SampName,
        season = Season,
        everEmMain = EverEmMain,
        riverTagged = RiverTagged,
        areaTagged = AreaTagged,
        riverKM = River.KM,
        date = Date,
        medianSampleDate = MedianDateSamp,
        fishNumber = Fishnumber,
        recap = Recap,
        yearOfStocking = MinOfYear.of.Stocking,
        estimatedYearOfStocking = Estimated.MinOfYear.of.Stocking,
        age = Age,
        observedLength = Length,
        observedLength2 = Length2,
        grLen,
        grOstrovskyLength,
        observedMass = Weight,
        predMass = PredMass,
        observedMass2 = Mass2,
        grMass,
        grOstrovskyMass = grOstrovskyMass,
        grInterval,
        cf = CF,
        mature = Mature,
        everMat = EverMat,
        #retro2Smolt = RETRO_.2..smolt,
        #retro3Smolt = RETRO_.3..smolt,
        #retroRes = RETRO_res,
        geneticSample = Genetic.Sample,
        readerID = ReaderID,
        aliveDead = Alive.Dead,
        antInstance = AntInstance,
        pass = Pass,
        effort = Effort,
        comment = Comment,
        quarterWidth = Quart.Width,
        quarterLength = Quart.Length
      ) |>
      mutate(
        drainage = str_to_lower(drainage),
        river = str_to_lower(river)
        #dateDate = (mdy(date)),
       # medianSampleDateDate = mdy(medianSampleDate)
      ) |>
      replace_negative_9999()
  )



#tmp=replace_negative_9999(cdNoWB_all_target)
