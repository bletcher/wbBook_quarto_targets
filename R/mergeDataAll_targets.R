tar_option_set(packages = c("tidyverse", "lubridate", "getWBData"))

electro <- tar_read(cdWB_electro_target)
wanding <- tar_read(cdWB_wanding0_target) %>% dplyr::select(-riverOrdered)
antenna <- tar_read(cdWB_antenna0_target) %>% dplyr::select(-riverOrdered)

w_e <- setdiff(names(wanding), names(electro))
a_e <- setdiff(names(antenna), names(electro))
a_w <- unique(c(w_e, a_e))

df <- data.frame(matrix(NA, dim(electro)[1],length(a_w)))
colnames(df) <- a_w

dataAll_target <-
  tar_plan(
      # cdWB_all = cdWB_electro_target %>% 
      # add_row(cdWB_wanding0_target %>% dplyr::select(-riverOrdered)) %>% 
      # add_row(cdWB_antenna0_target %>% dplyr::select(-riverOrdered)) %>%
      # mutate(
      #   riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell', "wb obear"),
      #                                labels = c("West Brook", "WB Jimmy", "WB Mitchell", "WB OBear"), ordered = T
      #   )
      # )
      
      
    cdWB_all_target = electro %>% 
      add_column(df) |> #add empty cols that are in wanding and antenna but not in electro
      add_row(wanding) |> 
      add_row(antenna) |> 
      mutate(
        riverOrdered = factor(river, levels = c('west brook', 'wb jimmy', 'wb mitchell', "wb obear"),
                              labels = c("West Brook", "WB Jimmy", "WB Mitchell", "WB OBear"), ordered = T
        )
      )
    
  )

