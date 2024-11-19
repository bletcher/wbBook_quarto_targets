library(targets)
library(tarchetypes)
library(quarto)
library(knitr)

library(getPrepareWBData)
library(getWBData)
# run when there is a change in the data and to force target re-run
#tar_invalidate(c("envDataWB_target", "cdWB_electro0_target", ...)
# tar_invalidate(ends_with("_target"))

# to convert markdown headers to quarto in-body chunk options
#convert_chunk_header("dataAll.qmd", output = identity)

# may need to do this if packages don't load
# install.packages("xxxx", dependencies=TRUE, repos='http://cran.rstudio.com/')
tar_make()


# for parallel runs
# may need to change 'dummy  ~ dnorm(0,1)' in the model code to force targets to restart.
#or tar_invalidate(c(  'modelCMR_ttt_ft_cohort_WB_flow_target',
#                   'modelCMR_ttt_ft_cohort_WB_flowByRiver_target'))
tar_make_future(workers = 2)

quarto::quarto_render(output_format = "html")


quarto::quarto_render("modelsCMR_NN_SIM.qmd", 
                      #cache_refresh = TRUE, # default is FALSE
                      output_format = "html")


quarto::quarto_render("getDataElectroDataRelease.qmd", output_format = "html")#, cache_refresh = TRUE)

quarto::quarto_render("getDataElectro.qmd", 
                      cache_refresh = TRUE,
                      output_format = "html")

# working chapter
quarto::quarto_render("getDataAntenna.qmd", 
                      #cache_refresh = TRUE, # default is FALSE
                      output_format = "html")


quarto::quarto_render("modelsCMR_Growth_NN_OB.qmd", 
                      #cache_refresh = TRUE, # default is FALSE
                      output_format = "html")


