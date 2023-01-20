library(targets)
library(tarchetypes)
library(quarto)
library(knitr)

# run when there is a change in the data
#tar_invalidate(c("envDataWB_Target", "cdWB_electro0_target", ...)

# to convert markdown headers to quarto in-body chunk options
#convert_chunk_header("dataAll.qmd", output = identity)

tar_make()
quarto::quarto_render(output_format = "html")


quarto::quarto_render("dataAll.qmd", output_format = "html")
quarto::quarto_render("getDataEnv.qmd", output_format = "html")

quarto::quarto_render("modelGrowthInMass.qmd", output_format = "html")
