library(targets)
library(tarchetypes)
library(quarto)

# run when there is a change in the data
#tar_invalidate(c("envDataWB_Target", "cdWB_electro0_target", ...)

tar_make()
quarto::quarto_render(output_format = "html")


quarto::quarto_render("dataAll.qmd", output_format = "html")
