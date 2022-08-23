library(targets)
library(quarto)

# run when there is a change in the data
#tar_invalidate(c("envDataWB_Target", "cdWB_electro0_target", ...)

tar_make()
quarto::quarto_render(output_format = "html")
