### Main Drake overiview
library(drake)
library(here)

r_outdated(r_args = list(stderr = here("tmp", "drakemistake.log")))
r_make()
