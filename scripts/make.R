### Main Drake overiview
library(drake)
library(here)

# r_outdated(r_args = list(stderr = here("tmp", "drakemistake.log")))
r_make()

# # Plot the plan
# source(file = here("scripts", "Plan.R"))
# conf = drake_config(plan = plan)
# vis_drake_graph(conf)
