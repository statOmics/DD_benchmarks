## Load top-level .Rprofile to activate renv
load_renv <- function() {
    renv_wd <- rprojroot::find_root(criterion = "renv.lock")
    old_wd <- setwd(renv_wd)
    on.exit(setwd(old_wd))
    source(".Rprofile")
}
load_renv()
rm(list = ls())
