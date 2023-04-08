options(
  repos = c(
    MPN = "https://mpn.metworx.com/snapshots/stable/2022-10-20",
    CRAN = "https://cran.rstudio.com"
  )
)

# Run renv/activate.R after setting repos to bootstrap from MPN (if needed)
source("renv/activate.R")


if(interactive()){
  message("repos set to: \n\t", paste0(unique(getOption('repos')), collapse = "\n\t"))
  message("library paths set to: \n\t", paste0(.libPaths(), collapse = "\n\t"))
  
}
