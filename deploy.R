


rm(list=ls())
appDir <- here::here("app")
print(appDir)
setwd(appDir)


rsconnect::deployApp(
  forceUpdate = TRUE,
  appName = "AnaerobicDigestion",
  appTitle = "AnaerobicDigestion",
  appDir = ".",
  account = "kylebarrett", 
  server = "shinyapps.io",
  contentCategory = "site",
  appFiles = c(
    list.files(appDir, recursive = TRUE)
  ),
  lint = FALSE
)
