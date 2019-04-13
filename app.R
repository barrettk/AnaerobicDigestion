
pkg <- c("shiny", "mrgsolve", "shinyAce", "shinydashboard", "dplyr", "knitr", "plyr", "tidyverse", 
         "wrapr", "extrafont", "polynom", "ggplot2", "shinyWidgets", "gridExtra", "rmarkdown", 
         "markdown", "sn", "rlang", "lattice", "reshape", "reshape2", "magrittr", "stats")
new.pkg <- pkg[!(pkg %in% installed.packages())]
if (length(new.pkg)) {
  install.packages(new.pkg)
}
#Mrgsolve
library(mrgsolve,quietly = T)
options(mrgsolve_mread_quiet=TRUE)
#UI Packages
suppressMessages(library(shiny,quietly = T))
suppressMessages(library(shinyAce,quietly = T))
suppressMessages(library(shinydashboard,quietly = T))
#suppressMessages(library(shinyjs,quietly=T))
#library(shiny.semantic)
library(shinyWidgets,quietly = T)
suppressMessages(library(rsconnect,quietly = T))
library(rmarkdown); library(markdown)
library(knitr)

#Server Packages
suppressMessages(library(plyr,quietly = T))
suppressMessages(library(dplyr,quietly = T))
suppressMessages(library(ggplot2,quietly = T))
suppressMessages(library(tidyverse,quietly = T))
suppressMessages(library(sn,quietly = T))
library(stats)
library(rlang)
library(wrapr)
library(lattice,quietly = T)
library(gridExtra,quietly = T)
library(extrafont,quietly = T)
library(reshape,quietly = T)
library(reshape2,quietly = T)
library(polynom,quietly = T)
library(magrittr,quietly = T)


# Compile the model
mod <- mrgsolve::mread("AnaerobicDigestionShinyV3.cpp")
#assign("mod",mod,envir = globalenv())
#source("Bacteria_kinetics.R",local = T)
source("Bacteria_kinetics.R",local=globalenv())
.model <- paste(mrgsolve:::code(mod), collapse='\n')
assign(".model",.model,envir = globalenv())
SensChoices <- c("Acidogen Conc.","Acetogen Conc.","Methanogen Conc.","Bacteroid Conc.","Head Space Ratio","Temperature","Number of Wells","WT % of Guar Gum","WT % of PEG-400")


# Server Code -------------------------------------------------------------

server <- function(input, output,session) {

  mod <- mrgsolve::mread_cache("AnaerobicDigestionShinyV3.cpp")
  
  # Kinetic Data Import From Bacteria_kinetics.R

  # Acetogens
  # Vmax
  coef_Aceto_vmax <- coef(fit_Aceto_vmax)
  Intercept_Aceto_vmax <- as.numeric(coef_Aceto_vmax[1])
  a_Aceto_vmax <- as.numeric(coef_Aceto_vmax[2])
  b_Aceto_vmax <- as.numeric(coef_Aceto_vmax[3])
  c_Aceto_vmax <- as.numeric(coef_Aceto_vmax[4])
  # Km
  coef_Aceto_km <- coef(fit_Aceto_km)
  Intercept_Aceto_km <- as.numeric(coef_Aceto_km[1])
  a_Aceto_km <- as.numeric(coef_Aceto_km[2])
  b_Aceto_km <- as.numeric(coef_Aceto_km[3])
  c_Aceto_km <- as.numeric(coef_Aceto_km[4])
  # Methanogens
  # Vmax
  coef_Meth_vmax <- coef(fit_Meth_vmax)
  Intercept_Meth_vmax <- as.numeric(coef_Meth_vmax[1])
  a_Meth_vmax <- as.numeric(coef_Meth_vmax[2])
  b_Meth_vmax <- as.numeric(coef_Meth_vmax[3])
  c_Meth_vmax <- as.numeric(coef_Meth_vmax[4])
  # Km
  coef_Meth_km <- coef(fit_Meth_km)
  Intercept_Meth_km <- as.numeric(coef_Meth_km[1])
  a_Meth_km <- as.numeric(coef_Meth_km[2])
  b_Meth_km <- as.numeric(coef_Meth_km[3])
  c_Meth_km <- as.numeric(coef_Meth_km[4])
  
  PegDist <- reactive({
    shiny::req(input[["PegMeanMW"]])
    meanVal <- as.numeric(input[["PegMeanMW"]])
    meanVal2 <- (meanVal-18)/44
    test <- sn::rsn(5000, (meanVal2*1.267), -1.3, alpha=10) %>% as.vector() %>% as.integer()
    test <- test[test>4 & test<10]
    test
  })
  
  output$PegPlot <- renderPlot({
    test <- PegDist() %>% as.vector() %>% as.integer()
    test2 <- test*44+18 %>% as.vector() %>% as.integer()
    df <- data.frame(test2,test); names(df) <- c("MW","N")
    df$N <- df$N %>% as.character() %>% as.factor()
    breaks <- c(238, 282, 326, 370, 414)
    #Plot
    hist2 <- ggplot(df,aes(x=MW)) + labs(x="Molecular Weight",y="Density",caption = paste0("Mean MW: ",round(mean(test2),2))) +  
      theme(text=element_text(family="Times New Roman", face="bold", size=14),
            plot.caption = element_text(color = "red", face = "italic",size=14)) + expand_limits(x = c(min(test2), max(test2))) +
      geom_vline(aes(xintercept=mean(test2)),color="blue", linetype="dashed", size=1) + geom_density() +
      scale_x_continuous(breaks=breaks) 
    hist2
  })
  
  pars <- reactive({
    shiny::req(PegDist())
    shiny::req(input[["Bact_ScaleFact_Aceto"]])
    shiny::req(input[["Bact_ScaleFact_Acido"]])
    shiny::req(input[["Bact_ScaleFact_Meth"]])
    shiny::req(input[["Bact_ScaleFact_Bact"]])
    shiny::req(input[["WT_Perc_Guar_IN"]])
    shiny::req(input[["WT_Perc_PEG_IN"]])
    shiny::req(input[["wells"]])
    shiny::req(input[["Temp"]])
    shiny::req(input[["Head_Space_VolRatio"]])
    
      test <- PegDist() %>% as.vector() %>% as.numeric()
      test2 <- test*44+18 %>% as.vector() %>% as.numeric()
      hist1 <- graphics::hist(test,breaks=4:9,plot=F)
      yA <- round(hist1$density,3); MW_PEG_In <- mean(test2)
    #if(length(test)>1){
      param(Bact_ScaleFact_Aceto = (as.numeric(input[["Bact_ScaleFact_Aceto"]])/1000),
            Bact_ScaleFact_Acido = (as.numeric(input[["Bact_ScaleFact_Acido"]])/1000),
            Bact_ScaleFact_Meth = (as.numeric(input[["Bact_ScaleFact_Meth"]])/1000),
            Bact_ScaleFact_Bact = (as.numeric(input[["Bact_ScaleFact_Bact"]])/1000),
            WT_Perc_Guar_IN = (as.numeric(input[["WT_Perc_Guar_IN"]]))/100,
            WT_Perc_PEG_IN = (as.numeric(input[["WT_Perc_PEG_IN"]]))/100,
            wells = as.numeric(input[["wells"]]),
            Temp = as.numeric(input[["Temp"]]),
            Head_Space_VolRatio = as.numeric(input[["Head_Space_VolRatio"]]),
            molFracPEG9 = yA[5], molFracPEG8 = yA[4], molFracPEG7 = yA[3],
            molFracPEG6 = yA[2], molFracPEG5 = yA[1], MW_PEG_In = MW_PEG_In,
            Intercept_Aceto_vmax=as.numeric(coef_Aceto_vmax[1]),a_Aceto_vmax=as.numeric(coef_Aceto_vmax[2]),
            b_Aceto_vmax=as.numeric(coef_Aceto_vmax[3]),c_Aceto_vmax=as.numeric(coef_Aceto_vmax[4]),
            Intercept_Aceto_km=as.numeric(coef_Aceto_km[1]),a_Aceto_km=as.numeric(coef_Aceto_km[2]),
            b_Aceto_km=as.numeric(coef_Aceto_km[3]),c_Aceto_km=as.numeric(coef_Aceto_km[4]),
            Intercept_Meth_vmax=as.numeric(coef_Meth_vmax[1]),a_Meth_vmax=as.numeric(coef_Meth_vmax[2]),
            b_Meth_vmax=as.numeric(coef_Meth_vmax[3]),c_Meth_vmax=as.numeric(coef_Meth_vmax[4]),
            Intercept_Meth_km=as.numeric(coef_Meth_km[1]),a_Meth_km=as.numeric(coef_Meth_km[2]),
            b_Meth_km=as.numeric(coef_Meth_km[3]),c_Meth_km=as.numeric(coef_Meth_km[4])
      )
    # }else{
    #   NULL
    # }
  })
  
  SensParam <- reactive({
    shiny::req(input[["simType"]])
    simType <- try({as.character(input[["simType"]])},silent = F)
    if(simType!="Normal"){
      shiny::req(input[["SensParam"]])
      SensParam <- as.character(input[["SensParam"]])
      if(SensParam=="Acidogen Conc."){
        SensParam2 <- "Bact_ScaleFact_Acido"
      }else if(SensParam=="Acetogen Conc."){
        SensParam2 <- "Bact_ScaleFact_Aceto"
      }else if(SensParam=="Methanogen Conc."){
        SensParam2 <- "Bact_ScaleFact_Meth"
      }else if(SensParam=="Bacteroid Conc."){
        SensParam2 <- "Bact_ScaleFact_Bact"
      }else if(SensParam=="Head Space Ratio"){
        SensParam2 <- "Head_Space_VolRatio"
      }else if(SensParam=="Temperature"){
        SensParam2 <- "Temp"
      }else if(SensParam=="Number of Wells"){
        SensParam2 <- "wells"
      }else if(SensParam=="WT % of Guar Gum"){
        SensParam2 <- "WT_Perc_Guar_IN"
      }else if(SensParam=="WT % of PEG-400"){
        SensParam2 <- "WT_Perc_PEG_IN"
      }
      return(SensParam2)
    }else{
      return(NULL)
    }
  })
  
  SensRangeR <- reactive({
    shiny::req(input[["simType"]])
    simType <- try({as.character(input[["simType"]])},silent = F)
    if(simType=="Sensitivity Analysis"){
      shiny::req(input[["SensParam"]])
      shiny::req(input[["SensRange"]])
      SensParam <- SensParam() %>% as.character()
      shiny::req(input[[SensParam]])
      MedVal <- as.numeric(input[[SensParam]])
      SensRange <- as.numeric(input[["SensRange"]])/100
      if(SensParam=="WT_Perc_Guar_IN"|SensParam=="WT_Perc_PEG_IN"){
        MedVal <- MedVal/100
      }else if(SensParam=="Bact_ScaleFact_Acido"|SensParam=="Bact_ScaleFact_Aceto"|SensParam=="Bact_ScaleFact_Meth"|SensParam=="Bact_ScaleFact_Bact"){
        MedVal <- MedVal/1000
      }else{
        MedVal <- MedVal
      }
      SensVals <- sort(c(MedVal,SensRange*MedVal+MedVal))
      SensVals
    }
  })
  
  omega_Kinetics <- reactive({
    shiny::req(input[["kinetic_var"]])
      kinetic_Var <- (as.numeric(input[["kinetic_var"]])/100)^2
      omegaMatrix <- omat(mod)@data$Bact_kinetics %>% as.matrix()
      omegaNames <- unlist(omat(mod)@labels[[1]]) %>% as.character()
      diag(omegaMatrix) <- kinetic_Var
      row.names(omegaMatrix) <- omegaNames
    return(omegaMatrix)
  })
  
  omega_Yields <- reactive({
    shiny::req(input[["yield_var"]])
      yield_Var <- (as.numeric(input[["yield_var"]])/100)^2
      omegaMatrix <- omat(mod)@data$Bact_yields %>% as.matrix()
      omegaNames <- unlist(omat(mod)@labels[[2]]) %>% as.character()
      diag(omegaMatrix) <- yield_Var
      row.names(omegaMatrix) <- omegaNames
    return(omegaMatrix)
  })
  
  TimeSolv <- reactive({
    shiny::req(input[["cutOff"]])
    cutOffTime <- as.numeric(input[["cutOff"]])
    #Times to solve equations
    if(cutOffTime<100 & cutOffTime>=80){
      addVal <- sort(unique(c(0.25,0.5,1:80,seq(80,cutOffTime,2.5))))
      T1 <- tgrid(0,cutOffTime,cutOffTime, add=addVal)
    }else if(cutOffTime<80){
      addVal <- sort(unique(c(0.25,0.5,1:cutOffTime)))
      T1 <- tgrid(0,cutOffTime,cutOffTime, add=addVal)
    }else if(cutOffTime>=100){
      addVal <- sort(unique(c(0.25,0.5,1:70,seq(72.5,82.5,by=2.5),seq(85,100,by=5),seq(100,cutOffTime,10))))
      T1 <- tgrid(0,cutOffTime,cutOffTime, add=addVal)
    }
    T1
  })
  
  out <- reactive({
    shiny::req(input[["cutOff"]])
    shiny::req(input[["nSim"]])
    shiny::req(input[["simType"]])
    shiny::req(TimeSolv())
    shiny::req(omega_Kinetics())
    shiny::req(omega_Yields())
    shiny::req(pars())
    cutOffTime <- as.numeric(input[["cutOff"]])
    nSim <- as.numeric(input[["nSim"]])
    simType <- as.character(input[["simType"]])
    T1 <- TimeSolv()
    # If Normal Simulation
    if(simType=="Normal"){
      outDat <- mod %>% param(pars()) %>% omat(Bact_kinetics=omega_Kinetics(),Bact_yields=omega_Yields()) %>% 
        mrgsim(nid=nSim,tgrid=T1,end=cutOffTime,atol = 1E-50,maxsteps=50000,hmax = 0.01)

    # If Sensitivity Analysis
    }else if(simType=="Sensitivity Analysis"){
      shiny::req(input[["SensParam"]])
      SensParam <- SensParam() %>% as.character()
      SensVals <- SensRangeR() %>% as.numeric()
      #idata_set function
      solveModel <- function(SensParam,T1,cutOffTime,nSim,SensVals){
        SensParam <- qc(.(SensParam)) 
        wrapr::let(
          c(SensParam2=substitute(SensParam)),{ 
            datafile2 <- rlang::expr(expand.ev(SensParam2=sort(rep(SensVals,nSim)))) %>% rlang::eval_tidy()
            rlang::expr(mod %>% idata_set(datafile2) %>% omat(Bact_kinetics=omega_Kinetics(),Bact_yields=omega_Yields()) %>% carry.out(SensParam2) %>% 
                          mrgsim(end=cutOffTime,tgrid=T1,atol = 1E-50,maxsteps=50000,hmax = 0.01)) %>% rlang::eval_tidy()})
      }
      outDat <- solveModel(SensParam,T1,cutOffTime,nSim,SensVals)
    }
    #assign("out",outDat,envir = globalenv())
    outDat
  }) # End output
  
  # observeEvent(c(input$resim1,input$resim2,input$resim3,input$resim4,
  #                input$resim5,input$resim6,input$resim7,input$simType),{
  #                  pars <- pars()
  #                  out <- out()
  #                  assign("out",out,envir = globalenv())
  #                },ignoreNULL = T, ignoreInit = T) # End observeEvent
  
  varyParam <- reactive({
    shiny::req(input[["simType"]])
    simType <- try({as.character(input[["simType"]])},silent = F)
    if(simType=="Sensitivity Analysis"){
      shiny::req(input[["SensParam"]])
      varyParam <- SensParam() %>% as.character()
      }else{
    out <- out()
    varyParam <- try({knobs(out)},silent=T) 
      }
    return(varyParam)
  })
  varyParamClass <- reactive({
    out <- out()
    varyParam <- varyParam() %>% as.character()
    varyParamClass <- try({
      class(varyParam) %>% as.character
      },silent=T)
    return(varyParamClass)
  })
  
  confInterval <- reactive({
    input[["confInterval"]] %>% as.character()
  })
  
  #Standard Deviation Function (For Table)
  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=T,
                        conf.interval=.95, .drop=TRUE) {
    library(plyr)
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)}
    # This does the summary. For each group's data frame, return a vector with N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,.fun = function(xx, col) {
      c(N = length2(xx[[col]], na.rm=na.rm),
        mean = mean(xx[[col]], na.rm=na.rm),
        sd = stats::sd(xx[[col]], na.rm=na.rm))},measurevar)
    datac <- plyr::rename(datac, c("mean" = measurevar)) # Rename the "mean" column
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    return(datac)
  }
 
    
    #### Render Component Data Table
    output$sum1 <- renderTable({
      shiny::req(out())
      shiny::req(input[["cutOff"]])
      shiny::req(input[["simType"]])
      simType <- as.character(input[["simType"]])
      out <- out()
      Inputs_Main <- c("H2O", "GUAR", "GUAR_Conc","Pressure_atm","AVG_PEG_MW")
      Outputs_Main <- c("H2O", "GUAR", "GUAR_Conc", "H2_GAS", "CO2_GAS", "CH4_GAS","H2_LIQ", "CO2_LIQ", "CH4_LIQ")
      varyParamClass <- varyParamClass() %>% as.character()
      varyParam <- varyParam() %>% as.character()#try({knobs(out)},silent=T)
      cutOffTime <- as.numeric(input[["cutOff"]]) #100 #hours (4.167 days)
      if(simType=="Sensitivity Analysis"){
        shiny::req(input[["SensParam"]])
        SumDat <- out %>% dplyr::select(ID,time,varyParam,eval(Inputs_Main),eval(Outputs_Main)) %>% as.data.frame
        if(varyParam=="WT_Perc_Guar_IN"|varyParam=="WT_Perc_PEG_IN"){
          SumDat[varyParam] <- round(SumDat[varyParam]*100,3)
        }
      }else{
        SumDat <- out %>% dplyr::select(ID,time,eval(Inputs_Main),eval(Outputs_Main)) %>% as.data.frame
      }
      SumDat <- SumDat %>% dplyr::filter(time==0 | time==cutOffTime) 
      for(i in 1:dim(SumDat)[1]){
        if(SumDat$time[i]==0){
          row.names(SumDat)[i] <- paste("Input", " ID =",SumDat$ID[i])
        }else if(SumDat$time[i]==cutOffTime){
          row.names(SumDat)[i] <- paste("Output", "ID =",SumDat$ID[i])}
        SumDat$TotalMol[i] <- SumDat$H2O[i]+SumDat$GUAR[i]+SumDat$H2_GAS[i]+SumDat$CO2_GAS[i]+
          SumDat$CH4_GAS[i]+SumDat$H2_LIQ[i]+SumDat$CO2_LIQ[i]+SumDat$CH4_LIQ[i]
      }
      SumDat$GUAR <- round(SumDat$GUAR/6600,3) #divide by n_bonds
      SumDat$H2_GAS <- round(SumDat$H2_GAS,3);SumDat$CO2_GAS <- round(SumDat$CO2_GAS,3)
      ModelErrorPerc <<- round((abs(SumDat$TotalMol[1]-SumDat$TotalMol[2])/SumDat$TotalMol[1])*100,3)
      SumDat2 <- SumDat %>% dplyr::select(everything(),-c(ID,time,H2_LIQ,CO2_LIQ,CH4_LIQ,TotalMol)) %>% 
        dplyr::rename("Guar Gum (mol)"=GUAR,"Guar Gum (g/L)"=GUAR_Conc, "H2 (mol-gas)"=H2_GAS,
                      "CO2 (mol-gas)"=CO2_GAS, "CH4 (mol-gas)"=CH4_GAS,"Total Pressure (atm)"=Pressure_atm,
                      "Average PEG MW (g/mol)"=AVG_PEG_MW)
      SumDat2$H2O <- round(SumDat2$H2O*18/1000,2); SumDat2 <- SumDat2 %>% dplyr::rename("Water (L)"=H2O)
      SumDat2$"Guar Gum (g/L)" <- round(SumDat2$"Guar Gum (g/L)",3)
      SumDat2$"Total Pressure (atm)" <- round(SumDat2$"Total Pressure (atm)",2)
      SumDat2$"Average PEG MW (g/mol)" <- round(SumDat2$"Average PEG MW (g/mol)",2)
      if(simType=="Normal" & n_distinct(out$ID)==1){
        row.names(SumDat2) <- c("Input","Output")}
      SumDat_Chem <- SumDat2
      SumDat_Chem
      },rownames = TRUE)
   
    #### Render Bacteria Data Table
    output$sum2 <- renderTable({
      shiny::req(out())
      shiny::req(input[["cutOff"]])
      shiny::req(input[["simType"]])
      simType <- as.character(input[["simType"]])
      out <- out()
      Inputs_Bact <- c("Conc_ACIDOGEN","Conc_ACETOGEN","Conc_METHANOGEN","Conc_BACTEROID")
      Outputs_Bact <- c("Conc_ACIDOGEN","Conc_ACETOGEN","Conc_METHANOGEN","Conc_BACTEROID")
      varyParamClass <- varyParamClass() %>% as.character()
      varyParam <- varyParam() %>% as.character()
      cutOffTime <- as.numeric(input[["cutOff"]]) #100 #hours (4.167 days)
      
      if(simType=="Sensitivity Analysis"){
        shiny::req(input[["SensParam"]])
        varyParam <- varyParam %>% as.character()
        SumDat <- out %>% dplyr::select(ID,time,varyParam,eval(Inputs_Bact),eval(Outputs_Bact)) %>% as.data.frame
        if(varyParam=="WT_Perc_Guar_IN"|varyParam=="WT_Perc_PEG_IN"){
          SumDat[varyParam] <- round(SumDat[varyParam]*100,3)
        }
      }else{
        SumDat <- out %>% dplyr::select(ID,time,eval(Inputs_Bact),eval(Outputs_Bact)) %>% as.data.frame
      }
      SumDat <- SumDat %>% dplyr::filter(time==0 | time==cutOffTime) 
      for(i in 1:dim(SumDat)[1]){
        if(SumDat$time[i]==0){
          row.names(SumDat)[i] <- paste("Input", " ID =",SumDat$ID[i])
        }else if(SumDat$time[i]==cutOffTime){
          row.names(SumDat)[i] <- paste("Output", "ID =",SumDat$ID[i])}}
      SumDat2 <- SumDat %>% dplyr::select(everything(),-c(ID,time)) 
      SumDat2$Conc_ACIDOGEN <- round(SumDat2$Conc_ACIDOGEN,2); SumDat2$Conc_ACETOGEN <- round(SumDat2$Conc_ACETOGEN,2)
      SumDat2$Conc_METHANOGEN <- round(SumDat2$Conc_METHANOGEN,2); SumDat2$Conc_BACTEROID <- round(SumDat2$Conc_BACTEROID,2)
      SumDat2 <- SumDat2 %>% dplyr::rename("Acidogens (g/L)"=Conc_ACIDOGEN, "Acetogens (g/L)"=Conc_ACETOGEN,
                                           "Methanogens (g/L)"=Conc_METHANOGEN,"Bacteroides (g/L)"=Conc_BACTEROID)
      if(simType=="Normal" & n_distinct(out$ID)==1){
        row.names(SumDat2) <- c("Input","Output")
        }
      SumDat_Bact <- SumDat2
      SumDat_Bact
      },rownames = TRUE)
    
    # Render Model Error
    output$modelError <- renderText({
      ErrorMessage <- paste("Model Error in Mole Balance: <b>",ModelErrorPerc,"% </b>")
      HTML(paste(ErrorMessage))
    })
    
 # }) #end observe
  
  # Plots --------------------------------------------------------------
  
  
  observe(priority=3,{
    shiny::req(out())
    shiny::req(input[["cutOff"]])
    shiny::req(input[["simType"]])
    simType <- as.character(input[["simType"]])
    out <- out() %>% as.data.frame()
    varyParamClass <- varyParamClass() %>% as.character()
    varyParam <- varyParam() %>% as.character()#try({knobs(out)},silent=T)
    cutOffTime <- as.numeric(input[["cutOff"]]) #100 #hours (4.167 days)
    
    #Plot Datasets
    if(simType=="Sensitivity Analysis"){
      shiny::req(input[["SensParam"]])
      varyParam <- SensParam() %>% as.character()
      MainProdDat <- out %>% dplyr::select(ID,time,varyParam,GUAR_Conc,AVG_PEG_MW,CH4_GAS, H2_GAS, CO2_GAS) %>%
        dplyr::filter(time<=cutOffTime)  %>% dplyr::rename("Guar Gum (g/L)"=GUAR_Conc, "AVG PEG MW (g/mol)"=AVG_PEG_MW,
                                                           "H2 (mol-gas)"=H2_GAS,"CO2 (mol-gas)"=CO2_GAS,"CH4 (mol-gas)"=CH4_GAS) %>% as.data.frame()
      IntermediateDat <- out %>% dplyr::select(ID,time,varyParam,GUAR_Conc,Conc_GLUCOSE,Conc_ETHANOL,Conc_PropAcid,Conc_ACETATE,CH4_GAS) %>%
        dplyr::filter(time<=cutOffTime)  %>% dplyr::rename("Guar Gum (g/L)"=GUAR_Conc,"CH4 (mol-gas)"=CH4_GAS, "Glucose (mol/L)"=Conc_GLUCOSE,
                                                           "Ethanol (mol/L)"=Conc_ETHANOL,"Propanoic Acid (mol/L)"=Conc_PropAcid,
                                                           "Acetate (mol/L)"=Conc_ACETATE) %>% as.data.frame()
      PEG_Dat <- out %>% dplyr::select(ID,time,varyParam,Conc_PEG9,Conc_PEG8,Conc_PEG7,Conc_PEG6,Conc_PEG5,Conc_PEG4,Conc_PEG3,Conc_DEG,Conc_EG,Conc_AcetHyde,Conc_ACETATE,CH4_GAS) %>%
        dplyr::filter(time<=cutOffTime)  %>% dplyr::rename("PEG-9 (mol/L)"=Conc_PEG9,"PEG-8 (mol/L)"=Conc_PEG8,"PEG-7 (mol/L)"=Conc_PEG7,
                                                           "PEG-6 (mol/L)"=Conc_PEG6,"PEG-5 (mol/L)"=Conc_PEG5,"PEG-4 (mol/L)"=Conc_PEG4,
                                                           "PEG-3 (mol/L)"=Conc_PEG3,"DEG (mol/L)"=Conc_DEG,"EG (mol/L)"=Conc_EG,
                                                           "Acetaldehyde (mol/L)"=Conc_AcetHyde, "Acetate (mol/L)"=Conc_ACETATE,
                                                           "CH4 (mol-gas)"=CH4_GAS) %>% as.data.frame()
      SystemDat <- out %>% dplyr::select(ID,time,varyParam,H2O,V_TOT,Pressure_atm,Temp2) %>% 
        dplyr::filter(time<=cutOffTime)  %>% dplyr::rename("Total Pressure (atm)"=Pressure_atm, 
                                                           "Liquid Volume (% Change)"=V_TOT,"Water (% Change)"=H2O,"Temperature (C)" = Temp2) %>% as.data.frame()
      SystemDat$"Liquid Volume (% Change)" <- (SystemDat$"Liquid Volume (% Change)"/SystemDat$"Liquid Volume (% Change)"[1])*100
      SystemDat$"Water (% Change)" <- (SystemDat$"Water (% Change)"/SystemDat$"Water (% Change)"[1])*100
      BacteriaDat <- out %>% dplyr::select(ID,time,varyParam,Conc_ACIDOGEN,Conc_ACETOGEN,Conc_METHANOGEN,Conc_BACTEROID) %>% 
        dplyr::filter(time<=cutOffTime)  %>% dplyr::rename("Acidogen Biomass (g/L)"=Conc_ACIDOGEN, "Acetogen Biomass (g/L)"=Conc_ACETOGEN,
                                                           "Methanogen Biomass (g/L)"=Conc_METHANOGEN,"Bacteroid Biomass (g/L)"=Conc_BACTEROID) %>% as.data.frame()
      if(varyParam=="WT_Perc_Guar_IN"|varyParam=="WT_Perc_PEG_IN"){
        MainProdDat[varyParam] <- round(MainProdDat[varyParam]*100,3)
        IntermediateDat[varyParam] <- round(IntermediateDat[varyParam]*100,3)
        PEG_Dat[varyParam] <- round(PEG_Dat[varyParam]*100,3)
        SystemDat[varyParam] <- round(SystemDat[varyParam]*100,3)
        BacteriaDat[varyParam] <- round(BacteriaDat[varyParam]*100,3)
      }
      MainProdDat <<- MainProdDat
      IntermediateDat <<- IntermediateDat
      PEG_Dat <<- PEG_Dat
      SystemDat <<- SystemDat
      BacteriaDat <<- BacteriaDat
    }else{
      MainProdDat <- out %>% dplyr::select(ID,time,GUAR_Conc,AVG_PEG_MW,CH4_GAS, H2_GAS, CO2_GAS) %>%
        dplyr::filter(time<=cutOffTime) %>% dplyr::rename("Guar Gum (g/L)"=GUAR_Conc, "AVG PEG MW (g/mol)"=AVG_PEG_MW,
                                                           "H2 (mol-gas)"=H2_GAS,"CO2 (mol-gas)"=CO2_GAS,"CH4 (mol-gas)"=CH4_GAS) %>% as.data.frame()
      IntermediateDat <- out %>% dplyr::select(ID,time,GUAR_Conc,Conc_GLUCOSE,Conc_ETHANOL,Conc_PropAcid,Conc_ACETATE,CH4_GAS) %>%
        dplyr::filter(time<=cutOffTime) %>% dplyr::rename("Guar Gum (g/L)"=GUAR_Conc,"CH4 (mol-gas)"=CH4_GAS, "Glucose (mol/L)"=Conc_GLUCOSE,
                                                           "Ethanol (mol/L)"=Conc_ETHANOL,"Propanoic Acid (mol/L)"=Conc_PropAcid,
                                                           "Acetate (mol/L)"=Conc_ACETATE) %>% as.data.frame()
      PEG_Dat <- out %>% dplyr::select(ID,time,Conc_PEG9,Conc_PEG8,Conc_PEG7,Conc_PEG6,Conc_PEG5,Conc_PEG4,Conc_PEG3,Conc_DEG,Conc_EG,Conc_AcetHyde,Conc_ACETATE,CH4_GAS) %>%
        dplyr::filter(time<=cutOffTime) %>% dplyr::rename("PEG-9 (mol/L)"=Conc_PEG9,"PEG-8 (mol/L)"=Conc_PEG8,"PEG-7 (mol/L)"=Conc_PEG7,
                                                           "PEG-6 (mol/L)"=Conc_PEG6,"PEG-5 (mol/L)"=Conc_PEG5,"PEG-4 (mol/L)"=Conc_PEG4,
                                                           "PEG-3 (mol/L)"=Conc_PEG3,"DEG (mol/L)"=Conc_DEG,"EG (mol/L)"=Conc_EG,
                                                           "Acetaldehyde (mol/L)"=Conc_AcetHyde, "Acetate (mol/L)"=Conc_ACETATE,
                                                           "CH4 (mol-gas)"=CH4_GAS) %>% as.data.frame()
      SystemDat <- out %>% dplyr::select(ID,time,H2O,V_TOT,Pressure_atm,Temp2) %>% 
        dplyr::filter(time<=cutOffTime)  %>% dplyr::rename("Total Pressure (atm)"=Pressure_atm, 
                                                           "Liquid Volume (% Change)"=V_TOT,"Water (% Change)"=H2O,"Temperature (C)" = Temp2) %>% as.data.frame()
      SystemDat$"Liquid Volume (% Change)" <- (SystemDat$"Liquid Volume (% Change)"/SystemDat$"Liquid Volume (% Change)"[1])*100
      SystemDat$"Water (% Change)" <- (SystemDat$"Water (% Change)"/SystemDat$"Water (% Change)"[1])*100
      BacteriaDat <- out %>% dplyr::select(ID,time,Conc_ACIDOGEN,Conc_ACETOGEN,Conc_METHANOGEN,Conc_BACTEROID) %>% 
        dplyr::filter(time<=cutOffTime)  %>% dplyr::rename("Acidogen Biomass (g/L)"=Conc_ACIDOGEN, "Acetogen Biomass (g/L)"=Conc_ACETOGEN,
                                                           "Methanogen Biomass (g/L)"=Conc_METHANOGEN,"Bacteroid Biomass (g/L)"=Conc_BACTEROID) %>% as.data.frame()
      MainProdDat <<- MainProdDat
      IntermediateDat <<- IntermediateDat
      PEG_Dat <<- PEG_Dat
      SystemDat <<- SystemDat
      BacteriaDat <<- BacteriaDat
    }
  }) #end observe
  
  ### Plot 1 (Reactor Properties)
  output$plot1 <- renderPlot({
    shiny::req(out())
    shiny::req(input[["cutOff"]])
    shiny::req(input[["simType"]])
    simType <- as.character(input[["simType"]])
    #Load Reactive Objects
    varyParamClass <- varyParamClass() %>% as.character()
    varyParam <- varyParam() %>% as.character()
    cutOffTime <- as.numeric(input[["cutOff"]]) #100 #hours (4.167 days)
    confInterval <- confInterval() %>% as.character()
    if(exists("SystemDat")){
    #Sensitivity
    if(simType=="Sensitivity Analysis"){
      shiny::req(input[["SensParam"]])
      varyParam <- varyParam() %>% as.character()
      colNames <- names(SystemDat)[4:length(SystemDat)]
      p1 <- vector("list",length = length(colNames)); names(p1) <- names(SystemDat)[4:length(SystemDat)]
      colScale <- scale_color_discrete(name = eval(varyParam))
      #Make Plot
      stdDevs <- yminVal <- ymaxVal <- yminE <- ymaxE <- p1
      p1 <- lapply(1:length(colNames),function(i){
        stdDevsDat <- summarySE(data=SystemDat,measurevar = colNames[i],groupvars = c("time",varyParam)) %>% as.data.frame()
        yminE <- stdDevsDat[,4] - stdDevsDat[,6]
        ymaxE <- stdDevsDat[,4] + stdDevsDat[,6]
        yminVal <- round(min(yminE)-min(yminE)*.01,15)
        ymaxVal <- round(max(ymaxE)+max(ymaxE)*.015,15)
        yLimits <- c(yminVal, ymaxVal)
        #Make Plot
        if(confInterval=="Error Bars"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,4], group=stdDevsDat[,varyParam])) + geom_line(aes_string(color=factor(stdDevsDat[,varyParam])),size=0.5) +
            geom_errorbar(color="black",aes(ymin=yminE,ymax=ymaxE),width=2.5) + geom_point(aes_string(color=factor(stdDevsDat[,varyParam])),shape=21,size=1.1,alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + colScale + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }else if(confInterval=="Confidence Band"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,4], group=stdDevsDat[,varyParam])) + geom_line(aes_string(color=factor(stdDevsDat[,varyParam])),size=0.5) +
            geom_ribbon(aes_string(ymin=yminE,ymax=ymaxE,fill=factor(stdDevsDat[,varyParam])),alpha=0.3,show.legend = FALSE) + geom_point(aes_string(color=factor(stdDevsDat[,varyParam])),shape=21,size=1.1,alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + colScale + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }
      })
      #Single Simulation 
    }else if(simType=="Normal" & n_distinct(SystemDat$ID)==1){
      colNames <- names(SystemDat)[3:length(SystemDat)]
      p1 <- vector("list",length = length(colNames)); names(p1) <- names(SystemDat)[3:length(SystemDat)]
      #Make Plot
      for(i in colNames){
        if(round(min(SystemDat[,i]))==0){
          yminVal <- 0
        }else{yminVal <- round(min(SystemDat[,i])-min(SystemDat[,i])*.01,15)}
        ymaxVal <- round(max(SystemDat[,i])+max(SystemDat[,i])*.015,15)
          p1[[i]] <- ggplot(SystemDat, aes_string(x=SystemDat$time, y = SystemDat[,i])) + 
            geom_point(color="red") + geom_hline(yintercept=0, size=0.6, color="black") + labs(x="Time (h)",y=i) + 
            geom_vline(xintercept=0, size=0.6, color="black") + scale_y_continuous(limits = c(yminVal, ymaxVal)) +
            theme(text=element_text(family="Times New Roman", face="bold", size=13))}
      #Multiple Simulations
    }else if(simType=="Normal" & n_distinct(SystemDat$ID)!=1){
      colNames <- names(SystemDat)[3:length(SystemDat)]
      p1 <- vector("list",length = length(colNames)); names(p1) <- colNames
      stdDevs <- yminVal <- ymaxVal <- yminE <- ymaxE <- p1
      p1 <- lapply(1:length(colNames),function(i){
        stdDevsDat <- summarySE(data=SystemDat,measurevar = colNames[i],groupvars = c("time")) %>% as.data.frame()
        yminE <- stdDevsDat[,3] - stdDevsDat[,5]
        ymaxE <- stdDevsDat[,3] + stdDevsDat[,5]
        yminVal <- round(min(yminE)-min(yminE)*.01,15)
        ymaxVal <- round(max(ymaxE)+max(ymaxE)*.015,15)
        yLimits <- c(yminVal, ymaxVal)
        #Make Plot
        if(confInterval=="Error Bars"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,3])) + geom_line(color="red",size=0.5) +
            geom_errorbar(color="black",aes(ymin=yminE,ymax=ymaxE),width=2.5) + geom_point(color="black",shape=21,size=1.1,fill="white",alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }else if(confInterval=="Confidence Band"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,3])) + geom_line(color="red",size=0.5) +
            geom_ribbon(aes(ymin=yminE,ymax=ymaxE),alpha=0.3) + geom_point(color="black",shape=21,size=1.1,fill="white",alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }
          })
    }
    p1 <<- p1; Plot1 <- grid.arrange(grobs=p1,ncol=2); assign("Plot1",Plot1,envir = globalenv()) #; print(Plot1)
    grid.arrange(Plot1)
    }
  },res=110)
  
  ### Plot 2 (Main Component Data)
  output$plot2 <- renderPlot({
    shiny::req(out())
    shiny::req(input[["cutOff"]])
    shiny::req(input[["simType"]])
    simType <- as.character(input[["simType"]])
    #Load Reactive Objects
    varyParamClass <- varyParamClass() %>% as.character()
    varyParam <- varyParam() %>% as.character()
    cutOffTime <- as.numeric(input[["cutOff"]]) #100 #hours (4.167 days)
    confInterval <- confInterval() %>% as.character()
    if(exists("MainProdDat")){
    #Sensitivity
    if(simType=="Sensitivity Analysis"){
      colNames <- names(MainProdDat)[4:length(MainProdDat)]
      p2 <- vector("list",length = length(colNames)); names(p2) <- names(MainProdDat)[4:length(MainProdDat)]
      colScale <- scale_color_discrete(name = eval(varyParam))
      #Make Plot
      stdDevs <- yminVal <- ymaxVal <- yminE <- ymaxE <- p2
      p2 <- lapply(1:length(colNames),function(i){
        stdDevsDat <- summarySE(data=MainProdDat,measurevar = colNames[i],groupvars = c("time",varyParam)) %>% as.data.frame()
        yminE <- stdDevsDat[,4] - stdDevsDat[,6]
        ymaxE <- stdDevsDat[,4] + stdDevsDat[,6]
        yminVal <- round(min(yminE)-min(yminE)*.01,15)
        ymaxVal <- round(max(ymaxE)+max(ymaxE)*.015,15)
        yLimits <- c(yminVal, ymaxVal)
        #Make Plot
        if(confInterval=="Error Bars"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,4], group=stdDevsDat[,varyParam])) + geom_line(aes_string(color=factor(stdDevsDat[,varyParam])),size=0.5) +
            geom_errorbar(color="black",aes(ymin=yminE,ymax=ymaxE),width=2.5) + geom_point(aes_string(color=factor(stdDevsDat[,varyParam])),shape=21,size=1.1,alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + colScale + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }else if(confInterval=="Confidence Band"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,4], group=stdDevsDat[,varyParam])) + geom_line(aes_string(color=factor(stdDevsDat[,varyParam])),size=0.5) +
            geom_ribbon(aes_string(ymin=yminE,ymax=ymaxE,fill=factor(stdDevsDat[,varyParam])),alpha=0.3,show.legend = FALSE) + geom_point(aes_string(color=factor(stdDevsDat[,varyParam])),shape=21,size=1.1,alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + colScale + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }
      })
      #Single Simulation 
    }else if(simType=="Normal" & n_distinct(MainProdDat$ID)==1){
      colNames <- names(MainProdDat)[3:length(MainProdDat)]
      p2 <- vector("list",length = length(colNames)); names(p2) <- names(MainProdDat)[3:length(MainProdDat)]
      #Make Plot
      for(i in colNames){
        if(round(min(MainProdDat[,i]))==0){
          yminVal <- 0
        }else{yminVal <- round(min(MainProdDat[,i])-min(MainProdDat[,i])*.01,15)}
        ymaxVal <- round(max(MainProdDat[,i])+max(MainProdDat[,i])*.015,15)
        p2[[i]] <- ggplot(MainProdDat, aes_string(x=MainProdDat$time, y = MainProdDat[,i])) + 
          geom_point(color="red") + geom_hline(yintercept=0, size=0.6, color="black") + labs(x="Time (h)",y=i) + 
          geom_vline(xintercept=0, size=0.6, color="black") + scale_y_continuous(limits = c(yminVal, ymaxVal)) +
          theme(text=element_text(family="Times New Roman", face="bold", size=13))}
      #Multiple Simulations
    }else if(simType=="Normal" & n_distinct(MainProdDat$ID)!=1){
      colNames <- names(MainProdDat)[3:length(MainProdDat)]
      p2 <- vector("list",length = length(colNames)); names(p2) <- colNames
      stdDevs <- yminVal <- ymaxVal <- yminE <- ymaxE <- p2
      p2 <- lapply(1:length(colNames),function(i){
        stdDevsDat <- summarySE(data=MainProdDat,measurevar = colNames[i],groupvars = c("time")) %>% as.data.frame()
        yminE <- stdDevsDat[,3] - stdDevsDat[,5]
        ymaxE <- stdDevsDat[,3] + stdDevsDat[,5]
        yminVal <- round(min(yminE)-min(yminE)*.01,15)
        ymaxVal <- round(max(ymaxE)+max(ymaxE)*.015,15)
        yLimits <- c(yminVal, ymaxVal)
        #Make Plot
        if(confInterval=="Error Bars"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,3])) + geom_line(color="red",size=0.5) +
            geom_errorbar(color="black",aes(ymin=yminE,ymax=ymaxE),width=2.5) + geom_point(color="black",shape=21,size=1.1,fill="white",alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }else if(confInterval=="Confidence Band"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,3])) + geom_line(color="red",size=0.5) +
            geom_ribbon(aes(ymin=yminE,ymax=ymaxE),alpha=0.3) + geom_point(color="black",shape=21,size=1.1,fill="white",alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }
          })
    }
    p2 <<- p2; Plot2 <- grid.arrange(grobs=p2,ncol=2); assign("Plot2",Plot2,envir = globalenv()) #; print(Plot2)
    grid.arrange(Plot2)
    }
  },res=110)
  
  ### Plot 3 (Intermediate Component Data for Guar Gum)
  output$plot3 <- renderPlot({
    shiny::req(out())
    shiny::req(input[["cutOff"]])
    shiny::req(input[["simType"]])
    simType <- as.character(input[["simType"]])
    #Load Reactive Objects
    varyParamClass <- varyParamClass() %>% as.character()
    varyParam <- varyParam() %>% as.character()
    cutOffTime <- as.numeric(input[["cutOff"]]) #100 #hours (4.167 days)
    confInterval <- confInterval() %>% as.character()
    if(exists("IntermediateDat")){
    #Sensitivity
    if(simType=="Sensitivity Analysis"){
      colNames <- names(IntermediateDat)[4:length(IntermediateDat)]
      p3 <- vector("list",length = length(colNames)); names(p3) <- names(IntermediateDat)[4:length(IntermediateDat)]
      colScale <- scale_color_discrete(name = eval(varyParam))
      #Make Plot
      stdDevs <- yminVal <- ymaxVal <- yminE <- ymaxE <- p3
      p3 <- lapply(1:length(colNames),function(i){
        stdDevsDat <- summarySE(data=IntermediateDat,measurevar = colNames[i],groupvars = c("time",varyParam)) %>% as.data.frame()
        yminE <- stdDevsDat[,4] - stdDevsDat[,6]
        ymaxE <- stdDevsDat[,4] + stdDevsDat[,6]
        yminVal <- round(min(yminE)-min(yminE)*.01,15)
        ymaxVal <- round(max(ymaxE)+max(ymaxE)*.015,15)
        yLimits <- c(yminVal, ymaxVal)
        #Make Plot
        if(confInterval=="Error Bars"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,4], group=stdDevsDat[,varyParam])) + geom_line(aes_string(color=factor(stdDevsDat[,varyParam])),size=0.5) +
            geom_errorbar(color="black",aes(ymin=yminE,ymax=ymaxE),width=2.5) + geom_point(aes_string(color=factor(stdDevsDat[,varyParam])),shape=21,size=1.1,alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + colScale + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }else if(confInterval=="Confidence Band"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,4], group=stdDevsDat[,varyParam])) + geom_line(aes_string(color=factor(stdDevsDat[,varyParam])),size=0.5) +
            geom_ribbon(aes_string(ymin=yminE,ymax=ymaxE,fill=factor(stdDevsDat[,varyParam])),alpha=0.3,show.legend = FALSE) + geom_point(aes_string(color=factor(stdDevsDat[,varyParam])),shape=21,size=1.1,alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + colScale + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }
      })
      #Single Simulation 
    }else if(simType=="Normal" & n_distinct(IntermediateDat$ID)==1){
      colNames <- names(IntermediateDat)[3:length(IntermediateDat)]
      p3 <- vector("list",length = length(colNames)); names(p3) <- names(IntermediateDat)[3:length(IntermediateDat)]
      #Make Plot
      for(i in colNames){
        if(round(min(IntermediateDat[,i]))==0){
          yminVal <- 0
        }else{yminVal <- round(min(IntermediateDat[,i])-min(IntermediateDat[,i])*.01,15)}
        ymaxVal <- round(max(IntermediateDat[,i])+max(IntermediateDat[,i])*.015,15)
        p3[[i]] <- ggplot(IntermediateDat, aes_string(x=IntermediateDat$time, y = IntermediateDat[,i])) + 
          geom_point(color="red") + geom_hline(yintercept=0, size=0.6, color="black") + labs(x="Time (h)",y=i) + 
          geom_vline(xintercept=0, size=0.6, color="black") + scale_y_continuous(limits = c(yminVal, ymaxVal)) +
          theme(text=element_text(family="Times New Roman", face="bold", size=13))}
      #Multiple Simulations
    }else if(simType=="Normal" & n_distinct(IntermediateDat$ID)!=1){
      colNames <- names(IntermediateDat)[3:length(IntermediateDat)]
      p3 <- vector("list",length = length(colNames)); names(p3) <- colNames
      stdDevs <- yminVal <- ymaxVal <- yminE <- ymaxE <- p3
      p3 <- lapply(1:length(colNames),function(i){
        stdDevsDat <- summarySE(data=IntermediateDat,measurevar = colNames[i],groupvars = c("time")) %>% as.data.frame()
        yminE <- stdDevsDat[,3] - stdDevsDat[,5]
        ymaxE <- stdDevsDat[,3] + stdDevsDat[,5]
        yminVal <- round(min(yminE)-min(yminE)*.01,15)
        ymaxVal <- round(max(ymaxE)+max(ymaxE)*.015,15)
        yLimits <- c(yminVal, ymaxVal)
        #Make Plot
        if(confInterval=="Error Bars"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,3])) + geom_line(color="red",size=0.5) +
            geom_errorbar(color="black",aes(ymin=yminE,ymax=ymaxE),width=2.5) + geom_point(color="black",shape=21,size=1.1,fill="white",alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }else if(confInterval=="Confidence Band"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,3])) + geom_line(color="red",size=0.5) +
            geom_ribbon(aes(ymin=yminE,ymax=ymaxE),alpha=0.3) + geom_point(color="black",shape=21,size=1.1,fill="white",alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }
        })
    }
    p3 <<- p3; Plot3 <- grid.arrange(grobs=p3,ncol=2); assign("Plot3",Plot3,envir = globalenv()) #; print(Plot3)
    grid.arrange(Plot3)
    }
  },res=110)
  
  ### Plot 4 (Intermediate Component Data for PEG)
  output$plot4 <- renderPlot({
    shiny::req(out())
    shiny::req(input[["cutOff"]])
    shiny::req(input[["simType"]])
    simType <- as.character(input[["simType"]])
    #Load Reactive Objects
    varyParamClass <- varyParamClass() %>% as.character()
    varyParam <- varyParam() %>% as.character()
    cutOffTime <- as.numeric(input[["cutOff"]]) #100 #hours (4.167 days)
    confInterval <- confInterval() %>% as.character()
    if(exists("PEG_Dat")){
    #Sensitivity
    if(simType=="Sensitivity Analysis"){
      colNames <- names(PEG_Dat)[4:length(PEG_Dat)]
      p4 <- vector("list",length = length(colNames)); names(p4) <- names(PEG_Dat)[4:length(PEG_Dat)]
      colScale <- scale_color_discrete(name = eval(varyParam))
      #Make Plot
      stdDevs <- yminVal <- ymaxVal <- yminE <- ymaxE <- p4
      p4 <- lapply(1:length(colNames),function(i){
        stdDevsDat <- summarySE(data=PEG_Dat,measurevar = colNames[i],groupvars = c("time",varyParam)) %>% as.data.frame()
        yminE <- stdDevsDat[,4] - stdDevsDat[,6]
        ymaxE <- stdDevsDat[,4] + stdDevsDat[,6]
        yminVal <- round(min(yminE)-min(yminE)*.01,15)
        ymaxVal <- round(max(ymaxE)+max(ymaxE)*.015,15)
        yLimits <- c(yminVal, ymaxVal)
        #Make Plot
        if(confInterval=="Error Bars"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,4], group=stdDevsDat[,varyParam])) + geom_line(aes_string(color=factor(stdDevsDat[,varyParam])),size=0.5) +
            geom_errorbar(color="black",aes(ymin=yminE,ymax=ymaxE),width=2.5) + geom_point(aes_string(color=factor(stdDevsDat[,varyParam])),shape=21,size=1.1,alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + colScale + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }else if(confInterval=="Confidence Band"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,4], group=stdDevsDat[,varyParam])) + geom_line(aes_string(color=factor(stdDevsDat[,varyParam])),size=0.5) +
            geom_ribbon(aes_string(ymin=yminE,ymax=ymaxE,fill=factor(stdDevsDat[,varyParam])),alpha=0.3,show.legend = FALSE) + geom_point(aes_string(color=factor(stdDevsDat[,varyParam])),shape=21,size=1.1,alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + colScale + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }
      })
      #Single Simulation 
    }else if(simType=="Normal" & n_distinct(PEG_Dat$ID)==1){
      colNames <- names(PEG_Dat)[3:length(PEG_Dat)]
      p4 <- vector("list",length = length(colNames)); names(p4) <- names(PEG_Dat)[3:length(PEG_Dat)]
      #Make Plot
      for(i in colNames){
        if(round(min(PEG_Dat[,i]))==0){
          yminVal <- 0
        }else{yminVal <- round(min(PEG_Dat[,i])-min(PEG_Dat[,i])*.01,15)}
        ymaxVal <- round(max(PEG_Dat[,i])+max(PEG_Dat[,i])*.015,15)
        p4[[i]] <- ggplot(PEG_Dat, aes_string(x=PEG_Dat$time, y = PEG_Dat[,i])) + 
          geom_point(color="red") + geom_hline(yintercept=0, size=0.6, color="black") + labs(x="Time (h)",y=i) + 
          geom_vline(xintercept=0, size=0.6, color="black") + scale_y_continuous(limits = c(yminVal, ymaxVal)) +
          theme(text=element_text(family="Times New Roman", face="bold", size=13))}
      #Multiple Simulations
    }else if(simType=="Normal" & n_distinct(PEG_Dat$ID)!=1){
        PEG_Dat1 <- PEG_Dat %>% dplyr::select(everything(),-c("CH4 (mol-gas)")) %>% tidyr::gather(key=key,value=value,-ID,-time)
        PEG_Dat2 <- MainProdDat %>% dplyr::select(ID,time,"AVG PEG MW (g/mol)")
        stdDevsDat <- summarySE(data=PEG_Dat1,measurevar = "value",groupvars = c("key","time")) %>% as.data.frame()
        stdDevsDat2 <- summarySE(data=PEG_Dat2,measurevar = "AVG PEG MW (g/mol)",groupvars = c("time")) %>% as.data.frame()
        yminE <- stdDevsDat[,4] - stdDevsDat[,6]; yminE2 <- stdDevsDat2[,3] - stdDevsDat2[,5]
        ymaxE <- stdDevsDat[,4] + stdDevsDat[,6]; ymaxE2 <- stdDevsDat2[,3] + stdDevsDat2[,5]
        yminVal <- round(min(yminE)-min(yminE)*.01,15); yminVal2 <- round(min(yminE2)-min(yminE2)*.01,15)
        ymaxVal <- round(max(ymaxE)+max(ymaxE)*.015,15); ymaxVal2 <- round(max(ymaxE2)+max(ymaxE2)*.015,15)
        yLimits <- c(yminVal, ymaxVal); yLimits2 <- c(yminVal2, ymaxVal2)
        plotPoints <- unique(PEG_Dat1$key)
        p4 <- vector("list",length = 2); names(p4) <- c("Intermediate Products","AVG MW Weight")
      #Make Plot
      if(confInterval=="Error Bars"){
        p4[[1]] <- ggplot(stdDevsDat, aes(y=value, x=time, group=key)) + scale_y_continuous(limits = yLimits) + 
          geom_errorbar(color="black",aes(ymin=yminE,ymax=ymaxE),width=2.5) + 
          geom_line(aes(color=key), size=1) + geom_hline(yintercept=0, size=0.6, color="black") + 
          labs(x="Time (h)",y="Concentration (mol/L)") + geom_vline(xintercept=0, size=0.6, color="black") + 
          theme(text=element_text(family="Times New Roman", face="bold", size=13)) +
          scale_colour_discrete(name="Polyethylene Glycol (PEG) \nIntermediate Product")
        p4[[2]] <- ggplot(stdDevsDat2, aes_string(y=stdDevsDat2$"AVG PEG MW (g/mol)", x=stdDevsDat2$time)) +
          geom_errorbar(color="black",aes(ymin=yminE2,ymax=ymaxE2),width=2.5) + scale_y_continuous(limits = yLimits2) +
          geom_line(size=1) + geom_hline(yintercept=0, size=0.6, color="black") + 
          labs(x="Time (h)",y="Average PEG MW") + geom_vline(xintercept=0, size=0.6, color="black") + 
          theme(text=element_text(family="Times New Roman", face="bold", size=13))
        
      }else if(confInterval=="Confidence Band"){
        p4[[1]] <- ggplot(stdDevsDat, aes(y=value, x=time, group=key)) + scale_y_continuous(limits = yLimits) + 
          geom_ribbon(aes_string(ymin=yminE,ymax=ymaxE,fill=factor(stdDevsDat$key)),alpha=0.3,show.legend = FALSE) +
          geom_line(aes(color=key), size=1) + geom_hline(yintercept=0, size=0.6, color="black") +  
          labs(x="Time (h)",y="Concentration (mol/L)") + geom_vline(xintercept=0, size=0.6, color="black") + 
          theme(text=element_text(family="Times New Roman", face="bold", size=13)) +
          scale_colour_discrete(name="Polyethylene Glycol (PEG) \nIntermediate Product")
        p4[[2]] <- ggplot(stdDevsDat2, aes_string(y=stdDevsDat2$"AVG PEG MW (g/mol)", x=stdDevsDat2$time)) + 
          geom_ribbon(aes_string(ymin=yminE2,ymax=ymaxE2),alpha=0.3,show.legend = FALSE) + scale_y_continuous(limits = yLimits2) +
          geom_line(size=1) + geom_hline(yintercept=0, size=0.6, color="black") +
          labs(x="Time (h)",y="Average PEG MW") + geom_vline(xintercept=0, size=0.6, color="black") + 
          theme(text=element_text(family="Times New Roman", face="bold", size=13))
      }
        #p4[[1]] + geom_point(data=subset(PEG_Dat3[PEG_Dat3$key==plotPoints[1],]), color='black', shape=18, size=2)
    }
    p4 <<- p4;
    if(simType=="Normal" & n_distinct(PEG_Dat$ID)!=1){
      Plot4 <- grid.arrange(grobs=p4,ncol=1,heights=c(2,1.5)); assign("Plot4",Plot4,envir = globalenv())
    }else{
      Plot4 <- grid.arrange(grobs=p4,ncol=2); assign("Plot4",Plot4,envir = globalenv())
    }
    grid.arrange(Plot4)
    }
  },res=110)
  
  
  
  ### Plot 5 (Bacteria Data)
  output$plot5 <- renderPlot({
    shiny::req(out())
    shiny::req(input[["cutOff"]])
    shiny::req(input[["simType"]])
    simType <- as.character(input[["simType"]])
    #Load Reactive Objects
    varyParamClass <- varyParamClass() %>% as.character()
    varyParam <- varyParam() %>% as.character()
    cutOffTime <- as.numeric(input[["cutOff"]]) #100 #hours (4.167 days)
    confInterval <- confInterval() %>% as.character()
    if(exists("BacteriaDat")){
    #Sensitivity
    if(simType=="Sensitivity Analysis"){
      colNames <- names(BacteriaDat)[4:length(BacteriaDat)]
      p5 <- vector("list",length = length(colNames)); names(p5) <- names(BacteriaDat)[4:length(BacteriaDat)]
      colScale <- scale_color_discrete(name = eval(varyParam))
      #Make Plot
      stdDevs <- yminVal <- ymaxVal <- yminE <- ymaxE <- p5
      p5 <- lapply(1:length(colNames),function(i){
        stdDevsDat <- summarySE(data=BacteriaDat,measurevar = colNames[i],groupvars = c("time",varyParam)) %>% as.data.frame()
        yminE <- stdDevsDat[,4] - stdDevsDat[,6]
        ymaxE <- stdDevsDat[,4] + stdDevsDat[,6]
        yminVal <- round(min(yminE)-min(yminE)*.01,15)
        ymaxVal <- round(max(ymaxE)+max(ymaxE)*.015,15)
        yLimits <- c(yminVal, ymaxVal)
        #Make Plot
        if(confInterval=="Error Bars"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,4], group=stdDevsDat[,varyParam])) + geom_line(aes_string(color=factor(stdDevsDat[,varyParam])),size=0.5) +
            geom_errorbar(color="black",aes(ymin=yminE,ymax=ymaxE),width=2.5) + geom_point(aes_string(color=factor(stdDevsDat[,varyParam])),shape=21,size=1.1,alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + colScale + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }else if(confInterval=="Confidence Band"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,4], group=stdDevsDat[,varyParam])) + geom_line(aes_string(color=factor(stdDevsDat[,varyParam])),size=0.5) +
            geom_ribbon(aes_string(ymin=yminE,ymax=ymaxE,fill=factor(stdDevsDat[,varyParam])),alpha=0.3,show.legend = FALSE) + geom_point(aes_string(color=factor(stdDevsDat[,varyParam])),shape=21,size=1.1,alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + colScale + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }
      })
      #Single Simulation 
    }else if(simType=="Normal" & n_distinct(BacteriaDat$ID)==1){
      colNames <- names(BacteriaDat)[3:length(BacteriaDat)]
      p5 <- vector("list",length = length(colNames)); names(p5) <- names(BacteriaDat)[3:length(BacteriaDat)]
      #Make Plot
      for(i in colNames){
        if(round(min(BacteriaDat[,i]))==0){
          yminVal <- 0
        }else{yminVal <- round(min(BacteriaDat[,i])-min(BacteriaDat[,i])*.01,15)}
        ymaxVal <- round(max(BacteriaDat[,i])+max(BacteriaDat[,i])*.015,15)
        p5[[i]] <- ggplot(BacteriaDat, aes_string(x=BacteriaDat$time, y = BacteriaDat[,i])) + 
          geom_point(color="red") + geom_hline(yintercept=0, size=0.6, color="black") + labs(x="Time (h)",y=i) + 
          geom_vline(xintercept=0, size=0.6, color="black") + scale_y_continuous(limits = c(yminVal, ymaxVal)) +
          theme(text=element_text(family="Times New Roman", face="bold", size=13))}
      #Multiple Simulations
    }else if(simType=="Normal" & n_distinct(BacteriaDat$ID)!=1){
      colNames <- names(BacteriaDat)[3:length(BacteriaDat)]
      p5 <- vector("list",length = length(colNames)); names(p5) <- colNames
      stdDevs <- yminVal <- ymaxVal <- yminE <- ymaxE <- p5
      p5 <- lapply(1:length(colNames),function(i){
        stdDevsDat <- summarySE(data=BacteriaDat,measurevar = colNames[i],groupvars = c("time")) %>% as.data.frame()
        yminE <- stdDevsDat[,3] - stdDevsDat[,5]
        ymaxE <- stdDevsDat[,3] + stdDevsDat[,5]
        yminVal <- round(min(yminE)-min(yminE)*.01,15)
        ymaxVal <- round(max(ymaxE)+max(ymaxE)*.015,15)
        yLimits <- c(yminVal, ymaxVal)
        #Make Plot
        if(confInterval=="Error Bars"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,3])) + geom_line(color="red",size=0.5) +
            geom_errorbar(color="black",aes(ymin=yminE,ymax=ymaxE),width=2.5) + geom_point(color="black",shape=21,size=1.1,fill="white",alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }else if(confInterval=="Confidence Band"){
          ggplot(stdDevsDat, aes_string(x=stdDevsDat$time, y = stdDevsDat[,3])) + geom_line(color="red",size=0.5) +
            geom_ribbon(aes(ymin=yminE,ymax=ymaxE),alpha=0.3) + geom_point(color="black",shape=21,size=1.1,fill="white",alpha=0.9) +
            scale_y_continuous(limits = yLimits) + geom_hline(yintercept=0, size=0.6, color="black") + 
            labs(x="Time (h)",y=colNames[i]) + geom_vline(xintercept=0, size=0.6, color="black") + 
            theme(text=element_text(family="Times New Roman", face="bold", size=13)) 
        }
        })
    }
    p5 <<- p5; Plot5 <- grid.arrange(grobs=p5,ncol=2); assign("Plot5",Plot5,envir = globalenv()) #; print(Plot5)
    grid.arrange(Plot5)
    }
  },res=110)
  
  
  #Render ReadMe UI
  rmdfiles <- c("MathModel.Rmd","GettingStarted.Rmd")
  sapply(rmdfiles, knit, quiet = T)
  output$markdownRM <- renderUI({
    fluidPage(
    #tags$head(HTML("<script type='text/x-mathjax-config'>MathJax.Hub.Config({ TeX: { equationNumbers: {autoNumber: 'all'} } });</script>")),
    #includeMarkdown("MathModel.Rmd")
    #includeMarkdown(rmarkdown::render("MathModel.Rmd",'html_document'))
    #HTML(markdown::markdownToHTML(knit('MathModel.Rmd', quiet = TRUE), fragment.only=TRUE))
    withMathJax(includeMarkdown("MathModel.Rmd"))
    #withMathJax(HTML(markdown::markdownToHTML(knit("MathModel.Rmd", quiet = TRUE), fragment.only=TRUE,title="Frack Off",header="Mathematical Model for Anaerobic Digestion")))
    )
    })
  #Render GettingStarted UI
  output$markdownGS <- renderUI({
    fluidPage(
      #tags$head(HTML("<script type='text/x-mathjax-config'>MathJax.Hub.Config({ TeX: { equationNumbers: {autoNumber: 'all'} } });</script>")),
      #includeMarkdown("MathModel.Rmd")
      #includeMarkdown(rmarkdown::render("MathModel.Rmd",'html_document'))
      #HTML(markdown::markdownToHTML(knit('MathModel.Rmd', quiet = TRUE), fragment.only=TRUE))
      withMathJax(includeMarkdown("GettingStarted.Rmd"))
      #withMathJax(HTML(markdown::markdownToHTML(knit("MathModel.Rmd", quiet = TRUE), fragment.only=TRUE,title="Frack Off",header="Mathematical Model for Anaerobic Digestion")))
    )
  })
  
  
} # End Server


# UI Code -----------------------------------------------------------------

ui <- dashboardPage(#theme = shinytheme("slate"),
  #shinythemes::themeSelector(),
  dashboardHeader(title = "Anaerobic Digestion"),
  dashboardSidebar(
    #hr(),
    sidebarMenu(id="tabs", 
                menuItem("Getting Started", tabName = "GettingStarted", icon = icon("dashboard"),selected=TRUE),
                menuItem("Model Parameters", icon = icon("balance-scale"),
                         shinyWidgets::dropdown(
                           shinydashboard::box(width = 12, status = "primary", solidHeader = FALSE,
                                               h2("Initial Bacteria Concentrations (mg/L)", align = "center"),
                                               numericInput("Bact_ScaleFact_Acido", label=p("Acidogens",style="color:#080808"), min=1,step=1,value=300,max=10000),
                                               numericInput("Bact_ScaleFact_Aceto", label=p("Acetogens",style="color:#080808"), min=1,step=1,value=3000,max=10000),
                                               numericInput("Bact_ScaleFact_Meth", label=p("Methanogens",style="color:#080808"), min=1,step=1,value=3000,max=10000),
                                               numericInput("Bact_ScaleFact_Bact", label=p("Bacteroides",style="color:#080808"), min=1,step=1,value=300,max=10000)
                           ),
                           label = "Bacteria", style = "stretch", size="sm",#up=TRUE,
                           status = "primary", width = "420px",
                           #tooltip = tooltipOptions(title = "Click to see options",placement = "bottom"),
                           animate = animateOptions(
                             enter = animations$fading_entrances$fadeInLeftBig,
                             exit = animations$fading_exits$fadeOutRightBig)
                         ),
                         shinyWidgets::dropdown(
                           shinydashboard::box(width = 12, status = "primary", solidHeader = FALSE,
                                               h2("Initial Weight Percents (%)", align = "center"),
                                               numericInput("WT_Perc_Guar_IN", 
                                                            label=p("Guar Gum",style="color:#080808"), min=0.1,step=.01,value=0.83,max=2),
                                               numericInput("WT_Perc_PEG_IN", 
                                                            label=p("PEG 400",style="color:#080808"), min=0.1,step=.01,value=0.4,max=2),
                                               shinyWidgets::dropdown(
                                                 shinydashboard::box(width = 12, status = "primary", solidHeader = FALSE,
                                                                     h2("PEG MW Distribution", align = "center"),
                                                                     plotOutput("PegPlot",width="100%",height="400px"),
                                                                     p("Note: The distribution was purposely skewed to more closely match",style="color:#080808"),
                                                                     p("observations in the literature.",style="color:#080808"),
                                                                     noUiSliderInput("PegMeanMW",label= p("Adjust Mean PEG MW",style="color:#080808"), 300,400,400,1,tooltips = F)
                                                 ),
                                                 label = "Adjust PEG MW Distribution", style = "stretch",size="sm",
                                                 status = "primary", width = "600px",
                                                 #tooltip = tooltipOptions(title = "Click to see options",placement = "bottom"),
                                                 animate = animateOptions(enter = "fadeInDown", exit = "fadeOutUp",duration = 0.8)
                                               )
                           ),
                           label = "Chemical Compounds", style = "stretch",size="sm", #up=TRUE,
                           status = "primary", width = "420px",
                           #tooltip = tooltipOptions(title = "Click to see options",placement = "bottom"),
                           animate = animateOptions(
                             enter = animations$fading_entrances$fadeInLeftBig,
                             exit = animations$fading_exits$fadeOutRightBig)
                         ),
                         shinyWidgets::dropdown(
                           shinydashboard::box(width = 12, status = "primary", solidHeader = FALSE,
                                               h2("System Properties", align = "center"),
                                               sliderInput("wells",label= p("Number of Wells to Pull From",style="color:#080808"), 1,55,55,1),
                                               sliderInput("Temp",label= p("Reactor Temperature (C)",style="color:#080808"), 20,40,30,1),
                                               numericInput("Head_Space_VolRatio",
                                                            label=p("Ratio of Headspace to Reactor Volume (L/L)",style="color:#080808"),
                                                            min=0.25,max=3,value=2,step=0.05)
                           ),
                           label = "System Properties", style = "stretch",size="sm", #up=TRUE,
                           status = "primary", width = "420px",
                           #tooltip = tooltipOptions(title = "Click to see options",placement = "bottom"),
                           animate = animateOptions(
                             enter = animations$fading_entrances$fadeInLeftBig,
                             exit = animations$fading_exits$fadeOutRightBig)
                         ),
                         shinyWidgets::dropdown(
                           shinydashboard::box(width = 12, status = "primary", solidHeader = FALSE,
                                               h2("Parameter Variability (CV %)", align = "center"),
                                               sliderInput("kinetic_var",label= p("Kinetic Rates",style="color:#080808"), 0,100,25,1),
                                               sliderInput("yield_var",label= p("Bacteria Yields",style="color:#080808"), 0,50,10,1),
                                               radioGroupButtons("confInterval",label=p("Confidence Interval (95%)",style="color:#080808"),justified = TRUE,
                                                                 checkIcon = list(yes = icon("ok", lib = "glyphicon")),
                                                                 choices = c("Confidence Band","Error Bars"),selected = "Confidence Band")
                           ),
                           label = "Variability", style = "stretch",size="sm", #up=TRUE,
                           status = "primary", width = "420px",
                           #tooltip = tooltipOptions(title = "Click to see options",placement = "bottom"),
                           animate = animateOptions(
                             enter = animations$fading_entrances$fadeInLeftBig,
                             exit = animations$fading_exits$fadeOutRightBig)
                         ),
                         shinyWidgets::dropdown(
                           shinydashboard::box(width = 12, status = "primary", solidHeader = FALSE,
                                               h2("Simulation Specification", align = "center"),
                                               radioGroupButtons("simType",label=p("Type of Simulation",style="color:#080808"),justified = TRUE,
                                                                 checkIcon = list(yes = icon("ok", lib = "glyphicon")),
                                                                 choices = c("Normal","Sensitivity Analysis"),selected = "Normal"),
                                               numericInput("nSim", label=p("Number of Simulations",style="color:#080808"), min=1,step=1,value=10,max=100),
                                               conditionalPanel("input.simType=='Sensitivity Analysis'",
                                                                h2("Sensitivity Analysis", align = "center"),
                                                                pickerInput(inputId = "SensParam",label = p("Choose Parameter to Vary",style="color:#080808"), 
                                                                            choices = SensChoices,options = list(size = 5,`style` = "btn-info"),selected = "WT % of Guar Gum"),
                                                                sliderInput("SensRange", p("Range of Parameter Variance:",style="color:#080808"), min = -150, max = 150, post="%",value = c(-50,50))
                                               )
                           ),
                           label = "Simulation", style = "stretch",size="sm", #up=TRUE,
                           status = "primary", width = "420px",
                           #tooltip = tooltipOptions(title = "Click to see options",placement = "bottom"),
                           animate = animateOptions(
                             enter = animations$fading_entrances$fadeInLeftBig,
                             exit = animations$fading_exits$fadeOutRightBig)
                         )
                ),
                menuItem("Plots",  icon = icon("line-chart"), 
                         menuSubItem("Main Components", tabName = "plotTab2", icon = icon("angle-right")), 
                         menuSubItem("Guar Gum Intermediates", tabName = "plotTab3", icon = icon("angle-right")),
                         menuSubItem("PEG Intermediates", tabName = "plotTab4", icon = icon("angle-right")),
                         menuSubItem("Bacteria Growth", tabName = "plotTab5", icon = icon("angle-right")),
                         menuSubItem("Reactor Properties", tabName = "plotTab1", icon = icon("angle-right")) 
                ),
                menuItem("Tables", icon = icon("table"), #tabName="tabtable",
                         menuSubItem("Main Component Data", tabName = "sumTab1", icon = icon("angle-right")),
                         menuSubItem("Bacteria Data", tabName = "sumTab2", icon = icon("angle-right"))
                ), hr(),
                menuItem("Mathematical Model", tabName = "readme", icon = icon("mortar-board")),
                menuItem("Background of Process", tabName = "Background", icon = icon("book")),
                menuItem("Codes",  icon = icon("file-text-o"),
                         menuSubItem("Model Code", tabName = "ModelCode", icon = icon("angle-right")),
                         #menuSubItem("ui.R", tabName = "ui", icon = icon("angle-right")),
                         #menuSubItem("server.R", tabName = "server", icon = icon("angle-right")),
                         menuSubItem("app.R", tabName = "app", icon = icon("angle-right"))
                )
    ),
    hr(),
    # conditionalPanel("input.tabs=='plot1' | input.tabs=='plot2' | input.tabs=='plot3' | input.tabs=='plot4'",
    sidebarMenu(
      menuItem("Data Truncation", icon = icon("chevron-circle-right"),
               fluidRow(column(1),
                        column(10,
                               sliderInput("cutOff",label="Set Time (h) to Truncate Data and Evaluate Output",
                                     25,300,100,5)
                 )))
      ),
    hr(),
    sidebarMenu(
      br(),
      div(img(src="FrackOff.jpg",height=301,width=180),style="text-align: center;")
    )   
    #)
  ), # End dashboardSidebar
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "GettingStarted",
              fluidPage(
                tags$head(HTML("<script type='text/x-mathjax-config'>MathJax.Hub.Config({ TeX: { equationNumbers: {autoNumber: 'all'} } });</script>")),
                withMathJax(),
                uiOutput('markdownGS')
              )
      ),
      # tabItem(tabName = "Input",
      #         ),
      tabItem(tabName = "plotTab1",
              fluidRow( 
                shinydashboard::box(width = 10, status = "primary",
                                    plotOutput("plot1",width="100%",height="800px")
                )#,
                # shinydashboard::box(width = 2, status = "primary",
                # actionBttn(inputId = "resim1",label = "Resimulate",style = "fill",color="primary",size = "lg",block = T)
              )#)
      ),
      tabItem(tabName = "plotTab2",
              fluidRow( 
                shinydashboard::box(width = 10, status = "primary",
                                    plotOutput("plot2",width="100%",height="800px")
                )#,
                # shinydashboard::box(width = 2, status = "primary",
                # actionBttn(inputId = "resim2",label = "Resimulate",style = "fill",color="primary",size = "lg",block = T)
              )#)
      ),         
      tabItem(tabName = "plotTab3",
              fluidRow( 
                shinydashboard::box(width = 10, status = "primary",
                                    plotOutput("plot3",width="100%",height="800px")
                )#,
                # shinydashboard::box(width = 2, status = "primary",
                # actionBttn(inputId = "resim3",label = "Resimulate",style = "fill",color="primary",size = "lg",block = T)
              )#)
      ),  
      tabItem(tabName = "plotTab4",
              fluidRow( 
                shinydashboard::box(width = 10, status = "primary",
                                    plotOutput("plot4",width="100%",height="1050px")
                )#,
                # shinydashboard::box(width = 2, status = "primary",
                # actionBttn(inputId = "resim4",label = "Resimulate",style = "fill",color="primary",size = "lg",block = T)
              )#)
      ),  
      tabItem(tabName = "plotTab5",
              fluidRow( 
                shinydashboard::box(width = 10, status = "primary",
                                    plotOutput("plot5",width="100%",height="800px")
                )#,
                # shinydashboard::box(width = 2, status = "primary",
                # actionBttn(inputId = "resim5",label = "Resimulate",style = "fill",color="primary",size = "lg",block = T)
              )#)
      ),         
      tabItem(tabName = "sumTab1",
              fluidRow(
              shinydashboard::box(width = 10, status = "primary", solidHeader = TRUE, title="Main Component Data",                
                                  tableOutput("sum1"),
                                  htmlOutput(outputId = "modelError")
              )#,
              # shinydashboard::box(width = 2, status = "primary",
              # actionBttn(inputId = "resim6",label = "Resimulate",style = "fill",color="primary",size = "lg",block = T)
              )#)
      ),
      tabItem(tabName = "sumTab2",
              fluidRow(
              shinydashboard::box(width = 10, status = "primary", solidHeader = TRUE, title="Bacteria Data",                
                                  tableOutput("sum2")#,
                                  #htmlOutput(outputId = "modelError")
              )#,
              # shinydashboard::box(width = 2, status = "primary",
              # actionBttn(inputId = "resim7",label = "Resimulate",style = "fill",color="primary",size = "lg",block = T)
              )#)
      ),
      # tabItem(tabName = "ui",
      #         shinydashboard::box(width = NULL, status = "primary", solidHeader = TRUE, title="ui.R",
      #                             downloadButton('downloadData2', 'Download'),
      #                             br(),br(),
      #                             pre(includeText("ui.R"))
      #         )
      # ),
      # tabItem(tabName = "server",
      #         shinydashboard::box(width = NULL, status = "primary", solidHeader = TRUE, title="server.R",
      #                             downloadButton('downloadData3', 'Download'),
      #                             br(),br(),
      #                             pre(includeText("server.R"))
      #         )
      # ),
      tabItem(tabName = "app",
              shinydashboard::box(width = NULL, status = "primary", solidHeader = TRUE, title="app.R",
                                  downloadButton('downloadData4', 'Download'),
                                  br(),br(),
                                  pre(includeText("app.R"))
              )
      ),
      tabItem(tabName = "ModelCode", 
              aceEditor("mod_code", .model,mode="r", height="1000px")
      ),
      tabItem(tabName = "readme",
              fluidPage(
                tags$head(HTML("<script type='text/x-mathjax-config'>MathJax.Hub.Config({ TeX: { equationNumbers: {autoNumber: 'all'} } });</script>")),
                withMathJax(), 
                uiOutput('markdownRM')
              )
      ),
      tabItem(tabName = "Background", 
              tags$iframe(style="height:820px; width:100%; scrolling=yes", 
                          src="AnaerobicDigestion_ReadMe.pdf")
              )
    )#,
    # conditionalPanel("input.tabs=='plotTab1' | input.tabs=='plotTab2' | input.tabs=='plotTab3' | input.tabs=='plotTab4'| input.tabs=='sumTab1' | input.tabs=='sumTab2'",
    #   tags$h2("Adjust parameters in dropdown menus:"),
    #   fluidPage(
    #    ) #End fluidPage
    # ) #End Conditional Panel
  ), # End dashboardBody
  tags$head(
    tags$script(src = "js/session.js"),
    tags$script(src = "js/modal_vid.js"),
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    tags$style(HTML(".main-sidebar { font-size: 16px; }"))
  )
    
) # End UI  

shinyApp(ui = ui, server = server)
#runApp(shinyApp(ui = ui, server = server),launch.browser=TRUE)
