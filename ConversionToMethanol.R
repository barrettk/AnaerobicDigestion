library(magrittr)

#### Summary of Data ####
out <- read.csv(paste0(getwd(),"/data/result.csv"))
Inputs <- out[out$time==0,] %>% dplyr::select(CO2_GAS,CO2_LIQ,H2_GAS,H2_LIQ,CH4_GAS,CH4_LIQ) %>%
  summary()

Outputs <- out[out$time==max(out$time),] %>% dplyr::select(CO2_GAS,CO2_LIQ,H2_GAS,H2_LIQ,CH4_GAS,CH4_LIQ) %>%
  summary()

Numextract <- function(string){
  unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
}
Inputs <- Numextract(Inputs[4,]) %>% as.data.frame(); names(Inputs) <- "Inlet"
Inputs$Component <- row.names(Inputs)
Outputs <- Numextract(Outputs[4,]) %>% as.data.frame(); names(Outputs) <- "Outlet"
Outputs$Component <- row.names(Outputs)

newDat <- full_join(Inputs,Outputs) %>% dplyr::select(Component,Inlet,Outlet)
newDat$Inlet %<>% as.character() %<>% as.numeric();#newDat
newDat$Outlet %<>% as.character() %<>% as.numeric();#newDat
CO2_Tot <- sum(newDat[1:2,2:3]$Outlet)
H2_Tot <- sum(newDat[3:4,2:3]$Outlet)
CH4_Tot <- sum(newDat[5:6,2:3]$Outlet)
sum <- data.frame(c("H2","CO2","CH4"),c(H2_Tot,CO2_Tot,CH4_Tot)); names(sum) <- c("Component","Value"); sum

#### Theoretical Yield of Methanol ####

# Carbon Dioxide/Hydrogen to Methanol
  # CO2 + 3H2 → CH3OH + H2O --> H2 limiting
  # Link: http://www.chemikinternational.com/wp-content/uploads/2014/02/8.pdf
out_Meth_H2 <- H2_Tot/3
out_CO2 <- CO2_Tot - H2_Tot/3
# Methane to Methanol
  # CH4 + (1/2)O2 → CH3OH , DelG = -111 kJ/mol --> CH4 limiting
  # Link: https://www.sciencedirect.com/science/article/pii/S1877705813000581
out_Meth_CH4 <- CH4_Tot
Methanol_Tot <- out_Meth_CH4 + out_Meth_H2
#Yield of Methanol per Methane
Methanol_Tot/out_Meth_CH4
#Yield of Methanol per Methane
Methanol_Tot/out_Meth_H2
#Total Methanol
Methanol_Tot
#CO2 Left
out_CO2

data.frame(Components=c("H2","CO2","CH4","Methanol"),Case1 = c(H2_Tot,CO2_Tot,CH4_Tot,0),Case2 =c(0,out_CO2,0,Methanol_Tot))

###Summary Tables###

head(SumDat_Chem); summary(SumDat_Bact)
  #Main Products
  SumDat_Chem$Flow <- rep(c("Inlet","Outlet"),dim(SumDat_Chem)[[1]]/2)
  summary(SumDat_Chem[SumDat_Chem$Flow=="Outlet",1:10] %>% round(5)) %>% kable('html',align = "c",caption="Main Outputs") %>% 
    kable_styling(bootstrap_options = "striped", full_width = F,font_size=12)
  #Bacteria
  SumDat_Bact$Flow <- rep(c("Inlet","Outlet"),dim(SumDat_Bact)[[1]]/2)
  tmp <- SumDat_Bact %>% gather(Bacteria,Value,-Flow) %>% dplyr::group_by(Bacteria) %>%
    dplyr::summarise(Initial=min(Value),Final=max(Value)) %>% as.data.frame() #%>%
  tmp2 <- cast(tmp,~Bacteria,value="Initial");tmp2$value <- "Input"; tmp2 %<>% dplyr::rename(Flow=value)
  tmp3 <- cast(tmp,~Bacteria,value="Final");tmp3$value <- "Output"; tmp3 %<>% dplyr::rename(Flow=value)
  rbind(tmp2,tmp3) %>% kable('html',align = "c",caption="Bacteria Biomass") %>% 
    kable_styling(bootstrap_options = "striped", full_width = F,font_size=12)
   # summary(SumDat_Chem[SumDat_Chem$Flow=="Outlet",1:10] %>% round(5)) %>% kable('latex',align = "c") %>% 
  #   kable_styling(bootstrap_options = "striped", full_width = F,font_size=12)
  # summary(SumDat_Chem[SumDat_Chem$Flow=="Outlet",1:10] %>% round(5)) %>% kable('markdown',align = "c") 
  

  library(summarytools)
  data <- SumDat_Chem[SumDat_Chem$Flow=="Outlet",1:10] %>% round(5)
  summarytools::descr(data)
  summarytools::descr(data, transpose = TRUE)  
  tmp <- dfSummary(data)
  
  