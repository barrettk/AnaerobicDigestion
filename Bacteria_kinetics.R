
try({
suppressMessages(library(reshape,quietly = T))
suppressMessages(library(dplyr,quietly = T))
suppressMessages(library(polynom,quietly = T))
suppressMessages(library(ggplot2,quietly = T))
#suppressMessages(library(cowplot,quietly = T))
},silent = T)
# 

Meth_dat <- read.csv("data/Methanogen_MSB_kinetics.csv",sep=",",row.names = NULL) %>% as.data.frame()
Aceto_dat <- read.csv("data/Acetogen_bakii_kinetics.csv",sep=",",row.names = NULL) %>% as.data.frame()

#sample size
sampl <- 500
font_size <- 4

# Acetogens ---------------------------------------------------------------

#Vmax
vmax_T4 <- round(rnorm(sampl,mean=Aceto_dat$vmax_Mean[Aceto_dat$Temperature==4],sd=Aceto_dat$vmax_SE[Aceto_dat$Temperature==4]/3),2)
vmax_T10 <- round(rnorm(sampl,mean=Aceto_dat$vmax_Mean[Aceto_dat$Temperature==10],sd=Aceto_dat$vmax_SE[Aceto_dat$Temperature==10]/3),2)
vmax_T15 <- round(rnorm(sampl,mean=Aceto_dat$vmax_Mean[Aceto_dat$Temperature==15],sd=Aceto_dat$vmax_SE[Aceto_dat$Temperature==15]/3),2)
vmax_T20 <- round(rnorm(sampl,mean=Aceto_dat$vmax_Mean[Aceto_dat$Temperature==20],sd=Aceto_dat$vmax_SE[Aceto_dat$Temperature==20]/3),2)
vmax_T25 <- round(rnorm(sampl,mean=Aceto_dat$vmax_Mean[Aceto_dat$Temperature==25],sd=Aceto_dat$vmax_SE[Aceto_dat$Temperature==25]/3),2)
vmax_T30 <- round(rnorm(sampl,mean=Aceto_dat$vmax_Mean[Aceto_dat$Temperature==30],sd=Aceto_dat$vmax_SE[Aceto_dat$Temperature==30]/3),2)
df_vmax <- data.frame(vmax_T4,vmax_T10,vmax_T15,vmax_T20,vmax_T25,vmax_T30)
DF_Aceto_vmax <- suppressMessages(melt(df_vmax)); DF_Aceto_vmax$variable <- as.character(DF_Aceto_vmax$variable)
DF_Aceto_vmax <- DF_Aceto_vmax %>% dplyr::rename(Temperature=variable, Vmax=value)
for(i in 1:dim(DF_Aceto_vmax)[1]){
  if(DF_Aceto_vmax$Temperature[i]=="vmax_T4") DF_Aceto_vmax$Temperature[i] <- as.numeric(4)
  if(DF_Aceto_vmax$Temperature[i]=="vmax_T10") DF_Aceto_vmax$Temperature[i] <- as.numeric(10)
  if(DF_Aceto_vmax$Temperature[i]=="vmax_T15") DF_Aceto_vmax$Temperature[i] <- as.numeric(15)
  if(DF_Aceto_vmax$Temperature[i]=="vmax_T20") DF_Aceto_vmax$Temperature[i] <- as.numeric(20)
  if(DF_Aceto_vmax$Temperature[i]=="vmax_T25") DF_Aceto_vmax$Temperature[i] <- as.numeric(25)
  if(DF_Aceto_vmax$Temperature[i]=="vmax_T30") DF_Aceto_vmax$Temperature[i] <- as.numeric(30)}
DF_Aceto_vmax$Temperature <- as.numeric(DF_Aceto_vmax$Temperature)
#Plot and Fit
fit_Aceto_vmax <<- lm(Vmax ~ poly(Temperature, 3,raw = T), DF_Aceto_vmax)
assign("fit_Aceto_vmax",fit_Aceto_vmax,envir = globalenv())
fit_Aceto_vmax_txt <- as.character(signif(as.polynomial(coef(fit_Aceto_vmax)), 3))
fit_Aceto_vmax_lbl <- paste(gsub("x", "~italic(x)", fit_Aceto_vmax_txt, fixed = TRUE),paste("italic(R)^2",
                       format(summary(fit_Aceto_vmax)$r.squared,digits = 2), sep = "~`=`~"),sep = "~~~~")
predict(fit_Aceto_vmax,data.frame(Temperature=30))
try({
Aceto_vmax_p <- ggplot(DF_Aceto_vmax,aes(Temperature,Vmax)) +
  geom_boxplot(aes(group = reshape::round_any(Temperature,5,floor)), outlier.alpha = 0.1) +
  scale_x_continuous(breaks = c(4,10,15,20,25,30)) + 
  stat_smooth(method="lm", se=TRUE, fill=NA,formula = y ~ poly(x, 3,raw = T), colour="red") +
  ggtitle("Vmax for Acetogens") + labs(caption = parse(text=fit_Aceto_vmax_lbl))
},silent=T)

#Km
km_T4 <- round(rnorm(sampl,mean=Aceto_dat$km_Mean[Aceto_dat$Temperature==4],sd=Aceto_dat$km_SE[Aceto_dat$Temperature==4]/3),2)
km_T10 <- round(rnorm(sampl,mean=Aceto_dat$km_Mean[Aceto_dat$Temperature==10],sd=Aceto_dat$km_SE[Aceto_dat$Temperature==10]/3),2)
km_T15 <- round(rnorm(sampl,mean=Aceto_dat$km_Mean[Aceto_dat$Temperature==15],sd=Aceto_dat$km_SE[Aceto_dat$Temperature==15]/3),2)
km_T20 <- round(rnorm(sampl,mean=Aceto_dat$km_Mean[Aceto_dat$Temperature==20],sd=Aceto_dat$km_SE[Aceto_dat$Temperature==20]/3),2)
km_T25 <- round(rnorm(sampl,mean=Aceto_dat$km_Mean[Aceto_dat$Temperature==25],sd=Aceto_dat$km_SE[Aceto_dat$Temperature==25]/3),2)
km_T30 <- round(rnorm(sampl,mean=Aceto_dat$km_Mean[Aceto_dat$Temperature==30],sd=Aceto_dat$km_SE[Aceto_dat$Temperature==30]/3),2)
df_km <- data.frame(km_T4,km_T10,km_T15,km_T20,km_T25,km_T30)
DF_Aceto_km <- suppressMessages(melt(df_km)); DF_Aceto_km$variable <- as.character(DF_Aceto_km$variable)
DF_Aceto_km <- DF_Aceto_km %>% dplyr::rename(Temperature=variable, Km=value)
for(i in 1:dim(DF_Aceto_km)[1]){
  if(DF_Aceto_km$Temperature[i]=="km_T4") DF_Aceto_km$Temperature[i] <- as.numeric(4)
  if(DF_Aceto_km$Temperature[i]=="km_T10") DF_Aceto_km$Temperature[i] <- as.numeric(10)
  if(DF_Aceto_km$Temperature[i]=="km_T15") DF_Aceto_km$Temperature[i] <- as.numeric(15)
  if(DF_Aceto_km$Temperature[i]=="km_T20") DF_Aceto_km$Temperature[i] <- as.numeric(20)
  if(DF_Aceto_km$Temperature[i]=="km_T25") DF_Aceto_km$Temperature[i] <- as.numeric(25)
  if(DF_Aceto_km$Temperature[i]=="km_T30") DF_Aceto_km$Temperature[i] <- as.numeric(30)}
DF_Aceto_km$Temperature <- as.numeric(DF_Aceto_km$Temperature)
#Plot and Fit
fit_Aceto_km <<- lm(Km ~ poly(Temperature, 3,raw = T), DF_Aceto_km)
assign("fit_Aceto_km",fit_Aceto_km,envir = globalenv())
fit_Aceto_km_txt <- as.character(signif(as.polynomial(coef(fit_Aceto_km)), 3))
fit_Aceto_km_lbl <- paste(gsub("x", "~italic(x)", fit_Aceto_km_txt, fixed = TRUE),paste("italic(R)^2",
                       format(summary(fit_Aceto_km)$r.squared,digits = 2), sep = "~`=`~"),sep = "~~~~")
try({
Aceto_km_p <- ggplot(DF_Aceto_km,aes(Temperature,Km)) +
  geom_boxplot(aes(group = reshape::round_any(Temperature,5,floor)), outlier.alpha = 0.1) +
  scale_x_continuous(breaks = c(4,10,15,20,25,30)) + 
  stat_smooth(method="lm", se=TRUE, fill=NA,formula = y ~ poly(x, 3,raw = T), colour="red") +
  ggtitle("Km for Acetogens") +labs(caption = parse(text=fit_Aceto_km_lbl))
},silent=T)

# Methanogens -------------------------------------------------------------

#vmax
vmax_T4 <- round(rnorm(sampl,mean=Meth_dat$vmax_Mean[Meth_dat$Temperature==4],sd=Meth_dat$vmax_SE[Meth_dat$Temperature==4]/3),2)
vmax_T10 <- round(rnorm(sampl,mean=Meth_dat$vmax_Mean[Meth_dat$Temperature==10],sd=Meth_dat$vmax_SE[Meth_dat$Temperature==10]/3),2)
vmax_T15 <- round(rnorm(sampl,mean=Meth_dat$vmax_Mean[Meth_dat$Temperature==15],sd=Meth_dat$vmax_SE[Meth_dat$Temperature==15]/3),2)
vmax_T20 <- round(rnorm(sampl,mean=Meth_dat$vmax_Mean[Meth_dat$Temperature==20],sd=Meth_dat$vmax_SE[Meth_dat$Temperature==20]/3),2)
vmax_T25 <- round(rnorm(sampl,mean=Meth_dat$vmax_Mean[Meth_dat$Temperature==25],sd=Meth_dat$vmax_SE[Meth_dat$Temperature==25]/3),2)
vmax_T30 <- round(rnorm(sampl,mean=Meth_dat$vmax_Mean[Meth_dat$Temperature==30],sd=Meth_dat$vmax_SE[Meth_dat$Temperature==30]/3),2)
df_vmax <- data.frame(vmax_T4,vmax_T10,vmax_T15,vmax_T20,vmax_T25,vmax_T30)
DF_Meth_vmax <- suppressMessages(melt(df_vmax)); DF_Meth_vmax$variable <- as.character(DF_Meth_vmax$variable)
DF_Meth_vmax <- DF_Meth_vmax %>% dplyr::rename(Temperature=variable, Vmax=value)
for(i in 1:dim(DF_Meth_vmax)[1]){
  if(DF_Meth_vmax$Temperature[i]=="vmax_T4") DF_Meth_vmax$Temperature[i] <- as.numeric(4)
  if(DF_Meth_vmax$Temperature[i]=="vmax_T10") DF_Meth_vmax$Temperature[i] <- as.numeric(10)
  if(DF_Meth_vmax$Temperature[i]=="vmax_T15") DF_Meth_vmax$Temperature[i] <- as.numeric(15)
  if(DF_Meth_vmax$Temperature[i]=="vmax_T20") DF_Meth_vmax$Temperature[i] <- as.numeric(20)
  if(DF_Meth_vmax$Temperature[i]=="vmax_T25") DF_Meth_vmax$Temperature[i] <- as.numeric(25)
  if(DF_Meth_vmax$Temperature[i]=="vmax_T30") DF_Meth_vmax$Temperature[i] <- as.numeric(30)}
DF_Meth_vmax$Temperature <- as.numeric(DF_Meth_vmax$Temperature)
#Plot and Fit
fit_Meth_vmax <<- lm(Vmax ~ poly(Temperature, 3,raw = T), DF_Meth_vmax)
assign("fit_Meth_vmax",fit_Meth_vmax,envir = globalenv())
fit_Meth_vmax_txt <- as.character(signif(as.polynomial(coef(fit_Meth_vmax)), 3))
fit_Meth_vmax_lbl <- paste(gsub("x", "~italic(x)", fit_Meth_vmax_txt, fixed = TRUE),paste("italic(R)^2",
                     format(summary(fit_Meth_vmax)$r.squared,digits = 2), sep = "~`=`~"),sep = "~~~~")
try({
Meth_vmax_p <- ggplot(DF_Meth_vmax,aes(Temperature,Vmax)) +
  geom_boxplot(aes(group = reshape::round_any(Temperature,5,floor)), outlier.alpha = 0.1) +
  scale_x_continuous(breaks = c(4,10,15,20,25,30)) + 
  stat_smooth(method="lm", se=TRUE, fill=NA,formula = y ~ poly(x, 3,raw = T), colour="red") +
  ggtitle("Vmax for Methanogens") +labs(caption = parse(text=fit_Meth_vmax_lbl))
},silent=T)

#km
km_T4 <- round(rnorm(sampl,mean=Meth_dat$km_Mean[Meth_dat$Temperature==4],sd=Meth_dat$km_SE[Meth_dat$Temperature==4]/3),2)
km_T10 <- round(rnorm(sampl,mean=Meth_dat$km_Mean[Meth_dat$Temperature==10],sd=Meth_dat$km_SE[Meth_dat$Temperature==10]/3),2)
km_T15 <- round(rnorm(sampl,mean=Meth_dat$km_Mean[Meth_dat$Temperature==15],sd=Meth_dat$km_SE[Meth_dat$Temperature==15]/3),2)
km_T20 <- round(rnorm(sampl,mean=Meth_dat$km_Mean[Meth_dat$Temperature==20],sd=Meth_dat$km_SE[Meth_dat$Temperature==20]/3),2)
km_T25 <- round(rnorm(sampl,mean=Meth_dat$km_Mean[Meth_dat$Temperature==25],sd=Meth_dat$km_SE[Meth_dat$Temperature==25]/3),2)
km_T30 <- round(rnorm(sampl,mean=Meth_dat$km_Mean[Meth_dat$Temperature==30],sd=Meth_dat$km_SE[Meth_dat$Temperature==30]/3),2)
df_km <- data.frame(km_T4,km_T10,km_T15,km_T20,km_T25,km_T30)
DF_Meth_km <- suppressMessages(melt(df_km)); DF_Meth_km$variable <- as.character(DF_Meth_km$variable)
DF_Meth_km <- DF_Meth_km %>% dplyr::rename(Temperature=variable, Km=value)
for(i in 1:dim(DF_Meth_km)[1]){
  if(DF_Meth_km$Temperature[i]=="km_T4") DF_Meth_km$Temperature[i] <- as.numeric(4)
  if(DF_Meth_km$Temperature[i]=="km_T10") DF_Meth_km$Temperature[i] <- as.numeric(10)
  if(DF_Meth_km$Temperature[i]=="km_T15") DF_Meth_km$Temperature[i] <- as.numeric(15)
  if(DF_Meth_km$Temperature[i]=="km_T20") DF_Meth_km$Temperature[i] <- as.numeric(20)
  if(DF_Meth_km$Temperature[i]=="km_T25") DF_Meth_km$Temperature[i] <- as.numeric(25)
  if(DF_Meth_km$Temperature[i]=="km_T30") DF_Meth_km$Temperature[i] <- as.numeric(30)}
DF_Meth_km$Temperature <- as.numeric(DF_Meth_km$Temperature)
#Plot and Fit
fit_Meth_km <<- lm(Km ~ poly(Temperature, 3,raw = T), DF_Meth_km)
assign("fit_Meth_km",fit_Meth_km,envir = globalenv())
fit_Meth_km_txt <- as.character(signif(as.polynomial(coef(fit_Meth_km)), 3))
fit_Meth_km_lbl <- paste(gsub("x", "~italic(x)", fit_Meth_km_txt, fixed = TRUE),paste("italic(R)^2",
                      format(summary(fit_Meth_km)$r.squared,digits = 2), sep = "~`=`~"),sep = "~~~~")
try({
Meth_km_p <- ggplot(DF_Meth_km,aes(Temperature,Km)) +
  geom_boxplot(aes(group = reshape::round_any(Temperature,5,floor)), outlier.alpha = 0.1) +
  scale_x_continuous(breaks = c(4,10,15,20,25,30)) + 
  stat_smooth(method="lm", se=TRUE, fill=NA,formula = y ~ poly(x, 3,raw = T), colour="red") +
  ggtitle("Km for Methanogens") +labs(caption = parse(text=fit_Meth_km_lbl))
},silent=T)
rm("i")
# Relavent Parameters
# Acetogens
# Vmax
coef_Aceto_vmax <<- coef(fit_Aceto_vmax)
Intercept_Aceto_vmax <<- as.numeric(coef_Aceto_vmax[1])
a_Aceto_vmax <<- as.numeric(coef_Aceto_vmax[2])
b_Aceto_vmax <<- as.numeric(coef_Aceto_vmax[3])
c_Aceto_vmax <<- as.numeric(coef_Aceto_vmax[4])
# Km
coef_Aceto_km <<- coef(fit_Aceto_km)
Intercept_Aceto_km <<- as.numeric(coef_Aceto_km[1])
a_Aceto_km <<- as.numeric(coef_Aceto_km[2])
b_Aceto_km <<- as.numeric(coef_Aceto_km[3])
c_Aceto_km <<- as.numeric(coef_Aceto_km[4])
# Methanogens
# Vmax
coef_Meth_vmax <<- coef(fit_Meth_vmax)
Intercept_Meth_vmax <<- as.numeric(coef_Meth_vmax[1])
a_Meth_vmax <<- as.numeric(coef_Meth_vmax[2])
b_Meth_vmax <<- as.numeric(coef_Meth_vmax[3])
c_Meth_vmax <<- as.numeric(coef_Meth_vmax[4])
# Km
coef_Meth_km <<- coef(fit_Meth_km)
Intercept_Meth_km <<- as.numeric(coef_Meth_km[1])
a_Meth_km <<- as.numeric(coef_Meth_km[2])
b_Meth_km <<- as.numeric(coef_Meth_km[3])
c_Meth_km <<- as.numeric(coef_Meth_km[4])

#plot all 4 graphs
#multiplot(Aceto_vmax_p,Aceto_km_p,Meth_vmax_p,Meth_km_p,cols=2)
