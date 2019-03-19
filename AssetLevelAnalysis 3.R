##############################################
## Bivariate Exit Dynamics (Multiple & Timing)
##############################################
## *) load packages & functions ----
rm(list=ls()) # remove workspace objects
computer <- "mac"
library(dplyr)
library(robustbase)
library(VineCopula)
library(MASS)
library(plot3D)
# library(rgl)
library(kdecopula)
library(ggplot2)
library(fclust)
library(latex2exp)
library(quantreg)
library(zoo)
library(survival)
library(lattice)
library(sgt)
# library(VGAM)
# library(Zelig)
# library(zeligverse)
# library(plm)
# library(lme4)
# library(censReg)
# library(xtable)
# require(KernelKnn)
###################################
## a) read & prepare data  -----
if(FALSE){
if(computer == "win") setwd("~/Data/Company Level/Gerd Export")
if(computer == "mac") setwd("~/Desktop/Company Level/Gerd Export")
# G <- read.csv2("20170511_qryExportFundPfc.csv")
# G <- read.csv2("20170522_qryExportFundPfc_v02.csv")
G <- read.csv2("170601_export.csv")
# G$FundInv_Quarter <- as.Date(G$FundInv_Quarter,format="%m/%d/%Y")
# G$Investment_Date <- as.Date(G$Investment_Date,format="%m/%d/%Y")
G$FundInv_Quarter <- as.Date(G$FundInv_Quarter,format="%d.%m.%y")
G$Investment_Date <- as.Date(G$Investment_Date,format="%d.%m.%y")
G$Exit_Date <- as.Date(G$Exit_Date, format="%d.%m.%y")
G$Pfc_Date1stInvestment <- as.Date(G$Pfc_Date1stInvestment, format="%d.%m.%")

G[,c(21:35)] <- apply(G[,c(21:35)],2,function(x)as.numeric(as.character(x)))

length(unique(as.factor(G$FundInv_FundId)))
length(unique(as.factor(G$Inv_EmittentId))) 
G$Fund_Emi_ID <- paste(G$FundInv_FundId,G$Inv_EmittentId,sep="_")
length(unique(as.factor(G$Fund_Emi_ID)))

G[G$Fund_Emi_ID =="265_4997" & G$FundInv_Quarter == "2016-12-31",c("Disposal","FMV_PeriodEnd")] <- c(11567173,0)

'
g.spz <- data.frame(G %>% group_by(Fund_Emi_ID) %>% summarise(
 FundInv_FundId = FundInv_FundId[1],
 Inv_EmittentId = Inv_EmittentId[1],
 uhp = sum(duplicated(Holding_Period)) - sum(is.na(Holding_Period))
))

g.spz <- g.spz[g.spz$uhp > 0, ]
write.csv2(g.spz, file="Duplicated Entries.csv")
View(G[G$Fund_Emi_ID == "28_5596",])
'

G <- data.frame(G %>% group_by(Fund_Emi_ID) %>% 
                   filter(!duplicated(Holding_Period)))

G.sum <- data.frame(G %>% group_by(Fund_Emi_ID) %>% summarise(
  MinDate = min(FundInv_Quarter,na.rm=T),
  MaxDate = max(FundInv_Quarter,na.rm=T),
  InvDate= Investment_Date[1],
  ExitDate= Exit_Date[1],
  CFdur= sum(FundInv_Proceed * Holding_Period ,na.rm=TRUE)/sum(FundInv_Proceed,na.rm=TRUE)
  - sum(FundInv_Addition * Holding_Period,na.rm=T)/sum(FundInv_Addition,na.rm=T),
  Industry = Industry[1],
  Industry2 = Industry_Level1[1],
  Sector = GICS_Sector[1],
  Vintage = Fund_VintageYear[1],
  Region = Fund_Region[1],
  Currency = Currency_Pfc[1],
  Type = Fund_InvestTypes[1],
  MOIC= (last(MOIC)) ))

table(G.sum$Industry,G.sum$Region)

G.sum$HoPi <- as.numeric(G.sum$ExitDate - G.sum$InvDate)/365.25
G.sum$GeoR <- (G.sum$MOIC)^(1/G.sum$HoPi)

# Load Public Data
if(computer == "win") setwd("C:/Users/christian.tausch/Dropbox/Project D")
if(computer == "mac") setwd("~/Dropbox/Project D")
trindex <- read.csv("TRindex.csv")
trindex$Date <- as.Date(trindex$Date)
hy_spreads <- read.csv("ML-HYOAS.csv")
hy_spreads$DATE <- as.Date(hy_spreads$DATE)
date_seq <- data.frame(DATE= seq(min((hy_spreads$DATE)), max((hy_spreads$DATE)),1))
hy_spreads <- merge(date_seq,hy_spreads,by="DATE",all.x=TRUE)
hy_spreads$ML_HYOAS <- na.locf(hy_spreads$BAMLH0A0HYM2)/100
trindex <- merge(trindex,hy_spreads[,c("DATE","ML_HYOAS")],by.x="Date",by.y="DATE",all.x= TRUE)
rm(date_seq,hy_spreads)

msci_world <- data.frame(Index= trindex$MSCI.World.Net.Return.Daily, Lag_Index= dplyr::lag(trindex$MSCI.World.Net.Return.Daily))
msci_world <- msci_world[-1,]
msci_world$Return_monthly <- msci_world$Index/msci_world$Lag_Index



G.sum <- merge(G.sum,trindex[,c(1,2)],by.x = "MaxDate",by.y="Date")
G.sum <- merge(G.sum,trindex[,c(1,2)],by.x = "MinDate",by.y="Date",suffixes = c(".Exit",".Entry"))
G.sum$MSCI.Multiple <- G.sum$MSCI.World.Net.Return.Daily.Exit/G.sum$MSCI.World.Net.Return.Daily.Entry
gdp <- read.csv2("Cum GDP April 2017.csv",dec=".")
G.sum$YearInvest <- as.numeric(format(G.sum$InvDate,"%Y"))
G.sum <- merge(G.sum,gdp[,c(1,4)],by.x = "YearInvest",by.y="Year")

sum(is.na(G.sum$HoPi))

by(G.sum,format(G.sum$InvDate,"%Y"),function(x) median(x$MOIC,na.rm=T))

quantile(G.sum$Geo,seq(0.95,1,0.002),na.rm=T)
## b-1) Timing: Create G.sucox
G.sucox <- G.sum
G.sucox$Ti <- ifelse(is.na(G.sum$HoPi), as.numeric(G.sum$MaxDate-G.sum$MinDate)/365.25 ,G.sum$HoPi)
G.sucox$Ev <- ifelse(is.na(G.sum$HoPi), 0,1)
G.sucox <- G.sucox[G.sucox$Ti > 0,]
G.sucox$FundAgeAtEntry <-  pmax(0.1, as.numeric(format(G.sucox$MinDate, "%Y"))-as.numeric(G.sucox$Vintage))

df_sector <- data.frame(table(G.sucox$Sector))
names(df_sector) <- c("Sector","Freq_G.sucox")
df_sector$CoxSector <- c("A_ConFinHeaMatUti","A_ConFinHeaMatUti","Energy","A_ConFinHeaMatUti","A_ConFinHeaMatUti",
                         "Industrials","Information Technology","A_ConFinHeaMatUti","A_ConFinHeaMatUti","Real Estate",
                         "Telecommunication Services","A_ConFinHeaMatUti")

G.sucox <- G.sucox[,1:24]
G.sucox <- merge(G.sucox,df_sector[,c("Sector","CoxSector")],by= "Sector",all.x = TRUE)
G.sucox$CoxSector <- as.factor(G.sucox$CoxSector)
levels(G.sucox$CoxSector)
## b-1) Timing: Create G2 (including all assets)
G2 <- merge(G,G.sucox[,c("Fund_Emi_ID","MinDate","MaxDate","MOIC","CoxSector")],by="Fund_Emi_ID",suffixes=c(".dyn",".end"),all.x = TRUE)
G2 <- G2[is.na(G2$Exit_Date) | !(G2$FundInv_Quarter > G2$Exit_Date),]
# G2 <- G2[!is.na(G2$Exit_Date),] # just realized exits

G2$Ev <- ifelse(is.na(G2$Exit_Date),0,
                ifelse(G2$FundInv_Quarter < G2$Exit_Date,0,
                       ifelse(G2$FundInv_Quarter ==  G2$Exit_Date,1,0)))
hpt <- 10
G2$Ev2 <- ifelse(G2$Holding_Period > hpt,NA,
                 ifelse(is.na(G2$Exit_Date) & G2$Holding_Period == hpt,1,
                        ifelse(is.na(G2$Exit_Date) & G2$Holding_Period < hpt,0,
                               ifelse(G2$FundInv_Quarter < G2$Exit_Date & G2$Holding_Period < hpt,0,
                                      ifelse(G2$FundInv_Quarter ==  G2$Exit_Date | G2$Holding_Period == hpt,1,NA)))))

G2$FundAge <- pmax(0,as.numeric(format(G2$FundInv_Quarter,"%Y")) - G2$Fund_VintageYear)
G2$FundAgeAI <- pmax(0,as.numeric(format(G2$Investment_Date,"%Y")) - G2$Fund_VintageYear)

G2$RLT <- ifelse(is.na(G2$Exit_Date),G2$MaxDate-G2$FundInv_Quarter,G2$Exit_Date-G2$FundInv_Quarter)
G2$RLT <- pmax(0.1,G2$RLT/365.25)

G2$T1 <- pmax(0.1,(G2$FundInv_Quarter - G2$Investment_Date)/365.25)
G2$T2 <- G2$T1 + 0.25

G2$default <- as.factor(ifelse(G2$MOIC.dyn < 0.5,"Yes","No"))
G2$default2 <- sqrt(ifelse(G2$MOIC.dyn < 1, 1/pmax(G2$MOIC.dyn,0.1),0))
G2$rocket <- as.factor(ifelse(G2$MOIC.dyn > 3, "Yes","No"))
G2$FLE <- as.factor(ifelse(G2$FundAge > 10,"Dead","Alive"))
fle.age <- 12
G2$FLE <- ifelse(G2$FundAge > fle.age, G2$FundAge-fle.age,0)
G2$FundAgeAI.trans <- pmin(7,G2$FundAgeAI)

mofu <- function(x){
  xx <- x - 3 ; k <- 0.5
  y <- (1/(1+exp(-xx*k)))
  y <- y + exp(-x/k)
  y <- y^(2)  
  return(y)
}
mofu(seq(0,4,0.5))
curve(mofu(x),0,15,ylim=c(0,1)) ; abline(v=0,col="red",lty=3)

# Gamma MOIC transformation
G2$Moic.trans <- mofu(G2$MOIC.dyn)
curve(dgamma(x,shape=2,scale=1),0,5,ylim=c(0,5)) ; abline(a=0,b=1,col="red") ; dgamma(0:10,shape=2,scale=1)
G2$Moic.trans2 <- dgamma(G2$MOIC.dyn,shape=2,scale=1)

colnames(G2)

#if(computer == "win") setwd("C:/Users/christian.tausch/Dropbox/Project D")
#if(computer == "mac") setwd("~/Dropbox/Project D")
#gdp <- read.csv2("Cum GDP April 2017.csv",dec=".")
G2$YearInvest <- as.numeric(format(G2$Investment_Date,"%Y"))
G2$YearQuarter <- as.numeric(format(G2$FundInv_Quarter,"%Y"))
G2 <- merge(G2,gdp[,c(1,4)],by.x = "YearQuarter",by.y="Year")
G2 <- merge(G2,gdp[,c(1,4)],by.x = "YearInvest",by.y="Year",suffixes = c("Quarter","Invest"))


G2 <- merge(G2,trindex[,c(1,2,7,27)],by.x = "FundInv_Quarter",by.y="Date")
G2 <- merge(G2,trindex[,c(1,2,7)],by.x = "Investment_Date",by.y="Date",suffixes = c(".Quarter",".Invest"))
G2 <- merge(G2,trindex[,c(1,2,7)],by.x = "Exit_Date",by.y="Date",all.x=TRUE)
names(G2)[names(G2) == "MSCI.World.Net.Return.Daily"] <- "MSCI.World.Net.Return.Daily.Exit"
names(G2)[names((G2)) == "NASDAQ.100.Total.Return"] <- "NASDAQ.100.Total.Return.Exit"
names(G2)[names(G2) == "ML_HYOAS"] <- "ML_HYOAS.quarter"


G2$MSCI.multi.hist <- G2$MSCI.World.Net.Return.Daily.Quarter/G2$MSCI.World.Net.Return.Daily.Invest
G2$MSCI.multi.future <- G2$MSCI.World.Net.Return.Daily.Exit/G2$MSCI.World.Net.Return.Daily.Quarter
G2$MSCI.geo.future <- G2$MSCI.multi.future^(1/(as.numeric(G2$Exit_Date - G2$FundInv_Quarter)/365.25))-1
## b-1) Timing: Cox Regression (coxph,basehaz)
# coxph & basehaz from survival package

#sun(G2$MSCI.geo.future)
res.NP5 <-  survival::coxph(Surv(time= T1, time2= T2, event= Ev) ~ 
                  strata(Fund_InvestTypes) + FundAgeAI + Moic.trans2 + ML_HYOAS.quarter + GLOQuarter, # + CoxSector, # GICS_Sector
                  data = G2)
#res.wb5 <-  survreg(Surv(time= T1, time2=T2, event = Ev) ~
#                      strata(Fund_InvestTypes) + FundAgeAI + MOIC.dyn + Moic.trans2 + GLOQuarter, 
#                    data = G2) # does not work with interval-censored data
summary(res.NP5)
cox.zph(res.NP5)
survfit(res.NP5)

bh.NP5 <- survival::basehaz(res.NP5,centered=FALSE) # Cumulative basehazard
bh.NP5$S.base <- exp(-1*bh.NP5$hazard)
# bh.NP5$h <- c(0,diff(bh.NP5$hazard))

# Estimate Weibull hazard with non-parametric simulation
np_sim <- function(n=100,strata="VC"){
  bh.NP5 <- bh.NP5[bh.NP5$strata==strata,]
  x <- list()
  for(i in 1:n){    # non-parametric life-time simulation
    pos <- tail(which(runif(1) < bh.NP5$S.base),1)
    if(length(pos) < 1) pos <- 1
    x[i] <- bh.NP5$time[pos]
  }
  return(unlist(x))
}
wb.bo <- fitdistr(np_sim(20000,"BO"),densfun="weibull")
shape= wb.bo$estimate["shape"] ; scale=wb.bo$estimate["scale"] ; bh.NP5$H.wb.bo <- (bh.NP5$time/scale)^shape
wb.vc <- fitdistr(np_sim(20000,"VC"),densfun="weibull")
shape= wb.vc$estimate["shape"] ; scale=wb.vc$estimate["scale"] ; bh.NP5$H.wb.vc <- (bh.NP5$time/scale)^shape

exp.bo <- fitdistr(np_sim(20000,"BO"),densfun="exponential")
exp.vc <- fitdistr(np_sim(20000,"VC"),densfun="exponential")

plot(bh.NP5$time,bh.NP5$S.base)
points(bh.NP5$time,exp(-bh.NP5$time*exp.bo$estimate),col="blue")
points(bh.NP5$time,exp(-bh.NP5$time*exp.vc$estimate),col="red")

table(G2$CoxSector)

ccp_df <- function(n=1,age=NA,cox_model=res.NP5){
  # create input dataframe
  df <-data.frame(Fund_InvestTypes= sample(c("BO","VC"),n,replace=TRUE),
                  AssetAge= sample(seq(1,6),n,replace=TRUE),
                  MOIC.dyn=rlnorm(n,0,0.5),
                  ML_HYOAS.quarter= rep(sample(sample(trindex$ML_HYOAS[!(is.na(trindex$ML_HYOAS))],1),1),n),
                  MSCI.Multiple.Exit=rlnorm(n,0,0.25),
                  Sector= sample(df_sector$Sector,n,replace=TRUE),
                  Fund_Region= sample(c("EU","US"),n,replace = TRUE)
  )
  
  df$RVPI <- pmax(0,df$MOIC.dyn - rlnorm(n,-2,0.6))
  
  df <- merge(df,df_sector[,colnames(df_sector) %in% c("Sector","CoxSector","Sector_BO_m","Sector_VC_m")], by="Sector",all.x=TRUE)
                                                      
  if(!is.na(age)) df$AssetAge <- rep(age,n)
  df$GLOQuarter <- sample(gdp$GLO,1)
  df$FundAgeAI <-  max(df$AssetAge) - df$AssetAge
  df$Holding_Period <- df$AssetAge
  df$Moic.trans2 <- dgamma(df$MOIC.dyn,shape=2,scale=1)
  df$CoxFactor <- exp(predict(cox_model,df))
  df$NAV <- rlnorm(n,meanlog = 10,sd=2)
  
  invisible(df)
}
mean.np5 <- function(t=0.5,strata="BO",CoxF=1){
  if(t<0.5) t=0.5
  df1 <- bh.NP5[bh.NP5$strata==strata,]
  df1$S <- exp(-CoxF*df1$hazard)
  a= diff(df1[df1$time > t,"time"])
  b= df1[df1$time > t & df1$t < df1[nrow(df1),"time"],"S"]
  c= df1[last(which(df1$time <= t)),"S"]
  x= sum(a*b)/c
  return(x)
}
quant.np5 <- function(t_start,t_end,strata="BO",CoxF=1){
  df1 <- bh.NP5[bh.NP5$strata==strata,]
  df1$S <- exp(-CoxF*df1$hazard)
  
  if(t_start > df1$time[1]){
    P_start <- min(df1[df1$time < t_start, "S"])
    P_end <- min(df1[df1$time < t_end, "S"])
  }else{
    P_start <- df1$S[1]
    P_end <- min(df1[df1$time < t_end, "S"])
  }
  
  y <- (1-P_end/P_start)
  
  return(y)
}
cox_plot <- function(df,time=0){
  n <- nrow(df)
  
  print(df)
  
  # plot it
  par(mfrow=c(1,1),mar=c(4.2,4.2,4,2))
  plot(x=0,y=0,xlim= c(-max(df$AssetAge),17),ylim=c(0,1), 
       main="Non-parametric Cox Regression (Cox.Reg: coxph, Cum.Haz: basehaz)",
       ylab="Survival Probability",xlab="Time")
  abline(h=seq(0,0.3,0.1),col="grey")
  abline(v=0,col="grey",lty=2)
  abline(v=time,col="red",lty=3,lwd=2)
  abline(v=-max(df$AssetAge),col="green")
  
  for(i in seq(1,nrow(df))){
    CoxFactor <- df$CoxFactor[i]
    AA <- df$AssetAge[i]
    
    if(df$Fund_InvestTypes[i] == "BO"){
      x.bo <- bh.NP5$time[bh.NP5$strata=="BO"]-AA
      y.bo.np <- exp(-CoxFactor*bh.NP5$hazard[bh.NP5$strata=="BO"])
      y.bo.wb <- exp(-CoxFactor*bh.NP5$H.wb.bo[bh.NP5$strata=="BO"])
      lines(x.bo, y.bo.np, type="l",lwd=2,col="blue")
      lines(x.bo, y.bo.wb, type="l",lty=3,lwd=2,col="blue")
      abline(v= time + mean.np5(t= time+AA,strata = "BO",CoxF=CoxFactor),lwd=1,col="blue")
      history <- which(x.bo < 0)
      lines(x.bo[history], y.bo.np[history], type="l",lwd=2,col="darkgrey")
      lines(x.bo[history], y.bo.wb[history], type="l",lty=3,lwd=2,col="darkgrey")
    }else{
      x.vc <- bh.NP5$time[bh.NP5$strata=="VC"]-AA
      y.vc.np <- exp(-CoxFactor*bh.NP5$hazard[bh.NP5$strata=="VC"])
      y.vc.wb <- exp(-CoxFactor*bh.NP5$H.wb.vc[bh.NP5$strata=="VC"])
      lines(x.vc, y.vc.np, type="l",lwd=2) # Nelson-Aalen Esimator for H
      lines(x.vc, y.vc.wb, type="l",lty=3,lwd=2)
      abline(v= time + mean.np5(t=time+AA,strata = "VC",CoxF=CoxFactor),lwd=1)
      history <- which(x.vc < 0)
      lines(x.vc[history], y.vc.np[history], type="l",lwd=2,col="darkgrey")
      lines(x.vc[history], y.vc.wb[history], type="l",lty=3,lwd=2,col="darkgrey")
    }
    
    
  }
  
  legend("topright",bty="n",legend = c("BO (non-para.)","VC (non-para.)","BO (weibull)","VC (weibull)"),
         col=c("blue","black"),lty=c(1,1,3,3),lwd=2,cex=1.2)
  legend("right",bty="n",legend=c("Fund Age","MOIC (current)","GDP (outlook)","Sector"),cex=1.1)
  
  
}
sim_cox <- function(t=0.5,strata="BO",CoxF=1,real_quantile=NA){
if(t < 0.5) t <- 0.5
if(t > 13) t <- 13

  df1 <- bh.NP5[bh.NP5$strata==strata,]
  df1$S <- exp(-CoxF*df1$hazard)
  c= df1[last(which(df1$time <= t)),"S"]
  
  if(is.na(real_quantile)){
    d=  runif(1) * c
  }else{
    d= real_quantile * c # realized quantile correction
  }
  
  if(d < min(df1$S)){
    rlt <- max(df1$time) - t
  }else{
    wd <- which(df1$S <= d)
    e= df1[first(wd),"time"]
    rlt <- max(0.250002,e - t)
  }
  
  return(rlt)
}
sim_cox(t=0,real_quantile = 0)


df1 <- ccp_df(5)
cox_plot(df1,0)
rm(df1)
## b-2) Multiple: Create G.p2c (just realized exits)
G.p2c <- data.frame(G2 %>% group_by(Fund_Emi_ID) %>% summarise(
  Proceeds_Exit= Total_Proceeds[FundInv_Quarter == Exit_Date[1]][1],
  Proceeds_Last= last(Total_Proceeds),
  FMV_Exit=FMV_PeriodEnd[FundInv_Quarter == Exit_Date[1]][1],
  FMV_Last=last(FMV_PeriodEnd),
  Cost_Exit= Total_Cost[FundInv_Quarter == Exit_Date[1]][1],
  Cost_Last= last(Total_Cost)
))
# View(G[G$Fund_Emi_ID == "103_2223",])

G.p2c <- merge(G2,G.p2c,by="Fund_Emi_ID",all.x=TRUE)
# G.p2c_BACKUP <- G.p2c
G.p2c <- G.p2c_BACKUP

G.p2c_realized <- data.frame(G.p2c %>% group_by(Fund_Emi_ID) %>% filter(FundInv_Quarter <= Exit_Date
                                                                        & Exit.Status == "realized"))
G.p2c_unrealized <- data.frame(G.p2c %>% group_by(Fund_Emi_ID) %>% filter(Exit.Status == "unrealized"))
G.p2c <- rbind(G.p2c_realized,G.p2c_unrealized)

sum(G.p2c$FundInv_Quarter > G.p2c$Exit_Date)
'
View(G.p2c[G.p2c$Fund_Emi_ID =="335_12824",])
'
G.p2c$Time2Exit <- as.numeric(G.p2c$Exit_Date - G.p2c$FundInv_Quarter)/365.25
G.p2c$Distance2Exit <- as.numeric(G.p2c$FundInv_Quarter - G.p2c$Investment_Date)/as.numeric(G.p2c$Exit_Date - G.p2c$Investment_Date)
summary(G.p2c$Distance2Exit)
G.p2c$Distance2Exit <- pmax(0,G.p2c$Distance2Exit)
G.p2c$DPI <- G.p2c$Total_Proceeds/G.p2c$Total_Cost
G.p2c$RVPI <- G.p2c$FMV_PeriodEnd/G.p2c$Total_Cost

G.p2c$FMV_C2come <- G.p2c$FMV_PeriodEnd + G.p2c$Cost_Exit - G.p2c$Total_Cost

G.p2c$P2Come1 <- G.p2c$Proceeds_Exit - G.p2c$Total_Proceeds + ifelse(G.p2c$FMV_Exit > 0, G.p2c$FMV_Exit,0)
G.p2c$P2Come2 <- G.p2c$Proceeds_Last - G.p2c$Total_Proceeds + ifelse(G.p2c$FMV_Last > 0, G.p2c$FMV_Last,0)
G.p2c$P2C.multi1 <- G.p2c$P2Come1 / G.p2c$FMV_C2come
G.p2c$P2C.multi2 <- G.p2c$P2Come2 / G.p2c$FMV_C2come
# Denominator Threshold to reduce variability 
floor.fmv <- 250000 # FMV + Cost 2 Come
G.p2c$P2C.multi1[!is.finite(G.p2c$P2C.multi1) | G.p2c$FMV_C2come < floor.fmv] <- NA
G.p2c$P2C.multi2[!is.finite(G.p2c$P2C.multi2) | G.p2c$FMV_C2come < floor.fmv] <- NA
quantile(G.p2c$P2C.multi1,seq(0,1,0.05),na.rm=TRUE)
quantile(G.p2c$P2C.multi1,c(seq(0,0.02,0.0025),seq(0.98,1,0.0025)),na.rm=TRUE)
quantile(G.p2c$P2C.multi2,seq(0,1,0.05),na.rm=TRUE)
quantile(G.p2c$P2C.multi2,c(seq(0,0.02,0.0025),seq(0.98,1,0.0025)),na.rm=TRUE)
# Negative FVM correction
table(G.p2c$P2C.multi1 < 0, G.p2c$FMV_Exit < 0) # negative FMVs
table(G.p2c$P2C.multi2 < 0, G.p2c$FMV_Last < 0)
'
View(G.p2c[G.p2c$Fund_Emi_ID =="393_13360",])
View(G.p2c[G.p2c$Fund_Emi_ID =="305_4619",])
'
G.p2c$FMV_Exit[G.p2c$FMV_Exit < 0] <- 0
G.p2c$FMV_Last[G.p2c$FMV_Last < 0] <- 0
# negative P2Come Multiple
levels(as.factor(G.p2c$Fund_Emi_ID[G.p2c$P2Come1 < 0])) # IDs with negative P2Come1
'
View(G.p2c[G.p2c$Fund_Emi_ID =="108_2126",]) # intermediate negative Proceeds
View(G.p2c[G.p2c$Fund_Emi_ID =="117_2092",]) # intermediate re-booking Proceeds
View(G.p2c[G.p2c$Fund_Emi_ID =="119_3906",]) # last Proceeds are negative
View(G.p2c[G.p2c$Fund_Emi_ID =="141_2052",]) # last Proceeds are negative = disposal
View(G.p2c[G.p2c$Fund_Emi_ID =="158_3090",]) # intermediate negative Proceeds
View(G.p2c[G.p2c$Fund_Emi_ID =="159_6340",]) # intermediate re-booking Proceeds
View(G.p2c[G.p2c$Fund_Emi_ID =="189_3875",]) # last Proceeds are negative
View(G.p2c[G.p2c$Fund_Emi_ID =="196_4609",]) # last Proceeds are negative
View(G.p2c[G.p2c$Fund_Emi_ID =="203_4531",]) # some minor intermediate negative Proceeds
View(G.p2c[G.p2c$Fund_Emi_ID =="203_4538",]) # some minor intermediate negative Proceeds
# Fund 203: some minor intermediate negative Proceeds ...
View(G.p2c[G.p2c$Fund_Emi_ID =="209_2571",]) # intermediate negative Proceeds
View(G.p2c[G.p2c$Fund_Emi_ID =="214_6392",]) # last Proceeds are negative
View(G.p2c[G.p2c$Fund_Emi_ID =="216_2585",]) # last Proceeds are negative
View(G.p2c[G.p2c$Fund_Emi_ID =="216_6703",]) # last Proceeds are negative
View(G.p2c[G.p2c$Fund_Emi_ID =="22_2208",]) # last Proceeds are negative
View(G.p2c[G.p2c$Fund_Emi_ID =="22_2274",]) # last Proceeds are negative
View(G.p2c[G.p2c$Fund_Emi_ID =="24_2230",]) # last Proceeds are negative
View(G.p2c[G.p2c$Fund_Emi_ID =="24_2398",]) # intermediate negative Proceeds
View(G.p2c[G.p2c$Fund_Emi_ID =="282_6382",]) # intermediate negative Proceeds

'
table(G.p2c$P2C.multi1 < 0, G.p2c$P2C.multi2 < 0)
table(G.p2c$P2C.multi1 < 0, G.p2c$P2Come1 < 0)
table(G.p2c$P2C.multi2 < 0, G.p2c$P2Come2 < 0)


G.p2c$P2C.multiA <- ifelse(G.p2c$P2C.multi1 < 0 & G.p2c$P2C.multi2 < 0, NA, G.p2c$P2C.multi1)
G.p2c$P2C.multiA <- ifelse(G.p2c$P2C.multiA < 0, G.p2c$P2C.multi2,G.p2c$P2C.multiA)
G.p2c$P2C.multiB <- ifelse(G.p2c$P2C.multi1 < 0 & G.p2c$P2C.multi2 < 0, 0, G.p2c$P2C.multi2)
G.p2c$P2C.multiB <- ifelse(G.p2c$P2C.multiB < 0, G.p2c$P2C.multi1,G.p2c$P2C.multiB)
cap <- 10000000
G.p2c$P2C.multi1[G.p2c$P2C.multi1 < 0] <- 0
G.p2c$P2C.multi1[G.p2c$P2C.multi1 > cap] <- cap
G.p2c$P2C.multi2[G.p2c$P2C.multi2 < 0] <- 0
G.p2c$P2C.multi2[G.p2c$P2C.multi2 > cap] <- cap
quantile(G.p2c$P2C.multi1,seq(0,1,0.05),na.rm=TRUE)
rm(cap,floor.fmv)

# merge GDP
# if(computer == "win") setwd("C:/Users/christian.tausch/Dropbox/Project D")
# if(computer == "mac") setwd("~/Dropbox/Project D")
# gdp <- read.csv2("Cum GDP April 2017.csv",dec=".")

'
G.p2c$YearInvest <- as.numeric(format(G.p2c$Investment_Date,"%Y"))
G.p2c$YearQuarter <- as.numeric(format(G.p2c$FundInv_Quarter,"%Y"))
# merge GDP
G.p2c <- merge(G.p2c,gdp[,c(1,4)],by.x = "YearQuarter",by.y="Year")
G.p2c <- merge(G.p2c,gdp[,c(1,4)],by.x = "YearInvest",by.y="Year",suffixes = c("Quarter","Invest"))
# merge Index
G.p2c <- merge(G.p2c,trindex[,c(1,2,27)],by.x = "FundInv_Quarter",by.y="Date")
G.p2c <- merge(G.p2c,trindex[,c(1,2,7)],by.x = "Investment_Date",by.y="Date",suffixes = c(".Quarter",".Invest"))
G.p2c <- merge(G.p2c,trindex[,c(1,2,7)],by.x = "Exit_Date",by.y="Date")
names(G.p2c)[names(G.p2c) == "MSCI.World.Net.Return.Daily"] <- "MSCI.World.Net.Return.Daily.Exit"
names(G.p2c)[names(G.p2c) == "NASDAQ.100.Total.Return"] <- "NASDAQ.100.Total.Return.Exit"
names(G.p2c)[names(G.p2c) == "ML_HYOAS"] <- "ML_HYOAS.quarter"
'

G.p2c <- G.p2c[order(G.p2c$Fund_Emi_ID,G.p2c$FundInv_Quarter),]
G.p2c$MSCI.Multiple <- G.p2c$MSCI.World.Net.Return.Daily.Quarter/G.p2c$MSCI.World.Net.Return.Daily.Invest
G.p2c$MSCI.Multiple.Exit <- G.p2c$MSCI.World.Net.Return.Daily.Exit/G.p2c$MSCI.World.Net.Return.Daily.Quarter


G.p2c$Fund_Age <- pmax(0,G.p2c$YearQuarter -G.p2c$Fund_VintageYear)
G.p2c$Fund_AgeAI <- pmax(0,G.p2c$YearInvest -G.p2c$Fund_VintageYear)

zombie_hp <- 10
G.p2c$ZombieStage <- pmax(zombie_hp,G.p2c$Holding_Period + G.p2c$Time2Exit) - zombie_hp
G.p2c$ZombieStage2 <- G.p2c$ZombieStage > 0
G.p2c$Time2Exit_ZS <- G.p2c$Time2Exit + G.p2c$ZombieStage
G.p2c$FundAgeAE <- as.numeric(format(G.p2c$Exit_Date,"%Y")) - G.p2c$Fund_VintageYear
G.p2c$AssetAgeAE <- as.numeric(format(G.p2c$Exit_Date,"%Y")) - as.numeric(format(G.p2c$Investment_Date,"%Y"))

G.p2c$RVPI_1 <- G.p2c$RVPI - 1
G.p2c$MSCI.Multiple.Exit_1 <- G.p2c$MSCI.Multiple.Exit - 1


# Merge: quarterly public equity returns
qper <- trindex[,1:2]
qper <- qper[format(qper$Date,"%m") %in% c("03","06","09","12"),]
qper$MSCI.qrtly.return <- c(diff(qper$MSCI.World.Net.Return.Daily) / qper$MSCI.World.Net.Return.Daily[1:(nrow(qper)-1)],NA)
G.p2c <- merge(G.p2c,qper[,c(1,3)],by.x = "FundInv_Quarter",by.y = "Date",all.x = TRUE)


# Sector separate for BO and VC
# G.p2c <- G.p2c[,1:111]
VC_m <- data.frame(table(G.p2c$GICS_Sector[G.p2c$Fund_InvestTypes =="VC"]))
names(VC_m) <- c("Sector","Freq_VC_m")
VC_m$Sector_VC_m <- c("IT_Con","IT_Con","Health_Rest","Health_Rest","Health_Rest","Industrials","IT_Con","Health_Rest","Health_Rest","Health_Rest","IT_Con","Health_Rest")
G.p2c <- merge(G.p2c,VC_m[,c(1,3)],by.x = "GICS_Sector",by.y="Sector",all.x = TRUE)
df_sector <- merge(df_sector,VC_m, by= "Sector",all.x= TRUE)

BO_m <- data.frame(table(G.p2c$GICS_Sector[G.p2c$Fund_InvestTypes =="BO"]))
names(BO_m) <- c("Sector","Freq_BO_m")
BO_m$Sector_BO_m <- c("Consumer","Consumer","Energy","Financials","Health Care","Industrials","Information Technology","Rest","Rest","Rest","Telecommunication Services","Rest")
G.p2c <- merge(G.p2c,BO_m[,c(1,3)],by.x = "GICS_Sector",by.y="Sector",all.x = TRUE)
df_sector <- merge(df_sector,BO_m, by= "Sector",all.x= TRUE)

rm(VC_m,BO_m)
# save.image("z_b2.RData")
}else{
  #### load workspace
  if(computer == "win") setwd("C:/Users/christian.tausch/Dropbox/Project D/3_R_Code_Pro_D")
  if(computer == "mac") setwd("~/Dropbox/Project D/3_R_Code_Pro_D")
  load("z_b2.RData")
}


rm(list=setdiff(ls(), c("G.p2c","msci_world","trindex","G.sum")))
# rm(list=setdiff(ls(), c("G.p2c","msci_world","ccp_df")))

dev.off()
par_default <- par()

# define useful functions
# Summarize Univariate distributioN
sun <- function(X,cdf=0,plot.it=TRUE,round.to=3,fit.sgt=FALSE,
                Xlim=range(X),Ylim=c(0,1),main=NA){
  options(warn = -1)
  X <- as.numeric(X)
  nobs <- length(X)
  NAs <- sum(is.na(X))
  Infs <- sum(!is.finite(X) & !is.na(X))
  if(NAs+Infs == length(X)){
    stop("No Numeric Input Data")
  }else{ 
    X <- X[!is.na(X) & is.finite(X)]
    Qs <- quantile(X,probs = c(0.05,0.1,0.25,0.5,0.75,0.9,0.95))
    meanX <- mean(X) ; sdX <- sd(X) ; maxX <- max(X) ; minX <- min(X)
    ECDF <- ecdf(X)
    P.cdf <- ECDF(cdf)
    names(P.cdf) <- paste("Pr(x<=",cdf,")",sep="")
    # SGT fit
    if(length(X) > 10000) fit.sgt <- FALSE
    if(fit.sgt){
      library(sgt)
      SGT <- sgt.mle(~X,start = list(mu=meanX,sigma=sdX,lambda=0,p=2,q=2))
    } 
    
    if(plot.it){
      par(mar=c(2.5,4.1,3.5,1.5))
      if(is.na(main)){
        m1 <- "Summarize Univariate distributioN"
      }else{
        m1 <- main
      }
      plot(ECDF,main=m1,cex=0.5,xlim=Xlim,ylim=Ylim)
      hist(X,add=TRUE,freq = F,border = "darkgrey",breaks=max(10,length(X)/20))
      rug(X,col="black")
      curve(dnorm(x,mean = meanX,sd=sdX),add=TRUE,col="blue",lwd=2)
      abline(v=meanX,col="blue",lwd=2)
      if(fit.sgt) curve(dsgt(x,mu=SGT$estimate[1],sigma=SGT$estimate[2],lambda=SGT$estimate[3],p=SGT$estimate[4],q=SGT$estimate[5] ),col="green",lty=2,add=TRUE)
      
      # Quantiles
      X0 <- c("25%"=minX,"50%"=NA,"75%"=maxX)
      for(i in c(0.25,0.5,0.75)){
        qy <- paste(i*100,"%",sep = "")
        segments(x0=X0[qy],y0=i,x1=Qs[qy], y1=i,lty=3,col="red",lwd=2)
        segments(x0=Qs[qy],y0=0,x1=Qs[qy], y1=i,lty=3,col="red",lwd=2)
      }
      # Legend
      position <- ifelse(meanX - minX > 0.5*(maxX - meanX), "topleft","right")
      if(fit.sgt){
        legend(position, inset = 0.02,bty="n",cex=0.7,legend = c("ECDF","mean/normal","Q25/50/75","SGT"),lty = c(1,1,3,2),col=c("black","blue","red","green"),pch=c(20,NA,NA,NA))
      }else{
        legend(position, inset = 0.02,bty="n",cex=0.7,legend = c("ECDF","mean/normal","Q25/50/75"),lty = c(1,1,3),col=c("black","blue","red"),pch=c(20,NA,NA))
      } 
    }
    out <- c(mean=meanX,sd=sdX,P.cdf,min=minX,Qs,max=maxX)
    out <- round(out,round.to)
    print(out)
    print(c(n=nobs,NAs=NAs,Infs=Infs))
    out <- c(out,n=nobs,NAs=NAs,Infs=Infs)
    if(fit.sgt) out.sgt <- c("SGT"=round(SGT$estimate,round.to),"SGT.pq"=round(prod(SGT$estimate[4:5]),round.to))
    if(fit.sgt) print(out.sgt)
    if(fit.sgt) out <- c(out,out.sgt)
    invisible(out)
  }
}
# logit to probability
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
# Color legend
legend.col <- function(col, lev){
  # https://aurelienmadouasse.wordpress.com/2012/01/13/legend-for-a-continuous-color-scale-in-r/
  
  opar <- par
  
  n <- length(col)
  
  bx <- par("usr")
  
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  
  xx <- rep(box.cx, each = 2)
  
  par(xpd = TRUE)
  for(i in 1:n){
    
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
    
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
}
# R2 calculation
r2 <- function(actual,predict){
  R2 <- 1 - (sum((actual-predict )^2)/sum((actual-mean(actual))^2))
  return(R2)
}
## x-0) MOIC Regression (Full Sample Benchmark Regression on MOIC) -----
if(FALSE){
  G.sum <- merge(G.sum,trindex[,c("Date","ML_HYOAS")],by.x = "InvDate",by.y="Date",all.x = TRUE)
  G.sum$ZombieStage <- pmax(10,G.sum$HoPi) - 10
  G.sum$MSCI.Multiple_1 <- G.sum$MSCI.Multiple - 1
  G.sum$Time2Exit_ZS <- G.sum$HoPi + G.sum$ZombieStage
  
  
  
  # Multivariate Regressions (Parametric)
  summary(lm(MOIC ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1,data= G.sum[G.sum$MOIC > 0,],subset = Type == "VC"))
  summary(robustbase::lmrob(MOIC ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1 ,data= G.sum,subset = Type == "VC",method="SMDM"))
  summary(censReg::censReg(MOIC ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1 ,data= G.sum[G.sum$Type == "VC",],left=0))
  
  summary(lm(MOIC ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1 ,data= G.sum,subset = Type == "BO"))
  summary(robustbase::lmrob(MOIC ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1 ,data= G.sum,subset = Type == "BO",method="SMDM"))
  summary(censReg::censReg(MOIC ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1 ,data= G.sum[G.sum$Type == "BO",],left=0),m=1)
  
  
  
  
  G.sum$MOIC_pois <- round(G.sum$MOIC*100,0)
  G.sum$HoPi_pois <- round(G.sum$HoPi*4,0)
  
  
  
  
  '
  summary(crch::trch(MOIC ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1,data= G.sum[G.sum$MOIC > 0.2,],subset = Type == "VC",
  dist="student",link.scale = "identity",left=0.2))
  
  summary(fm <- flexmix::flexmix(MOIC_pois ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1,data= G.sum[complete.cases(G.sum) & G.sum$Type == "VC",],
  k=2,model = list(FLXMRglm(MOIC_pois~ ., family = "poisson"), 
  FLXMRglm(MOIC_pois ~ ML_HYOAS, family = "poisson"))))
  fm <- flexmix::stepFlexmix(MOIC_pois ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1, model = FLXMRglm(family="Gamma"), 
  data = G.sum[complete.cases(G.sum) & G.sum$Type == "VC",], k = 2, nrep = 10)
  summary(fm)
  summary(refit(fm))
  flexmix::parameters(fm, component = 1, model = 1)
  xy <- G.sum[complete.cases(G.sum) & G.sum$Type == "VC",]
  summary(mixtools::poisregmixEM(y=as.numeric(xy$MOIC_pois), as.matrix(xy[,c("HoPi","ML_HYOAS","MSCI.Multiple_1")]),addintercept=TRUE, k=3))
  pgm_reg <- poisson.glm.mix::pois.glm.mix(reference= xy[,c("ML_HYOAS","MSCI.Multiple_1")],response= xy[,c("HoPi_pois","MOIC_pois")],
  mnr=5,L=c(1,1),Kmax = 2,Kmin = 1,max.iter=1000,m1=3, m2=3, t1=3, t2=3, msplit=3, tsplit=3,m=2)
  xy$cluster <- pgm_reg$est.sel.mod.icl$clust
  xy$MOIC_gauss <- LambertW::Gaussianize(xy$MOIC,type="hh")
  sun(xy$MOIC_gauss)
  sun(xy$MOIC)
  '
  
  # https://stats.idre.ucla.edu/r/dae/zip/
  zi_reg <- pscl::zeroinfl(MOIC_pois ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1 | HoPi ,dist = "negbin", data=G.sum[complete.cases(G.sum) & G.sum$Type == "VC",])
  summary(zi_reg)
  coef(zi_reg)
  predict(zi_reg,G.sum[1:10,])

  
  plot(zi_reg$y,zi_reg$fitted.values,xlim=c(0,1000),ylim=c(0,1000))

  logit2prob(-2.347 )
  
  hurdle_reg <- pscl::hurdle(MOIC_pois ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1 | HoPi ,dist = "negbin", data=G.sum[complete.cases(G.sum) & G.sum$Type == "VC",])
  summary(hurdle_reg)
  coef(hurdle_reg)
  predict(hurdle_reg,G.sum[1:10,])
  
  glm.nb_reg <- MASS::glm.nb(MOIC_pois ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1, data=G.sum[complete.cases(G.sum) & G.sum$Type == "VC",])
  summary(glm.nb_reg)
  coef(glm.nb_reg)
  glm.nb_reg$theta ; glm.nb_reg$SE.theta
  sun(glm.nb_reg$residuals)
  predict(glm.nb_reg,G.sum[1:10,])
  par(mfrow=c(2,1))
  set.seed(99)
  sun(rnegbin(length(glm.nb_reg$fitted.values),mu=glm.nb_reg$fitted.values,0.4025)/100,Xlim = c(0,25))
  sun(G.sum$MOIC[complete.cases(G.sum) & G.sum$Type == "VC"],Xlim = c(0,25))
  
  
  
  glm_ln_vc <- glm(MOIC_pois ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1, data=G.sum[complete.cases(G.sum) & G.sum$Type == "VC",],
                   family=gaussian(link="identity"))
  summary(glm_ln_vc)
  
  
  
  '
  x <- model.matrix( ~ .-1,data=G.sum[,c("MOIC","HoPi","ZombieStage","ML_HYOAS","MSCI.Multiple_1")])
  y <- x[,1]
  x <- x[,-1]
  # Ridge: L2 Regularization
  cv.ridge <- glmnet::cv.glmnet(x, y, family="gaussian"", alpha=0, parallel=TRUE, standardize=TRUE,intercept=TRUE)
  # LASSO: L1 Regularization
  cv.lasso <- glmnet::cv.glmnet(x, y, family="gaussian", alpha=1, parallel=TRUE, standardize=TRUE,intercept=TRUE)
  plot(cv.ridge)
  plot(cv.lasso)
  coef(cv.ridge)
  coef(cv.lasso)
  '
  
  # Procedure to test several regression approches (iter-fold-cross-validation)
  knn_cv <- function(iter=10, max_y=6, df= G.sum,
                     ia_depth= 1,n_minObsNode=10, shrnk = 0.001, n_tree_multi= 10, # bgm
                     method_kknn= "mahalanobis",  # knn
                     weights_kknn= "triangular", k_percent= 40,    # knn
                     do_legend= FALSE){
    
    # which independent parameters to use?
    
    df <- df[complete.cases(df),c("MOIC","HoPi","ZombieStage","ML_HYOAS","MSCI.Multiple_1","Type")]
    linear_formula <- formula("MOIC ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1")
    # df <- df[complete.cases(df),c("MOIC","Time2Exit_ZS","ML_HYOAS","MSCI.Multiple_1","Type")]
    # linear_formula <- formula("MOIC ~ Time2Exit_ZS + ML_HYOAS + MSCI.Multiple_1")
    
    
    
    par(mfrow=c(2,1),mar=c(2,2,2,1),oma=c(3,3,1,1))
    out_list <- list()
    do_regressions <- TRUE
    cex_para <- 0.6
    corr_method <- "kendall"
    cor_tobit <- cor_ols <- cor_robust <- cor_bgm <- cor_knn <- 0
    
    
    
    for(type in c("BO","VC")){
      plot(c(0,max_y),c(0,max_y),type="l",lyt=1,xlab=NA, ylab=NA,
           xlim=c(-max_y,max_y),ylim=c(0,max_y),main=type)
      lines(-c(0,max_y),c(0,max_y))
      abline(v=0)
      
      # subset input by type
      df1 <- df[df$Type == type,-ncol(df)]
      all_rows <- seq(1,nrow(df1),1)
      cv_buckets <- split(sample(all_rows), letters[seq(1,iter)])
      cv_buckets <- split(sample(all_rows), seq(1,iter))
      
      measure <- list()
      for(a in seq(1,iter)){
        print(paste(type,a))
        
        # split df in train & test
        # train_rows <- sample(all_rows,round(length(all_rows)*train_fold),replace = FALSE)
        # test_rows <- all_rows[!(all_rows %in% train_rows)]
        test_rows <- unlist(cv_buckets[a])
        train_rows <- all_rows[!(all_rows %in% test_rows)]
        df_train <- df1[train_rows,]
        df_test <- df1[test_rows,]
        x_train <- model.matrix( ~ .-1,data=df_train)
        x_test <- model.matrix( ~ .-1,data=df_test)
        df_test <- data.frame(x_test)
        
        
        # parametric + boosting + k-nearest-neighbors
        if(do_regressions){
          x_para <- - df_test$MOIC
          
          # Tobit
          tobit_reg <- censReg::censReg(linear_formula,data=df_train,left=0)
          predict_tobit <- function(reg=tobit_reg,df_in=df_test){
            df_test2 <- df_in[,-1] # remove dependent variable
            df_test2$Intercept <- 1 # add intercept
            df_test2 <- df_test2[,c(ncol(df_test2),seq(1,ncol(df_test2)-1))] # intercept as 1. column
            prediction <- list()
            for(i in seq(nrow(df_test2))){
              prediction[i] <- sum(reg$estimate[1:5] * df_test2[i,])
            }
            prediction <- as.numeric(unlist(prediction))
            return(prediction)
          }
          y_tobit <- predict_tobit()
          cor_tobit <- cor_tobit + cor(-x_para,y_tobit,method = corr_method)
          print(paste("Corr",corr_method,"tobit:",round(cor_tobit,4)))
          
          # OLS
          ols_reg <- lm(linear_formula,data=df_train)
          y_ols <- predict(ols_reg,df_test)
          cor_ols <- cor_ols + cor(-x_para,y_ols,method = corr_method)
          print(paste("Corr",corr_method,"ols:", round(cor_ols,4)))
          
          # robust
          rob_reg <- robustbase::lmrob(linear_formula,data=df_train)
          y_rob <- predict(rob_reg,df_test)
          cor_robust <- cor_robust + cor(-x_para,y_rob,method = corr_method)
          print(paste("Corr",corr_method,"robust:",round(cor_robust,4)))
          
          # Boosting
          n_trees <-  n_tree_multi * nrow(df_train)
          gbm_reg <- gbm::gbm(linear_formula,
                              distribution="tdist",
                              n.trees= n_trees,
                              n.minobsinnode=n_minObsNode,shrinkage = shrnk,interaction.depth= ia_depth,
                              # weights= dgamma(df_train$MOIC,shape = 2,scale = 1),bag.fraction=1,
                              data=df_train)
          y_gbm <- predict(gbm_reg,df_test,n.trees=n_trees)
          cor_bgm <- cor_bgm + cor(-x_para,y_gbm,method = corr_method)
          print(paste("Corr",corr_method,"bgm:", round(cor_bgm,4)))
          
          # k-nearest means
          k_min <- round(nrow(df_train)/100)
          #knn.fit <- FNN::knn.reg(train= x_train,test= x_test,y= x_train[,1],k=k_min*i)
          y_knn <- KernelKnn::KernelKnn(data=x_train[,-1], 
                                        TEST_data= x_test[,-1],
                                        y= x_train[,1],
                                        k= k_min*k_percent,
                                        h=1, 
                                        method= method_kknn, 
                                        weights_function= weights_kknn,
                                        regression = TRUE)
          cor_knn <- cor_knn + cor(df_test$MOIC,y_knn,method = corr_method)
          print(paste("Corr",corr_method,"knn:", round(cor_knn,4)))
          '
          # Generalized Additive Model
          gam_reg <- gam::gam(MOIC ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1,data=df_train)
          y_gam <- predict(gam_reg,df_test)
          cor_gam <- cor_gam + cor(-x_para,y_gam,method = corr_method)
          print(paste("Corr",corr_method,"gam:", round(cor_gam,4)))
          # random forest
          rF_reg <- randomForest::randomForest(MOIC ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1,
          data=df_train,ntree= n_trees,maxnodes=10)
          y_rF <- predict(rF_reg,df_test)
          points(x_para,y_rF,pch=20,cex=cex_para*0.7,col="orange")
          print(paste("Corr rForest:", cor(-x_para,y_rF,method= corr_method)))
          '
        }
        
        # Regression Results (Predictions)
        df_test$knn_pred <- y_knn
        df_test$gbm_pred <- y_gbm
        df_test$tob_pred <- y_tobit
        df_test$rob_pred <- y_rob
        df_test$ols_pred <- y_ols
        out_list[[paste("Type",type,sep="")]][[paste("Run",a)]] <- df_test
        
        
        if(a > 0){
          points(x_para,y_tobit,pch=20,cex=cex_para,col="gray75")
          points(x_para,y_ols,pch=20,cex=cex_para*0.9,col="gray50")
          points(x_para,y_rob,pch=20,cex=cex_para*0.9*0.9,col="gray25")
          points(-x_para,y_gbm,pch=20,cex=cex_para,col="palevioletred2")
          points(-x_para, y_knn,pch=20,cex=cex_para*0.9,col="peachpuff")
          # points(-x_para,y_gam,pch=20,cex=cex_para*0.9*0.9,col="lavenderblush1")
        }
        
        }
      if(do_legend){
        legend(x=0.1,y=max_y+0.4,bty="n",col=c(NA,"peachpuff","palevioletred2","lavenderblush1"),legend=c("non-linear",paste("knn ",k_percent,"%",sep=""),"boosting"),pch=c(NA,20,20),cex=1.5)
        legend(x=-2.1,y=max_y+0.4,bty="n",legend=c("linear","OLS","Robust","Tobit"),col=c(NA,"gray50","gray25","gray75"),pch=c(NA,20,20,20),cex=1.5)
      }
      }
    
    mtext("Predictions",side=2,outer=TRUE,line=1,cex=1.4)
    mtext("Observations",side=1,outer=TRUE,line=1,cex=1.4)
    
    cum_cor <- c(tobit= cor_tobit,ols= cor_ols,robust= cor_robust,bgm= cor_bgm, knn= cor_knn)
    cum_cor <-  cum_cor/(2*iter)
    # return(cum_cor)
    invisible(out_list)
    }
  # setEPS() ; postscript("Parametric Vs NonParametric Regression Prediction.eps", width = 5, height = 5, family = "Helvetica",pointsize=5)
  set.seed(99)
  system.time( knn_output <- knn_cv(2) )
  # dev.off()
  # Analyze k-nearest-neighbors output
  knn_analyzer <- function(knn_out= knn_output){
    outlist <- list()
    par(mfrow=c(2,1),mar=c(2,2,2,1),oma=c(3,3,1,1))
    for(type in c("BO","VC")){
      # Extract & Prepare Data
      knn_type <- do.call(rbind,knn_out[[paste("Type",type,sep="")]])
      knn_type <- knn_type[complete.cases(knn_type),]
      print(paste(type,nrow(knn_type),"predictions"))
      df_cor <- knn_type[,colnames(knn_type) %in% "MOIC" | grepl("pred",colnames(knn_type))]
      
      # Plot Residuals
      plot(c(-1,5),c(0,0),type="l",lty=3,lwd=2,ylim=c(-5,80),ylab=NA,xlab=NA)
      i <- 1
      for(reg in colnames(df_cor)[-1]){
        points(knn_type[,reg], knn_type$MOIC-knn_type[,reg], ylab=NA,pch=20,col=1+i,cex=0.7)
        i <- i + 1
      }
      mtext(type,side=2,line=2.5)
      legend("topright",bty="n",legend=c(paste(substring(colnames(df_cor)[-1],1,3))),col=seq(2,length(colnames(df_cor)[-1])+1,1),pch=20)
      
      # Goodness of Fit
      for(come in c("pearson","spearman","kendall")){
        correl <- round(cor(df_cor,method=come)[-1,1],4)
        outlist[[paste("Type",type,sep="_")]][[paste("cor",come,sep="_")]] <- correl
      }
      
    }
    mtext("Predictions",side=1,outer=TRUE,line=1) ; mtext("Residuals",side=2,outer=TRUE,line=1)
    return(outlist)
  }
  knn_analyzer()
  
  
  
  sun(knn_bo$MOIC-knn_bo$knn_pred,main="kNN Residuals (BO)")
  sun(knn_vc$MOIC-knn_vc$knn_pred,main="kNN Residuals (BO)")
  
  gam.object <- gam::gam(MOIC ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1 ,data= G.sum,subset = Type == "BO")
  summary(gam.object) ; sun(gam.object$residuals) ; sun(gam.object$fitted.values)
  
  n_trees <- 10*G.sum[complete.cases(G.sum) & G.sum$Type == "BO",]
  system.time(gbm.object <- gbm::gbm(MOIC ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1 , 
                                     distribution="tdist",n.trees=n_trees,
                                     data= G.sum[complete.cases(G.sum) & G.sum$Type == "BO",]))
  summary(gbm.object)
  gbm.object
  sun(predict(gbm.object,G.sum,n.trees=n_trees))
  
  system.time(rF.object <- randomForest::randomForest(MOIC ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1 , 
                                                      n.trees=n_trees,
                                                      data= G.sum[complete.cases(G.sum) & G.sum$Type == "BO",]))
  sun(predict(rF.object,G.sum))
  
  
  
  
  
  '
  df_X <- G.sum[,c("MOIC","HoPi","ZombieStage","ML_HYOAS","MSCI.Multiple_1")]
  df_X <- df_X[complete.cases(df_X),]
  df_X1 <- df_X[1:500,]
  df_X2 <- df_X[501:nrow(df_X),]
  knn_fit <- FNN::knn.reg(train=df_X1, test = df_X2,  y= df_X1$MOIC,k=1)
  plot(df_X2$MOIC,knn_fit$pred) ; abline(a=0,b=1,col="red")
  sun(df_X2$MOIC - knn_fit$pred)
  
  Kknn_fit <- KernelKnn::KernelKnn(data=df_X1, TEST_data = df_X2,  y= df_X1$MOIC,k=50,h=1,method="euclidean",regression = TRUE)
  plot(df_X2$MOIC,Kknn_fit) ; abline(a=0,b=1,col="red")
  '
  
}
## x-1) Multiple Regression
if(FALSE){
  # outlier down-weighting function
  odw1 <- function(x){
    out <- dgamma(x, shape=2,scale=1) + 0.1
    return(out)
  }
  curve(odw1,to=10)
  
  
  
  formula_multiple_bo <- formula("P2C.multi1 ~ Fund_Region + Sector_BO_m + Holding_Period + Time2Exit + ZombieStage + I(RVPI-1) + ML_HYOAS.quarter + I(MSCI.Multiple.Exit-1)")
  formula_multiple_vc <- formula("P2C.multi1 ~ Fund_Region + Sector_VC_m + Holding_Period + Time2Exit + ZombieStage + I(RVPI-1) + ML_HYOAS.quarter + I(MSCI.Multiple.Exit-1)")
  formula_multiple_bo2 <- formula("P2C.multi1 ~ Holding_Period + Time2Exit + ZombieStage2 +  I(RVPI-1) + ML_HYOAS.quarter + I(MSCI.Multiple.Exit-1)")
  formula_multiple_vc2 <- formula("P2C.multi1 ~ Holding_Period + Time2Exit + ZombieStage +  I(RVPI-1)  + ML_HYOAS.quarter + I(MSCI.Multiple.Exit-1)")
  
  re <- list()
  p2c.threshold <- 0
  
  # OLS
  re$p2c.ols.bo <- lm(formula_multiple_bo2, data=bo_ss)
  re$p2c.ols.vc <- lm(formula_multiple_vc2, data=vc_ss)
  re$p2c.ols_trunc.bo <- lm(formula_multiple_bo2, data=bo_ss, subset= P2C.multi1 > p2c.threshold)
  re$p2c.ols_trunc.vc <- lm(formula_multiple_vc2, data=vc_ss, subset= P2C.multi1 > p2c.threshold)
  re$p2c.ols_odw.bo <- lm(formula_multiple_bo2, data=bo_ss, weights = odw1(P2C.multi1))
  re$p2c.ols_odw.vc <- lm(formula_multiple_vc2, data=vc_ss, weights = odw1(P2C.multi1))
  
  
  # ROBUST (outliers)
  re$p2c.rob.bo <- robustbase::lmrob(formula_multiple_bo2, data=bo_ss, method = "SMDM")
  re$p2c.rob.vc <- robustbase::lmrob(formula_multiple_vc2, data=vc_ss, method = "SMDM")
  re$p2c.rob_trunc.bo <- robustbase::lmrob(formula_multiple_bo2, data=bo_ss, subset= P2C.multi1 > p2c.threshold,method = "SMDM")
  re$p2c.rob_trunc.vc <- robustbase::lmrob(formula_multiple_vc2, data=vc_ss, subset= P2C.multi1 > p2c.threshold,method = "SMDM")
  
  # censReg TOBIT (censoring + possibly random effects)
  re$p2c.censReg.bo <- censReg::censReg(formula_multiple_bo2, left=0,data=bo_ss)
  re$p2c.censReg.vc <- censReg::censReg(formula_multiple_vc2, left=0,data=vc_ss)
  re$p2c.censRegX.bo <- censReg::censReg(formula_multiple_bo2, left=0,right=4,data=bo_ss)
  re$p2c.censRegX.vc <- censReg::censReg(formula_multiple_vc2, left=0,right=4,data=vc_ss)
  
  predict(re$p2c.censReg.bo,bo_ss)
  
  
  
  # Ridge & Lasso
  for(type in c("BO","VC")){
    if(type == "BO"){df_ss1 <- bo_ss}else{df_ss1 <- vc_ss}
    
    df_ss2 <- df_ss1[,c("P2C.multi1","Holding_Period","Time2Exit","ZombieStage","RVPI","ML_HYOAS.quarter","MSCI.Multiple.Exit")]
    df_ss2$MSCI_ <- df_ss2$MSCI.Multiple.Exit - 1
    # df_ss2$RVPI_ <- df_ss2$RVPI - 1
    df_ss2 <- df_ss2[,c("P2C.multi1","Holding_Period","Time2Exit","ZombieStage","ML_HYOAS.quarter","MSCI_")]
    x <- model.matrix(P2C.multi1 ~ .-1,data=df_ss2)
    y <- df_ss2[,1]
    
    # Ridge: L2 Regularization
    cv.ridge <- glmnet::cv.glmnet(x, y, family='gaussian', alpha=0, parallel=TRUE, standardize=TRUE,intercept=TRUE,weights=odw1(df_ss2$P2C.multi1))
    # LASSO: L1 Regularization
    cv.lasso <- glmnet::cv.glmnet(x, y, family='gaussian', alpha=1, parallel=TRUE, standardize=TRUE,intercept=TRUE,weights=odw1(df_ss2$P2C.multi1))
    
    re[[paste("p2c.ridge",type,sep=".")]] <- cv.ridge
    re[[paste("p2c.lasso",type,sep=".")]] <- cv.lasso
    
    rm(type,x,y,cv.ridge,cv.lasso,df_ss1,df_ss2)
  }
  
  
  
  
  re.bag <- function(split=10, runs=10, trunc_low= 0.15, type = "VC"){
    geelist <- list()
    censReglist <- list()
    lmroblist <- list()
    lmroblist_trunc <- list()
    lmerlist <- list()
    
    for(j in 1:runs){
      print(paste(">>>> Run",j))
      
      if(type== "VC"){  df_ss <- vc_ss}
      if(type== "BO"){  df_ss <- bo_ss}
      
      # randomize & split
      ids <- levels(as.factor(df_ss$Fund_Emi_ID))
      ids <- sample(ids)
      stepsize <- round(length(ids)/split,0)
      steps <- seq(1,length(ids),stepsize)
      steps <- steps[1:split]
      
      for(i in steps){
        df_ss1 <- df_ss
        if(i == steps[split]){
          end <- length(ids)
        }else{
          end <- i+stepsize-1
        }
        df_ss1 <- df_ss1[df_ss1$Fund_Emi_ID %in% ids[i:end],]
        # Do Regressions
        if(sum(df_ss1$ZombieStage != 0) > 0 & sum(df_ss1$P2C.multi1 <= 0) > 0){
          # print(i)
          p_df_ss <- plm::pdata.frame(df_ss1,c("Fund_Emi_ID","FundInv_Quarter"))
          
          if(type == "BO"){
            p2c.lmrob <- robustbase::lmrob(formula_multiple_bo2, data=df_ss1, method = "SMDM",weights = log(FMV_C2come))
            p2c.lmrob_trunc <- robustbase::lmrob(formula_multiple_bo2, data=df_ss1[df_ss1$P2C.multi1 > trunc_low,], method = "SMDM")
            p2c.lmer <- lme4::lmer(formula_multiple_bo_lmer, data=df_ss1,weights= log(FMV_C2come))
            p2c.geese <- geepack::geeglm(formula_multiple_bo2, id= Fund_Emi_ID, data=df_ss1,corstr= "ar1")
            p2c.censReg <- censReg::censReg(formula_multiple_bo2, left=0, data=p_df_ss)
          }else{ # VC
            p2c.lmrob <- robustbase::lmrob(formula_multiple_vc2, data=df_ss1, method = "SMDM",weights = log(FMV_C2come))
            p2c.lmrob_trunc <- robustbase::lmrob(formula_multiple_vc2, data=df_ss1[df_ss1$P2C.multi1 > trunc_low,], method = "SMDM")
            p2c.lmer <- lme4::lmer(formula_multiple_vc_lmer, data=df_ss1,weights= log(FMV_C2come))
            p2c.geese <- geepack::geeglm(formula_multiple_vc2, id= Fund_Emi_ID, data=df_ss1,corstr= "ar1")
            p2c.censReg <- censReg::censReg(formula_multiple_vc2, left=0, data=p_df_ss)
          }
          # Fill Lists
          geelist[[paste(i,j,sep="_")]] <- c(p2c.geese$coefficients, gamma= as.numeric(p2c.geese$geese$gamma), alpha= as.numeric(p2c.geese$geese$alpha))
          censReglist[[paste(i,j,sep="_")]] <- p2c.censReg$estimate
          lmroblist[[paste(i,j,sep="_")]] <- p2c.lmrob$coefficients
          lmroblist_trunc[[paste(i,j,sep="_")]] <- p2c.lmrob_trunc$coefficients
          lmerlist[[paste(i,j,sep="_")]] <- c(colMeans(coef(p2c.lmer)$Fund_Emi_ID),sd_intercept= sd(unlist(coef(p2c.lmer)$Fund_Emi_ID["(Intercept)"])))
        }else{print(paste("Exclude",i))}
      }
    }
    
    out_frame <- function(inlist){
      models <- data.frame(do.call(rbind,inlist))
      models <- models[complete.cases(models),]
      bag_reg <- data.frame(Mean= colMeans(models),
                            Median= apply(models,2,median),
                            SD= apply(models,2,sd))
      bag_reg$t_value <- bag_reg$Mean/bag_reg$SD
      return(bag_reg)
    }
    
    return(list(geeglm= out_frame(geelist),
                lmer= out_frame(lmerlist),
                censReg= out_frame(censReglist),
                lmrob= out_frame(lmroblist),
                lmrob_trunc= out_frame(lmroblist_trunc)))
  }
  # system.time( bag1 <- re.bag() )
  
  
  
  # Ensemble Learning Approach
  re.bag2 <- function(runs=25, sample_size=200, type = "BO",downweighting= FALSE){
    censReglist <- list()
    censRegXlist <- list()
    lmroblist <- list()
    olslist <- list()
    ridge_list <- list()
    lasso_list <- list()
    
    if(type== "VC"){
      df_ss <- vc_ss
      formula_reg <- formula("P2C.multi1 ~ Holding_Period + Time2Exit + ZombieStage + I(RVPI-1) + ML_HYOAS.quarter + I(MSCI.Multiple.Exit-1)")
    }
    if(type== "BO"){  
      df_ss <- bo_ss
      # formula_reg <- formula("P2C.multi1 ~ Holding_Period + Time2Exit + ZombieStage + ML_HYOAS.quarter + I(MSCI.Multiple.Exit-1)")
      formula_reg <- formula("P2C.multi1 ~ Time2Exit + I(RVPI-1) + ML_HYOAS.quarter + I(MSCI.Multiple.Exit-1)")
    }
    
    
    for(j in 1:runs){
      # randomized samples (option: outlier down-weighting within ID groups)
      if(downweighting){
        df_ss1 <- as.data.frame(df_ss %>% group_by(Fund_Emi_ID) %>% sample_n(sample_size, replace=TRUE, weight= odw1(P2C.multi1) ))
      }else{
        df_ss1 <- as.data.frame(df_ss %>% group_by(Fund_Emi_ID) %>% sample_n(sample_size, replace=TRUE))
      }
      # df_ss1 <- as.data.frame(df_ss %>% group_by(Fund_Emi_ID) %>% sample_n(sample_size, replace=TRUE, weight= log(FMV_C2come)  ))
      # df_ss1 <- as.data.frame(df_ss %>% sample_n(sample_size*250,replace=TRUE,weight= odw1(P2C.multi1)))
      # df_ss1 <- as.data.frame(df_ss %>% sample_n(sample_size*250,replace=TRUE,weight= log(FMV_C2come)))
      
      # LASSO & Ridge Data Preparation
      df_ss2 <- df_ss1[,c("P2C.multi1","Holding_Period","Time2Exit","ZombieStage","RVPI","ML_HYOAS.quarter","MSCI.Multiple.Exit")]
      df_ss2$MSCI_ <- df_ss2$MSCI.Multiple.Exit - 1
      df_ss2$RVPI_ <- df_ss2$RVPI
      df_ss2 <- df_ss2[,c("P2C.multi1","Holding_Period","Time2Exit","ZombieStage","RVPI_","ML_HYOAS.quarter","MSCI_")]
      x <- as.matrix(df_ss2[,2:ncol(df_ss2)])
      y <- df_ss2[,1]
      
      if(sum(df_ss1$ZombieStage != 0) > 0 & sum(df_ss1$P2C.multi1 <= 0) > 0){
        # Do Regressions
        p2c.ols <- lm(formula_reg, data=df_ss1)
        p2c.lmrob <- robustbase::lmrob(formula_reg, data=df_ss1, method = "SMDM")
        p2c.censReg <- censReg::censReg(formula_reg, left=0, right=Inf, data=df_ss1)
        p2c.censRegX <- censReg::censReg(formula_reg, left=0.5, right=4, data=df_ss1)
        # Ridge: L2 Regularization
        cv.ridge <- glmnet::cv.glmnet(x, y, family='gaussian', alpha=0, parallel=TRUE, standardize=TRUE)
        # LASSO: L1 Regularization
        cv.lasso <- glmnet::cv.glmnet(x, y, family='gaussian', alpha=1, parallel=TRUE, standardize=TRUE)
        
        # Fill Lists
        censReglist[[paste("run",j,sep="_")]] <- p2c.censReg$estimate
        censRegXlist[[paste("run",j,sep="_")]] <- p2c.censRegX$estimate
        lmroblist[[paste("run",j,sep="_")]] <- p2c.lmrob$coefficients
        olslist[[paste("run",j,sep="_")]] <- p2c.ols$coefficients
        ridge_list[[paste("run",j,sep="_")]] <- t(as.matrix(coef(cv.ridge, s=cv.ridge$lambda.min)))
        lasso_list[[paste("run",j,sep="_")]] <- t(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)))
        
      }else{print(paste("Exclude",i))}
    }
    
    out_frame <- function(inlist){
      models <- data.frame(do.call(rbind,inlist))
      models <- models[complete.cases(models),]
      bag_reg <- data.frame(Mean= colMeans(models),
                            Median= apply(models,2,median),
                            SD= apply(models,2,sd))
      bag_reg$'Mean/SD' <- bag_reg$Mean/bag_reg$SD
      bag_reg <- round(bag_reg,4)
      return(bag_reg)
    }
    
    return(list(censReg= out_frame(censReglist),
                censRegX= out_frame(censRegXlist),
                lasso= out_frame(lasso_list),
                ridge= out_frame(ridge_list),
                lmrob= out_frame(lmroblist),
                ols= out_frame(olslist)))
  }
  # system.time( bo.bag <- re.bag2(type = "BO",runs=1000,sample_size=1) )
  # saveRDS(bo.bag ,"z_bo.bag3_ss1.rds")
  bo.bag <- readRDS("z_bo.bag3_ss1.rds")
  # system.time( vc.bag <- re.bag2(type = "VC",runs=1000,sample_size=1) )
  # saveRDS(vc.bag ,"z_vc.bag3_ss1.rds")
  vc.bag <- readRDS("z_vc.bag3_ss1.rds")
  
  re.bag.coefs <- function(modelz,models_bag= c("censReg","lmrob","ols"),print_it=FALSE){
    all_row_names <- c("Intercept","Holding_Period","Time_To_Exit","Zombie_Stage","RVPI_1","High_Yield_Spread","MSCI_1")
    len_df <- min(sapply(modelz,nrow))
    mean_list <- list()
    for(i in models_bag){
      df <- modelz[[i]]
      df <- df[1:len_df,]
      
      if(len_df == 7){ # VC
        row_names <- all_row_names
        mean_list[[i]] <- df$Mean
      }else{ # BO
        row_names <- c("Intercept","Time_To_Exit","RVPI_1","High_Yield_Spread","MSCI_1")
        mean_list[[i]] <- c(df$Mean[1],0,df$Mean[2],0,df$Mean[3:5])
      }
      
      rownames(df) <- row_names
      modelz[[i]] <- df
    }
    if(print_it) print(modelz[names(modelz) %in% models_bag])
    
    bag_average <- data.frame(do.call(rbind,mean_list))
    colnames(bag_average) <- all_row_names
    bag_average <- colMeans(bag_average)
    
    return(bag_average)
  }
  re.bag.coefs(bo.bag,print_it=TRUE)
  re.bag.coefs(vc.bag,print_it=TRUE)
  
  
  rm(formula_multiple_bo,formula_multiple_vc,formula_multiple_bo2,formula_multiple_vc2,p2c.threshold)
}
## x-2) Multiple Prediction
if(FALSE){
  pred.bag2 <- function(modelz,df){
    out2 <- list()
    all_models <- c("censReg","lmrob","ols")
    all_models_bagging <- c(all_models,"bagging")
    
    for(j in all_models_bagging){
      if(j == "bagging"){
        len <- length(modelz$ols$Mean)
        model <- re.bag.coefs(modelz,all_models)
      }else{
        model <- modelz[[j]]$Mean
      }
      
      if(length(model) < 7){
        print(paste("5 Parameter",j))
        betas <- c(model[1],0,model[2],0,model[3:5])
      }else{
        print(paste("7 Parameter",j))
        betas <- model[1:7]
      }
      
      out <- list()
      
      for(i in 1:nrow(df)){
        X_vector <- c(1,
                      df$Holding_Period[i],
                      df$Time2Exit[i],
                      df$ZombieStage[i],
                      (df$RVPI[i]-1),
                      df$ML_HYOAS.quarter[i],
                      (df$MSCI.Multiple.Exit[i]-1))
        
        out[i] <- sum(betas * X_vector)
      }
      
      out2[[j]] <- unlist(out)
    }
    
    out3 <- data.frame(Obs=(df$P2C.multi1),
                       Pred_Bag= out2[["bagging"]],
                       Pred_Tobit= out2[["censReg"]],
                       Pred_Rob= out2[["lmrob"]],
                       Pred_OLS= out2[["ols"]])
    
    out3[,c("Resi_Bag","Resi_Tobit","Resi_Rob","Resi_OLS")] <- out3$Obs - out3[,2:5]
    # out3$Pred_Tobit_zm <- out3$Pred_Tobit + mean(out3$Resi_Tobit)
    # out3$Pred_Rob_zm <- out3$Pred_Rob + mean(out3$Resi_Robust)
    # out3$Pred_OLS_zm <- out3$Pred_OLS + mean(out3$Resi_OLS)
    
    return(out3)
  }
  system.time(pre_bo <- pred.bag2(modelz= bo.bag,df= bo_ss))
  system.time(pre_vc <- pred.bag2(modelz= vc.bag,df= vc_ss))
  
  
  bo_ss[,c("Pred_Multiple","Resi_Multiple")] <- pre_bo[,c("Pred_Bag","Resi_Bag")]
  vc_ss[,c("Pred_Multiple","Resi_Multiple")] <- pre_vc[,c("Pred_Bag","Resi_Bag")]
  
  
  # realized quantile multiple
  rqm <- function(df){
    for(i in 1:nrow(df)){
      region <- df$Fund_Region[i]
      fund_type <- df$Fund_InvestTypes[i]
      
      realized <- df$P2C.multi1[i]
      estimate <- df$Pred_Multiple[i]
      residualz <-df$Resi_Multiple
      
      np_distribution <- estimate + residualz
      
      real_quantile <- sum(realized > np_distribution)/length(np_distribution)
      
      df[i,"RQ"] <- real_quantile
    }
    return(df$RQ)
  }
  system.time(bo_ss$RQ_Multiple <- rqm(bo_ss))
  system.time(vc_ss$RQ_Multiple <- rqm(vc_ss))
  
}
## b) Timing Regression & Prediction --------
if(FALSE){
  cox.model <- function(df.bo,df.vc){
    
    bo_vc <- rbind(df.bo, df.vc)
    # Exit Events are truncated by Filtering (!!!) <-> correction of this error
    d = as.data.frame(aggregate(bo_vc$FundInv_Quarter,by=list(bo_vc$Fund_Emi_ID),max))
    colnames(d) <- c("Fund_Emi_ID","MaxDateID")
    bo_vc <- merge(bo_vc,d,by="Fund_Emi_ID",all.x=TRUE)
    bo_vc$EvLast <- as.numeric(bo_vc$FundInv_Quarter == bo_vc$MaxDateID & !is.na(bo_vc$P2C.multi1))
    print(sum(levels(bo_vc$Fund_Emi_ID)))
    print(sum(bo_vc$EvLast))

    # internal covariates can NOT be used in Cox regression
    bo_vc$RVPI.trans2 <- dgamma(bo_vc$RVPI,shape=2,scale=1)
    bo_vc$DPI <- pmax(0,bo_vc$MOIC.dyn - bo_vc$RVPI)
    
    # sun(bo_vc$ML_HYOAS.quarter)
    # sun(bo_vc$FundAgeAI)
    re.Cox <-  survival::coxph(survival::Surv(time= T1, time2= T2, event= Ev) ~ 
                                 strata(Fund_InvestTypes) +
                                 FundAgeAI +
                                 # GLOQuarter +
                                 ML_HYOAS.quarter,
                                 # MOICtrans2,
                                 # RVPI.trans2 + DPI,
                                data = bo_vc[bo_vc$Fund_InvestTypes %in% c("BO","VC"),], ties="efron")
    print(summary(re.Cox))
    # cox.zph(re.Cox) ; survfit(re.Cox)
    
    # Cumulative basehazard
    bh.Cox <- survival::basehaz(re.Cox,centered=FALSE) # Covariates == 0
    bh.Cox$S.base <- exp(-1*bh.Cox$hazard)
    bh.Cox.Mean <- survival::basehaz(re.Cox,centered=TRUE) # Covariates == 0
    bh.Cox.Mean$S.base <- exp(-1*bh.Cox.Mean$hazard)
    
    out <- list(Regression= re.Cox, Base_Hazard= bh.Cox,Base_Hazard_Mean=bh.Cox.Mean, Data_Frame= bo_vc)
    
    return(out)
  }
  system.time( cox_model <- cox.model(bo_ss,vc_ss) )
  cox_fac <- function(df,co_mo){
    out <- list()
    for(i in 1:nrow(df)){
      out[i] <-   sum(co_mo$Regression$coefficients * df[i, c("FundAgeAI", "ML_HYOAS.quarter", "Moic.trans2")])
    }
    return(as.numeric(unlist(out)))
  }
  cox.it <- function(cox_model){
    re.Cox <- cox_model$Regression
    bh.Cox <- cox_model$Base_Hazard
    bo_vc <- cox_model$Data_Frame
    
    mean.Cox <- function(t=0.5,strata="BO",CoxF=1){
      if(t<0.5) t=0.5
      df1 <- bh.Cox[bh.Cox$strata==strata,]
      df1$S <- exp(-CoxF*df1$hazard)
      a= diff(df1[df1$time > t,"time"])
      b= df1[df1$time > t & df1$t < df1[nrow(df1),"time"],"S"]
      c= df1[last(which(df1$time <= t)),"S"]
      x= sum(a*b)/c
      return(x)
    }
    quant.Cox <- function(t_start,t_end,strata="BO",CoxF=1){
      df1 <- bh.Cox[bh.Cox$strata==strata,]
      df1$S <- exp(-CoxF*df1$hazard)
      
      if(t_start > df1$time[1]){
        P_start <- min(df1[df1$time < t_start, "S"])
        P_end <- min(df1[df1$time < t_end, "S"])
      }else{
        P_start <- df1$S[1]
        P_end <- min(df1[df1$time < t_end, "S"])
      }
      
      y <- (1-P_end/P_start)
      
      return(y)
    }
    
    
    bo_vc$RLT <- as.numeric(bo_vc$Exit_Date - bo_vc$FundInv_Quarter)/365.25
    bo_vc$predict.cox.np <- cox_fac(bo_vc,cox_model)
    
    for(i in 1:nrow(bo_vc)){
      asset.age <- bo_vc[i,"Holding_Period"]
      type <- as.character(bo_vc[i,"Fund_InvestTypes"])
      coxfac <- exp(bo_vc[i,"predict.cox.np"])
      rml <- bo_vc$RLT[i]
      if(!any(is.na(c(asset.age,type,coxfac)))){
        bo_vc[i,"ERLT"] <- mean.Cox(t=asset.age,
                                    strata=type,
                                    CoxF=coxfac)
        bo_vc[i,"RQ_RLT"] <- quant.Cox(t_start = asset.age,
                                       t_end = asset.age + rml,
                                       strata=type,
                                       CoxF=coxfac)
      }
    }
    
    bo_vc$Resi_Timing <- bo_vc$RLT - bo_vc$ERLT
    
    # bo1 <- bo_vc[bo_vc$Fund_InvestTypes == "BO",]
    # vc1 <- bo_vc[bo_vc$Fund_InvestTypes == "VC",]
    # return(list(BO= bo1,VC= vc1))
    return(bo_vc)
    
  }
  system.time( bo_vc <- cox.it(cox_model) )
  cox_weibull <- function(df){
    # Estimate Weibull hazard with non-parametric simulation
    np_sim <- function(n=100,strata){
      df <- df[df$strata==strata,]
      x <- list()
      for(i in 1:n){    # non-parametric life-time simulation
        pos <- tail(which(runif(1) < df$S.base),1)
        if(length(pos) < 1) pos <- 1
        x[i] <- df$time[pos]
      }
      return(unlist(x))
    }
    
    wb.bo <- MASS::fitdistr(np_sim(30000,"BO"),densfun="weibull")
    wb.vc <- MASS::fitdistr(np_sim(30000,"VC"),densfun="weibull")
    
    
    out <- list()
    out$BO <- c(wb.bo$estimate["scale"], wb.bo$estimate["shape"])
    out$VC <- c(wb.vc$estimate["scale"], wb.vc$estimate["shape"])
    
    return(out)
  }
  system.time( cox_wb_zero  <- cox_weibull(cox_model$Base_Hazard) )
  system.time( cox_wb_mean  <- cox_weibull(cox_model$Base_Hazard_Mean) )
  S_wb <- function(x,type,cox_wb1){
    y=exp(-(x/cox_wb1[[type]]["scale"])^cox_wb1[[type]]["shape"])
    return(as.numeric(y))}
  plot_cox.weibull <- function(eps=FALSE,df_in=cox_model$Base_Hazard_Mean){
    if(eps){
      setEPS() ; postscript("Timing Cox Weibull Base Survival Function.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=5)
    }
    
    cox_wb= cox_weibull(df_in)  
    bo_col <- "royalblue1" ; vc_col <- "maroon1"
    
    df <- df_in
    par(mfrow=c(1,1),mar=c(4.5,4.5,2,2),cex=1.5)
    plot(df$time[df$strata == "BO"],df$S.base[df$strata == "BO"],lty=3,type="l",xlim=c(0,20),ylim=c(0,1),col=bo_col,xlab="Time",ylab="Survival Function")
    curve(S_wb(x,"BO",cox_wb),add=TRUE,col=bo_col)
    lines(df$time[df$strata == "VC"],df$S.base[df$strata == "VC"],lty=3,col=vc_col)
    curve(S_wb(x,"VC",cox_wb),add=TRUE,col=vc_col)
    abline(h=0,v=0,col="grey",lty=2)
    legend("topright",bty="n",legend=c("BO",
                                       paste("Scale",round(cox_wb$BO["scale"],2)),
                                       paste("Shape",round(cox_wb$BO["shape"],2)),
                                       "VC",
                                       paste("Scale",round(cox_wb$VC["scale"],2)),
                                       paste("Shape",round(cox_wb$VC["shape"],2))),
           col=c(bo_col,NA,NA,vc_col,NA,NA),lty=c(1,NA,NA,1,NA,NA))
    
    if(eps){ dev.off() }
  }
  plot_cox.weibull(df_in=cox_model$Base_Hazard_Mean)
  survival::cox.zph(cox_model$Regression)
  summary(cox_model$Regression)
  
  quantile_weibull <- function(t_now=0,T_exit=0,cofa=0,type="BO"){
    q_bo <- S_wb(T_exit,"BO",cox_wb_zero)^exp(cofa) / S_wb(t_now,"BO",cox_wb_zero)^exp(cofa)
    q_vc <- S_wb(T_exit,"VC",cox_wb_zero)^exp(cofa) / S_wb(t_now,"VC",cox_wb_zero)^exp(cofa)
    out <- ifelse(type == "BO",q_bo,q_vc)
    out[out < 0 | out > 1] <- NA
    return(out)
  }

  bo_vc$RQ_RLT.weibull <- quantile_weibull(bo_vc$Holding_Period,
                                           as.numeric((bo_vc$MaxDate-bo_vc$Investment_Date)/365.25),
                                           cox_fac(bo_vc,cox_model),
                                           bo_vc$Fund_InvestTypes)
  
  Weibull_Roseblatt <- function(eps=FALSE){
    if(eps){
      setEPS() ; postscript("Timing_Rosenblatt.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=5)
    }
    bo_col <- "royalblue1" ; vc_col <- "maroon1"
    
    par(mar=c(3,3,2,1),mfrow=c(1,2),oma=c(2.5,2,1,1))
    for(type in c("BO","VC")){
      hi_col <- ifelse(type=="BO",bo_col,vc_col)
      hist(bo_vc$RQ_RLT.weibull[bo_vc$Fund_InvestTypes == type],freq = FALSE,
           main=type,border=hi_col,lty=3,ylab=NA,xlab=NA,ylim=c(0,1.5),xlim=c(0,1))
      abline(a=0,b=1,col="darkgray",lwd=1)
      wb_ecdf <- ecdf(bo_vc$RQ_RLT.weibull[bo_vc$Fund_InvestTypes == type])
      curve(wb_ecdf,add=TRUE,col="darkgreen",lwd=1)
      legend("topright",bty="n",cex=1.3,col=c("darkgreen","darkgray"),legend = c("empirical","theoretical"),lty=1)
      abline(h=1,col="black",lty=2)
    }
    mtext(latex2exp::TeX('$\\S^{(wb)}(d_i) (S^{(wb)}(t_i))^{-1}'),side=1,outer= TRUE,line=1,cex=1.7)
    mtext("ECDF or Density",side=2,outer= TRUE,line=0.5,cex=1.5)
    
    
    if(eps){ dev.off() }
  }
  Weibull_Roseblatt()
  
  sim_wb <- function(u= 0.9, type="BO",cofa=0,cox_wb=cox_wb_zero){
    shape <- cox_wb[[type]]["shape"]
    scale <- cox_wb[[type]]["scale"]
    o1 <- scale * (- log(u^(1/exp(cofa))))^(1/shape)
    return(as.numeric(o1))
  }
  # Check if quantile_weibull(...) and sim_wb(...) are consistent
  quantile_weibull(0,sim_wb(0.8,type="BO",cofa=-0.5),-0.5,"BO")
}
if(FALSE){
  # non-parametric Cox simulation
  sim_cox2 <- function(t=0.5,strata="BO",CoxF=1,real_quantile=NA,bh.NP5=cox_model$Base_Hazard){
    if(t < 0.5) t <- 0.5
    if(t > 15) t <- 15
    
    df1 <- bh.NP5[bh.NP5$strata==strata,]
    df1$S <- exp(-CoxF*df1$hazard)
    c= df1[last(which(df1$time <= t)),"S"]
    
    if(is.na(real_quantile)){
      d=  runif(1) * c
    }else{
      d= real_quantile * c # realized quantile correction
    }
    
    if(d < min(df1$S)){
      rlt <- max(0.250002,max(df1$time) - t)
    }else{
      wd <- which(df1$S <= d)
      e= df1[first(wd),"time"]
      rlt <- max(0.250002,e - t)
    }
    
    return(rlt)
  }
  
}
if(FALSE){
  # Cox Regression with time-invariant covariates using G.sum
  G.sum$Timing <- pmax(0.25,as.numeric(G.sum$MaxDate - G.sum$MinDate)/365.25)
  G.sum$Event  <- ifelse(is.na(G.sum$HoPi),0,1)
  G.sum$FundAge <- pmax(0,as.numeric(format(G.sum$MinDate,"%Y")) - G.sum$Vintage)
  
  re.Cox <-  survival::coxph(survival::Surv(time= Timing,event= Event) ~ 
                               strata(Type) + FundAge,
                              data = G.sum, ties="efron")
  print(summary(re.Cox))
  re.Cox.np <- eha::coxreg(survival::Surv(time= Timing,event= Event) ~ 
                            strata(Type) + FundAge,
                            method="efron",data = G.sum)
  print(summary(re.Cox.np)) ; plot(re.Cox.np)
  re.Cox.wb <- eha::phreg(survival::Surv(time= Timing,event= Event) ~ 
                            strata(Type) + FundAge,
                          dist="weibull",data = G.sum)
  print(summary(re.Cox.wb)) ; plot(re.Cox.wb)
  re.WB <- survival::survreg(survival::Surv(time= Timing,event= Event) ~ 
                                   strata(Type) + FundAge,
                                   dist="weibull",data = G.sum)
  print(summary(re.WB))
  re.lolo <- survival::survreg(survival::Surv(time= Timing,event= Event) ~ 
                              strata(Type) + FundAge,
                            dist="loglogistic",data = G.sum)
  print(summary(re.lolo))
  
  
  S_aft <- function(t){
    mue_aft <- exp(2.0986-0.4928-0.0394*7)
    scale_aft <- 0.611
    
    S <- exp(-(t/mue_aft)^(1/scale_aft))
    return(S)
  }
  S_aft(0:15)
  
  
  
}
if(TRUE){
cox.model2 <- function(G.p2c){
  # delete young entries (bias correction?)
  # G.p2c <- G.p2c[G.p2c$Investment_Date < as.Date("2010-01-01"),]
  G.p2c <- G.p2c[order(G.p2c$Fund_Emi_ID,G.p2c$FundInv_Quarter),]
  
  # equidistant time-grid
  G.p2c$T1_r <- round(G.p2c$T1*4,0)/4
  G.p2c$T2_r <- round(G.p2c$T2*4,0)/4
  
  outlist <- list()
  
  for(Fund_Types in c("BO","VC")){
    re.Cox <-  survival::coxph(survival::Surv(time= T1_r, time2= T2_r, event= Ev) ~ 
                                 FundAgeAI +
                                 # GLOQuarter +
                                 ML_HYOAS.quarter +
                                 # MSCI.multi.hist +
                                 MSCI.qrtly.return,
                               data = G.p2c[G.p2c$Fund_InvestTypes %in% Fund_Types,], ties="efron")
    
    # Cumulative basehazard
    bh.Cox <- survival::basehaz(re.Cox,centered=FALSE) # Covariates == 0
    bh.Cox$S.base <- exp(-1*bh.Cox$hazard)
    bh.Cox$hazard_rate <- c(0,diff(bh.Cox$hazard))
    
    bh.Cox.Mean <- survival::basehaz(re.Cox,centered=TRUE) # Covariates == 0
    bh.Cox.Mean$S.base <- exp(-1*bh.Cox.Mean$hazard)
    bh.Cox.Mean$hazard_rate <- c(0,diff(bh.Cox.Mean$hazard))
    
    outlist[[Fund_Types]] <- list(Regression= re.Cox, BaseHaze_Zero= bh.Cox, BaseHaze_Mean= bh.Cox.Mean)
  }
  

  '
  par(mfrow=c(2,2))
  plot(bh.Cox[bh.Cox$strata =="BO",c(2,4)],main="BO zero",ylim=c(0,1)) ; abline(h=0,v=0,col="grey")
  plot(bh.Cox.Mean[bh.Cox.Mean$strata =="BO",c(2,4)],main="BO mean",ylim=c(0,1)) ; abline(h=0,v=0,col="grey")
  plot(bh.Cox[bh.Cox$strata =="VC",c(2,4)],main="VC zero",ylim=c(0,1)) ; abline(h=0,v=0,col="grey")
  plot(bh.Cox.Mean[bh.Cox.Mean$strata =="VC",c(2,4)],main="VC mean",ylim=c(0,1)) ; abline(h=0,v=0,col="grey")
  par(par_default)
  '
  
  out <- outlist
  return(out)
}
system.time( cox_model <- cox.model2(G.p2c) )
summary(cox_model$BO$Regression)
summary(cox_model$VC$Regression)



public_data <- trindex[,c("Date","MSCI.World.Net.Return.Daily","ML_HYOAS")]
public_data$MSCI_monthly_return <- c(NA,diff(public_data$MSCI.World.Net.Return.Daily) / public_data$MSCI.World.Net.Return.Daily[1:(nrow(public_data)-1)])
as.Date(levels(as.factor(G.p2c$FundInv_Quarter))) %in% public_data$Date


Li_WBcox <- function(from_Date, to_Date,
                          beta_MSCI, beta_HYS, beta_FundAge,
                          scale_wb, shape_wb,
                          FundAge,
                          public_data){
  from_Date <- as.Date(from_Date)
  if(is.na(to_Date)){
    to_Date <- as.Date("2016-12-31")
    Censoring <- TRUE
  }else{
    to_Date <- as.Date(to_Date)
    Censoring <- FALSE
  }
  public_data <- public_data[public_data$Date >= from_Date & public_data$Date <= to_Date,]
  public_data$t <- as.numeric(public_data$Date - public_data$Date[1])/365.25
  public_data$FundAge <- FundAge
  
  public_data$base_haze_rate_wb_exact <- c(0,diff( (public_data$t/scale_wb)^shape_wb ))
  
  public_data$expBX <- exp( as.matrix(public_data[,c("MSCI_monthly_return","ML_HYOAS","FundAge")]) %*% 
                              c(beta_MSCI,beta_HYS,beta_FundAge) )
  
  public_data$Cum_Haze_exact <- cumsum(as.numeric(public_data$base_haze_rate_wb_exact) * public_data$expBX)
  public_data$Surv_WB <- exp(-public_data$Cum_Haze_exact)

  if(Censoring){
    out <- public_data$Surv_WB[nrow(public_data)]
  }else{
    out <- - diff(public_data$Surv_WB)
    out <- out[length(out)]
  }
  return(out)
}
Li_WBcox("2016-01-31","2016-12-31",1,-2,0.05,2,0.5,3,public_data)
Li_WBcox("2015-01-31",NA,1,-2,0.05,2,0.5,3,public_data)



LoLi_WBcox <- function(par, data){
  x <- par  
  loglikelihood <- list()
  for(i in 1:nrow(data)){
    L <- Li_WBcox(from_Date = data$InvDate[i],
                  to_Date = data$ExitDate[i],
                  x[1], x[2], x[3], x[4], x[5],
                  FundAge = max(0, data$YearInvest[i] - data$Vintage[i]),
                  public_data = public_data)
    
    loglikelihood[i] <- log(L)
  }
return(sum(as.numeric(loglikelihood)))
}
LoLi_WBcox2 <- function(par, data){
  x <- c(0,0,0,par)
  
  loglikelihood <- list()
  for(i in 1:nrow(data)){
    L <- Li_WBcox(from_Date = data$InvDate[i],
                  to_Date = data$ExitDate[i],
                  x[1], x[2], x[3], x[4], x[5],
                  FundAge = max(0, data$YearInvest[i] - data$Vintage[i]),
                  public_data = public_data)
    
    loglikelihood[i] <- log(L)
  }
  return(sum(as.numeric(loglikelihood)))
}

G.sum$Timing <- ifelse(is.na(G.sum$HoPi),
                       as.numeric(as.Date("2016-12-31") - G.sum$InvDate)/365.25,
                       as.numeric(G.sum$ExitDate - G.sum$InvDate)/365.25)
G.sum$RightCensored <- ifelse(is.na(G.sum$HoPi),1,0)
data_LoLi <- G.sum[G.sum$InvDate < G.sum$ExitDate | is.na(G.sum$ExitDate),]
# data_LoLi <- data_LoLi[data_LoLi$Timing > 0,] # our model can handle zero time for censored events
sun(data_LoLi$Timing)

table(is.na(data_LoLi$HoPi),data_LoLi$Type) # exit vs right censored events

system.time( loli  <- LoLi_WBcox(c(-0.6,-2,0.05,2,0.5),data_LoLi) )
# sun(loli) ; sum(loli) ; which(is.infinite(loli))

Para_Cox_Weibull <- list()
for(fund_type in c("BO","VC")){
  data_sample <- data_LoLi[data_LoLi$Type == fund_type,]

  system.time( Para_Cox_Weibull$fit[[fund_type]] <- optim(par= c(1,-2,0.05,2,0.5), 
                            fn= LoLi_WBcox, 
                            data= data_sample,
                            method= "L-BFGS-B", 
                            lower=c(rep(-Inf,3),0,0),
                            upper=c(Inf,Inf,Inf,Inf,Inf),
                            control=list("fnscale"=-1), # makes it a min problem
                            hessian=TRUE))
  
  system.time( Para_Cox_Weibull$fit2[[fund_type]] <- optim(par= c(2,0.5), 
                             fn= LoLi_WBcox2,
                             data= data_sample,
                             method= "L-BFGS-B", 
                             lower=c(0,0),
                             upper=c(Inf,Inf),
                             control=list("fnscale"=-1), # makes it a min problem
                             hessian=TRUE))
}


# saveRDS(Para_Cox_Weibull, "z_Para_Cox_Weibull.rds")
Para_Cox_Weibull <- readRDS("z_Para_Cox_Weibull.rds")


SEandAIC <- function(fit){
  # 95% Conf_Intervals
  # https://stats.stackexchange.com/questions/27033/in-r-given-an-output-from-optim-with-a-hessian-matrix-how-to-calculate-paramet
  fisher_info<-solve(-fit$hessian) # postive hessian since we minimize -log(Likelihood) ??
  prop_sigma<-sqrt(diag(fisher_info))
  
  result <- data.frame(MLE= fit$par, SE= prop_sigma)
  result$t_value <- result$MLE / result$SE
  result$p_value <- pnorm(-abs(result$t_value))
  print(round(result,4))
  #prop_sigma<-diag(prop_sigma)
  #upper<-fit$par+1.96*prop_sigma
  #lower<-fit$par-1.96*prop_sigma
  #interval<-data.frame(value=fit$par, upper=upper, lower=lower)
  
  fit$Coefs <- result
  fit$AIC <-  2 * length(par) -  2 * fit$value
  # fit$BIC <- log(n) * length(par) -  2 * fit$value
  
  return(fit)
}

Para_Cox_Weibull$fit$BO <- SEandAIC(Para_Cox_Weibull$fit$BO)
Para_Cox_Weibull$fit$VC <- SEandAIC(Para_Cox_Weibull$fit$VC)
Para_Cox_Weibull$fit2$BO <- SEandAIC(Para_Cox_Weibull$fit2$BO)
Para_Cox_Weibull$fit2$VC <- SEandAIC(Para_Cox_Weibull$fit2$VC)

Para_Cox_Weibull$fit$VC$AIC / Para_Cox_Weibull$fit2$VC$AIC
Para_Cox_Weibull$fit$BO$AIC / Para_Cox_Weibull$fit2$BO$AIC

fitWeib <- survreg(Surv(Timing, RightCensored) ~ 1, dist="weibull", data=data_LoLi)
summary(fitWeib)


plot_cox.weibull <- function(eps=FALSE,para_in = Para_Cox_Weibull, non_para_in= cox_model){
  if(eps){
    setEPS() ; postscript("Timing Cox Weibull Base Survival Function 2.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=5)
  }
  bo_col <- "royalblue1" ; vc_col <- "maroon1"
  
  sc_bo <- para_in$BO$par[4]
  sh_bo <- para_in$BO$par[5]
  sc_vc <- para_in$VC$par[4]
  sh_vc <- para_in$VC$par[5]
  sc_bo * gamma(1+1/sh_bo) # mean Weibull
  sc_vc * gamma(1+1/sh_vc) # mean Weibull

  
  par(mfrow=c(1,1),mar=c(4.5,4.5,2,2),cex=1.5)
  curve(1- pweibull(x, scale = sc_bo, shape = sh_bo),0,20,col=bo_col,xlab="Time",ylab="Survival Function")
  curve(1- pweibull(x, scale = sc_vc, shape = sh_vc),0,20,col=vc_col,add=TRUE)
  abline(h=0,v=0,col="grey",lty=2)
  legend("topright",bty="n",legend=c("BO",
                                     paste("Scale",round(sc_bo,2)),
                                     paste("Shape",round(sh_bo,2)),
                                     "VC",
                                     paste("Scale",round(sc_vc,2)),
                                     paste("Shape",round(sh_vc,2))),
         col=c(bo_col,NA,NA,vc_col,NA,NA),lty=c(1,NA,NA,1,NA,NA))
  
  
  lines(non_para_in$BO$BaseHaze_Zero$time, non_para_in$BO$BaseHaze_Zero$S.base,
        col=bo_col, lty=2, type = "s")
  lines(non_para_in$VC$BaseHaze_Zero$time, non_para_in$VC$BaseHaze_Zero$S.base,
        col=vc_col, lty=2, type = "s")
  
  
  if(eps){ dev.off() }
}
plot_cox.weibull()



S_WBcox <- function(from_Date, to_Date,
                     beta_MSCI, beta_HYS, beta_FundAge,
                     scale_wb, shape_wb,
                     FundAge,
                     public_data){
  from_Date <- as.Date(from_Date)
  if(is.na(to_Date)){
    to_Date <- as.Date("2016-12-31")
    Censoring <- TRUE
  }else{
    to_Date <- as.Date(to_Date)
    Censoring <- FALSE
  }
  public_data <- public_data[public_data$Date >= from_Date & public_data$Date <= to_Date,]
  public_data$t <- as.numeric(public_data$Date - public_data$Date[1])/365.25
  public_data$FundAge <- FundAge
  
  public_data$base_haze_rate_wb_exact <- c(0,diff( (public_data$t/scale_wb)^shape_wb ))
  
  public_data$expBX <- exp( as.matrix(public_data[,c("MSCI_monthly_return","ML_HYOAS","FundAge")]) %*% 
                              c(beta_MSCI,beta_HYS,beta_FundAge) )
  
  public_data$Cum_Haze_exact <- cumsum(as.numeric(public_data$base_haze_rate_wb_exact) * public_data$expBX)
  public_data$Surv_WB <- exp(-public_data$Cum_Haze_exact)
  
  if(Censoring){
    out <- public_data$Surv_WB[nrow(public_data)]  * runif(1) # basically NA
  }else{
    out <- public_data$Surv_WB[nrow(public_data)]
  }
  return(out)
}
Rosenblatt_WBcox <- function(FundType){
  
  data <- data_LoLi[data_LoLi$Type == FundType,]
  x <- Para_Cox_Weibull[[FundType]]$par
  
  outlist <- list()
  for(i in 1:nrow(data)){
    S <- S_WBcox(from_Date = data$InvDate[i],
                  to_Date = data$ExitDate[i],
                  x[1], x[2], x[3], x[4], x[5],
                  FundAge = max(0, data$YearInvest[i] - data$Vintage[i]),
                  public_data = public_data)
    
    outlist[i] <- S
  }
  return(as.numeric(outlist))
}
Weibull_Roseblatt <- function(eps=FALSE){
  if(eps){
    setEPS() ; postscript("Timing_Rosenblatt.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=5)
  }
  bo_col <- "royalblue1" ; vc_col <- "maroon1"
  
  par(mar=c(3,3,2,1),mfrow=c(1,2),oma=c(2.5,2,1,1))
  for(type in c("BO","VC")){
    robla <- Rosenblatt_WBcox(type)
    
    hi_col <- ifelse(type=="BO",bo_col,vc_col)
    hist(robla,freq = FALSE,
         main=type,border=hi_col,lty=3,ylab=NA,xlab=NA,ylim=c(0,1.5),xlim=c(0,1))
    abline(a=0,b=1,col="darkgray",lwd=1)
    wb_ecdf <- ecdf(robla)
    curve(wb_ecdf,add=TRUE,col="darkgreen",lwd=1)
    legend("topright",bty="n",cex=1.3,col=c("darkgreen","darkgray"),legend = c("empirical","theoretical"),lty=1)
    abline(h=1,col="black",lty=2)
  }
  mtext(latex2exp::TeX('$\\S_{Rosenblatt}(t_i,T_i)'),side=1,outer= TRUE,line=1,cex=1.7)
  mtext("ECDF or Density",side=2,outer= TRUE,line=0.5,cex=1.5)
  
  
  if(eps){ dev.off() }
}
Weibull_Roseblatt()



}
## c) Multiple ---------
if(TRUE){ # Data Preparation for Multiple (250/250 subset)
  G.p2c <- G.p2c[!is.na(G.p2c$P2C.multi1),] # JUST EXITED INVESTMENTS
  G.p2c  <- G.p2c[G.p2c$Fund_InvestTypes %in% c("BO","VC") & G.p2c$Fund_Region %in% c("US","EU","Asia"),]
  regression_variables <- c("P2C.multi1","RVPI_1","Fund_Region","GICS_Sector","Holding_Period","Time2Exit","Time2Exit_ZS","ZombieStage","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
  G.p2c <- G.p2c[complete.cases(G.p2c[,regression_variables]),]
  
  # delete young entries (bias correction?)
  G.p2c <- G.p2c[G.p2c$Investment_Date < as.Date("2010-01-01"),]
  G.p2c <- G.p2c[order(G.p2c$Fund_Emi_ID,G.p2c$FundInv_Quarter),]
  sun(G.p2c$P2C.multi1,main="Multiple Regression: Dependent Variable")
  
  # craeating subset
  bo_ss <- G.p2c[G.p2c$Fund_InvestTypes =="BO",]
  if(TRUE){
    set.seed(99)
    bo_ss <- bo_ss[bo_ss$Fund_Emi_ID %in% sample(levels(as.factor(bo_ss$Fund_Emi_ID)),250),]
    if(nrow(bo_ss) != 4729) rm(bo_ss)
  }
  cor(bo_ss[,c("P2C.multi1","ZombieStage","Time2Exit","Holding_Period","RVPI","ML_HYOAS.quarter","MSCI.Multiple.Exit")],method="pearson")
  id_bo <- unique(bo_ss$Fund_Emi_ID)
  
  vc_ss <- G.p2c[G.p2c$Fund_InvestTypes =="VC",]
  if(TRUE){
    set.seed(99)
    vc_ss <- vc_ss[vc_ss$Fund_Emi_ID %in% sample(levels(as.factor(vc_ss$Fund_Emi_ID)),250),]
    if(nrow(vc_ss) != 4421) rm(vc_ss)
  }
  cor(vc_ss[,c("P2C.multi1","ZombieStage","Time2Exit","Holding_Period","RVPI","ML_HYOAS.quarter","MSCI.Multiple.Exit")],method="pearson")
  id_vc <- unique(vc_ss$Fund_Emi_ID)
  
  length(id_bo) + length(id_vc)
  
  sun(bo_ss$P2C.multi1[bo_ss$FundInv_Quarter == bo_ss$Investment_Date])
  sun(vc_ss$P2C.multi1[vc_ss$FundInv_Quarter == vc_ss$Investment_Date])
}
if(TRUE){
  # Create data.frame for test procedure
  creat_reg_df <- function(df= G.p2c,ID_BO= id_bo,ID_VC= id_vc,reg_vars= regression_variables,stochastic=TRUE){
    reg_vars2 <- reg_vars[!(reg_vars %in% c("Time2Exit","ZombieStage","GICS_Sector","Fund_Region"))]
    # df2 <- df[df$Fund_Emi_ID %in% c(ID_BO,ID_VC),]
    df2 <- df[df$Investment_Date < as.Date("2010-01-01"),]
    df2 <- df2[df2$Fund_InvestTypes %in% c("BO","VC"),c(reg_vars2,
                                                        "Fund_Emi_ID","Fund_InvestTypes",
                                                        "FundInv_Quarter","Investment_Date")]
    
    if(stochastic){
      dt <- data.table::as.data.table(df2)
      dt2 <- dt[,.SD[sample(.N,1)], by=Fund_Emi_ID]
      df2 <- as.data.frame(dt2)
      # df2 <- as.data.frame(df2 %>% group_by(Fund_Emi_ID) %>% sample_n(1, replace=TRUE))
    }else{
      df2 <- as.data.frame(df2 %>% group_by(Fund_Emi_ID) %>% filter(FundInv_Quarter == Investment_Date))
    }
    DF <- list()
    DF[["BO"]] <- df2[df2$Fund_InvestTypes == "BO",reg_vars2]
    DF[["VC"]] <- df2[df2$Fund_InvestTypes == "VC",reg_vars2]
    return(DF)
  }
  system.time( df2 <- creat_reg_df(stochastic=TRUE) )
  nrow(df2$VC) ; nrow(df2$BO)
  # View(df2$VC)
  # binomial-gamma hurdle model
  # http://seananderson.ca/2014/05/18/gamma-hurdle.html
  # https://www.researchgate.net/post/Are_there_generalizations_of_zero-inflated_negative_binomial_and_hurdle_modeling_that_address_continuous_ie_non-count_variables
  # https://www.ibm.com/developerworks/library/ba-optimR-john-nash/
  
  iterativ_bbgh <- function(n=10,type= "BO", METH= "ucminf", benchmark_optimx=FALSE){
    if(n == 1){
      predictors_1 <- c("Time2Exit_ZS","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
      predictors_0 <- c("Time2Exit_ZS","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
      predictors_2 <- c("ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
    }else{
      # predictors_1 <- c("RVPI_1","Holding_Period","Time2Exit_ZS","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
      # predictors_1 <- c("Time2Exit_ZS","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
      predictors_2 <- c("RVPI_1","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
      predictors_1 <- c("Holding_Period","Time2Exit_ZS","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
      predictors_0 <- c("RVPI_1","Holding_Period","Time2Exit_ZS","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
    }
    response1 <- "P2C.multi1"
    hurdle0 <- 0.1
    hurdle2 <- 5
    
    # Function generates start values for Likelihood Maximation
    bbgammahu <- function(df){
      df$hurdle0 <- ifelse(df[,response1] > hurdle0, 1, 0)
      df$hurdle2 <- ifelse(df[,response1] > hurdle2, 1, 0)
      
      formula0 <- as.formula(paste("hurdle0", paste(predictors_0, collapse=" + "), sep=" ~ "))
      formula0_IO <- as.formula(paste("hurdle0", 1, sep=" ~ "))
      formula1 <- as.formula(paste(response1, paste(predictors_1, collapse=" + "), sep=" ~ "))
      formula1_IO <- as.formula(paste(response1, 1, sep=" ~ "))
      formula2 <- as.formula(paste("hurdle2", paste(predictors_2, collapse=" + "), sep=" ~ "))
      formula2_IO <- as.formula(paste("hurdle2", 1, sep=" ~ "))
      
      glms <- list()
      glms$m0 <- glm(formula0, data = df, family = binomial(link = logit))
      glms$m1 <- glm(formula1, data = subset(df, hurdle0 == 1 & hurdle2 == 0), family = Gamma(link = log))
      # glms$m1 <- glm(formula1, data = subset(df, hurdle0 == 1), family = Gamma(link = log))
      glms$m2 <- glm(formula2, data = df, family = binomial(link = logit))
      glms$m2_IO <- glm(formula2_IO, data = df, family = binomial(link = logit))
      glms$m1_IO <- glm(formula1_IO, data = subset(df, hurdle0 == 1 & hurdle2 == 0), family = Gamma(link = log))
      # glms$m1_IO <- glm(formula1_IO, data = subset(df, hurdle0 == 1), family = Gamma(link = log))
      glms$m0_IO <- glm(formula0_IO, data = df, family = binomial(link = logit))
      
      out <- list()
      for(modelz in names(glms)[1:3]){
        m_IO <- glms[[paste(modelz,"IO",sep = "_")]]
        m <- glms[[modelz]]
        out[[modelz]] <- c(Intercept_IO= as.numeric(coef(m_IO)), 
                           Theta_IO= summary(m_IO)$dispersion, 
                           AIC_IO= summary(m_IO)$aic,
                           coef(m), 
                           Theta= summary(m)$dispersion,
                           AIC= summary(m)$aic)
      }
      
      return(out)
    }
    
    out <- list() ; out_nb2 <- list()
    for(i in seq(1,n)){
      print(i)
      # create right input data
      if(n == 1){ stoch_df = FALSE }else{ stoch_df = TRUE }
      df_bovc <- creat_reg_df(stochastic = stoch_df) # create random input
      df_bovc <- df_bovc[[type]]
      # if(n == 1){ df_bovc = df_bovc[,!(colnames(df_bovc) %in% c("RVPI_1","Holding_Period"))] }

      
# PART 1 - Negative Binomial GLM
      count_depvar <- paste("I(round(",response1," * ",100,",0))",sep="")
      predictors_nb2 <- c("ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
      formula_NB2 <- as.formula(paste(count_depvar, paste(predictors_0, collapse=" + "), sep=" ~ "))
      
      error_status_glm.nb <- try( nb2_reg <- MASS::glm.nb(formula_NB2,data=df_bovc) )
      if(class(error_status_glm.nb) == "try-error"){
        print(paste("NB2: Error occured: Run ",i))
        create_output <- FALSE
        # break
      }else{
        if(nb2_reg$theta > 100){
          print(paste("NB2: Theta too high: Run ",i))
          create_output <- FALSE
        }else{
          out_nb2[[paste("r",i,sep="")]] <- c(nb2_reg$coefficients, Theta= nb2_reg$theta, AIC= nb2_reg$aic)
          create_output <- TRUE
        }
      }
      
      
# PART 2 - Hurdle Model
      # m_1
      if(create_output){
        out_123 <- bbgammahu(df_bovc)
        # truncdist::dtrunc
        trunc_dist <- function(x, spec="gamma", a = hurdle0, b = hurdle2, ...){
          if (a >= b) 
            stop("argument a is greater than or equal to b")
          tt <- rep(0, length(x))
          g <- get(paste("d", spec, sep = ""), mode = "function")
          G <- get(paste("p", spec, sep = ""), mode = "function")
          G.a <- G(a, ...)
          G.b <- G(b, ...)
          '
          if (G.a == G.b) {
            print(data.frame(G_a= G.a, G_b= G.b))
            stop(">> Trunction interval is not inside the domain of the density function")
          }
          '
          tt[x >= a & x <= b] <- g(x[x >= a & x <= b], ...)/(G(b, ...) - 
                                                               G(a, ...))
          return(tt)
        }
        # define likelihood functions
        ml_gamma_trunc <- function(par, df= df_bovc, 
                                   response= response1, predictors= predictors_1, 
                                   low= hurdle0, high= hurdle2){
          df1 <- df[df[,response] > hurdle0 & df[,response] < hurdle2, ]
          df1 <- df1[complete.cases(df1),]
          df1$Intercept <- 1
          Y <- df1[,response]
          X <- as.matrix(df1[,c("Intercept", predictors)])
          
          shape1 <- par[length(par)]
          b <- par[-length(par)]
          
          Xb <- X %*% b
          mu <- exp(Xb)
          scale1 <- mu / shape1
          
          dtg <- trunc_dist(Y, shape=shape1,scale=scale1)
          #LikFun <- prod(dtg)
          #return(-LikFun)
          LogLikFun <- sum(log(dtg))
          return(-LogLikFun)
        }
        ml_gamma_trunc_IO <- function(par, df= df_bovc, 
                                      response= response1, predictors= predictors_1, 
                                      low= hurdle0, high= hurdle2){
          df1 <- df[df[,response] > hurdle0 & df[,response] < hurdle2, ]
          df1 <- df1[complete.cases(df1),]
          df1$Intercept <- 1
          Y <- df1[,response]
          X <- as.matrix(df1[,c("Intercept")])
          
          shape <- par[length(par)]
          b <- par[-length(par)]
          
          Xb <- X %*% b
          mu <- exp(Xb)
          scale <- mu / shape
          
          dtg <- trunc_dist(Y, spec="gamma",shape=shape,scale=scale,a=low,b=high)
          #LikFun <- prod(dtg)
          #return(-LikFun)
          LogLikFun <- sum(log(dtg))
          return(-LogLikFun)
        }
        # Optimize m1
        for(j in seq(1,1)){
          start <- out_123$m1[!(names(out_123$m1) %in% c("Intercept_IO","Theta_IO","AIC_IO","AIC"))]
          m1_result <- optimx::optimx(start, ml_gamma_trunc, method = METH)
          # m1_result <- Rsolnp::solnp(start, ml_gamma_trunc)
          
          # compare solnp_result to other opimization methods
          if(benchmark_optimx){ 
            optim_result <- optimx::optimx(start, ml_gamma_trunc,
                                           control = list(all.methods=TRUE, save.failures=TRUE))
            print(optim_result)
          }
          
          # report non-convergence
          if(m1_result$convcode != 0){
            print(paste("*** !!!! --- No convergence m0: Run", i))
          }

          # store optimzation results
          end_par <- unlist(m1_result[1:(length(predictors_1) + 2)])
          out_123$m1[!(names(out_123$m1) %in% c("Intercept_IO","Theta_IO","AIC_IO","AIC"))] <- end_par
          out_123$m1["AIC"] <- 2 * (length(predictors_1) + 2) - 2 * (-m1_result$value)

          
          # intercept only model update
          start_IO <- out_123$m1[c("Intercept_IO","Theta_IO")]
          m1_result_IO <- optimx::optimx(start_IO,ml_gamma_trunc_IO, method = METH)
          # m1_result_IO <- Rsolnp::solnp(start_IO, ml_gamma_trunc_IO)
          
          end_par_IO <- unlist(m1_result_IO[1:2])
          out_123$m1[c("Intercept_IO","Theta_IO")] <- end_par_IO
          out_123$m1["AIC_IO"] <- 2 * 2 - 2 * (-m1_result_IO$value)
        }
        
      }
      # m_0 and m_2
      if(create_output & is.finite(hurdle2)){
        # likelihood function full parameter model
        ml_hurdles02 <- function(par, df= df_bovc, 
                                 response= response1, 
                                 predictors_lower = predictors_0, 
                                 predictors_upper = predictors_2, 
                                 low= hurdle0, high= hurdle2){
          
          df$hurdle0 <- ifelse(df[,response1] > hurdle0, 1, 0)
          df$hurdle2 <- ifelse(df[,response1] > hurdle2, 1, 0)
          df <- df[complete.cases(df),]
          df$Intercept <- 1
          
          Y <- df[,response]
          
          X_0 <- as.matrix(df[,c("Intercept", predictors_lower)])
          b_0 <- par[1 : (1 + length(predictors_lower))]
          X_2 <- as.matrix(df[,c("Intercept", predictors_upper)])
          b_2 <- par[(length(predictors_lower) + 1) + (1 : (1 + length(predictors_upper)))]
          
          # cf. Barnett Dobson "8.4.3 Adjacent categories logit model", p. 160
          
          f_0 <- 1 / (1 + exp(-(X_0 %*% b_0))) # P(Y > h_0)
          f_0 <- 1 - f_0                       # P(Y <= h_0)
          f_2 <- 1 / (1 + exp(-(X_2 %*% b_2))) # P(Y > h_2)
          
          likelihood <- ifelse(Y <= low, f_0,
                               ifelse(Y > high, f_2,
                                      1 - f_0 - f_2))
          
          # check 1 - f_0 - f_2 >= 0
          if(min(likelihood) < 0){
            print("Issue m0 & m2:  1 - f_0 - f_2 >= 0")
            likelihood <- pmax(0.0001,likelihood)
          }
          
          return(-sum(log(likelihood)))
        }
        ml_hurdles02_IO <- function(par, df= df_bovc, 
                                 response= response1, 
                                 low= hurdle0, high= hurdle2){
          
          df$hurdle0 <- ifelse(df[,response1] > hurdle0, 1, 0)
          df$hurdle2 <- ifelse(df[,response1] > hurdle2, 1, 0)
          df <- df[complete.cases(df),]
          df$Intercept <- 1
          
          Y <- df[,response]
          
          X_0 <- as.matrix(df[,c("Intercept")])
          b_0 <- par[1]
          X_2 <- as.matrix(df[,c("Intercept")])
          b_2 <- par[2]
          
          # cf. Barnett Dobson "8.4.3 Adjacent categories logit model", p. 160
          
          f_0 <- 1 / (1 + exp(-(X_0 %*% b_0))) # P(Y > h_0)
          f_0 <- 1 - f_0                       # P(Y <= h_0)
          f_2 <- 1 / (1 + exp(-(X_2 %*% b_2))) # P(Y > h_2)
          
          likelihood <- ifelse(Y <= low, f_0,
                               ifelse(Y > high, f_2,
                                      1 - f_0 - f_2))
          
          # check 1 - f_0 - f_2 >= 0
          if(min(likelihood) < 0){
            print("Issue m0 & m2:  1 - f_0 - f_2 >= 0")
            likelihood <- pmax(0.0001,likelihood)
          }
          
          return(-sum(log(likelihood)))
        }
        
        # Optimize m0 and m2
        for(j in seq(1,1)){
          # A) Full Parameter Model
          par_hurdles02 <- c(out_123$m0[!(names(out_123$m0) %in% c("Intercept_IO","Theta_IO","AIC_IO","Theta","AIC"))],
                             out_123$m2[!(names(out_123$m2) %in% c("Intercept_IO","Theta_IO","AIC_IO","Theta","AIC"))])
          
          m02_result <- optimx::optimx(par_hurdles02, ml_hurdles02, method = METH)

          # compare solnp_result to other opimization methods
          if(benchmark_optimx){ 
            optim_result_02 <- optimx::optimx(par_hurdles02, ml_hurdles02,
                                           control = list(all.methods=TRUE, save.failures=TRUE))
            print(optim_result_02)
          }
          
          # report non-convergence
          if(m02_result$convcode != 0){
            print(paste("*** !!!! --- No convergence m0 & m2: Run", i))
          }
          
          # store optimzation results
          end_par_m0 <- unlist(m02_result[1:(1 + length(predictors_0))])
          out_123$m0[!(names(out_123$m0) %in% c("Intercept_IO","Theta_IO","AIC_IO","Theta","AIC"))] <- end_par_m0
          
          end_par_m2 <- unlist(m02_result[(1 + length(predictors_0)) + (1 : (1 + length(predictors_2)))])
          out_123$m2[!(names(out_123$m2) %in% c("Intercept_IO","Theta_IO","AIC_IO","Theta","AIC"))] <- end_par_m2

          AIC_m02 <- 2 * length(par_hurdles02) - 2 * (-m02_result$value)
          out_123$m2["AIC"] <- out_123$m0["AIC"] <- AIC_m02
          out_123$m2["Theta"] <- out_123$m0["Theta"] <- NA
          
          # B) Intercept Only Model
          par_hurdles02_IO <- c(out_123$m0["Intercept_IO"],
                                out_123$m2["Intercept_IO"])
          
          m02_result_IO <- optimx::optimx(par_hurdles02_IO, ml_hurdles02_IO, method = METH)
          
          # report non-convergence
          if(m02_result_IO$convcode != 0){
            print(paste("*** !!!! --- No convergence m0 & m2 IO: Run", i))
          }
          
          # store optimzation results
          end_par_m0_IO <- unlist(m02_result_IO[1])
          out_123$m0["Intercept_IO"] <- end_par_m0_IO
          
          end_par_m2_IO <- unlist(m02_result_IO[2])
          out_123$m2["Intercept_IO"] <- end_par_m2_IO
          
          AIC_m02_IO <- 2 * length(par_hurdles02_IO) - 2 * (-m02_result_IO$value)
          out_123$m2["AIC_IO"] <- out_123$m0["AIC_IO"] <- AIC_m02_IO
          out_123$m2["Theta_IO"] <- out_123$m0["Theta_IO"] <- NA
        }
        
      }
      # store result
      if(create_output){
      out[[paste("r",i,sep = "")]] <- out_123
      }
    }
    
    # Prep Output: NB2
    out_nb2 <- data.frame(do.call(rbind,out_nb2))
    
    
    # Prep Output: Hurdle Model (h1)
    out2 <- do.call(abind::abind, c(
      lapply(out, function(l){
        l2 <- as.data.frame(do.call(rbind, lapply(lapply(l, unlist), "[",
                                                  unique(unlist(c(sapply(l,names)))))))
        names(l2) <- unique(unlist(c(sapply(l,names))))
        l2
      })
      , along=3))
    
    return(list(NB2= out_nb2, HM=out2))
  }
  prep_MultiReg_output <- function(input){
    in_hm <- input$HM
    in_nb <- input$NB2
    
    out_hm <- list() ; out_nb <- list()
    for(statics in c("mean","median","sd")){
      out_nb[[statics]] <- data.frame(apply(in_nb, 2, statics))
      out_hm[[statics]] <- data.frame(apply(in_hm,  1:2, statics))
    }
    
    prepare_out <- function(out){
      out <- data.frame(do.call(cbind,out))
      colnames(out) <- c("Mean","Median","SD")
      out$t_value <- out$Mean / out$SD
      out$p_value <- pnorm(-abs(out$t_value))
      return(round(out,3))
    }
    
    df_M <- data.frame(cbind(t(out_hm$mean), 
                             pnorm(-abs( t(out_hm$mean) / t(out_hm$sd) ))
                             ))
    colnames(df_M)[4:6] <- paste("p_value",0:2,sep="")
    df_M <- round(df_M,3)

    OUTPUT <- list(NB= prepare_out(out_nb), HM= df_M)
    return((OUTPUT))
  }
  
  set.seed(98) ; system.time(out_bbgh_bo <- iterativ_bbgh(3,"BO"))
  set.seed(99) ; system.time(out_bbgh_bo_1 <- iterativ_bbgh(100,"BO"))
  prep_MultiReg_output(out_bbgh_bo)
  prep_MultiReg_output(out_bbgh_bo_1)
  
  set.seed(98) ; system.time(out_bbgh_vc <- iterativ_bbgh(500,"VC"))
  set.seed(99) ; system.time(out_bbgh_vc_1 <- iterativ_bbgh(100,"VC"))
  prep_MultiReg_output(out_bbgh_vc)
  prep_MultiReg_output(out_bbgh_vc_1)
  
  
  # out_bbgh <- out_bbgh_bo
  # df_M <- data.frame(cbind(t(out_bbgh$mean), t(out_bbgh$mean) / t(out_bbgh$sd)))
  # colnames(df_M)[4:6] <- paste("z",1:3,sep="")
  # write.csv(round(df_M,3),"Hurdle_Gamma_GLM.csv")
  
  predict_dgh <- function(df_in= bo_vc, model_in= list(BO=out_bbgh_bo, VC=out_bbgh_vc),
                          plot_it=FALSE, eps=FALSE,
                          sim_out= FALSE, sim_in= bo_vc[sample(seq(nrow(bo_vc)),1),], runo=runif(nrow(bo_vc))){
    mobo <- data.frame(t(model_in$BO$mean))
    movc <- data.frame(t(model_in$VC$mean))
    h0 <- 0.3 ; h2 <- 3
    
    ecdf_bo_h0 <- ecdf(df_in$P2C.multi1[df_in$Fund_InvestTypes == "BO" & df_in$P2C.multi1 < h0])
    ecdf_bo_h2 <- ecdf(df_in$P2C.multi1[df_in$Fund_InvestTypes == "BO" & df_in$P2C.multi1 > h2])
    ecdf_vc_h0 <- ecdf(df_in$P2C.multi1[df_in$Fund_InvestTypes == "VC" & df_in$P2C.multi1 < h0])
    ecdf_vc_h2 <- ecdf(df_in$P2C.multi1[df_in$Fund_InvestTypes == "VC" & df_in$P2C.multi1 > h2])
    
    output <- list()
    df_in <- if(sim_out)  sim_in else df_in
    iter <-   seq(nrow(df_in))
    for(i in iter){
      real_multiple <- df_in$P2C.multi1[i]
      
      type= as.character(df_in$Fund_InvestTypes[i])
      Xi <- c(Intercept= 1,
              RVPI_1= df_in$RVPI_1[i],
              HP= df_in$Holding_Period[i],
              T2E= df_in$Time2Exit_ZS[i],
              HYS= df_in$ML_HYOAS.quarter[i],
              PEM_1= df_in$MSCI.Multiple.Exit_1[i])
      b0 <- if(type=="BO") mobo[4:9,"m0"] else movc[4:9,"m0"]
      b1 <- if(type=="BO") mobo[4:9,"m1"] else movc[4:9,"m1"]
      b2 <- if(type=="BO") mobo[4:9,"m2"] else c(movc[1,"m2"],rep(0,5))
      
      p0 <- 1 - plogis(sum(Xi * b0)) # P(Y < hurdle0)
      p2 <- plogis(sum(Xi * b2))    # P(Y > hurdle2)
      p1 <- max(0.00999, 1 - p0 - p2)    # P(Y <= hurdle2 & Y >= hurdle0)
      mu1 <- exp(sum(Xi * b1))
      shape1 <- ifelse(type=="BO",mobo["Theta","m1"],movc["Theta","m1"])
      scale1 <- mu1/shape1
      
      if(plot_it){
        print(Xi)
        P <- c(P0= p0,P1= p1,P2= p2)
        print(P)
      }
      compound_CDF <- function(Y,tyype){
        ecdf_h0 <- if(tyype == "BO") ecdf_bo_h0 else ecdf_vc_h0
        ecdf_h2 <- if(tyype == "BO") ecdf_bo_h2 else ecdf_vc_h2
        out <- list()
        for(i in 1:length(Y)){
          value <- NA
          if(Y[i] < h0) value <- ecdf_h0(Y[i]) * p0
          if(Y[i] > h2) value <- ecdf_h2(Y[i]) * p2 + (p0 + p1)
          if((Y[i] <= h2) & (Y[i] >= h0)){value <- truncdist::ptrunc(Y[i], spec="gamma",shape=shape1,scale=scale1,a=h0,b=h2) * p1 + p0} 
          out[i] <- value
        }
        return(unlist(out))
      }
      inverse_CDF <- function(prob){
        inverse <- function(f, lower = -0.0001, upper = 100){
          function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
        }
        inv_fun <- inverse(function(x){compound_CDF(x,type)})
        out <- as.numeric(inv_fun(prob))
        out <- max(0,out)
        return(out)
      }
      
      if(sim_out){
        output[i] <-  inverse_CDF(runo[i])
      }else{
        output[i] <- compound_CDF(real_multiple,type)
      }
      
      if(plot_it){
        if(eps){
          setEPS() ; postscript("Compound_CDF_DGH.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=5)
        }
        # BO
        par(cex=1.3)
        curve(ecdf_bo_h0,0,h0,xlim=c(0,10),ylim=c(0,1.1),ylab="CDF",xlab="Multiple (Y)",col="red")
        abline(v=c(h0,h2),col="gray",lty=2) ; abline(h=c(0,1),col="gray")
        curve(ecdf_bo_h2,add=TRUE,col="green",h2,10)
        curve(truncdist::ptrunc(x, spec="gamma",shape=shape1,scale=scale1,a=h0,b=h2),h0,h2,col="blue",add=TRUE)
        
        seq_0 <- seq(0,h0,0.01) ; seq_1 <- seq(h0,h2,0.01) ; seq_2 <- seq(h2,10,0.01)
        points(seq_0,compound_CDF(seq_0,type),col="red",type="l",lwd=2)
        points(seq_1,compound_CDF(seq_1,type),col="blue",type="l",lwd=2)
        points(seq_2,compound_CDF(seq_2,type),col="green",type="l",lwd=2)
        legend("bottomright",bty="n",col=c("green","blue","red"),lty=1,cex=1.2,lwd=2,
               legend=c(latex2exp::TeX('$\\Y > h_{2} $'),
                        latex2exp::TeX('$\\h_{0} \\leq Y \\leq h_{2} $'),
                        latex2exp::TeX('$\\Y < h_{0} $')
               ))
        legend(x=6,y=0.5,bty="n",legend=paste(names(Xi[-1]),round(Xi[-1],2)),cex=0.6,text.col="darkgrey")
        text(x=0,y=1.08,round(p0,3),cex=0.8,col= "red")
        text(x=1.7,y=1.08,round(p1,3),cex=0.8, col= "blue")
        text(x=4,y=1.08,round(p2,3),cex=0.8, col="green")
        
        if(eps){ dev.off() }
      }
      
    }
    output <- as.numeric(unlist(output))
    return(output)
  }
  # bo_vc$RQ_Multiple.dgh <- predict_dgh()
  predict_dgh(plot_it = TRUE,sim_out=TRUE,eps=FALSE)
  
  DGH_Roseblatt <- function(eps=FALSE){
    if(eps){
      setEPS() ; postscript("Hurdle_Rosenblatt.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=5)
    }
    bo_col <- "royalblue1" ; vc_col <- "maroon1"
    
    par(mar=c(3,3,2,1),mfrow=c(1,2),oma=c(3,2,1,1),cex=1.2)
    for(type in c("BO","VC")){
      hi_col <- ifelse(type=="BO",bo_col,vc_col)
      hist(bo_vc$RQ_Multiple.dgh[bo_vc$Fund_InvestTypes == type],freq = FALSE,
           main=type,border=hi_col,lty=3,ylab=NA,xlab=NA,ylim=c(0,1.7),xlim=c(0,1))
      abline(a=0,b=1,col="darkgray",lwd=1,lty=2)
      wb_ecdf <- ecdf(bo_vc$RQ_Multiple.dgh[bo_vc$Fund_InvestTypes == type])
      curve(wb_ecdf,add=TRUE,col="forestgreen",lwd=1)
      legend("topright",bty="n",cex=1.3,col=c("forestgreen","darkgray"),legend = c("empirical","theoretical"),lty=c(1,2))
      abline(h=1,col="black",lty=2)
    }
    mtext(latex2exp::TeX('$\\F_{0,1,2}^{(DGH)}(Y_i)}'),side=1,outer= TRUE,line=1,cex=1.7)
    mtext("ECDF or Density",side=2,outer= TRUE,line=0.5,cex=1.5)
    
    
    if(eps){ dev.off() }
  }
  DGH_Roseblatt()
  
  
  if(FALSE){
    shape <- theta <- df_M["Theta_IO","m1"]
    mu <- exp(df_M["X.Intercept.","m1"])
    scale <- mu / shape
    
    seq_x <- seq(0,8,0.01) ; bo_col <- "royalblue1"
    
    par(par_default)
    # setEPS() ; postscript("Truncated_Gamma_PDF.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=5)
    par(cex=1.3)
    plot(seq_x,truncdist::ptrunc(seq_x, spec="gamma",shape=shape,scale=scale,a=0.3,b=3),lwd=1,
         col= bo_col, type="l",ylab="CDF / PDF",xlab="Multiple")
    lines(seq_x,truncdist::ptrunc(seq_x, spec="gamma",shape=shape,scale=scale,a=0.3,b=100),lty=3,col=bo_col,lwd=2)
    abline(v=c(0.3,3),col="gray",lty=2)
    lines(seq_x,truncdist::dtrunc(seq_x, spec="gamma",shape=shape,scale=scale,a=0.3,b=3),lwd=1, col= bo_col, type="l",ylab="PDF",xlab="Multiple")
    lines(seq_x,truncdist::dtrunc(seq_x, spec="gamma",shape=shape,scale=scale,a=0,b=100),lty=3,col=bo_col,lwd=2)
    legend("right",bty="n",legend=c("truncated gamma","gamma"),lty=c(1,3),col=c(bo_col),lwd=c(1,2))
    # dev.off()
  }
}

if(FALSE){
  # rbind to summary function
  rbind2summary <- function(df_rbind){
    df_rbind <- df_rbind[complete.cases(df_rbind),]
    
    df_rbind <- data.frame(t(do.call(rbind,list(Mean= colMeans(df_rbind), 
                                                Median= apply(df_rbind, 2, median),
                                                P_negative= apply(df_rbind,2,function(x) sum(x<=0)/length(x)), 
                                                SD= apply(df_rbind, 2, sd)))))
    df_rbind[,"Mean/SD"] <- df_rbind$Mean / df_rbind$SD
    return(df_rbind)
  }

  # Procedure to test several regression approches (iter-fold-cross-validation)
  moa_cv <- function(iterations=10,cv_fold=10, max_y=10,reg_vars= regression_variables,do_legend= TRUE,eps=FALSE){
    # ia_depth= 3 ; n_minObsNode=30 ; shrnk = 0.001 ; n_tree_multi= 10 # bgm
    # method_kknn= "mahalanobis" ; weights_kknn= "triangular" ; k_percent= 75    # knn
    
    # Plot Basis
    if(eps){
      setEPS() ; postscript("Multiple CrossValidation.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=5)
    }
    par(mfrow=c(1,2),mar=c(2,2,2,1),oma=c(3,3,1,1))
    pred_list <- list() ; coef_list <- list()
    do_regressions <- TRUE
    cex_para <- 1.25
    
    for(iter in seq(iterations)){
      print(iter)
      df= creat_reg_df()

      for(type in c("BO","VC")){
        if(iter==1){
          colines <- "blue"
          plot(c(0,max_y),c(0,max_y),type="l",lyt=1,xlab=NA, ylab=NA,
               xlim=c(-max_y,max_y),ylim=c(-max_y*0.5,max_y*1.5),main=type,col=colines)
          lines(-c(0,max_y),c(0,max_y),col=colines)
          abline(v=0,h=0,col=colines,lty=2)
        }
        
        # subset input by type & create CV buckets
        df1 <- df[[type]]
        all_rows <- seq(1,nrow(df1),1)
        cv_buckets <- split(sample(all_rows), letters[seq(1,cv_fold)])
        cv_buckets <- split(sample(all_rows), seq(1,cv_fold))
        
        # Regression Formula
        if(TRUE){
          if(type == "BO"){
            reg_vars2 <- reg_vars[!(reg_vars %in% c("Holding_Period","RVPI_1","Time2Exit","ZombieStage","GICS_Sector","Fund_Region"))]
          }else{
            reg_vars2 <- reg_vars[!(reg_vars %in% c("ML_HYOAS.quarter","RVPI_1","Time2Exit","ZombieStage","GICS_Sector","Fund_Region"))]
          }
          reg_vars2 <- reg_vars[!(reg_vars %in% c("Time2Exit","ZombieStage","GICS_Sector","Fund_Region"))]
          
          linear_formula <- as.formula(paste(reg_vars2[1], paste(reg_vars2[-1], collapse=" + "), sep=" ~ "))
          cound_factor <- 100
          count_depvar <- paste("I(round(",reg_vars2[1]," * ",cound_factor,",0))",sep="")
          count_depvar_truncated <- paste("I(pmax(0,round(",reg_vars2[1]," * ",cound_factor,",0)-0))",sep="")
          count_formula <- as.formula(paste(count_depvar, paste(reg_vars2[-1], collapse=" + "), sep=" ~ "))
          count_formula_IO <- as.formula(paste(count_depvar, 1, sep=" ~ "))
          
          # zero inflated NB2
          znb_count_part <- paste(count_depvar_truncated, paste(reg_vars2[-1], collapse=" + "), sep=" ~ ")
          znb_zero_part <- paste(reg_vars2[-1], collapse=" + ")
          zero_inf_formula_truncated <- as.formula(paste(znb_count_part, znb_zero_part, sep = "|"))
          zero_inf_formula_IO_truncated <- as.formula(paste(count_depvar_truncated, 1, sep=" ~ "))
          
          # Log-Normal glm
          gln_depvar <- paste("I(pmax(",reg_vars2[1],", 0.01))",sep="")
          gln_formula <- as.formula(paste(gln_depvar, paste(reg_vars2[-1], collapse=" + "), sep=" ~ "))
          gln_formula_IO <- as.formula(paste(gln_depvar, 1, sep=" ~ "))
        }
        
        # Cross Validation Scheme
        for(a in seq(1,cv_fold)){
        # Split df1 in train & test set
          if(do_regressions){
            test_rows <- unlist(cv_buckets[a])
            train_rows <- all_rows[!(all_rows %in% test_rows)]
            df_train <- df1[train_rows,]
            df_test <- df1[test_rows,]
            x_train <- model.matrix( ~ .-1,data=df_train)
            x_test <- model.matrix( ~ .-1,data=df_test)
            df_test <- data.frame(x_test)
          }
        # Estimation of linear + GLM.nb + non-linear models
          if(do_regressions){
            # glm (gamma)
            gga_reg_IO <- glm(gln_formula_IO, family=Gamma(link = "log"), data=df_train)
            try_gamma <- try(gga_reg <- glm(gln_formula, family=Gamma(link = "log"), data=df_train))
            if(class(try_gamma) == "try-error"){
              print(paste("Gamma: Error occured at",iter,type,a))
              break
            }
            y_gga <- exp(predict(gga_reg,df_test,type="link"))
            
            # negative binomial regression
            # https://stats.idre.ucla.edu/r/dae/negative-binomial-regression/
            glm.nb_reg_IO <- MASS::glm.nb(count_formula_IO,data=df_train)
            error_status_glm.nb <- try(glm.nb_reg <- MASS::glm.nb(count_formula,data=df_train))
            if(class(error_status_glm.nb) == "try-error"){
              print(paste("NB2: Error occured at",iter,type,a))
              break
            }
            y_nb <- exp(predict(glm.nb_reg,df_test)) / cound_factor

            # zero-inflation
            # https://stats.idre.ucla.edu/r/dae/zip/
            zi.nb_reg_IO <- pscl::zeroinfl(zero_inf_formula_IO_truncated ,dist = "negbin",data=df_train)
            zi.nb_reg <- pscl::zeroinfl(zero_inf_formula_truncated ,dist = "negbin",data=df_train)
            y_zi.nb <- predict(zi.nb_reg,df_test) / cound_factor

            # hurdle
            hu.nb_reg_IO <- pscl::hurdle(count_formula_IO ,data=df_train,
                                      dist = "negbin", zero.dist= "binomial")
            hu.nb_reg <- pscl::hurdle(count_formula,data=df_train,
                                      dist = "negbin", zero.dist= "binomial")
            y_hu.nb <- predict(hu.nb_reg,df_test) / cound_factor




            
            # Tobit
            tobit_reg <- VGAM::vglm(linear_formula,VGAM::tobit(Lower = 0),data=df_train)
            y_tobit <- VGAM::predict(tobit_reg,df_test,type="response")
            # OLS
            ols_reg <- lm(linear_formula,data=df_train)
            y_ols <- predict(ols_reg,df_test)
            # robust
            rob_reg <- robustbase::lmrob(linear_formula,data=df_train)
            y_rob <- predict(rob_reg,df_test)


            

            
            '
            tobit_reg <- censReg::censReg(linear_formula,data=df_train,left=0)
            predict_tobit <- function(reg=tobit_reg,df_in=df_test){
            df_test2 <- df_in[,-1] # remove dependent variable
            df_test2$Intercept <- 1 # add intercept
            df_test2 <- df_test2[,c(ncol(df_test2),seq(1,ncol(df_test2)-1))] # intercept as 1. column
            prediction <- list()
            for(i in seq(nrow(df_test2))){
            prediction[i] <- sum(reg$estimate[1:5] * df_test2[i,])
            }
            prediction <- as.numeric(unlist(prediction))
            return(prediction)
            }
            y_tobit <- predict_tobit()
            coef_tobit <- coef(tobit_reg)[1:length(reg_vars2)]
            

            # glm (log-normal)
            gln_reg <- glm(gln_formula,data=df_train, family=gaussian(link="log"))
            y_gln <- exp(predict(gln_reg,df_test,type="link"))
            coef_gln <- coef(gln_reg)
            

            # Generalized Additive Model
            gam_reg <- gam::gam(count_formula, family = quasipoisson(link = "log"), data=df_train)
            y_gam <- predict(gam_reg,df_test)
            print(summary(gam_reg))


            # Boosting
            n_trees <-  n_tree_multi * nrow(df_train)
            gbm_reg <- gbm::gbm(linear_formula,
                                distribution="tdist",
                                n.trees= n_trees,
                                n.minobsinnode=n_minObsNode,shrinkage = shrnk,interaction.depth= ia_depth,
                                # weights= dgamma(df_train[,1],shape = 2,scale = 1),bag.fraction=1,
                                data=df_train)
            y_gbm <- predict(gbm_reg,df_test,n.trees=n_trees)
            
            # k-nearest means
            k_min <- round(nrow(df_train)/100)
            #knn.fit <- FNN::knn.reg(train= x_train,test= x_test,y= x_train[,1],k=k_min*i)
            y_knn <- KernelKnn::KernelKnn(data=x_train[,-1], 
                                          TEST_data= x_test[,-1],
                                          y= x_train[,1],
                                          k= k_min*k_percent,
                                          h=1, 
                                          method= method_kknn, 
                                          weights_function= weights_kknn,
                                          regression = TRUE)
            
            # SVM
            svm_reg <- e1071::svm(linear_formula,data=df_train)
            y_svm <- predict(svm_reg,df_test)
            
            # random forest
            rF_reg <- randomForest::randomForest(MOIC ~ HoPi + ZombieStage + ML_HYOAS + MSCI.Multiple_1,
            data=df_train,ntree= n_trees,maxnodes=10)
            y_rF <- predict(rF_reg,df_test)
            points(x_para,y_rF,pch=20,cex=cex_para*0.7,col="orange")
            '
          }
        # Write coefs to list
          if(do_regressions){
            # linear (OLS/ROB/TOB)
            coef_list[["linear"]][[paste("Type",type,sep="")]][[paste("CV",a,iter)]] <- data.frame(
              do.call(rbind,list(TOB= coef(tobit_reg)[-2],OLS= coef(ols_reg),ROB= coef(rob_reg))))
            
            # NB2
            coef_list[["NB2"]][["InterceptOnly"]][[paste("Type",type,sep="")]][[paste("CV",a,iter)]] <- data.frame(
              do.call(rbind,list(GNB= c(coef(glm.nb_reg_IO), 
                                        Theta= glm.nb_reg_IO$theta, 
                                        TwoLogLik= glm.nb_reg_IO$twologlik,
                                        AIC=AIC(glm.nb_reg_IO)))))
            coef_list[["NB2"]][["Regression"]][[paste("Type",type,sep="")]][[paste("CV",a,iter)]] <- data.frame(
              do.call(rbind,list(GNB= c(coef(glm.nb_reg),
                                        Theta= glm.nb_reg$theta, 
                                        TwoLogLik= glm.nb_reg$twologlik,
                                        AIC=AIC(glm.nb_reg)))))
            
            # ZNB2
            coef_list[["ZNB2"]][["InterceptOnly"]][[paste("Type",type,sep="")]][[paste("CV",a,iter)]] <- data.frame(
              do.call(rbind,list(ZNB= c(coef(zi.nb_reg_IO), 
                                        Theta= zi.nb_reg_IO$theta, 
                                        Log_Lik= zi.nb_reg_IO$loglik),
                                 AIC = AIC(zi.nb_reg_IO))))
            coef_list[["ZNB2"]][["Regression"]][[paste("Type",type,sep="")]][[paste("CV",a,iter)]] <- data.frame(
              do.call(rbind,list(ZNB= c(coef(zi.nb_reg),
                                        Theta= zi.nb_reg$theta,
                                        Log_Lik= zi.nb_reg$loglik, 
                                        AIC=AIC(zi.nb_reg)))))
            
            # HNB2
            coef_list[["HNB2"]][["InterceptOnly"]][[paste("Type",type,sep="")]][[paste("CV",a,iter)]] <- data.frame(
              do.call(rbind,list(HNB= c(coef(hu.nb_reg_IO),
                                        Theta= hu.nb_reg_IO$theta,
                                        Log_Lik= hu.nb_reg_IO$loglik, 
                                        AIC= AIC(hu.nb_reg_IO)))))
            coef_list[["HNB2"]][["Regression"]][[paste("Type",type,sep="")]][[paste("CV",a,iter)]] <- data.frame(
              do.call(rbind,list(HNB= c(coef(hu.nb_reg),
                                        Theta= hu.nb_reg$theta,
                                        Log_Lik= hu.nb_reg$loglik, 
                                        AIC= AIC(hu.nb_reg)))))
            
            # Gamma
            coef_list[["Gamma"]][["InterceptOnly"]][[paste("Type",type,sep="")]][[paste("CV",a,iter)]] <- data.frame(
              do.call(rbind,list(HNB= c(coef(gga_reg_IO),
                                        Theta= summary(gga_reg_IO)$dispersion,
                                        AIC= gga_reg_IO$aic))))
            coef_list[["Gamma"]][["Regression"]][[paste("Type",type,sep="")]][[paste("CV",a,iter)]] <- data.frame(
              do.call(rbind,list(HNB= c(coef(gga_reg),
                                        Theta= summary(gga_reg)$dispersion,
                                        AIC= gga_reg$aic))))
            
          }
        # Save regression results to pred_list
          if(do_regressions){
            df_test$ols_pred <- y_ols
            df_test$rob_pred <- y_rob
            df_test$tob_pred <- y_tobit
            
            df_test$gnb_pred <- y_nb
            df_test$znb_pred <- y_zi.nb
            df_test$hnb_pred <- y_hu.nb
            
            # df_test$gln_pred <- y_gln
            df_test$gga_pred <- y_gga
            
            # df_test$knn_pred <- y_knn
            # df_test$gbm_pred <- y_gbm
            # df_test$svm_pred <- y_svm
            pred_list[[paste("Type",type,sep="")]][[paste("CV",a,iter)]] <- df_test
          }
        # Plot points
          if(a > 0 & iter == 1){
            x_para <- - df_test[,1]
            
            points(x_para,y_tobit,pch=20,cex=cex_para,col="gray75")
            points(x_para,y_ols,pch=20,cex=cex_para*0.9,col="gray50")
            points(x_para,y_rob,pch=20,cex=cex_para*0.9*0.9,col="gray25")
            
            points(-x_para, y_zi.nb,pch=20,cex=cex_para,col="palevioletred2")
            points(-x_para, y_nb,pch=20,cex=cex_para*0.9,col="peachpuff")
            points(-x_para, y_gga,pch=20,cex=cex_para*0.9*0.9,col="blue")
            # points(-x_para, y_hu.nb,pch=20,cex=cex_para*0.9*0.9,col="indianred3")
            
            # points(-x_para,y_gbm,pch=20,cex=cex_para,col="palevioletred2")
            # points(-x_para, y_knn,pch=20,cex=cex_para*0.9,col="peachpuff")
            # points(-x_para, y_svm,pch=20,cex=cex_para*0.9*0.9,col="indianred3")
            # points(-x_para,y_gam,pch=20,cex=cex_para*0.9*0.9,col="lavenderblush1")
          }
        # Add legends
          if(do_legend & a == 1 & iter == 1){
            # legend("topright",bty="n",col=c(NA,"peachpuff","palevioletred2","indianred3"),legend=c("non-linear",paste("KNN ",k_percent,"%",sep=""),"GBM","SVM"),pch=c(NA,20,20,20),cex=1)
            legend("topright",bty="n",col=c(NA,"peachpuff","palevioletred2","indianred3"),legend=c("count-data","GNB","ZNB","HNB"),pch=c(NA,20,20,20),cex=1)
            legend("topleft",bty="n",legend=c("linear","OLS","ROB","TOB"),col=c(NA,"gray50","gray25","gray75"),pch=c(NA,20,20,20),cex=1)
          }
          }
      }
  }
    
    
    mtext("Predictions",side=2,outer=TRUE,line=1,cex=1.4)
    mtext("Observations",side=1,outer=TRUE,line=1,cex=1.4)
    if(eps){ dev.off() }
    
    # Prepare Cross-Validation Coefs of Regressions
    for(type in c("TypeBO","TypeVC")){
      coef_list[["linear"]][[type]] <- data.frame(do.call(rbind,coef_list[["linear"]][[type]]))
      for(reg in c("Regression","InterceptOnly")){
        for(modelz in c("NB2","ZNB2","HNB2","Gamma")){
          coef_list[[modelz]][[reg]][[type]] <- data.frame(do.call(rbind,coef_list[[modelz]][[reg]][[type]]))
        }
      }
    }
    
    invisible(list(Coef= coef_list, 
                   Pred= pred_list))
}
#  set.seed(99) ; system.time( moa_output <- moa_cv(10) )
  
  # Analyze output
  moa_analyzer <- function(moa_out= moa_output,eps=FALSE){
    outlist <- list()
    df_out <- list()
    if(eps){
      setEPS() ; postscript("Multiple Residuals CrossValidation.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=5)
    }
    par(mfrow=c(1,2),mar=c(2,2,2,1),oma=c(3,3,1,1))
    for(type in c("BO","VC")){
      # Calculate Average Regression Coefficients (OLS/ROB/TOB)
      df_coef <- moa_out$Coef$linear[[paste("Type",type,sep="")]]
      avg_coef <- rbind2summary(df_coef)

      # Calculate Average Regression Coefficients (Count - GNB)
      avg_coef_glm <- list()
      for(modelz in c("NB2","ZNB2","HNB2","Gamma")){
        df_model_IO <- moa_out$Coef[[modelz]]$InterceptOnly[[paste("Type",type,sep="")]]
        avg_model_IO <- rbind2summary(df_model_IO)
        df_model_full <- moa_out$Coef[[modelz]]$Regression[[paste("Type",type,sep="")]]
        avg_model <- data.frame(rbind(avg_model_IO,
                                      rbind2summary(df_model_full)))
        if(modelz == "NB2"){
          avg_model["Pseudo R^2","Mean"] <- 1  - (avg_model["TwoLogLik1","Mean"]/avg_model["TwoLogLik","Mean"])
        }
        if(modelz %in% c("ZNB2","HNB2")){
          avg_model["Pseudo R^2","Mean"] <- 1  - (avg_model["Log_Lik1","Mean"]/avg_model["Log_Lik","Mean"])
        }
        
        avg_coef_glm[[modelz]] <- avg_model
      }
      

      # Extract & Prepare Data
      moa_type <- do.call(rbind,moa_out$Pred[[paste("Type",type,sep="")]])
      moa_type <- moa_type[complete.cases(moa_type),]
      print(paste(type,nrow(moa_type),"predictions"))
      df_cor <- moa_type[,colnames(moa_type) %in% colnames(moa_type)[1] | grepl("pred",colnames(moa_type))]
      df_cor[df_cor < 0] <- 0
      
      # Plot Residuals
      plot(c(-5,10),c(0,0),type="l",lty=2,ylim=c(-10,30),ylab=NA,xlab=NA,col="blue",main=type)
      i <- 1
      colist <- c("gray50","gray25","gray75","peachpuff","palevioletred2","indianred3")
      rows2plot <- seq(1,min(1000,nrow(moa_type)))
      for(reg in colnames(df_cor)[-1][1:6]){
        points(moa_type[rows2plot,reg],(moa_type[rows2plot,1]-moa_type[rows2plot,reg]), ylab=NA,pch=20,col=colist[i],cex=1.25*(0.95^(i-1)))
        i <- i + 1
      }
      # mtext(type,side=2,line=2.5)
      abline(a=0,b=-1,col="blue") ; abline(v=0,col="blue",lty=2)
      legend("topright",bty="n",legend=c(paste(toupper(substring(colnames(df_cor)[-1][1:6],1,3)))),
             col=colist,pch=20,cex=1.3)
      
      # Correlations
      for(come in c("pearson","spearman","kendall")){
        correl <- cor(df_cor,method=come)[-1,1]
        outlist[[paste("Type",type,sep="_")]][[paste("cor_pred",come,sep="_")]] <- correl
        df_cor2 <- -(df_cor-df_cor[,1])
        coresi <- diag(cor(df_cor[,-1],df_cor2[,-1],method=come))
        outlist[[paste("Type",type,sep="_")]][[paste("cor_resi",come,sep="_")]] <- coresi
      }
      
      # Mean Squared Error & Mean Absolute Deviation
      mse <- list() ; mad <- list() ; bias <- list()
      for(cona in colnames(df_cor)[-1]){
        mse[[cona]] <- (mean((df_cor[,1]-df_cor[,cona])^2))
        mad[[cona]] <- mean(abs(df_cor[,1]-df_cor[,cona]))
        bias[[cona]] <- mean(df_cor[,1] - df_cor[,cona])
      }
      outlist[[paste("Type",type,sep="_")]][[paste("MAE")]] <- unlist(mad)
      outlist[[paste("Type",type,sep="_")]][[paste("MSE")]] <- unlist(mse)
      outlist[[paste("Type",type,sep="_")]][[paste("BIAS")]] <- unlist(bias)
      
      df_out$Coef[[type]] <- avg_coef
      df_out$GLM[[type]] <- avg_coef_glm
      df_out$Perf[[type]] <- data.frame(do.call(rbind, outlist[[paste("Type",type,sep="_")]]))
      df_out$Perf[[type]]["VAR",] <- df_out$Perf[[type]]["MSE",] - (df_out$Perf[[type]]["BIAS",])^2
      colnames(df_out$Perf[[type]]) <- toupper(substring(colnames(df_out$Perf[[type]]),1,3))
    }
    mtext("Predictions",side=1,outer=TRUE,line=1,cex=1.3) ; mtext("Residuals",side=2,outer=TRUE,line=1,cex=1.3)
    if(eps){ dev.off()  }
    return(df_out)
  }
#  system.time( moa_ana <- moa_analyzer() )
#  lapply(moa_ana$GLM$BO,function(x)round(x,3))
#  lapply(moa_ana$GLM$VC,function(x)round(x,3))
  
  
  # function to compare fit of intercept only model and to simultaneously determine THETA
  fit_gnb_intercept_only <- function(N=10,simulate=FALSE,plot_it=TRUE,eps=FALSE){
    Q_data <- list()
    overdispersion <- list() ; loli <- list()
    for(i in 1:N){
      print(i)
      df <- creat_reg_df()
      for(type in c("BO","VC")){
        df1 <- df[[type]]
        countz <- 100
        observations <- pmax(0,round(df1$P2C.multi1*countz,0)-0)
        n <- length(observations)
        Q_seq <- seq(0,1,0.05)
        df_quantiles <- data.frame(Emp_Obs= quantile(observations,Q_seq))
        
        # estimate mu & size (NB2)
        nb_obs <- MASS::fitdistr(observations,densfun="negative binomial")
        mu_obs <- nb_obs$estimate["mu"]
        size_obs <- nb_obs$estimate["size"]
        loli[[paste("T",type,sep="_")]][["NB2"]][[paste("Run",i,sep="_")]] <- c(nb_obs$estimate, AIC= AIC(nb_obs))
        overdispersion[[paste("T",type,sep="_")]][[paste("Run",i,sep="_")]] <- size_obs
        df_quantiles[,"NB2_PMF"] <- qnbinom(Q_seq,mu=mu_obs,size = size_obs)

        # test other distributions
        observations2 <- pmax(0.01,df1$P2C.multi1)
        dtb_list <- c("gamma","weibull","lognormal")
        for(dtb in dtb_list){
          error <- try(fit_dtb <- MASS::fitdistr(observations2,densfun=dtb))
          if(class(error) != "try-error"){
            loli[[paste("T",type,sep="_")]][[dtb]][[paste("Run",i,sep="_")]] <- c(fit_dtb$estimate, AIC= AIC(fit_dtb))
          }
        }
        
        if(simulate){
          nb_sim <- rnbinom(n,mu=mu_obs,size = size_obs)
          df_quantiles[,"NB2_sim"] <- quantile(nb_sim,Q_seq)
          for(lepto_para in seq(0.5,5,0.5)){
            gamma_noise <- rgamma(n,lepto_para,lepto_para)
            nb_sim <- rnbinom(n,mu=mu_obs,size = size_obs * gamma_noise)
            df_quantiles[,paste("LP",lepto_para,sep="_")] <- quantile(nb_sim,Q_seq)
          }
        }
        
        df_quantiles$Prob <- Q_seq
        df_quantiles <- df_quantiles[-nrow(df_quantiles),]
        Q_data[[paste("T",type,sep="_")]][[paste("Run",i,sep="_")]] <- df_quantiles
      }
    }
    
    out <- list() # calculate mean quantiles
    out$BO <- as.data.frame(apply(abind::abind(Q_data$T_BO,along = 3),  1:2, mean))
    out$VC <- as.data.frame(apply(abind::abind(Q_data$T_VC,along = 3),  1:2, mean))
    
    for(type in c("BO","VC")){
      for(dtb in c("NB2",dtb_list)){
        loli[[paste("T",type,sep="_")]][[dtb]] <- rbind2summary(data.frame(do.call(rbind,
                                                    loli[[paste("T",type,sep="_")]][[dtb]])))
      }
    }

    
    
    if(plot_it){
      if(eps){
        setEPS() ; postscript("Multiple negative binomial intercept only fit.eps", width = 5, height = 2, family = "Helvetica",pointsize=5)
      }
      
      mean_THETA_bo <- as.numeric(1/mean(overdispersion[[paste("T","BO",sep="_")]]))
      mean_THETA_vc <- as.numeric(1/mean(overdispersion[[paste("T","VC",sep="_")]]))
      
      par(mfrow=c(1,2),mar=c(4,4,2,2),cex=1.3)
      XLIM <- 7
      bo_col <- "royalblue1" ; vc_col <- "maroon1"
      
      plot(out$BO$Emp_Obs/countz,
           out$BO$Prob,
           ylim=c(0,1),xlim=c(0,XLIM),col=bo_col,main="BO",
           xlab="Multiple in %",ylab="CDF",type="b",pch=0,cex=0.7)
      points(out$BO$NB2_PMF/countz,
             out$BO$Prob,
             col=bo_col,type="l")
      curve_p_dtb <- function(type){
        dtb <- "gamma"
        curve(pgamma(x, loli[[paste("T",type,sep="_")]][[dtb]][1,"Mean"],
                     loli[[paste("T",type,sep="_")]][[dtb]][2,"Mean"]),0,XLIM,lty=2,add=TRUE)
        dtb <- "weibull"
        curve(pweibull(x, loli[[paste("T",type,sep="_")]][[dtb]][1,"Mean"],
                       loli[[paste("T",type,sep="_")]][[dtb]][2,"Mean"]),0,XLIM,lty=3,add=TRUE)
        dtb <- "lognormal"
        curve(plnorm(x, loli[[paste("T",type,sep="_")]][[dtb]][1,"Mean"],
                     loli[[paste("T",type,sep="_")]][[dtb]][2,"Mean"]),0,XLIM,lty=4,add=TRUE)
      }
      curve_p_dtb("BO")

      legend("bottomright",bty="n",legend=c("empirical","GNB (intercept only)"),col=bo_col,pch=c(0,NA),lty=c(NA,1),cex=0.9)
      legend("right",bty="n",legend= latex2exp::TeX(paste("$\\Theta_{BO} = $",mean_THETA_bo,sep="")))

      plot(out$VC$Emp_Obs/countz,
           out$VC$Prob,
           ylim=c(0,1),xlim=c(0,XLIM),col=vc_col,main="VC",
           xlab="Multiple in %",ylab="CDF",type="b",pch=0,cex=0.7)
      points(out$BO$NB2_PMF/countz,
             out$BO$Prob,
             col=vc_col,type="l")
      curve_p_dtb("VC")
      
      legend("bottomright",bty="n",legend=c("empirical","GNB (intercept only)"),col=vc_col,pch=c(0,NA),lty=c(NA,1),cex=0.9)
      legend("right",bty="n",legend= latex2exp::TeX(paste("$\\Theta_{VC} = $",mean_THETA_vc,sep="")))
      
      
      if(eps){ 
        dev.off()
      }
    }
    
    invisible(loli)
  }
#  set.seed(99) ; system.time( out <- fit_gnb_intercept_only(10) )
  
  
  # preapare moa_ana for .csv export for Lyx
  create_table_csvs <- function(ma= moa_ana,write.it=FALSE){
    df_cc <- do.call(rbind,ma$NB2)
    df_cc$model <- "negative_binomial"
    df_c <- do.call(rbind,ma$Coef)
    df_c$model <- "ols_tob_rob"
    bagging_coefs <- rbind(df_c,df_cc)
    bagging_coefs[,1:5] <- apply(bagging_coefs[,1:5],2,function(x)round(x,2))
    if(write.it){
      write.csv(bagging_coefs,"LM_Coefs.csv")
      write.csv(round(rbind(ma$Perf$BO,ma$Perf$VC),2),"LM_Perf.csv")
      saveRDS(ma ,"z_CV_Output.rds")
    }
  }
  # create_table_csvs()
}
if(FALSE){
  moa_ana1 <- readRDS("z_CV_Output.rds")
  
  predict_MoaCV <- function(df = bo_vc, cvpara= moa_ana1){
    b_bo <- c(Intercept          = cvpara$NB2$BO$Mean[1],
              Holding_Period     = 0,
              Time2Exit_ZS       = cvpara$NB2$BO$Mean[2],
              High_Yield_Spread  = cvpara$NB2$BO$Mean[3],
              MSCI_Multiple      = cvpara$NB2$BO$Mean[4])
    size_bo <- cvpara$THETA_BO_sampling
    b_vc <- c(Intercept          = cvpara$NB2$VC$Mean[1],
              Holding_Period     = cvpara$NB2$VC$Mean[2],
              Time2Exit_ZS       = cvpara$NB2$VC$Mean[3],
              High_Yield_Spread  = 0,
              MSCI_Multiple      = cvpara$NB2$VC$Mean[4])
    size_vc <- cvpara$THETA_VC_sampling
    
    out_list <- list()
    for(i in seq(nrow(df))){
      X_i <- c(Intercept          = 1,
               Holding_Period     = df$Holding_Period[i],
               Time2Exit_ZS       = df$Time2Exit[i] + df$ZombieStage[i],
               High_Yield_Spread  = df$ML_HYOAS.quarter[i],
               MSCI_Multiple      = df$MSCI.Multiple.Exit[i] - 1)
      type <- df$Fund_InvestTypes[i]
      if(type == "BO"){
        b_i <- b_bo
      }else{
        b_i <- b_vc
      }
      
      mu_i <- exp(sum(X_i * b_i))
      size_i <- ifelse(type == "BO", size_bo, size_vc)
      
      out_list[[i]] <- c(Mu= mu_i, Size= size_i)
    }
    
    out <- data.frame(do.call(rbind,out_list))
    
    return(out)
  }
  RQ_MoaCV <- function(df = bo_vc, cvpara= moa_ana){
    pred <- predict_MoaCV(bo_vc,moa_ana)
    realized <- round(df$P2C.multi1 * 100,0)
    RQ <- as.numeric(pnbinom(realized, mu= pred$Mu, size= pred$Size))
    
    out <- data.frame(RQ_Multiple.nb2 = RQ)
    
    return(out)
  }
  
  bo_vc$RQ_Multiple.nb2 <- RQ_MoaCV()
}
## e) Data Summary & Vizualisation (for publication) ---------
if(FALSE){
  length(levels(as.factor(G.p2c$FundInv_FundId))) # of realized funds
  length(levels(as.factor(G.p2c$Fund_Emi_ID))) # of exited companies
  
  # 3D density chart: Bivariate Input Variable (MOIC vs. Holding Period)
  # setEPS() ; postscript("Bivariate Random Variable.eps", width = 5, height = 5, family = "Helvetica",pointsize=7)
  # png(filename = "Bivariate Random Variable.png", width = 1500,height = 1000,pointsize = 18)
  plot_bivariate_input <- function(G.sum,id_bo,id_vc){
    par(mfrow=c(1,2),mar=c(1,2,1,2),cex=0.7,oma=c(1,1,5,1))
    for(fund_type in c("BO","VC")){
      kernel_grid <- 100
      timing_correction <- TRUE
      
      if(fund_type == "BO"){
        ID_set <- id_bo
      }else{
        ID_set <- id_vc
      }
      
      sub_df2 <- G.sum[G.sum$Fund_Emi_ID %in% ID_set,c("MOIC","HoPi")]
      aCDP <- as.copuladata(pobs(sub_df2))
      aCDP <- aCDP[complete.cases(aCDP),]
      
      d3 <- kde2d(aCDP$MOIC, aCDP$HoPi, n=kernel_grid)
      
      rotate <- 180 - 45 - 90 - 90 + 45
      persp3D(x= d3$x, y=d3$y, z=d3$z,phi = -90+100 +45, theta = 270+45+rotate,colkey = FALSE,
              bty="u", col.panel=NULL,col.axis="darkgray",col.grid="darkgray",lwd=2,resfac=1,
              image=F,xlab="Holding Period",ylab="Multiple",zlab="Density",
              main=paste(fund_type))
      
      rm(aCDP,sub_df2)
    }
    # mtext("Timing vs Multiple: 2-dimensional Rosenblatt transformation (realized quantiles)",side=3,line=2.5,outer=TRUE,cex=1.5)
    par(par_default)
  }
  plot_bivariate_input(G.sum,id_bo,id_vc)
  #dev.off()
  
  
  
  
  
  # png(filename = "Bagging_Residuals_BOVC.png", width = 700,height = 800,pointsize = 16)
  # pdf("Bagging_Residuals_2.pdf", width = 5, height = 4, family = "Helvetica",pointsize=4)
  # setEPS() ; postscript("Bagging_Residuals_3.eps", width = 5, height = 4, family = "Helvetica",pointsize=6)
  par(mfrow=c(2,1))
  sun(pre_bo$Resi_Bag,main="BO: Bagging Residuals",Xlim = c(-5,15))
  #sun(pre_bo$Real-pre_bo$Pred_Rob,main="BO: Robust Residuals",Xlim = c(-5,65))
  #sun(pre_bo$Real-pre_bo$Pred_OLS,main="BO: OLS Residuals",Xlim = c(-5,65))
  ## VC
  sun(pre_vc$Resi_Bag,main="VC: Bagging Residuals",Xlim = c(-5,15))
  #sun(pre_vc$Real-pre_vc$Pred_Rob,main="VC: Robust Residuals",Xlim = c(-5,50))
  #sun(pre_vc$Real-pre_vc$Pred_OLS,main="VC: OLS Residuals",Xlim = c(-5,50))
  # dev.off()
  par(par_default)
  
  
  
  
  # Scatterplot (Observations x Predictions)
  # png(filename = "Bagging_Scatter_Obs_Predict.png", width = 350,height = 175,pointsize = 5,res=100)
  # pdf("Bagging_Scatter_Obs_Predict2.pdf", width = 5, height = 2.5, family = "Helvetica",pointsize=8)
  # setEPS() ; postscript("Bagging_Scatter_Obs_Predict2.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=7)
  plot_scatter_multiple <- function(pre_bo,pre_vc){
    par(mfrow=c(1,2),mar=c(2,2,3,1),oma=c(2,2,0,1))
    plot(pre_bo$Obs,pre_bo$Pred_Tobit + mean(pre_bo$Resi_Tobit),col="blue",pch=20,cex=0.8,
         ylim=c(0,5),xlim=c(0,20),main="BO",ylab=NA,xlab=NA)
    points(pre_bo$Obs,pre_bo$Pred_OLS + mean(pre_bo$Resi_OLS),col="green",pch=20,cex=0.8)
    points(pre_bo$Obs,pre_bo$Pred_Rob + mean(pre_bo$Resi_Rob),col="red",pch=20,cex=0.8)
    abline(a=0,b=1)
    legend("topright",bty="n",col=c("blue","red","green"),pch=20,legend=c("Tobit","Robust","OLS"))
    plot(pre_vc$Obs,pre_vc$Pred_Tobit + mean(pre_vc$Resi_Tobit),col="blue",pch=20,cex=0.8,
         ylim=c(0,5),xlim=c(0,20),main="VC",ylab=NA,xlab=NA)
    points(pre_vc$Obs,pre_vc$Pred_Rob + mean(pre_vc$Resi_Rob),col="red",pch=20,cex=0.8)
    points(pre_vc$Obs,pre_vc$Pred_OLS + mean(pre_vc$Resi_OLS),col="green",pch=20,cex=0.8)
    abline(a=0,b=1)
    legend("topright",bty="n",col=c("blue","red","green"),pch=20,legend=c("Tobit","Robust","OLS"))
    mtext("Observations",side=1,outer = TRUE,line = 0.5,cex=1.1)
    mtext("Prediction + Bias Correction",side=2,outer = TRUE,line = 0.5,cex=1.1)
  }
  plot_scatter_multiple(pre_bo,pre_vc)
  par(par_default)
  # dev.off()
  
  
  
  
  # Categorical Stratification (Residual Mean)
  table(bo_ss$GICS_Sector,bo_ss$Fund_Region)
  (aggregate(pre_bo$Resi_Bag, list(bo_ss$GICS_Sector,bo_ss$Fund_Region), function(x) round(mean(x),4)))
  table(vc_ss$GICS_Sector,vc_ss$Fund_Region)
  (aggregate(pre_vc$Resi_Bag, list(vc_ss$GICS_Sector,vc_ss$Fund_Region), function(x) round(mean(x),4)))
  
  
  # P(x==0) vs mu in NB2
  par(mar=c(4.5,4.5,1,1))
  seq_x <- seq(0,400) ; Y <- 10
  plot(seq_x/100,pnbinom(Y,mu=seq_x,size=moa_ana1$THETA_BO_sampling),ylim=c(0,1),type="l",col="royalblue1",
       ylab=latex2exp::TeX(paste("$\\textit{P}(\\hat{Y} $\\leq$ 0.1) $",sep="")),
       xlab= latex2exp::TeX(paste("$\\mu $",sep=""))) 
  lines(seq_x/100,pnbinom(Y,mu=seq_x,size=moa_ana1$THETA_VC_sampling),ylim=c(0,1),col="maroon1")  
  abline(h=c(0,0.2,0.4),v=0,col="grey",lty=3)
  legend("topright",bty="n",legend = c("BO","VC"),col=c("royalblue1","maroon1"),lty=1)
  
  
  ######  BiVariate
  par(mfrow=c(3,2),mar=c(1,2,1,2),cex=0.7,oma=c(1,0,2,0))
  sun(bo_vc$RQ_RLT.weibull[bo_vc$Fund_InvestTypes == "BO"], Xlim = c(0,1), Ylim = c(0,2),main="BO: Timing (Weibull)") # in realized set there are just few late observations
  sun(bo_vc$RQ_RLT.weibull[bo_vc$Fund_InvestTypes == "VC"], Xlim = c(0,1), Ylim = c(0,2),main="VC: Timing (Weibull)") # in realized set there are just few late observations
  sun(bo_vc$RQ_Multiple.nb2[bo_vc$Fund_InvestTypes == "BO"], Xlim = c(0,1), Ylim = c(0,2),main="BO: Multiple (NB2)") # in realized set there are just few late observations
  sun(bo_vc$RQ_Multiple.nb2[bo_vc$Fund_InvestTypes == "VC"], Xlim = c(0,1), Ylim = c(0,2),main="VC: Multiple (NB2)") # in realized set there are just few late observations
  

  # 3D density chart (Residual Dependency via "realized quantiles" two-dimensional Rosenblatt transformation)
  # we estimated non-parametric distribution functions for multiple and timing
  
  # png(filename = "Bivariate Residual Density (Type Region).png", width = 1500,height = 1000,pointsize = 18)
  # setEPS() ; postscript("Bivariate_Model_Dependencies.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=7)
  plot_bivariate_output <- function(df_in= bo_vc,kernel_grid= 100,timing_correction= TRUE){
  par(mfrow=c(1,2),mar=c(1,2,1,2),cex=0.7,oma=c(1,0,2,0))
    
    for(timing_correction in c(TRUE)){
      for(fund_type in c("BO","VC")){
        
        bo_vc <- as.data.frame(df_in %>% group_by(Fund_Emi_ID) %>% sample_n(1, replace=TRUE))
        # bo_vc <- df_in[df_in$FundInv_Quarter == df_in$Investment_Date,]
        
        multi_var <- "RQ_Multiple.dgh"  # "RQ_Multiple.nb2"
        
        sub_df2 <- bo_vc[complete.cases(bo_vc[,c("RQ_RLT.weibull",multi_var)])
                         & (is.finite(bo_vc$RQ_RLT.weibull)) &
                           bo_vc$Fund_InvestTypes == fund_type,
                         c("RQ_RLT.weibull",multi_var)]
        aCDP <- as.copuladata(pobs(sub_df2))
        sub_df2$RQ_RLT.weibull_ascopula <- aCDP$RQ_RLT.weibull
        sub_df2[,paste(multi_var,"ascopula",sep = "_")] <- aCDP[,multi_var]
        
        
        if(timing_correction){
          d3 <- MASS::kde2d(sub_df2$RQ_RLT.weibull_ascopula,sub_df2[,paste(multi_var,"ascopula",sep = "_")],n=kernel_grid)
          label <- "uniformized"
        }else{
          d3 <- MASS::kde2d(sub_df2$RQ_RLT.weibull,sub_df2[,multi_var],n=kernel_grid)
          label <- "predicted"
        }
        
        rotate <- 180 - 45 - 90 - 90 +45
        persp3D(x= d3$x, y=d3$y, z=d3$z,phi = -90+100 +45, theta = 270+45+rotate,colkey = TRUE,
                bty="u", col.panel=NULL,col.axis="darkgray",col.grid="darkgray",lwd=2,resfac=1,
                image=F,xlab="Timing",ylab="Multiple",zlab="Density",
                main=paste(fund_type,": ",label,sep=""))
      }
    }
    
    # mtext("Bivariate Rosenblatt Transformation of Marginal Estimation Models",side=3,line=2.5,outer=TRUE,cex=1.5)
    par(par_default)
  }
  plot_bivariate_output(bo_vc)
  # dev.off()
  par(par_default)
  
}

## f) Simulation Engine --------
create_random_input <- function(n=1,age=NA,cox_mo=cox_model,hp_shift=0){
  # create input dataframe
  df <-data.frame(Fund_InvestTypes= sample(c("BO","VC"),n,replace=TRUE),
                  AssetAge= sample(seq(1,6),n,replace=TRUE),
                  MOIC.dyn=rlnorm(n,0,0.5),
                  ML_HYOAS.quarter= rep(sample(sample(trindex$ML_HYOAS[!(is.na(trindex$ML_HYOAS))],1),1),n),
                  MSCI.Multiple.Exit=rlnorm(n,0,0.25),
                  # Sector= sample(df_sector$Sector,n,replace=TRUE),
                  Fund_Region= sample(c("EU","US"),n,replace = TRUE)
  )
  df$AssetAge <- 0
  
  df$RVPI <- pmax(0.1,df$MOIC.dyn - rlnorm(n,-2,0.6))
  df$RVPI_1 <- df$RVPI - 1
  df$MSCI.Multiple.Exit_1 <- df$MSCI.Multiple.Exit - 1
  
  # df <- merge(df,df_sector[,colnames(df_sector) %in% c("Sector","CoxSector","Sector_BO_m","Sector_VC_m")], by="Sector",all.x=TRUE)
  
  if(!is.na(age)) df$AssetAge <- rep(age,n)
  # df$GLOQuarter <- sample(gdp$GLO,1)
  df$FundAgeAI <-  max(df$AssetAge) - df$AssetAge
  df$Holding_Period <- df$AssetAge + hp_shift
  df$Moic.trans2 <- dgamma(df$MOIC.dyn,shape=2,scale=1)
  df$CoxFactor <- exp(rowSums(data.frame(mapply(`*`,df[,c("FundAgeAI","ML_HYOAS.quarter","Moic.trans2")],cox_model$Regression$coefficients))))
  df$NAV <- 1 # rlnorm(n,meanlog = 10,sd=2)
  
  invisible(df)
}
esitmation <- function(iteration=1,companies=5,bo_vc_frame=bo_vc){
  startz <- Sys.time()
  # 0.) create input data
  df_input <- create_random_input(companies)

  # Simulate
  # df_out <- data.frame(Timing=NA,Multiple=NA,ID= NA,Iteration=NA)
  df_out <- list()
  for(i in seq(iteration)){
    # Macro Variables for Iteration
    msci_multiple <- cumprod(sample(x= msci_world$Return_monthly,size= 400,replace = TRUE))
    df_input$CoxFactor <- exp(rowSums(data.frame(mapply(`*`,
                                      df_input[,c("FundAgeAI","ML_HYOAS.quarter","Moic.trans2")],
                                      cox_model$Regression$coefficients))))
    # Prepapre bi_dependence
    multi_var <- "RQ_Multiple.dgh" # "RQ_Multiple.nb2"
    bi_dep_list <- list()
    for(type in c("BO","VC")){
      bo_vc_frame2 <- as.data.frame(bo_vc_frame %>% group_by(Fund_Emi_ID) %>% sample_n(1, replace=TRUE))
      bi_dependence <- bo_vc_frame2[bo_vc_frame2$Fund_InvestTypes == type, c("RQ_RLT.weibull",multi_var)]
      bi_dependence <- data.frame(pobs(bi_dependence))
      bi_dep_list[[type]] <- bi_dependence[complete.cases(bi_dependence),]
    }
  
    for(id in seq(nrow(df_input))){
      df_in <- df_input[id,]
      type <- as.character(df_in$Fund_InvestTypes)
      region <- as.character(df_in$Fund_Region)
      nav <- df_in$NAV
      label <- paste(region,type,sep="_")
      
      # 1.) Estimating Dependence
      bi_dep <- bi_dep_list[[type]][sample(nrow(bi_dep_list[[type]]),1),]
      
      # 2.) Estimate Timing
      # timing <- sim_cox2(t=df_in$Holding_Period, strata=type, CoxF=df_in$CoxFactor,real_quantile=bi_dep$RQ_RLT.weibull)
      if(FALSE){
        U_in <- bi_dep$RQ_RLT.weibull * S_wb(df_in$Holding_Period,type,cox_wb_zero)^(df_in$CoxFactor)
        timing <- sim_wb(U_in,type, log(df_in$CoxFactor), cox_wb_zero)
      }else{ # intercept only weibull
        U_in <- bi_dep$RQ_RLT.weibull * S_wb(df_in$Holding_Period,type,cox_wb_mean)^(1)
        timing <- sim_wb(U_in,type, 0, cox_wb_mean)
      }
      df_in$Time2Exit <- timing
      zombie_hp <- 10
      df_in$ZombieStage <- pmax(zombie_hp,df_in$Holding_Period + df_in$Time2Exit) - zombie_hp
      df_in$Time2Exit_ZS <- df_in$Time2Exit + df_in$ZombieStage
      
      pos_msci <- max(1,min(400,round(timing*12,0)))
      df_in$MSCI.Multiple.Exit <- msci_multiple[pos_msci]
      df_in$MSCI.Multiple.Exit_1 <- df_in$MSCI.Multiple.Exit - 1
      
      # 3. Estimate Multiple
      if(FALSE){
        if(type == "BO"){
          bag_model <- bo.bag
        }else{
          bag_model <- vc.bag
        }
        
        betas <- re.bag.coefs(bag_model)
        X_vector <- c(1,
                      df_in$Holding_Period,
                      df_in$Time2Exit,
                      df_in$ZombieStage,
                      (df_in$RVPI-1),
                      df_in$ML_HYOAS.quarter,
                      (df_in$MSCI.Multiple.Exit-1))
        
        multiple <- sum(betas * X_vector)
        residualz <- bo_vc_frame$Resi_Multiple[bo_vc_frame$Fund_InvestTypes == type]
        np_distribution <- multiple + residualz
        np_distribution <- as.numeric(unlist(np_distribution))
        np_distribution <- pmax(0,np_distribution) # limited-liability floor
        realized_multiple <- as.numeric(quantile(np_distribution,bi_dep[,multi_var]))
      }else{
        realized_multiple <- predict_dgh(sim_in= df_in,sim_out= TRUE, runo=bi_dep[,multi_var])
      }
      
      
      # Create Output for this Iteration
      df_out[[paste(i,id,sep = "_")]] <- c(timing,realized_multiple,id,i,nav)
    }
  }
  
  df_out <- data.frame(do.call(rbind,df_out))
  names(df_out) <- c("Timing","Multiple","ID","Iteration","NAV")
  
  # Result
  print(">>>> Input:")
  # df_input$GLOQuarter <- gdp_glo # initial GDP estimate
  print(df_input[,c("Fund_InvestTypes","Fund_Region","FundAgeAI","Holding_Period","MOIC.dyn","ML_HYOAS.quarter","NAV")])
  
  if(iteration > 1){
    out <- data.frame(MULTIPLE=mean(df_out$Multiple) ,TIMING=mean(df_out$Timing))
    rownames(out) <- c("Mean")
    
    print(out)
    print(Sys.time() - startz)
    invisible(df_out)
  }else{
    print("Multiple: summary statistics")
    par(par_default)
    sun(np_distribution, main="Estimated Multiple Distribution")
    abline(v=realized_multiple,col="orange",lwd=3)
    print(">>>>>") ; print(">>>>>")
    
    out <- data.frame(timing=c(bi_dep$RQ_Timing,timing),
                      muliple=c(bi_dep$RQ_Multiple.nb2,realized_multiple))
    rownames(out) <- c("Quantile","Value")
    return(out)
  }
  
  
  
}
plot_heatmap <- function(x1,y1,n=10){
  df <- data.frame(X= x1, Y= y1)
  df <- df[complete.cases(df) & !is.infinite(df$X) & !is.infinite(df$Y),]
  x1 <- df$X
  y1 <- df$Y
  
  dens3d <- kde2d(x1,y1,n)
  de3d <- expand.grid(dens3d$x,dens3d$y)
  de3d$Z <- as.vector(dens3d$z)
  colr <- rev(heat.colors(10))
  de3d$col <- cut(de3d$Z, breaks= 10, labels=  colr)
  # levelplot(Z ~ Var1*Var2, data=de3d)
  plot(de3d$Var1,de3d$Var2,pch=".")
  for(cool in levels(de3d$col)){
    df_uni2 <- de3d[de3d$col == cool,]
    points(df_uni2[,1], df_uni2[,2], pch= 15, cex=3,col= cool)
  }
  points(x1,y1,pch=20,cex=0.33,col="grey")
  legend.col(col= colr, lev= de3d$Z)
}
plot_heatmap(bo_vc$RQ_RLT.weibull[bo_vc$Fund_InvestTypes == "VC"],bo_vc$RQ_Multiple.dgh[bo_vc$Fund_InvestTypes == "VC"])
scatterhist = function(df){
  # https://www.r-bloggers.com/example-8-41-scatterplot-with-marginal-histograms/
  
  par(par_default)
  par(oma=c(2,2,4,4))
  
  x= df[,1]
  xlab= names(df)[1]
  y= df[,2]
  ylab= names(df)[2]
  
  zones=matrix(c(2,4,1,3), ncol=2, byrow=TRUE)
  layout(zones)
  
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  
  # Chart 1
  par(mar=c(3,3,1,1))
  plot(x,y)
  abline(h=mean(y),v=mean(x),col="red",lwd=1)
  
  # Chart 2
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  
  # Chart 3
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  
  # Chart 4
  par(mar=c(3,3,1,1))
  df_uni <- as.data.frame(pobs(df))
  make_heatmap <- TRUE
  if(make_heatmap){
    plot_heatmap(df_uni$Timing,df_uni$Multiple)
  }else{
    k_cluster <- 3
    plot(df_uni,pch=".")
    df_uni$Cluster <- kmeans(df_uni,k_cluster)$cluster
    for(cool in 1:k_cluster){
      df_uni2 <- df_uni[df_uni$Cluster == cool,]
      points(df_uni2[,1], df_uni2[,2], col= cool+1)
    }
  }
  
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE)
  mtext(ylab, side=2, line=1, outer=TRUE)
  mtext("Exit Simulation",side=3,line=-3,cex=2,outer=TRUE)
}
cf_paths <- function(df,eps=FALSE){
  df$CashFlow <- df$NAV * df$Multiple
  nav_start <- sum(df$NAV[df$Iteration==1])
  agg_cf <- aggregate(df$CashFlow, by=list(Category=df$Iteration),FUN=sum)
  agg_cf$CF2NAV <- agg_cf$x / nav_start
  agg_cf$CF2NAV_rounded <- floor(agg_cf$CF2NAV)
  max_cf <- max(agg_cf$x)/nav_start
  
  if(eps){
    setEPS() ; postscript("CF_Simulation_3.eps", width = 5, height = 3.5, family = "Helvetica",pointsize=7)
  }
  par(mar=c(4.5,4.5,4,2),mfrow=c(1,1))
  limit_x <- max(df$Timing)+3.5
  plot(0,1, type="s", lwd=2,col="blue", ylim=c(-1,max_cf) ,xlim=c(-1,limit_x),
       xlab="Time",ylab="Multiple on NAV(t)",main="Cash Flow Paths (Fund/Portfolio)")
  for(iter in unique(df$Iteration)){
    df_iter <- df[df$Iteration == iter, c("Timing","CashFlow")]
    df_iter <- rbind(df_iter,c(0,0))
    df_iter <- df_iter[order(df_iter$Timing),]
    df_iter$CumCF <- cumsum(df_iter$CashFlow)
    lines(df_iter$Timing,df_iter$CumCF/nav_start,col="darkgrey",type="s")
    points(df_iter$Timing,df_iter$CumCF/nav_start,col="darkgrey",pch=20,cex=1)
    
    cf2nav <- round(sum(df_iter$CashFlow)/nav_start,2)
  }
  # Final Fund/Portfolio CashFlow-Multiple on NAV
  points(y=agg_cf$CF2NAV,x=rep(max(df$Timing)+1,length(agg_cf$CF2NAV)),col="red",pch=18,cex=1)
  agg_cf_unq <- agg_cf[!(duplicated(agg_cf$CF2NAV_rounded)),]
  for(j in seq(nrow(agg_cf_unq))){
    text(max(df$Timing)+1.75,agg_cf_unq$CF2NAV,round(agg_cf_unq$CF2NAV,1),cex=0.8)
  }
  
  abline(h=1,col="blue",lwd=2)
  points(0,1,col="blue",pch=15)
  abline(h=0,v=0,col="black",lty=3)
  abline(h=mean(agg_cf$x)/nav_start,col="red",lty=2,lwd=1)
  
  # Inverse of ECDF = Quantile Function
  x <- seq(0,1,0.01)
  y <- quantile(agg_cf$CF2NAV,x)
  lines(-x,y,col="green",lwd=2) # ECDF
  abline(v=-1,h=-1,lty=3)
  point_size <- 1.25
  points(-0.75,quantile(agg_cf$CF2NAV,0.75),pch=19,col="green",cex=point_size)
  points(-0.5,quantile(agg_cf$CF2NAV,0.5),pch=19,col="green",cex=point_size)
  points(-0.25,quantile(agg_cf$CF2NAV,0.25),pch=19,col="green",cex=point_size)


  agg_timing <- aggregate(df$Timing, by=list(Category=df$Iteration),FUN=max)
  timing_ecdf1 <- ecdf(agg_timing$x)
  curve(-timing_ecdf1(x),from=0,to=max(df$Timing),col="green",lwd=2,add=TRUE)
  points(quantile(agg_timing$x,0.25),-0.25,pch=19,col="green",cex=point_size)
  points(quantile(agg_timing$x,0.5),-0.5,pch=19,col="green",cex=point_size)
  points(quantile(agg_timing$x,0.75),-0.75,pch=19,col="green",cex=point_size)

  timing_ecdf2 <- ecdf(df$Timing)
  curve(-timing_ecdf2(x),from=0,to=max(df$Timing),col="orange",lwd=2,add=TRUE)
  points(quantile(df$Timing,0.25),-0.25,pch=19,col="orange",cex=1.5)
  points(quantile(df$Timing,0.5),-0.5,pch=19,col="orange",cex=1.5)
  points(quantile(df$Timing,0.75),-0.75,pch=19,col="orange",cex=1.5)

  if(eps){dev.off()}
}

# Start Simulation
set.seed(99)
df_out1 <- esitmation(100,10)
png(filename = "CF Simulation 3.png",width = 480*1.5, height = 480, units = "px", pointsize = 14)
cf_paths(df_out1)
dev.off()

scatterhist(df_out1[df_out1$ID ==1, 1:2])
reshape(df_out1,v.names=c("Timing","Multiple"),timevar="ID",idvar="Iteration",drop="NAV",direction="wide")

#df_out2 <- df_out1[df_out1$Iteration == 233,]
#View(df_out2)
#sum(df_out2$Multiple * df_out2$NAV)/sum(df_out2$NAV)
