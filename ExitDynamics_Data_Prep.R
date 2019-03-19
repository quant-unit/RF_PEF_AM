####################################
#### Exit Dynamics: Data Preparation
####################################
'
This code creates the input file used within ExitDynamics_Timing.R & ExitDynamics_Multiple.R
Output is stored as ExitDynamics_Data_V0.RData in the data.wd
'
## 0) Load (& install) required Packages -------
rm(list=ls()) # remove workspace objects
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

##### ROOT -------------
if(getwd() == "/Users/christausch/Dropbox/Project D/3_R_Pro_D/Exit_Dynamics/Code"){
  computer <- "mac"
}
if(getwd() == "C:/Users/christian.tausch/Dropbox/Project D/3_R_Pro_D/Exit_Dynamics/Code"){
  computer <- "win"
}


## 1) read & prepare data (old) -----

if(computer == "win") root <- "/Users/christian.tausch/Dropbox/Project D/3_R_Pro_D/Exit_Dynamics"
if(computer == "mac") root <- "~/Dropbox/Project D/3_R_Pro_D/Exit_Dynamics"

code.wd <- paste(root, "Code", sep="/")
data.wd <- paste(root, "Data", sep="/")
eps.wd <-  paste(root, "EPS_output", sep="/")
csv.wd <-  paste(root, "CSV_output", sep="/")

setwd(code.wd)
source("Useful_Functions.R")

if(FALSE){
  if(computer == "win") setwd("~/Data/Company Level/Gerd Export")
  if(computer == "mac") setwd("~/Desktop/Company Level/Gerd Export")
  G <- read.csv2("170601_export.csv")
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
  
  G <- data.frame(G %>% group_by(Fund_Emi_ID) %>% 
                    filter(!duplicated(Holding_Period)))
  
  g.sum <- data.frame(G %>% group_by(Fund_Emi_ID) %>% summarise(
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
  
  table(g.sum$Industry,g.sum$Region)
  
  g.sum$HoPi <- as.numeric(g.sum$ExitDate - g.sum$InvDate)/365.25
  g.sum$GeoR <- (g.sum$MOIC)^(1/g.sum$HoPi)
  
  # Load Public Data
  setwd(data.wd)
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
  
  
  
  g.sum <- merge(g.sum,trindex[,c(1,2)],by.x = "MaxDate",by.y="Date")
  g.sum <- merge(g.sum,trindex[,c(1,2)],by.x = "MinDate",by.y="Date",suffixes = c(".Exit",".Entry"))
  g.sum$MSCI.Multiple <- g.sum$MSCI.World.Net.Return.Daily.Exit/g.sum$MSCI.World.Net.Return.Daily.Entry
  gdp <- read.csv2("Cum GDP April 2017.csv",dec=".")
  g.sum$YearInvest <- as.numeric(format(g.sum$InvDate,"%Y"))
  g.sum <- merge(g.sum,gdp[,c(1,4)],by.x = "YearInvest",by.y="Year")
  
  sum(is.na(g.sum$HoPi))
  
  by(g.sum,format(g.sum$InvDate,"%Y"),function(x) median(x$MOIC,na.rm=T))
  
  quantile(g.sum$Geo,seq(0.95,1,0.002),na.rm=T)
  ## b-1) Timing: Create G.sucox
  G.sucox <- g.sum
  G.sucox$Ti <- ifelse(is.na(g.sum$HoPi), as.numeric(g.sum$MaxDate-g.sum$MinDate)/365.25 ,g.sum$HoPi)
  G.sucox$Ev <- ifelse(is.na(g.sum$HoPi), 0,1)
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
  # G.p2c <- G.p2c_BACKUP
  
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
  sun(G.p2c$Distance2Exit)
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
  # G.p2c$P2C.multi1[!is.finite(G.p2c$P2C.multi1) | G.p2c$FMV_C2come < floor.fmv] <- NA
  # G.p2c$P2C.multi2[!is.finite(G.p2c$P2C.multi2) | G.p2c$FMV_C2come < floor.fmv] <- NA
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
  setwd(data.wd)
  load("z_b2.RData")
}


## 2) Prep g.p2c ------
colnames(G.p2c)

keep_cols_G.p2c <- c("Fund_Emi_ID","FundInv_FundId","Inv_EmittentId",
                     "FundInv_Quarter","Investment_Date","Exit_Date",
                     "Exit.Status","Holding_Period",
                     "Fund_VintageYear","Fund_InvestTypes","Fund_Region",
                     "GICS_Sector","Country_Pfc",
                     "MOIC.dyn","MOIC.end","DPI","RVPI_1",
                     "FundAge","FundAgeAI",
                     "T1","T2","Ev","Ev2",
                     "P2C.multi1",
                     "Distance2Exit", # weighting for Multiple Regression
                     "Time2Exit","ZombieStage","Time2Exit_ZS",
                     "MSCI.Multiple.Exit_1","MSCI.qrtly.return","ML_HYOAS.quarter"
                     )

g.p2c2 <- G.p2c[,keep_cols_G.p2c]
g.p2c2 <- g.p2c2[base::order(g.p2c2$Fund_Emi_ID, g.p2c2$FundInv_Quarter),]

g.p2c2$MOIC.end <- ifelse(g.p2c2$Exit.Status == "unrealized", NA, g.p2c2$MOIC.end) 
g.p2c2$P2C.multi1 <- ifelse(g.p2c2$Time2Exit <= 0.3, NA, g.p2c2$P2C.multi1) # last Multiple observations are irrelevant 


g.p2c2$Fund_ID <- as.factor(sapply(strsplit(g.p2c2$Fund_Emi_ID, "_"), `[`, 1))
g.p2c2$Company_ID <- as.factor(sapply(strsplit(g.p2c2$Fund_Emi_ID, "_"), `[`, 2))


## 3) Prep g.sum ------
g.sum <- G.sum
g.sum$Timing <- ifelse(is.na(g.sum$HoPi),
                       as.numeric(as.Date("2016-12-31") - g.sum$InvDate)/365.25,
                       as.numeric(g.sum$ExitDate - g.sum$InvDate)/365.25)
g.sum$RightCensored <- ifelse(is.na(g.sum$HoPi),1,0)
g.sum$Fund_ID <- as.factor(sapply(strsplit(g.sum$Fund_Emi_ID, "_"), `[`, 1))
g.sum$Company_ID <- as.factor(sapply(strsplit(g.sum$Fund_Emi_ID, "_"), `[`, 2))

## 4.a) Final Filter: Fund of Funds -----
# Prep g.sum for parametric Cox
g.sum <- g.sum[g.sum$InvDate < g.sum$ExitDate | is.na(g.sum$ExitDate), ]
# our model can handle zero time for right censored events
# g.sum <- g.sum[g.sum$Timing > 0,]

sun(table(g.sum$Fund_ID)) # why so many investments for some funds? fund of funds !?!?
fund_ids <- names(table(g.sum$Fund_ID)[table(g.sum$Fund_ID) < 100]) # filter out FoFs
g.sum <- g.sum[g.sum$Fund_ID %in% fund_ids, ]

# filter also g.p2c2
g.p2c2 <- g.p2c2[g.p2c2$FundInv_FundId %in% fund_ids, ]



## 4.b) Final Filter: unrelevant categorical information -------
colnames(g.p2c2)
# g.p2c2 <- g.p2c2[, !(colnames(g.p2c2) %in% c("Fund_Region", "GICS_Sector", "Country_Pfc", "Distance2Exit"))]

colnames(g.sum)
g.sum <- g.sum[, !(colnames(g.sum) %in% c("Industry", "Industry2", "Sector", "Currency"))]


## 5) Prep public.data ----
public.data <- trindex[,c("Date","MSCI.World.Net.Return.Daily","ML_HYOAS")]
public.data$MSCI_monthly_return <- c(NA,diff(public.data$MSCI.World.Net.Return.Daily) / public.data$MSCI.World.Net.Return.Daily[1:(nrow(public.data)-1)])
as.Date(levels(as.factor(g.p2c2$FundInv_Quarter))) %in% public.data$Date

# create output ------
rm(list=setdiff(ls(), c("g.p2c2","public.data","g.sum","data.wd")))

dev.off() # reset graphical parameters
par_default <- par()

setwd(data.wd)
rm("data.wd")

# write.csv2(g.sum, "Static_Multiple_Timing2.csv")

if(FALSE){
  save.image("ExitDynamics_Data_V2.RData")
}