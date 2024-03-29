###########################
## Exit Dynamics: Main File
###########################
## 0) Prologue -------
rm(list=ls()) # remove workspace objects
root <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
setwd('..')
root <- getwd()

wd <- list()
wd$code <- file.path(root, "Exit_Dynamics", "Code")
#wd$data <- file.path(root, "Exit_Dynamics", "Code", "Data")
wd$data.in <- file.path(root, "Exit_Dynamics", "Code", "data_in")
wd$data.out <- file.path(root, "Exit_Dynamics", "Code", "data_out")
wd$eps <-  file.path(root, "Exit_Dynamics", "Code", "eps_out")
wd$eps2 <-  file.path(root, "Exit_Dynamics", "Code", "eps_out_small")
#wd$csv <-  file.path(root, "Exit_Dynamics", "Code", "CSV_output")
rm(root)
dir.create(wd$eps)
dir.create(wd$eps2)

setwd(wd$code)
source("Useful_Functions.R")

# packages used
library(dplyr)
library(xtable)
library(latex2exp)
library(gamlss)
library(gamlss.tr)
library(VineCopula)
library(copula)

create.eps <- FALSE
create.csv <- FALSE
## 00) Load & Reduce Input Data ------
setwd(wd$data.out)
load("ExitDynamics_Data_V3.RData")


g.p2c2 <- g.p2c2[g.p2c2$Fund_InvestTypes %in% c("BO", "VC"), ]
colnames(g.p2c2)

# How many syndicated/club deals?
g.syn <- g.sum[g.sum$Company_ID %in% g.sum$Company_ID[duplicated(g.sum$Company_ID)], ]
g.syn$SynID <- paste(g.syn$Company_ID, g.syn$InvDate, g.syn$ExitDate)
length(g.syn$SynID[duplicated(g.syn$SynID)]) / nrow(g.syn)

g.sum <- g.sum[g.sum$Type %in% c("BO", "VC"), 
               c("Type", "Fund_ID", "InvDate", "ExitDate", "YearInvest", "Vintage", "MOIC")]
g.sum$ExitYear <- as.factor(format(g.sum$ExitDate, "%Y"))

table(g.sum$Type)
table(g.p2c2$Fund_InvestTypes)

## 000) Summary Statistics --------
sun(as.numeric(g.sum$ExitDate - g.sum$InvDate)/365.25)
length(names(table(g.sum$Fund_ID[g.sum$Type == "BO"])))

list(
  BO_Funds = length(levels(as.factor(as.character(g.sum$Fund_ID[g.sum$Type == "BO"])))),
  VC_Funds = length(levels(as.factor(as.character(g.sum$Fund_ID[g.sum$Type == "VC"])))),
  BO_VC_Funds = length(levels(as.factor(as.character(g.sum$Fund_ID[g.sum$Type %in% c("BO", "VC")])))),
  Censoring_Companies = table(is.na(g.sum$ExitDate),g.sum$Type) # exit vs right censored events
)

df.entries <- merge(
  x = data.frame(table(g.sum$YearInvest[g.sum$Type == "BO"]), row.names = TRUE),
  y = data.frame(table(g.sum$YearInvest[g.sum$Type == "VC"]), row.names = TRUE),
  by = "row.names"
)
colnames(df.entries) <- c("Year", "Entries.BO", "Entries.VC")

df.exits <- merge(
  x = data.frame(summary(g.sum$ExitYear[g.sum$Type == "BO"])),
  y = data.frame(summary(g.sum$ExitYear[g.sum$Type == "VC"])),
  by = "row.names", all = TRUE
)
colnames(df.exits) <- c("Year", "Exits.BO", "Exits.VC")
df.enex <- merge(df.entries, df.exits, by = "Year", all = TRUE)
rownames(df.enex) <- df.enex$Year
df.enex$Year <- NULL
xtable::xtable(df.enex)

## 1) Timing ------
par(par_default)
setwd(wd$code)
source("ExitDynamics_Timing.R")

# Non-Parametric Cox Regression (partial likelihood)
np.cox <- timing$NonParaCoxRegression(g.p2c2)

# create RDS & .csv output
setwd(wd$data.out)
if (FALSE) {
  # Parametric Cox Regression
  wb.cox <- timing$ParaCoxRegression(public.data, g.sum)
  saveRDS(wb.cox, "Timing_Para_Cox_Weibull6.rds")
} else {
  wb.cox <- readRDS("Timing_Para_Cox_Weibull6.rds")
}
timing$CreateOutput(wb.cox, create.csv)

# Create Plots
if(TRUE) {
  if(create.eps){
    setwd(wd$eps)
  }
  
  timing$PlotCoxWeibullSurvival(eps = create.eps,
                                para_in = wb.cox, 
                                non_para_in = np.cox)
  
  uni.CDF.timing <- timing$Weibull_Roseblatt(eps= create.eps, 
                                             data_LoLi = g.sum, 
                                             public.data = public.data, 
                                             Para_Cox_Weibull = wb.cox)
  
  
  sun(uni.CDF.timing$BO, Ylim = c(0,2))
  sun(uni.CDF.timing$VC, Ylim = c(0,2))
}

## 2) Multiple -------
par(par_default)
setwd(wd$code)
source("ExitDynamics_Copula2.R")

## 3) Empirical Plot & Visualize mMPP framework -----
empirical.timing.multiple <- function(size = "big"){
  dy.bo <- ecdf(G.p2c$P2C.multi1[G.p2c$Fund_InvestTypes == "BO"])
  st.bo <- ecdf(g.sum$MOIC[g.sum$Type == "BO"])
  dy.vc <- ecdf(G.p2c$P2C.multi1[G.p2c$Fund_InvestTypes == "VC"])
  st.vc <- ecdf(g.sum$MOIC[g.sum$Type == "VC"])
  
  if(create.eps){
    old.wd <- getwd()
    
    if(size == "big"){
      setwd(wd$eps)
      setEPS()
      postscript("Empirical-Timing-Multiple.eps", width = 5, height = 2,
                 family = "Helvetica", pointsize = 2)
    } else {
      setwd(wd$eps2)
      setEPS()
      postscript("Empirical-Timing-Multiple_small.eps", 
                 width = 3.25, height = 3.25/2, 
                 family = "Helvetica", pointsize = 1)
    }
    
    if(size == "bigger"){
      setwd(wd$eps)
      setEPS()
      postscript("Empirical-Timing-Multiple_bigger.eps", width = 6, height = 4.5/2,
                 family = "Helvetica", pointsize = 3)
    }
    
  }
  cex.par <- 1.3
  if(size == "small") cex.par <- 1
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2, 1), cex = cex.par)
  plot(np.cox$BO$KaMe, main = "Timing", xlab = "Timing", ylab = "Survival Function", conf.int = FALSE)
  lines(np.cox$VC$KaMe, xlab = "Timing", ylab = "Survival Function", conf.int = FALSE, lty=2)
  legend("topright", bty = "n", legend = c("BO", "VC"), lty = c(1,2))
  
  curve(st.bo, 0, 10, ylim = c(0, 1), xlab = "Multiple", ylab = "ECDF", main = "Multiple")
  # curve(dy.bo, 0, 10, add = TRUE, lty = 2)
  curve(st.vc, 0, 10, ylim = c(0, 1), xlab = "Multiple", ylab = "ECDF", main = "VC", add = TRUE, lty=2)
  # curve(dy.vc, 0, 10, add = TRUE, lty = 2)
  legend("bottomright", bty = "n", legend = c("BO", "VC"), lty = c(1,2))
  
  if(create.eps){ 
    dev.off() 
    setwd(old.wd)
  }
  
}
empirical.timing.multiple("big")
empirical.timing.multiple("small")
empirical.timing.multiple("bigger")

empirical.multiple <- function() {
  st.vc <- ecdf(g.sum$MOIC[g.sum$Type == "VC"])
  
  par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2, 1), cex = 1.3)
  plot(st.vc, xlim = c(0, 10),verticals=TRUE, do.points=FALSE, lwd = 2,
       xlab = "Multiple on invested capital (MOIC)", ylab = "ECDF", main = "Exit Multiple: VC deals")
  # abline(v=c(0), col = "grey", lty = 2)
  points(0, st.vc(0), lwd = 2)
  print(st.vc(0))
}
empirical.multiple()

mMMP_chart <- function(size = "big"){
  
  mpp <- data.frame(Time = c(1.564, 6.404), CashFlow = c(-0.8, 1.3))
  mpp3 <- data.frame(Time = c(3.323), CashFlow = c(-0.4))
  mpp2 <- data.frame(Time = c(4.38, 8.3434), CashFlow = c(-1.2, 0.954))
  smargin <- 0.5
  
  
  if(create.eps){
    old.wd <- getwd()

    if(size == "big"){
      setwd(wd$eps)
      setEPS()
      postscript("MPP-Visualization.eps", 
                 width = 4.5, height = 2, 
                 family = "Helvetica", pointsize = 2)
    } else {
      setwd(wd$eps2)
      setEPS()
      postscript("MPP-Visualization.eps", 
                 width = 3, height = 2, 
                 family = "Helvetica", pointsize = 2)
    }

  }
  #par(par_default)
  #par(mar = c(4.5,4.5,2,1))
  
  plot(mpp, type = "h", xlim = c(0, 10), ylim = c(-2, 2), ylab = "Cash Flow", xlab = "Time (in Quarters)", 
       cex.lab = 0.8, xaxs = "i", cex.axis = 0.8, at = seq(-10,10,1))
  points(mpp, pch = 22, cex = 2, lwd = 1)
  text(x = mpp$Time,y = c(smargin, -smargin), c(latex2exp::TeX('$\\c_1'),latex2exp::TeX('$\\d_1')), cex=1)
  text(x = mpp$Time,y = mpp$CashFlow - c(smargin, -smargin), c(latex2exp::TeX('$\\-C_1'),latex2exp::TeX('$\\D_1')), cex=1)
  
  abline(h = 0, v = 0)
  abline(v = seq(20), lty = 2, col= "darkgrey")
  
  lines(mpp2, type = "h")
  points(mpp2, pch = 21, cex = 2, lwd = 1)
  text(x = mpp2$Time,y = c(smargin, -smargin), c(latex2exp::TeX('$\\c_3'),latex2exp::TeX('$\\d_3')), cex=1)
  text(x = mpp2$Time,y = mpp2$CashFlow - c(smargin, -smargin), c(latex2exp::TeX('$\\-C_3'),latex2exp::TeX('$\\D_3')), cex=1)
  
  lines(mpp3, type = "h")
  points(mpp3, pch = 25, cex = 2, lwd = 1)
  text(x = mpp3$Time,y = smargin, c(latex2exp::TeX('$\\c_2')), cex=1)
  text(x = mpp3$Time,y = mpp3$CashFlow - c(smargin), c(latex2exp::TeX('$\\-C_2')), cex=1)
  
  
  
  if(create.eps){ 
    dev.off() 
    setwd(old.wd)
  }
}
mMMP_chart("big")
mMMP_chart("small")

mMMP_chart2 <- function(CEX = 1, meth = "reduced-form", size = "big"){
  
  mpp <- data.frame(Time = c(1.564, 6.404), CashFlow = c(-0.8, 1.3))
  mpp3 <- data.frame(Time = c(3.323), CashFlow = c(-0.4))
  mpp2 <- data.frame(Time = c(4.38, 8.3434), CashFlow = c(-1.2, 0.954))
  smargin <- 0.5
  
  
  if(create.eps){
    old.wd <- getwd()
    
    if(size == "big"){
      setwd(wd$eps)
      setEPS()
      postscript("MPP-Visualization2.eps", 
                 width = 4.5, height = 2, 
                 family = "Helvetica", pointsize = 2)
    } else {
      setwd(wd$eps2)
      setEPS()
      postscript("MPP-Visualization2.eps", 
                 width = 3, height = 2, 
                 family = "Helvetica", pointsize = 2)
    }
  }
  #par(par_default)
  #par(mar = c(4.5,4.5,2,1), cex = CEX)
  
  plot(mpp, type = "h", xlim = c(0, 10), ylim = c(-2, 2), ylab = "Cash Flow", xlab = "Time (in Quarters)", 
       cex.lab = 0.8, xaxs = "i", cex.axis = 0.8, at = seq(-10,10,1))
  points(mpp, pch = 22, cex = 2, lwd = 1)
  text(x = mpp$Time,y = c(smargin, -smargin), c(latex2exp::TeX('$\\c_1'),latex2exp::TeX('$\\d_1')), cex=1)
  text(x = mpp$Time,y = mpp$CashFlow - c(smargin, -smargin), c(latex2exp::TeX('$\\-C_1'),latex2exp::TeX('$\\D_1')), cex=1)
  
  abline(h = 0, v = 0)
  abline(v = seq(20), lty = 2, col= "darkgrey")
  
  if(meth != "reduced-form"){
    lines(mpp2, type = "h")
    points(mpp2, pch = 21, cex = 2, lwd = 1)
    text(x = mpp2$Time,y = c(smargin, -smargin), c(latex2exp::TeX('$\\c_3'),latex2exp::TeX('$\\d_3')), cex=1)
    text(x = mpp2$Time,y = mpp2$CashFlow - c(smargin, -smargin), c(latex2exp::TeX('$\\C_3'),latex2exp::TeX('$\\D_3')), cex=1)
    
    lines(mpp3, type = "h")
    points(mpp3, pch = 25, cex = 2, lwd = 1)
    text(x = mpp3$Time,y = smargin, c(latex2exp::TeX('$\\c_2')), cex=1)
    text(x = mpp3$Time,y = mpp3$CashFlow - c(smargin), c(latex2exp::TeX('$\\C_2')), cex=1)
  } else {
    valueprocess <- data.frame(Time = 4.5, CashFlow = 0.9)
    value.col <- "black"
    lines(valueprocess, type = "h", col = value.col)
    points(valueprocess, pch = 22, cex = 2, lwd = 1, col = value.col)
    text(x = valueprocess$Time,y = c(-smargin), c(latex2exp::TeX('$\\t')), cex = 1, col = value.col)
    text(x = valueprocess$Time,y = valueprocess$CashFlow + c(smargin), c(latex2exp::TeX('$\\V_1(t)')), cex = 1, col = value.col)
    
    
  }

  
  
  
  if(create.eps){ 
    dev.off() 
    setwd(old.wd)
  }
}
mMMP_chart2(size = "big")
mMMP_chart2(size = "small")

if(FALSE) {
  create.eps <- FALSE
  old.wd <- getwd()
  setwd(wd$eps)
  setEPS()
  postscript("MPP-Visualization-Combi.eps", 
             width = 5, height = 2, 
             family = "Helvetica", pointsize = 3)
  par(mfrow=c(1,2), mar = c(4.5, 4.5, 1.0, 1.0), cex=1.2)
  mMMP_chart()
  mMMP_chart2()
  dev.off() 
  setwd(old.wd)  
}


if(FALSE){
  setwd(wd$eps)
  witdth <- 250
  png("MPP1.png", witdth*2.2, witdth, pointsize = 10)
  mMMP_chart2(1.3, "other")
  dev.off()
  
  png("MPP2.png", witdth*2.2, witdth, pointsize = 10)
  mMMP_chart2(1.3)
  dev.off()
}

## 4) Simulate ------
par(par_default)
setwd(wd$code)
source("ExitDynamics_Bivariate.R")

set.seed(98)
system.time(x <- simulate$Final.Simulator(5000, N = c(3, 15), c(2, 8), "VC"))
setwd(wd$data.out)
saveRDS(x, "test_simulation_result_newPubSce.RDS")
x <- readRDS("test_simulation_result_newPubSce.RDS")


# Table: Final Multiples
quants <- c(0.005, 0.01, 0.05, 0.5)
xtable::xtable(
  data.frame(
    N3_A2 = quantile(x$Output$N_3_Age_2$Joe180, quants),
    N3_A8 = quantile(x$Output$N_3_Age_8$Joe180, quants),
    N15_A2 = quantile(x$Output$N_15_Age_2$Joe180, quants),
    N15_A8 = quantile(x$Output$N_15_Age_8$Joe180, quants)
  ), digits = 4,
  label = c("tab:final_multiples"),
  caption = c("CFaR with infinite horizon (i.e., final fund multiple) for VC funds from table XXX. 
              VC funds with more deals and longer realized holding periods are gauged safer by our model.")
)



# create EPS
if(create.eps) {
  old.wd <- getwd()
  setwd(wd$eps)
  setEPS()
  postscript("Simulation-Combi.eps", 
             width = 5, height = 3.5, 
             family = "Helvetica", pointsize = 1)
  
  par(mfrow=c(2,2), mar = c(4.5, 4.5, 1.0, 1.0), cex=1.2)
  
  df.3companies <- simulate$cf2quantiles(x$df.out$N_3_Age_2$Joe180, 3, 2)
  df.15companies <- simulate$cf2quantiles(x$df.out$N_15_Age_2$Joe180, 15, 2)
  df.3companies <- simulate$cf2quantiles(x$df.out$N_3_Age_8$Joe180, 3, 8)
  df.15companies <- simulate$cf2quantiles(x$df.out$N_15_Age_8$Joe180, 15, 8)
  
  dev.off() 
  setwd(old.wd)  
}


for(pfs in names(x$Output)){
  df.summary <- data.frame(sapply(x$Output[[pfs]], sun))
  x$CSV[[pfs]] <- df.summary
  
  if(create.csv){
    setwd(wd$data.out)
    write.csv(df.summary, file = paste("PoFo_",pfs,".csv",sep = ""), na = "", row.names = TRUE, dec = ".")
  }
}

x$Output

x$CSV
## X) Simulate bigger ------
par(par_default)
setwd(wd$code)
source("ExitDynamics_Bivariate.R")

pofo.sizes <- c(seq(1,15), 30, 60, 100, 150, 250, 500)
ages <- c(1, 5, 9)

set.seed(98)
x <- list()
for (type in c("BO", "VC")) {
  x[[type]] <- simulate$Final.Simulator(iterations = 5000, 
                                        N = pofo.sizes, 
                                        deal.ages = ages, 
                                        Type = type, 
                                        fixed.pub.sce = FALSE,
                                        # use.hist.average = TRUE,
                                        base.on.same.vintage = FALSE)
}
setwd(wd$data.out)
saveRDS(x, "test_simulation_result_newPubSce_BIG.RDS")
x <- readRDS("test_simulation_result_newPubSce_BIG.RDS")

# public scenarios for BO and VC differ
tail(x$BO$PublicScenario)
tail(x$VC$PublicScenario)

plot.quantiles <- function(fund.type = "BO", age=5, quantiles = c(0.01, 0.05, 0.1, 0.5, 0.9)) {
  l <- list()
  for (n in pofo.sizes) {
    y <- x[[fund.type]]$Output[[paste0("N_", n ,"_Age_", age)]]$Joe180
    y <- quantile(y, quantiles)
    df <- data.frame(N = n)
    for (i in 1:length(y)) {
      col <- paste("Q", names(y)[i])
      df[, col] <- y[i]
    }
    l[[paste(n)]] <- df
  }
  df <- data.frame(do.call(rbind, l))
  
  chart.name <- paste0(fund.type, " (deal age: ", age, ")")
  
  # plot
  plot(df$N, df[, 2], type="l", ylim=c(-1, 5), 
       ylab="Multiple", 
       #        main = chart.name,
       xlab="Deals in Portfolio"
       )
  for (i in 2:(length(quantiles))) {
    points(df$N, df[, i + 1], type="l", col=i)
  }
  abline(h=c(1,2), col="grey", lty=2)
  legend("bottomright", bty = "n", cex=0.5,
         legend = paste(sort(quantiles, TRUE), "quantile"), col = length(quantiles):1, lty = 1)
  
  legend("topright", bty="n", legend = chart.name)
}

# create EPS
if(create.eps) {
  old.wd <- getwd()
  setwd(wd$eps)
  setEPS()
  postscript("Simulation-Bigger-Random-Public.eps", 
             width = 5, height = 3.5, 
             family = "Helvetica", pointsize = 1)
  
  par(mfrow=c(2,2), mar = c(4.5, 4.5, 1.0, 1.0), cex=1.2)
  
  plot.quantiles("BO", age=1)
  plot.quantiles("BO", age=9)
  plot.quantiles("VC", age=1)
  plot.quantiles("VC", age=9)
  
  dev.off() 
  setwd(old.wd)  
}

## 5) Epilogue -----
sessionInfo()
