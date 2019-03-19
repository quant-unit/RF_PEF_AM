###########################
## Exit Dynamics: Main File
###########################
## 0) Prologue -------

rm(list=ls()) # remove workspace objects
setwd("C:/Users/christian.tausch")
setwd("/Users/christausch/")
root <- file.path(getwd(), "Dropbox", "Project D", "3_R_Pro_D")

wd <- list()
wd$code <- file.path(root, "Exit_Dynamics", "Code")
wd$data <- file.path(root, "Exit_Dynamics", "Code", "Data")
wd$eps <-  file.path(root, "Exit_Dynamics", "Code", "EPS_output")
wd$eps2 <-  file.path(root, "Exit_Dynamics", "Code", "EPS_output_small")
wd$csv <-  file.path(root, "Exit_Dynamics", "Code", "CSV_output")
rm(root)
dir.create(wd$eps2)

setwd(wd$code)
source("Useful_Functions.R")
source("ExitDynamics_Packages.R")

setwd(wd$data)
load("ExitDynamics_Data_V1.RData")

create.eps <- FALSE
create.csv <- FALSE
## 00) Reduce Input Data ------
g.p2c2 <- g.p2c2[g.p2c2$Fund_InvestTypes %in% c("BO", "VC"), ]
colnames(g.p2c2)

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

## 1) Timing   ------
par(par_default)
setwd(wd$code)
source("ExitDynamics_Timing.R")

# Non-Parametric Cox Regression (partial likelihood)
np.cox <- timing$NonParaCoxRegression(g.p2c2)


# Parametric Cox Regression
setwd(wd$data)
system.time(
  # wb.cox <- timing$ParaCoxRegression(public.data, g.sum)
  # saveRDS(wb.cox, "z_Para_Cox_Weibull5.rds")
  wb.cox <- readRDS("z_Para_Cox_Weibull5.rds")
)


# create .csv output
setwd(wd$csv)
timing$CreateOutput(wb.cox, create.csv)


# Create Plots
if(FALSE){
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
source("ExitDynamics_Multiple.R")
source("ExitDynamics_Copula2.R")

coef2part <- copula2part$coef2part

if (FALSE) {
  # Summary Statistics (How many observations, i.e. company IDs)
  list(First_Entry_Date = min(G.p2c$Investment_Date),
       BO_stochastic = nrow(multiple$creat_reg_df("BO", stochastic = TRUE)),
       BO_first = nrow(multiple$creat_reg_df("BO", stochastic = FALSE)),
       VC_stochastic = nrow(multiple$creat_reg_df("VC", stochastic = TRUE)),
       VC_first = nrow(multiple$creat_reg_df("VC", stochastic = FALSE)) )
  
  
  # Run Multiple Regressions
  if(FALSE){
    multi.reg <- list()
    
    for(type in c("BO","VC")){
      for(iter in c(1, 500)){
        set.seed(99)
        multi.reg[[paste(type, iter, sep="_")]] <- multiple$MLE(iter, type)
      }
    }
    
    setwd(wd$data)
    # saveRDS(multi.reg, "multi.reg2.RDS")
  } else {
    setwd(wd$data)
    multi.reg <- readRDS("multi.reg2.RDS")
  }
  
  # NB2 Model
  sapply(multi.reg$BO_500$NB2, sun)
  sapply(multi.reg$VC_500$NB2, sun)
  
  # Create .csv Output
  setwd(wd$csv)
  multi.reg.sum <- multiple$CreateOutput(multi.reg, create.csv)
  
  
  # Plot color compound CDF
  setwd(wd$eps)
  set.seed(97)
  multiple$predict_dgh(sim_out = TRUE, plot_it = TRUE, eps = create.eps, size = "big")
  setwd(wd$eps2)
  set.seed(97)
  multiple$predict_dgh(sim_out = TRUE, plot_it = TRUE, eps = create.eps, size = "small")
  
}

if(FALSE){
  
  # Plot dependent variables (Multiple)
  setwd(wd$eps)
  multiple$plot.multi.dist(create.eps)
  
  
  system.time(
    uni.CDF.multiple <- multiple$predict_dgh(sim_out = FALSE, plot_it = FALSE, eps = FALSE)
  )
  sun(uni.CDF.multiple, Ylim = c(0,2))
  
  multiple$DGH_Roseblatt(G.p2c, uni.CDF.multiple)
  
  plot.ecdf(uni.CDF.multiple)
  abline(v = c(0,1), a = 0, b = 1, col = "grey", lty = 2)
  points(seq(0,1,0.001), quantile(uni.CDF.multiple, seq(0,1,0.001)), col ="blue", pch = 20, cex = 0.7)
  
  # curve(truncdist::dtrunc(x, spec="gamma" ,a=0.1 ,b=5 ,shape=1.304 , scale= exp(1.405)/1.304), 0,6)
  
}


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
system.time(x <- simulate$Final.Simulator(5000, N = c(3, 15), "VC"))
setwd(wd$data)
# saveRDS(x, "test_simulation_result_newPubSce.RDS")
x <- readRDS("test_simulation_result_newPubSce.RDS")

tail(x$PublicScenario)
x$Output$N_3$independent
x$df.out$N_3$independent
cf2quantiles <- function(list.in, no.companies) {
  q.seq <- seq(0,1,0.01)
  df.N3 <- data.frame(do.call(rbind, list.in)) / no.companies
  df.N3 <- data.frame(apply(df.N3, 2, function(x){quantile(x, q.seq)}))
  
  if(create.eps){
    par.old <- par()
    old.wd <- getwd()
    setwd(wd$eps)
    setEPS()
    postscript(paste("Simulation_ECDF_",no.companies , "_companies.eps",sep=""), 
               width = 3, height = 2.5, family = "Helvetica", pointsize = 3)
  }
  #par(cex = 1.5, mar = c(4.5, 4.5, 1, 1))
  plot(df.N3[, paste("X", 5, sep="")], q.seq, xlim = c(0,5), type = "l",
       # main = paste(no.companies,"company VC fund"),
       xlab = "Cumulative cash flow per current value", 
       ylab = "ECDF")
  lines(df.N3[, paste("X", 10, sep="")], q.seq, lty = 2)
  abline(h = seq(0,1,0.1), col= "grey", lty = 3) ; abline(v=seq(0,10,1), col ="grey", lty = 3)
  legend("bottomright", bty = "n", legend = c("Horizon:", "5 years", "10 years"), 
         lty = c(NA, 1,2 ))
  invisible(df.N3)
  
  if(create.eps) {
    dev.off()
    setwd(old.wd)
    par(par.old)
  }
}
df.3companies <- cf2quantiles(x$df.out$N_3$independent, 3)
df.15companies <- cf2quantiles(x$df.out$N_15$independent, 15)

if(FALSE) {
  create.eps <- FALSE
  old.wd <- getwd()
  setwd(wd$eps)
  setEPS()
  postscript("Simulation-Combi.eps", 
             width = 5, height = 2, 
             family = "Helvetica", pointsize = 1)
  
  par(mfrow=c(1,2), mar = c(4.5, 4.5, 1.0, 1.0), cex=1.2)
  
  df.3companies <- cf2quantiles(x$df.out$N_3$independent, 3)
  df.15companies <- cf2quantiles(x$df.out$N_15$independent, 15)
  
  dev.off() 
  setwd(old.wd)  
}

df.3companies["10%", c("X5")]
df.15companies["10%", c("X5")]
df.3companies["90%", c("X10")]
df.15companies["90%", c("X10")]

for(pfs in names(x$Output)){
  df.summary <- data.frame(sapply(x$Output[[pfs]], sun))
  x$CSV[[pfs]] <- df.summary
  
  if(create.csv){
    setwd(wd$csv)
    write.csv(df.summary, file = paste("PoFo_",pfs,".csv",sep = ""), na = "", row.names = TRUE, dec = ".")
  }
}

x$Output

x$CSV
## 5) Epilogue -----
sessionInfo()
