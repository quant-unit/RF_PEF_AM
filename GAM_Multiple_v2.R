###########################################
## Generalized Additive Models for Multiple
###########################################
# source data ------
rm(list=ls()) # remove workspace objects

x <- unlist(strsplit(getwd(), "/"))
root <- paste(x[1:3], collapse = "/")

setwd(file.path(root, "Dropbox", "Project D", "3_R_Pro_D", "Exit_Dynamics", "Code"))
source("ExitDynamics_Data_Prep.R")

if(FALSE){
  saveRDS(g.p2c2, "g_p2c2.RDS")
  saveRDS(g.sum, "g_sum.RDS")
  saveRDS(public.data, "public_data.RDS")
  saveRDS(par_default, "par_default.RDS")
}


############# ROOT -------------
set.root <- function(){
  roo <- list()
  roo$mac <- "/Users/christausch"
  roo$win <- "C:/Users/christian.tausch"
  roo$mac.test <- try(setwd(roo$mac), silent = TRUE)
  roo$win.test <- try(setwd(roo$win), silent = TRUE)
  
  if(roo$mac.test == roo$mac){
    root <- roo$mac
    data.root <- "/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Data"
    computer <- "mac"
  } 
  if(roo$win.test == roo$win){
    root <- roo$win
    data.root <- "C:/Users/christian.tausch/iCloudDrive/Data"
    computer <- "win"
  } 
  
  print(computer)
  
  return(list(root = root, data.root = data.root, computer = computer))
}
root <- set.root()$root
data.root <- set.root()$data.root
root <- set.root()$root
data.root <- set.root()$data.root

# load data ------
# Private Data
# source("ExitDynamics_Data_Prep.R")
if(FALSE){
  g.p2c2 <- readRDS("g_p2c2.RDS")
  g.sum <- readRDS("g_sum.RDS")
  public.data <- readRDS("public_data.RDS")
  par_default <- readRDS("par_default.RDS")
}



# Fama French public Equity Market Data
setwd(file.path(data.root, "FamaFrench"))

df.fafr <- readRDS("df.famafrench.indices.RDS")
li.5fac <- readRDS("list.fivefactor.RDS")
li.USindustries <- readRDS("list_US_industries.RDS")

# Pastor Stambaugh 2003 data
setwd(file.path(data.root, "PS_2003"))
df.ps03 <- readRDS("Liq_PS2003.RDS")

# Spread Data
setwd(file.path(data.root, "Fred"))
df.spreads <- readRDS("df.spreads.RDS")


# create data.frame for historical simulation -----------
li.5fac2 <- lapply(li.5fac, function(x){
  x$Mkt <- x$Mkt.RF + x$RF
  return(x)
})

df.5fac <- data.frame(do.call(cbind, li.5fac2))
not.use <- c(grep(".Index", colnames(df.5fac)), grep(".X", colnames(df.5fac)), grep(".Date", colnames(df.5fac)))
df.5fac <- cbind(Date = df.5fac$Asia_Pacific_ex_Japan_5_Factors.Date, df.5fac[, -not.use])

df.hist <- merge(df.5fac, df.ps03[,c("Date", "Traded.Liq..LIQ_V.")], by = "Date", all.x = TRUE)
colnames(df.hist)[colnames(df.hist) == "Traded.Liq..LIQ_V."] <- "Traded.Liq.PS03"
df.hist$Traded.Liq.PS03 <- df.hist$Traded.Liq.PS03 * 100


df.industry <- li.USindustries$Value_Weighted_Returns$`10_Industry`
df.industry <- df.industry[, colnames(df.industry) != "X"]
df.industry <- df.industry[, -grep(".Index", colnames(df.industry))]
colnames(df.industry) <- paste(colnames(df.industry), "Indst", sep = ".")

df.hist <- merge(df.hist, df.industry, by.x = "Date", by.y = "Date.Indst", all.x = TRUE)

df.hist <- merge(df.hist, df.spreads[, c(1,3,4)], by = "Date", all.x = TRUE)

saveRDS(df.hist, "df_hist.RDS")

#setwd(file.path("C:", "Users", "christian.tausch", "analytics", "exit_simulation_neu", "data_input))
#write.csv(df.hist, "public_market_data.csv", row.names = FALSE)


rm(df.5fac, df.industry)
# data preparation ----------
# df.e2e <- multiple$creat_reg_df("VC", stochastic = FALSE)

G.p2c <- g.p2c2 # data set used for timing regression (load directly from ExitDynamics_Data_Prep.R)

G.p2c <- G.p2c[!is.na(G.p2c$P2C.multi1),] # JUST EXITED INVESTMENTS
regression_variables <- c("P2C.multi1","RVPI_1","Holding_Period","Time2Exit","Time2Exit_ZS","ZombieStage","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
G.p2c <- G.p2c[complete.cases(G.p2c[,regression_variables]),]

# delete young entries (bias correction?)
G.p2c <- G.p2c[G.p2c$Investment_Date < as.Date("2010-01-01"), ]

# filter by RVPI
G.p2c <- G.p2c[G.p2c$RVPI_1 > - 0.9, ]
G.p2c <- G.p2c[order(G.p2c$Fund_Emi_ID,G.p2c$FundInv_Quarter), ]


df.all <- G.p2c[G.p2c$Fund_InvestTypes %in% c("BO", "VC", "DIS"), ]
df.all <- df.all[complete.cases(df.all), ]
df.all$Fund_Emi_ID <- as.factor(df.all$Fund_Emi_ID)

table(df.all$Fund_Region)
rm(G.p2c, g.p2c2, g.sum)

# Merge with Fama French Data
public.multiple <- function(df.main = df.all, df.index = df.fafr){
  df.main$Index.now <- NA
  df.main$Index.exit <- NA
  
  for(i in 1:nrow(df.main)){
    region <- df.main$Fund_Region[i]
    if(region == "Asia") region <- "AP"
    type <- df.main$Fund_InvestTypes[i]
    date.now <- df.main$FundInv_Quarter[i]
    date.exit <- df.main$Exit_Date[i]
    
    df.main$Index.now[i] <- df.index[df.index$Date ==  date.now, paste(type, region, "Index", sep = ".")]
    df.main$Index.exit[i] <- df.index[df.index$Date ==  date.exit, paste(type, region, "Index", sep = ".")]
  }
  
  multiple <- df.main$Index.exit / df.main$Index.now
  return(multiple)
}
public.multiple2 <- function(df.main = df.all, li.index = li.5fac, df.ps = df.ps03){
  df.main[, c("SMB.Index", "HML.Index", "RMW.Index", "CMA.Index", "RF.Index", "Mkt", "PS03")] <- NA

  for(i in 1:nrow(df.main)){
    region <- df.main$Fund_Region[i]
    # type <- df.main$Fund_InvestTypes[i]
    date.now <- df.main$FundInv_Quarter[i]
    date.exit <- df.main$Exit_Date[i]
    
    if(region == "US") df.index <- li.index$North_America_5_Factors
    if(region == "Asia") df.index <- li.index$Asia_Pacific_ex_Japan_5_Factors
    if(region == "EU") df.index <- li.index$Europe_5_Factors
    
    for(cols in c("SMB.Index", "HML.Index", "RMW.Index", "CMA.Index", "RF.Index", "Mkt")){
      df.main[i, paste(cols, "Multi", sep=".")] <- df.index[df.index$Date ==  date.exit, cols] / df.index[df.index$Date ==  date.now, cols]
    }
    # PS 2003 Traded Liquidity Factor
    df.main[i, paste("PS03", "Multi", sep=".")] <- df.ps$Traded.Liq..LIQ_V..Index[df.ps$Date == date.exit] / df.ps$Traded.Liq..LIQ_V..Index[df.ps$Date == date.now]
    
  }
  
  df.out <- df.main[, paste(c("SMB.Index", "HML.Index", "RMW.Index", "CMA.Index", "RF.Index", "Mkt", "PS03"), "Multi", sep=".")]
  return(df.out)
}
public.multiple3 <- function(df.main = df.all, li.index = li.5fac, li.Industry = li.USindustries){
  # df.main = df.all ; li.index = li.5fac ; li.Industry = li.USindustries
  df.market <- li.index$North_America_5_Factors
  # df.index <- li.USindustries$Equal_Weighted_Returns$`10_Industry`
  df.index <- li.USindustries$Value_Weighted_Returns$`10_Industry`
  
  
  df.IndustryMapping <- data.frame(GICS  = levels(df.main$GICS_Sector), 
             Index = paste(c("Durbl", "NoDur", "Enrgy", "Other", "Hlth", "Manuf", "HiTec", "Other", "Other", "Other", "Telcm", "Utils"),
                           "Index", sep = "."))
  print(df.IndustryMapping)
  
  df.main$Index.now <- NA
  df.main$Index.exit <- NA
  df.main$Market.now <- NA
  df.main$Market.exit <- NA
  
  for(i in 1:nrow(df.main)){
    industry <- df.main$GICS_Sector[i]
    index.col <- df.IndustryMapping[df.IndustryMapping$GICS == industry, "Index"]
    index.col <- as.character(index.col)
    
    date.now <- df.main$FundInv_Quarter[i]
    date.exit <- df.main$Exit_Date[i]
    
    df.main$Index.now[i] <- df.index[df.index$Date ==  date.now, index.col]
    df.main$Index.exit[i] <- df.index[df.index$Date ==  date.exit, index.col]
    df.main$Market.now[i] <- df.market[df.market$Date == date.now, "Mkt"]
    df.main$Market.exit[i] <- df.market[df.market$Date == date.exit, "Mkt"]
  }
  
  multiple.index <- df.main$Index.exit / df.main$Index.now
  multiple.market <- df.main$Market.exit / df.main$Market.now
  spread <- multiple.index - multiple.market

  return(spread)
}


system.time(df.multis <- public.multiple2(df.all, li.5fac, df.ps03))
df.all <- cbind(df.all, df.multis)
# system.time(df.all$FaFr.Multi <- public.multiple(df.all, df.fafr))
system.time(df.all$IndustryMinusMarket <- public.multiple3())



# merge spreads
df.all1 <- merge(df.all, df.spreads, by.x = "FundInv_Quarter", by.y = "Date", all.x = TRUE)

df.all1$HYS_mixed <- ifelse(df.all1$Fund_Region == "Asia", 
                            (df.all1$HYS_US + df.all1$HYS_EU) / 2,
                            ifelse(df.all1$Fund_Region == "US",
                                   df.all1$HYS_US, df.all1$HYS_EU))
df.all1$HYS_mixed <- df.all1$HYS_mixed / 100

# saveRDS(df.all1, "df_all1.RDS")
# read public index return ------------
library(readxl)
path.data <- file.path(data.root, "PublicPrivate")
setwd(path.data)
file.name <- "2018-10-08 Update-2d1882-91cc0a.xlsx"

df <- list()
for(sheet in c("ca_index_100", "usd2fx", "public_indices")) {
  df.data <- read_xlsx(file.name, sheet = sheet)
  df.data <- as.data.frame(df.data)
  quarterly.subset <- sapply(df.data$Date, function(x){
    substr(x, 6, 7) %in% c("03", "06", "09", "12")
  })
  df.data <- df.data[quarterly.subset, ]
  df.data$Date <- as.Date(df.data$Date)
  
  if(sheet != "ca_index_100") {
    y <- apply(df.data[, colnames(df.data) != "Date"], 2, function(x) {
      x <- as.numeric(x)
      c(diff(x), NA) / x
    })
    dates <- df.data$Date[-1]
    df.data <- data.frame(y[-nrow(y),])
    df.data$Date <- dates
    rm(y)
  }
  
  df[[sheet]] <- df.data
  rm(df.data)
}
fill.dates <- as.Date(df$public_indices$Date[df$public_indices$Date < min(df$ca_index_100$Date) | df$public_indices$Date > max(df$ca_index_100$Date)])
df.ca100x <- df$ca_index_100[1:length(fill.dates), ]
df.ca100x$Date <- fill.dates
ncol.ca <- ncol(df.ca100x)
for(i in 1:length(fill.dates)) df.ca100x[i, 2:ncol.ca] <- rep(0, (ncol.ca-1))
df$ca_index_100 <- rbind(df$ca_index_100, df.ca100x)
df$ca_index_100$RE <- df$ca_index_100$RE_FUNDS
df$ca_index_100$Infrastructure <- df$ca_index_100$INF_FUNDS
df$ca_index_100$NATRES <- df$ca_index_100$NATRES_FUNDS
df$ca_index_100$Date <- as.Date(df$ca_index_100$Date)
df$ca_index_100 <- df$ca_index_100[order(df$ca_index_100$Date), ]

# View(df$public_indices)
# View(df$usd2fx)
# View(df$ca_index_100)
colnames(df$public_indices)

to_return <- function(vector) {
  x <- exp(cumsum(log(1 + vector)))
  return(x)
}
return_from_to <- function(am_type, am_region, currency, start.date, end.date) {
  start.date <- as.Date(start.date)
  end.date <- as.Date(end.date)
  
  # Public Indices
  public.dates <- (df$public_indices$Date >= start.date) & (df$public_indices$Date <= end.date)
  
  msci.world.multiple <- to_return(df$public_indices[public.dates, "MSCI.World.Net.Return.Daily"])
  if(is.na(am_region) | am_region %in% c("US", "NA")) {
    msci.region.multiple <- to_return(df$public_indices[public.dates, "MSCI.North.America.Net.Return.Daily"])
  }
  if(am_region %in% c("AP", "Asia", "Australasia")) {
    msci.region.multiple <- to_return(df$public_indices[public.dates, "MSCI.AC.Asia.Pacific.Net.Return.Daily"])
  }
  if(am_region %in% c("EU","Europe")) {
    msci.region.multiple <- to_return(df$public_indices[public.dates, "MSCI.Europe.Net.Return.Daily"])
  }
  if(am_region %in% c("ROW","Africa","Americas","Diversified Multi-Regional", "Middle East & Israel")) {
    msci.region.multiple <- msci.world.multiple
  }
  if(am_region %in% c("GLO")) {
    msci.region.multiple <- msci.world.multiple
  }
  russell.2000.multiple <- to_return(df$public_indices[public.dates, "Russell.2000.Total.Return.Index"])
  nasdaq.100.multiple <- to_return(df$public_indices[public.dates, "NASDAQ.100.Total.Return"])
  sovg.bonds.multiple <- to_return(df$public_indices[public.dates, "S.P.Global.Developed.Sovereign.Bond.Index.Total.Return"])
  
  # Preqin
  private.dates <- (df$ca_index_100$Date >= start.date) & (df$ca_index_100$Date <= end.date)
  private.returns <- df$ca_index_100[private.dates, am_type]
  private.multiple <- to_return(private.returns)
  
  if(currency == "USD") {
    currency.returns <- rep(0, length(private.returns))
  } else {
    currency.dates <- (df$usd2fx$Date >= start.date) & (df$usd2fx$Date <= end.date)
    currency.returns <- df$usd2fx[currency.dates, paste(currency, "USD", sep = ".")]
  }
  currency.multiple <- to_return(currency.returns)
  
  return(list(Transaction.Quarter = df$public_indices$Date[public.dates],
              msci.world.multiple = msci.world.multiple - sovg.bonds.multiple,
              msci.region.multiple = msci.region.multiple - msci.world.multiple,
              russell.2000.multiple = russell.2000.multiple - msci.world.multiple,
              nasdaq.100.multiple = nasdaq.100.multiple - msci.world.multiple,
              sovg.bonds.multiple = sovg.bonds.multiple,
              private.multiple = private.multiple - msci.world.multiple,
              currency.multiple = currency.multiple))
}
return_from_to("BO", "NA", "USD", "2001-03-31", "2018-09-30")


extract.CA.factor <- function(am_type, start.date, end.date) {
  start.date <- as.Date(start.date)
  end.date <- as.Date(end.date)
  
  if(am_type == "DIS") am_type <- "DD"
  
  # CA
  private.dates <- (df$ca_index_100$Date >= start.date) & (df$ca_index_100$Date <= end.date)
  private.returns <- df$ca_index_100[private.dates, am_type]
  private.multiple <- exp(sum(log(1 + private.returns)))
  return(private.multiple)
}

ca <- apply(df.all1, 1, function(x){
  extract.CA.factor(x[["Fund_InvestTypes"]], x[["FundInv_Quarter"]], x[["Exit_Date"]])
})

df.all1$CA_private <- ca

# saveRDS(df.all1, "df_all1.RDS")


# load packages, util functions, define truncated family -----
library(gamlss)
library(gamlss.tr)

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}



# Create Truncated Generalized Gamma Family
gamlss.tr::gen.trun(par = 0.1, family = GG, name = "tr", type = "left")

gamlss.tr::gen.trun
gamlss.tr::trun.q

# Fit Truncated GGtr by stepGAIC(k) ---------
# https://rdrr.io/cran/gamlss/man/stepGAIC.html

# Fixed Parameters
nC <- detectCores() # no of cores for parallel stepGAICAll.A
default.threshold <- 0.1

create.FORM <- function(fund.type = "undefined"){
  predictors <- c("Mkt.Multi", 
                  # "IndustryMinusMarket", # problematic for VC
                  "I(Mkt.Multi - 1)", 
                  "I(SMB.Index.Multi - 1)",                       
                  "I(CMA.Index.Multi - 1)", 
                  "I(RMW.Index.Multi - 1)", 
                  "I(PS03.Multi - 1)",
                  # "Holding_Period",
                  "HYS_mixed", 
                  "Time2Exit_ZS", 
                  "RVPI_1")
  if(fund.type != "VC"){
    predictors <- c(predictors, "IndustryMinusMarket")
  }
  linear.terms <- paste(predictors, collapse = "+")
  if(FALSE){
    local.pb.spline.method <- "ML"
    pb.terms <- paste(paste("pb(", predictors,", method='",local.pb.spline.method, "')",sep=""), collapse = "+")
  } else {
    df <- 2
    pb.terms <- paste(paste("pb(", predictors,", df=",df, ")",sep=""), collapse = "+")
  }
  all.terms <- paste(linear.terms, "+", pb.terms)
  FORM = as.formula(paste("~", all.terms, sep = ""))
  
  return(FORM)
}



# For Loop 0
reg.family <- "GGtr"
# For Loop 1
fund.type <- "BO"
# For Loop 2
past.holding.period <- 0
# For Loop


li <- list()
fund.typeZ <- c("BO", "VC" ,"DIS")
fund.typeZ <- "VC"
holding.periodZ <- c(0,5,10)
k <- list(AIC = 2, BIC = NA, SqrtBIC = NA)
k["BIC"] <- NULL
k["AIC"] <- NULL
# k["SqrtBIC"] <- NULL
system.time(
for(fund.type in fund.typeZ){
  for(holding.period in holding.periodZ){
    reg.data <- df.all1[df.all1$Fund_InvestTypes == fund.type & 
                          df.all1$Holding_Period == holding.period &
                          df.all1$P2C.multi1 > default.threshold, ]
  no.funds <- length(levels(as.factor(reg.data$FundInv_FundId)))
  no.obs <- nrow(reg.data)
  k["BIC"] <- log(nrow(reg.data))
  k["SqrtBIC"] <- sqrt(log(nrow(reg.data)))
    for(k2use in names(k)){
      print(paste(fund.type, holding.period, k2use))

      # Start with Intercept Only Model
      m1 <- try(gamlss::gamlss(P2C.multi1 ~ 1, family = GGtr, data = reg.data))
      m2 <- try(gamlss::stepGAICAll.A(m1,
                                  k = unlist(k[k2use]),
                                  parallel = "multicore", ncpus = nC,
                                  scope = list(lower= ~ 1, upper = create.FORM(fund.type))))
      drop.bic <- try(dropterm(m2, k = unlist(k["BIC"]), parallel="multicore",  ncpus=nC))
      li[[fund.type]][[as.character(holding.period)]][[k2use]] <- list(m1 = m1, m2 = m2, drop.bic = drop.bic, no.funds = no.funds, no.obs = no.obs)
    }
  }
}
)

model <- li$VC$`0`$SqrtBIC$m2
summary(model)
Rsq(model)
term.plot(model, "mu")
term.plot(model, "sigma")


# gamlss::term.plot(li$DIS$`0`$SqrtBIC$m2, "mu")

getwd()
#setwd(file.path(root, "Dropbox", "Project D", "3_R_Pro_D", "Exit_Dynamics", "Code", "Data"))
#saveRDS(li, "GAMLSS_list.RDS")
# li2 <- readRDS("GAMLSS_list.RDS")


# VC
# sigma:  P2C.multi1 ~ pb(Time2Exit_ZS, method = "ML") + HYS_mixed + I(RMW.Index.Multi -  1) 
# mu: P2C.multi1 ~ pb(I(Mkt.Multi - 1), method = "ML") + pb(I(PS03.Multi -  1), method = "ML") + Time2Exit_ZS 

#summary(m1)
#summary(m2)
#term.plot(m2, "mu")
#term.plot(m2, "sigma")






#######################################
## OLD APPROACHES
#################
# Fit Truncated Generalized Gamma for BO and VC ---------
# Truncated Generalized Gamma Model (nu = 0 --> log-normal)
gam_models <- list()

# Buy Out: Fama French Factors
system.time( gam_models$BO$GGtr <- gamlss(P2C.multi1 ~
    # MSCI.Multiple.Exit_1
    I(Mkt.Multi - 1) + IndustryMinusMarket + I(SMB.Index.Multi - 1) + I(CMA.Index.Multi - 1) + I(RMW.Index.Multi - 1) + HYS_mixed + pb(Time2Exit_ZS),
    sigma.formula = ~ RVPI_1  + I(PS03.Multi - 1) + pb(Time2Exit_ZS) + (Holding_Period) + I(SMB.Index.Multi - 1) + I(CMA.Index.Multi - 1) + I(RMW.Index.Multi - 1),
    nu.formula = ~ Mkt.Multi + Holding_Period + 0, # explainable skewness (intercept == 0)
    family = GGtr,
    data = df.all1[df.all1$Fund_InvestTypes == "BO" & df.all1$P2C.multi1 > 0.1, ]) )
summary(gam_models$BO$GGtr)


'
# Buy Out: Industry and Tailored FamaFrench Portfolios
system.time( gam_models$BO$GGtr2 <- gamlss(P2C.multi1 ~
        I(FaFr.Multi - 1) + IndustryMinusMarket + HYS_mixed + pb(Time2Exit_ZS),
        # sigma.formula = ~ HYS_mixed + pb(Time2Exit_ZS),
        # nu.formula = ~ Mkt.Multi + Holding_Period + 0, # explainable skewness (intercept == 0)
        family = GGtr,
        data = df.all1[df.all1$Fund_InvestTypes == "BO" & df.all1$P2C.multi1 > 0.1, ]) )
summary(gam_models$BO$GGtr2)
'

# Venture Capital: Fama French Factors
system.time( gam_models$VC$GGtr <- gamlss(P2C.multi1 ~
    # MSCI.Multiple.Exit_1
    I(Mkt.Multi - 1) + I(PS03.Multi - 1) + I(SMB.Index.Multi - 1) + I(CMA.Index.Multi - 1) + I(RMW.Index.Multi - 1) + HYS_mixed + Time2Exit_ZS,
    sigma.formula = ~ I(SMB.Index.Multi - 1) + I(CMA.Index.Multi - 1) + I(RMW.Index.Multi - 1)
    + (RVPI_1) + pb(Time2Exit_ZS) + Holding_Period + IndustryMinusMarket,
    nu.formula = ~ Mkt.Multi, # unexplainable skewness (intercept != 0)
    family = GGtr,
    data = df.all1[df.all1$Fund_InvestTypes == "VC" & df.all1$P2C.multi1 > 0.1, ]) )
summary(gam_models$VC$GGtr)

# DD: Fama French Factors
system.time( gam_models$DIS$GGtr <- gamlss(P2C.multi1 ~
       # MSCI.Multiple.Exit_1
     0 + I(Mkt.Multi - 1) + IndustryMinusMarket + I(PS03.Multi - 1) + I(RMW.Index.Multi - 1),
     # + HYS_mixed + I(SMB.Index.Multi - 1) + Time2Exit_ZS,
     sigma.formula = ~ 0 + I(RMW.Index.Multi - 1)
     + Holding_Period + IndustryMinusMarket,
     nu.formula = ~ HYS_mixed + 0, # unexplainable skewness (intercept != 0)
     family = GGtr,
     data = df.all1[df.all1$Fund_InvestTypes == "DIS" & df.all1$P2C.multi1 > 0.1, ]) )
summary(gam_models$DIS$GGtr)


# All: Fama French Factors
system.time( gam_models$ALL$GGtr <- gamlss(P2C.multi1 ~
      # MSCI.Multiple.Exit_1
      I(Mkt.Multi - 1) + I(SMB.Index.Multi - 1) + I(CMA.Index.Multi - 1) + I(RMW.Index.Multi - 1) + HYS_mixed + Time2Exit_ZS,
    sigma.formula = ~ I(PS03.Multi - 1) + I(CMA.Index.Multi - 1) + I(RMW.Index.Multi - 1)
    + (RVPI_1) + pb(Time2Exit_ZS) + Holding_Period + IndustryMinusMarket,
    nu.formula = ~ Mkt.Multi + I(PS03.Multi - 1), # unexplainable skewness (intercept != 0)
    family = GGtr,
    data = df.all1[df.all1$P2C.multi1 > 0.1, ]) )
summary(gam_models$ALL$GGtr)


# select BO or VC model
ggtr.model <- gam_models$VC$GGtr
ggtr.model <- gam_models$BO$GGtr

summary(ggtr.model)
Rsq(ggtr.model, "both")
#term.plot(ggtr.model, "mu")
#term.plot(ggtr.model, "sigma")
#term.plot(ggtr.model, "nu")

ggtr.model$mu.coefficients

PEM_1 <- 0.5
SMB_1 <- CMA_1 <- RMW_1 <- PS03 <- 0
HYS <- 0.05
T2E_ZS <- 5
HP <- 5
RVPI_1 <- 0

mu <- exp(sum(ggtr.model$mu.coefficients * c(1, PEM_1, SMB_1, CMA_1, RMW_1,  HYS, T2E_ZS)))
si <- exp(sum(ggtr.model$sigma.coefficients * c(1, RVPI_1, PS03, T2E_ZS, HP, SMB_1, CMA_1, RMW_1)))
nu <- sum(ggtr.model$nu.coefficients * c(1, HP))
print(c(Mu = mu, Sigma = si, NU = nu))

curve(dGGtr(x, mu = mu, sigma = si, nu = nu), from = 0.1, to = 10, col = "red")

curve(qGGtr(x, mu = mu, sigma = si, nu = nu), from = 0.0001, to = 0.99, col = "red")
curve(qGGtr(x, mu = mu, sigma = 1.3793288, nu = nu), from = 0.0001, to = 0.99, add = TRUE)
abline(h = 0.1) ; abline(h = 1:5, col = "grey", lty = 3) ; abline(v = seq(0, 1, 0.1), col = "grey", lty = 3)


'
# Truncated Generalized Gammma (Random Effect)
system.time( gamlss.h1b <- gamlss(P2C.multi1 ~ 
          (MSCI.Multiple.Exit_1) + pb(ML_HYOAS.quarter) + pb(Time2Exit_ZS) + 
          random(Fund_Emi_ID, 2),
          # re(random = ~ 1 | Fund_Emi_ID),
          nu.formula = ~  MSCI.Multiple.Exit_1,
          sigma.formula = ~ pb(RVPI_1) + (ML_HYOAS.quarter) + pb(Time2Exit_ZS) + pb(Holding_Period),
          family = GGtr,
          data = df.all[df.all$P2C.multi1 > 0.1, ]))
summary(gamlss.h1b)
Rsq(gamlss.h1b, "both")
term.plot(gamlss.h1b, "mu")
term.plot(gamlss.h1b, "sigma")
sun(gamlss.h1b$mu.coefSmo[[1]]$coef) # random coef (add to intercept term)
'

# Fit Binomial Logit Hurdle for BO and VC --------
# (i.e. logistic regression)

# Buy Out
system.time( gam_models$BO$BI <- gamlss(I(P2C.multi1 > 0.1) ~ 
                                   # random(Fund_Emi_ID) + 
                                   # pb(MSCI.Multiple.Exit_1) +
                                   I(Mkt.Multi - 1) +  I(CMA.Index.Multi - 1) + I(RMW.Index.Multi - 1) + 
                                   IndustryMinusMarket + # I(SMB.Index.Multi - 1) +
                                   HYS_mixed +
                                   RVPI_1 +
                                   # Time2Exit_ZS + # not significant for BO
                                   Holding_Period,
                      family = BI(mu.link = "logit"), 
                      data = df.all1[df.all1$Fund_InvestTypes == "BO", ]) )
summary(gam_models$BO$BI)
gamlss::Rsq(gam_models$BO$BI)

system.time( gam_models$VC$BI <- gamlss(I(P2C.multi1 > 0.1) ~ 
                                  # random(Fund_Emi_ID) + 
                                  # pb(MSCI.Multiple.Exit_1) +
                                  I(Mkt.Multi - 1) + I(PS03.Multi - 1) +  I(CMA.Index.Multi - 1) + I(RMW.Index.Multi - 1) + 
                                  IndustryMinusMarket + # I(SMB.Index.Multi - 1) +
                                  HYS_mixed +
                                  RVPI_1 +
                                  Time2Exit_ZS +
                                  Holding_Period,
                                family = BI(mu.link = "logit"), 
                                data = df.all1[df.all1$Fund_InvestTypes == "VC", ]) )
summary(gam_models$VC$BI)

system.time( gam_models$DIS$BI <- gamlss(I(P2C.multi1 > 0.1) ~ 
                            # random(Fund_Emi_ID) + 
                            # pb(MSCI.Multiple.Exit_1) +
                            I(Mkt.Multi - 1) + I(PS03.Multi - 1) + I(RMW.Index.Multi - 1) + 
                            # IndustryMinusMarket + I(SMB.Index.Multi - 1) + I(CMA.Index.Multi - 1) +
                            # HYS_mixed + Holding_Period
                            RVPI_1 +
                            Time2Exit_ZS,
                          family = BI(mu.link = "logit"), 
                          data = df.all1[df.all1$Fund_InvestTypes == "DIS", ]) )
summary(gam_models$DIS$BI)

system.time( gam_models$ALL$BI <- gamlss(I(P2C.multi1 > 0.1) ~ 
                                          # random(Fund_Emi_ID) + 
                                          # pb(MSCI.Multiple.Exit_1) + I(PS03.Multi - 1) +
                                          I(Mkt.Multi - 1)  + I(CMA.Index.Multi - 1) + I(RMW.Index.Multi - 1) + 
                                          IndustryMinusMarket + I(SMB.Index.Multi - 1) +
                                          HYS_mixed +
                                          RVPI_1 +
                                          Time2Exit_ZS +
                                          Holding_Period,
                                        family = BI(mu.link = "logit"), 
                                        data = df.all1) )
summary(gam_models$ALL$BI)

# select logistic model
bi.model <- gam_models$VC$BI
bi.model <- gam_models$BO$BI

summary(bi.model)
Rsq(bi.model, "both")
# term.plot(bi.model, "mu")
# sun(bi.model$mu.coefSmo[[1]]$coef)

logit2prob(1.386691 + 0.38)



'
system.time( gamlss.h0b <- gamlss(I(P2C.multi1 > 0.01) ~ re(
                                    random = ~ 1 | Fund_Emi_ID,
                                    correlation = corAR1()
                                    ) +
                                    # pb(MSCI.Multiple.Exit_1) +
                                    pb(ML_HYOAS.quarter) +
                                    pb(RVPI_1) +
                                    pb(Time2Exit_ZS) +
                                    pb(Holding_Period),
                                  family = BI(mu.link = "logit"), 
                                  data = df.all) )

?nlme::lme

summary(gamlss.h0b)
Rsq(gamlss.h0b, "both")
term.plot(gamlss.h0b, "mu")
str(gamlss.h0)
head(gamlss.h0$mu.s, 25)
sun(gamlss.h0$mu.coefSmo[[1]]$coef) # random coef (add to intercept term)
'

# Summary --------
to.df <- function(x, col.name = "Coef"){
  y <- data.frame(x)
  colnames(y) <- col.name
  y$Variable <- rownames(y)
  y
}

df.sum <- merge(
  merge(
    merge(
      to.df(gam_models$BO$BI$mu.coefficients, "BO.BI.mu"),
      to.df(gam_models$BO$GGtr$mu.coefficients, "BO.GG.mu"),
      by = "Variable", all = TRUE
    ),
    merge(
      to.df(gam_models$VC$BI$mu.coefficients, "VC.BI.mu"),
      to.df(gam_models$VC$GGtr$mu.coefficients, "VC.GG.mu"),
      by = "Variable", all = TRUE
    ),
    by = "Variable", all = TRUE
  ),
  merge(
    merge(
      to.df(gam_models$BO$GGtr$sigma.coefficients, "BO.GG.sigma"),
      to.df(gam_models$BO$GGtr$nu.coefficients, "BO.GG.nu"),
      by = "Variable", all = TRUE
    ),
    merge(
      to.df(gam_models$VC$GGtr$sigma.coefficients, "VC.GG.sigma"),
      to.df(gam_models$VC$GGtr$nu.coefficients, "VC.GG.nu"),
      by = "Variable", all = TRUE
    ),
    by = "Variable", all = TRUE
  ),
  by = "Variable", all = TRUE
)

View(df.sum)

make_summary <- function() {
  
  list.out <- list()
  for(type in names(gam_models)){
    list.type <- gam_models[[type]]
    for(model in names(list.type)){
      list.model <- gam_models[[type]][[model]]
      if(model == "GGtr") {
        coef.names <- c("mu.coefficients", "sigma.coefficients", "nu.coefficients")
      }
      if(model == "BI") {
        coef.names <- "mu.coefficients"
      }
      for(coef in coef.names){
        para <- gam_models[[type]][[model]][[coef]]
        para <- data.frame(para)
        colnames(para) <- "Coef"
        para$VariableName <- rownames(para)
        coef.new <- strsplit(coef, "[.]")[[1]][1]
        type.new <- ifelse(type == "DIS", "DD", type)
        model.new <- ifelse(model == "GGtr", "GG", model)
        para$ID <- paste(type.new, model.new, coef.new, sep = "_")
        list.out[[paste(type, model, coef, sep = "_")]] <- para
      }
    }
  }
  df.out <- do.call(rbind, list.out)
  rownames(df.out) <- NULL
  df.out <- reshape(data = df.out, v.names = "Coef", direction = "wide", timevar = "VariableName", idvar = "ID")
  names(df.out) <- gsub("Coef.", "", names(df.out))
  df.out[is.na(df.out)] <- 0
  
  df.out <- df.out[, c("ID", sort(colnames(df.out[2:ncol(df.out)])))]
  
  return(df.out)
}
df.sum2 <- make_summary()
View(df.sum2)
strsplit("mu.coefficients", "[.]")[[1]][1]


#setwd(file.path("C:", "Users", "christian.tausch", "analytics", "exit_simulation_neu", "data_input"))
#write.csv(df.sum2, "GG_regression_params2.csv", row.names = FALSE)


