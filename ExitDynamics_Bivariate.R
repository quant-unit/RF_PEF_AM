###########################
## Exit Dynamics: Bivariate
###########################
## x) Create new environment  ------
simulate <- new.env()

## A) Input Data ------
simulate$create.sim.input <- function(no_companies = 10, type = "VC", 
                                      today = "2016-12-31", do.random = FALSE){
  if(do.random){
    df <- data.frame(CompanyAge =  rpois(no_companies, 3),
                     Type = rep(type, no_companies),
                     TVPI = rlnorm(no_companies,0,0.5))
    df$RVPI <- df$TVPI * runif(no_companies, 0.5, 1)
  }else{
    df <- data.frame(CompanyAge = rep(5, no_companies), # rpois(no_companies, 3),
                     Type = rep(type, no_companies),
                     TVPI = rep(1.5, no_companies))
    df$RVPI <- df$TVPI - 0.5
  }
  
  ApproxCompanyStartDate <- as.Date(today) - 365 * df$CompanyAge
  pos <- sapply(ApproxCompanyStartDate, function(y){
    which.min(abs(y - public.data$Date))
  })  
  df$CompanyStartDate <- public.data$Date[pos]
  
  df$FundAge <- max(df$CompanyAge)
  df$FundAgeAtEntry <- df$FundAge - df$CompanyAge
  return(df)
}
simulate$create.sim.input()

## B) Public Scenario  -----
setwd(wd$data)
# df.hist <- readRDS("df_hist.RDS")
set.seed(99)
public.data$Global.CMA <- rnorm(nrow(public.data),0,0.02)
simulate$create.public.scenario <- function(periods = 12 * 20, upto = "2016-12-31"){
  # empirical past
  upto <- as.Date(upto)
  indices <- c("Date", "ML_HYOAS", "MSCI_monthly_return", "Global.CMA")
  em.past <- public.data[public.data$Date <= upto, indices]
  em.past$MSCI.Multiple.Exit_1 <- 0
  em.past$ML_HYOAS.quarter <- 0.05
  em.past$CMA.Multiple.Exit_1 <- 0
  
  # future scenario
  public.scenario <- public.data[complete.cases(public.data), indices]

  public.scenario <- public.scenario[sample(nrow(public.scenario), periods, replace = TRUE), ]
  public.scenario$Date <- seq.Date(as.Date(upto) + 1, by = "month", length.out = periods) - 1
  
  public.scenario$MSCI.Multiple.Exit_1 <- exp(cumsum(log(1 + public.scenario$MSCI_monthly_return))) - 1
  public.scenario$CMA.Multiple.Exit_1 <- exp(cumsum(log(1 + public.scenario$Global.CMA))) - 1
  
  public.scenario$ML_HYOAS.quarter <- public.scenario$ML_HYOAS[1]
  
  # combine past and future
  out <- rbind(em.past, public.scenario[-1, ])

  return(out)
}
tail(simulate$create.public.scenario())



## 0) Bivariate / Dependence / Copula -------
simulate$TM.copula <- function(n, cluster, joe.para = 1.10){
  myCop.clayton <- copula::archmCopula(family = "clayton", dim = 2, param = 1)
  U.bivariate <- data.frame(copula::rCopula(n, myCop.clayton))
  colnames(U.bivariate) <- c("Timing", "Multiple")
  
  if(cluster == "bottom-left"){
    out <- U.bivariate
  }
  if(cluster == "top-left"){
    U.bivariate$Multiple <- 1 - U.bivariate$Multiple
    out <- U.bivariate
  }
  if(cluster == "top-right"){
    U.bivariate$Multiple <- 1 - U.bivariate$Multiple
    U.bivariate$Timing <- 1 - U.bivariate$Timing
    out <- U.bivariate
  }
  if(cluster == "bottom-right"){
    U.bivariate$Timing <- 1 - U.bivariate$Timing
    out <- U.bivariate
  }
  if(cluster == "independent"){
    out <- data.frame(Multiple = runif(n),
                      Timing = runif(n))
  }
  if(cluster == "t.copula"){
    tcop <- copula::ellipCopula(family = "t", dim = 2, dispstr = "toep", param = 0, df = 1)
    U.bivariate <- data.frame(copula::rCopula(n, tcop))
    colnames(U.bivariate) <- c("Timing", "Multiple")
    out <- U.bivariate
  }
  if(cluster == "Joe180") {
    U.bivariate <- data.frame(VineCopula::BiCopSim(n, par = joe.para, family = 16))
    colnames(U.bivariate) <- c("Timing", "Multiple")
    out <- U.bivariate
  }
  
  # plot(out, xlab = "Timing", ylab = "Multiple")
  return(out)
}

simulate$TM.copula(1000, "Joe180", sample(c(2,4,3), 1000, replace = TRUE))

# TODO: we need rotated copulas as "stress scenarios"
tcop <- copula::ellipCopula(family = "t", dim = 2, dispstr = "toep", param = 0, df = 1)
plot(copula::rCopula(2000, tcop))

# 180-degree rotated Joe Copula
plot(VineCopula::BiCopSim(1000, par = 1.10, family = 16))

## 1) Cox Weibull Simulation Function -----
simulate$Simulate.WBcox <- function(UniForm, # from Copula
                           FundAgeEntry = 5, # at Entry
                           CompanyStartDate = "2016-10-31", 
                           today = "2016-12-31",
                           Type = "BO",
                           para.in = wb.cox,
                           df.timing){
  if(Type == "BO"){
    para_fit = para.in$BO$Full$par
  }
  if(Type == "VC"){
    para_fit = para.in$VC$Full$par
  }
  
  periods <- nrow(df.timing)

  # Betas
  beta_MSCI = para_fit["MSCI"]
  beta_HYS = para_fit["HYS"] 
  beta_FundAge = para_fit["FundAge"]
  scale_wb = para_fit["Scale"]
  shape_wb = para_fit["Shape"]
  
  today <- as.Date(today)
  CompanyStartDate <- as.Date(CompanyStartDate)
  df.timing <- df.timing[df.timing$Date >= CompanyStartDate, ]
  
  # df.timing$Time2Exit
  df.timing$d <- seq(1, nrow(df.timing)) / 12 
  df.timing$Time2Exit <- pmax(0, df.timing$d - df.timing$d[df.timing$Date == today])
  df.timing$FundAgeEntry <- FundAgeEntry
  
  # df.timing$Date <- seq.Date(from = as.Date("2016-12-31") + 1, length.out = periods, by = "month") - 1
  
  df.timing$base_haze_rate_wb_exact <- c(0,diff( (df.timing$d / scale_wb)^shape_wb ))
  
  df.timing$expBX <- exp( as.matrix(df.timing[,c("MSCI_monthly_return","ML_HYOAS","FundAgeEntry")]) %*% 
                              c(beta_MSCI, beta_HYS, beta_FundAge) )
  
  df.timing$Cum_Haze_exact <- cumsum(as.numeric(df.timing$base_haze_rate_wb_exact) * df.timing$expBX)
  df.timing$Surv_WB <- exp(-df.timing$Cum_Haze_exact)
  
  Future.Uniform <- UniForm * df.timing$Surv_WB[df.timing$Date == today]
  sim_id <- which(abs(df.timing$Surv_WB - Future.Uniform) == min(abs(df.timing$Surv_WB - Future.Uniform)))
  out <- df.timing[sim_id, ]
  
  # plot(df.timing$d, df.timing$Surv_WB) ; abline(v = out$d, col = "red")
  # print(head(df.timing, 10))
  return(out)
}

simulate$Simulate.WBcox(UniForm = 0.99, df.timing = simulate$create.public.scenario())


simulate$Timing.Simulator <- function(df.pr, df.pu, U.timing = runif(nrow(df.pr))){
  timing.list <- list()
  for(i in seq(nrow(df.pr))){
    Timing.sim <- simulate$Simulate.WBcox(
                   UniForm = U.timing[i],
                   FundAgeEntry = df.pr$FundAgeAtEntry[i], 
                   CompanyStartDate = df.pr$CompanyStartDate[i],
                   Type = df.pr$Type[i],
                   df.timing = df.pu)
    timing.list[[i]] <- Timing.sim[,c("MSCI.Multiple.Exit_1", "CMA.Multiple.Exit_1", "ML_HYOAS.quarter", "d", "Time2Exit", "Date")]
  }
  timing.list <- data.frame(do.call(rbind, timing.list))
  df.pr <- cbind(df.pr, timing.list)
  
  # rename for multiple predicion
  df.pr$Fund_InvestTypes <- df.pr$Type
  df.pr$RVPI_1 <- df.pr$RVPI - 1
  df.pr$Holding_Period <- df.pr$CompanyAge
  df.pr$Time2Exit_ZS <- df.pr$Time2Exit + pmax(0, df.pr$d - 10)
  
  # order by exit date
  df.pr <- df.pr[order(df.pr$Date), ]

  return(df.pr)
}

simulate$Timing.Simulator(simulate$create.sim.input(), simulate$create.public.scenario())


## 2) Simulate from Cox Weibull x Two-part Hurdle Model -------
simulate$df.out2cumcf <- function(df.out){
  df.cf <- data.frame(Time = seq(0,15,0.25), CF = 0)
  df.cf2 <- data.frame(Time = round(df.out$Time2Exit*4)/4, CF = df.out$Multiple)
  df.cf2$Time <- ifelse(df.cf2$Time > 15, 15, df.cf2$Time)
  df.cf <- rbind(df.cf, df.cf2)
  df.cf <- aggregate(CF ~ Time, df.cf, sum)
  df.cf$CumCF <- cumsum(df.cf$CF)
  cum.cf <- df.cf$CumCF
  names(cum.cf) <- df.cf$Time
  return(cum.cf)
}

simulate$Final.Simulator <- function(iterations, N, Type, fixed.pub.sce = TRUE){
  all_out <- list()
  # fixed public scenario
  if(fixed.pub.sce){
    pub.sce <- simulate$create.public.scenario()
    all_out$PublicScenario <- pub.sce
    print(tail(pub.sce))
  }

  for(element.N in N){
    print(paste(element.N, "company portfolio"))
    # fixed private portfolio
    df.companies <- simulate$create.sim.input(element.N, Type)
    all_out$PF[[paste("N", element.N, sep ="_")]] <- df.companies
    print(df.companies)
    
    corner.list <- list()
    df.out.list2 <- list()
    for(corner in c("independent", "Joe180")){
      print(corner)
      final.multiple <- list()
      df.out.list1 <- list()
      for(i in seq(iterations)){
        # Dynamic Public Scenario
        if(!fixed.pub.sce){
          pub.sce <- simulate$create.public.scenario()
        }
        
        # Copula
        joe.para <- ifelse(df.companies$Type == "VC", coef2part$VC$joe.para, coef2part$BO$joe.para)
        U.copula <- simulate$TM.copula(nrow(df.companies), corner, joe.para)

        # Timing 
        df.out <- simulate$Timing.Simulator(df.pr = df.companies,
                                   df.pu = pub.sce,
                                   U.timing = U.copula$Timing)
        
        # Multiple
        if(FALSE) {
          df.out$Multiple <- multiple$predict_dgh(sim_out = TRUE,
                                                  U.multiple = U.copula$Multiple,
                                                  sim_in = df.out)
        }
        df.out$Multiple <- copula2part$predict.2part(U.multiple = U.copula, sim_in = df.out)
        
        # Result
        # round_df(df.out, 3)
        
        # Cash Flow on Current NAV
        final.multiple[[i]] <- sum(df.out$RVPI / sum(df.out$RVPI) * df.out$Multiple)
        df.out.list1[[i]] <- simulate$df.out2cumcf(df.out)
      }
      
      final.multiple <- unlist(final.multiple)
      corner.list[[corner]] <- final.multiple
      df.out.list2[[corner]]  <- df.out.list1
      
    }
    all_out$Output[[paste("N", element.N, sep ="_")]] <- corner.list
    all_out$df.out[[paste("N", element.N, sep ="_")]] <- df.out.list2
  }
  
  return(all_out)
}
x <- simulate$Final.Simulator(5, c(3, 15), "VC")


## x) Attach new environment -----
while("simulate" %in% search())
  detach("simulate")
attach(simulate)