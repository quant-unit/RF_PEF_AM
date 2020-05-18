###########################
## Exit Dynamics: AC copula
## inference for margins ##
###########################
## x) Create new environment  ------
copula2part <- new.env()

## Clayton Copula <<<<< ---------
ClayCop <- function(u,v,para) {
    (u^(-para) + v^(-para) - 1)^(-1/para)
}
ClayCop(0.9,0.9,4)

ClayCopDerv <- function(u,v,para) {
  (para+1)*((u*v)^(-(para+1)))*(u^(-para) + v^(-para) - 1)^(-(2*para+1)/para)
}
ClayCopDerv(0.2,0.2,4)

ClayCopFPD <- function(u, v, para) {
  (u^(-para)+v^(-para)-1)^(1/para-1)*u^(-para-1)
}

## Weibull Survival <<< -------

weibull.Surv <- function(start, end, fund.type, fund.age, return.Survival) {
  para_fit <- wb.cox[[fund.type]]$Full$par
  
  beta_MSCI = para_fit["MSCI"]
  beta_HYS = para_fit["HYS"] 
  beta_FundAge = para_fit["FundAge"]
  scale_wb = para_fit["Scale"]
  shape_wb = para_fit["Shape"]

    Surv <- timing$Li_WBcox(from_Date = start, 
                  to_Date = end,
                  beta_MSCI = para_fit["MSCI"], 
                  beta_HYS = para_fit["HYS"] , 
                  beta_FundAge = para_fit["FundAge"],
                  scale_wb = para_fit["Scale"], 
                  shape_wb = para_fit["Shape"],
                  FundAge = fund.age,
                  public.input = public.data,
                  return.Surv = return.Survival)
  return(Surv)
}
weibull.Surv("2002-12-31", "2011-03-31", "BO", 4, TRUE)
weibull.Surv("2002-12-31", "2011-03-31", "BO", 4, FALSE)


## Multiple - Data Prep -----------
if(TRUE) {
  G.p2c <- g.p2c2 # data set used for timing regression
  
  G.p2c <- G.p2c[!is.na(G.p2c$P2C.multi1),] # JUST EXITED INVESTMENTS
  G.p2c  <- G.p2c[G.p2c$Fund_InvestTypes %in% c("BO","VC"),]
  regression_variables <- c("P2C.multi1","RVPI_1","Holding_Period","Time2Exit","Time2Exit_ZS","ZombieStage","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
  G.p2c <- G.p2c[complete.cases(G.p2c[,regression_variables]),]
  
  # delete young entries (bias correction)
  G.p2c <- G.p2c[G.p2c$Investment_Date < as.Date("2010-01-01"), ]
  # filter by RVPI
  G.p2c <- G.p2c[G.p2c$RVPI_1 > -0.9, ]
  G.p2c <- G.p2c[order(G.p2c$Fund_Emi_ID,G.p2c$FundInv_Quarter), ]
  
  summary(G.p2c$P2C.multi1[G.p2c$P2C.multi1 > 0])
  
} else {
  setwd(wd$data)
  G.p2c <- readRDS("df_all1.RDS")
  # G.p2c <- G.p2c[G.p2c$Fund_Region == "US", ]
}


sample.subset <- function(fund.type, df) {
  df <- df[df$Fund_InvestTypes %in% c(fund.type), ]
  dt <- data.table::as.data.table(df)
  dt2 <- dt[,.SD[sample(.N,1)], by=Company_ID]
  df2 <- as.data.frame(dt2)
  df2 <- df2[, !(colnames(df2) %in% c("Ev", "Ev2", "MOIC.dyn", "MOIC.end"))]
  return(df2)
}

get.first.obs <- function(fund.type, df) {
  df2 <- df[df$Fund_InvestTypes == fund.type, ]
  df2 <- as.data.frame(df2 %>% group_by(Fund_Emi_ID) %>% filter(FundInv_Quarter == min(FundInv_Quarter)))
  return(df2)
}

# Summary Statistics (How many observations, i.e. company IDs)
list(First_Entry_Date = min(G.p2c$Investment_Date),
     BO_stochastic = nrow(sample.subset("BO", G.p2c)),
     BO_first = nrow(get.first.obs("BO", G.p2c)),
     VC_stochastic = nrow(sample.subset("VC", G.p2c)),
     VC_first = nrow(get.first.obs("VC", G.p2c)) )
## Fit Joint Model <<<< ------------
library(gamlss)
library(VineCopula)

hurdle0 <- 0 # truncation for zero hurdle unnecessary
gamlss.tr::gen.trun(par = hurdle0, family = GA, name = "tr", type = "left")
gamlss.tr::gen.trun(par = hurdle0, family = GG, name = "tr", type = "left")
gamlss.tr::gen.trun(par = hurdle0, family = BCTo, name = "tr", type = "left")

f <- function(x) x * dGAtr(x, mu = 3, sigma = 2)
integrate(f, 0.000001, 1000)

estimate.bi.model <- function(f.type) {
  no.iterations <- 1000
  final.result <- list()
  
  response1 <- "P2C.multi1"

  predictors_0 <- c("RVPI_1", "Holding_Period", "ML_HYOAS.quarter", "MSCI.Multiple.Exit_1")
  predictors_1mu <- c("Time2Exit", "Holding_Period", "ML_HYOAS.quarter", "MSCI.Multiple.Exit_1")
  predictors_1sigma <- c("Holding_Period", "Time2Exit")

  formula0 <- as.formula(paste("hurdle0", paste(predictors_0, collapse=" + "), sep=" ~ "))
  formula0_IO <- as.formula(paste("hurdle0", 1, sep=" ~ "))
  formula1 <- as.formula(paste(response1, paste(predictors_1mu, collapse=" + "), sep=" ~ "))
  formula1sigma <- as.formula(paste("", paste(predictors_1sigma, collapse=" + "), sep=" ~ "))
  formula1_IO <- as.formula(paste(response1, 1, sep=" ~ "))
  
  G.p2c$hurdle0 <- ifelse(G.p2c[,response1] > hurdle0, 1, 0)
  
  full.entry2exit.list <- list()
  if(TRUE) {
    m0.full <- gamlss(formula0, data = G.p2c[G.p2c$Fund_InvestTypes == f.type, !(colnames(G.p2c) %in% c("Ev", "Ev2", "MOIC.dyn", "MOIC.end"))], family = BI(mu.link = logit))
    m1.full <- gamlss(formula1, sigma.formula = formula1sigma, 
                      data = subset(G.p2c[, !(colnames(G.p2c) %in% c("Ev", "Ev2", "MOIC.dyn", "MOIC.end"))], 
                                    hurdle0 == 1 & Fund_InvestTypes == f.type), family = GAtr)
    
    m0.entryexit <- gamlss(formula0, data = G.p2c[G.p2c$Holding_Period == 0 &  G.p2c$Fund_InvestTypes == f.type, !(colnames(G.p2c) %in% c("Ev", "Ev2", "MOIC.dyn", "MOIC.end"))], family = BI(mu.link = logit))
    m1.entryexit <- gamlss(formula1, sigma.formula = formula1sigma, 
                           data = subset(G.p2c[G.p2c$Holding_Period == 0, !(colnames(G.p2c) %in% c("Ev", "Ev2", "MOIC.dyn", "MOIC.end"))], 
                                         hurdle0 == 1 & Fund_InvestTypes == f.type), family = GAtr)
    
    full.entry2exit.list$m0.full <- m0.full
    full.entry2exit.list$m1.full <- m1.full
    full.entry2exit.list$m0.entryexit <- m0.entryexit
    full.entry2exit.list$m1.entryexit <- m1.entryexit
  }

  set.seed(99)
  iter.list <- list()
  for(i in 1:no.iterations) {
    print(i)
    df <- sample.subset(f.type, G.p2c)
    # df$hurdle0 <- ifelse(df[,response1] > hurdle0, 1, 0)
    
    glms <- list()
    glms.out <- list()
    # Full Covariate Model
    glms$m0 <- gamlss(formula0, data = df, family = BI(mu.link = logit), control = gamlss.control(trace = FALSE))
    try.m1 <- try(gamlss(formula1, sigma.formula = formula1sigma, data = subset(df, hurdle0 == 1), family = GAtr, control = gamlss.control(trace = FALSE)))
    
    if(class(try.m1) == "try-error") {
      next
    }
    glms$m1 <- try.m1
    
    glms.out$m0$mu.coefficients <- glms$m0$mu.coefficients
    glms.out$m0$fitness <- c(AIC = glms$m0$aic, BIC = glms$m0$sbc)
    
    glms.out$m1$mu.coefficients <- glms$m1$mu.coefficients
    glms.out$m1$sigma.coefficients <- glms$m1$sigma.coefficients
    glms.out$m1$fitness <- c(AIC = glms$m1$aic, BIC = glms$m1$sbc)
    glms.out$m1$Rsq <- Rsq(try.m1, type = "both")
    
    # Intercept Only
    glms$m0_IO <- gamlss(formula0_IO, data = df, family = BI(mu.link = logit), control = gamlss.control(trace = FALSE))
    glms$m1_IO <- gamlss(formula1_IO, data = subset(df, hurdle0 == 1), family = GAtr, control = gamlss.control(trace = FALSE))
    
    glms.out$m0_IO$fitness <- c(AIC = glms$m0_IO$aic, BIC = glms$m0_IO$sbc)
    glms.out$m1_IO$fitness <- c(AIC = glms$m1_IO$aic, BIC = glms$m1_IO$sbc)
    
    
    # predict quantiles
    df <- df[df$hurdle0 == 1, ]
    
    Q.multiple <- pGAtr(q = df$P2C.multi1, 
                      mu = exp(predict(glms$m1_IO, what = "mu")), 
                      # nu = glms$m1$nu.coefficients["(Intercept)"]
                      sigma = exp(predict(glms$m1_IO, what = "sigma")))
    d.multiple <- dGAtr(x = df$P2C.multi1, 
                      mu = exp(predict(glms$m1_IO, what = "mu")), 
                      # nu = glms$m1$nu.coefficients["(Intercept)"]
                      sigma = exp(predict(glms$m1_IO, what = "sigma")))
    
    glms.out$m1$Q.multiple <- Q.multiple
    
    Q.timing <- apply(df, 1, function(x) {
      S1 <- weibull.Surv(x[["Investment_Date"]], x[["Exit_Date"]], 
                         f.type, as.numeric(x[["FundAgeAI"]]), TRUE)
      S2 <- weibull.Surv(x[["Investment_Date"]], x[["FundInv_Quarter"]], 
                         f.type, as.numeric(x[["FundAgeAI"]]), TRUE)
      return(S1/S2)
    })
    glms.out$t$Q.timing <- Q.timing
    
    
    Q.timing_previousQuarter <- apply(df, 1, function(x) {
      beforeExitDate <- base::as.Date(cut(as.Date(x[["Exit_Date"]]), "quarter")) - 1
      S1 <- weibull.Surv(x[["Investment_Date"]], beforeExitDate, 
                         f.type, as.numeric(x[["FundAgeAI"]]), TRUE)
      S2 <- weibull.Surv(x[["Investment_Date"]], x[["FundInv_Quarter"]], 
                         f.type, as.numeric(x[["FundAgeAI"]]), TRUE)
      return(S1/S2)
    })

    d.survival <- apply(df, 1, function(x) {
      weibull.Surv(x[["Investment_Date"]], x[["Exit_Date"]], 
                   f.type, as.numeric(x[["FundAgeAI"]]), FALSE)  })
    
    psObs <- function(x) {
      x <- ifelse(x == 0, 0.0001, x)
      x <- ifelse(x == 1, 0.9999, x)
      return(x)
    }
    
    plot(Q.timing, Q.multiple)
    
    # continuous Copula
    loliCON <- function(para) {
      # -sum(log(d.survival * d.multiple * ClayCopDerv(pobs(Q.multiple), pobs(Q.timing),para)))
      -sum(log(d.survival * d.multiple * BiCopPDF(pobs(Q.multiple), pobs(Q.timing), par = para, family = 16)))
      }
    # plot(BiCopSim(2500, family = 16, par = 4))

    glms.out$m1t$claytonCON <- optimize(loliCON, interval = c(1.001, 10))
    glms.out$m1t$claytonCON
    if(is.na(glms.out$m1t$claytonCON$objective)) {
      print("next CON")
      next
    }

    # discrete Copula (does not work)
    if(FALSE) {
      loliDIS <- function(para) {
        -sum(log(d.multiple * 
                   (ClayCopFPD(pobs(Q.multiple), pobs(Q.timing),para) - ClayCopFPD(pobs(Q.multiple), pobs(Q.timing_previousQuarter),para))
        ))
      }
      summary(Q.timing)
      summary(Q.timing_previousQuarter)
      summary(ClayCopFPD(pobs(Q.multiple), pobs(1-Q.timing),2) - ClayCopFPD(pobs(Q.multiple), pobs(1-Q.timing_previousQuarter),2))
      loliDIS(0.1)
      glms.out$m1t$claytonDIS <- optimize(loliDIS, interval = c(0,10))
      glms.out$m1t$claytonDIS
      if(is.na(glms.out$m1t$claytonDIS$objective)) {
        print("next DIS")
        next
      }
    }

    
    iter.list[[paste("i",i,sep="")]] <- glms.out
  }
  
  # Output Preparation
  if(TRUE) {
    
    coefs <- list()
    
    # Logit
    x0 <- sapply(iter.list, function(x){
      x$m0$mu.coefficients
    })
    coefs$m0$mu$obs <- x0
    coefs$m0$mu$Mean <- colMeans(data.frame(t(x0)))
    coefs$m0$mu$SD <- apply(data.frame(t(x0)), 2, sd)
    
    # Gamma
    for(coef in c("mu", "sigma")) {
      obs <- sapply(iter.list, function(x){
        x$m1[[paste(coef, "coefficients", sep = ".")]]
      })
      coefs$m1[[coef]][["obs"]] <- obs
      coefs$m1[[coef]][["Mean"]] <- colMeans(data.frame(t(obs))) # mu
      coefs$m1[[coef]][["SD"]] <- apply(data.frame(t(obs)), 2, sd)
    }
    
    Rsquared <- sapply(iter.list, function(x){
      unlist(x$m1$Rsq)
    })
    
    coefs$m1$Rsq <- list(Mean = colMeans(data.frame(t(Rsquared))),
                         SD = apply(data.frame(t(Rsquared)), 2, sd))
    
    # Bivariate
    U.Timing <- obs <- sapply(iter.list, function(x){
      x$t$Q.timing
    })
    U.Multiple <- obs <- sapply(iter.list, function(x){
      x$m1$Q.multiple
    })
    CopDF <- data.frame(Timing = unlist(U.Timing), Multiple = unlist(U.Multiple))
    
    cor(Q.timing, Q.multiple, method = "kendall")
    plot(Q.timing, Q.multiple)
    
    ccp <- sapply(iter.list, function(x){
      x$m1t$claytonCON$minimum
    })
    coefs$copCON$Mean <- mean(ccp)
    coefs$copCON$SD <- sd(ccp)
    
    # AIC and BIC
    get.fitness <- function(model) {
      y <- sapply(iter.list, function(x){
        x[[model]]$fitness
      })
      
      Mean <- colMeans(data.frame(t(y)))
      SD <- apply(data.frame(t(y)), 2, sd)
      return(list(c(Model = model, Var = "Mean", Mean), 
                  c(Model = model, Var = "SD", SD)))
    }
    
    baic <- list()
    for(model in c("m0", "m0_IO", "m1", "m1_IO")) {
      baic[[model]] <- get.fitness(model)
    }
    
    df.baic <- data.frame(rbind_all(combine(baic)))
    
    final.result <- list(coefs = coefs, BAIC = df.baic, CopDF = CopDF,
                         Iter = list(no.iterations = no.iterations, no.converged = length(iter.list)),
                         full.entry2exit = full.entry2exit.list)
    return(final.result)
  }
}
if(FALSE) {
  system.time(bo.GA <- estimate.bi.model("BO"))
  system.time(vc.GA <- estimate.bi.model("VC"))
  setwd(wd$data.out)
  saveRDS(bo.GA, "multiple_bo.GA20200508.RDS")
  saveRDS(vc.GA, "multiple_vc.GA20200508.RDS")
} else {
  setwd(wd$data.out)
  bo.GA <- readRDS("multiple_bo.GA20200508.RDS")
  vc.GA <- readRDS("multiple_vc.GA20200508.RDS")
}

summary(bo.GA$full.entry2exit$m1.full)
summary(bo.GA$full.entry2exit$m1.entryexit)
summary(vc.GA$full.entry2exit$m1.full)
summary(vc.GA$full.entry2exit$m1.entryexit)

bo.GA$coefs$copCON
vc.GA$coefs$copCON

plot(BiCopSim(2500, family = 16, par = 1.095231))


plot.pairs <- function(df, do_eps = FALSE, name = "") {
  old_par <- par()
  if(do_eps) {
    old.wd <- getwd()
    setwd(wd$eps)
    setEPS()
    postscript(paste("dependency", name, ".eps", sep = ""), width = 3, height = 3,
               family = "Helvetica", pointsize = 3)
  }

  par(mfrow = c(2, 2), mar = c(4.1, 4.1, 1, 1), oma = c(0,0,1,0), cex = 1.3)
  
  # df 1
  hist(df$Multiple, main = "", xlab = "Multiple", probability = TRUE)
  abline(h = 1, col = "darkgrey", lwd = 2)
  plot(df[1:1500, ], cex = 0.5)
  BiCopMetaContour(df$Timing, df$Multiple, 
                   # levels = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5),
                   levels = c(0.2,0.9, 1.3),
                   margins = "unif",
                   xlab = "Timing", ylab = "Multiple")
  hist(df$Timing, main = "", xlab = "Timing", probability = TRUE)
  abline(h = 1, col = "darkgrey", lwd = 2)
  
  mtext(name, line = -0.2, cex = 1.5, font = 2, outer = TRUE)
  
  if(do_eps) {
    dev.off()
    setwd(old.wd)
  }
  par(old_par)
}
# plot.pairs(bo.GA$CopDF, name = "BO") 
# plot.pairs(vc.GA$CopDF, name = "VC")

# View Results <<<<<<<< ---------
coefs <- bo.GA$coefs
coefs <- vc.GA$coefs
round(coefs$m0$mu$Mean ,3)
round(coefs$m0$mu$SD ,3)
round(coefs$m1$mu$Mean ,3)
round(coefs$m1$mu$SD ,3)
round(coefs$m1$sigma$Mean, 3)
round(coefs$m1$sigma$SD, 3)

vc.GA$BAIC
bo.GA$BAIC

make.coef2part <- function() {
  coef2part <- list()
  coef2part$BO$m0mu <- bo.GA$coefs$m0$mu$Mean
  coef2part[["BO"]][["m1mu"]] <- bo.GA$coefs$m1$mu$Mean
  coef2part[["BO"]][["m1sigma"]] <- bo.GA$coefs$m1$sigma$Mean
  coef2part$BO$joe.para <- bo.GA$coefs$copCON$Mean
  
  coef2part[["VC"]][["m0mu"]] <- vc.GA$coefs$m0$mu$Mean
  coef2part[["VC"]][["m1mu"]] <- vc.GA$coefs$m1$mu$Mean
  coef2part[["VC"]][["m1sigma"]] <- vc.GA$coefs$m1$sigma$Mean
  coef2part$VC$joe.para <- vc.GA$coefs$copCON$Mean
  
  return(coef2part)
}

copula2part$coef2part <- make.coef2part()
