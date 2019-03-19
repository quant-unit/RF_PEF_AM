########################
## Exit Dynamics: Timing
########################
## x) Create new environment  ------
timing <- new.env()

## 1) non-parametric (partial likelihood) -------
timing$NonParaCoxRegression <- function(G.p2c2){

  G.p2c2 <- G.p2c2[order(G.p2c2$Fund_Emi_ID,G.p2c2$FundInv_Quarter), ]
  
  # equidistant time-grid
  G.p2c2$T1_r <- round(G.p2c2$T1 * 4, 0) / 4
  G.p2c2$T2_r <- round(G.p2c2$T2 * 4, 0) / 4
  
  outlist <- list()
  
  for(Fund_Types in c("BO", "VC")){
    
    SurvData <- G.p2c2[G.p2c2$Fund_InvestTypes %in% Fund_Types,]
    attach(SurvData)
    intcen.Surv <- survival::Surv(time= T1_r, time2= T2_r, event= Ev)
    detach(SurvData)
    rm(SurvData)
    
    kaplan.meier <- survival::survfit(intcen.Surv ~ 1)

    # Cox Model
    re.Cox <-  survival::coxph(intcen.Surv ~ 
                                 # GLOQuarter +
                                 ML_HYOAS.quarter +
                                 # MSCI.multi.hist +
                                 MSCI.qrtly.return +
                               FundAgeAI,
                               data = G.p2c2[G.p2c2$Fund_InvestTypes %in% Fund_Types,], ties="efron")
    
    # Cumulative basehazard
    bh.Cox <- survival::basehaz(re.Cox,centered=FALSE) # Covariates == 0
    bh.Cox$S.base <- exp(-1*bh.Cox$hazard)
    bh.Cox$hazard_rate <- c(0,diff(bh.Cox$hazard))
    
    bh.Cox.Mean <- survival::basehaz(re.Cox,centered=TRUE) # Covariates == 0
    bh.Cox.Mean$S.base <- exp(-1*bh.Cox.Mean$hazard)
    bh.Cox.Mean$hazard_rate <- c(0,diff(bh.Cox.Mean$hazard))
    
    outlist[[Fund_Types]] <- list(Regression= re.Cox, BaseHaze_Zero= bh.Cox, 
                                  BaseHaze_Mean= bh.Cox.Mean, KaMe = kaplan.meier)
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



## 2) parametric Weibull (full likelihood) -------
# Likelihhod function for one observation
timing$Li_WBcox <- function(from_Date, to_Date,
                            beta_MSCI, beta_HYS, beta_FundAge,
                            scale_wb, shape_wb,
                            FundAge,
                            public.input,
                            return.Surv=FALSE){
  
  from_Date <- as.Date(from_Date)
  
  if(is.na(to_Date)){
    to_Date <- as.Date("2016-12-31")
    Censoring <- TRUE
  } else {
    to_Date <- as.Date(to_Date)
    Censoring <- FALSE
  }
  
  puda <- public.input[public.input$Date >= from_Date & public.input$Date <= to_Date,]
  puda$t <- as.numeric(puda$Date - puda$Date[1])/365.25
  puda$FundAge <- FundAge
  
  # Weibull: Cumulative Hazard Function
  puda$base_haze_rate_wb_exact <- c(0,diff( (puda$t/scale_wb)^shape_wb ) )
  
  # Gompertz: Cumulative Hazard Function
  # puda$base_haze_rate_wb_exact <- c(0,diff(scale_wb/shape_wb * (exp(puda$t * shape_wb) - 1)  ) )
  
  # Loglogistic: Cumulative Hazard Function
  # puda$base_haze_rate_wb_exact <- c(0,diff( log(1 + exp(scale_wb) * puda$t^shape_wb ) ) )
  
  
  puda$expBX <- exp( as.matrix(puda[,c("MSCI_monthly_return","ML_HYOAS","FundAge")]) %*% 
                       c(beta_MSCI,beta_HYS,beta_FundAge) )
  
  puda$Cum_Haze_exact <- cumsum(as.numeric(puda$base_haze_rate_wb_exact) * puda$expBX)
  
  puda$Surv_WB <- exp( - puda$Cum_Haze_exact)
  
  
  if(Censoring | return.Surv){
    out <- puda$Surv_WB[nrow(puda)]
  } else {
    out <- - diff(puda$Surv_WB)
    out <- sum(tail(out, 3)) # integrate over last three months
  }
  return(out)
}

timing$ParaCoxRegression <- function(public.data, data_LoLi){
  
  # Full Parameter Log Likelihood
  LoLi_WBcox <- function(par, data, min.it = TRUE){

    full_para <- c("MSCI", "HYS", "FundAge", "Scale", "Shape")
    paras <- full_para %in% names(par)
    x <- rep(0, length(full_para))
    names(x) <- full_para
    x[paras] <- par
    
    
    fund_ids <- base::unique(data$Fund_ID) 
    loglikelihood_outer <- list()

    for(j in fund_ids){
      sub.data <- data[data$Fund_ID == j, ]
      loglikelihood_inner <- list()
      
      
      for(i in 1:nrow(sub.data)){
        L <- timing$Li_WBcox(from_Date = sub.data$InvDate[i],
                      to_Date = sub.data$ExitDate[i],
                      x[1], x[2], x[3], # betas
                      x[4], x[5],       # shape & scale
                      FundAge = max(0, sub.data$YearInvest[i] - sub.data$Vintage[i]),
                      public.input = public.data)
        
        
        if(is.na(L)){
          print(paste(">> Fund",j, "Company", i))
          print(x)
        }
        
        if(L == 0){
          print(paste("Fund",j, "Company", i))
          print(round(x,3))
        }
        
        loglikelihood_inner[i] <- log(L)
      }
      
      loglikelihood_outer[[paste("r",j,sep = "")]] <- sum(as.numeric(loglikelihood_inner))
    }
    
    out <- sum(as.numeric(loglikelihood_outer))
    
    if(min.it){ # converts it to a minimization problem
      out <- - out
    }

    return(out)
  }
  
  
  # Create Start Values and Bounds for all fits
  fit.list <- list()
  fit.list[["WB_only"]] <- list(Start.Para = c(Scale = 2, Shape = 0.5))
  fit.list[["Full"]] <- list(Start.Para = c(MSCI = 1, HYS = -2, FundAge = 0.05, Scale = 2, Shape = 0.5))
  fit.list[["HYS_Age"]] <- list(Start.Para = c(HYS = -2, FundAge = 0.05, Scale = 2, Shape = 0.5))
  fit.list[["MSCI_Age"]] <- list(Start.Para = c(MSCI = 1, FundAge = 0.05, Scale = 2, Shape = 0.5))
  fit.list[["Age"]] <- list(Start.Para = c(FundAge = 0.05, Scale = 2, Shape = 0.5))
  
  for(fits in names(fit.list)){
    len <- length(fit.list[[fits]]$Start.Para)
    fit.list[[fits]]$Bounds.Low <- c(rep(-Inf, len - 2), 0, 0)
    fit.list[[fits]]$Bounds.Upp <- c(rep(Inf, len))
  }

  
  # Optimize Maximum Likelihood
  fund.types <- base::unique(data_LoLi$Type)
  print(table(data_LoLi$Type))
  optimX_result <- list()
  
  for(fund_type in fund.types){
    
      data_sample <- data_LoLi[data_LoLi$Type == fund_type, ]
      
      for(fits in names(fit.list)){
        optimX_fit <- optimx::optimx(par = fit.list[[fits]][["Start.Para"]], 
                                      fn = LoLi_WBcox,
                                      data = data_sample,
                                      min.it = TRUE, # makes it a min problem
                                      method = "nlminb",
                                      lower = fit.list[[fits]][["Bounds.Low"]],
                                      upper = fit.list[[fits]][["Bounds.Upp"]],
                                      control = list(all.methods = FALSE, 
                                                     save.failures = TRUE),
                                      hessian = TRUE)
        
        optimX_prep <- list()
        optimX_prep[["summary"]] <- optimX_fit
        
        optimX_prep[["value"]] <- - optimX_fit$value # maximum log likelihood estimate
        
        par <- as.numeric(coef(optimX_fit))
        names(par) <- colnames(coef(optimX_fit))
        optimX_prep[["par"]] <- par
        
        
        optimX_prep["hessian"] <- attr(optimX_fit, "details")[, "nhatend"]

        # optimX_prep["gradient"] <- attr(optimX_fit, "details")[, "ngatend"]
        # optimX_prep["hessian_eigenvalues"] <- attr(optimX_fit, "details")[, "hev"]
        
        optimX_result[[fund_type]][[fits]] <- optimX_prep
        
        
        
        
        # old method optim
        
        '
        optim_result[[fits]][[fund_type]] <- optim(par = fit.list[[fits]][["Start.Para"]],
        fn = LoLi_WBcox, 
        data = data_sample,
        mode = fits,
        method = "L-BFGS-B", 
        lower = fit.list[[fits]][["Bounds.Low"]],
        upper = fit.list[[fits]][["Bounds.Upp"]],
        hessian = TRUE)
        
        '
      }
      
    }

  
  # Calculate SEs, p-values and AICs
  SEandAIC <- function(fit){
    # 95% Conf_Intervals
    # https://stats.stackexchange.com/questions/27033/in-r-given-an-output-from-optim-with-a-hessian-matrix-how-to-calculate-paramet
    
    fisher_info <- solve(fit$hessian) # sign issue (maximization vs. minimization)

    prop_sigma <- sqrt(diag(fisher_info))
    
    result <- data.frame(MLE = fit$par, SE = prop_sigma)
    result$t_value <- result$MLE / result$SE
    result$p_value <- pnorm(-abs(result$t_value))

    rownames(result) <- names(fit$par)
    
    print(round(result,4))
    #prop_sigma<-diag(prop_sigma)
    #upper<-fit$par+1.96*prop_sigma
    #lower<-fit$par-1.96*prop_sigma
    #interval<-data.frame(value=fit$par, upper=upper, lower=lower)
    
    fit$Coefs <- result
    fit$AIC <-  2 * length(fit$par) -  2 * fit$value
    # fit$BIC <- log(n) * length(fit$par) -  2 * fit$value
    
    return(fit)
  }
  
  for(fit.names in names(fit.list)){
    for(fund_type in fund.types){
      print(paste(">>>", fund_type, fit.names))
      
      optimX_result[[fund_type]][[fit.names]] <- SEandAIC(optimX_result[[fund_type]][[fit.names]])
    }
  }

  
  return(optimX_result)
}



# fitWeib <- survreg(Surv(Timing, RightCensored) ~ 1, dist="weibull", data=data_LoLi[data_LoLi$Timing > 0,])
# summary(fitWeib)

## 3) Plot: Survival Function & Rosenblatt Transform --------
timing$PlotCoxWeibullSurvival <- function(eps = FALSE,
                                   para_in,
                                   non_para_in){
  if(eps){
    setEPS() ; postscript("Timing Cox Weibull Base Survival Function 2.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=5)
  }
  bo_col <- "royalblue1" ; vc_col <- "maroon1"
  
  sc_bo <- para_in$BO$Full$par[4]
  sh_bo <- para_in$BO$Full$par[5]
  sc_vc <- para_in$VC$Full$par[4]
  sh_vc <- para_in$VC$Full$par[5]
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


timing$Weibull_Roseblatt <- function(eps = FALSE,
                              data_LoLi,
                              public.data,
                              Para_Cox_Weibull){
# Define Functions  
  S_WBcox <- function(from_Date, to_Date,
                      beta_MSCI, beta_HYS, beta_FundAge,
                      scale_wb, shape_wb,
                      FundAge, # at Entry
                      public.data){
    from_Date <- as.Date(from_Date)
    if(is.na(to_Date)){
      to_Date <- as.Date("2016-12-31")
      Censoring <- TRUE
    }else{
      to_Date <- as.Date(to_Date)
      Censoring <- FALSE
    }
    public.data <- public.data[public.data$Date >= from_Date & public.data$Date <= to_Date,]
    public.data$t <- as.numeric(public.data$Date - public.data$Date[1])/365.25
    public.data$FundAge <- FundAge
    
    public.data$base_haze_rate_wb_exact <- c(0,diff( (public.data$t/scale_wb)^shape_wb ))
    
    public.data$expBX <- exp( as.matrix(public.data[,c("MSCI_monthly_return","ML_HYOAS","FundAge")]) %*% 
                                c(beta_MSCI,beta_HYS,beta_FundAge) )
    
    public.data$Cum_Haze_exact <- cumsum(as.numeric(public.data$base_haze_rate_wb_exact) * public.data$expBX)
    public.data$Surv_WB <- exp(-public.data$Cum_Haze_exact)
    
    if(Censoring){
      out <- public.data$Surv_WB[nrow(public.data)]  * runif(1) # basically NA
    }else{
      out <- public.data$Surv_WB[nrow(public.data)]
    }
    return(out)
  }
  
  Rosenblatt_WBcox <- function(FundType){
    
    data <- data_LoLi[data_LoLi$Type == FundType,]
    x <- Para_Cox_Weibull[[FundType]]$Full$par
    
    outlist <- list()
    for(i in 1:nrow(data)){
      S <- S_WBcox(from_Date = data$InvDate[i],
                   to_Date = data$ExitDate[i],
                   x[1], x[2], x[3], x[4], x[5],
                   FundAge = max(0, data$YearInvest[i] - data$Vintage[i]),
                   public.data = public.data)
      
      outlist[i] <- S
    }
    return(as.numeric(outlist))
  }
  
# Create Plot  
  if(eps){
    setEPS() ; postscript("Timing_Rosenblatt.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=5)
  }
  bo_col <- "royalblue1" ; vc_col <- "maroon1"
  
  out <- list()
  par(mar=c(3,3,2,1),mfrow=c(1,2),oma=c(2.5,2,1,1))
  for(type in c("BO","VC")){
    robla <- Rosenblatt_WBcox(type)
    out[[type]] <- robla
    
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
  
  return(out)
}


## 4) CSV: Regression Coefs & Parameters -------
timing$CreateOutput <- function(para.estim, do.csv){
  iterlist <- list()
  for(types in names(para.estim)){
    print(types)
    for(modelz in names(para.estim[[1]])){
      
      df1 <- data.frame(do.call(rbind, list(
        data.frame(Variable = rownames(para.estim[[types]][[modelz]]$Coefs),
                   Stats = "Mean",
                   Value = para.estim[[types]][[modelz]]$Coefs$MLE),
        data.frame(Variable = rownames(para.estim[[types]][[modelz]]$Coefs),
                   Stats = "SD",
                   Value = para.estim[[types]][[modelz]]$Coefs$SE),
        data.frame(Variable = "AIC", Stats = "Mean", Value = para.estim[[types]][[modelz]]$AIC) ) ) )
      df1$Type <- types
      df1$Model <- modelz
      
      iterlist[[paste(types, modelz, sep="!")]] <- df1
    }
  }
  
  df <- data.frame(do.call(rbind, iterlist))
  
  # Reshape to wide format
  df2 <- reshape(data = df, direction = "wide",
                 v.names = "Value", idvar = c("Type", "Variable", "Stats"), timevar = "Model")
  df2 <- df2[order(df2$Type, factor(df2$Variable, levels = c("MSCI", "HYS", "FundAge", "Scale", "Shape", "AIC"))), ]
  rownames(df2) <- NULL
  df2 <- round_df(df2)
  
  
  # Formatting SD (in parentheses)
  for(colz in names(para.estim[[1]])){
    colz <- paste("Value", colz, sep = ".")
    df2[, colz] <- ifelse(df2$Stats == "Mean", df2[, colz],
                          ifelse(is.na(df2[, colz]), NA,
                                 paste("(", df2[, colz], ")", sep = "")))
  }
  
  
  # create csv
  if (do.csv) {
    csv.dir <- paste(Sys.Date(), "Timing", "Parameters.csv", sep = "_")
    write.csv2(df2, file = csv.dir, na = "", row.names = FALSE)
  } else {
    print(df2)
  }
  
  invisible(df2)
  
  
}


timing$para_coefs_csv <- function(in_list, to.csv = FALSE){
  ff = function(x){ 
    if (class(x) == "list") 
      lapply(x, ff) 
    else if(class(x) == "data.frame") 
      TRUE
    else
      NULL
  }
  lnames = names(unlist(lapply(in_list, ff)))
  varnames = strsplit(lnames, split = ".", fixed = TRUE)
  
  list_df <- lapply(seq_along(varnames), function(i){
    names <- varnames[[i]]
    df <- in_list[[names]]
    df <- round_df(df)
    
    df$Model <- names[2]
    df$Type <- names[1]
    df$Variable <- rownames(df)
    rownames(df) <- NULL
    df <- df[c(5,6,1,2,3,4)]
    return(df)
  })
  
  out_df <- as.data.frame(do.call(rbind, list_df))

  return(out_df)
}
timing$para_AIC_csv <- function(in_list, to.csv = FALSE){
  rlist_list.flatten <- function (x, use.names = TRUE, classes = "ANY") 
  {
    len <- sum(rapply(x, function(x) 1L, classes = classes))
    y <- vector("list", len)
    i <- 0L
    items <- rapply(x, function(x) {
      i <<- i + 1L
      y[[i]] <<- x
      TRUE
    }, classes = classes)
    if (use.names && !is.null(nm <- names(items))) 
      names(y) <- nm
    y
  }
  
  x <- rlist_list.flatten(in_list)
  AIC_BO <- as.data.frame(do.call(rbind, x[grepl("BO", names(x)) & grepl("AIC", names(x))]))
  AIC_VC <- as.data.frame(do.call(rbind, x[grepl("VC", names(x)) & grepl("AIC", names(x))]))
  
  df_AIC <- data.frame(Model = names(wb.cox[[1]]),
                        BO = AIC_BO[,1], 
                        VC = AIC_VC[,1])
  
  df_AIC <- df_AIC[order(df_AIC$BO),]
  df_AIC <- round_df(df_AIC, 2)

  return(df_AIC)
}
timing$non_para_coefs_csv <- function(df.in, create.csv = FALSE){
  sum.bo <- summary(np.cox$BO$Regression)
  sum.bo <- data.frame(sum.bo$coefficients)
  sum.bo$Type <- "BO"
  sum.bo$Variable <- row.names(sum.bo)
  
  sum.vc <- summary(np.cox$VC$Regression)
  sum.vc <- data.frame(sum.vc$coefficients)
  sum.vc$Type <- "VC"
  sum.vc$Variable <- row.names(sum.vc)
  
  df_out <- rbind(sum.bo, sum.vc)
  row.names(df_out) <- NULL
  df_out <- df_out[,c(7,6,1:5)]
  df_out <- round_df(df_out)
  
  if(create.csv){
    write.csv(df_out, "NonParaCoxSummary.csv")
  }
  
  return(df_out)
}
timing$TimingFormatedCSV <- function(para.in, non.para.in, make.csv = FALSE){
  para.coef <- timing$para_coefs_csv(para.in)
  para.aic <- timing$para_AIC_csv(para.in)
  non.para.coef <- timing$non_para_coefs_csv(non.para.in)
  
  re.shaped <- list()
  for(type in c("BO", "VC")){
    v.nameZ <- c("MLE", "SE")
    for(v.name in v.nameZ){
      drop.name <- v.nameZ[v.nameZ != v.name]
      df2 <- reshape(data = para.coef[grep(type, para.coef$Model), ], direction = "wide", 
                     v.names = v.name, idvar = "Variable", timevar = "Model", 
                     drop = c(drop.name, "t_value", "p_value"))
      
      # add (parametric) AICs
      if(v.name == "MLE"){
        AICs <- para.aic[c(1,3,2), type]
        # df2[c(1 + nrow(df2)), ] <- c("AIC", AICs)
      }
      
      colnames(df2) <- gsub(paste(v.name, type, sep = "."), "Coef", colnames(df2))
      df2$M_S <- v.name
      df2$Type <- type
      
      # merge non-parametric frame
      np.coef <- ifelse(v.name == "MLE", "coef", "se.coef.")
      df.np <- non.para.coef[non.para.coef$Type == type, c("Variable", np.coef)]
      colnames(df.np) <- c("Variable", "Coef_np")
      df.np$Variable <- c("HYS", "MSCI", "FundAge")
      df2 <- merge(df2, df.np, by = "Variable", all.x = TRUE)
      
      # brackets to SEs
      if(v.name == "SE"){
        to.sd <- function(x){
          if(is.numeric(x)){
            y <- paste("(",x,")",sep = "")
          } else {
            y <- x
          }
          return(y)
        }
        df2[] <- lapply(df2, to.sd)
        
        to.na <- function(x){
          if(x == "(NA)"){
            y <- NA
          } else {
            y <- x
          }
          return(y)
        }
        df2 <- as.data.frame(apply(df2, c(1,2),  to.na))
      }
      
      
      
      re.shaped[[paste(type, v.name, sep="_")]] <- df2
    }
  }
  
  df3 <- data.frame(do.call(rbind, re.shaped))
  df3 <- df3[order(df3$Type, 
                   factor(df3$Variable, 
                          levels = c("FundAge","MSCI","HYS","Scale","Shape","AIC"))
  ), ]
  
  df3 <- df3[, c("Type", "Variable", "M_S", "Coef_fit", "Coef_fit3", "Coef_fit2", "Coef_np")]
  rownames(df3) <- NULL
  
  if (make.csv) {
    csv.dir <- paste(Sys.Date(), "Timing", "non_and_para.csv", sep = "_")
    write.csv2(df3, file = csv.dir, na = "", row.names = FALSE)
    
    # txt.dir <- paste(Sys.Date(), "Timing", reg_names, "iter.txt", sep = "_")
    # print(xtable::xtable(df3), file = txt.dir)
  } else {
    print(df3)
  }
  
  invisible(df3)
}


## x) Attach new environment -----
#while("timing" %in% search())
#  detach("timing")
#attach(timing)