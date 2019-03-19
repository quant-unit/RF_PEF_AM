##########################
## Exit Dynamics: Multiple
##########################
## x) Create new environment  ------
multiple <- new.env()

## 1) Prep Data (G.p2c, regression_variables)  -----

G.p2c <- g.p2c2 # data set used for timing regression

G.p2c <- G.p2c[!is.na(G.p2c$P2C.multi1),] # JUST EXITED INVESTMENTS
G.p2c  <- G.p2c[G.p2c$Fund_InvestTypes %in% c("BO","VC"),]
regression_variables <- c("P2C.multi1","RVPI_1","Holding_Period","Time2Exit","Time2Exit_ZS","ZombieStage","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
G.p2c <- G.p2c[complete.cases(G.p2c[,regression_variables]),]

# delete young entries (bias correction?)
G.p2c <- G.p2c[G.p2c$Investment_Date < as.Date("2010-01-01"), ]
# filter by RVPI
G.p2c <- G.p2c[G.p2c$RVPI_1 > -0.9, ]
G.p2c <- G.p2c[order(G.p2c$Fund_Emi_ID,G.p2c$FundInv_Quarter), ]


multiple$data <- G.p2c

# Plot Dependent Variables (Entry-to-Exit vs longitudinal, BO vs VC)
multiple$plot.multi.dist <- function(eps = FALSE, static_df = g.sum, dynamic_df = G.p2c){
  if(eps){
    old.wd <- getwd()
    setwd(wd$eps)
    setEPS() ; postscript("Multiple_Dependent_Variable.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=5)
  }
  
  multi.bo.dyn <- dynamic_df$P2C.multi1[dynamic_df$Fund_InvestTypes == "BO"]
  multi.vc.dyn <- dynamic_df$P2C.multi1[dynamic_df$Fund_InvestTypes == "VC"]
  multi.bo.sta <- static_df$MOIC[static_df$Type == "BO"]
  multi.vc.sta <- static_df$MOIC[static_df$Type == "VC"]
  
  bins_by <- 0.1
  bins <- seq(0,100,bins_by)
  hist.bins <- hist(multi.bo.dyn, bins, plot=FALSE)[[1]]
  hist.bins <- hist.bins[-length(hist.bins)]
  hist.bo.dyn <- hist(multi.bo.dyn, bins, plot=FALSE)[["counts"]]
  hist.vc.dyn <- hist(multi.vc.dyn, bins, plot=FALSE)[["counts"]]
  hist.bo.sta <- hist(multi.bo.sta, bins, plot=FALSE)[["counts"]]
  hist.vc.sta <- hist(multi.vc.sta, bins, plot=FALSE)[["counts"]]
  
  ecdf.multi.bo.dyn <- ecdf(multi.bo.dyn)
  ecdf.multi.vc.dyn <- ecdf(multi.vc.dyn)
  ecdf.multi.bo.sta <- ecdf(multi.bo.sta)
  ecdf.multi.vc.sta <- ecdf(multi.vc.sta)
  
  if(FALSE){
    par(mfrow = c(2, 2))
    sun(multi.bo.dyn, main="BO - dynamic",Xlim=c(0,10))
    sun(multi.vc.dyn, main="VC - dynamic",Xlim=c(0,10))
    sun(multi.bo.sta, main="BO - static",Xlim=c(0,10))
    sun(multi.vc.sta, main="VC - static",Xlim=c(0,10))
  }
  
  bo_col <- "royalblue1"
  vc_col <- "maroon1"
  
  curve(ecdf.multi.bo.dyn, 0, 5, col = bo_col, lty = 3, ylim = c(0, 1),
        ylab = "CDF", xlab = "Multiple")
  curve(ecdf.multi.vc.dyn, 0, 10, col = vc_col, add = TRUE, lty = 3)
  curve(ecdf.multi.bo.sta, 0, 10, col = bo_col, add = TRUE)
  curve(ecdf.multi.vc.sta, 0, 10, col = vc_col, add = TRUE)
  
  abline(h = c(0,1), v = 0, col = "grey")
  abline(v = 1, col = "darkgrey", lty = 2)
  
  legend("right",bty = "n", legend = c("BO - static", "VC - static", "BO - dynamic", "VC - dynamic"), lty = c(1,1,3,3), col = c(bo_col, vc_col))
  
  lines(hist.bins, 1 - hist.vc.sta / sum(hist.vc.sta), col = vc_col, type = "s", lty = 2)
  lines(hist.bins, 1 - hist.bo.sta / sum(hist.bo.sta), col = bo_col, type = "s", lty = 2)
  
  
  '
  plot(hist.bins, hist.bo.sta / sum(hist.bo.sta), col = bo_col, type = "h", xlim = c(0, 5),
       ylim = c(-0.06, 0.06),
       # ylim = c(-max(hist.vc.sta / sum(hist.vc.sta)), max(hist.bo.sta / sum(hist.bo.sta))),
       xlab = "Multiple Observations", ylab = "Relative Frequency")
  abline(h = 0, col = "grey", lty=3)
  lines(hist.bins, - hist.vc.sta / sum(hist.vc.sta), col = vc_col, type = "h")
  '
  
  
  
  
  if(eps){ 
    dev.off() 
    setwd(old.wd)
  }
  
  par(par_default)
}



## 2) Define Functions -----

# Create data.frame for test procedure
multiple$creat_reg_df <- function(f.type, df = G.p2c, reg_vars = regression_variables, 
                                  stochastic = TRUE, private.source = "AssetMetrix"){
  if(private.source == "AssetMetrix"){
    reg_vars2 <- reg_vars[!(reg_vars %in% c("Time2Exit","ZombieStage"))]
    df2 <- df[df$Investment_Date < as.Date("2010-01-01"),]
    df2 <- df2[df2$Fund_InvestTypes %in% c(f.type),c(reg_vars2,
                                                     "Fund_Emi_ID", "Company_ID",
                                                     "Fund_InvestTypes",
                                                     "FundInv_Quarter","Investment_Date")]
    if(stochastic){
      dt <- data.table::as.data.table(df2)
      dt2 <- dt[,.SD[sample(.N,1)], by=Company_ID]
      df2 <- as.data.frame(dt2)
      # df2 <- as.data.frame(df2 %>% group_by(Fund_Emi_ID) %>% sample_n(1, replace=TRUE))
    }else{
      # df2 <- as.data.frame(df2 %>% group_by(Fund_Emi_ID) %>% filter(FundInv_Quarter == Investment_Date))
      df2 <- as.data.frame(df2 %>% group_by(Fund_Emi_ID) %>% filter(FundInv_Quarter == min(FundInv_Quarter)))
    }
    
    DF <- df2[,reg_vars2]
    
  } else {
    if(private.source == "Preqin"){
      DF <- super$creat_reg_df.preqin(super$stochastic.assigner()$B)
    }
  }
  
  return(DF)
}

# system.time( print(head(multiple$creat_reg_df("BO",stochastic = TRUE))) )


# binomial-gamma hurdle model
# http://seananderson.ca/2014/05/18/gamma-hurdle.html
# https://www.researchgate.net/post/Are_there_generalizations_of_zero-inflated_negative_binomial_and_hurdle_modeling_that_address_continuous_ie_non-count_variables
# https://www.ibm.com/developerworks/library/ba-optimR-john-nash/


# TODO: implement weighting in MLE (use Distance-to-Exit)


multiple$MLE <- function(n = 10, type= "BO", METH = "ucminf", benchmark_optimx = FALSE, 
                         generate.private = "AssetMetrix"){
  
  if(n == 1 | generate.private != "AssetMetrix"){
    predictors_1 <- c("Time2Exit_ZS","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
    predictors_0 <- c("Time2Exit_ZS","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
    predictors_2 <- c("ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
  } else {
    predictors_2 <- c("RVPI_1","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
    predictors_1 <- c("Holding_Period","Time2Exit_ZS","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
    predictors_0 <- c("RVPI_1","Holding_Period","Time2Exit_ZS","ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
  }
  
  response1 <- "P2C.multi1"
  hurdle0 <- 0.1
  hurdle2 <- 5
  
  formula0 <- as.formula(paste("hurdle0", paste(predictors_0, collapse=" + "), sep=" ~ "))
  formula0_IO <- as.formula(paste("hurdle0", 1, sep=" ~ "))
  formula1 <- as.formula(paste(response1, paste(predictors_1, collapse=" + "), sep=" ~ "))
  formula1_IO <- as.formula(paste(response1, 1, sep=" ~ "))
  formula2 <- as.formula(paste("hurdle2", paste(predictors_2, collapse=" + "), sep=" ~ "))
  formula2_IO <- as.formula(paste("hurdle2", 1, sep=" ~ "))
  

# Function to generate start values for Likelihood Maximation
  estimate.start.values <- function(df){
    df$hurdle0 <- ifelse(df[,response1] > hurdle0, 1, 0)
    df$hurdle2 <- ifelse(df[,response1] > hurdle2, 1, 0)
    
    glms <- list()
    # Full Covariate Model
    glms$m0 <- glm(formula0, data = df, family = binomial(link = logit))
    glms$m1 <- glm(formula1, data = subset(df, hurdle0 == 1 & hurdle2 == 0), family = Gamma(link = log))
    glms$m2 <- glm(formula2, data = df, family = binomial(link = logit))
    
    # Intercept Only
    glms$m0_IO <- glm(formula0_IO, data = df, family = binomial(link = logit))
    glms$m1_IO <- glm(formula1_IO, data = subset(df, hurdle0 == 1 & hurdle2 == 0), family = Gamma(link = log))
    glms$m2_IO <- glm(formula2_IO, data = df, family = binomial(link = logit))
    
    esv <- list()
    for(modelz in c("m0", "m1", "m2")){
      m_IO <- glms[[paste(modelz,"IO",sep = "_")]]
      m <- glms[[modelz]]
      
      esv[[modelz]][["Full"]][["Coef"]] <- coef(m)
      esv[[modelz]][["Full"]][["Theta"]] <- c(Theta = summary(m)$dispersion)
      esv[[modelz]][["Full"]][["AIC"]] <- c(AIC = summary(m)$aic)
      
      esv[[modelz]][["IO"]][["Coef"]] <- c(Intercept_IO = as.numeric(coef(m_IO)))
      esv[[modelz]][["IO"]][["Theta"]] <- c(Theta_IO = summary(m_IO)$dispersion)
      esv[[modelz]][["IO"]][["AIC"]] <- c(AIC_IO = summary(m_IO)$aic)
    }
    
    return(esv)
  }
  
# Resampling Procedure
  out_nb2 <- list()
  out_nb2_IO <- list()
  hurdle.estimates <- list()
  hurdle.estimates.n1 <- list()
  huest <- list()
  
  for(i in seq(1, n)){
    print(i)
    
# create correct input data
    if(n == 1){ stoch_df = FALSE }else{ stoch_df = TRUE }
    df_bovc <- multiple$creat_reg_df(f.type = type, stochastic = stoch_df, 
                                     private.source = generate.private) # create random input
    

# PART 1 - Negative Binomial GLM
    create_output <- TRUE
    do.NB2 <- TRUE
    if(do.NB2){
      count_depvar <- paste("I(round(",response1," * ",100,",0))",sep="")
      predictors_nb2 <- c("ML_HYOAS.quarter","MSCI.Multiple.Exit_1")
      formula_NB2 <- as.formula(paste(count_depvar, paste(predictors_nb2, collapse=" + "), sep=" ~ "))
      formula_NB2_IO <- as.formula(paste(count_depvar, 1, sep=" ~ "))
      
      # glm.tweedie <- glm(formula_NB2, family = statmod::tweedie(var.power=1.5,link.power=0), data=df_bovc)
      # print(summary(glm.tweedie))
      
      error_status_glm.nb <- try( nb2_reg <- MASS::glm.nb(formula_NB2, data = df_bovc) )
      
      if(class(error_status_glm.nb)[1] == "try-error"){
        print(paste("NB2: Error occured: Run ",i))
        create_output <- FALSE
        # break
      }else{
        if(nb2_reg$theta > 100){
          print(paste("NB2: Theta too high: Run ",i))
          create_output <- FALSE
        }else{
          out_nb2[[paste("r",i,sep="")]] <- c(nb2_reg$coefficients, 
                                              Theta= nb2_reg$theta, AIC= nb2_reg$aic)
          
          # NB2 intercept only model
          nb2_reg_IO <- MASS::glm.nb(formula_NB2_IO,data=df_bovc)
          out_nb2_IO[[paste("r",i,sep="")]] <- c(nb2_reg_IO$coefficients, 
                                                 Theta= nb2_reg_IO$theta, AIC= nb2_reg_IO$aic)
          
          if(n == 1){
            out_nb2_n1 <- nb2_reg
            out_nb2_n1_IO <- nb2_reg_IO
          }
          
        }
      }
    }

# PART 2 - Hurdle Model
    # m_1
    if(create_output){
      # Estimate Start Values for all models m_0, m_1, m_2
      esv_123 <- estimate.start.values(df_bovc)

      # truncdist::dtrunc
      trunc_dist <- function(x, spec = "gamma", a = hurdle0, b = hurdle2, ...){
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
      ml_trunc1 <- function(par, df = df_bovc, do.IO = FALSE){
        
        df1 <- df[df[, response1] > hurdle0 & df[, response1] < hurdle2, ]
        df1 <- df1[complete.cases(df1),]
        
        Y <- df1[, response1]
        
        if(do.IO){
          form1 <- formula1_IO
        } else {
          form1 <- formula1
        }
        
        mf_1 <- model.frame(formula = form1, data = df1)
        mt_1 <- attr(mf_1, "terms")
        X_1 <- model.matrix(mt_1, mf_1)
        b_1 <- par[1 : ncol(X_1)]
        
        shape1 <- par[ncol(X_1) + 1]
        

        Xb <- X_1 %*% b_1
        mu <- exp(Xb)
        scale1 <- mu / shape1
        
        dtg <- trunc_dist(Y, shape = shape1, scale = scale1)

        neg.log.lik <- - sum(log(dtg))
        
        return(neg.log.lik)
      }

      # Optimize m1
      for(j in seq(1,1)){
        start_m1 <- c(esv_123$m1$Full$Coef,  esv_123$m1$Full$Theta)

        m1_result <- optimx::optimx(start_m1, ml_trunc1, method = METH, do.IO = FALSE)

        # compare various optimx results
        if(benchmark_optimx){ 
          optim_result <- optimx::optimx(start, ml_trunc1,
                                         control = list(all.methods=TRUE, save.failures=TRUE))
          print(optim_result)
        }
        
        # report non-convergence
        if(m1_result$convcode != 0){
          print(paste("*** !!!! --- No convergence m0: Run", i))
        }
        
        # store optimzation results
        huest$m1$Full$Coef[[paste("r", i, sep = "")]] <- unlist(m1_result[1 : (length(start_m1) - 1)])
        huest$m1$Full$Theta[[paste("r", i, sep = "")]] <- unlist(m1_result[length(start_m1)])
        huest$m1$Full$AIC[[paste("r", i, sep = "")]] <- 2 * length(start_m1) - 2 * (-m1_result$value)
        
        
        # store SEs, t- & p-values, if n == 1
        if(n == 1){
          Hessian2Summary <- function(hezzian, parz){
            fisher_info<-solve(--hezzian) # postive hessian since we minimize -log(Likelihood) ??
            prop_sigma<-sqrt(diag(fisher_info))
            
            result <- data.frame(MLE = unlist(parz), SE = prop_sigma)
            result$t_value <- result$MLE / result$SE
            result$p_value <- pnorm(-abs(result$t_value))
            return(result)
          }
          
          hurdle.estimates.n1$m1$Full <- Hessian2Summary(attr(m1_result, "details")[1,"nhatend"][[1]],
                                                c(huest$m1$Full$Coef[[paste("r", i, sep = "")]], 
                                                  Shape_Gamma = huest$m1$Full$Theta[[paste("r", i, sep = "")]]))
          hurdle.estimates.n1$m1$Full["AIC", 1] <- huest$m1$Full$AIC[[paste("r", i, sep = "")]]
          hurdle.estimates.n1$m1$Full[, "Model"] <- "m1"
        }
        
        
        # intercept only model update
        # start_IO <- out_123$m1[c("Intercept_IO","Theta_IO")]
        start_m1_IO <- c(esv_123$m1$IO["Coef"], esv_123$m1$IO["Theta"])
        
        m1_result_IO <- optimx::optimx(start_m1_IO, ml_trunc1, method = METH, do.IO = TRUE)

        # store optimzation results
        huest$m1$IO$Coef[[paste("r", i, sep = "")]] <- unlist(m1_result_IO[1:(length(start_m1_IO) - 1)])
        huest$m1$IO$Theta[[paste("r", i, sep = "")]] <- unlist(m1_result_IO[length(start_m1_IO)])
        huest$m1$IO$AIC[[paste("r", i, sep = "")]] <- 2 * length(start_m1_IO) - 2 * (-m1_result_IO$value)

        
        # store SEs, t- & p-values, if n == 1
        if(n == 1){
          hurdle.estimates.n1$m1$IO <- Hessian2Summary(attr(m1_result_IO, "details")[1,"nhatend"][[1]],
                                              c(Intercept = huest$m1$IO$Coef[[paste("r", i, sep = "")]], 
                                                Shape_Gamma = huest$m1$IO$Theta[[paste("r", i, sep = "")]]))
          hurdle.estimates.n1$m1$IO["AIC", 1] <- huest$m1$IO$AIC[[paste("r", i, sep = "")]]
          hurdle.estimates.n1$m1$IO[, "Model"] <- paste("m1","IO",sep = "_")
        }
      }
      
  }
    # m_0 and m_2 (joint likelihood approach)
    if(create_output ){
      
      # likelihood function for hurdles h_0 and h_2
      ml_hurdles02 <- function(par, df = df_bovc, do.IO = FALSE){
        
        df$hurdle0 <- ifelse(df[, response1] > hurdle0, 1, 0)
        df$hurdle2 <- ifelse(df[, response1] > hurdle2, 1, 0)
        
        Y <- df[, response1]
        
        if(do.IO){
          form0 <- formula0_IO
          form2 <- formula2_IO
        } else {
          form0 <- formula0
          form2 <- formula2
        }

        mf_0 <- model.frame(formula = form0, data = df)
        mt_0 <- attr(mf_0, "terms")
        X_0 <- model.matrix(mt_0, mf_0)
        b_0 <- par[1 : ncol(X_0)]

        
        mf_2 <- model.frame(formula = form2, data = df)
        mt_2 <- attr(mf_2, "terms")
        X_2 <- model.matrix(mt_2, mf_2)
        b_2 <- par[ncol(X_0) + (1 : ncol(X_2))]

        
        # cf. Agresti (2010) 4.1 "adjacent-categories logit models"
        # P(Y < h_0)
        f_0 <- 1 - (1 / (1 + exp(-(X_0 %*% b_0))))
        # P(Y > h_2)
        f_2 <- (1 - f_0) * (1 / (1 + exp(-(X_2 %*% b_2))))

        
        likelihood <- ifelse(Y <= hurdle0, f_0,
                             ifelse(Y > hurdle2, f_2,
                                    1 - f_0 - f_2))
        
        # check 1 - f_0 - f_2 >= 0
        if(min(likelihood) < 0){
          print("Issue m0 & m2:  (1 - f_0 - f_2) >= 0")
          likelihood <- pmax(0.0001,likelihood)
        }
        
        neg.log.lik <- -sum(log(likelihood))

        return(neg.log.lik)
      }

      # Optimize m0 and m2
      for(j in seq(1,1)){
        
        # A) Full Parameter Model
        if(TRUE){
          start_hurdles02 <- c(esv_123$m0$Full$Coef, esv_123$m2$Full$Coef)

          m02_result <- optimx::optimx(start_hurdles02, ml_hurdles02, method = METH, do.IO = FALSE)
          
          # compare solnp_result to other opimization methods
          if(benchmark_optimx){
            optim_result_02 <- optimx::optimx(par_hurdles02, ml_hurdles02,
                                              control = list(all.methods = TRUE, save.failures = TRUE))
            print(optim_result_02)
          }
          
          # report non-convergence
          if(m02_result$convcode != 0){
            print(paste("*** !!!! --- No convergence m0 & m2: Run", i))
          }
          
          # store optimzation results
          huest$m0$Full$Coef[[paste("r", i, sep = "")]] <- unlist(m02_result[1 : length(esv_123$m0$Full$Coef)])
          huest$m2$Full$Coef[[paste("r", i, sep = "")]] <- unlist(m02_result[length(esv_123$m0$Full$Coef) + (1 : length(esv_123$m2$Full$Coef))])
          
          AIC_m02 <- 2 * length(start_hurdles02) - 2 * (-m02_result$value)
          huest$m0$Full$AIC[[paste("r", i, sep = "")]] <- AIC_m02
          huest$m2$Full$AIC[[paste("r", i, sep = "")]] <- AIC_m02
          huest$m0$Full$Theta[[paste("r", i, sep = "")]] <- NA
          huest$m2$Full$Theta[[paste("r", i, sep = "")]] <- NA
          
          if(n == 1){
            # m0 & m2
            hurdle.estimates.n1$m02$Full <- Hessian2Summary(attr(m02_result, "details")[1,"nhatend"][[1]],
                                                   c(huest$m0$Full$Coef[[paste("r", i, sep = "")]], 
                                                     huest$m2$Full$Coef[[paste("r", i, sep = "")]]))
            hurdle.estimates.n1$m02$Full["AIC",1] <- AIC_m02
            hurdle.estimates.n1$m02$Full[,"Model"] <- c(rep("m0",length(huest$m0$Full$Coef[[paste("r", i, sep = "")]])), 
                                                        rep("m2",length(huest$m2$Full$Coef[[paste("r", i, sep = "")]])), "m02")
          }
        }
        
        # B) Intercept Only Model
        if(TRUE){
          start_hurdles02_IO <- c(esv_123$m0$IO["Coef"], esv_123$m2$IO["Coef"])
          
          m02_result_IO <- optimx::optimx(start_hurdles02_IO, ml_hurdles02, method = METH, do.IO = TRUE)
          
          # report non-convergence
          if(m02_result_IO$convcode != 0){
            print(paste("*** !!!! --- No convergence m0 & m2 IO: Run", i))
          }
          
          # store optimzation results
          huest$m0$IO$Coef[[paste("r", i, sep = "")]] <- unlist(m02_result_IO[1])
          huest$m2$IO$Coef[[paste("r", i, sep = "")]] <- unlist(m02_result_IO[2])
          
          AIC_m02_IO <- 2 * length(start_hurdles02_IO) - 2 * (-m02_result_IO$value)
          huest$m0$IO$AIC[[paste("r", i, sep = "")]] <- AIC_m02_IO
          huest$m2$IO$AIC[[paste("r", i, sep = "")]] <- AIC_m02_IO
          huest$m0$IO$Theta[[paste("r", i, sep = "")]] <- NA
          huest$m2$IO$Theta[[paste("r", i, sep = "")]] <- NA
          
          
          if(n == 1){
            # m0 & m2
            hurdle.estimates.n1$m02$IO <- Hessian2Summary(attr(m02_result_IO, "details")[1,"nhatend"][[1]],
                                                 c(Intercept_0 = huest$m0$IO$Coef[[paste("r", i, sep = "")]], 
                                                   Intercept_2 = huest$m2$IO$Coef[[paste("r", i, sep = "")]]))
            hurdle.estimates.n1$m02$IO["AIC",1] <- AIC_m02_IO
            hurdle.estimates.n1$m02$IO[,"Model"] <- paste(c(rep("m0", 1),rep("m2", 1),"m02"),"IO",sep = "_")
          }
        }

      }
      
    }
}
  
# Prep Output: NB2
  if(do.NB2){
    out_nb2 <- data.frame(do.call(rbind,out_nb2))
    out_nb2_IO <- data.frame(do.call(rbind,out_nb2_IO))
    colnames(out_nb2_IO) <- paste(colnames(out_nb2_IO),"IO",sep="_")
    
    if(n == 1){
      prep_nb2_n1 <- function(nb2_reg, tag){
        sum.nb2.n1 <- summary(nb2_reg)
        df <- data.frame(coef(sum.nb2.n1))
        colnames(df) <- c("MLE", "SE", "t_value", "p_value")
        df["Theta", 1:2] <- c(nb2_reg$theta,nb2_reg$SE.theta)
        df["Theta", 3] <- df["Theta", 1] / df["Theta", 2]
        df["Theta", 4] <- pnorm(-abs(df["Theta", 3]))
        df["AIC", 1] <- AIC(nb2_reg)
        df$Model <- tag # 
        return(df)
      }
      
      # Full NB2 n == 1 model
      df.nb2.n1 <- prep_nb2_n1(out_nb2_n1, "NB2")
      # IO NB2 n == 1 model
      df.nb2.n1_IO <- prep_nb2_n1(out_nb2_n1_IO, "NB2_IO")
      
      # merge
      df.nb2.n1 <- data.frame(rbind(df.nb2.n1, df.nb2.n1_IO))
    }
  }

# Prep Output: Hurdle Model
  for(model in c("m0", "m1", "m2")){
      huest[[model]]$Full$Coef <- data.frame(do.call(rbind, huest[[model]][["Full"]]$Coef))
      huest[[model]]$Full$Coef$Shape <- huest[[model]][["Full"]]$Theta
      huest[[model]]$Full$Coef$AIC <- huest[[model]][["Full"]]$AIC
      huest[[model]]$Full <- huest[[model]][["Full"]]$Coef
      
      huest[[model]]$IO <- data.frame(Intercept = huest[[model]]$IO$Coef,
                                      Shape = huest[[model]]$IO$Theta,
                                      AIC = huest[[model]]$IO$AIC)
  }

# Define Final Output
  if(do.NB2){
    if(n == 1){
      final_out <- list(NB2 = out_nb2,  NB2_n1 = df.nb2.n1,  
                        HM = huest, HMn1 = hurdle.estimates.n1)
    } else {
      final_out <- list(NB2 = out_nb2, NB2_IO = out_nb2_IO, 
                        HM = huest)
    }
  } else {
    if(n == 1){
      final_out <- list(HM = huest, HMn1 = hurdle.estimates.n1)
    } else {
      final_out <- list(HM = huest)
    }
  }

  return(final_out)
}



multiple$CreateOutput <- function(list_in, make.csv = FALSE){
# new core function
  create.stats <- function(input, reg_name){
    
    modelz <- c("m0", "m1", "m2")
    
    # detect Inf_ids 
    Inf_ids <- list()
    for(model in modelz){
      Inf_ids[[paste("Full", model, "AIC")]] <- which(is.infinite(input$HM[[model]]$Full$AIC))
      Inf_ids[[paste("IO", model, "AIC")]] <- which(is.infinite(input$HM[[model]]$IO$AIC)) 
    }
    Inf_id <- unlist(Inf_ids)
    print(Inf_id)

    
    
    # prep HM
    hm.models <- list()
    for(model in modelz){
      
      para.stats.list <- list()
      for(para in c("Full", "IO")){
        df.data <- input$HM[[model]][[para]]
        
        if(length(Inf_id) > 0){
          df.data <- df.data[- Inf_id, ]
        }
        
        for(stats in c("mean", "sd")){
          def <- data.frame(
                     Value = apply(df.data, 2, stats),
                     Model = model, 
                     Para = para,
                     Stats = stats)
          def$Variable <- rownames(def)
          
          # Correct Variable Names
          def$Variable <- ifelse(grepl("Intercept",def$Variable), "Intercept", def$Variable)
          def$Variable <- ifelse(grepl("ML_HYOAS",def$Variable), "HYS", def$Variable)
          def$Variable <- ifelse(grepl("MSCI",def$Variable), "PEM-1", def$Variable)
          def$Variable <- ifelse(grepl("AIC",def$Variable), "AIC", def$Variable)
          def$Variable <- ifelse(grepl("Time2Exit_ZS",def$Variable), "T2E_ZS", def$Variable)
          def$Variable <- ifelse(grepl("Theta",def$Variable), "Theta", def$Variable)
          def$Variable <- ifelse(grepl("Shape",def$Variable), "Shape", def$Variable)
          def$Variable <- ifelse(grepl("RVPI",def$Variable), "RVPI-1", def$Variable)
          
          def$Variable <- ifelse(grepl("Holding_Period",def$Variable), "HP", def$Variable)
        
          para.stats.list[[paste(para, stats, sep = "_")]] <- def
        }

      }
      def2 <- data.frame(do.call(rbind, para.stats.list))
      def2[, model] <- def2$Value
      rownames(def2) <- NULL
      def2 <- def2[, c("Variable", "Stats", "Para", model)]
      hm.models[[model]] <- def2
    }
    
    hm.out <- Reduce(function(x, y) merge(x, y, all=T, 
                                                   by = c("Variable", "Stats", "Para")), 
                              hm.models, accumulate=F)
    hm.out <- round_df(hm.out)
    
    # Formatting SD (in parentheses)
    for(colz in modelz){
      hm.out[, colz] <- ifelse(hm.out$Stats == "mean", hm.out[, colz],
                               ifelse(is.na(hm.out[, colz]), NA,
                               paste("(", hm.out[, colz], ")", sep = "")))
    }
    
    hm.out <- hm.out[order(hm.out$Para, 
                           factor(hm.out$Variable, levels = c("Intercept","HP","T2E_ZS","RVPI-1","PEM-1","HYS","Shape","Scale (NB2)", "AIC")),
                           hm.out$Stats), ]
    return(hm.out)
  }
  
# define core function
  create.function <- function(input, reg_name){
    # Static or Dynamic Model
    n <- nrow(input$NB2)
    
    if (n > 1) {
      # Filter Inf data
      Inf_nb <- which(is.infinite(input$NB2[, "AIC"]))
      Inf_nb_IO <- which(is.infinite(input$NB2_IO[, "AIC_IO"]))
      Inf_hm <- which(is.infinite(input$HM[2, 11, ]))
      Inf_IDs <- c(Inf_nb, Inf_hm, Inf_nb_IO)
      
      if (length(Inf_IDs) > 0) {
        in_hm <- input$HM[ , ,-Inf_IDs]
        in_nb <- input$NB2[-Inf_IDs, ]
        # If n > 1: we need to calculate the intercept only NB2 statistics
        in_nb_IO <- input$NB2_IO[-Inf_IDs, ]
      } else {
        in_hm <- input$HM
        in_nb <- input$NB2
        # If n > 1: we need to calculate the intercept only NB2 statistics
        in_nb_IO <- input$NB2_IO
      }
      
      
    } else {
      in_hm <- input$HM
      in_nb <- input$NB2
    }
    
    # in_hm <- test.tph$Reg1$HM
    
    
    out_hm <- list() ; out_nb <- list() ; out_nb_IO <- list()
    for(statics in c("mean","median","sd")){
      out_nb[[statics]] <- data.frame(apply(in_nb, 2, statics))
      out_hm[[statics]] <- data.frame(apply(in_hm$Full,  2, statics))
      if(n > 1){
        out_nb_IO[[statics]] <- data.frame(apply(in_nb_IO,  2, statics))
      }
    }
    print(out_hm)
    
    # modified hurdle model output
    out_hm$t_value <- out_hm$mean / out_hm$sd
    out_hm$p_value <- apply(out_hm$t_value, c(1, 2), function(x) pnorm(-abs(x))) 
    modified_out_hm <- lapply(out_hm, function(x){
      df <- as.data.frame(t(x))
      df <- reshape(df, idvar = "Variable", ids = rownames(df), 
                    times = colnames(df), timevar = "Model",
                    varying = list(colnames(df)), direction = "long")
      return(df)
    })
    modified_out_hm <- Reduce(function(x, y) merge(x, y, all=T, 
                                                   by=c("Variable","Model")), 
                              modified_out_hm, accumulate=F)
    modified_out_hm <- data.frame(modified_out_hm)
    colnames(modified_out_hm)[3:7] <- c("Mean","Median","SD","t_value","p_value")
    modified_out_hm <- modified_out_hm[order(modified_out_hm$Model, modified_out_hm$Variable),]

    modified_out_hm <- round_df(modified_out_hm)
    
    
    # NB2 Model (resampling)
    prepare_out <- function(out){
      out <- data.frame(do.call(cbind,out))
      colnames(out) <- c("Mean","Median","SD")
      out$t_value <- out$Mean / out$SD
      out$p_value <- pnorm(-abs(out$t_value))
      out$Variable <- rownames(out)
      out$Model <- "NB2"
      return(round_df(out))
    }
    modified_out_nb <- prepare_out(out_nb)
    if(n > 1){
      modified_out_nb_IO <- prepare_out(out_nb_IO)
    }
    
    # Format long output
    format_long_resampling <- function(outp){
      outp$IO <- ifelse(grepl("IO",outp$Variable), "IO", "Full")
      
      outp$Variable <- ifelse(grepl("Intercept",outp$Variable), "Intercept", outp$Variable)
      outp$Variable <- ifelse(grepl("ML_HYOAS",outp$Variable), "HYS", outp$Variable)
      outp$Variable <- ifelse(grepl("MSCI",outp$Variable), "PEM_1", outp$Variable)
      outp$Variable <- ifelse(grepl("AIC",outp$Variable), "AIC", outp$Variable)
      outp$Variable <- ifelse(grepl("Time2Exit_ZS",outp$Variable), "T2E_ZS", outp$Variable)
      outp$Variable <- ifelse(grepl("Theta",outp$Variable), "Theta", outp$Variable)
      outp$Variable <- ifelse(grepl("Holding_Period",outp$Variable), "HP", outp$Variable)
      
      outp <- outp[complete.cases(outp),]
      # outp <- outp[,c("Model","IO","Variable","Mean","Median","SD","t_value","p_value")]
      outp <- outp[,c("Model","IO","Variable","Mean","SD","t_value","p_value")]
      outp <- outp[order(outp$Model, outp$IO, outp$Variable), ]
      rownames(outp) <- NULL
      
      return(outp)
    }
    format_long_n1 <- function(outp){
      
      outp$Variable <- rownames(outp)
      outp$IO <- ifelse(grepl("IO",outp$Model), "IO", "Full")
      
      outp$Variable <- ifelse(grepl("Intercept",outp$Variable), "Intercept", outp$Variable)
      outp$Variable <- ifelse(grepl("ML_HYOAS",outp$Variable), "HYS", outp$Variable)
      outp$Variable <- ifelse(grepl("MSCI",outp$Variable), "PEM_1", outp$Variable)
      outp$Variable <- ifelse(grepl("AIC",outp$Variable), "AIC", outp$Variable)
      outp$Variable <- ifelse(grepl("Time2Exit_ZS",outp$Variable), "T2E_ZS", outp$Variable)
      outp$Variable <- ifelse(grepl("Theta",outp$Variable), "Theta", outp$Variable)
      outp$Variable <- ifelse(grepl("Holding_Period",outp$Variable), "HP", outp$Variable)
      
      for(model in c("m1","m02","NB2")){
        outp$Model <- ifelse(grepl(model,outp$Model), model, outp$Model)
      }
      outp$Model <- ifelse(outp$Model == "m0_IO", "m0", outp$Model)
      outp$Model <- ifelse(outp$Model == "m2_IO", "m2", outp$Model)
      
      outp <- outp[,c("Model","IO","Variable","MLE","SE","t_value","p_value")]
      outp <- outp[order(outp$Model, outp$IO, outp$Variable), ]
      rownames(outp) <- NULL
      
      return(outp)
    }
    
    
    # http://hydroecology.net/iterating-through-lists-of-lists-of-lists/
    magicCSV.round = function(startlist, to.csv = FALSE, targetdir=NA){
      # use ff function to get all list name paths
      ff = function(x){ 
        if (class(x) == "list") 
          lapply(x, ff) 
        else if(class(x) == "data.frame") 
          TRUE
        else
          NULL
      }
      lnames = names(unlist(lapply(startlist, ff)))
      varnames = strsplit(lnames, split = ".", fixed = TRUE)
      
      # create .csv paths
      fnames = file.path(file.path(targetdir), paste0(gsub(".", "/", 
                                                           lnames, fixed = TRUE), ".csv"))
      
      # create directories
      if(to.csv){
        dirnames = unlist(lapply(fnames, dirname))
        suppressWarnings(lapply(dirnames, dir.create, recursive = TRUE))
      }
      
      # apply round_df function to each dataframe
      list_df <- lapply(seq_along(varnames), function(i){
        df <- startlist[[varnames[[i]]]]
        df <- round_df(df)
        if(to.csv){
          write.csv(df, file = fnames[[i]])
        }
        return(df)
      })
      
      out_df <- as.data.frame(do.call(rbind,list_df))
      return(out_df)
    }
    
    
    
    if("HMn1" %in% names(input)){
      OUTPUT <- list(NB = input$NB2_n1, HM = input$HMn1)
      OUTPUT2 <- magicCSV.round(OUTPUT)
      OUT <- format_long_n1(OUTPUT2)
    }else{
      OUTPUT <- as.data.frame(do.call(rbind,list(modified_out_nb ,modified_out_hm, modified_out_nb_IO)))
      OUT <- format_long_resampling(OUTPUT)
    }
    
    # Re-Name Theta correctly
    OUT$Variable <- ifelse(OUT$Variable == "Theta",
                           ifelse(OUT$Model == "m1", "Shape (Gamma)", "Scale (NB2)"),
                           OUT$Variable)
    
    return(OUT)
  }

# nice formating function
  MS_format <- function(df){
    is.resampling <- any(colnames(df) %in% "Mean")
    if(is.resampling){
      v.nameZ <- c("Mean", "SD")
    }else{
      v.nameZ <- c("MLE", "SE")
    }

    re.shaped <- list()
    for(IO_type in c("Full", "IO")){
      for(v.name in v.nameZ){
        drop.name <- v.nameZ[v.nameZ != v.name]
        df2 <- reshape(data = df[df$IO == IO_type,], direction = "wide", 
                       v.names = v.name, idvar = "Variable", timevar = "Model", 
                       drop = c(drop.name, "t_value", "p_value"))
        colnames(df2) <- gsub(v.name, "Coef", colnames(df2))
        df2$M_S <- v.name
        
        if(v.name %in% c("SD", "SE")){
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
        
        re.shaped[[paste(IO_type, v.name, sep="_")]] <- df2
      }
    }
    df3 <- data.frame(do.call(rbind, re.shaped))
    df3 <- df3[order(df3$IO, 
                     factor(df3$Variable, levels = c("Intercept","HP","T2E_ZS","RVPI_1","PEM_1","HYS","Shape (Gamma)","Scale (NB2)", "AIC"))
    ), ]
    
    if(!is.resampling){
      m02AIC.full <- df3$Coef.m02[df3$Variable == "AIC" & df3$M_S == "MLE" & df3$IO == "Full"]
      m02AIC.IO <- df3$Coef.m02[df3$Variable == "AIC" & df3$M_S == "MLE" & df3$IO == "IO"]
      
      df3$Coef.m0[df3$Variable == "AIC" & df3$M_S == "MLE" & df3$IO == "Full"] <- m02AIC.full
      df3$Coef.m0[df3$Variable == "AIC" & df3$M_S == "MLE" & df3$IO == "IO"] <- m02AIC.IO
      df3$Coef.m2[df3$Variable == "AIC" & df3$M_S == "MLE" & df3$IO == "Full"] <- m02AIC.full
      df3$Coef.m2[df3$Variable == "AIC" & df3$M_S == "MLE" & df3$IO == "IO"] <- m02AIC.IO
    }
    
    df3 <- df3[, c("IO", "Variable", "M_S", paste("Coef",c("m0", "m1","m2", "NB2"), sep = "."))]
    rownames(df3) <- NULL
    
    return(df3)
  }

# iterate over list
  out_list <- list()
  for(reg_names in names(list_in)){
    input1 <- list_in[[reg_names]]
    
    OUT <- create.stats(input1, reg_names)
    # OUT <- create.function(input1, reg_names)
    out_list[[reg_names]] <- OUT

    # OUT <- MS_format(OUT)
    
    if (make.csv) {
      csv.dir <- paste(Sys.Date(), "Multiple", reg_names, "iter.csv", sep = "_")
      write.csv2(OUT, file = csv.dir, na = "", row.names = FALSE)
    } else {
      print(paste(">>> ---- <<<", reg_names))
      print(OUT)
    }
    
  }
  
  invisible(out_list)
}




multiple$predict_dgh <- function(df_in= G.p2c, 
                                 U.multiple = NA,
                                 model_in= list(BO = multi.reg.sum$BO_500, 
                                                VC= multi.reg.sum$VC_500),
                        sim_in = G.p2c[sample(seq(nrow(G.p2c)),1),],
                        sim_out = FALSE, plot_it = FALSE, eps = FALSE, size = "big"){
  model_in$BO <- model_in$BO[model_in$BO$Para == "Full" & model_in$BO$Stats == "mean", ]
  model_in$VC <- model_in$VC[model_in$VC$Para == "Full" & model_in$VC$Stats == "mean", ]

  
  h0 <- 0.1 ; h2 <- 5
  
  ecdf_bo_h0 <- ecdf(df_in$P2C.multi1[df_in$Fund_InvestTypes == "BO" & df_in$P2C.multi1 < h0])
  ecdf_bo_h2 <- ecdf(df_in$P2C.multi1[df_in$Fund_InvestTypes == "BO" & df_in$P2C.multi1 > h2])
  ecdf_vc_h0 <- ecdf(df_in$P2C.multi1[df_in$Fund_InvestTypes == "VC" & df_in$P2C.multi1 < h0])
  ecdf_vc_h2 <- ecdf(df_in$P2C.multi1[df_in$Fund_InvestTypes == "VC" & df_in$P2C.multi1 > h2])
  
  # linear estimates & gamma shapa parameter
  linear.prediction <- function(type1, x.vec){
    mo.frame <- model_in[[type1]]

    out <- list()
    for(m in c("m0", "m1", "m2")){
      df <- mo.frame[, c("Variable", m)]
      colnames(df) <- c("Variable", "Beta")
      df <- merge(df, data.frame(Variable = names(x.vec), X = x.vec), by = "Variable", all.y = TRUE)

      betaX <- sum(as.numeric(df$Beta) * df$X, na.rm = TRUE)
      out[[m]] <- betaX
    }
    
    out[["shape"]] <- as.numeric(mo.frame[mo.frame$Variable == "Shape", "m1"])
    return(out)
  }

  
  output <- list()
  df_in <- if(sim_out)  sim_in else df_in
  runo <- nrow(df_in)
  if(is.na(U.multiple)){
    U.multiple <- runif(runo)
  }
  iter <-   seq(runo)
  for(i in iter){
    real_multiple <- df_in$P2C.multi1[i]
    
    type= as.character(df_in$Fund_InvestTypes[i])
    Xi <- c(Intercept = 1,
            'RVPI-1' = df_in$RVPI_1[i],
            HP = df_in$Holding_Period[i],
            T2E_ZS = df_in$Time2Exit_ZS[i],
            HYS = df_in$ML_HYOAS.quarter[i],
            'PEM-1' = df_in$MSCI.Multiple.Exit_1[i])
    
    betaXshape <- linear.prediction(type, Xi)

    p0 <- 1 - plogis( betaXshape[["m0"]] ) # P(Y < hurdle0)
    p2 <- (1 - p0) * plogis( betaXshape[["m2"]] )    # P(Y > hurdle2)
    p1 <- 1 - p0 - p2    # P(Y <= hurdle2 & Y >= hurdle0)
    mu1 <- exp( betaXshape[["m1"]] )
    shape1 <- betaXshape[["shape"]]
    scale1 <- mu1/shape1
    
    if(plot_it){
      print(Xi)
      P <- c(P0= p0,P1= p1,P2= p2)
      print(P)
    }
    compound_CDF <- function(Y, tyype){
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
      output[i] <-  inverse_CDF(U.multiple[i])
    }else{
      output[i] <- compound_CDF(real_multiple,type)
    }
    
    if(plot_it){
      if(eps){
        setEPS()
        if(size == "small"){
          postscript("Compound_CDF_DGH_small.eps", 
                     width = 3.5, height = 2.5, 
                     family = "Helvetica", pointsize = 2)
        } else {
          postscript("Compound_CDF_DGH.eps", 
                     width = 5, height = 2.5, 
                     family = "Helvetica",pointsize = 5)
        }
        
      }
      # BO
      par(cex=1.3)
      curve(ecdf_bo_h0, 0, h0, xlim = c(0, 10), ylim = c(0, 1), ylab = "CDF", xlab = "Multiple (Y)", col = "red")
      abline(v = c(h0, h2), col = "gray", lty = 2) ; abline(h = c(0, 1), col = "gray")
      curve(ecdf_bo_h2, add = TRUE, col = "green", h2, 10)
      curve(truncdist::ptrunc(x, spec = "gamma", shape = shape1, scale = scale1, a = h0, b = h2), h0, h2, col = "blue", add = TRUE)
      
      seq_0 <- seq(0,h0,0.01) ; seq_1 <- seq(h0,h2,0.01) ; seq_2 <- seq(h2,10,0.01)
      points(seq_0,compound_CDF(seq_0,type),col="red",type="l",lwd=2)
      points(seq_1,compound_CDF(seq_1,type),col="blue",type="l",lwd=2)
      points(seq_2,compound_CDF(seq_2,type),col="green",type="l",lwd=2)
      leg.cex <- 1.2
      if(size == "small") leg.cex <- 1
      legend("bottomright",bty="n",col=c("green","blue","red"), lty = 1,cex = leg.cex, lwd = 2,
             legend=c(latex2exp::TeX('$\\Y > h_{2} $'),
                      latex2exp::TeX('$\\h_{0} < Y \\leq h_{2} $'),
                      latex2exp::TeX('$\\Y \\leq h_{0} $')
             ))
      #legend(x=6,y=0.5,bty="n",legend=paste(names(Xi[-1]),round(Xi[-1],2)),cex=0.6,text.col="darkgrey")
      # text(x=0,y=1.08,round(p0,3),cex=0.8,col= "red")
      # text(x=1.7,y=1.08,round(p1,3),cex=0.8, col= "blue")
      # text(x=6,y=1.08,round(p2,3),cex=0.8, col="green")
      
      if(eps){ dev.off() }
    }
    
  }
  output <- as.numeric(unlist(output))
  return(output)
}



multiple$DGH_Roseblatt <- function(df, uni, eps = FALSE){
  if(eps){
    setEPS() ; postscript("Hurdle_Rosenblatt.eps", width = 5, height = 2.5, family = "Helvetica",pointsize=5)
  }
  bo_col <- "royalblue1" ; vc_col <- "maroon1"
  
  df$RQ_Multiple.dgh <- uni
  
  par(mar=c(3,3,2,1),mfrow=c(1,2),oma=c(3,2,1,1),cex=1.2)
  for(type in c("BO","VC")){
    hi_col <- ifelse(type=="BO",bo_col,vc_col)
    hist(df$RQ_Multiple.dgh[df$Fund_InvestTypes == type],freq = FALSE,
         main=type,border=hi_col,lty=3,ylab=NA,xlab=NA,ylim=c(0,1.7),xlim=c(0,1))
    abline(a=0,b=1,col="darkgray",lwd=1,lty=2)
    wb_ecdf <- ecdf(df$RQ_Multiple.dgh[df$Fund_InvestTypes == type])
    curve(wb_ecdf,add=TRUE,col="forestgreen",lwd=1)
    legend("topright",bty="n",cex=1.3,col=c("forestgreen","darkgray"),legend = c("empirical","theoretical"),lty=c(1,2))
    abline(h=1,col="black",lty=2)
  }
  mtext(latex2exp::TeX('$\\F_{0,1,2}^{(TPH)}(Y_i)}'),side=1,outer= TRUE,line=1,cex=1.7)
  mtext("ECDF or Density",side=2,outer= TRUE,line=0.5,cex=1.5)
  
  
  if(eps){ dev.off() }
}



## x) Attach new environment -----
#while("multiple" %in% search())
#  detach("multiple")
#attach(multiple)