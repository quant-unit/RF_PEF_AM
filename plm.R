## plm 
require(plm)

#View(G.p2c)
response <- "P2C.multi1.discounted"
re <- "re(random=~1|Fund_Emi_ID, correlation=corARMA(p=1, q=0), opt='optim')"
predictors <- c(re, "RVPI_1","Holding_Period","ML_HYOAS.quarter", "MSCI.Multiple.Exit_1", "I(CMA.Index.Multi-1)")
formula <- as.formula(paste(response, paste(predictors, collapse=" + "), sep=" ~ "))
predictors_1sigma <- c(re, "Holding_Period", "Time2Exit")
formula1sigma <- as.formula(paste("", paste(predictors_1sigma, collapse=" + "), sep=" ~ "))

f.type <- "BO"
rm(m1.full)
m1.full <- gamlss(P2C.multi1.discounted ~ MSCI.Multiple.Exit_1,
                  #+ re(random=~MSCI.Multiple.Exit_1|Fund_Emi_ID, correlation=corARMA(p=1,q=0),opt='optim'), 
                  #sigma.formula = formula1sigma, 
                  data = subset(G.p2c[, !(colnames(G.p2c) %in% c("Ev", "Ev2", "MOIC.dyn", "MOIC.end"))], 
                                P2C.multi1.discounted > 0 & Fund_InvestTypes == f.type), family = GA)
summary(m1.full)


# PLM
E <- pdata.frame(G.p2c[G.p2c$Fund_InvestTypes == f.type & G.p2c[, response] > 0, ], 
                 index=c("Fund_Emi_ID","FundInv_Quarter"), drop.index=TRUE, row.names=TRUE)
m <- plm(formula, data = E, model = "within")
summary(m)
?re
?random


