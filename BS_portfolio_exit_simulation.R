##################### ### ### ### ### ###
## Simulate B&S Exits
##################### ### ### ### ### ###
# Prep Data for simulate$Timing.Simulator() --------

df.bs <- g.p2c2[g.p2c2$FundInv_Quarter == "2016-12-31", ]
# df.bs <- g.p2c2[g.p2c2$FundInv_Quarter == "2000-12-31", ]
df.bs$CompanyAge <- as.numeric(df.bs$FundInv_Quarter - df.bs$Investment_Date) / 365.25
df.bs$Type <- df.bs$Fund_InvestTypes
df.bs$TVPI <- df.bs$MOIC.dyn
df.bs$RVPI <- df.bs$RVPI_1 + 1
df.bs$CompanyStartDate <- df.bs$Investment_Date
df.bs$FundAgeAtEntry <- df.bs$FundAgeAI
df.bs <- df.bs[, colnames(df.bs) %in% colnames(simulate$create.sim.input())]


# df.bs <- simulate$Timing.Simulator(df.bs, simulate$create.public.scenario())

# Simulation function  ------
BS.Simulator <- function(iterations){
  all_out <- list()

  # fixed private portfolio
  df.companies <- df.bs
  all_out$PF <- df.companies

  copula.selection <- c("independent", "t.copula",
                        "bottom-left", "bottom-right",
                        "top-left", "top-right")
  copula.selection <- "independent"
  
  df.finalMultiple <- data.frame(Iter = seq(iterations))
  df.finalMultiple[, copula.selection] <- NA
  
  for(i in seq(iterations)){
      if(i %% 10 == 0) print(i)
      
      # public scenario
      pub.sce <- simulate$create.public.scenario()
      
      for(corner in copula.selection){
        # Copula
        U.copula <- simulate$TM.copula(nrow(df.companies), corner)
        
        # Timing 
        df.out <- simulate$Timing.Simulator(df.pr = df.companies,
                                            df.pu = pub.sce,
                                            U.timing = U.copula$Timing)
        
        # Multiple
        df.out$Multiple <- multiple$predict_dgh(sim_out = TRUE,
                                                U.multiple = U.copula$Multiple,
                                                sim_in = df.out)
        
        # Cash Flow on Current NAV
        df.finalMultiple[i, corner] <- sum(df.out$RVPI / sum(df.out$RVPI) * df.out$Multiple)
      }
      
    }
  
  all_out$Output <- df.finalMultiple
  
  
  return(all_out)
}

system.time(bssim1 <- BS.Simulator(500))

sapply(bssim1$Output, sun)
sapply(bssim2$Output, sun)

# Backtest with Preqin Perf Analyst Data ---------
root <- "C:/Users/christian.tausch/Dropbox/Project D/3_R_Pro_D/"
setwd(file.path(root, "Persistence", "Data"))
ppa <- read.csv2("Preqin_PerfAna_2016Q4.csv")
# setwd(file.path(root, "Superposition", "Data"))
# cat.merge <- read.csv2("MergeCatTypesV05.csv")
# ppa <- merge(ppa, cat.merge[ ,c(1,2)], by = )

ppa1 <- ppa[ppa$Status == "Liquidated" & ppa$Type == "Fund of Funds", ]
ppa1$Net.Multiple..X. <- as.numeric(as.character(ppa1$Net.Multiple..X.))
sun(ppa1$Net.Multiple..X.)

pf.list <- list()
for(i in 1:1000){
  ppa1$ID <- sample(rep(seq(1, 33), 3), 99)
  pf.list[[i]] <- aggregate(ppa1$Net.Multiple..X., by = list(ppa1$ID), FUN = mean)$x
}
sun(unlist(pf.list))


sun(ppa1$Net.Multiple..X.)



