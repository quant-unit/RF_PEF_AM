################################
## Exit Dynamics: w/ Preqin Data
################################
## 0) Prologue -------
rm(list=ls()) # remove workspace objects
computer <- "mac"

if(computer == "win") root <- file.path("C:/Users/christian.tausch", "Dropbox", "Project D", "3_R_Pro_D")
if(computer == "mac") root <- file.path("~", "Dropbox", "Project D", "3_R_Pro_D")

wd <- list()
wd$code <- file.path(root, "Exit_Dynamics", "Code")
wd$data <- file.path(root, "Exit_Dynamics", "Code", "Data")
wd$super$code <- file.path(root, "Superposition", "Code")
rm(root,computer)


setwd(wd$code)
source("Useful_Functions.R")
source("ExitDynamics_Packages.R")

setwd(wd$data)
load("ExitDynamics_Data_V1.RData")


## X) Fund Level Data --------
setwd(wd$super$code)
source("Superposition_Main.R")

set.seed(99)
preq.al <- super$creat_reg_df.preqin(super$stochastic.assigner()$B)
head(preq.al)
## X) Timing   ------
setwd(wd$code)
source("ExitDynamics_Timing.R")

system.time(
  wb.cox2 <- timing$ParaCoxRegression(public.data, preq.al[preq.al$InvDate != preq.al$ExitDate, ][1:50, ])
)


## X) Multiple --------
setwd(wd$code)
source("ExitDynamics_Multiple.R")
system.time( multi.test <- multiple$MLE(5, super$fund.type, generate.private = "Preqin"))
multi.test.sum <- multiple$CreateOutput(list(A = multi.test))

