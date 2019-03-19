

library(gamlss)

x <- seq(0,1,0.1)
pGG(x, mu = 1, sigma = 1, nu = 1)
pGA(x, mu = 1, sigma = 1)

setwd("C:/Users/christian.tausch/Dropbox/Project D/3_R_Pro_D/Exit_Dynamics/Code/CSV_output")

GA.list <- readRDS("GA.list20181023.RDS")

round(GA.list$VC$VC$coefs$m0$mu$Mean, 3)
round(GA.list$VC$VC$coefs$m0$mu$SD, 3)
round(GA.list$VC$VC$coefs$m1$mu$Mean, 3)
round(GA.list$VC$VC$coefs$m1$mu$SD, 3)
round(GA.list$VC$VC$coefs$m1$sigma$Mean, 3)
round(GA.list$VC$VC$coefs$m1$sigma$SD, 3)

GA.list$VC$VC$Iter
