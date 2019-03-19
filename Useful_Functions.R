########################
## Some Useful Functions
########################
## define useful functions  -------
# Summarize Univariate distributioN
sun <- function(X,cdf=0,plot.it=TRUE,round.to=3,fit.sgt=FALSE,
                Xlim=range(X),Ylim=c(0,1),main=NA){
  options(warn = -1)
  X <- as.numeric(X)
  nobs <- length(X)
  NAs <- sum(is.na(X))
  Infs <- sum(!is.finite(X) & !is.na(X))
  if(NAs+Infs == length(X)){
    stop("No Numeric Input Data")
  }else{ 
    X <- X[!is.na(X) & is.finite(X)]
    Qs <- quantile(X,probs = c(0.05,0.1,0.25,0.5,0.75,0.9,0.95))
    meanX <- mean(X) ; sdX <- sd(X) ; maxX <- max(X) ; minX <- min(X)
    ECDF <- ecdf(X)
    P.cdf <- ECDF(cdf)
    names(P.cdf) <- paste("Pr(x<=",cdf,")",sep="")
    # SGT fit
    if(length(X) > 10000) fit.sgt <- FALSE
    if(fit.sgt){
      library(sgt)
      SGT <- sgt.mle(~X,start = list(mu=meanX,sigma=sdX,lambda=0,p=2,q=2))
    } 
    
    if(plot.it){
      par(mar=c(2.5,4.1,3.5,1.5))
      if(is.na(main)){
        m1 <- "Summarize Univariate distributioN"
      }else{
        m1 <- main
      }
      plot(ECDF,main=m1,cex=0.5,xlim=Xlim,ylim=Ylim)
      hist(X,add=TRUE,freq = F,border = "darkgrey",breaks=max(10,length(X)/20))
      rug(X,col="black")
      curve(dnorm(x,mean = meanX,sd=sdX),add=TRUE,col="blue",lwd=2)
      abline(v=meanX,col="blue",lwd=2)
      if(fit.sgt) curve(dsgt(x,mu=SGT$estimate[1],sigma=SGT$estimate[2],lambda=SGT$estimate[3],p=SGT$estimate[4],q=SGT$estimate[5] ),col="green",lty=2,add=TRUE)
      
      # Quantiles
      X0 <- c("25%"=minX,"50%"=NA,"75%"=maxX)
      for(i in c(0.25,0.5,0.75)){
        qy <- paste(i*100,"%",sep = "")
        segments(x0=X0[qy],y0=i,x1=Qs[qy], y1=i,lty=3,col="red",lwd=2)
        segments(x0=Qs[qy],y0=0,x1=Qs[qy], y1=i,lty=3,col="red",lwd=2)
      }
      # Legend
      position <- ifelse(meanX - minX > 0.5*(maxX - meanX), "topleft","right")
      if(fit.sgt){
        legend(position, inset = 0.02,bty="n",cex=0.7,legend = c("ECDF","mean/normal","Q25/50/75","SGT"),lty = c(1,1,3,2),col=c("black","blue","red","green"),pch=c(20,NA,NA,NA))
      }else{
        legend(position, inset = 0.02,bty="n",cex=0.7,legend = c("ECDF","mean/normal","Q25/50/75"),lty = c(1,1,3),col=c("black","blue","red"),pch=c(20,NA,NA))
      } 
    }
    out <- c(mean=meanX,sd=sdX,P.cdf,min=minX,Qs,max=maxX)
    out <- round(out,round.to)
    print(out)
    print(c(n=nobs,NAs=NAs,Infs=Infs))
    out <- c(out,n=nobs,NAs=NAs,Infs=Infs)
    if(fit.sgt) out.sgt <- c("SGT"=round(SGT$estimate,round.to),"SGT.pq"=round(prod(SGT$estimate[4:5]),round.to))
    if(fit.sgt) print(out.sgt)
    if(fit.sgt) out <- c(out,out.sgt)
    invisible(out)
  }
}

round_df <- function(df, digits=3) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  return(df)
}

'
# logit to probability
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
# Color legend
legend.col <- function(col, lev){
  # https://aurelienmadouasse.wordpress.com/2012/01/13/legend-for-a-continuous-color-scale-in-r/
  
  opar <- par
  
  n <- length(col)
  
  bx <- par("usr")
  
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  
  xx <- rep(box.cx, each = 2)
  
  par(xpd = TRUE)
  for(i in 1:n){
    
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
    
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
}
# R2 calculation
r2 <- function(actual,predict){
  R2 <- 1 - (sum((actual-predict )^2)/sum((actual-mean(actual))^2))
  return(R2)
}
'