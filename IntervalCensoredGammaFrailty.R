############################
# Modeling interval-censored, clustered cow udder quarter infection times through the shared gamma frailty model
############################
# Goethals et al (2009)
# Electronic Supplementary Material
# https://link.springer.com/article/10.1198/jabes.2009.0001
############################

fund_ids <- names(table(g.sum.para.cox$Fund_ID)[table(g.sum.para.cox$Fund_ID) < 15]) # filter out FoFs
g.sum.para.cox <- g.sum.para.cox[g.sum.para.cox$Fund_ID %in% fund_ids, ]
g.sum.para.cox <- g.sum.para.cox[g.sum.para.cox$Timing > 0.25, ]

cluster <- g.sum.para.cox$Fund_ID
upper <- g.sum.para.cox$Timing
lower <- upper - 0.25
fail <- ifelse(g.sum.para.cox$RightCensored == 1, 0, 1)
X1 <- pmax(0, g.sum.para.cox$YearInvest - g.sum.para.cox$Vintage)
X2 <- g.sum.para.cox$MSCI.Multiple
X3 <- g.sum.para.cox$GLO

#Calculate the number of clusters.
clusternames <- levels(as.factor(cluster))
ncluster <- length(clusternames)


# Create a data set with the variables cluster, the lower bound,
#the upper bound, the censoring indicator and the covariates.
datasetint <- as.matrix(cbind(cluster,lower,upper,fail,X1,X2,X3))

ncovar = 3


# create subsets for right-censored and interval-censored observations
cendata <- datasetint[datasetint[,4]==0,]
intdata <- datasetint[datasetint[,4]==1,]

# Create a list of signs that corresponds to the n_ik
# (here restricted to 4 events)
signs <- list(1,c(1,-1))
for(i in 3:20) signs[[i]]<-kronecker(signs[[i-1]],c(1,-1))

# Function to calculate the loglikelihood per cluster
CalcLogLikClust <-function(i,x)
{
  theta<-exp(x[1])   
  lambda<-exp(x[2])	
  gamma<-x[3]
  beta<-x[4:(3+ncovar)]
  if(ncovar==1) #univariate case
  {	
    cenX<-cendata[cendata[,1]==clusternames[i],5]
    if (length(cenX)==0)Ci<-0
    else {Ci<-sum(x[2]*as.vector(cendata[cendata[,1]==clusternames[i],3])^x[3]*exp(cenX*beta))}
    intL<-intdata[intdata[,1]==clusternames[i],2]
    
    # if there are no events in that cluster
    nevents <- length(intL)
    crossprod <- 1
    if(nevents>0)
    {
      intX<-intdata[intdata[,1]==clusternames[i],5]
      intRster <- x[2]*(intdata[intdata[,1]==clusternames[i],3]^x[3])*exp(intX*beta)
      intLster <- x[2]*(intL^x[3])*exp(intX*beta)
      crossprod<-c(exp(intLster[1]),exp(intRster[1]))
      if(nevents>1)
      {
        for(ik in 2:nevents)
        {
          crossprod <- kronecker(crossprod, c(exp(intLster[ik]),exp(intRster[ik]))); 
        }
      }
    }
  }
  else #multivariate
  {
    #cenR<-cendata[cendata[,1]==clusternames[i],3]
    cenX<-cendata[cendata[,1]==clusternames[i],5:(4+ncovar)]
    if (length(cendata[cendata[,1]==clusternames[i],5])==0)Ci<-0
    else 
    {
      Ci<-sum(x[2]*(as.vector(cendata[cendata[,1]==clusternames[i],3])^x[3])*exp(cenX%*%beta))
    }
    
    # if there are no events in that cluster
    crossprod <- 1
    intL<-intdata[intdata[,1]==clusternames[i],2]
    nevents <- length(intL)
    if(nevents>0)
    {
      intX<-intdata[intdata[,1]==clusternames[i],5:(4+ncovar)]
      expiXb <- exp(intX%*%beta)		
      intRster <- x[2]*(intdata[intdata[,1]==clusternames[i],3]^x[3])*expiXb
      intLster <- x[2]*(intL^x[3])*expiXb
      crossprod<-c(exp(intLster[1]),exp(intRster[1]))
      if(nevents>1)
      {
        for(ik in 2:nevents)
        {
          crossprod <- kronecker(crossprod, c(exp(intLster[ik]),exp(intRster[ik]))); 
        }
      }
    }
    
  }	
  
  # Loglikelihood for 1 cluster
  log(1/(theta^(1/theta))*sum((1/((sum(lambda*
                                         (as.vector(cendata[cendata[,1]==clusternames[i],3])^gamma)
                                       *exp(cenX*beta))+1/theta+log(crossprod))^(1/theta)))*signs[[nevents+1]]))
}

# Calculate full marginal loglikelihood (formula 5)
CalcLogLik <- function(x)
{
  -sum(sapply(1:ncluster,CalcLogLikClust,x=x))
}

# Maximising the full marginal loglikelihood to obtain parameter estimates
init<-c(1,1,1,1)
print(results <- nlm(CalcLogLik,init,print.level=2, hessian=TRUE))

# Calculate covariance matrix
covmatr<-solve(results$hessian)
