####################
## Generalized Gamma
####################
library(gamlss)

# Create uniforms for truncated distributions
U <- 0
truncation_threshold <- 0.1
prob_of_default <- pGG(truncation_threshold)
while(U < prob_of_default){
  U <- runif(1)
}
U

# Create quantile function for Generalized Gamma distribution
# qGG
qGA
qNO
pGG

quantile.f.Gamma <- function (p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\\n", ""))
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\\n", ""))
  q <- qgamma(p, shape = 1/sigma^2, scale = mu * sigma^2, lower.tail = lower.tail, log.p = log.p)
  q
}
quantile.f.Normal <- function (p, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\\n", ""))
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\\n", ""))
  q <- qnorm(p, mean = mu, sd = sigma, lower.tail = lower.tail)
  q
}
quantile.f.GG <- function (p, mu = 1, sigma = 0.5, nu = 2, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu < 0)) 
    stop(paste("mu must be positive", "\\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\\n", ""))
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\\n", ""))
  if (length(nu) > 1) {
    p <- ifelse(nu > 0, p, 1 - p)
    z <- ifelse(abs(nu) > 1e-06, 
                quantile.f.Gamma(p, mu = 1, sigma = sigma * abs(nu)), 
                quantile.f.Normal(p, mu = log(mu), sigma = sigma))
    y <- ifelse(abs(nu) > 1e-06, 
                mu * z^(1/nu), 
                exp(z))
  }
  else if (abs(nu) > 1e-06) {
    p <- if (nu > 0) 
      p
    else 1 - p
    z <- quantile.f.Gamma(p, mu = 1, sigma = sigma * abs(nu))
    y <- mu * z^(1/nu)
  }
  else {
    z <- quantile.f.Normal(p, mu = log(mu), sigma = sigma)
    y <- exp(z)
  }
  y
}
CDF_GG <- function (q, mu = 1, sigma = 0.5, nu = 2, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\\n", ""))
  if (any(q < 0)) 
    stop(paste("q must be positive", "\\n", ""))
  z <- (q/mu)^nu
  if (length(nu) > 1) 
    cdf <- ifelse(abs(nu) > 1e-06, 
                  pGA(z, mu = 1, sigma = sigma * abs(nu), lower.tail = (nu < 0) - lower.tail, log.p = log.p), 
                  pNO(log(q), mu = log(mu), sigma = sigma))
                  # pNO(log(z), mu = log(mu), sigma = sigma))
  else if (abs(nu) > 1e-06) 
    cdf <- pGA(z, mu = 1, sigma = sigma * abs(nu), lower.tail = (nu < 0) - lower.tail, log.p = log.p)
  else cdf <- pNO(log(q), mu = log(mu), sigma = sigma)
  cdf
}
CDF_Gamma <- function (q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\\n", ""))
  if (any(q < 0)) 
    stop(paste("y must be positive", "\\n", ""))
  cdf <- pgamma(q, shape = 1/sigma^2, scale = mu * sigma^2, 
                lower.tail = lower.tail, log.p = log.p)
  cdf
}
CDF_Normal <- function (q, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\\n", ""))
  cdf <- pnorm(q, mean = mu, sd = sigma, lower.tail = lower.tail, log.p = log.p)
  cdf
}

CDF_GG(0.5)

nu <- 0
q <- 0.1
CDF_GG(q, 1, 0.5, nu)
CDF_GG(c(q, q), c(1,1), c(0.5,0.5), c(nu, nu))


quantile.f.GG(0.25)


#############################
## Binomial (link == "logit")
#############################

BI
pBI


##############
## GG flexsurv
##############

flexsurv::qgengamma

flexsurv_qgengamma <- function (p, mu = 0, sigma = 1, Q, lower.tail = TRUE, log.p = FALSE) {
  d <- dbase("gengamma", lower.tail = lower.tail, log = log.p, 
             p = p, mu = mu, sigma = sigma, Q = Q)
  for (i in seq_along(d)) assign(names(d)[i], d[[i]])
  p[Q < 0] <- 1 - p[Q < 0]
  ret[ind] <- numeric(sum(ind))
  ret[ind][Q == 0] <- qlnorm(p[Q == 0], mu[Q == 0], 1/sigma[Q ==  0]^2)
  qn0 <- Q != 0
  p <- p[qn0]
  mu <- mu[qn0]
  sigma <- sigma[qn0]
  Q <- Q[qn0]
  ret[ind][qn0] <- exp(mu + sigma * (log(Q^2 * qgamma(p, 1/Q^2, 1))/Q))
  ret
}


