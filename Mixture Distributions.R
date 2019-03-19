## mixdist & flexmix -------
library(mixdist)
data(pike65)
data(pikepar)
fitpike1 <- mix(pike65, pikepar, "lnorm", constr = mixconstr(consigma = "CCV"), emsteps = 3)
summary(fitpike1)
plot(fitpike1)

library(flexmix)

## JM::weibull.frailty --------
# https://rdrr.io/github/drizopoulos/JM/src/R/weibull.frailty.R

weibull.frailty <-
  function (formula = formula(data), data = parent.frame(), id = "id", subset, na.action, init, 
            control = list()) {
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "subset", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
    Terms <- if (missing(data)) terms(formula) else terms(formula, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    if (NROW(m) == 0) 
      stop("No (non-missing) observations.\n")
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) 
      stop("Response must be a survival object.\n")
    logT <- log(Y[, 1])
    d <- Y[, 2]
    id <- if (is.character(id) && length(id) == 1) {
      if (missing(data) || !id %in% names(data))
        stop("'id' not a 'data'.\n")
      nam.id <- id
      dd <- if (missing(na.action)) na.omit(data) else na.action(data)
      as.vector(unclass(factor(dd[[id]])))
    } else {
      as.vector(unclass(factor(id)))
    }
    attr(Terms, "intercept") <- 1
    X <- model.matrix(Terms, m)[, -1, drop = FALSE]
    type <- attr(Y, "type")
    if (type != "right") 
      stop("weibull.frailty() currently supports only right-censored data.\n")
    if (missing(init)) 
      init <- NULL
    out <- weibull.frailty.fit(logT, d, X, id, init, control)
    out$y <- Y
    out$x <- X
    out$id <- id
    out$nam.id <- nam.id
    out$terms <- Terms
    out$data <- if (missing(data)) m else data
    out$call <- call
    class(out) <- "weibull.frailty"
    out
  }

weibull.frailty.fit <-
  function (logT, d, X, id, init.thetas, control = list()) {
    Time <- exp(logT)
    p <- ncol(X)
    N <- nrow(X)
    n <- length(unique(id))
    D <- as.vector(tapply(d, id, sum))
    if (!ncol(X))
      X <- as.matrix(rep(0, N))
    Xd <- colSums(X * d)
    sum.d <- sum(d)
    sum.dlogT <- sum(d * logT)
    con <- list(optimizer = "optim", parscale = NULL, maxit = 500, numeriDeriv = "cd", eps.Hes = 1e-03)
    con[names(control)] <- control
    clnams <- colnames(X)
    dimnames(X) <- names(logT) <- names(Time) <- names(d) <- names(id) <- names(Xd) <- NULL 
    fn <- function (thetas) {
      betas <- thetas[seq_len(p)]
      scale <- exp(thetas[p + 1])
      shape <- exp(thetas[p + 2])
      var.fr <- exp(thetas[p + 3])
      theta <- 1 / var.fr
      eta <- if (p > 0) c(X %*% betas) else rep(0, N)
      log.lambda0 <- log(shape * scale) + (shape - 1) * logT
      Lambda0 <- scale * Time^shape
      P1 <- sum(d * (log.lambda0 + eta))
      P2 <- n * (theta * log(theta) - lgamma(theta))
      P3 <- sum(lgamma(D + theta) - (D + theta) * log(theta + c(tapply(Lambda0 * exp(eta), id, sum))))
      - (P1 + P2 + P3)
    }
    gr <- function (thetas) {
      betas <- thetas[seq_len(p)]
      scale <- exp(thetas[p + 1])
      shape <- exp(thetas[p + 2])
      var.fr <- exp(thetas[p + 3])
      theta <- 1 / var.fr
      eta <- if (p > 0) c(X %*% betas) else rep(0, N)
      exp.eta <- exp(eta)
      log.lambda0 <- log(shape * scale) + (shape - 1) * logT
      Lambda0 <- scale * Time^shape
      Lambda0.eta <- Lambda0 * exp(eta)
      mat.id <- cbind(Lambda0.eta, Time^shape * exp(eta), logT * Lambda0.eta, Lambda0.eta * X)
      P <- rowsum(mat.id, id, FALSE)
      Lambda0.eta <- P[, 1]
      theta.Lambda0.eta <- theta + Lambda0.eta
      log.theta.Lambda0.eta <- log(theta.Lambda0.eta)
      X.Lambda0.eta <- P[, seq(4, ncol(P)), drop = FALSE]
      sc.betas <- - c(Xd - colSums((D + theta) * X.Lambda0.eta / theta.Lambda0.eta))
      sc.scale <- - scale * c(sum.d / scale - sum((D + theta) * P[, 2] / theta.Lambda0.eta))
      sc.shape <- - shape * (sum.d / shape + sum.dlogT - sum((D + theta) * P[, 3] / theta.Lambda0.eta))
      sc.var.fr <- theta * sum(log(theta) + 1 - digamma(theta) + digamma(D + theta) - 
                                 log.theta.Lambda0.eta - (D + theta) / theta.Lambda0.eta)
      if (p > 0) c(sc.betas, sc.scale, sc.shape, sc.var.fr) else c(sc.scale, sc.shape, sc.var.fr)
    }
    if (is.null(init.thetas) || length(init.thetas) != p + 3)
      init.thetas <- rep(0.01, p + 3)
    names(init.thetas) <- NULL
    psc <- if (is.null(xx <- con$parscale)) rep(c(1, 0.1), c(p, 3)) else xx
    opt <- if (con$optimizer == "optim") {
      optim(init.thetas, fn, gr, method = "BFGS", control = list(maxit = con$maxit, parscale = psc))
    } else {
      nlminb(init.thetas, fn, gr, scale = psc, control = list(iter.max = con$maxit))
    }
    H <- if (con$numeriDeriv == "cd") cd.vec(opt$par, gr, eps = con$eps.Hes) else fd.vec(opt$par, gr, eps = con$eps.Hes)
    if (any(is.na(H) | !is.finite(H))) {
      warning("infinite or missing values in Hessian at convergence.\n")
    } else {
      ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
      if (!all(ev >= -1e-06 * abs(ev[1]))) 
        warning("Hessian matrix at convergence is not positive definite.\n")
    }
    betas <- opt$par[seq_len(p)]
    names(betas) <- clnams
    scale <- exp(opt$par[p + 1])
    shape <- exp(opt$par[p + 2])
    var.fr <- exp(opt$par[p + 3])
    list(coefficients = list(betas = betas, scale = scale, shape = shape, var.frailty = var.fr), hessian = H, 
         logLik = -opt[[2]], control = con)
  }

## parfm ---------
library(parfm)

parfm <- 
function (formula, cluster = NULL, strata = NULL, data, inip = NULL, 
          iniFpar = NULL, dist = c("weibull", "inweibull", "frechet", 
                                   "exponential", "gompertz", "loglogistic", "lognormal", 
                                   "logskewnormal"), frailty = c("none", "gamma", "ingau", 
                                                                 "possta", "lognormal", "loglogistic"), method = "nlminb", 
          maxit = 500, Fparscale = 1, showtime = FALSE, correct = 0) 
{
  if (missing(data)) {
    data <- eval(parse(text = paste("data.frame(", paste(all.vars(formula), 
                                                         collapse = ", "), ")")))
  }
  dist <- tolower(dist)
  dist <- match.arg(dist)
  frailty <- tolower(frailty)
  frailty <- match.arg(frailty)
  if (frailty == "none" && !is.null(cluster)) {
    warning(paste0("With frailty='none' the cluster variable '", 
                   cluster, "' is not used!"))
  }
  if (frailty == "none" && !is.null(iniFpar)) {
    warning("With frailty='none' the argument 'iniFpar' is not used!")
  }
  if (frailty == "possta") {
    if (10^correct == Inf || 10^-correct == 0) {
      stop("'correct' is too large!")
    }
    if (10^correct == 0 || 10^-correct == Inf) {
      stop("'correct' is too small!")
    }
  }
  else if (correct != 0) {
    warning(paste0("'correct' has no effect when 'frailty = ", 
                   frailty, "'"))
  }
  obsdata <- NULL
  if (length(formula[[2]]) == 3) {
    obsdata$time <- eval(formula[[2]][[2]], envir = data)
    obsdata$event <- eval(formula[[2]][[3]], envir = data)
  }
  else if (length(formula[[2]]) == 4) {
    obsdata$trunc <- eval(formula[[2]][[2]], envir = data)
    obsdata$time <- eval(formula[[2]][[3]], envir = data)
    obsdata$event <- eval(formula[[2]][[4]], envir = data)
  }
  if (!all(levels(as.factor(obsdata$event)) %in% 0:1)) {
    stop(paste("The status indicator 'event' in the Surv object", 
               "in the left-hand side of the formula object", "must be either 0 (no event) or 1 (event)."))
  }
  obsdata$x <- as.data.frame(model.matrix(formula, data = data))
  if (is.null(cluster)) {
    if (frailty != "none") {
      stop(paste("if you specify a frailty distribution,\n", 
                 "then you have to specify the cluster variable as well"))
    }
    else {
      obsdata$cluster <- rep(1, nrow(data))
    }
    obsdata$ncl <- 1
    obsdata$di <- sum(obsdata$event)
  }
  else {
    if (!cluster %in% names(data)) {
      stop(paste0("object '", cluster, "' not found"))
    }
    obsdata$cluster <- eval(as.name(cluster), envir = data)
    obsdata$ncl <- length(levels(as.factor(obsdata$cluster)))
    obsdata$di <- aggregate(obsdata$event, by = list(obsdata$cluster), 
                            FUN = sum)[, , drop = FALSE]
    cnames <- obsdata$di[, 1]
    obsdata$di <- as.vector(obsdata$di[, 2])
    names(obsdata$di) <- cnames
  }
  if (is.null(strata)) {
    obsdata$strata <- rep(1, length(obsdata$time))
    obsdata$nstr <- 1
    obsdata$dq <- sum(obsdata$event)
  }
  else {
    if (!strata %in% names(data)) {
      stop(paste0("object '", strata, "' not found"))
    }
    obsdata$strata <- eval(as.name(strata), envir = data)
    obsdata$nstr <- length(levels(as.factor(obsdata$strata)))
    obsdata$dq <- aggregate(obsdata$event, by = list(obsdata$strata), 
                            FUN = sum)[, , drop = FALSE]
    snames <- obsdata$dq[, 1]
    obsdata$dq <- as.vector(obsdata$dq[, 2])
    names(obsdata$dq) <- snames
  }
  if (!is.null(cluster) && !is.null(strata)) {
    obsdata$dqi <- xtabs(x ~ Group.1 + Group.2, data = aggregate(obsdata$event, 
                                                                 by = list(obsdata$cluster, obsdata$strata), FUN = sum))
    dimnames(obsdata$dqi) <- list(cluster = dimnames(obsdata$dqi)[[1]], 
                                  strata = dimnames(obsdata$dqi)[[2]])
  }
  else if (!is.null(cluster)) {
    obsdata$dqi <- obsdata$di
  }
  else if (!is.null(strata)) {
    obsdata$dqi <- obsdata$dq
  }
  else {
    obsdata$dqi <- sum(obsdata$event)
  }
  if (frailty == "none") {
    nFpar <- 0
  }
  else if (frailty %in% c("gamma", "ingau", "possta", "lognormal")) {
    nFpar <- 1
  }
  obsdata$nFpar <- nFpar
  if (dist == "exponential") {
    nBpar <- 1
  }
  else if (dist %in% c("weibull", "inweibull", "frechet", "gompertz", 
                       "lognormal", "loglogistic")) {
    nBpar <- 2
  }
  else if (dist %in% c("logskewnormal")) {
    nBpar <- 3
  }
  obsdata$nBpar <- nBpar
  nRpar <- ncol(obsdata$x) - 1
  obsdata$nRpar <- nRpar
  if (!is.null(inip)) {
    if (length(inip) != nBpar * obsdata$nstr + nRpar) {
      stop(paste("number of initial parameters 'inip' must be", 
                 nBpar * obsdata$nstr + nRpar))
    }
    p.init <- inip
    if (dist %in% c("exponential", "weibull", "inweibull", 
                    "frechet", "gompertz")) {
      if (any(p.init[1:obsdata$nstr] <= 0)) {
        stop(paste("with that baseline, the 1st parameter has to be > 0"))
      }
      p.init[1:obsdata$nstr] <- log(p.init[1:obsdata$nstr])
    }
    if (dist %in% c("weibull", "inweibull", "frechet", "gompertz", 
                    "lognormal", "loglogistic", "logskewnormal")) {
      if (any(p.init[obsdata$nstr + 1:obsdata$nstr] <= 
              0)) {
        stop(paste("with that baseline, the 2nd parameter has to be > 0"))
      }
      p.init[obsdata$nstr + 1:obsdata$nstr] <- log(p.init[obsdata$nstr + 
                                                            1:obsdata$nstr])
    }
  }
  else {
    sink("NUL")
    inires <- optimx(par = rep(0, nRpar + nBpar), fn = optMloglikelihood, 
                     method = method, obs = obsdata, dist = dist, frailty = "none", 
                     correct = correct, hessian = FALSE, control = list(maxit = maxit, 
                                                                        starttests = FALSE, dowarn = FALSE))
    sink()
    p.init <- inires[1:(nRpar + nBpar)]
    rm(inires)
  }
  if (frailty == "none") {
    pars <- NULL
  }
  else if (frailty %in% c("gamma", "ingau", "lognormal")) {
    if (is.null(iniFpar)) {
      iniFpar <- 1
    }
    else if (iniFpar <= 0) {
      stop("initial heterogeneity parameter (theta) has to be > 0")
    }
    pars <- log(iniFpar)
  }
  else if (frailty == "possta") {
    if (is.null(iniFpar)) {
      iniFpar <- 0.5
    }
    else if (iniFpar <= 0 || iniFpar >= 1) {
      stop("initial heterogeneity parameter (nu) must lie in (0, 1)")
    }
    pars <- log(-log(iniFpar))
  }
  pars <- c(pars, unlist(p.init))
  res <- NULL
  sink("NUL")
  todo <- expression({
    res <- optimx(par = pars, fn = optMloglikelihood, method = method, 
                  obs = obsdata, dist = dist, frailty = frailty, correct = correct, 
                  hessian = FALSE, control = list(maxit = maxit, starttests = FALSE, 
                                                  dowarn = FALSE))
  })
  if (showtime) {
    extime <- system.time(eval(todo))[1]
  }
  else {
    eval(todo)
    extime <- NULL
  }
  sink()
  if (res$convcode > 0) {
    warning("optimisation procedure did not converge,\n              conv = ", 
            bquote(.(res$convergence)), ": see ?optimx for details")
  }
  it <- res$niter
  lL <- -res$value
  if (frailty == "possta") {
    lL <- lL + correct * log(10) * obsdata$ncl
  }
  estim_par <- as.numeric(res[1:(nFpar + nBpar * obsdata$nstr + 
                                   nRpar)])
  if (frailty %in% c("gamma", "ingau")) {
    theta <- exp(estim_par[1:nFpar])
    sigma2 <- NULL
    nu <- NULL
  }
  else if (frailty == "lognormal") {
    theta <- NULL
    sigma2 <- exp(estim_par[1:nFpar])
    nu <- NULL
  }
  else if (frailty == "possta") {
    theta <- NULL
    sigma2 <- NULL
    nu <- exp(-exp(estim_par[1:nFpar]))
  }
  else if (frailty == "none") {
    theta <- NULL
    sigma2 <- NULL
    nu <- NULL
  }
  if (dist == "exponential") {
    lambda <- exp(estim_par[nFpar + 1:obsdata$nstr])
    ESTIMATE <- c(lambda = lambda)
  }
  else if (dist %in% c("weibull", "inweibull", "frechet")) {
    rho <- exp(estim_par[nFpar + 1:obsdata$nstr])
    lambda <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(rho = rho, lambda = lambda)
  }
  else if (dist == "gompertz") {
    gamma <- exp(estim_par[nFpar + 1:obsdata$nstr])
    lambda <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(gamma = gamma, lambda = lambda)
  }
  else if (dist == "lognormal") {
    mu <- estim_par[nFpar + 1:obsdata$nstr]
    sigma <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(mu = mu, sigma = sigma)
  }
  else if (dist == "loglogistic") {
    alpha <- estim_par[nFpar + 1:obsdata$nstr]
    kappa <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(alpha = alpha, kappa = kappa)
  }
  else if (dist == "logskewnormal") {
    xi <- estim_par[nFpar + 1:obsdata$nstr]
    omega <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    alpha <- estim_par[nFpar + 2 * obsdata$nstr + 1:obsdata$nstr]
    ESTIMATE <- c(xi = xi, omega = omega, alpha = alpha)
  }
  if (nRpar == 0) {
    beta <- NULL
  }
  else {
    beta <- estim_par[-(1:(nFpar + nBpar * obsdata$nstr))]
    names(beta) <- paste("beta", names(obsdata$x), sep = ".")[-1]
  }
  ESTIMATE <- c(theta = theta, sigma2 = sigma2, nu = nu, ESTIMATE, 
                beta = beta)
  resHessian <- optimHess(par = ESTIMATE, fn = Mloglikelihood, 
                          obs = obsdata, dist = dist, frailty = frailty, correct = correct, 
                          transform = FALSE)
  var <- try(diag(solve(resHessian)), silent = TRUE)
  if (class(var) == "try-error" | any(is.nan(var))) {
    warning(var[1])
    STDERR <- rep(NA, nFpar + nBpar * obsdata$nstr + nRpar)
    PVAL <- rep(NA, nFpar + nBpar * obsdata$nstr + nRpar)
  }
  else {
    if (any(var <= 0)) {
      warning(paste("negative variances have been replaced by NAs\n", 
                    "Please, try other initial values", "or another optimisation method"))
    }
    if (frailty %in% c("gamma", "ingau")) {
      seTheta <- sapply(1:nFpar, function(x) {
        ifelse(var[x] > 0, sqrt(var[x]), NA)
      })
      seSigma2 <- seNu <- NULL
    }
    else if (frailty == "lognormal") {
      seSigma2 <- sapply(1:nFpar, function(x) {
        ifelse(var[x] > 0, sqrt(var[x]), NA)
      })
      seTheta <- seNu <- NULL
    }
    else if (frailty == "possta") {
      seNu <- sapply(1:nFpar, function(x) {
        ifelse(var[x] > 0, sqrt(var[x]), NA)
      })
      seTheta <- seSigma2 <- NULL
    }
    if (dist == "exponential") {
      seLambda <- sapply(1:obsdata$nstr, function(x) {
        ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), 
               NA)
      })
      STDERR <- c(seLambda = seLambda)
    }
    else if (dist %in% c("weibull", "inweibull", "frechet")) {
      seRho <- sapply(1:obsdata$nstr, function(x) {
        ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), 
               NA)
      })
      seLambda <- sapply(1:obsdata$nstr, function(x) {
        ifelse(var[nFpar + obsdata$nstr + x] > 0, sqrt(var[nFpar + 
                                                             obsdata$nstr + x]), NA)
      })
      STDERR <- c(seRho = seRho, seLambda = seLambda)
    }
    else if (dist == "gompertz") {
      seGamma <- sapply(1:obsdata$nstr, function(x) {
        ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), 
               NA)
      })
      seLambda <- sapply(1:obsdata$nstr, function(x) {
        ifelse(var[nFpar + obsdata$nstr + x] > 0, sqrt(var[nFpar + 
                                                             obsdata$nstr + x]), NA)
      })
      STDERR <- c(seGamma = seGamma, seLambda = seLambda)
    }
    else if (dist == "lognormal") {
      seMu <- sapply(1:obsdata$nstr, function(x) {
        ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), 
               NA)
      })
      seSigma <- sapply(1:obsdata$nstr, function(x) {
        ifelse(var[nFpar + obsdata$nstr + x] > 0, sqrt(var[nFpar + 
                                                             obsdata$nstr + x]), NA)
      })
      STDERR <- c(seMu = seMu, seSigma = seSigma)
    }
    else if (dist == "loglogistic") {
      seAlpha <- sapply(1:obsdata$nstr, function(x) {
        ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), 
               NA)
      })
      seKappa <- sapply(1:obsdata$nstr, function(x) {
        ifelse(var[nFpar + obsdata$nstr + x] > 0, sqrt(var[nFpar + 
                                                             obsdata$nstr + x]), NA)
      })
      STDERR <- c(seAlpha = seAlpha, seKappa = seKappa)
    }
    else if (dist == "logskewnormal") {
      seXi <- sapply(1:obsdata$nstr, function(x) {
        ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), 
               NA)
      })
      seOmega <- sapply(1:obsdata$nstr, function(x) {
        ifelse(var[nFpar + obsdata$nstr + x] > 0, sqrt(var[nFpar + 
                                                             obsdata$nstr + x]), NA)
      })
      seAlpha <- sapply(1:obsdata$nstr, function(x) {
        ifelse(var[nFpar + 2 * obsdata$nstr + x] > 0, 
               sqrt(var[nFpar + 2 * obsdata$nstr + x]), NA)
      })
      STDERR <- c(seXi = seXi, seOmega = seOmega, seAlpha = seAlpha)
    }
    if (nRpar == 0) {
      seBeta <- NULL
    }
    else {
      seBeta <- numeric(nRpar)
      varBeta <- var[-(1:(nFpar + nBpar * obsdata$nstr))]
      for (i in 1:nRpar) {
        seBeta[i] <- ifelse(varBeta[i] > 0, sqrt(varBeta[i]), 
                            NA)
      }
      PVAL <- c(rep(NA, nFpar + nBpar * obsdata$nstr), 
                2 * pnorm(q = -abs(beta/seBeta)))
    }
    STDERR <- c(STDERR, se.beta = seBeta)
    if (frailty != "none") {
      STDERR <- c(se.theta = seTheta, se.sigma2 = seSigma2, 
                  se.nu = seNu, STDERR)
    }
  }
  resmodel <- cbind(ESTIMATE = ESTIMATE, SE = STDERR)
  rownames(resmodel) <- gsub("beta.", "", rownames(resmodel))
  if (nRpar > 0) 
    resmodel <- cbind(resmodel, `p-val` = PVAL)
  class(resmodel) <- c("parfm", class(resmodel))
  Call <- match.call()
  if (!match("formula", names(Call), nomatch = 0)) 
    stop("A formula argument is required")
  Terms <- terms(formula, data = data)
  attributes(resmodel) <- c(attributes(resmodel), list(call = Call, 
                                                       convergence = res$convergence, it = it, extime = extime, 
                                                       nobs = nrow(data), shared = (nrow(data) > obsdata$ncl), 
                                                       loglik = lL, dist = dist, cumhaz = attributes(Mloglikelihood(p = estim_par, 
                                                                                                                    obs = obsdata, dist = dist, frailty = frailty, correct = correct))$cumhaz, 
                                                       cumhazT = attributes(Mloglikelihood(p = estim_par, obs = obsdata, 
                                                                                           dist = dist, frailty = frailty, correct = correct))$cumhazT, 
                                                       di = obsdata$di, dq = obsdata$dq, dqi = obsdata$dqi, 
                                                       frailty = frailty, clustname = cluster, stratname = strata, 
                                                       correct = correct, formula = as.character(Call[match("formula", 
                                                                                                            names(Call), nomatch = 0)]), terms = attr(Terms, 
                                                                                                                                                      "term.labels"), FisherI = resHessian))
  if (frailty != "none") {
    names(attr(resmodel, "cumhaz")) <- names(attr(resmodel, 
                                                  "di")) <- unique(obsdata$cluster)
  }
  if (showtime) 
    cat("\nExecution time:", extime, "second(s) \n")
  return(resmodel)
}