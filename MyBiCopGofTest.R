library(mvtnorm)
library(VineCopula)
#-------------------------------------
MyBiCopGofTest<-function (u1, u2, family, par = 0, par2 = 0, method = "white", 
          max.df = 30, B = 100, obj = NULL) 
{
  if (method == "White") 
    method <- "white"
  if (method == "Kendall") 
    method <- "kendall"
  args <- preproc(c(as.list(environment()), call = match.call()), 
                  check_u, remove_nas, check_nobs, check_if_01, extract_from_BiCop, 
                  na.txt = " Only complete observations are used.")
  list2env(args, environment())
  allfams <- c(0:10,
               13, 14, 16:20,
               23, 24, 26:30, 33, 34, 36:40,
               104, 114, 124, 134, 204, 214, 224, 234)
  tawns <- which(allfams > 100)
  if (!(family %in% allfams[-tawns])) 
    stop("Copula family not implemented.")
  if (par != 0) 
    BiCopCheck(family, par, par2)
  if (family == 2 && method == "kendall") 
    stop("The goodness-of-fit test based on Kendall's process is not ", 
         "\n             implemented for the t-copula.")
  if (family %in% c(7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 
                    30, 37, 38, 39, 40) && method == "white") 
    stop("The goodness-of-fit test based on White's information matrix ", 
         "equality is not implemented for the BB copulas.")
  T <- length(u1)
  if (method == "white") {
    if (family == 2) {
      if (par == 0) {
        pars <- BiCopEst(u1, u2, family = family, method = "mle", 
                         max.df = max.df)
        theta <- pars$par
        nu <- pars$par2
      }
      else {
        theta <- par
        nu <- par2
      }
    }
    else {
      nu <- 0
      theta <- BiCopEst(u1, u2, family = family, method = "mle")$par
    }
    if (family == 2) {
      Dprime <- matrix(0, 3, T)
      Vt <- array(0, dim = c(3, 3, T))
      Bt <- array(0, dim = c(2, 2, T))
      grad <- c(0, 0)
      gradD <- gradDtcopula(u1, u2, theta, nu)
      for (t in 1:T) {
        H <- hesseTcopula(u1[t], u2[t], theta, nu)
        Hprime <- as.vector(H[lower.tri(H, diag = TRUE)])
        C <- OPGtcopula(u1[t], u2[t], family, theta, 
                        nu)
        Cprime <- as.vector(C[lower.tri(C, diag = TRUE)])
        Dprime[, t] <- Hprime + Cprime
        Bt[, , t] <- H
        tmp <- Dprime[, t] - gradD %*% solve(Bt[, , t]) %*% 
          grad
        Vt[, , t] <- (tmp) %*% t(tmp)
      }
      D <- apply(Dprime, 1, mean)
      V0 <- apply(Vt, c(1, 2), mean)
    }
    else {
      b <- BiCopPDF(u1, u2, family, theta, nu)
      d <- BiCopDeriv2(u1, u2, family, theta, nu, deriv = "par")/b
      D <- mean(d)
      eps <- 1e-04
      b_eps1 <- BiCopPDF(u1, u2, family, theta - eps, nu, 
                         check.pars = FALSE)
      d_eps1 <- BiCopDeriv(u1, u2, family, theta - eps, 
                           nu, deriv = "par", check.pars = FALSE)/b_eps1
      gradD_1 <- mean(d_eps1)
      b_eps2 <- BiCopPDF(u1, u2, family, theta + eps, nu, 
                         check.pars = FALSE)
      d_eps2 <- BiCopDeriv(u1, u2, family, theta + eps, 
                           nu, deriv = "par", check.pars = FALSE)/b_eps2
      gradD_2 <- mean(d_eps2)
      gradD <- (gradD_2 - gradD_1)/(2 * eps)
      tmp1 <- BiCopDeriv(u1, u2, family, theta, nu, deriv = "par")
      tmp2 <- tmp1/b^2
      tmp3 <- -tmp2 + d
      H <- mean(tmp3)
      Vt <- (d - gradD/H * tmp1/b)^2
      V0 <- mean(Vt)
    }
    if (family == 2) {
      handle <- try(solve(V0), TRUE)
      if (is.null(dim(handle))) 
        handle <- ginv(V0)
      test <- T * (t(D) %*% handle %*% D)
      pvalue <- 1 - pchisq(test, df = length(D))
    }
    else {
      test <- T * D * solve(V0) * D
      pvalue <- 1 - pchisq(test, df = 1)
    }
    if (B == 0) {
      out <- list(statistic = test, p.value = pvalue)
    }
    else {
      test_boot <- bootWhite(family, theta, nu, B, N = length(u1))
      test_boot <- test_boot[!is.na(test_boot)]
      pvalue <- mean(test_boot >= as.numeric(test))
      out <- list(statistic = test, p.value = pvalue)
    }
  }
  else if (method == "IR") {
    if (family == 2) {
      if (par == 0) {
        pars <- BiCopEst(u1, u2, family = family, method = "mle", 
                         max.df = max.df)
        theta <- pars$par
        nu <- pars$par2
      }
      else {
        theta <- par
        nu <- par2
      }
    }
    else {
      nu <- 0
      theta <- BiCopEst(u1, u2, family = family, method = "mle")$par
    }
    if (family == 2) {
      grad <- c(0, 0)
      rho_teil <- f_rho(u1, u2, theta, nu)
      nu_teil <- f_nu(u1, u2, theta, nu)
      rho_nu_teil <- f_rho_nu(u1, u2, theta, nu)
      H <- matrix(c(rho_teil, rho_nu_teil, rho_nu_teil, 
                    nu_teil), 2, 2)
      grad[1] <- BiCopDeriv(u1, u2, family = family, par = theta, 
                            par2 = nu, deriv = "par", log = TRUE, check.pars = FALSE)
      grad[2] <- BiCopDeriv(u1, u2, family = family, par = theta, 
                            par2 = nu, deriv = "par2", log = TRUE, check.pars = FALSE)
      C <- grad %*% t(grad)
    }
    else {
      d <- rep(0, T)
      for (t in 1:T) {
        b <- BiCopPDF(u1[t], u2[t], family, theta, nu, 
                      check.pars = FALSE)
        d[t] <- BiCopDeriv2(u1[t], u2[t], family = family, 
                            par = theta, par2 = nu, deriv = "par", check.pars = FALSE)/b
      }
      H <- mean(d)
      C <- BiCopDeriv(u1, u2, family = family, par = theta, 
                      par2 = nu, deriv = "par", log = TRUE, check.pars = FALSE)
    }
    Phi <- -solve(H) %*% C
    IR <- trace(Phi)/dim(H)[1]
    if (B == 0) {
      out <- list(IR = IR, p.value = NULL)
    }
    else {
      IR_boot <- boot.IR(family, theta, nu, B, length(u1))
      sigma2 <- var(IR_boot)
      IR_new <- ((IR - 1)/sqrt(sigma2))^2
      IR_boot <- ((IR_boot - 1)/sqrt(sigma2))^2
      p.value <- mean(IR_boot >= IR_new)
      out <- list(IR = IR, p.value = p.value)
    }
  }
  else if (method == "kendall") {
    if (family %in% c(13, 14, 16, 17, 18, 19, 20)) {
      u1 <- 1 - u1
      u2 <- 1 - u2
      family <- family - 10
    }
    else if (family %in% c(23, 24, 26, 27, 28, 29, 30)) {
      u1 <- 1 - u1
      family <- family - 20
    }
    else if (family %in% c(33, 34, 36, 37, 38, 39, 40)) {
      u2 <- 1 - u2
      family <- family - 30
    }
    param <- suppressWarnings({
      BiCopEst(u1, u2, family = family)
    })
    ostat <- obs.stat(u1, u2, family, param)
    if (B == 0) {
      sn.obs <- ostat$Sn
      tn.obs <- ostat$Tn
      out <- list(Sn = sn.obs, Tn = tn.obs)
    }
    else {
      numError <- 0
      sn.boot <- rep(0, B)
      tn.boot <- rep(0, B)
      for (i in 1:B) {
        ax <- try({
          tmp <- boot.stat(u1, u2, family, param)
        }, silent = TRUE)
        if (inherits(ax, "try-error")) {
          sn.boot[i] <- NA
          tn.boot[i] <- NA
          numError <- numError + 1
        }
        else {
          sn.boot[i] <- tmp$sn
          tn.boot[i] <- tmp$tn
        }
      }
      #if (numError > 0) {
      #  warning(paste("In the calculation of the p-values for the copula\ngoodness-of-fit test based on Kendall's process\nerrors occured in", 
      #                numError, "of", B, "bootstraps which were suppressed.\nThe erroneous bootstraps were deleted. Note that this may cause erroneous p-values.\nThis may be an indicator for copula misspecification.\nMost probably Kendall's tau is close to zero.\nConsider an independence test first."))
      #}
      sn.boot <- sn.boot[!is.na(sn.boot)]
      tn.boot <- tn.boot[!is.na(tn.boot)]
      B.sn <- length(sn.boot)
      B.tn <- length(tn.boot)
      if(B.sn!=B.tn){stop("Number of bootstraps in CvM and KS based stat are not same !!!")}
      sn.obs <- ostat$Sn
      tn.obs <- ostat$Tn
      pv.sn <- sapply(sn.obs, function(x) (1/B.sn) * length(which(sn.boot[1:B.sn] >= 
                                                                    x)))
      pv.tn <- sapply(tn.obs, function(x) (1/B.tn) * length(which(tn.boot[1:B.tn] >= 
                                                                    x)))
      out <- list(p.value.CvM = pv.sn, p.value.KS = pv.tn, 
                  statistic.CvM = sn.obs, statistic.KS = tn.obs,
                  B_success=B.sn)
    }
  }
  else {
    stop("Method not implemented")
  }
  return(out)
}



###########################

# sub functions for the calculation of the Hessian matrix of the t-copula

f_rho <- function(u1, u2, par, par2) {
  a <- .C("diff2lPDF_rho_tCopula",
          as.double(u1),
          as.double(u2),
          as.integer(length(u1)),
          as.double(c(par, par2)),
          as.integer(2),
          as.double(rep(0, length(u1))),
          PACKAGE = "VineCopula")[[6]]
  
  return(sum(a))
}

f_nu <- function(u1, u2, par, par2) {
  a <- .C("diff2lPDF_nu_tCopula_new",
          as.double(u1),
          as.double(u2),
          as.integer(length(u1)),
          as.double(c(par, par2)),
          as.integer(2),
          as.double(rep(0, length(u1))),
          PACKAGE = "VineCopula")[[6]]
  
  return(sum(a))
}

f_rho_nu <- function(u1, u2, par, par2) {
  a <- .C("diff2lPDF_rho_nu_tCopula_new",
          as.double(u1),
          as.double(u2),
          as.integer(length(u1)),
          as.double(c(par, par2)),
          as.integer(2),
          as.double(rep(0, length(u1))),
          PACKAGE = "VineCopula")[[6]]
  
  return(sum(a))
}


#####################

# sub functions for the Kendall GOF


boot.stat <- function(u, v, fam, param) {
  
  n <- length(u)
  t <- seq(1, n)/(n + 1e-04)
  kt <- rep(0, n)
  
  # estimate paramemter for different copula family from (u,v)
  # param <- suppressWarnings({
  #     BiCopEst(u, v, family = fam)
  # })
  # calulate k(t) and kn(t) of bootstrap sample data
  if (fam == 1) {
    # normal
    sam <- BiCopSim(n, 1, param$par, param$par2)  # generate data for the simulation of K(t)
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # parameter estimation of sample data
    sim <- BiCopSim(10000, 1, sam.par$par, sam.par$par2)  # generate data for the simulation of theo. K(t)
    cormat <- matrix(c(1, param$par, param$par, 1), 2, 2)
    dcop <- rep(0, 10000)
    for (i in 1:10000) dcop[i] <- pmvnorm(upper = c(qnorm(sim[i, 1]),
                                                    qnorm(sim[i, 2])),
                                          corr = cormat)
    kt <- sapply(t,
                 function(x) (1/10000) * length(which(dcop[1:10000] <= x)))  # simulate K(t) of sample data
  }
  if (fam == 2) {
    # t
    sam <- BiCopSim(n, fam, param$par, param$par2)  # generate data for the simulation of K(t)
    sam.par <- suppressWarnings({
      BiCopEst(sam[, 1], sam[, 2], family = fam)
    })  # parameter estimation of sample data
    sim <- BiCopSim(10000, fam, sam.par$par, sam.par$par2)  # generate data for the simulation of theo. K(t)
    
    # par2 muss auf einen Integer gesetzt werden f?r mvtnorm
    param$par2 <- round(param$par2)
    
    cormat <- matrix(c(1, param$par, param$par, 1), 2, 2)
    dcop <- rep(0, 10000)
    for (i in 1:10000) dcop[i] <- pmvt(upper = c(qt(sim[i, 1], df = param$par2),
                                                 qt(sim[i, 2], df = param$par2)),
                                       corr = cormat,
                                       df = param$par2)
    kt <- sapply(t, function(x) (1/10000) * length(which(dcop[1:10000] <= x)))  # simulate K(t) of sample data
  } else if (fam == 3) {
    # Clayton
    sam <- BiCopSim(n, 3, param$par)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)$par  # estimate parameter of sample data
    kt <- t + t * (1 - t^sam.par)/sam.par
  } else if (fam == 4) {
    # gumbel
    sam <- BiCopSim(n, 4, param$par)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)$par  # estimate parameter of sample data
    kt <- t - t * log(t)/(sam.par)
  } else if (fam == 5) {
    # frank
    sam <- BiCopSim(n, 5, param$par)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)$par  # estimate parameter of sample data
    kt <- t + log((1 - exp(-sam.par))/(1 - exp(-sam.par * t))) * (1 - exp(-sam.par * t))/(sam.par * exp(-sam.par * t))
  } else if (fam == 6) {
    sam <- BiCopSim(n, 6, param$par)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)$par  # estimate parameter of sample data
    kt <- t - (log(1 - (1 - t)^sam.par) * (1 - (1 - t))^sam.par)/(sam.par * (1 - t)^(sam.par - 1))
  } else if (fam == 7) {
    # BB1
    sam <- BiCopSim(n, 7, param$par, param$par2)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # estimate parameter of sample data
    theta <- sam.par$par
    delta <- sam.par$par2
    kt <- t + 1/(theta * delta) * (t^(-theta) - 1)/(t^(-1 - theta))
  } else if (fam == 8) {
    # BB6
    sam <- BiCopSim(n, 8, param$par, param$par2)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # estimate parameter of sample data
    theta <- sam.par$par
    delta <- sam.par$par2
    kt <- t + log(-(1 - t)^theta + 1) * (1 - t - (1 - t)^(-theta) + (1 - t)^(-theta) * t)/(delta * theta)
  } else if (fam == 9) {
    # BB7
    sam <- BiCopSim(n, 9, param$par, param$par2)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # estimate parameter of sample data
    theta <- sam.par$par
    delta <- sam.par$par2
    kt <- t + 1/(theta * delta) * ((1 - (1 - t)^theta)^(-delta) - 1)/((1 - t)^(theta - 1) * (1 - (1 - t)^theta)^(-delta - 1))
  } else if (fam == 10) {
    # BB8
    sam <- BiCopSim(n, 10, param$par, param$par2)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # estimate parameter of sample data
    theta <- sam.par$par
    delta <- sam.par$par2
    kt <- t + log(((1 - t * delta)^theta - 1)/((1 - delta)^theta - 1)) * (1 - t * delta - (1 - t * delta)^(-theta) + (1 - t * delta)^(-theta) * t * delta)/(theta * delta)
  }
  
  # calculate emp. Kn
  w <- rep(0, n)
  w[1:n] <- mapply(function(x, y) (1/n) * length(which(x > sam[, 1] & y > sam[, 2])),
                   sam[, 1],
                   sam[, 2])
  w <- sort(w)
  kn <- rep(0, n)
  kn <- sapply(t, function(x) (1/n) * length(which(w[1:n] <= x)))
  
  # calculate test statistic Sn
  Sn1 <- 0
  Sn2 <- 0
  for (j in 1:(n - 1)) {
    Sn1 <- Sn1 + ((kn[j])^2 * (kt[j + 1] - kt[j]))
    Sn2 <- Sn2 + (kn[j]) * ((kt[j + 1])^2 - (kt[j])^2)
  }
  sn <- n/3 + n * Sn1 - n * Sn2
  # calculation of test statistics Tn
  tm <- matrix(0, n - 1, 2)
  # mit i=0
  for (j in 1:(n - 1)) {
    tm[j, 1] <- abs(kn[j] - kt[j])
  }
  # mit i=1
  for (j in 1:(n - 1)) {
    tm[j, 2] <- abs(kn[j] - kt[j + 1])
  }
  tn <- max(tm) * sqrt(n)
  
  sn <- sort(sn)  # vector of ordered statistic Sn
  tn <- sort(tn)  # vector of ordered statistic Tn
  out <- list(sn = sn, tn = tn)
}



obs.stat <- function(u, v, fam, param) {
  
  n <- length(u)
  t <- seq(1, n)/(n + 1e-04)
  kt <- rep(0, n)
  
  # estimate paramemter for different copula family from (u,v)
  # param <- suppressWarnings({
  #     BiCopEst(u, v, family = fam)
  # })
  
  # calculate observed K(t) of (u,v)
  kt.obs <- rep(0, n)
  if (fam == 1) {
    sim <- BiCopSim(10000, 1, param$par)  # generate data for the simulation of K(t)
    cormat <- matrix(c(1, param$par, param$par, 1), 2, 2)
    dcop <- rep(0, 10000)
    # TODO: for-loop as apply
    for (i in 1:10000) dcop[i] <- pmvnorm(upper = c(qnorm(sim[i, 1]),
                                                    qnorm(sim[i, 2])),
                                          corr = cormat)
    kt.obs <- sapply(t,
                     function(x) (1/10000) * length(which(dcop[1:10000] <= x)))  # simulate K(t) of sample data
  } else if (fam == 2) {
    sim <- BiCopSim(10000, 2, param$par, param$par2)  # generate data for the simulation of K(t)
    cormat <- matrix(c(1, param$par, param$par, 1), 2, 2)
    dcop <- rep(0, 10000)
    # TODO: for-loop as apply
    for (i in 1:10000) dcop[i] <- pmvt(upper = c(qt(sim[i, 1], df = param$par2),
                                                 qt(sim[i, 2], df = param$par2)),
                                       corr = cormat,
                                       df = param$par2)
    kt.obs <- sapply(t,
                     function(x) (1/10000) * length(which(dcop[1:10000] <= x)))  # simulate K(t) of sample data
  } else if (fam == 3) {
    kt.obs <- t + t * (1 - t^param$par)/param$par
  } else if (fam == 4) {
    kt.obs <- t - t * log(t)/(param$par)
  } else if (fam == 5) {
    kt.obs <- t + log((1 - exp(-param$par))/(1 - exp(-param$par * t))) * (1 - exp(-param$par * t))/(param$par * exp(-param$par * t))
  } else if (fam == 6) {
    kt.obs <- t - (log(1 - (1 - t)^param$par) * (1 - (1 - t))^param$par)/(param$par * (1 - t)^(param$par - 1))
  } else if (fam == 7) {
    theta <- param$par
    delta <- param$par2
    kt.obs <- t + 1/(theta * delta) * (t^(-theta) - 1)/(t^(-1 - theta))
  } else if (fam == 8) {
    theta <- param$par
    delta <- param$par2
    kt.obs <- t + log(-(1 - t)^theta + 1) * (1 - t - (1 - t)^(-theta) + (1 - t)^(-theta) * t)/(delta * theta)
  } else if (fam == 9) {
    theta <- param$par
    delta <- param$par2
    kt.obs <- t + 1/(theta * delta) * ((1 - (1 - t)^theta)^(-delta) - 1)/((1 - t)^(theta - 1) * (1 - (1 - t)^theta)^(-delta - 1))
  } else if (fam == 10) {
    theta <- param$par
    delta <- param$par2
    kt.obs <- t + log(((1 - t * delta)^theta - 1)/((1 - delta)^theta - 1)) * (1 - t * delta - (1 - t * delta)^(-theta) + (1 - t * delta)^(-theta) * t * delta)/(theta * delta)
  }
  # calculation of observed Kn
  w <- rep(0, n)
  w[1:n] <- mapply(function(x, y) (1/n) * length(which(x > u & y > v)), u, v)
  
  w <- sort(w)
  kn.obs <- rep(0, n)
  kn.obs <- sapply(t, function(x) (1/n) * length(which(w[1:n] <= x)))
  
  # calculation of observed value Sn
  Sn1 <- 0
  Sn2 <- 0
  for (j in 1:(n - 1)) {
    Sn1 <- Sn1 + ((kn.obs[j])^2 * (kt.obs[j + 1] - kt.obs[j]))
    Sn2 <- Sn2 + (kn.obs[j]) * ((kt.obs[j + 1])^2 - (kt.obs[j])^2)
  }
  Sn <- n/3 + n * Sn1 - n * Sn2  # observed Sn
  
  # calculation of observed Tn
  tn.obs <- matrix(0, n - 1, 2)
  # mit i=0
  for (j in 1:(n - 1)) {
    tn.obs[j, 1] <- abs(kn.obs[j] - kt.obs[j])
  }
  # mit i=1
  for (j in 1:(n - 1)) {
    tn.obs[j, 2] <- abs(kn.obs[j] - kt.obs[j + 1])
  }
  Tn <- max(tn.obs) * sqrt(n)
  out <- list(Sn = Sn, Tn = Tn)
  return(out)
}


############################


# boot.IR
#
# bootstrap for IR
#
# @param family copula family
# @param theta first copula parameter
# @param nu second copula parameter
# @param B number of bootstraps
# @param n Number of observations
#
# @return IR vector of test statistics
#
# @author Ulf Schepsmeier
#

boot.IR <- function(family, theta, nu, B, n) {
  # theta und nu sind die geschaetzten Parameter
  IR <- rep(0, B)
  # TODO: for-loop as apply
  for (i in 1:B) {
    sam <- BiCopSim(n, family, theta, nu)
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = family)  # parameter estimation of sample data
    if (family == 2) {
      theta2 <- sam.par[1]
      nu2 <- sam.par[2]
      grad <- c(0, 0)
      rho_teil <- f_rho(sam[, 1], sam[, 2], theta2, nu2)
      nu_teil <- f_nu(sam[, 1], sam[, 2], theta2, nu2)
      rho_nu_teil <- f_rho_nu(sam[, 1], sam[, 2], theta2, nu2)
      H <- matrix(c(rho_teil, rho_nu_teil, rho_nu_teil, nu_teil), 2, 2)  # Hesse matrix
      grad[1] <- BiCopDeriv(sam[, 1],
                            sam[, 2],
                            family = family,
                            par = theta2,
                            par2 = nu2,
                            deriv = "par",
                            log = TRUE,
                            check.pars = FALSE)
      grad[2] <- BiCopDeriv(sam[, 1],
                            sam[, 2],
                            family = family,
                            par = theta2,
                            par2 = nu2,
                            deriv = "par2",
                            log = TRUE,
                            check.pars = FALSE)
      C <- grad %*% t(grad)
    } else {
      theta2 <- sam.par
      nu2 <- 0
      d <- rep(0, T)
      # TODO: for-loop as apply
      for (t in 1:T) {
        b <- BiCopPDF(sam[t, 1], sam[t, 2], family, theta2, nu2)
        d[t] <- BiCopDeriv2(sam[t, 1],
                            sam[t, 2],
                            family = family,
                            par = theta2,
                            par2 = nu2,
                            deriv = "par",
                            check.pars = FALSE)/b
      }
      H <- mean(d)
      C <- BiCopDeriv(sam[, 1],
                      sam[, 2],
                      family = family,
                      par = theta2,
                      par2 = nu2,
                      deriv = "par",
                      log = TRUE,
                      check.pars = FALSE)
    }
    Phi <- -solve(H) %*% C
    IR[i] <- trace(Phi)/dim(H)[1]
  }
  
  return(IR)
}



## sub-functions

# hesseTcopula
#
# This small function calculates the Hessian matrix for the t-copula
#
# @param u1 first copula argument
# @param u2 second copula argument
# @param theta first copula parameter
# @param nu second copula parameter
#
# @return H Hesse matrix for the t-copula
#
# @author Ulf Schepsmeier
#

hesseTcopula <- function(u1, u2, theta, nu){
  rho_teil <- f_rho(u1, u2, theta, nu)
  nu_teil <- f_nu(u1, u2, theta, nu)
  rho_nu_teil <- f_rho_nu(u1, u2, theta, nu)
  H <- matrix(c(rho_teil, rho_nu_teil, rho_nu_teil, nu_teil), 2, 2)
}



# OPGtcopula
#
# This small function calculates the outer product of gradient for the t-copula
#
# @param u1 first copula argument
# @param u2 second copula argument
# @param family copula family (here Student's t copula = 2)
# @param theta first copula parameter
# @param nu second copula parameter
#
# @return C outer product of gradient
#
# @author Ulf Schepsmeier
#

OPGtcopula <- function(u1, u2, family, theta, nu){
  grad <- numeric(2)
  # gradient
  grad[1] <- mean(BiCopDeriv(u1,
                             u2,
                             family = family,
                             par = theta,
                             par2 = nu,
                             deriv = "par",
                             log = TRUE,
                             check.pars = FALSE))
  grad[2] <- mean(BiCopDeriv(u1,
                             u2,
                             family = family,
                             par = theta,
                             par2 = nu,
                             deriv = "par2",
                             log = TRUE,
                             check.pars = FALSE))
  
  ## outer product of gradient
  C <- grad %*% t(grad)
}


# gradDtcopula
#
# derivative of D (i.e. gradD) for the t-copula
#
# @param u1 first copula argument
# @param u2 second copula argument
# @param theta first copula parameter
# @param nu second copula parameter
#
# @return gradD gradient of D
#
# @author Ulf Schepsmeier
#

gradDtcopula <- function(u1, u2, theta, nu){
  eps <- 0.001
  H_theta_eps_plus <- hesseTcopula(u1, u2, theta+eps, nu)
  H_theta_eps_minus <- hesseTcopula(u1, u2, theta-eps, nu)
  H_nu_eps_plus <- hesseTcopula(u1, u2, theta, nu+eps)
  H_nu_eps_minus <- hesseTcopula(u1, u2, theta, nu-eps)
  C_theta_eps_plus <- OPGtcopula(u1, u2, family=2, theta=(theta+eps), nu=nu)
  C_theta_eps_minus <- OPGtcopula(u1, u2, family=2, theta=(theta-eps), nu=nu)
  C_nu_eps_plus <- OPGtcopula(u1, u2, family=2, theta=theta, nu=(nu+eps))
  C_nu_eps_minus <- OPGtcopula(u1, u2, family=2, theta=theta, nu=(nu-eps))
  
  Hprime_theta_eps_plus <- as.vector(H_theta_eps_plus[lower.tri(H_theta_eps_plus, diag = TRUE)])
  Hprime_theta_eps_minus <- as.vector(H_theta_eps_minus[lower.tri(H_theta_eps_minus, diag = TRUE)])
  Hprime_nu_eps_plus <- as.vector(H_nu_eps_plus[lower.tri(H_nu_eps_plus, diag = TRUE)])
  Hprime_nu_eps_minus <- as.vector(H_nu_eps_minus[lower.tri(H_nu_eps_minus, diag = TRUE)])
  Cprime_theta_eps_plus <- as.vector(C_theta_eps_plus[lower.tri(C_theta_eps_plus, diag = TRUE)])
  Cprime_theta_eps_minus <- as.vector(C_theta_eps_minus[lower.tri(C_theta_eps_minus, diag = TRUE)])
  Cprime_nu_eps_plus <- as.vector(C_nu_eps_plus[lower.tri(C_nu_eps_plus, diag = TRUE)])
  Cprime_nu_eps_minus <- as.vector(C_nu_eps_minus[lower.tri(C_nu_eps_minus, diag = TRUE)])
  
  Dprime_theta_eps_plus <- Hprime_theta_eps_plus + Cprime_theta_eps_plus
  Dprime_theta_eps_minus <- Hprime_theta_eps_minus + Cprime_theta_eps_minus
  Dprime_nu_eps_plus <- Hprime_nu_eps_plus + Cprime_nu_eps_plus
  Dprime_nu_eps_minus <- Hprime_nu_eps_minus + Cprime_nu_eps_minus
  
  gradD_theta <- (Dprime_theta_eps_minus - Dprime_theta_eps_plus)/(2*eps)
  gradD_nu <- (Dprime_nu_eps_minus - Dprime_nu_eps_plus)/(2*eps)
  
  gradD <- cbind(gradD_theta, gradD_nu)
}





# bootWhite
#
# This small function provides the code to calculated bootstrapped p-values
# for the White test.
#
# @param family copula family
# @param theta first copula parameter
# @param nu second copula parameter
# @param B number of bootstraps
# @param N number of observations
#
# @return testStat
#
# @author Ulf Schepsmeier
#

bootWhite <- function(family, theta, nu, B, N){
  testStat <- rep(0, B)
  numError <- 0
  # TODO: for-loop as apply
  for(i in 1:B){
    ax <- try({
      sam <- BiCopSim(N, family, theta, nu)
      sam.par <- BiCopEst(sam[, 1], sam[, 2], family = family)  # parameter estimation of sample data
      testStat[i] <- BiCopGofTest(u1 = sam[,1], u2 = sam[,2], family = family,
                                  par = sam.par$par, par2 = sam.par$par2, method = "White",
                                  B = 0)$statistic
    }, silent = TRUE)
    if(inherits(ax, "try-error")){
      testStat[i] <- NA
      numError <- numError + 1
    }
  }
  if(numError > 0){
    warning(paste("In the calculation of the p-values for the copula
goodness-of-fit test based on White's test
errors occured in ", numError, "of", B, "bootsrtaps which were suppressed.
The erroneous bootstraps were deleted. Note that this may cause erroneous p-values.
This may be an indicator for copula misspecification.
Most probably Kendall's tau is close to zero.
Consider an independence test first."))
  }
  
  return(testStat)
}