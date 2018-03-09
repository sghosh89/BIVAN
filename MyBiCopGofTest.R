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
      if (numError > 0) {
        warning(paste("In the calculation of the p-values for the copula\ngoodness-of-fit test based on Kendall's process\nerrors occured in", 
                      numError, "of", B, "bootstraps which were suppressed.\nThe erroneous bootstraps were deleted. Note that this may cause erroneous p-values.\nThis may be an indicator for copula misspecification.\nMost probably Kendall's tau is close to zero.\nConsider an independence test first."))
      }
      sn.boot <- sn.boot[!is.na(sn.boot)]
      tn.boot <- tn.boot[!is.na(tn.boot)]
      B.sn <- length(sn.boot)
      B.tn <- length(tn.boot)
      sn.obs <- ostat$Sn
      tn.obs <- ostat$Tn
      pv.sn <- sapply(sn.obs, function(x) (1/B.sn) * length(which(sn.boot[1:B.sn] >= 
                                                                    x)))
      pv.tn <- sapply(tn.obs, function(x) (1/B.tn) * length(which(tn.boot[1:B.tn] >= 
                                                                    x)))
      out <- list(p.value.CvM = pv.sn, p.value.KS = pv.tn, 
                  statistic.CvM = sn.obs, statistic.KS = tn.obs)
    }
  }
  else {
    stop("Method not implemented")
  }
  return(out)
}