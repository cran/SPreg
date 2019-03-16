skewProbit <- function(formula, data = list(), penalty = "Jeffrey", initial = NULL, cvtCov = TRUE, delta0 = 3, level = 0.95){
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  
  est <- skewProbit.fit(y, x, penalty, initial, cvtCov, delta0, level)
  est$call <- match.call()
  est$formula <- formula
  class(est) <- "skewProbit"
  est
}

print.skewProbit <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

summary.skewProbit <- function(object, ...){
  
  TAB <- cbind(Estimate = coef(object),
               StdErr = object$stderr,
               lower = object$lower,
               upper = object$upper,
               zscore = object$zscore,
               p.value = object$pval
  )
  colnames(TAB)[3:4] <- c(paste("L", (1-object$level)/2, "%", sep = ""),
                          paste("U", 1-(1-object$level)/2, "%", sep = ""))
  res <- list(call = object$call, coefficients = TAB)
  class(res) <- "summary.skewProbit"
  res
  
}

print.summary.skewProbit <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
}

skewProbit.fit <- function(y, x, penalty = "Jeffrey", initial = NULL, cvtCov = TRUE, delta0 = 3, level = 0.95){
  
  ####
  n <- length(y)
  if (is.null(n1 <- nrow(x))) stop("'x' must be a matrix")
  if(n1 == 0L) stop("0 (non-NA) cases")
  
  if(dim(x)[1] != n) stop("The number of subjects in y and x are different")
  
  if(is.null(initial)){
    initial <- coef(glm(y ~ x, family = binomial(link = "probit")))
    initial <- initial[!is.na(initial)]
  } else{
    if(length(initial) != dim(x)[2]) stop("initial guess should be the same length of covariate matrix x")
  }
  
  if(cvtCov){
    idx_cvt <- which(apply(x, 2, function(x) length(table(x))) >= n/2)
    for(i in idx_cvt){
      x[,i] <- scale(x[,i], TRUE, TRUE)
    }
  }
  
  p <- ncol(x)
  
  u_p <- 1-(1-level)/2
  l_p <- (1-level)/2
  ####
  if(penalty == "Naive"){
    delta0.tmp <- delta0 
    while(1){
      fit_naive <- try(ucminf(c(initial, delta0.tmp), fn = loglik, X = x, y = y, hessian = 2), silent=TRUE)
      if(class(fit_naive) != "try-error"){
        if(fit_naive$convergence == 1 || fit_naive$convergence == 2 || fit_naive$convergence == 4) break
      }
      delta0.tmp <- runif(1, 1, 10)
    }
    Nest <- fit_naive$par
    Ninv <- try(suppressWarnings(sqrt(diag(fit_naive$invhessian))), silent = TRUE)
    if("try-error" %in% class(Ninv)){ Nse <- rep(NA, p+1) } else {
      Nse <- Ninv
    }	
    coef <- cbind(Nest, Nse, Nest/Nse, 2*(1-pnorm(abs(Nest/Nse))), Nest + qnorm(l_p)*Nse,  Nest + qnorm(u_p)*Nse)
  } else if(penalty == "Jeffrey"){
    delta0.tmp <- delta0
    while(1){
      fit_Jeff <- try(ucminf(c(initial, delta0.tmp), fn = Jloglikp, X = x, y = y, hessian = 2), silent=TRUE)
      if(class(fit_Jeff) != "try-error"){						
        if(fit_Jeff$convergence == 1 || fit_Jeff$convergence == 2 || fit_Jeff$convergence == 4) break
      }
      delta0.tmp <- runif(1, 1, 10)
    }
    Jest <- fit_Jeff$par
    Jinv <- try(suppressWarnings(sqrt(diag(fit_Jeff$invhessian))), silent = TRUE)
    if("try-error" %in% class(Jinv)){ Jse <- rep(NA, p+1) } else {
      Jse <- Jinv
    }
    coef <- cbind(Jest, Jse, Jest/Jse, 2*(1-pnorm(abs(Jest/Jse))), Jest + qnorm(l_p)*Jse, Jest + qnorm(u_p)*Jse)
  } else if(penalty == "Cauchy"){
    delta0.tmp <- delta0
    while(1){
      fit_Cauchy <- try(ucminf(c(initial, delta0.tmp), fn = Cloglikp, X = x, y = y, hessian = 2), silent=TRUE)
      if(class(fit_Cauchy) != "try-error"){	
        if(fit_Cauchy$convergence == 1 || fit_Cauchy$convergence == 2 || fit_Cauchy$convergence == 4) break
      }
      delta0.tmp <- runif(1, 1, 10)
    }	
    Cest <- fit_Cauchy$par
    Cinv <- try(suppressWarnings(sqrt(diag(fit_Cauchy$invhessian))), silent = TRUE)
    if("try-error" %in% class(Cinv)){ Cse <- rep(NA, p+1) } else {
      Cse <- Cinv
    }
    coef <- cbind(Cest, Cse, Cest/Cse, 2*(1-pnorm(abs(Cest/Cse))), Cest + qnorm(l_p)*Cse, Cest + qnorm(u_p)*Cse)	
  }
  rownames(coef)[1] <- "(Intercept)"
  rownames(coef)[2:ncol(x)] <- colnames(x)[-1]
  rownames(coef)[nrow(coef)] <- "delta"
  re <- list(coefficients = coef[,1],
             stderr = coef[,2],
             zscore = coef[,3],
             pval = coef[,4],
             lower = coef[,5],
             upper = coef[,6],
             level = level)
  return(re)
  
  
}


##############################################################################################################

loglik <- function(paras, X, y){
  delta <- paras[length(paras)]
  eta1 <- as.vector(X %*% paras[-length(paras)])
  mu1 <- psn(as.vector(eta1), alpha = delta)
  mu1[which(mu1==0)] <- min(mu1[mu1 != 0])
  mu1[which(mu1==1)] <- max(mu1[mu1 != 1])
  re <- sum(y*log(mu1) + (1-y)*log(1-mu1))
  return(-re)
}

## Faster than derlik1 below
derlik <- function(paras, X, y){
  n <- length(y)
  delta <- paras[length(paras)]
  eta1 <- as.vector(X %*% paras[-length(paras)])
  mu1 <- psn(as.vector(eta1), alpha = delta)
  mu1[which(mu1==0)] <- min(mu1[mu1 != 0])
  mu1[which(mu1==1)] <- max(mu1[mu1 != 1])
  
  re <- rep(0, length(paras))
  for(i in 1:n){
    re <- re + ((y[i] - mu1[i])/(mu1[i]*(1-mu1[i])))*
      c(-2*dnorm(eta1[i])*pnorm(delta*eta1[i])*X[i,], exp(-(1+delta^2)*eta1[i]^2/2)/(pi*(1+delta^2)))
  }
  return(re)
}


Jloglikp <- function(paras, X, y){
  inf.mat=matrix(0, nrow=length(paras), ncol=length(paras))
  delta=paras[length(paras)]
  eta1 <- as.vector(X %*% paras[-length(paras)])
  mu1 <- psn(as.vector(eta1), alpha = delta)
  mu1 <- psn(as.vector(eta1), alpha = paras[length(paras)])
  #mu1[which(mu1==0)] <- min(mu1[mu1 != 0])
  mu1[which(mu1==0)] <- 10e-10
  #mu1[which(mu1==1)] <- max(mu1[mu1 != 1])
  mu1[which(mu1==1)] <- 1-10e-10
  term0= dnorm(eta1)*pnorm(delta*eta1)
  term1=mu1*(1-mu1)
  term2=term0*term0/term1
  term3=exp(-0.5*eta1^2*(1+delta^2))
  term4=term3*term3
  inf.mat.b <- 4*t(term2*X) %*% X
  inf.mat.bd <- -2*colSums((term0*term3/term1)*X)/(pi*(1+delta^2))
  inf.mat.d <- sum(term4/term1)/(pi*(1+delta^2))^2
  inf.mat[-length(paras), -length(paras)] <- inf.mat.b
  inf.mat[length(paras), length(paras)] <- inf.mat.d
  inf.mat[length(paras), -length(paras)] <- inf.mat.bd
  inf.mat[-length(paras), length(paras)] <- inf.mat.bd	
  if(det(inf.mat) < 0) qnty1=0 else qnty1=0.5*log(det(inf.mat))
  #######
  re <- sum(y*log(mu1) + (1-y)*log(1-mu1)) + qnty1
  return(-re)
}

Cloglikp <- function(paras, X, y){
  delta=paras[length(paras)]
  eta1 <- as.vector(X %*% paras[-length(paras)])
  mu1 <- psn(as.vector(eta1), alpha = delta)
  mu1[which(mu1==0)] <- min(mu1[mu1 != 0])
  mu1[which(mu1==1)] <- max(mu1[mu1 != 1])
  re <- sum(y*log(mu1) + (1-y)*log(1-mu1)) - sum(log(1+paras^2/2.5^2))
  return(-re)
}
