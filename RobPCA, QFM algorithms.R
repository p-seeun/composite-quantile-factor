library(Matrix)
library(expm)

loss.fun <-function(tau, vec){
  loss = sapply(vec, function(x){ return(ifelse( x>0, tau*x, (tau-1)*x)) }) 
  return( mean(loss) )
}

# QFM estimation by Chen et al.(2021)
qfm_est <- function(X, r, tol, tau) {
  Te <- nrow(X)
  N <- ncol(X)
  
  F0 <- matrix(rnorm(Te * r), nrow = Te)
  
  obj <- 0
  newobj <- 1
  
  while (abs(newobj - obj) > tol) {
    lambda0 <- matrix(0, nrow=N, ncol=r)
    for (i in 1:N) {
      lambda0[i,]<-rq(X[,i] ~ F0-1, tau = tau)$coeff
    }
    obj <- loss.fun(tau, as.vector(X - F0 %*% t(lambda0)))
    
    F1 <- matrix(0,nrow = Te, ncol = r)
    for (j in 1:Te) {
      F1[j,] <- rq(X[j,] ~ lambda0-1, tau=tau)$coeff
    }
    F0<-F1  
    
    lambda1 <- matrix(nrow = N, ncol = r)
    for (i in 1:N) {
      lambda1[i,]<-rq(X[,i] ~ F1-1, tau = tau)$coeff
    }
    newobj <- (loss.fun(tau, as.vector(X - F1 %*% t(lambda1))))
  }
  
  Fhat <- F1
  Lhat <- lambda1
  
  sigmaF <- t(Fhat) %*% Fhat / Te
  sigmaA <- t(Lhat) %*% Lhat / N
  
  dum1 <- sqrtm(sigmaF) %*% sigmaA %*% sqrtm(sigmaF)
  dum2 <- eigen(dum1)$vectors
  
  R<- solve(sqrtm(sigmaF), dum2)
  Fhat <- Fhat %*% R
  Lhat <- Lhat %*% t(solve(R))
  
  return(list(fmat = Fhat, lmat = Lhat))
}

# Robust PCA estimation by He et al.(2022)
robpca_est <- function(X, r) {
  Te <- nrow(X)
  N <- ncol(X)
  K <- matrix(0, nrow = N, ncol = N)
  
  for (p in 2:Te) {
    for (q in 1:(p-1)) {
      xp <- X[p,]
      xq <- X[q,]
      K <- K + ((xp-xq) %*% t(xp-xq)) / sum((xp-xq)^2)
    }
  }
  K <- K*2 / N / (N - 1)
  dum <- eigen(K)
  Lhat <- dum$vectors[, 1:r] * sqrt(N)
  Fhat <- X %*% Lhat / N
  return(list(lmat = Lhat, fmat = Fhat))
}



