source("CQFM algorithm.R")
source("RobPCA, QFM algorithms.R")

X = readRDS("./Real data analysis/Fred-MD/FRED-MD.rds")

N=ncol(X)
Te=nrow(X)

## Estimate number of factors
set.seed(1)

# CQFM_unif
CQFM_unif = cqfm_est(X, r=15, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(1,1,1,1,1,1,1))
lmat = CQFM_unif$lmat
rhats = lapply(lmat, function(L){ x= eigen(t(L)%*%L/N)$values; return(sum( x > min(N,Te)^(-1/3)*x[1] ))}) 
rhat = max(unlist(rhats))
rhat # 7

#QFM01
QFM01 = qfm_est(X, 15, tol=1e-5, tau=0.01)
lmat = QFM01$lmat
dum = t(lmat) %*% lmat / N
rhat <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1])  
rhat # =1

QFM1 = qfm_est(X, 15, tol=1e-5, tau=0.1)
lmat = QFM1$lmat
dum = t(lmat) %*% lmat / N
rhat <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1]) 
rhat # =2

QFM3 = qfm_est(X, 15, tol=1e-5, tau=0.3)
lmat = QFM3$lmat
dum = t(lmat) %*% lmat / N
rhat <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1]) 
rhat # =6

QFM5 = qfm_est(X, 15, tol=1e-5, tau=0.5)
lmat = QFM5$lmat
dum = t(lmat) %*% lmat / N
rhat <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1]) 
rhat # =6

QFM7 = qfm_est(X, 15, tol=1e-5, tau=0.7)
lmat = QFM7$lmat
dum = t(lmat) %*% lmat / N
rhat <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1]) 
rhat # =7

QFM9 = qfm_est(X, 15, tol=1e-5, tau=0.9)
lmat = QFM9$lmat
dum = t(lmat) %*% lmat / N
rhat <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1]) 
rhat # =2

QFM99 = qfm_est(X, 15, tol=1e-5, tau=0.99)
lmat = QFM99$lmat
dum = t(lmat) %*% lmat / N
rhat <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1]) #log(N*Te/(N+Te))*(N+Te)/N/Te 
rhat # =1

# CQFM_low
CQFM_low = cqfm_est(X, r=15, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(2,2,1,1,1,1,1))
lmat = CQFM_low$lmat
rhats = lapply(lmat, function(L){ x= eigen(t(L)%*%L/N)$values; return(sum( x > min(N,Te)^(-1/3)*x[1] ))}) 
rhat = max(unlist(rhats))
rhat # 7

# CQFM_med
CQFM_med = cqfm_est(X, r=15, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(1,1,2,2,2,1,1))
lmat = CQFM_med$lmat
rhats = lapply(lmat, function(L){ x= eigen(t(L)%*%L/N)$values; return(sum( x > min(N,Te)^(-1/3)*x[1] ))}) 
rhat = max(unlist(rhats))
rhat # 7

# CQFM_high
CQFM_high = cqfm_est(X, r=15, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(1,1,1,1,1,2,2))
lmat = CQFM_high$lmat
rhats = lapply(lmat, function(L){ x= eigen(t(L)%*%L/N)$values; return(sum( x > min(N,Te)^(-1/3)*x[1] ))}) 
rhat = max(unlist(rhats))
rhat # 7

# Mean FM
icp <- function(r){
  svd = svd(X) 
  lambda = t(as.matrix(sqrt(N)*t(svd$v)[1:r,])) 
  Factor = (X %*% lambda)/N 
  IC = log( sum((X - Factor%*%t(lambda))^2) /(N*Te) )+ r * (N+Te)/(N*Te) * log(N*Te/(N+Te) ) #icp1
  return(IC)
}
icp.results = c()
for(r in 2:10){
  icp.results[r] = icp(r)
}
rhat = which.min(icp.results)
rhat # 8

# Huang(2023)
tau.vec = c(0.01,0.1,0.3,0.5,0.7,0.9,0.99)
icp <- function(r){
  cqrfactor = cqrfactor(Z, r=r, tau=tau.vec)
  Factor = cqrfactor$fmat
  Lambda = cqrfactor$lmat
  loss = 0
  for(k in 1:length(tau.vec)){
    loss = loss + loss.fun(tau.vec[k], Z - Factor %*% t(Lambda)-cqrfactor$quantiles[k])
  }
  IC = log(loss) + r*(N+Te)/(N*Te)*log( N*Te/(N+Te) )
  return(IC)
}

icp.results = c()
for(r in 1:10){
  icp.results[r] = icp(r)
  print(r)
}
rhat = which.min(icp.results)
rhat # 6
