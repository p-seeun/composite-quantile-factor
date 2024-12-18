#Location-shift model: t(2)

source("CQFM algorithm.R")
source("RobPCA, QFM algorithms.R")
library(cqrfactor)

rsquare <- function(a, B){
  lm<- summary(lm(a ~ B))
  return(lm$adj.r.squared)
}

## Rob.PCA and CQFM
RobPCA = list()
Huang = list()
CQFM = list()

for(n in 1:5){
  if(n==1){i=1; j=1}
  if(n==2){i=2; j=1}
  if(n==3){i=2; j=2}
  if(n==4){i=2; j=3}
  if(n==5){i=3; j=3}
  
  # Set the parameters
  N = 50*i
  Te = 50*j
  rep = 100
  tol = 1e-5
  b = c(0.8, 0.5, 0.2)
  
  # Initialize vectors to store results
  R2_RobPCA = matrix(NA, nrow = rep, ncol = 3)
  R2_Huang = matrix(NA, nrow = rep, ncol = 3)
  R2_CQFM = matrix(NA, nrow = rep, ncol = 3)
  
  for (s in 1:rep) {
    set.seed(s)
    
    # Generate data
    F <- matrix(0, nrow = Te, ncol = 3)
    u <- matrix(rnorm(Te * 3), nrow = Te, ncol = 3)
    F[1, ] <- u[1, ]
    for (h in 2:Te) {
      F[h, 1] <- b[1] * F[h - 1, 1] + u[h, 1]
      F[h, 2] <- b[2] * F[h - 1, 2] + u[h, 2]
      F[h, 3] <- b[3] * F[h - 1, 3] + u[h, 3]
    }
    Lambda <- matrix(rnorm(N * 3), nrow = N, ncol = 3)
    
    e =  matrix(rt(Te * N, 2), nrow = Te, ncol = N)
    
    X <- F %*% t(Lambda) + e
    weight = dt(qt(c(0.1,0.3,0.5,0.7,0.9),df=2),df=2)
    
    # Estimate factors
    Fhat_RobPCA = robpca_est(X, r=3)$fmat
    Fhat_Huang = cqrfactor(X, r=3, tau=c(0.1,0.3,0.5,0.7,0.9), tol=1e-4)$fmat
    Fhat_CQFM = cqfm_est(X, r=4, tau.vec=c(0.1,0.3,0.5,0.7,0.9), tol=1e-4, weight=weight)$fmat
    
    r2_RobPCA = apply(F, 2, function(x) rsquare(x, Fhat_RobPCA))
    r2_Huang = apply(F, 2, function(x) rsquare(x, Fhat_Huang))
    r2_CQFM = apply(F, 2, function(x) rsquare(x, Fhat_CQFM))
    
    # Store results
    R2_RobPCA[s, ] <- r2_RobPCA
    R2_Huang[s, ] <- r2_Huang
    R2_CQFM[s, ] <- r2_CQFM
    
    print(s)
  }
  Mean_R2_RobPCA = colMeans(R2_RobPCA)
  Mean_R2_Huang = colMeans(R2_Huang)
  Mean_R2_CQFM = colMeans(R2_CQFM)
  
  RobPCA[[n]] = Mean_R2_RobPCA
  Huang[[n]] = Mean_R2_Huang
  CQFM[[n]] = Mean_R2_CQFM
  print(n)
}

## QFMs 
QFM0.1 = list()
QFM0.3 = list()
QFM0.5 = list()
QFM0.7 = list()
QFM0.9 = list()

for(n in 1:5){
  if(n==1){i=1; j=1}
  if(n==2){i=2; j=1}
  if(n==3){i=2; j=2}
  if(n==4){i=2; j=3}
  if(n==5){i=3; j=3}
  
  # Set the parameters
  N = 50*i
  Te = 50*j
  rep = 100
  tol = 1e-5
  b = c(0.8, 0.5, 0.2)
  
  # Initialize vectors to store results
  R2_QFM0.1 = matrix(NA, nrow = rep, ncol = 3)
  R2_QFM0.3 = matrix(NA, nrow = rep, ncol = 3)
  R2_QFM0.5 = matrix(NA, nrow = rep, ncol = 3)
  R2_QFM0.7 = matrix(NA, nrow = rep, ncol = 3)
  R2_QFM0.9 = matrix(NA, nrow = rep, ncol = 3)
  
  for (s in 1:rep) {
    set.seed(s)
    
    # Generate data
    F <- matrix(0, nrow = Te, ncol = 3)
    u <- matrix(rnorm(Te * 3), nrow = Te, ncol = 3)
    F[1, ] <- u[1, ]
    for (h in 2:Te) {
      F[h, 1] <- b[1] * F[h - 1, 1] + u[h, 1]
      F[h, 2] <- b[2] * F[h - 1, 2] + u[h, 2]
      F[h, 3] <- b[3] * F[h - 1, 3] + u[h, 3]
    }
    Lambda <- matrix(rnorm(N * 3), nrow = N, ncol = 3)
    
    e =  matrix(rt(Te * N, 2), nrow = Te, ncol = N)
    
    X <- F %*% t(Lambda) + e
    
    # Estimate factors
    Fhat_QFM0.1 = qfm_est(X, 4, tol, tau=0.1)$fmat
    Fhat_QFM0.3 = qfm_est(X, 4, tol, tau=0.3)$fmat
    Fhat_QFM0.5 = qfm_est(X, 3, tol, tau=0.5)$fmat
    Fhat_QFM0.7 = qfm_est(X, 4, tol, tau=0.7)$fmat
    Fhat_QFM0.9 = qfm_est(X, 4, tol, tau=0.9)$fmat
    
    r2_QFM0.1 = apply(F, 2, function(x) rsquare(x, Fhat_QFM0.1))
    r2_QFM0.3 = apply(F, 2, function(x) rsquare(x, Fhat_QFM0.3))
    r2_QFM0.5 = apply(F, 2, function(x) rsquare(x, Fhat_QFM0.5))
    r2_QFM0.7 = apply(F, 2, function(x) rsquare(x, Fhat_QFM0.7))
    r2_QFM0.9 = apply(F, 2, function(x) rsquare(x, Fhat_QFM0.9))
    
    # Store results
    
    R2_QFM0.1[s, ] <- r2_QFM0.1
    R2_QFM0.3[s, ] <- r2_QFM0.3
    R2_QFM0.5[s, ] <- r2_QFM0.5
    R2_QFM0.7[s, ] <- r2_QFM0.7
    R2_QFM0.9[s, ] <- r2_QFM0.9
    
    print(s)
  }
  Mean_R2_QFM0.1 = colMeans(R2_QFM0.1)
  Mean_R2_QFM0.3 = colMeans(R2_QFM0.3)
  Mean_R2_QFM0.5 = colMeans(R2_QFM0.5)
  Mean_R2_QFM0.7 = colMeans(R2_QFM0.7)
  Mean_R2_QFM0.9 = colMeans(R2_QFM0.9)
  
  QFM0.1[[n]] = Mean_R2_QFM0.1
  QFM0.3[[n]] = Mean_R2_QFM0.3
  QFM0.5[[n]] = Mean_R2_QFM0.5
  QFM0.7[[n]] = Mean_R2_QFM0.7
  QFM0.9[[n]] = Mean_R2_QFM0.9
}

# RobPCA[[n]], QFM0.1[[n]], ..., QFM0.9[[n]], Huang[[n]], CQFM[[n]] include results for n-th combination of (N,T).
Table1_t = matrix(nrow=8, ncol=0)
for(n in 1:5){
  table = rbind(RobPCA[[n]], QFM0.1[[n]], QFM0.3[[n]], QFM0.5[[n]], QFM0.7[[n]], QFM0.9[[n]], Huang[[n]], CQFM[[n]])
  Table1_t = cbind(Table1_t, table)
}

Table1_t # first table in Table 1
