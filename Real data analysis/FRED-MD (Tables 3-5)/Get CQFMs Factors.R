source("CQFM algorithm.R")
source("RobPCA, QFM algorithms.R")

X = readRDS("./Real data analysis/Fred-MD/FRED-MD.rds")

## CQFM_unif
f.CQFM_unif = list()
for(st in 1:(nrow(X)-120)){

  Win = X[st:(st+120-1),]
  N = ncol(Win)
  Te = nrow(Win)
  
  set.seed(st)
  f.CQFM_unif[[st]] = cqfm_est(Win, r=7, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(1,1,1,1,1,1,1))$fmat
  
  print(st)
}

#saveRDS(list(f.CQFM=f.CQFM_unif), "./Real data analysis/Fred-MD/CQFM_unif.rds")

## CQFM_low
f.CQFM_low = list()
for(st in 1:(nrow(X)-120)){
  
  Win = X[st:(st+120-1),]
  N = ncol(Win)
  Te = nrow(Win)
  
  set.seed(st)
  f.CQFM_low[[st]] = cqfm_est(Win, r=7, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(2,2,1,1,1,1,1))$fmat
  
  print(st)
}

#saveRDS(list(f.CQFM=f.CQFM_low), "./Real data analysis/Fred-MD/CQFM_low.rds")

## CQFM_med
f.CQFM_med = list()
for(st in 1:(nrow(X)-120)){
  
  Win = X[st:(st+120-1),]
  N = ncol(Win)
  Te = nrow(Win)
  
  set.seed(st)
  f.CQFM_med[[st]] = cqfm_est(Win, r=7, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(1,1,2,2,2,1,1))$fmat
  
  print(st)
}

#saveRDS(list(f.CQFM=f.CQFM_med), "./Real data analysis/Fred-MD/CQFM_med.rds")

## CQFM_high
f.CQFM_high = list()
for(st in 1:(nrow(X)-120)){
  
  Win = X[st:(st+120-1),]
  N = ncol(Win)
  Te = nrow(Win)
  
  set.seed(st)
  f.CQFM_high[[st]] = cqfm_est(Win, r=7, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(1,1,1,1,1,2,2))$fmat
  
  print(st)
}

#saveRDS(list(f.CQFM=f.CQFM_high, "./Real data analysis/Fred-MD/CQFM_high.rds")
