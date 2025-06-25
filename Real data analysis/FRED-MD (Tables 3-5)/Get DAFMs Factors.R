source("DAFM algorithm.R")
source("RobPCA, QFM algorithms.R")

X = readRDS("./Real data analysis/Fred-MD/FRED-MD.rds")

## DAFM_unif
f.DAFM_unif = list()
for(st in 1:(nrow(X)-120)){

  Win = X[st:(st+120-1),]
  N = ncol(Win)
  Te = nrow(Win)
  
  set.seed(st)
  f.DAFM_unif[[st]] = dafm_est(Win, r=7, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(1,1,1,1,1,1,1))$fmat
  
  print(st)
}

#saveRDS(list(f.DAFM=f.DAFM_unif), "./Real data analysis/Fred-MD/DAFM_unif.rds")

## DAFM_low
f.DAFM_low = list()
for(st in 1:(nrow(X)-120)){
  
  Win = X[st:(st+120-1),]
  N = ncol(Win)
  Te = nrow(Win)
  
  set.seed(st)
  f.DAFM_low[[st]] = dafm_est(Win, r=7, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(2,2,1,1,1,1,1))$fmat
  
  print(st)
}

#saveRDS(list(f.DAFM=f.DAFM_low), "./Real data analysis/Fred-MD/DAFM_low.rds")

## DAFM_med
f.DAFM_med = list()
for(st in 1:(nrow(X)-120)){
  
  Win = X[st:(st+120-1),]
  N = ncol(Win)
  Te = nrow(Win)
  
  set.seed(st)
  f.DAFM_med[[st]] = dafm_est(Win, r=7, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(1,1,2,2,2,1,1))$fmat
  
  print(st)
}

#saveRDS(list(f.DAFM=f.DAFM_med), "./Real data analysis/Fred-MD/DAFM_med.rds")

## DAFM_high
f.DAFM_high = list()
for(st in 1:(nrow(X)-120)){
  
  Win = X[st:(st+120-1),]
  N = ncol(Win)
  Te = nrow(Win)
  
  set.seed(st)
  f.DAFM_high[[st]] = dafm_est(Win, r=7, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(1,1,1,1,1,2,2))$fmat
  
  print(st)
}

#saveRDS(list(f.DAFM=f.DAFM_high, "./Real data analysis/Fred-MD/DAFM_high.rds")
