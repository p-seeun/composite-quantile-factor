source("DAFM algorithm.R")
source("RobPCA, QFM algorithms.R")

X = readRDS("./Real data analysis/Fred-MD/FRED-MD.rds")

f.QFM01 = list()
f.QFM1 = list()
f.QFM3 = list()
f.QFM5 = list()
f.QFM7 = list()
f.QFM9 = list()
f.QFM99 = list()

# st: a start point of rolling window with size 120
for(st in 1:(nrow(X)-120)){
  Win = X[st:(st+120-1),]
  N = ncol(Win)
  Te = nrow(Win)
  
  set.seed(st)
  f.QFM01[[st]] = qfm_est(Win, r=1, tol=1e-5, tau=0.01)$fmat
  
  set.seed(st)
  f.QFM1[[st]] = qfm_est(Win, r=2, tol=1e-5, tau=0.1)$fmat
  
  set.seed(st)
  f.QFM3[[st]] = qfm_est(Win, r=6, tol=1e-5, tau=0.3)$fmat
  
  set.seed(st)
  f.QFM5[[st]] = qfm_est(Win, r=6, tol=1e-5, tau=0.5)$fmat
  
  set.seed(st)
  f.QFM7[[st]] = qfm_est(Win, r=7, tol=1e-5, tau=0.7)$fmat
  
  set.seed(st)
  f.QFM9[[st]] = qfm_est(Win, r=2, tol=1e-5, tau=0.9)$fmat

  set.seed(st)
  f.QFM99[[st]] = qfm_est(Win, r=1, tol=1e-5, tau=0.99)$fmat
  
  print(st)
}

#saveRDS(list(f.QFM01=f.QFM01, f.QFM1=f.QFM1, f.QFM3=f.QFM3, f.QFM5=f.QFM5, f.QFM7=f.QFM7, f.QFM9=f.QFM9, f.QFM99=f.QFM99),
#        "./Real data analysis/Fred-MD/QFM Factors.rds")
