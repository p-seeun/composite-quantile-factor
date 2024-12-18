library(cqrfactor)

X = readRDS("./Real data analysis/Fred-MD/FRED-MD.rds")

## Huang
f.Huang = list()
for(st in 1:(nrow(X)-120)){

  Win = X[st:(st+120-1),]
  N = ncol(Win)
  Te = nrow(Win)
  
  set.seed(st)
  f.Huang[[st]] = cqrfactor(Win, r=6, tau=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99))$fmat
  
  print(st)
}

saveRDS(list(f.Huang=f.Huang), "./Real data analysis/Fred-MD/Huang Factors.rds")

