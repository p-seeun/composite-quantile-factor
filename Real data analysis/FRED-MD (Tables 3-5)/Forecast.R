source("./Real data analysis/Fred-MD/Forecast functions.R")
library(readr) 
library(tseries)

# Differentiate unemployment rate
dum <- read_csv("./Real data analysis/Fred-MD/2020-03.csv")
dum$sasdate[c(493,733)] #1999 Dec - 2019 Dec
Y = dum$UNRATE[494:733]-dum$UNRATE[493:732]

# Get estimated factors from FRED-MD
X = readRDS("./Real data analysis/Fred-MD/FRED-MD.rds")
QFMs = readRDS("./Real data analysis/Fred-MD/QFM Factors.rds")
Huang = readRDS("./Real data analysis/Fred-MD/Huang Factors.rds")
CQFM_unif = readRDS("./Real data analysis/Fred-MD/CQFM_unif.rds")
CQFM_low= readRDS("./Real data analysis/Fred-MD/CQFM_low.rds")
CQFM_med = readRDS("./Real data analysis/Fred-MD/CQFM_med.rds")
CQFM_high = readRDS("./Real data analysis/Fred-MD/CQFM_high.rds")

pmax = 5
for(h in 1:6){
  lags = dum$UNRATE[494:733] - dum$UNRATE[(494-h):(733-h)]
  Y.true = lags[(120+h):240] + dum$UNRATE[(733-120):(733-h)]
  
  # AR
  pred = c()
  for(st in 1:(120-h+1)){
    Z.win = lags[st:(st+119)]
    Y.win = Y[st:(st+119)]
    pred[st] = Frcst_ar(Z.win, Y.win, h, pmax)$pred
  }
  pred = pred + dum$UNRATE[(733-120):(733-h)] 
  ar = sqrt(mean((pred-Y.true)^2)) 
  
  # AR + PCA
  pred = c()
  for(st in 1:(120-h+1)){
    Z.win = lags[st:(st+119)]
    Y.win = Y[st:(st+119)]
    X.win = X[st:(st+119),]
    svd = svd(X.win); N = ncol(X.win)
    lamb = t(as.matrix(sqrt(N)*t(svd$v)[1:8,]))
    FPCA = (X.win %*% lamb)/N
    pred[st] = Frcst(Z.win,Y.win, h, pmax, Fact=FPCA)$pred
  }
  pred = pred + dum$UNRATE[(733-120):(733-h)]
  pca = sqrt(mean((pred-Y.true)^2)) 
  
  # AR + QFM(tau) / tau=0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99
  pred01 = c(); pred1 = c(); pred3 = c(); pred5 = c() 
  pred7 = c(); pred9 = c(); pred99 = c(); 
  
  for(st in 1:(120-h+1)){
    Z.win = lags[st:(st+119)]
    Y.win = Y[st:(st+119)]
    pred01[st] = Frcst(Z.win, Y.win, h, pmax, Fact=cbind(QFMs$f.QFM01[[st]]))$pred
    pred1[st] = Frcst(Z.win, Y.win, h, pmax, Fact=cbind(QFMs$f.QFM1[[st]]))$pred
    pred3[st] = Frcst(Z.win, Y.win, h, pmax, Fact=cbind(QFMs$f.QFM3[[st]]))$pred
    pred5[st] = Frcst(Z.win, Y.win, h, pmax, Fact=cbind(QFMs$f.QFM5[[st]]))$pred
    pred7[st] = Frcst(Z.win, Y.win, h, pmax, Fact=cbind(QFMs$f.QFM7[[st]]))$pred
    pred9[st] = Frcst(Z.win, Y.win, h, pmax, Fact=cbind(QFMs$f.QFM9[[st]]))$pred
    pred99[st] = Frcst(Z.win, Y.win, h, pmax, Fact=cbind(QFMs$f.QFM99[[st]]))$pred
  }
  pred01 = pred01 + dum$UNRATE[(733-120):(733-h)]
  pred1 = pred1 + dum$UNRATE[(733-120):(733-h)]
  pred3 = pred3 + dum$UNRATE[(733-120):(733-h)]
  pred5 = pred5 + dum$UNRATE[(733-120):(733-h)]
  pred7 = pred7 + dum$UNRATE[(733-120):(733-h)]
  pred9 = pred9 + dum$UNRATE[(733-120):(733-h)]
  pred99 = pred99 + dum$UNRATE[(733-120):(733-h)]
  
  qfm01 = sqrt(mean((pred01-Y.true)^2)) 
  qfm1 = sqrt(mean((pred1-Y.true)^2))
  qfm3 = sqrt(mean((pred3-Y.true)^2)) 
  qfm5 = sqrt(mean((pred5-Y.true)^2)) 
  qfm7 = sqrt(mean((pred7-Y.true)^2)) 
  qfm9 = sqrt(mean((pred9-Y.true)^2)) 
  qfm99 = sqrt(mean((pred99-Y.true)^2)) 
  
  pred_Huang = c()
  for(st in 1:(120-h+1)){
    Z.win = lags[st:(st+119)]
    Y.win = Y[st:(st+119)]
    pred_Huang[st] = Frcst(Z.win,Y.win, h, pmax, Fact=cbind(Huang$f.Huang[[st]]))$pred
    }
  pred_Huang = pred_Huang + dum$UNRATE[(733-120):(733-h)]
  huang = sqrt(mean((pred_Huang-Y.true)^2)) 
  
  # AR + CQFM_{weight} / weight = unif, low, med, high
  pred_unif = c(); pred_low = c(); pred_med = c(); pred_high = c();
  for(st in 1:(120-h+1)){
    Z.win = lags[st:(st+119)]
    Y.win = Y[st:(st+119)]
    pred_unif[st] = Frcst(Z.win,Y.win, h, pmax, Fact=cbind(CQFM_unif$f.CQFM[[st]]))$pred
    pred_low[st] = Frcst(Z.win,Y.win, h, pmax, Fact=cbind(CQFM_low$f.CQFM[[st]]))$pred
    pred_med[st] = Frcst(Z.win,Y.win, h, pmax, Fact=cbind(CQFM_med$f.CQFM[[st]]))$pred
    pred_high[st] = Frcst(Z.win,Y.win, h, pmax, Fact=cbind(CQFM_high$f.CQFM[[st]]))$pred
  }
  pred_unif = pred_unif + dum$UNRATE[(733-120):(733-h)]
  pred_low = pred_low + dum$UNRATE[(733-120):(733-h)]
  pred_med = pred_med + dum$UNRATE[(733-120):(733-h)]
  pred_high = pred_high + dum$UNRATE[(733-120):(733-h)]
  
  cqfm_unif = sqrt(mean((pred_unif-Y.true)^2)) 
  cqfm_low = sqrt(mean((pred_low-Y.true)^2)) 
  cqfm_med = sqrt(mean((pred_med-Y.true)^2)) 
  cqfm_high = sqrt(mean((pred_high-Y.true)^2)) 
  
  print(h)
  print(list(ar=(ar/ar)^2, pca=(pca/ar)^2, qfm01=(qfm01/ar)^2, qfm1=(qfm1/ar)^2, 
             qfm3=(qfm3/ar)^2, qfm5=(qfm5/ar)^2, qfm7=(qfm7/ar)^2, qfm9=(qfm9/ar)^2, qfm99=(qfm99/ar)^2, huang=(huang/ar)^2, 
             cqfm_unif=(cqfm_unif/ar)^2, cqfm_low=(cqfm_low/ar)^2, cqfm_med=(cqfm_med/ar)^2, cqfm_high=(cqfm_high/ar)^2))
}

