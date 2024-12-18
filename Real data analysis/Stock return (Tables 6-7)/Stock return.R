source("CQFM algorithm.R")
source("RobPCA, QFM algorithms.R")
library(cqrfactor)

## 1. Get monthly stock return data of firms with share code 12 and 14.
library(tidyverse)
data = read.csv("./Real data analysis/Monthly stock return/CRSP rawdata(monthly).csv")

data2 = data[data$SHRCD %in% c(12,14),]
tib = as_tibble(data2[,c("date","RET","PERMNO")])
tib = tib %>% spread(key = PERMNO, value = RET)
X = as.matrix(tib)[,-1]

Y= matrix(nrow=nrow(X),ncol=ncol(X))
for(i in 1:nrow(X)){
  for(j in 1:ncol(X)){
    Y[i,j]=as.numeric(X[i,j])
  }
}
colnames(Y) <- colnames(X)

# Remove firms with NA
na.index = which(apply(Y, 2, function(x){sum(is.na(x))}!=0))
Z = Y[,-na.index]
colnames(Z) <- colnames(Y)[-na.index]
Z = scale(Z, center=T, scale=F)

dim(Z) # 240(months) x 233(firms)

## 2. Estimate factor number
N = ncol(Z)
Te = nrow(Z)

set.seed(1)
# CQFM_unif
CQFM_unif = cqfm_est(Z, r=15, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(1,1,1,1,1,1,1))
lmat = CQFM_unif$lmat
rhats = lapply(lmat, function(L){ x= eigen(t(L)%*%L/N)$values; return(sum( x > min(N,Te)^(-1/3)*x[1] ))}) 
rhat = max(unlist(rhats))
rhat # 2

#QFM01
QFM01 = qfm_est(Z, 15, tol=1e-5, tau=0.01)
lmat = QFM01$lmat
dum = t(lmat) %*% lmat / N
rhat <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1]) 
rhat # =1

#QFM1
QFM1 = qfm_est(Z, 15, tol=1e-5, tau=0.1)
lmat = QFM1$lmat
dum = t(lmat) %*% lmat / N
rhat <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1])
rhat # =1

#QFM3
QFM3 = qfm_est(Z, 15, tol=1e-5, tau=0.3)
lmat = QFM3$lmat
dum = t(lmat) %*% lmat / N
rhat <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1]) 
rhat # =2

#QFM5
QFM5 = qfm_est(Z, 15, tol=1e-5, tau=0.5)
lmat = QFM5$lmat
dum = t(lmat) %*% lmat / N
rhat <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1])
rhat # =2

#QFM7
QFM7 = qfm_est(Z, 15, tol=1e-5, tau=0.7)
lmat = QFM7$lmat
dum = t(lmat) %*% lmat / N
rhat <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1]) #log(N*Te/(N+Te))*(N+Te)/N/Te 
rhat # =2

#QFM
QFM9 = qfm_est(Z, 15, tol=1e-5, tau=0.9)
lmat = QFM9$lmat
dum = t(lmat) %*% lmat / N
rhat <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1]) #log(N*Te/(N+Te))*(N+Te)/N/Te 
rhat # =1

#QFM99
set.seed(1)
QFM99 = qfm_est(Z, 15, tol=1e-5, tau=0.99)
lmat = QFM99$lmat
dum = t(lmat) %*% lmat / N
rhat <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1]) #log(N*Te/(N+Te))*(N+Te)/N/Te 
rhat # =1

#CQFM_low
CQFM_low = cqfm_est(Z, r=15, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(2,2,1,1,1,1,1))
lmat = CQFM_low$lmat
rhats = lapply(lmat, function(L){ x= eigen(t(L)%*%L/N)$values; return(sum( x > min(N,Te)^(-1/3)*x[1] ))}) 
rhat = max(unlist(rhats))
rhat # 2

#CQFM_med
CQFM_med = cqfm_est(Z, r=15, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(1,1,2,2,2,1,1))
lmat = CQFM_med$lmat
rhats = lapply(lmat, function(L){ x= eigen(t(L)%*%L/N)$values; return(sum( x > min(N,Te)^(-1/3)*x[1] ))}) 
rhat = max(unlist(rhats))
rhat # 2

#CQFM_high
CQFM_high = cqfm_est(Z, r=15, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), tol=1e-5, weight=c(1,1,1,1,1,2,2))
lmat = CQFM_high$lmat
rhats = lapply(lmat, function(L){ x= eigen(t(L)%*%L/N)$values; return(sum( x > min(N,Te)^(-1/3)*x[1] ))}) 
rhat = max(unlist(rhats))
rhat # 2

#Mean FM
icp<-function(r){
  svd = svd(Z) 
  lambda = t(as.matrix(sqrt(N)*t(svd$v)[1:r,])) 
  Factor = (Z %*% lambda)/N 
  IC = log(sum((Z-Factor%*%t(lambda))^2)/(N*Te))+r*(N+Te)/(N*Te)*log(min(N,Te)) 
  return(IC)
}

icp.results = c()
for(r in 2:10){
  icp.results[r] = icp(r)
}
rhat = which.min(icp.results)
rhat # 8


## 3. Compute CIV
DATA = read.csv("./Real data analysis/Monthly stock return/CRSP rawdata(daily).csv")

DATA2 = DATA[ DATA$SHRCD %in% c(12,14) ,] 
TIB = as_tibble( DATA2[ , c("date","RET","PERMNO")] )
TIB = TIB %>% spread( key = PERMNO, value = RET )
TIB = type.convert( TIB, as.is =TRUE )
TIB$YM = format(as.Date(TIB$date), '%Y%m')
Months = unique(TIB$YM) #2004-01 - 2023-12

df1 = as.data.frame(TIB)
df2 = df1[,colnames(df1) %in% c(colnames(Z),"YM","date")] # Match firms to that of monthly data Z
df3 = df2[complete.cases(df2),] # 233 firms matched, Remove rows with NA
df4 = subset(df3, select=-c(YM,date))

dim(df4) # 5032(days) x 233(firms)

# Get Residuals
Resid = matrix(nrow = nrow(Z), ncol=ncol(Z))
for(t in 1:length(Months)){
  index = (df3$YM == Months[t])
  m.data = as.matrix( df4[index, ] )
  m.data.sc = scale(m.data, scale=F)
  svd = svd(m.data.sc)
  N = ncol(m.data); Te = ncol(m.data)
  lamb = t(as.matrix(sqrt(N)*t(svd$v)[1:8,])) 
  f.mean = (m.data.sc %*% lamb)/N
  for(i in 1:N){
    Resid[t,i] = sd( lm(m.data[,i]~f.mean)$residual )
  }
  print(t)
}

# Get CIVs
CIV = apply(Resid,1,mean)

# Get Factors
set.seed(1)
f.CQFM_unif = cqfm_est(X=Z, r=2, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), weight=c(1,1,1,1,1,1,1), tol=1e-5)$fmat
f.CQFM_low = cqfm_est(X=Z, r=2, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), weight=c(2,2,1,1,1,1,1), tol=1e-5)$fmat
f.CQFM_med = cqfm_est(X=Z, r=2, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), weight=c(1,1,2,2,2,1,1), tol=1e-5)$fmat
f.CQFM_high = cqfm_est(X=Z, r=2, tau.vec=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99), weight=c(1,1,1,1,1,2,2), tol=1e-5)$fmat

f.QFM01 = qfm_est(X=Z, r=1, tau=0.01, tol=1e-5)$fmat
f.QFM1 = qfm_est(X=Z, r=1, tau=0.1, tol=1e-5)$fmat
f.QFM3 = qfm_est(X=Z, r=1, tau=0.3, tol=1e-5)$fmat
f.QFM5 = qfm_est(X=Z, r=3, tau=0.5, tol=1e-5)$fmat
f.QFM7 = qfm_est(X=Z, r=2, tau=0.7, tol=1e-5)$fmat
f.QFM9 = qfm_est(X=Z, r=1, tau=0.9, tol=1e-5)$fmat
f.QFM99 = qfm_est(X=Z, r=1, tau=0.99, tol=1e-5)$fmat

f.Huang = cqrfactor(Z, r=2, tau=c(0.01,0.1,0.3,0.5,0.7,0.9,0.99))$fmat

svd = svd(Z)
N = ncol(Z)
lamb = t(as.matrix(sqrt(N)*t(svd$v)[1:8,]))
f.Mean = (Z %*% lamb)/N

# R^2 values from CIV factor (Total)
summary(lm(CIV~f.CQFM_unif))$adj.r.squared
summary(lm(CIV~f.QFM01))$adj.r.squared
summary(lm(CIV~f.QFM1))$adj.r.squared
summary(lm(CIV~f.QFM3))$adj.r.squared
summary(lm(CIV~f.QFM5))$adj.r.squared
summary(lm(CIV~f.QFM7))$adj.r.squared
summary(lm(CIV~f.QFM9))$adj.r.squared
summary(lm(CIV~f.QFM99))$adj.r.squared
summary(lm(CIV~f.Mean))$adj.r.squared

summary(lm(CIV~f.CQFM_low))$adj.r.squared
summary(lm(CIV~f.CQFM_med))$adj.r.squared
summary(lm(CIV~f.CQFM_high))$adj.r.squared

summary(lm(CIV~f.Huang))$adj.r.squared

# Get index for each share code
new = unique( DATA2 %>% select(PERMNO, SHRCD) %>% unite(new, PERMNO, SHRCD) )$new
NumCd.Tot = matrix(as.numeric(unlist(str_split(new,"_"))), ncol=2, byrow=T)
NumCd = NumCd.Tot[ NumCd.Tot[,1] %in% colnames(Z),]
Cmp12 = NumCd[NumCd[,2]==12,1]
Cmp14 = NumCd[NumCd[,2]==14,1]
index12 = colnames(Z) %in% Cmp12 
index14 = colnames(Z) %in% Cmp14 

# R^2 values from CIV factor (Share code 12)
CIV12 = apply(Resid[,index12],1,mean)
summary(lm(CIV12 ~ f.CQFM_unif))$adj.r.squared

summary(lm(CIV12 ~ f.QFM01))$adj.r.squared
summary(lm(CIV12 ~ f.QFM1))$adj.r.squared
summary(lm(CIV12 ~ f.QFM3))$adj.r.squared
summary(lm(CIV12 ~ f.QFM5))$adj.r.squared
summary(lm(CIV12 ~ f.QFM7))$adj.r.squared
summary(lm(CIV12 ~ f.QFM9))$adj.r.squared
summary(lm(CIV12 ~ f.QFM99))$adj.r.squared

summary(lm(CIV12 ~ f.Mean))$adj.r.squared

summary(lm(CIV12 ~ f.CQFM_low))
summary(lm(CIV12 ~ f.CQFM_med))
summary(lm(CIV12 ~ f.CQFM_high))

summary(lm(CIV12~f.Huang))$adj.r.squared

# R^2 values from CIV factor (Share code 14)
CIV14 = apply(Resid[,index14],1,mean)
summary(lm(CIV14 ~ f.CQFM_unif))$adj.r.squared

summary(lm(CIV14 ~ f.QFM01))$adj.r.squared
summary(lm(CIV14 ~ f.QFM1))$adj.r.squared
summary(lm(CIV14 ~ f.QFM3))$adj.r.squared
summary(lm(CIV14 ~ f.QFM5))$adj.r.squared
summary(lm(CIV14 ~ f.QFM7))$adj.r.squared
summary(lm(CIV14 ~ f.QFM9))$adj.r.squared
summary(lm(CIV14 ~ f.QFM99))$adj.r.squared

summary(lm(CIV14 ~ f.Mean))$adj.r.squared

summary(lm(CIV14 ~ f.CQFM_low))$adj.r.squared
summary(lm(CIV14 ~ f.CQFM_med))$adj.r.squared
summary(lm(CIV14 ~ f.CQFM_high))$adj.r.squared

summary(lm(CIV14~f.Huang))$adj.r.squared
