library(cqrReg)
library(expm)
library(doParallel)

# Modified cqr.admm (for step 3) 
cqr.admm_mod = function (X.star, y.star, tau.vec, rho, maxit, toler) { 
  # X.star : Nk x p, y.star: Nk x 1, tau.vec : k x 1
  if (missing(maxit)) {
    maxit = 200
  }
  if (missing(toler)) {
    toler = 0.001
  }
  if (missing(rho)) {
    rho = 0.4 
  }
  k = length(tau.vec)
  N = nrow(X.star)/k
  p = ncol(X.star)
  
  tau.star = rep(tau.vec, each=N)
  beta = solve(t(X.star) %*% X.star, t(X.star) %*% y.star) # beta: p x 1
  
  betah = CQRADMMCPP(X.star, y.star, beta, toler, maxit, tau.star, rho, p)
  return(betah)
}

loss.fun <-function(tau, vec){
  loss = sapply(vec, function(x){ return(ifelse( x>0, tau*x, (tau-1)*x)) }) 
  return( mean(loss) )
}


# Main algorithm
dafm_est = function(X, r, tau.vec, tol, k.star=1, weight){ # X : Te*N 
  
  Te = nrow(X)
  N = ncol(X)
  K = length(tau.vec)
  maxiter = 200
  
  ##### Step 1. Initial value F^(0)
  fac = matrix(rnorm(Te*r), nrow=Te)
  
  cl <- makeCluster(detectCores())  # Use available cores
  clusterExport(cl, c("cqr.admm_mod","CQRADMMCPP")) 
  registerDoParallel(cl)
  
  preobj = 0
  for(iter in 1:maxiter){
    ##### Step 2. Update Lambda
    if(r==1){
      lamb <- lapply(tau.vec, function(tau) {
        as.matrix(apply(X, 2, function(x) rq(x ~ fac - 1, tau = tau)$coef),nrow=N) 
      })
    }
    if(r!=1){
    lamb <- lapply(tau.vec, function(tau) {
      t(apply(X, 2, function(x) rq(x ~ fac - 1, tau = tau)$coef)) #N x r
    })
    }
    fac <- foreach(t = 1:Te, .combine = rbind) %dopar% {
      X.st <- rep(X[t, ], K)*rep(weight,each=N)
      lamb.st <- diag(rep(weight,each=N))%*%do.call(rbind, lamb)
      t(cqr.admm_mod(X.star = lamb.st, y.star = X.st, tau.vec = tau.vec))
    }
    #### Step 4. Iterate steps 2,3 until convergence
    obj = 0
    for(t in 1:Te){
      for(k in 1:K){
        obj = obj + loss.fun(tau.vec[k], X[t,]-lamb[[k]]%*%fac[t,])
      }
    }
    obj = obj / K / Te
    if((preobj!=0) && (abs(obj-preobj) < tol)){
      print(c(iter,"converged",obj))
      break
    }
    print(c(iter, abs(obj-preobj), obj))
    preobj = obj
  }
  #### Step 5. Normalize F, L_1, ..., L_K
  dum = t(fac) %*% fac / Te
  M = sqrtm(dum) %*% (t(lamb[[k.star]]) %*% lamb[[k.star]]/N) %*% sqrtm(dum) #tau_{k.star}에 대해 normalize
  H = solve(sqrtm(dum), eigen(M)$vectors)
  
  fac.norm = fac %*% H
  H.inv = solve(H)
  lamb.norm <- lapply(lamb, function(l){ l %*% t(H.inv) })
  
  stopCluster(cl)
  
  return(list(fmat=fac.norm, lmat=lamb.norm))
}

