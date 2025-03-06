#####This file contains the main functions we used to implement our RMLC method#####
#####including functions for Algorithm 1, 2, and 3#####


#####Robust loss functions used in the paper#####

#Tukey biweight loss function
tukey_loss=function(x, M){
  y = x^2/2-x^4/(2*M^2)+x^6/(6*M^4)
  y = y*(abs(x)<M)+M^2/6*(abs(x)>=M)
  return(y)
}

#first derivative of the Tukey biweight loss function
tukey_psi=function(x, M){
  y = x*(1-(x/M)^2)^2;
  y = y*(abs(x)<M)+0*(abs(x)>=M)
  return(y)
}

#Huber's loss function
huber_loss=function(x, M){
  y = x^2/2
  y = y*(abs(x)<M)+(M*(abs(x)-0.5*M))*(abs(x)>=M)
  return(y)
}

#first derivative of the Huber's loss function
huber_psi=function(x, M){
  y = x
  y = y*(abs(x)<M)+M*sign(x)*(abs(x)>=M)
  return(y)
}




#####Algorithm 1: multivariate robust lasso with a given Omega#####

MVRlasso=function(Y, X, type=c("square","huber","tukey"), M, Omega, lambda, weight=1, tau1=1, max_iters = 100, w = 10,
              backtrack = TRUE, stepsizeShrink = 0.5, convergence_in = 1e-8){
  
  
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)
  Ymean=colMeans(Y)
  Xmean=colMeans(X)
  Y=scale(Y,scale=FALSE)
  X=scale(X,scale=FALSE)
  
  if(type=="square"){
  #objective functions squared loss
  f <- function(beta){sum(((X%*%beta-Y)%*%pracma::sqrtm(Omega)$B)^2)/n}
  gradf <- function(beta){2*t(X)%*%(X%*%beta-Y)%*%Omega/n}
  } else if(type=="huber"){
  #objective functions huber loss
  f <- function(beta){2*sum(huber_loss((X%*%beta-Y)%*%pracma::sqrtm(Omega)$B,M))/n}
  gradf <- function(beta){2*t(X)%*%huber_psi((X%*%beta-Y)%*%pracma::sqrtm(Omega)$B,M)%*%pracma::sqrtm(Omega)$B/n}
  } else {
  #objective functions tukey loss
  f <- function(beta){2*sum(tukey_loss((X%*%beta-Y)%*%pracma::sqrtm(Omega)$B,M))/n}
  gradf <- function(beta){2*t(X)%*%tukey_psi((X%*%beta-Y)%*%pracma::sqrtm(Omega)$B,M)%*%pracma::sqrtm(Omega)$B/n}
  }
  
  ###non-smooth part, lasso penalty, proximal soft thresholding
  g <- function(beta) {lambda*sum(abs(beta))}
  proxg <- function(beta, tau) { sign(beta)*(sapply(abs(beta) - tau*(lambda*weight),
                                                    FUN=function(x) {max(x,0)})) }
  
  beta0 <- mrce(Y=Y, X=X, lam2=lambda, method="fixed.omega", omega=Omega)$Bhat 
  if(max(beta0)> 1e+10){
    beta0 <- matrix(0,p,q)
  }
  
  sol=apgmatrix(f, gradf, g, proxg, beta0, tau1=tau1, max_iters = max_iters, w = w,
                backtrack = backtrack, stepsizeShrink = stepsizeShrink,
                convergence_in = convergence_in)
  
  mu=Ymean-t(sol$x)%*%Xmean
  
  return(list(beta =sol$x, mu=mu, obj=sol$objective, nit=length(sol$objective),taus=sol$taus))
}


#gradients of huber and tukey loss function

gradOmega_huber=function(E,Omega,M){
  S=pracma::sqrtm(Omega)$B
  q=dim(Omega)[1]
  n=dim(E)[1]
  A=solve(S%x%diag(1,q)+diag(1,q)%x%S)
  B=rep(NA,q*q)
  for (i in 1:(q*q)){
    C=matrix(A[i,],q,q,byrow=TRUE)
    B[i]=sum(diag(t(E)%*%huber_psi(E%*%S,M)%*%(C+t(C))/2))
  }
  return(matrix(B,q,q)/n)
}

gradOmega_tukey=function(E,Omega,M){
  S=pracma::sqrtm(Omega)$B
  q=dim(Omega)[1]
  n=dim(E)[1]
  A=solve(S%x%diag(1,q)+diag(1,q)%x%S)
  B=rep(NA,q*q)
  for (i in 1:(q*q)){
    C=matrix(A[i,],q,q,byrow=TRUE)
    B[i]=sum(diag(t(E)%*%tukey_psi(E%*%S,M)%*%(C+t(C))/2))
  }
  return(matrix(B,q,q)/n)
}


#####Algorithm 2: robust graphical lasso with a given Beta#####

MVRlasso.Omega=function(Y, X, type=c("square","huber","tukey"), M, Omega,beta, mu, lambda, diag=FALSE, weight=1, tau1=1, max_iters = 100, w = 10,
                  backtrack = TRUE, stepsizeShrink = 0.5, convergence_in = 1e-8){
  
  
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)
  
  E=Y-X%*%beta-rep(1,n)%*%t(mu)

  if(type=="square"){
    #objective functions squared loss
    f <- function(Omega){sum((E%*%pracma::sqrtm(Omega)$B)^2)/n-log(det(Omega))}
    gradf <- function(Omega){t(E)%*%E/(n)-solve(Omega)}
  } else if(type=="huber"){
    #objective functions huber loss
    f <- function(Omega){2*sum(huber_loss(E%*%pracma::sqrtm(Omega)$B,M))/n-log(det(Omega))}
    gradf <- function(Omega){2*gradOmega_huber(E,Omega,M)-solve(Omega)}
  } else {
    #objective functions tukey loss
    f <- function(Omega){2*sum(tukey_loss(E%*%pracma::sqrtm(Omega)$B,M))/n-log(det(Omega))}
    gradf <- function(Omega){2*gradOmega_tukey(E,Omega,M)-solve(Omega)}
  }
  
  ###non-smooth part, lasso penalty, proximal soft thresholding
  if(diag==TRUE){
  g <- function(Omega) {lambda*sum(abs(Omega))}
  proxg <- function(Omega, tau) { sign(Omega)*(sapply(abs(Omega) - tau*(lambda*weight),
                                                    FUN=function(x) {max(x,0)})) }
  }else if(diag==FALSE){
  g <- function(Omega) {Omegano=Omega-diag(diag(Omega))
    lambda*sum(abs(Omegano))}
  proxg <- function(Omega, tau) {Omegano=Omega-diag(diag(Omega))
    Omegano=sign(Omegano)*(sapply(abs(Omegano) - tau*(lambda*weight),
                                                        FUN=function(x) {max(x,0)}))
    diag(Omegano)=diag(Omega)
    return(Omegano)}
  }
  
  Omegaini <- glasso(t(E)%*%E/n,lambda,penalize.diagonal = FALSE)$wi
  
  sol=apgmatrix.Omega(f, gradf, g, proxg, Omegaini, tau1=tau1, max_iters = max_iters, w = w,
                backtrack = backtrack, stepsizeShrink = stepsizeShrink,
                convergence_in = convergence_in)
  
  return(list(Omega =sol$x, obj=sol$objective, nit=length(sol$objective),taus=sol$taus))
}




#####Algorithm 3: the faster approximation that performs the iteration once#####

# PF: Omega will be input
MVRLC=function(Y, X, type=c("square","huber","tukey"), M1, M2, Omega,lambda1,lambda2, diag=FALSE, weight=1, weighto=1, tau1=1, max_iters = 50, w = 10,
               backtrack = TRUE, stepsizeShrink = 0.5, convergence_in = 1e-6){
  
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)
  
  
#  betaupdatef=MVRlasso(Y, X, type=type, M=M2, Omega=diag(q), lambda=lambda2, weight=weight, tau1=tau1, max_iters = max_iters, w = w,
#                       backtrack = backtrack, stepsizeShrink = stepsizeShrink, convergence_in = convergence_in)
   betaupdatef=MVRlasso(Y, X, type=type, M=M2, Omega=Omega, lambda=lambda2, weight=weight, tau1=tau1, max_iters = max_iters, w = w,
                       backtrack = backtrack, stepsizeShrink = stepsizeShrink, convergence_in = convergence_in)

  omegaupdatef=MVRlasso.Omega(Y, X, type=type, M=M1, beta=betaupdatef$beta, mu=betaupdatef$mu, lambda=lambda1, diag=diag, weight=weighto, tau1=tau1, max_iters = max_iters, w = w,
                                backtrack = backtrack, stepsizeShrink = stepsizeShrink, convergence_in = convergence_in)

  betaupdatef=MVRlasso(Y, X, type=type, M=M2, Omega=omegaupdatef$Omega, lambda=lambda2, weight=weight, tau1=tau1, max_iters = max_iters, w = w,
                       backtrack = backtrack, stepsizeShrink = stepsizeShrink, convergence_in = convergence_in)
  
  
  return(list(beta =betaupdatef$beta, mu=betaupdatef$mu, Omega=omegaupdatef$Omega, taus2=betaupdatef$taus, taus1=omegaupdatef$taus))
}



#####Algorithm S.1: performs the full outer iteration#####

MVRLC.complete=function(Y, X, type=c("square","huber","tukey"), M1, M2, Omega,lambda1,lambda2, diag=FALSE, weight=1, weighto=1, tau1=1, max_iters = 50, max_iters_out=50, w = 10,
                        backtrack = TRUE, stepsizeShrink = 0.5, convergence_in = 1e-6, convergence_out=1e-3){
  
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)
  
  betaridge=solve(t(X)%*%X+lambda2*diag(p))%*%t(X)%*%Y
  betaridgenorm=sum(abs(betaridge))
  
  betai=MVRlasso(Y, X, type=type, M=M2, Omega=Omega, lambda=lambda2, weight=weight, tau1=tau1, max_iters = max_iters, w = w,
                 backtrack = backtrack, stepsizeShrink = stepsizeShrink, convergence_in = convergence_in)
  
  betaupdate=betai$beta
  muupdate=betai$mu
  
  betanorm=NULL
  
  for(i in 1:max_iters_out){
    
    omegaupdatef=MVRlasso.Omega(Y, X, type=type, M=M1, beta=betaupdate, mu=muupdate, lambda=lambda1, diag=diag, weight=weighto, tau1=tau1, max_iters = max_iters, w = w,
                                backtrack = backtrack, stepsizeShrink = stepsizeShrink, convergence_in = convergence_in)
    omegaupdate=omegaupdatef$Omega
    
    betaold=betaupdate
    
    betaupdatef=MVRlasso(Y, X, type=type, M=M2, Omega=omegaupdate, lambda=lambda2, weight=weight, tau1=tau1, max_iters = max_iters, w = w,
                         backtrack = backtrack, stepsizeShrink = stepsizeShrink, convergence_in = convergence_in)
    
    betaupdate=betaupdatef$beta
    muupdate=betaupdatef$mu
    
    betaupnorm=sum(abs(betaupdate-betaold))
    betanorm=c(betanorm,betaupnorm)
    
    if (betaupnorm < convergence_out*betaridgenorm) {
      break
    }
    
  }
  
  if(i==max_iters_out){
    warning("Outer loop reaches the maximun number of iteration")
  }
  
  return(list(beta =betaupdatef$beta, mu=betaupdatef$mu, Omega=omegaupdatef$Omega, nit=i, betanorm=betanorm, taus2=betaupdatef$taus, taus1=omegaupdatef$taus))
}



#####fucntion to select tuning parameters using TMSPE#####

MVRLC.cv=function(X, Y, K, trim=0.8, Omega,lambda1.min=0.01, lambda2.min=0.01, lambda1.max=1, lambda2.max=1, nlambda1=5, nlambda2=5, type=c("square","huber","tukey"), M1, M2, diag=FALSE, weight=1, weighto=1, tau1=1, max_iters = 50, max_iters_out=10, w = 10,
                  backtrack = TRUE, stepsizeShrink = 0.5, convergence_in = 1e-6, convergence_out=1e-3){
  q=ncol(Y)
  p=ncol(X)
  Data=cbind(Y,X)
  #Randomly shuffle the data
  sData<-Data[sample(nrow(Data)),]
  
  #Create K equally size folds
  folds <- cut(seq(1,nrow(sData)),breaks=K,labels=FALSE)
  testData=list()
  trainData=list()
  
  #Perform K fold cross validation
  for(i in 1:K){
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData[[i]] <- sData[testIndexes, ]
    trainData[[i]] <- sData[-testIndexes, ]
  }
  
  lambda1seq=seq(lambda1.min,lambda1.max,length.out=nlambda1)
  lambda2seq=seq(lambda2.min,lambda2.max,length.out=nlambda2)
  lambda=expand.grid(lambda1=lambda1seq,lambda2=lambda2seq)
  
  criteria=vector()
  sparsity=vector()

  for (k in 1:nrow(lambda)){
    lambda1=lambda[k,1]
    lambda2=lambda[k,2]
    
    cv.error=NULL
    cv.sparsity=NULL
    for (i in 1:K){
      train=trainData[[i]]
      test=testData[[i]]
      sol=MVRLC(train[,1:q], train[,-(1:q)], type=type, Omega,M1=M1, M2=M2, lambda1,lambda2, diag=diag, weight=weight, weighto=weighto, tau1=tau1, max_iters = max_iters, w = w,
                backtrack = backtrack, stepsizeShrink = stepsizeShrink, convergence_in = convergence_in)
      
      spars <- sum(sol$beta!=0)
      if (spars==0){warning(paste("Solution for lambda1=",lambda1,"and lambda2=",
                                             lambda2,"fully sparse!"))}

            ###calculate prediction error
      error=test[,1:q]-test[,-(1:q)]%*%sol$beta-rep(1,nrow(test))%*%t(sol$mu)
      errorsq=sort(error^2)[1:as.integer(trim*(length(error^2)))]
      if (spars==0){errorsq <- NA} # don't allow for full sparsity
      cv.error=c(cv.error, sum(errorsq))
      cv.sparsity=c(cv.sparsity,spars)
    }  
    criteria[k]=mean(cv.error)
    sparsity[k]=mean(cv.sparsity)
    #plot(sparsity,criteria)
  }
  mink=which.min(criteria)
  lambda1=lambda[mink,1]
  lambda2=lambda[mink,2]
  
  return(list(error.CV =criteria, lambda1=lambda1, lambda2=lambda2))
}


#####main function to perform APG algorithm for muti response#####
apgmatrix=function (f, gradf, g, proxg, x0, tau1, max_iters = 100, w = 10, 
                    backtrack = TRUE, stepsizeShrink = 0.5, convergence_in = 1e-8) 
{
  residual <- double(max_iters)
  normalizedResid <- double(max_iters)
  taus <- double(max_iters)
  fVals <- double(max_iters)
  objective <- double(max_iters + 1)
  totalBacktracks <- 0
  backtrackCount <- 0
  x1 <- x0
  xprev <- x1
  d1 <- x1
  f1 <- f(d1)
  fVals[1] <- f1
  gradf1 <- gradf(d1)
  
  maxResidual <- -Inf
  minObjectiveValue <- Inf
  objective[1] <- f1 + g(x0)
  for (i in 1:max_iters) {
    x0 <- x1 + (i/(i+3))*(x1 - xprev)
    xprev <- x1
    gradf0 <- matrix(gradf1)
    tau0 <- tau1
    x1hat <- x0 - tau0 * c(gradf0)
    x1 <- proxg(x1hat, tau0)
    Dx <- matrix(x1 - x0)
    d1 <- x1
    f1 <- f(d1)
    if (backtrack) {
      M <- max(fVals[max(i - w, 1):max(i - 1, 1)])
      backtrackCount <- 0
      prop <- (f1 - 1e-12 > M + t(Dx) %*% gradf0 + 0.5 * 
                 (norm(Dx, "f")^2)/tau0) && (backtrackCount < 
                                               20)
      while (prop) {
        tau0 <- tau0 * stepsizeShrink
        x1hat <- x0 - tau0 * c(gradf0)
        x1 <- proxg(x1hat, tau0)
        d1 <- x1
        f1 <- f(d1)
        Dx <- matrix(x1 - x0)
        backtrackCount <- backtrackCount + 1
        prop <- (f1 - 1e-12 > M + t(Dx) %*% gradf0 + 
                   0.5 * (norm(Dx, "f")^2)/tau0) && (backtrackCount < 
                                                       20)
      }
      totalBacktracks <- totalBacktracks + backtrackCount
    }
    taus[i] <- tau0
    residual[i] <- norm(Dx, "f")/tau0
    maxResidual <- max(maxResidual, residual[i])
    normalizer <- max(norm(gradf0, "f"), norm(as.matrix(x1 - 
                                                          x1hat), "f")/tau0) + convergence_in
    normalizedResid[i] <- residual[i]/normalizer
    fVals[i] <- f1
    objective[i + 1] <- f1 + g(x1)
    newObjectiveValue <- objective[i + 1]
    
    if (newObjectiveValue < minObjectiveValue) {
      bestObjectiveIterate <- x1
      minObjectiveValue <- min(minObjectiveValue, newObjectiveValue)
    }
    gradf1 <- gradf(d1)
    Dg <- matrix(gradf1 + (x1hat - x0)/tau0)
    dotprod <- t(Dx) %*% Dg
    tau_s <- as.vector(norm(Dx, "f")^2/dotprod)
    tau_m <- as.vector(dotprod/norm(Dg, "f")^2)
    tau_m <- max(tau_m, 0)
    if (abs(dotprod) < convergence_in) 
      break
    if (2 * tau_m > tau_s) {
      tau1 <- tau_m
    }
    else {
      tau1 <- tau_s - 0.5 * tau_m
    }
    if ((tau1 <= 0) || is.infinite(tau1) || is.nan(tau1)) {
      tau1 <- tau0 * 1.5
    }
  }
  if(i==max_iters){
    warning("reach the maximun number of iteration")
  }
  return(list(x = bestObjectiveIterate, objective = objective[1:(i + 
                                                                   1)], fVals = fVals[1:i], totalBacktracks = totalBacktracks, 
              residual = residual[1:i], taus = taus[1:i]))
}



#####main function to perform apg for solving Omega#####
apgmatrix.Omega=function (f, gradf, g, proxg, x0, tau1, max_iters = 100, w = 10, 
                          backtrack = TRUE, stepsizeShrink = 0.5, convergence_in = 1e-8) 
{
  residual <- double(max_iters)
  normalizedResid <- double(max_iters)
  taus <- double(max_iters)
  fVals <- double(max_iters)
  objective <- double(max_iters + 1)
  totalBacktracks <- 0
  backtrackCount <- 0
  x1 <- x0
  d1 <- x1
  f1 <- f(d1)
  fVals[1] <- f1
  gradf1 <- gradf(d1)
  
  maxResidual <- -Inf
  minObjectiveValue <- Inf
  objective[1] <- f1 + g(x0)
  for (i in 1:max_iters) {
    x0 <- x1 
    gradf0 <- matrix(gradf1)
    tau0 <- tau1
    x1hat <- x0 - tau0 * c(gradf0)
    x1 <- proxg(x1hat, tau0)
    while (det(x1)<=0){
      if(i==1){
        warning("initial tau is too large")
      }
      tau0=tau0 * stepsizeShrink
      x1hat <- x0 - tau0 * c(gradf0)
      x1 <- proxg(x1hat, tau0)
    }
    Dx <- matrix(x1 - x0)
    d1 <- x1
    f1 <- f(d1)
    if (backtrack) {
      M0 <- max(fVals[max(i - w, 1):max(i - 1, 1)])
      backtrackCount <- 0
      prop <- (f1 - 1e-12 > M0 + t(Dx) %*% gradf0 + 0.5 * 
                 (norm(Dx, "f")^2)/tau0) && (backtrackCount < 
                                               20)
      while (prop) {
        tau0 <- tau0 * stepsizeShrink
        x1hat <- x0 - tau0 * c(gradf0)
        x1 <- proxg(x1hat, tau0)
        d1 <- x1
        while(det(d1)<=0){
          tau0=tau0 * stepsizeShrink
          x1hat <- x0 - tau0 * c(gradf0)
          x1 <- proxg(x1hat, tau0)
          d1 <- x1
        }
        f1 <- f(d1)
        Dx <- matrix(x1 - x0)
        backtrackCount <- backtrackCount + 1
        prop <- (f1 - 1e-12 > M0 + t(Dx) %*% gradf0 + 
                   0.5 * (norm(Dx, "f")^2)/tau0) && (backtrackCount < 
                                                       20)
      }
      totalBacktracks <- totalBacktracks + backtrackCount
    }
    taus[i] <- tau0
    residual[i] <- norm(Dx, "f")/tau0
    maxResidual <- max(maxResidual, residual[i])
    normalizer <- max(norm(gradf0, "f"), norm(as.matrix(x1 - 
                                                          x1hat), "f")/tau0) + convergence_in
    normalizedResid[i] <- residual[i]/normalizer
    fVals[i] <- f1
    objective[i + 1] <- f1 + g(x1)
    newObjectiveValue <- objective[i + 1]
    
    if (newObjectiveValue < minObjectiveValue) {
      bestObjectiveIterate <- x1
      minObjectiveValue <- min(minObjectiveValue, newObjectiveValue)
    }
    gradf1 <- gradf(d1)
    Dg <- matrix(gradf1 + (x1hat - x0)/tau0)
    dotprod <- t(Dx) %*% Dg
    tau_s <- as.vector(norm(Dx, "f")^2/dotprod)
    tau_m <- as.vector(dotprod/norm(Dg, "f")^2)
    tau_m <- max(tau_m, 0)
    if (abs(dotprod) < convergence_in) 
      break
    if (2 * tau_m > tau_s) {
      tau1 <- tau_m
    } else {
      tau1 <- tau_s - 0.5 * tau_m
    } 
    if ((tau1 <= 0) || is.infinite(tau1) || is.nan(tau1)) {
      tau1 <- tau0 * 1.5
    }
  }
  if(i==max_iters){
    warning("reach the maximun number of iteration")
  }
  return(list(x = bestObjectiveIterate, objective = objective[1:(i + 
                                                                   1)], fVals = fVals[1:i], totalBacktracks = totalBacktracks, 
              residual = residual[1:i], taus = taus[1:i]))
}


#####multivariate robust MM-ridge#####

MVRLC.ridge=function(Y, X, type=c("square","huber","tukey"), Omega,M1, M2, lambda1=0,lambda2=0,lambdaridge=0.1, diag=FALSE, weight=1, weighto=1, cons1=5, cons2=1, tau1=1, max_iters = 50, w = 10,
                     backtrack = TRUE, stepsizeShrink = 0.5, convergence_in = 1e-6){
  
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)
  
  betaupdatef=MVRridge(Y, X, type=type, M=M2, Omega=Omega, lambda=lambda2, lambdaridge=lambdaridge, weight=weight, tau1=tau1, max_iters = max_iters, w = w,
                       backtrack = backtrack, stepsizeShrink = stepsizeShrink, convergence_in = convergence_in)
  
  omegaupdatef=suppressWarnings(MVRlasso.Omega(Y, X, type=type, M=M1, beta=betaupdatef$beta, mu=betaupdatef$mu, lambda=lambda1, diag=diag, weight=weighto, tau1=tau1, max_iters = max_iters, w = w,
                                               backtrack = backtrack, stepsizeShrink = stepsizeShrink, convergence_in = convergence_in))
  
  betaupdatef=MVRridge(Y, X, type=type, M=M2, Omega=omegaupdatef$Omega, lambda=lambda2, lambdaridge=lambdaridge, weight=weight, tau1=tau1, max_iters = max_iters, w = w,
                       backtrack = backtrack, stepsizeShrink = stepsizeShrink, convergence_in = convergence_in)
  
  
  return(list(Omega=omegaupdatef$Omega*cons1, beta =betaupdatef$beta*cons2, taus1=omegaupdatef$taus, taus2=betaupdatef$taus))
}


MVRridge=function(Y, X, type=c("square","huber","tukey"), M, Omega, lambda=0, lambdaridge, weight=1, tau1=1, max_iters = 100, w = 10,
                  backtrack = TRUE, stepsizeShrink = 0.5, convergence_in = 1e-8){
  
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)
  Ymean=colMeans(Y)
  Xmean=colMeans(X)
  Y=scale(Y,scale=FALSE)
  X=scale(X,scale=FALSE)
  
  if(type=="square"){
    #objective functions squared loss
    f <- function(beta){sum(((X%*%beta-Y)%*%pracma::sqrtm(Omega)$B)^2)/n + lambdaridge*sum(beta^2)}
    gradf <- function(beta){2*t(X)%*%(X%*%beta-Y)%*%Omega/n+2*lambdaridge*beta}
  } else if(type=="huber"){
    #objective functions huber loss
    f <- function(beta){2*sum(huber_loss((X%*%beta-Y)%*%pracma::sqrtm(Omega)$B,M))/n + lambdaridge*sum(beta^2)}
    gradf <- function(beta){2*t(X)%*%huber_psi((X%*%beta-Y)%*%pracma::sqrtm(Omega)$B,M)%*%pracma::sqrtm(Omega)$B/n+2*lambdaridge*beta}
  } else {
    #objective functions tukey loss
    f <- function(beta){2*sum(tukey_loss((X%*%beta-Y)%*%pracma::sqrtm(Omega)$B,M))/n + lambdaridge*sum(beta^2)}
    gradf <- function(beta){2*t(X)%*%tukey_psi((X%*%beta-Y)%*%pracma::sqrtm(Omega)$B,M)%*%pracma::sqrtm(Omega)$B/n+2*lambdaridge*beta}
  }
  
  ###non-smooth part, lasso penalty, proximal soft thresholding
  g <- function(beta) {lambda*sum(abs(beta))}
  proxg <- function(beta, tau) { sign(beta)*(sapply(abs(beta) - tau*(lambda*weight),
                                                    FUN=function(x) {max(x,0)})) }
  
  beta0 <- mrce(Y=Y, X=X, lam2=lambda, method="fixed.omega", omega=Omega)$Bhat 
  if(max(beta0)> 1e+10){
    beta0 <- matrix(0,p,q)
  }
  
  sol=apgmatrix(f, gradf, g, proxg, beta0, tau1=tau1, max_iters = max_iters, w = w,
                backtrack = backtrack, stepsizeShrink = stepsizeShrink,
                convergence_in = convergence_in)
  
  mu=Ymean-t(sol$x)%*%Xmean
  
  return(list(beta =sol$x, mu=mu, obj=sol$objective, nit=length(sol$objective),taus=sol$taus))
}

