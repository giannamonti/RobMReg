# library(pracma) # sqrtm() needs to be taken from here!!!
# library(MethylCapSig)
# library(MRCKF)
# source("MRCD.R") # from https://wis.kuleuven.be/statdatascience/robust/software
# # robust Lasso
# source("APG_Peter.R")
# library(MRCE)
# library(robustHD) # for sparseLTS
# library(mvoutlier)

### Knockoff Functions
computeLCDnew = function(B){
  # B is a matrix here with K columns!
  p  <- nrow(B)
  orig <- 1:(p/2)
  kf <- (p/2+1):p
  oB <- abs(B[orig,])
  kB <- abs(B[kf,])
  W <- oB-kB
  return(W)
}

evalsim <- function(listres,listtrue){
  nsim <- length(listtrue)
  fdr <- tpr <- rep(NA,nsim)
  for (i in 1:nsim){
    fdr[i] <- sum((listres[[i]]%in%listtrue[[i]])==FALSE)/length(listres[[i]])
    tpr[i] <- sum(listres[[i]]%in%listtrue[[i]])/length(listtrue[[i]])
  }
  list(fdr=fdr,tpr=tpr)
}

evalsim2 <- function(listres,listtrue){
  # for the e-values
  nsim <- length(listtrue)
  fdr <- tpr <- rep(NA,nsim)
  for (i in 1:nsim){
    fdr[i] <- sum((listres[[i]]%in%listtrue[[i]])==FALSE)/length(listres[[i]]) # FDR
    tpr[i] <- sum(listres[[i]]%in%listtrue[[i]])/length(listtrue[[i]]) # TPR
  }
  list(fdr=fdr,tpr=tpr)
}

evalues <- function(X,Y,fit.cv,typeset,Omega,M=100,q=0.2){
  # Generate knockoff matrix:
  p <- ncol(X)
  K <- ncol(Y)
  E <- array(NA,dim=c(M,p,K))
  all_rej <- c()
  for(m in 1:M){
    if(m %% 20 == 1){
      cat(sprintf("Implementing the %d-th run.\n", m))
    }
    if (typeset=="square"){
      Xkf <- create.second_order_PF(X,robust=FALSE,shrink=TRUE)
      Xaug <- cbind(X, Xkf)
      XX <- as.matrix(scale(Xaug))
    }
    else{
      Xkf <- create.second_order_PF(X,robust=TRUE,shrink=TRUE)
      Xaug <- cbind(X, Xkf)
      XX <- as.matrix(scale(Xaug,apply(Xaug,2,median),apply(Xaug,2,mad)))
    }
    fit <- MVRLC(Y, XX, type=typeset, M1=4.685, M2=4.685,Omega=Omega, 
                 lambda1=fit.cv$lambda1,lambda2=fit.cv$lambda2, 
                 weight=1, weighto=1,max_iters = 200)
    Wstat <- computeLCDnew(fit$beta)
    #tau <- knockoff.threshold(as.vector(Wstat), fdr = q/2 , offset = 0) ## GM PF
    #E[m,,] <- (as.vector(Wstat) >= tau) / (1 + sum(as.vector(Wstat) <= -tau))
    for (j in 1:K){
      tau <- knockoff.threshold(as.vector(Wstat[,j]), fdr = q/2, offset = 0) ## GM PF
      E[m,,j] <- (as.vector(Wstat[,j]) >= tau) / (1 + sum(as.vector(Wstat[,j]) <= -tau))
    }
  }
  ## run e-BH
  Em <- apply(E,c(2,3),mean)
  E_ord <- order(as.vector(Em), decreasing = TRUE)
  Es <- sort(Em, decreasing = TRUE)
  comp <- Es >= (1 / q / (1:(p*K)))
  plot(Es,col=comp+2)
  points((1 / q / (1:(p*K))),col=7,pch=3)
  #comp <- E >= ((p*K) / q / (1:(p*K)))
  id <- max(which(comp > 0))
  if(id > 0){
    rej <- E_ord[1:id]
  }else{
    rej <- NULL
  }
  return(sort(rej))
}
  

evalues.filter <- function(X,Y,fit.cv,var.screen,typeset,Omega,M=100,q=0.2){
  # var.screen are the screen variables to be used for X
  # Generate knockoff matrix:
  X <- X[,var.screen]
  p <- ncol(X)
  K <- ncol(Y)
  E <- array(NA,dim=c(M,p,K))
  all_rej <- c()
  for(m in 1:M){
    if(m %% 20 == 1){
      cat(sprintf("Implementing the %d-th run.\n", m))
    }
    if (typeset=="square"){
      Xkf <- create.second_order_PF(X,robust=FALSE,shrink=TRUE)
      Xaug <- cbind(X, Xkf)
      XX <- as.matrix(scale(Xaug))
    }
    else{
      Xkf <- create.second_order_PF(X,robust=TRUE,shrink=TRUE)
      Xaug <- cbind(X, Xkf)
      XX <- as.matrix(scale(Xaug,apply(Xaug,2,median),apply(Xaug,2,mad)))
    }
    fit <- MVRLC(Y, XX, type=typeset, M1=4.685, M2=4.685,Omega=Omega, 
                 lambda1=fit.cv$lambda1,lambda2=fit.cv$lambda2, 
                 weight=1, weighto=1,max_iters = 200)
    Wstat <- computeLCDnew(fit$beta)
    #tau <- knockoff.threshold(as.vector(Wstat), fdr = q/2 , offset = 0) ## GM PF
    #E[m,,] <- (as.vector(Wstat) >= tau) / (1 + sum(as.vector(Wstat) <= -tau))
    for (j in 1:K){
      #tau <- knockoff.threshold(as.vector(Wstat[,j]), fdr = q, offset = 0) ## GM PF
      tau <- knockoff.threshold(as.vector(Wstat[,j]), fdr = q/2, offset = 0) ## GM PF
      E[m,,j] <- (as.vector(Wstat[,j]) >= tau) / (1 + sum(as.vector(Wstat[,j]) <= -tau))
    }
  }
  ## run e-BH
  Em <- apply(E,c(2,3),mean)
  E_ord <- order(as.vector(Em), decreasing = TRUE)
  Es <- sort(Em, decreasing = TRUE)
  comp <- Es >= (1 / q / (1:(p*K)))
  plot(Es,col=comp+2)
  points((1 / q / (1:(p*K))),col=7,pch=3)
  #comp <- E >= ((p*K) / q / (1:(p*K)))
  id <- max(which(comp > 0))
  if(id > 0){
    rej <- E_ord[1:id]
  }else{
    rej <- NULL
  }
  return(sort(rej))
}


runmethodlasso <- function(Y,XX,X,typeset="square",trimset=1){
    if (typeset=="square"){
      all.residuals <- matrix(NA,nrow(Y),ncol(Y))
      for (i in 1:ncol(Y)){
        lam <- cv.glmnet(X,Y[,i])$lambda.1se
        res <- glmnet::glmnet(X,Y[,i],lambda=lam)
        all.residuals[,i] <- Y[,i]-predict(res,X)
      }      
      Omega <- solve(cov(all.residuals))
      #Omega <- diag(ncol(Y))
    }
    else{
      #all.residuals <- matrix(NA,nrow(Y),ncol(Y))
      #for (i in 1:ncol(Y)){
      #  all.residuals[,i] <- sparseLTS(XX,Y[,i])$residuals
      #}
      #Omega <- solve(covMcd(all.residuals)$cov)
      #Omega <- diag(ncol(Y))
      #require(mvoutlier)
      out <- pcout(cbind(Y,X))$wfinal01
      all.residuals <- matrix(NA,sum(out==1),ncol(Y))
      for (i in 1:ncol(Y)){
        lam <- cv.glmnet(X[out==1,],Y[out==1,i])$lambda.1se
        res <- glmnet(X[out==1,],Y[out==1,i],lambda=lam)
        all.residuals[,i] <- Y[out==1,i]-predict(res,X[out==1,])
      }      
      Omega <- solve(cov(all.residuals))
    }
    fit.cv <- MVRLC.cv(XX, Y, K=5, trim=trimset, Omega=Omega,lambda1.min=0.2, lambda1.max=2, 
                       lambda2.min=0.005,  lambda2.max=0.5, nlambda1=2, nlambda2=10, 
                       type=typeset, M1=4.685, M2=4.685, diag=FALSE, 
                       weight=1, weighto=1)
    #plot(fit.cv$error.CV,main=typeset)
    fit <- MVRLC(Y, XX, type=typeset, Omega=Omega,M1=4.685, M2=4.685, 
                 lambda1=fit.cv$lambda1,lambda2=fit.cv$lambda2, 
                 weight=1, weighto=1)
  Wstat <- computeLCDnew(fit$beta)
  tstat <- knockoff.threshold(as.vector(Wstat), fdr = 0.2, offset = 0)
  tpstat <- knockoff.threshold(as.vector(Wstat), fdr = 0.2, offset = 1)
  Sstat <- which(Wstat >= tstat)
  Spstat <- which(Wstat >= tpstat)
  #var.evalues <- evalues(X=X,Y=Y,fit.cv=fit.cv,typeset=typeset,Omega=Omega,M=20,q=0.2)
  list(Sstat=Sstat,Spstat=Spstat,tstat=tstat,tpstat=tpstat,Wstat=Wstat,beta=fit$beta)
      # var.evalues=var.evalues)
}

runmethodlassoE <- function(Y,X,typeset="square",trimset=1){
  if (typeset=="square"){
    all.residuals <- matrix(NA,nrow(Y),ncol(Y))
    for (i in 1:ncol(Y)){
      lam <- cv.glmnet(X,Y[,i])$lambda.1se
      res <- glmnet::glmnet(X,Y[,i],lambda=lam)
      all.residuals[,i] <- Y[,i]-predict(res,X)
    }      
    Omega <- solve(cov(all.residuals))
    Xkf <- create.second_order_PF(X,robust=FALSE,shrink=TRUE)
  }
  else{
    out <- pcout(cbind(Y,X))$wfinal01
    all.residuals <- matrix(NA,sum(out==1),ncol(Y))
    for (i in 1:ncol(Y)){
      lam <- cv.glmnet(X[out==1,],Y[out==1,i])$lambda.1se
      res <- glmnet::glmnet(X[out==1,],Y[out==1,i],lambda=lam)
      all.residuals[,i] <- Y[out==1,i]-predict(res,X[out==1,])
    }      
    Omega <- solve(cov(all.residuals))
    Xkf <- create.second_order_PF(X,robust=TRUE,shrink=TRUE)
  }
  Xaug <- cbind(X, Xkf)
  if (typeset=="square"){
    XX <- as.matrix(scale(Xaug))
  }
  else{
    XX <- as.matrix(scale(Xaug,apply(Xaug,2,median),apply(Xaug,2,mad)))
  }
  fit.cv <- MVRLC.cv(XX, Y, K=5, trim=trimset, Omega=Omega,lambda1.min=0.2, lambda1.max=2, 
                     lambda2.min=0.005,  lambda2.max=0.5, nlambda1=2, nlambda2=10, 
                     type=typeset, M1=4.685, M2=4.685, diag=FALSE, 
                     weight=1, weighto=1,max_iters = 200)
  #plot(fit.cv$error.CV,main=typeset)
  fit <- MVRLC(Y, XX, type=typeset, Omega=Omega,M1=4.685, M2=4.685, 
               lambda1=fit.cv$lambda1,lambda2=fit.cv$lambda2, 
               weight=1, weighto=1,max_iters = 200)
  Wstat <- computeLCDnew(fit$beta)
  tstat1 <- knockoff.threshold(as.vector(Wstat), fdr = 0.2, offset = 0)
  Sstat1 <- which(Wstat >= tstat1)
  tpstat1 <- knockoff.threshold(as.vector(Wstat), fdr = 0.2, offset = 1)
  Spstat1 <- which(Wstat >= tpstat1)
  Sstat <- Spstat <- NULL
  for (j in 1:ncol(Y)){
    tstat <- knockoff.threshold(as.vector(Wstat[,j]), fdr = 0.2, offset = 0)
    Sstat <- c(Sstat,which(Wstat[,j] >= tstat)+(j-1)*ncol(X))
    tpstat <- knockoff.threshold(as.vector(Wstat[,j]), fdr = 0.2, offset = 1)
    Spstat <- c(Sstat,which(Wstat[,j] >= tpstat)+(j-1)*ncol(X))
  }
  var.evalues <- evalues(X=X,Y=Y,fit.cv=fit.cv,typeset=typeset,Omega=Omega,M=100,q=0.2)
  list(Sstat=Sstat,Spstat=Spstat,tstat=tstat,tpstat=tpstat,
       Sstat1=Sstat1,Spstat1=Spstat1,tstat1=tstat1,tpstat1=tpstat1,
       Wstat=Wstat,beta=fit$beta, var.evalues=var.evalues)
}

runmethodlassoEfilter <- function(Y,X,typeset="square",trimset=1){
  if (typeset=="square"){
    res <- cv.glmnet(X,Y,family="mgaussian")
    #res <- cv.glmnet(scale(X),scale(Y),standardize=FALSE,family="mgaussian")
    coefm <- matrix(NA,nrow=ncol(X),ncol=ncol(Y))
    for (i in 1:ncol(Y)){
      coefm[,i] <- as.vector(coef(res)[[i]][-1])
    }      
    var.screen <- which(apply(coefm,1,sum)!=0)
    #var.screen <- which(as.vector(coefm)!=0)
    all.residuals <- Y-drop(predict(res,newx=X))
    Omega <- solve(cov(all.residuals))
    Xkf <- create.second_order_PF(X[,var.screen],robust=FALSE,shrink=TRUE)
  }
  else{
    out <- pcout(cbind(Y,X))$wfinal01
    res <- cv.glmnet(X[out==1,],Y[out==1,],family="mgaussian")
    #res <- cv.glmnet(scale(X,apply(X,2,median),apply(X,2,mad)),scale(Y,apply(Y,2,median),apply(Y,2,mad)),
    #                  standardize=FALSE,family="mgaussian")
    coefm <- matrix(NA,nrow=ncol(X),ncol=ncol(Y))
    for (i in 1:ncol(Y)){
      coefm[,i] <- as.vector(coef(res)[[i]][-1])
    }      
    var.screen <- which(apply(coefm,1,sum)!=0)
    #var.screen <- which(as.vector(coefm)!=0)
    all.residuals <- Y[out==1,]-drop(predict(res,newx=X[out==1,]))
    Omega <- solve(cov(all.residuals))
    Xkf <- create.second_order_PF(X[,var.screen],robust=TRUE,shrink=TRUE)
  }
  Xaug <- cbind(X[,var.screen], Xkf)
  if (typeset=="square"){
    XX <- as.matrix(scale(Xaug))
  }
  else{
    XX <- as.matrix(scale(Xaug,apply(Xaug,2,median),apply(Xaug,2,mad)))
  }
  fit.cv <- MVRLC.cv(XX, Y, K=5, trim=trimset, Omega=Omega,lambda1.min=0.2, lambda1.max=2, 
                     lambda2.min=0.005,  lambda2.max=0.5, nlambda1=2, nlambda2=10, 
                     type=typeset, M1=4.685, M2=4.685, diag=FALSE, 
                     weight=1, weighto=1,max_iters = 200)
  #plot(fit.cv$error.CV,main=typeset)
  fit <- MVRLC(Y, XX, type=typeset, Omega=Omega,M1=4.685, M2=4.685, 
               lambda1=fit.cv$lambda1,lambda2=fit.cv$lambda2, 
               weight=1, weighto=1,max_iters=200)
  Wstat <- computeLCDnew(fit$beta)
  tstat1 <- knockoff.threshold(as.vector(Wstat), fdr = 0.2, offset = 0)
  Sstat1 <- which(Wstat >= tstat1)
  tpstat1 <- knockoff.threshold(as.vector(Wstat), fdr = 0.2, offset = 1)
  Spstat1 <- which(Wstat >= tpstat1)
  Sstat <- Spstat <- NULL
  for (j in 1:ncol(Y)){
    tstat <- knockoff.threshold(as.vector(Wstat[,j]), fdr = 0.2, offset = 0)
    Sstat <- c(Sstat,which(Wstat[,j] >= tstat)+(j-1)*ncol(X))
    tpstat <- knockoff.threshold(as.vector(Wstat[,j]), fdr = 0.2, offset = 1)
    Spstat <- c(Sstat,which(Wstat[,j] >= tpstat)+(j-1)*ncol(X))
  }
  var.evalues <- evalues.filter(X=X,Y=Y,fit.cv=fit.cv,var.screen,typeset=typeset,Omega=Omega,M=100,q=0.2)
  list(Sstat=Sstat,Spstat=Spstat,tstat=tstat,tpstat=tpstat,
       Sstat1=Sstat1,Spstat1=Spstat1,tstat1=tstat1,tpstat1=tpstat1,
       Wstat=Wstat,beta=fit$beta, var.evalues=var.evalues,var.screen=var.screen)
}

runmethodlassoEfast <- function(Y,X,typeset="square",trimset=1,M=50,q=0.2,filter=FALSE){
  # M replications for e-values
  var.screen <- NULL
  if (typeset=="square"){
    if (filter==TRUE){
      res <- cv.glmnet(X,Y,family="mgaussian")
      coefm <- matrix(NA,nrow=ncol(X),ncol=ncol(Y))
      for (i in 1:ncol(Y)){
        coefm[,i] <- as.vector(coef(res)[[i]][-1])
      }      
      var.screen <- which(apply(coefm,1,sum)!=0)
      all.residuals <- Y-drop(predict(res,newx=X))
      X <- X[,var.screen]
    }
    else{
      all.residuals <- matrix(NA,nrow(Y),ncol(Y))
      for (i in 1:ncol(Y)){
        lam <- cv.glmnet(X,Y[,i])$lambda.1se
        res <- glmnet::glmnet(X,Y[,i],lambda=lam)
        all.residuals[,i] <- Y[,i]-predict(res,X)
      }
    }
    Omega <- solve(cov(all.residuals))
    mu <- colMeans(X)
    Sigma <- matrix(as.numeric(corpcor::cov.shrink(X, verbose = F)), nrow = ncol(X))
    #Xkf <- create.second_order_PF(X,robust=FALSE,shrink=TRUE)
  }
  else{
    out <- pcout(cbind(Y,X))$wfinal01
    if (filter==TRUE){
      res <- cv.glmnet(X[out==1,],Y[out==1,],family="mgaussian")
      coefm <- matrix(NA,nrow=ncol(X),ncol=ncol(Y))
      for (i in 1:ncol(Y)){
        coefm[,i] <- as.vector(coef(res)[[i]][-1])
      }      
      var.screen <- which(apply(coefm,1,sum)!=0)
      all.residuals <- Y[out==1,]-drop(predict(res,newx=X[out==1,]))
      X <- X[,var.screen]
    }
    else{
      all.residuals <- matrix(NA,sum(out==1),ncol(Y))
      for (i in 1:ncol(Y)){
        lam <- cv.glmnet(X[out==1,],Y[out==1,i])$lambda.1se
        res <- glmnet::glmnet(X[out==1,],Y[out==1,i],lambda=lam)
        all.residuals[,i] <- Y[out==1,i]-predict(res,X[out==1,])
      }    
    }
    Omega <- solve(cov(all.residuals))
    res <- covOGK(X,sigmamu=robustbase::scaleTau2)
    mu <- res$center
    Sigma <- res$cov
    #Xkf <- create.second_order_PF(X,robust=TRUE,shrink=TRUE)
  }
  method <- "asdp"
  if (nrow(Sigma) <= 500){method <- "sdp"}
  diag_s <- diag(switch(method, sdp = create.solve_sdp(Sigma), asdp = create.solve_asdp(Sigma)))
  if (is.null(dim(diag_s))) {
    diag_s = diag(diag_s, length(diag_s))
  }
  SigmaInv_s <- solve(Sigma, diag_s)
  mu_k <- X - sweep(X, 2, mu, "-") %*% SigmaInv_s
  Sigma_k <- 2 * diag_s - diag_s %*% SigmaInv_s
  cholSigma_k <- chol(Sigma_k) # this matrix will be used later!
  ####
  Xkf <- mu_k + matrix(rnorm(ncol(X) * nrow(X)), nrow(X)) %*% cholSigma_k
  Xaug <- cbind(X, Xkf)
  if (typeset=="square"){
    XX <- as.matrix(scale(Xaug))
  }
  else{
    XX <- as.matrix(scale(Xaug,apply(Xaug,2,median),apply(Xaug,2,mad)))
  }
  fit.cv <- MVRLC.cv(XX, Y, K=5, trim=trimset, Omega=Omega,lambda1.min=0.2, lambda1.max=2, 
                     lambda2.min=0.005,  lambda2.max=0.5, nlambda1=2, nlambda2=10, 
                     type=typeset, M1=4.685, M2=4.685, diag=FALSE, 
                     weight=1, weighto=1,max_iters = 200)
  #plot(fit.cv$error.CV,main=typeset)
  fit <- MVRLC(Y, XX, type=typeset, Omega=Omega,M1=4.685, M2=4.685, 
               lambda1=fit.cv$lambda1,lambda2=fit.cv$lambda2, 
               weight=1, weighto=1,max_iters = 200)
  Wstat <- computeLCDnew(fit$beta)
  tstat1 <- knockoff.threshold(as.vector(Wstat), fdr = 0.2, offset = 0)
  Sstat1 <- which(Wstat >= tstat1)
  Sstat <- NULL
  for (j in 1:ncol(Y)){
    tstat <- knockoff.threshold(as.vector(Wstat[,j]), fdr = 0.2, offset = 0)
    Sstat <- c(Sstat,which(Wstat[,j] >= tstat)+(j-1)*ncol(X))
  }
  #var.evalues <- evalues(X=X,Y=Y,fit.cv=fit.cv,typeset=typeset,Omega=Omega,M=100,q=0.2)
  #evalues <- function(X,Y,fit.cv,typeset,Omega,M=100,q=0.2){
    # Generate knockoff matrix:
    p <- ncol(X)
    K <- ncol(Y)
    E <- array(NA,dim=c(M,p,K))
    all_rej <- c()
    for(m in 1:M){
      if(m %% 20 == 1){
        cat(sprintf("Implementing the %d-th run.\n", m))
      }
      Xkf <- mu_k + matrix(rnorm(ncol(X) * nrow(X)), nrow(X)) %*% cholSigma_k
      Xaug <- cbind(X, Xkf)
      if (typeset=="square"){
        XX <- as.matrix(scale(Xaug))
      }
      else{
        XX <- as.matrix(scale(Xaug,apply(Xaug,2,median),apply(Xaug,2,mad)))
      }
      fit <- MVRLC(Y, XX, type=typeset, M1=4.685, M2=4.685,Omega=Omega, 
                   lambda1=fit.cv$lambda1,lambda2=fit.cv$lambda2, 
                   weight=1, weighto=1,max_iters = 200)
      Wstat <- computeLCDnew(fit$beta)
      #tau <- knockoff.threshold(as.vector(Wstat), fdr = q/2 , offset = 0) ## GM PF
      #E[m,,] <- (as.vector(Wstat) >= tau) / (1 + sum(as.vector(Wstat) <= -tau))
      for (j in 1:K){
        tau <- knockoff.threshold(as.vector(Wstat[,j]), fdr = q/2, offset = 0) ## GM PF
        #tau <- knockoff.threshold(as.vector(Wstat[,j]), fdr = q/2, offset = 0) ## GM PF
        E[m,,j] <- (as.vector(Wstat[,j]) >= tau) / (1 + sum(as.vector(Wstat[,j]) <= -tau))
      }
    }
    ## run e-BH
    Em <- apply(E,c(2,3),mean)
    E_ord <- order(as.vector(Em), decreasing = TRUE)
    Es <- sort(Em, decreasing = TRUE)
    comp <- Es >= (1 / q / (1:(p*K)))
    #plot(Es,col=comp+2)
    #points((1 / q / (1:(p*K))),col=7,pch=3)
    #comp <- E >= ((p*K) / q / (1:(p*K)))
    id <- max(which(comp > 0))
    if(id > 0){
      rej <- E_ord[1:id]
    }else{
      rej <- NULL
    }
    var.evalues <- sort(rej)
  
  list(Sstat=Sstat,tstat=tstat,Sstat1=Sstat1,tstat1=tstat1,
       Wstat=Wstat,beta=fit$beta, var.evalues=var.evalues,var.screen=var.screen)
}


gendata <- function(n=100,p=200,K=2,outX=FALSE,outE=FALSE,propout=0.1,facsigma=10,facerror=10,sp=0.1,errorRho = .7){
#  n        # number of samples
#  p        # number of features (without reference)
  pstar = p + 1  #number of features (including reference)
#  K = 2          # number of responses
  A = 3          # signal amplitude
#  sp = .05       # sparsity percentage
#  errorRho = .7  # error rho for E
  
  # Generate covariance of W
  sigma = matrix(0, nrow = pstar, ncol = pstar)
  for(i in 1:pstar){
    for(j in 1:pstar){
      sigma[i,j] = .5^(abs(i-j))
    }
  }
  
  # generate W (latent log-normal data)
  means = rep(1, pstar)
  W = mvlognormal(n, means, Sigma = diag(sigma), R = sigma)
  colnames(W) = paste("Taxa", 1:ncol(W))
  
  if (outX==TRUE & propout>0){
    n1 <- round(propout*n) # number of outliers
   # Wout = mvlognormal(n1, means*10, Sigma = diag(sigma), R = sigma)
    Wout = mvlognormal(n1, means, Sigma = facsigma*diag(sigma), R = diag(diag(sigma)))
  }
  
  # generate B (coefficient matrix)
  Aset = seq(-A,A)
  if(length(which(Aset==0))!=0){
    Aset = Aset[-which(Aset==0)]
  }
  B = matrix(0, nrow = p, ncol = K)
  totalSig = p*K
  nSig = floor(totalSig*sp)
  signalEntries = sort(sample(1:totalSig, size = nSig, replace = F))
  Bset = sample(Aset, nSig, replace = T)
  B[signalEntries] = Bset
  
  # generate error matrix E
  errorCovariance = matrix(0,K, K)
  for(i in 1:K){
    for(j in 1:K){
      errorCovariance[i,j] = errorRho^abs(i-j)
    }
  }
  
  E = matrix(0, nrow = n, ncol = K)
  for(i in 1:n){
    ei = rnorm(K, mean = rep(0, K), errorCovariance)
    E[i,] = ei
  }
  
  if (outE==TRUE & propout>0){
    n1 <- round(propout*n)
    Eout <- matrix(0, nrow = n1, ncol = K)
    for(i in 1:n1){
      eiout = rnorm(K, mean = rep(0, K), facerror*errorCovariance)
      #eiout = rnorm(K, mean = rep(facerror, K), errorCovariance)
      Eout[i,] = eiout
    }
    E[1:n1,] <- Eout
  }
  
  # generate Z (compositional matrix)
  Z = acomp(W)
  Z = as.matrix(Z)
  
  if (outX==TRUE & propout>0){
    Zout <- as.matrix(acomp(Wout))
  }
  
  # generate X (alr-transformation matrix)
  ref = Z[,pstar]
  X = matrix(0, nrow = n, ncol = p)
  for(t in 1:ncol(X)){
    X[,t] = log(Z[,t]/ref)
  }
  
  if (outX==TRUE & propout>0){
    refout = Zout[,pstar]
    Xout = matrix(0, nrow = n1, ncol = p)
    for(t in 1:ncol(X)){
      Xout[,t] = log(Zout[,t]/refout)
    }
  }
  
  # generate Y (response matrix)
  Y = X%*%B + E
  colnames(Y) = paste("Response", 1:ncol(Y))
  if (outX==TRUE & propout>0){
    X[1:n1,] <- Xout
  }
  
  list(X=X,Y=Y,signalEntries=signalEntries)
}

evalsim <- function(listres,listtrue){
  nsim <- length(listres)
  fdr <- tpr <- rep(NA,nsim)
  for (i in 1:nsim){
    if (is.numeric(listres[[i]])){
      fdr[i] <- sum((listres[[i]]%in%listtrue[[i]])==FALSE)/length(listres[[i]])
      tpr[i] <- sum(listres[[i]]%in%listtrue[[i]])/length(listtrue[[i]])
    }
  }
  if (any(is.na(c(fdr,tpr)))){
    warning("Nas in fdr/tpr!")
  }
  list(fdr=fdr,tpr=tpr)
}

evalsimF <- function(listres,listtrue,listvar.screen,p=200,K=5){
  # for the filter procedure
  nsim <- length(listres)
  fdr <- tpr <- rep(NA,nsim)
  for (i in 1:nsim){
    varsel <- listvar.screen[[i]]
    for (j in 1:(K-1)){
      varsel <- c(varsel,listvar.screen[[i]]+j*p)
    }
    varsel <- sort(varsel)
    listresi <- varsel[listres[[i]]]
    if (is.numeric(listresi)){
      fdr[i] <- sum((listresi%in%listtrue[[i]])==FALSE)/length(listresi) # FDR
      tpr[i] <- sum(listresi%in%listtrue[[i]])/length(listtrue[[i]]) # TPR
    }
  }
  if (any(is.na(c(fdr,tpr)))){
    warning("Nas in fdr/tpr!")
  }
  list(fdr=fdr,tpr=tpr)
}


create.second_order_PF <- function (X, method = c("asdp", "equi", "sdp"), shrink = TRUE, robust=TRUE) 
{
  method = match.arg(method)
  if (!shrink) {
    if (robust){
      require(rrcov)
      res <- CovRobust(X)
      mu <- res@center
      Sigma <- res@cov
    }
    else{
      mu <- colMeans(X)
      Sigma <- cov(X)
    }
    if (!is_posdef(Sigma)) {
      shrink = TRUE
    }
  }
  if (shrink) {
    if (robust){
      #res <- mrcd(t(X),h=0.75*nrow(X))
      #mu <- res$mu
      #Sigma <- res$cov
      res <- covOGK(X,sigmamu=robustbase::scaleTau2)
      mu <- res$center
      Sigma <- res$cov
    }
    else{
      mu <- colMeans(X)
      if (!requireNamespace("corpcor", quietly = T)) 
        stop("corpcor is not installed", call. = F)
      Sigma = tryCatch({
        suppressWarnings(matrix(as.numeric(corpcor::cov.shrink(X, 
                verbose = F)), nrow = ncol(X)))
      }, warning = function(w) {
      }, error = function(e) {
        stop("SVD failed in the shrinkage estimation of the covariance matrix. Try upgrading R to version >= 3.3.0")
      }, finally = {
      })
    }
  }
  create.gaussian(X, mu, Sigma, method = method)
}

