########################################################################
# 1. generate data
########################################################################
# (1) gen predictors: x
genX <- function(n, p, rho){
  CovMatrix <- outer(1:p, 1:p, function(x,y) {rho^abs(x-y)})
  x <- mvrnorm(n, rep(0,p), CovMatrix)
  x
}

# (2) gen random error: eps
genE <- function(n, error='normal'){
  if (error=='normal') eps <- rnorm(n, mean=0, sd=1)
  if (error=='t') eps <- rt(n,df=3)
  if (error=='chisq') eps <- rchisq(n, df=3)
  eps
}

# (3) gen simulations
genSim <- function(n, p, rho, error='normal',mode="iid"){
  x <- genX(n=n, p=p, rho=rho)
#   x1 <- runif(n,-5,5)
#   x2 <- runif(n,0,1)
#   x3 <- runif(n,-2,2)
#   x4 <- runif(n,-2,2)
#   x5 <- rnorm(n,sd=2)
#   x <- cbind(x1,x2,x3,x4,x5)
  eps <- genE(n=n, error=error)
  SN <- 10
  y <- x[,1]^2+2*exp(-5*x[,2]^2)+x[,3]*x[,4]+x[,5]
  sigma <- sqrt(var(y)/var(eps)/10)
  if(mode=="iid"){
    y <- y +as.vector(sigma)*eps
  }else{
    y <- y + as.vector(sigma)*(1+x[,1])*eps
  }
  dat <- data.frame(y=y,x1=x[,1],x2=x[,2],
                    x3=x[,3],x4=x[,4],x5=x[,5])
  dat
}

############################################
# 2.Compute precision indexes
############################################
RMSE <- function(y, yfit){
  if(is.matrix(y)==T){
    rmse <- sqrt(apply((y-yfit)^2,2,mean))
  }else{
    rmse <- sqrt(mean((y-yfit)^2))
  }
  rmse
}

MAD <- function(y, yfit){
  if(is.matrix(y)==T){
    mad <- apply(abs(y-yfit),2,mean)
  }else{
    mad <- mean(abs(y-yfit))
  }
  mad
}

############################################
# 3.expectile regression forest
############################################
qrf <- opt_QRF(Xtrain, Ytrain) 
X<-Xtrain
Y<-Ytrain

Mse <- matrix(NA,ncol(Xtrain),1)
opt_QRF <- function(X, Y){
  qrf <- list()
  mse <- numeric()
  for(j in 1:ncol(X)){
    qrf[[j]] <- quantregForest(x=X, y=Y, nodesize=5,mtry=j,ntree=1000)
    mse[j] <- mean(qrf[[j]]$mse)
    Mse[j,]<-mse[j]
  }
  #opt_mtry <- which.min(mse)
  qrf[[which.min(mse)]]
#   mse <- cumsum(qrf$mse)/cumsum(1:length(qrf$mse))
#   opt_ntree <-NA
#   for(i in 1:(length(mse)-1)){
#     if(abs(mse[i+1]-mse[i])<10^-4){
#       opt_ntree <- i
#       break
#     }
#   }
#   if(is.na(opt_ntree)) opt_ntree <- length(mse)
#   qrf <- quantregForest(x=X, y=Y, nodesize=5,
#                         mtry=opt_mtry,ntree=opt_ntree)
#   qrf
}

# opt_QRF <- function(X, Y){
#   #qrf <- list()
#   mse <- numeric()
#   for(j in 1:ncol(X)){
#     qrf <- quantregForest(x=X, y=Y, nodesize=5,mtry=j,ntree=100)
#     mse[j] <- mean(qrf$mse)
#   }
#   
#   qrf <- list()
#   msee <- numeric()
#   for(i in 2:20){
#     qrf[[i]] <- quantregForest(x=X, y=Y, nodesize=i,
#                                mtry=which.min(mse),ntree=100)
#     msee[i] <- mean(qrf[[i]]$mse)
#   }
#   #opt_mtry <- which.min(mse)
#   qrf[[which.min(na.omit(msee))+1]]
# }


expectile_beta <- function(z,tau){
  z <- na.omit(z)
  z <- sort(z)
  S <- length(z)
  tmp1 <- tmp2 <- matrix(NA,S+1,S)
  for(k in 1:(S+1)){
    for(s in 1:S){
      tmp1[k,s] <- (1-tau)*z[s]*(s<=(k-1))+tau*z[s]*(s>=k)
      tmp2[k,s] <- (1-tau)*(s<=(k-1))+tau*(s>=k)
    }
  }
  beta <- apply(tmp1,1,sum)/apply(tmp2,1,sum)
  for(k in 1:(S+1)){
    if(k==1){
      if(beta[k]<=z[k]) betas <- beta[k]
    }else if(k==S+1){
      if(beta[k]>=z[k-1]) betas <- beta[k]
    }else{
      if(beta[k]>=z[k-1]&beta[k]<=z[k]) betas <- beta[k]
    } 
  }
  betas
}


predict.QRF <- function(object,newdata=NULL,what=c(0.1,0.5,0.9),method="expectile"){
  class(object) <- "randomForest"
  predictNodes <- attr(predict(object,newdata=newdata,nodes=TRUE),"nodes")
  rownames(predictNodes) <- NULL
  valuesPredict <- 0*predictNodes
  ntree <- ncol(object[["valuesNodes"]])
  for (tree in 1:ntree){
    valuesPredict[,tree] <- object[["valuesNodes"]][predictNodes[,tree],tree]  
  }
  
  if(method=="mean") result <- apply(valuesPredict,1,mean,na.rm=TRUE)
  #result_MRF <- apply(valuesPredict,1,mode,na.rm=TRUE)
  if(method=="quantile") result <- t(apply(valuesPredict,1,quantile, what,na.rm=TRUE))
  if(method=="expectile"){
    result <- matrix(NA,nrow(valuesPredict),length(what))
    for(i in 1:length(what)){
      result[,i] <- apply(valuesPredict,1,expectile_beta, what[i]) 
    }
  }
  result
}

############################################
# 4.expectile regression
############################################
exptreg <- function(X, y, tau,  max.iter, tol){
  # X is the model matrix
  # y is the response vector of observed proportion
  # maxIter is the maximum number of iterations
  # tol is a convergence criterion
  X <- cbind(1, X) # add constant
  b <- bLast <- rep(0, ncol(X)) # initialize
  it <- 1 # iteration index
  while (it <= max.iter){
    
    ypred <- c(X %*% b)
    w <- as.vector(tau *(y>= ypred) + (1-tau)* (y<ypred))
    #w <- diag(w)
    b <- solve(t(X)%*%diag(w)%*%X)%*%(t(X)%*%diag(w)%*%y)
    b <- as.vector(b)
    #b <- lsfit(X, y, w, intercept=FALSE)$coef#lsfit() weighted Least Squard
    if (max(abs(b - bLast)/(abs(bLast) + 0.01*tol)) < tol) break
    bLast <- b
    it <- it + 1 # increment index
  }
  if (it > max.iter) warning('maximum iterations exceeded')
  
  b
}

############################################
# 5.expectile regression boosting tree
############################################
cv_boost <- function(form, data, tau){
  minobs <- c(2:10)
  error <- numeric()
  for(i in 1:length(minobs)){
    fit <- erboost(form, data, distribution=list(name="expectile",alpha=tau),
                   n.trees=3000, n.minobsinnode = minobs[i],cv.folds = 10)
    best.iter <- erboost.perf(fit,method="cv")
    error[i] <- mean(fit$cv.error[1:best.iter])
  }
  minobs <- minobs[which.min(error)]
  minobs
}


predictTree <- function(y,taus,method="mean"){
  if(method=="mean"){
    yhat <- mean(y)
  }
  if(method=="quantile"){
    yhat <- quantile(y, taus)
  }
  if(method=="expectile"){
    yhat <- expectile_beta(y, taus)
  }
  yhat
}


quantreg_tree <- function(fit,dat_train,dat_test,Index,tau,method="quantile"){
  pred_train <- predict(fit,dat_train)
  pred_test <- predict(fit,dat_test)
  freq <- table(pred_train)
  yhat <- c()
  yhat_test <- c()
  for(i in 1:length(freq)){
    ind1 <- which(pred_train==names(freq)[i])
    y <- dat_train[ind1,Index]
    yhat[ind1] <- predictTree(y,tau,method)
    #yhat[ind1] <- quantile(y,tau)
    
    ind2 <- which(pred_test==names(freq)[i])
    yhat_test[ind2] <- predictTree(y,tau,method)
    #yhat_test[ind2] <- quantile(y,tau)
  }
  ans <- list()
  ans$yhat <- yhat
  ans$yhat_test <- yhat_test
  ans
}
CVgroup <- function(K,datasize,random=T){
  cvlist <- list()
  n <- rep(1:K,each=ceiling(datasize/K))[1:datasize]
  if(random==T){
    n <- sample(n,datasize)
  }
  x <- 1:K
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[n==x])  #dataseq中随机生成k个随机有序数据列
  return(cvlist)
}


cv_Tree <- function(data, Index, K){
  cvlist <- CVgroup(K = K,datasize = nrow(data))
  obsinnode <- c(2:10)
  error <- matrix(NA,length(obsinnode),K)
  for(i in 1:length(obsinnode)){
    for(j in 1:K){
      train <- data[-cvlist[[j]],]
      test <- data[cvlist[[j]],]
      form <- as.formula(paste(names(train)[1],"~."))
      fit <- rpart(form, train, control=list(minbucket =obsinnode[i])) 
      tmp <- predict(fit,test)
      error[i,j] <- RMSE(test[,Index], tmp)
    } 
  }
  ind <- which.min(apply(error,1,mean))
  minobsinnode <- obsinnode[ind]
  minobsinnode
}

# (2) gen random parameter
genParm <- function(l, p){
  beta <- runif(l,-1,1)
  pl <- c()
  ul <- vl <- list()
  for(i in 1:l){
    pl[i] <- min(ceiling(1.5+rexp(1,0.5)),p)
    ul[[i]] <- rnorm(pl[i])
    dl <- runif(pl[i],0.1,2)^2
    dl <- diag(dl)
    Ul <- qr.Q(qr(matrix(rnorm(pl[i]*pl[i]),pl[i],pl[i])))
    vl[[i]] <- Ul%*%dl%*%t(Ul)
  }
  ans <- list()
  ans$beta <- beta
  ans$pl <- pl
  ans$ul <- ul
  ans$vl <- vl
  ans
}

# (3) gen simulations
genSim <- function(n, p, l, error='normal',mode="iid"){
  x <- matrix(rnorm(n*p),n,p)
  eps <- genE(n=n, error=error)
  parms <- genParm(l, p)
  beta <- parms$beta
  ul <- parms$ul
  pl <- parms$pl
  vl <- parms$vl
  #SN <- 10
  y <- c()
  for(i in 1:n){
    gl <- c()
    for(j in 1:length(pl)){
      zl <- x[i,1:pl[j]]
      gl[j] <- exp(-0.5*(t(zl-ul[[j]])%*%vl[[j]]%*%(zl-ul[[j]])))
    }
    y[i] <- t(beta)%*%gl
  }
  sigma <- sqrt(var(y)/var(eps)/10)
  if(mode=="iid"){
    y <- y +as.vector(sigma)*eps
  }else{
    y <- y + as.vector(sigma)*(1+x[,1])*eps
  }
  dat <- data.frame(y=y,x=x)
  dat
}
