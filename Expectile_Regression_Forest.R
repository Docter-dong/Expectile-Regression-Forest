########################################################
# Description:
# 1.Purpose: Expectile Regression Forest
# 2.Author: Chao Cai, Haotian Dong
# 3.Founded: July, 2021.
# 4.Revised: April, 2022.
########################################################

########################################################
# 1.Genetate Data
########################################################
#(1) generate random error: eps
genE <- function(n, error='normal'){
  if (error=='normal') eps <- rnorm(n, mean=0, sd=1)
  if (error=='t') eps <- rt(n,df=3)
  if (error=='chisq') eps <- rchisq(n, df=3)
  eps
}

# (2) generate parameter
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
  colnames(dat) <- c("y",paste("x",1:p,sep=""))
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

############################################
# 3.Compute expectile value
############################################
expectile_value <- function(z,tau){
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

############################################
# 4.Optimal Random Forest
############################################
opt_RF <- function(X, Y){
  require(quantregForest) 
  ntrees <- c(1:15)
  nodesizes <- c(2:15)
  mtrys <- 1:ncol(X)
  parm<- matrix(0,nrow=length(ntrees)*length(nodesizes)*length(mtrys),ncol=3)
  parm[,1] <- rep(ntrees,length(nodesizes)*length(mtrys))
  parm[,2] <- rep(mtrys,each=length(ntrees)*length(nodesizes))
  parm[,3] <- rep(rep(nodesizes,each=length(ntrees)),length(mtrys))
  mse <- c()
  for(i in 1:nrow(parm)){
    rf <- quantregForest(x=X, y=Y, nodesize=parm[i,3],mtry=parm[i,2],
                         ntree=100*parm[i,1])
    mse[i]<- mean(rf$mse)
  }
  opt_parm <- parm[which.min(mse),]
  rf <- quantregForest(x=X, y=Y, nodesize=opt_parm[3],
                        mtry=opt_parm[2],ntree=100*opt_parm[1])
  rf
}

############################################
# 5.Predict of Expectile Regression Forest
############################################
predict.ERF <- function(object,newdata=NULL,what=c(0.1,0.5,0.9)){
  class(object) <- "randomForest"
  predictNodes <- attr(predict(object,newdata=newdata,nodes=TRUE),"nodes")
  rownames(predictNodes) <- NULL
  valuesPredict <- 0*predictNodes
  ntree <- ncol(object[["valuesNodes"]])
  for (tree in 1:ntree){
    valuesPredict[,tree] <- object[["valuesNodes"]][predictNodes[,tree],tree]  
  }
  
  result <- matrix(NA,nrow(valuesPredict),length(what))
  for(i in 1:length(what)){
    result[,i] <- apply(valuesPredict,1,expectile_value, what[i]) 
  }
  result
}

#############################################################
# 6. Importance of Expectile regression forest
#############################################################
Importance.ERF <- function(object,X,what){
  predict_value <- predict.ERF(object,X,what)
  Import <- matrix(NA,ncol(X),length(what))
  for(j in 1:ncol(X)){
    Imp <- matrix(NA,100,length(what))
    for(i in 1:100){
      nX <- X
      nX[,j] <- sample(nX[,j],length(nX[,j]))
      pre_value_adjust <- predict.ERF(object,nX,what)
      Imp[i,] <- mean((pre_value-pre_value_adjust)^2)
    }
    Import[j,] <- apply(Imp,2,mean)
  }
  ImportSD <- Import/apply(Import,2,sum)#compute weight
  Imp <- cbind(Import,ImportSD)
  colnames(Imp) <- c(paste("Import values_",what,sep=""),
                     paste("Import weight_",what,sep=""))
  rownames(Imp) <- colnames(X)
  Imp
}
#############################################################
# 7. Partial dependence of Expectile regression forest
#############################################################
partial_dependence.ERF <- function(object,X,what,Varnames){#Varnames is the name of the variable
  ind <- grep(Varnames,colnames(X))[1]
  x <- sort(unique(X[,ind]))
  if(length(x)>50){
    x <- seq(min(X[,ind]),max(X[,ind]),length=50)
  }
  P <- matrix(NA,ncol(x),length(what))
  for(i in 1:length(x)){
    xx <- X
    xx[ind] <- rep(x[i],nrow(X))
    yhat <- predict.ERF(object,xx,what)
    P[i,] <- apply(yhat,2,mean)
  }
  P <- cbind(x,P)
  colnames(P) <- c(Varnames,paste("partial_dependence_",what),sep="")
  P
}

##################################################
# 8.Cross-validation of Expectile regression tree
##################################################
#(1) data grouping
CVgroup <- function(K=10,datasize,random=T){
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

#(2)Cross-validation
cv_Tree <- function(data, K=10){
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
      error[i,j] <- RMSE(test[,1], tmp)
    } 
  }
  ind <- which.min(apply(error,1,mean))
  minobsinnode <- obsinnode[ind]
  minobsinnode
}
############################################
# 9.Fit Optimal Expectile regression tree
############################################
opt_tree <- function(dat){
  form <- as.formula(paste(names(dat)[1],"~."))
  ops <- cv_Tree(dat, K=10)
  fit <- rpart(form, dat,control=list(minbucket =ops))
  fit
}

############################################
# 10.predict of regression tree
############################################
predict.tree <- function(fit,Xtest,taus){
  #required packages:rpart; fit is fitted result of tree
  freq <- table(fit$where)
  pred_test <- predict(fit,Xtest)
  yhat <- matrix(NA,nrow(Xtest),length(taus))
  for(i in 1:length(freq)){
    ind1 <- which(fit$where==as.numeric(names(freq))[i])
    y_mu <- round(mean(fit$y[ind1]),3)
    y_Rl <- fit$y[ind1]#y in each leaf
    ind2 <- which(round(pred_test,3)==y_mu)
    for(i in 1:length(taus)){
      yhat[ind2,i] <- rep(expectile_value(y_Rl,taus[i]),length(ind2))
    }
  }  
  yhat
}




library(rpart)
library(quantregForest)
dat <- genSim(n=200, p=10, l=20, error='normal',mode="iid")
X <- as.matrix(dat[,-1])
Y <- dat[,1]
rf <- opt_RF(X, Y)
ERF.value <- predict.ERF(rf,newdata=X,what=c(0.1,0.5,0.9))

tree <- opt_tree(dat)
ERT.value <- predict.tree(tree,as.data.frame(X),taus=c(0.1,0.5,0.9))
