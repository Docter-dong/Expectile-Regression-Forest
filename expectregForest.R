setwd("D:/expectile regression forest/R")
rm(list = ls())

library(quantregForest)
library(MASS)
library(quantreg)
library(erboost)
# library(openxlsx)
source("sub-forest.R")
source('ernn_L1L2.R')
library("expectreg")
library("rpart")


#dat <- read.csv("Folds5x2_pp1.csv")
#dat <- cbind(dat[,5],dat[,1:4])
#data(airquality)
#dat <- airquality[ !apply(is.na(airquality), 1,any), ]
#dat <- dat[,c(1:4)]

#dat<-read.table("machine.data",sep=",")
#dat <- dat[,-c(1:2)]
#dat <- cbind(dat[,8],dat[,1:7])
#dat[,4] <- as.character(dat[,4])
#dat[,4] <- as.numeric(dat[,4])
#dat <- na.omit(dat)
#
#
#data("Boston", package = "MASS")
#dat <- cbind(Boston[,14],Boston[,c(1,3,5:13)])
#
#dat<- read.csv("winequality-red.csv",sep=";")
#dat <- cbind(dat[,12],dat[,1:11])
  
what <- c(0.05,0.5,0.95)
M <- 100
mad_ERF_test <- rmse_ERF_test <- mad_ER_test <- rmse_ER_test <- 
  mad_ERF_train <- rmse_ERF_train <- mad_ER_train <- rmse_ER_train <- 
  mad_ERGB_test <- rmse_ERGB_test <- mad_ERGB_train <- rmse_ERGB_train <-
  mad_ERT_test <- rmse_ERT_test <- mad_ERT_train <- rmse_ERT_train <-
  mad_ERNN_test <- rmse_ERNN_test <- mad_ERNN_train <- rmse_ERNN_train <-
  mad_ERG_test <- rmse_ERG_test <- mad_ERG_train <- rmse_ERG_train <-
  matrix(NA,M,length(what))
for(m in 1:M){
  ## generate data ##
  set.seed(1+m)
  dat <- genSim(n=1000, p=10, l=20, error='normal',mode="iid")
  #dat <- genSim(n=1000, p=10, l=20, error='t',mode="iid")
  #dat <- genSim(n=1000, p=10, l=20, error='chisq',mode="iid")
  #dat <- genSim(n=1000, p=10, l=20, error='normal',mode="nid")
  #dat <- genSim(n=1000, p=10, l=20, error='t',mode="nid")
  #dat <- genSim(n=1000, p=10, l=20, error='chisq',mode="nid")
  colnames(dat) <- c("y",paste("x",1:(ncol(dat)-1),sep=""))
  ind <- sample(1:nrow(dat),floor(0.7*nrow(dat)))
  dat_train <- dat[ind,]
  dat_test <- dat[-ind,]
  #dat_train <- genSim(n=1000, p=5, rho=0.5, error='normal',mode="iid")
  Xtrain <- as.matrix(dat_train[,-1])
  Ytrain <- dat_train[,1]
  #dat_test <- genSim(n=500, p=5, rho=0.5, error='normal',mode="iid")
  Xtest <- as.matrix(dat_test[,-1])
  Ytest <- dat_test[,1]
  ## compute expectile Regression Forests ##
  qrf <- opt_QRF(Xtrain, Ytrain)  
  res_train <- predict.QRF(object=qrf,newdata=Xtrain,what,method="expectile")
  res_test <- predict.QRF(object=qrf,newdata=Xtest,what,method="expectile")
  ## compute expectile Regression ##
  y_train <- matrix(NA,nrow(Xtrain),length(what))
  y_test <- matrix(NA,nrow(Xtest),length(what))
  for(i in 1:length(what)){
    coefs <- exptreg(X=Xtrain, y=Ytrain, tau=what[i], max.iter=100, tol=1e-8)
    y_train[,i] <- cbind(rep(1,nrow(Xtrain)),Xtrain)%*%coefs
    y_test[,i] <- cbind(rep(1,nrow(Xtest)),Xtest)%*%coefs
  }
  explaws <- expectreg.ls(y~rb(x1,"pspline")+rb(x2,"pspline")+
                            rb(x3,"pspline")+rb(x4,"pspline"),
                          data=dat_train,smooth="aic",
                          expectiles=c(0.05,0.5,0.95))
  y_train_erg <- predict(explaws,dat_train)$fitted
  y_test_erg <- predict(explaws,dat_test)$fitted
  ## compute expectile Regression ##
  y_train_ernn <- matrix(NA,nrow(Xtrain),length(what))
  y_test_ernn <- matrix(NA,nrow(Xtest),length(what))
  for(i in 1:length(what)){
    ernn <- ernn.fit(Xtrain, as.matrix(Ytrain), n.hidden = 3, tau=what[i],iter.max = 5000, n.trials = 2, 
                     eps = 10^-2, penalty = 0)
    y_train_ernn[,i] <- ernn.predict(Xtrain, ernn)
    y_test_ernn[,i] <- ernn.predict(Xtest, ernn)
 }
  
  
  ytr <- cbind(Ytrain,Ytrain,Ytrain)
  yte <- cbind(Ytest,Ytest,Ytest)
  mad_ER_train[m,] <- MAD(ytr, y_train)
  mad_ERF_train[m,] <- MAD(ytr, res_train)
  mad_ER_test[m,] <- MAD(yte, y_test)
  mad_ERF_test[m,] <- MAD(yte, res_test)
  
  rmse_ER_train[m,] <- RMSE(ytr, y_train)
  rmse_ERF_train[m,] <- RMSE(ytr, res_train)
  rmse_ER_test[m,] <- RMSE(yte, y_test)
  rmse_ERF_test[m,] <- RMSE(yte, res_test)
  
  mad_ERG_train[m,] <- MAD(ytr, y_train_erg)
  mad_ERG_test[m,] <- MAD(yte, y_test_erg)
  
  rmse_ERG_train[m,] <- RMSE(ytr, y_train_erg)
  rmse_ERG_test[m,] <- RMSE(yte, y_test_erg)
  
  mad_ERNN_train[m,] <- MAD(ytr, y_train_ernn)
  mad_ERNN_test[m,] <- MAD(yte, y_test_ernn)
  
  rmse_ERNN_train[m,] <- RMSE(ytr, y_train_ernn)
  rmse_ERNN_test[m,] <- RMSE(yte, y_test_ernn)
  ## compute expectile Regression boosting tree ##
  form <- as.formula(paste(names(dat_train)[1],"~.",sep=""))
  ops <- cv_Tree(dat_train, Index=colnames(dat_train)[1], K=10)
  fit <- rpart(form, dat_train,control=list(minbucket =ops))
  
  for(i in 1:length(what)){
    parms <- cv_boost(form, dat_train, tau=what[i])
    model <- erboost(form, data = dat_train,distribution=list(name="expectile",alpha=what[i]),
                     n.trees=3000, n.minobsinnode = parms,cv.folds = 10)
    best.iter <- erboost.perf(model, method = "cv")
    y_insample <- predict.erboost(model,dat_train,best.iter)
    y_outsample <- predict.erboost(model,dat_test,best.iter)
    mad_ERGB_train[m,i] <- MAD(Ytrain,y_insample)
    mad_ERGB_test[m,i] <- MAD(Ytest,y_outsample)
    rmse_ERGB_train[m,i] <- RMSE(Ytrain,y_insample)
    rmse_ERGB_test[m,i] <- RMSE(Ytest,y_outsample)
    
    res <- quantreg_tree(fit,dat_train,dat_test,Index=colnames(dat_train)[1],
                         what[i],method="expectile")
    mad_ERT_train[m,i] <- MAD(dat_train[,1],res$yhat)
    mad_ERT_test[m,i] <- MAD(dat_test[,1],res$yhat_test)
    rmse_ERT_train[m,i] <- RMSE(dat_train[,1],res$yhat)
    rmse_ERT_test[m,i] <- RMSE(dat_test[,1],res$yhat_test)
  }
  print(m)
}

