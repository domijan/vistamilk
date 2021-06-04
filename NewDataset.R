rm(list = ls())
set.seed(1979)
library(tidyverse)
library(dendextend)
library(glmnet) # for the ridge regression
library(pls)
library(ranger)
library(bartMachine)
library(e1071)
library(BKPC)
library(kernlab)
library(caret) # For ten fold crossvalidation
library(effects)
library (stats)
library(lme4)
library(brnn)
library(gbm)
library(haven)

calcRMSE <- function(y, yhat){sqrt(mean((yhat - y)^2))}

mixed.data <- read_sas("dataset_for_mixed_model.sas7bdat")
mixed.data%>% glimpse()
y <- mixed.data$delta_bcs
x <- mixed.data %>% select(PM5:PM752)
# glimpse(x)
matplot(t(x), type = "l")

ys <- as.numeric(y)
a <- mean(ys, na.rm = TRUE)
b <- sd(ys, na.rm = TRUE)

ys <- ys-a
ys <- ys/b


ys %>% as.data.frame() %>% ggplot(aes(x = ys))+ geom_histogram()
# hist(ys)

me <- apply(x,2,mean)
sd <- apply(x,2,sd)

x.sc <- scale(x,me,sd)
matplot(t(x.sc), type = "l")


indexcor <- abs(cor(cbind(ys,x.sc), use = "pairwise.complete.obs"))[1, -1]
# matplot(t(indexcor), type = "l", col = 9, lty = 1)




# clustering the columns (variables) rather than observations
D2 <- as.matrix(1-abs(cor(x.sc)))

# image(D2)


h2<- hclust(as.dist(D2), method = "ward.D2")
dnew <- as.dendrogram(h2)



d2 <- color_branches(dnew,k=15, col = 2:16) 

plot(d2)

# plot(reorder(d2, D2, method = "OLO")) 

labl.cols <- cutree(d2,15) # 15 variables

nonz <- matrix(0,1,15)


for (j in unique(labl.cols)){
  nonz[1, j] <- which(indexcor==max(indexcor[labl.cols==j]))
  }


kfunc <-  rbfdot(sigma = 0.0005)
# kfunc <- laplacedot(sigma = 0.0001)
# kfunc <- anovadot(sigma = 1, degree = 1)

# id_te <- sample(nrow(x), 20000)
# ys_te <- ys[id_te]
# te <- x.sc[id_te, ]
# ys_tr <- ys[-id_te]
# tr <- x.sc[-id_te, ]
# 
# 
# 
# Ktrain <- kernelMatrix(kfunc, tr)
# image(Ktrain[1:1000, 1:1000])
# 
# # Ktest <- kernelMatrix(kfunc, te, tr)



# kpcTrain <- kPCA(Ktrain)


# kpcTest <- predict(kpcTrain, Ktest)



dat <- cbind(ys, x.sc)

dat <- as.data.frame(dat)
p <- dat %>% ggplot(aes(x= ys)) + geom_histogram() 
print(p)
dat <- dat %>% na.omit()
dim(dat)
nonzy <- nonz

save(dat, nonzy, file="temp.RData")
rm(list=ls())
load(file="temp.RData")
pls.fit <- plsr(ys ~., data=dat, scale= TRUE, validation = "none")

N <- 2
grid <- 10^seq(-3, 3, length = 20)
RMSEFULL <- vector(mode = "list", length = 20)
t1 <- vector(mode = "list", length = 1)
names(t1)<- "delta_bcs"
flds <- createFolds(1:20000, k = 5, list = TRUE, returnTrain = FALSE)


# pls.pc <-c(4, 6, 5, 4, 3, 6, 15, 4, 10, 15,8, 4,8,7)

lm.fit.SG <- vector(mode = "list", length = N)
nonneg.fit.SG<-  vector(mode = "list", length = N)
rf.fit.SG <- vector(mode = "list", length = N)
lasso.fit.SG <- vector(mode = "list", length = N)




  
  set.seed(1979)
  ttt <- replicate(N, sample(nrow(dat), 20000))
  
  
  RMSE <- matrix(0, N, 20)
  predictions.te <- matrix(0, nrow(dat)-20000, 15)
  predictions.tv <- matrix(0, 20000, 15)
  
  
  
  for (i in 1:N){
    
    
    dat.tr.full <- dat[ttt[,i],]
    dat.te <- dat[-ttt[,i],]
    
    for(l in 1:5){
      dat.tr <- dat.tr.full[-flds[[l]],]
      dat.tv <- dat.tr.full[flds[[l]],]
      
      j <- 1
      x <- model.matrix(ys ~. , data = dat.tr)
      y1 <- dat.tr[,1]
      
      set.seed(1951)
      lasso.fit <- glmnet(x,y1,alpha=1, lambda = grid) # for lasso
      cv.out <- cv.glmnet(x,y1,alpha=1, nfolds = 5)
      # print(cv.out$lambda.min)
      # plot(lasso.fit)
      set.seed(1951)
      lasso.fit <- glmnet(x,y1,alpha=1, 
                          lambda = cv.out$lambda.min)
      
      
      lasso.pred <- predict(lasso.fit, newx = model.matrix(ys ~. , data = dat.tv))
      # plot(dat.tv[,1], lasso.pred,  main = "lasso")
      # abline(a = 0, b = 1)
      # coef(lasso.fit)
      predictions.tv[flds[[l]] ,j] <- lasso.pred 
      j <- j + 1
      
      #############################################
      
      set.seed(1951)
      en.fit <- glmnet(x,y1,alpha=0.5, lambda = grid) # for lasso
      
      # plot(lasso.fit)
      
      cv.out <- cv.glmnet(x,y1,alpha=0.5, nfolds = 5)
      # print(cv.out$lambda.min)
      
      set.seed(1951)
      en.fit <- glmnet(x,y1,alpha=0.5, 
                       lambda = cv.out$lambda.min)
      
      
      en.pred <- predict(en.fit, newx = model.matrix(ys ~. , data = dat.tv))
      # plot(dat.tv[,1], en.pred,  main = "lasso")
      # abline(a = 0, b = 1)
      # coef(en.fit)
      predictions.tv[flds[[l]] ,j] <- en.pred 
      j <- j + 1
      
      
      # cbind(coef(lasso.fit), coef(en.fit))
      ####################################################
      set.seed(1951)
      pc <- prcomp(dat.tr[, 2:ncol(dat.tr)], scale = TRUE)
      
      pca.lm.fit <- lm(dat.tr[, 1]~pc$x[,1:50])
      # summary(pca.lm.fit)
      # pairs(pc$x[,1:5])
      x.tv <- predict(pc, dat.tv[, 2:ncol(dat.tr)])
      # pca.pred.poor <- as.matrix(cbind(rep(1, 20000), pc$x[,1:5])) %*% as.matrix(pca.lm.fit$coef[1:6])
      pca.pred <- as.matrix(cbind(rep(1, nrow(x.tv)), x.tv[,1:50])) %*% as.matrix(pca.lm.fit$coef[1:51])
      
      predictions.tv[flds[[l]] ,j] <- pca.pred 
      j <- j + 1
      # plot(dat.tv[,1], pca.pred,  main = "PCA")
      # abline(a = 0, b = 1)
      
      #############################################################
      
      set.seed(1951)
      
      save.image(file="temp.RData")
      rm(list=ls())
      load(file="temp.RData")
      pls.fit <- plsr(ys ~., data=dat.tr, scale= TRUE, validation = "none")
      
      bl <- 1:100
      for (i in 1:100) predict(pls.fit, dat.tv[,-1], ncomp = 5)
      plot(1:100, predict(pls.fit, dat.tv[,-1], ncomp = 5))
      pls.pred <- predict(pls.fit, dat.tv[,-1], ncomp = 55)
      calcRMSE(dat.tv[,1], pls.pred)
      plot(dat.tv[,1], pls.pred,  main = "PCA")
      abline(a = 0, b = 1)
      predictions.tv[flds[[l]] ,j] <- pls.pred 
      j <- j+ 1
      
      
      #############################################################
      set.seed(1951)
      lm.fit <- lm(dat.tr[, 1]~ as.matrix(dat.tr[, nonzy+1]) )
      x <- model.matrix(ys ~. , data = dat.tv[, c(1, nonzy+1)])
      lm.pred <- x %*% as.matrix(lm.fit$coef)
      
      
      predictions.tv[flds[[l]] ,j] <- lm.pred
      j <- j+ 1
      
      #############################################################
      set.seed(1951) 
      
      rf.fit <- ranger(ys ~ ., data = dat.tr, importance = "impurity")
      # plot(rf.fit$variable.importance)
      pred.rf <- predict(rf.fit , data = dat.tv)
      
      predictions.tv[flds[[l]] ,j] <- pred.rf$predictions 
      j <- j+1
      
      ####################################################
      set.seed(1951)      
      rf.vs.fit <- ranger(ys ~ ., data = dat.tr, importance = "impurity", regularization.factor = 0.2, regularization.usedepth=FALSE)
      # plot(rf.fit$variable.importance)
      pred.vs.rf <- predict(rf.vs.fit , data = dat.tv)
      
      predictions.tv[flds[[l]] ,j] <- pred.vs.rf$predictions
      j <- j+1
      # plot(dat.tv[,1], pred.rf$predictions, col = as.numeric(dat.tv$class))
      # abline(0,1)
      
      ####################################################
      set.seed(1951)
      bart.fit <- bartMachine(as.data.frame(dat.tr[,-1]), dat.tr[,1], verbose = FALSE)
      # summary(bart.fit)
      pred.bart <- predict(bart.fit, as.data.frame(dat.tv[,-1]))
      predictions.tv[flds[[l]] ,j] <- pred.bart
      j <-j+1
      # plot(dat.tv[,1], pred.bart, col = as.numeric(dat.tv$class), main = "BART")
      # abline(0,1)
      
      ####################################################
      set.seed(1951)        
      Ktrain <- kernelMatrix(kfunc, as.matrix(dat.tr[,2:ncol(dat.tr)]))
      Kvalid <- kernelMatrix(kfunc, as.matrix(dat.tv[,2:ncol(dat.tr)]), as.matrix(dat.tr[,2:ncol(dat.tr)]))
      kpcTrain <- kPCA(Ktrain)
      
      
      # plot the data projection on the principal components
      
      # pairs(cbind(dat.tr[,1], kpcTrain$KPCs[ , 1 : 6]), col = dat.tr$class)
      kpcValid <- predict(kpcTrain, Kvalid)
      # pairs(kpcTest[ , 1 : 6], col = dat.tv$class)
      kpca.lm.fit <- lm(dat.tr[, 1]~  kpcTrain$KPCs[ , 1 : 6] )
      
      kpca.pred <- cbind(rep(1,nrow(kpcValid)), kpcValid[ , 1 : 6]) %*% as.matrix(kpca.lm.fit$coef)
      # plot(dat.tv[,1], kpca.pred , col = as.numeric(dat.tv$class), main = "KPCs")
      # abline(0,1)
      
      predictions.tv[flds[[l]] ,j] <- kpca.pred
      j <- j+ 1
      
      
      ####################################################
      set.seed(1951)
      
      rf.kpc.fit <- ranger(y = dat.tr[, 1], x =  kpcTrain$KPCs[ , 1 : 6], importance = "impurity")
      # plot(rf.fit$variable.importance)
      pred.rf.kpc <- predict(rf.kpc.fit , data = kpcValid[ , 1 : 6])
      predictions.tv[flds[[l]] ,j] <- pred.rf.kpc$predictions
      j <- j+1
      
      
      ####################################################
      set.seed(1951)        
      rf.k.fit <- ranger(y = dat.tr[, 1], x =  Ktrain, importance = "impurity")
      # plot(rf.fit$variable.importance)
      pred.rf.k <- predict(rf.k.fit , data = Kvalid)
      predictions.tv[flds[[l]] ,j] <- pred.rf.k$predictions
      j <- j+1
      
      
      ####################################################
      set.seed(1951)        
      svm.fit <- svm(ys~.,dat.tr)
      
      #Predict using SVM regression
      pred.svm <- predict(svm.fit, dat.tv)
      predictions.tv[flds[[l]] ,j] <- pred.svm
      j <- j+1
      # plot(dat.tv[,1], predYsvm,  main = "RF")
      # abline(0,1)
      
      
      ####################################################
      set.seed(1951)
      ppr.fit <- ppr(dat.tr[,-1], dat.tr[, 1], nterms = 2, max.terms = 5)
      pred.ppr <- ppr.fit %>% predict(dat.tv[,-1])               
      
      predictions.tv[flds[[l]] ,j] <- pred.ppr
      j <- j+1
      # plot(dat.tv[,1], pred.ppr,  main = "RF")
      # abline(0,1)
      ####################################################
      set.seed(1951)        
      brnn.fit <- brnn(as.matrix(dat.tr[, -1]), dat.tr[, 1])
      ###prediction
      pred.brnn <- brnn.fit %>% predict(dat.tv[,-1])                          #predict y hat
      predictions.tv[flds[[l]] ,j] <- pred.brnn
      j <- j+1
      # plot(dat.tv[,1], pred.brnn,  main = "BRNN")
      # abline(0,1)
      
      ####################################################
      set.seed(1951)       
      gbm.fit <- gbm(ys~., data = dat.tr, distribution = "gaussian", cv.folds = 5, n.trees = 200)
      best.iter <- gbm.perf(gbm.fit, method = "cv")
      
      pred.gbm <- gbm.fit %>% predict(dat.tv[,-1],  n.trees = best.iter)               
      predictions.tv[flds[[l]] ,j] <- pred.gbm
      j <- j+1
      # plot(dat.tv[,1], pred.gbm,  main = "GBM")
      # abline(0,1)
      
    }
    
    
    colnames(predictions.tv) <- c("lasso", "elastic net", "pca", "pls",  "LM+15", "RF", "RF+VS", "Bart", "kPCA","RF+kPCA", "RF+kernel", "svm", "ppr", "brnn", "gbm")
    # pairs(cbind(dat.tr.full[,1], predictions.tv)) # correlated predictions
    
    
    ####################################################
    set.seed(1951)    
    lm.fit.SG[[i]][[yind]]<- lm(dat.tr.full[,1]~predictions.tv - 1) 
    ####################################################
    set.seed(1951)    
    nonneg.fit.SG[[i]][[yind]] <- glmnet(predictions.tv,dat.tr.full[,1],alpha=1, lambda = 0, lower.limits = 0, intercept = FALSE) 
    ####################################################
    set.seed(1951) 
    rf.fit.SG[[i]][[yind]] <- ranger(y = dat.tr.full[, 1], x = cbind( predictions.tv, dat.tr.full[, -1]), importance = "impurity")
    ####################################################
    set.seed(1951)    
    cv.out <- cv.glmnet(predictions.tv,dat.tr.full[,1],alpha=1, lambda = grid,  intercept = FALSE)
    # plot(cv.out)
    lasso.fit.SG[[i]][[yind]] <- glmnet(predictions.tv,dat.tr.full[,1],alpha=1, lambda = cv.out$lambda.min,  intercept = FALSE)
    
    # plot(coef(nonneg.fit.SG[[i]][[yind]]), coef(lasso.fit.SG[[i]][[yind]]))
    # sum(coef(lasso.fit.SG[[i]][[yind]]))
    
    # pred.stack.lm <- as.matrix( predictions.tv) %*% as.matrix(lm.fit.SG[[i]][[yind]]$coef)
    
    # pred.stack.nonneg <- predict(nonneg.fit.SG[[i]][[yind]], predictions.te)
    
    
    # pred.stack.rf <- predict(lm.fit.SG[[i]][[yind]], data = cbind(predictions.te, dat.te[, -1]))$prediction
    # 
    #  pred.stack.lasso <- predict(lasso.fit.SG[[i]][[yind]], predictions.te)
    
    # plot(dat.tr.full[, 1], as.matrix( predictions.tv) %*% as.matrix(lm.fit.SG[[i]][[yind]]$coef))
    #   abline(a = 0, b = 1)  
    
    #   plot(dat.tr.full[, 1], predict(nonneg.fit.SG[[i]][[yind]], predictions.tv))
    # abline(a = 0, b = 1)  
    
    #    plot(dat.tr.full[, 1], predict(rf.fit.SG[[i]][[yind]], data = cbind(predictions.tv, dat.tr.full[, -1]))$prediction)
    # abline(a = 0, b = 1)  
    
    # plot(coef(lasso.fit.SG[[i]][[yind]]), coef(nonneg.fit.SG[[i]][[yind]]))
    # plot(dat.tr.full[, 1], predict(lasso.fit.SG[[i]][[yind]], predictions.tv))
    # abline(a = 0, b = 1)
    ##-----------------------------------------------------------
    j <- 1
    ####################################################
    set.seed(1951) 
    x <- model.matrix(ys ~. , data = dat.tr.full)
    y1 <- dat.tr.full[,1]
    lasso.fit <- glmnet(x,y1,alpha=1, lambda = grid) # for lasso
    cv.out <- cv.glmnet(x,y1,alpha=1)
    # print(cv.out$lambda.min)
    
    lasso.fit <- glmnet(x,y1,alpha=1, 
                        lambda = cv.out$lambda.min)
    
    
    lasso.pred <- predict(lasso.fit, newx = model.matrix(ys ~. , data = dat.te))
    RMSE[i,j] <- calcRMSE(dat.te[,1], lasso.pred)
    predictions.te[,j] <- lasso.pred 
    j <- j + 1
    
    #############################################
    
    set.seed(1951)   
    en.fit <- glmnet(x,y1,alpha=0.5, lambda = grid) # for lasso
    
    # plot(lasso.fit)
    
    cv.out <- cv.glmnet(x,y1,alpha=0.5)
    # print(cv.out$lambda.min)
    
    en.fit <- glmnet(x,y1,alpha=0.5, 
                     lambda = cv.out$lambda.min)
    
    
    en.pred <- predict(en.fit, newx = model.matrix(ys ~. , data = dat.te))
    RMSE[i,j] <- calcRMSE(dat.te[,1], en.pred)
    predictions.te[,j] <- en.pred 
    j <- j + 1
    
    ####################################################
    
    set.seed(1951)    
    pc <- prcomp(dat.tr.full[, 2:ncol(dat.tr.full)], scale = TRUE)
    
    pca.lm.fit<- lm(dat.tr.full[, 1]~pc$x[,1:5])
    
    x.te <- predict(pc, dat.te[, 2:ncol(dat.tr.full)])
    # pca.pred.poor <- as.matrix(cbind(rep(1, 20000), pc$x[,1:5])) %*% as.matrix(pca.lm.fit$coef[1:6])
    pca.pred <- as.matrix(cbind(rep(1, nrow(x.te)), x.te[,1:5])) %*% as.matrix(pca.lm.fit$coef[1:6])
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pca.pred)
    predictions.te[,j] <- pca.pred
    j <- j + 1
    # plot(dat.te[,1], pca.pred, col = as.numeric(dat.te$class), main = "PCA")
    # abline(a = 0, b = 1)
    
    #############################################################
    
    set.seed(1951)   
    
    pls.fit <- plsr(ys ~., data=dat.tr.full, scale= TRUE, validation = "CV")
    
    
    pls.pred <- predict(pls.fit, dat.te[,-1], ncomp = pls.pc[yind])
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pls.pred)
    predictions.te[,j] <- pls.pred 
    j <- j+ 1
    
    
    #############################################################
    set.seed(1951)     
    lm.fit <- lm(dat.tr.full[, 1]~ as.matrix(dat.tr.full[, nonzy+1]) )
    x <- model.matrix(ys ~. , data = dat.te[, c(1, nonzy+1)])
    lm.pred <- x %*% as.matrix(lm.fit$coef)
    
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], lm.pred)
    predictions.te[,j] <- lm.pred 
    j <- j+ 1
    
    #############################################################
    
    
    set.seed(1951)   
    rf.fit <- ranger(ys ~ ., data = dat.tr.full, importance = "impurity")
    # plot(rf.fit$variable.importance)
    pred.rf <- predict(rf.fit , data = dat.te)
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pred.rf$predictions)
    predictions.te[,j] <- pred.rf$predictions 
    j <- j+1
    #############################################################
    
    
    set.seed(1951)   
    
    
    rf.vs.fit <- ranger(ys ~ ., data = dat.tr.full, importance = "impurity",regularization.factor = 0.2, regularization.usedepth=FALSE)
    # plot(rf.fit$variable.importance)
    pred.rf.vs <- predict(rf.vs.fit , data = dat.te)
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pred.rf.vs$predictions)
    predictions.te[,j] <- pred.rf.vs$predictions
    j <- j+1
    # plot(dat.te[,1], pred.rf$predictions, col = as.numeric(dat.te$class))
    # abline(0,1)
    
    #############################################################
    
    
    set.seed(1951)   
    
    
    bart.fit <- bartMachine(as.data.frame(dat.tr.full[,-1]), dat.tr.full[,1], verbose = FALSE)
    # summary(bart.fit)
    pred.bart <- predict(bart.fit, as.data.frame(dat.te[,-1]))
    RMSE[i,j] <- calcRMSE(dat.te[,1], pred.bart)
    predictions.te[,j] <- pred.bart
    j <- j+1
    # plot(dat.te[,1], pred.bart, col = as.numeric(dat.te$class), main = "BART")
    # abline(0,1)
    
    #############################################################
    
    
    set.seed(1951)   
    
    Ktrain <- kernelMatrix(kfunc, as.matrix(dat.tr.full[,2:ncol(dat.tr.full)]))
    Ktest <- kernelMatrix(kfunc, as.matrix(dat.te[,2:ncol(dat.tr.full)]), as.matrix(dat.tr.full[,2:ncol(dat.tr.full)]))
    kpcTrain <- kPCA(Ktrain)
    
    
    # plot the data projection on the principal components
    
    # pairs(cbind(dat.tr.full[,1], kpcTrain$KPCs[ , 1 : 6]), col = dat.tr.full$class)
    kpcTest <- predict(kpcTrain, Ktest)
    # pairs(kpcTest[ , 1 : 6], col = dat.te$class)
    kpca.lm.fit <- lm(dat.tr.full[, 1]~  kpcTrain$KPCs[ , 1 : 6] )
    
    kpca.pred <- cbind(rep(1,nrow(kpcTest)), kpcTest[ , 1 : 6]) %*% as.matrix(kpca.lm.fit$coef)
    # plot(dat.te[,1], kpca.pred , col = as.numeric(dat.te$class), main = "KPCs")
    # abline(0,1)
    
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], kpca.pred)
    predictions.te[,j] <- kpca.pred
    j <- j+ 1
    
    
    
    #############################################################
    set.seed(1951)   
    
    rf.kpca.fit <- ranger(y = dat.tr.full[, 1], x =  kpcTrain$KPCs[ , 1 : 6], importance = "impurity")
    # plot(rf.fit$variable.importance)
    pred.rf.kpca <- predict(rf.kpca.fit , data = kpcTest[ , 1 : 6])
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pred.rf.kpca$predictions)
    predictions.te[,j] <- pred.rf.kpca$predictions
    j <- j+1
    
    #############################################################
    set.seed(1951)    
    
    rf.k.fit <- ranger(y = dat.tr.full[, 1], x =  Ktrain, importance = "impurity")
    # plot(rf.fit$variable.importance)
    pred.rf.k <- predict(rf.k.fit , data = Ktest)
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pred.rf.k$predictions)
    predictions.te[,j] <- pred.rf.k$predictions
    j <- j+1
    
    
    #############################################################
    set.seed(1951) 
    svm.fit <- svm(ys~.,dat.tr.full)
    
    #Predict using SVM regression
    pred.svm <- predict(svm.fit, dat.te)
    RMSE[i,j] <- calcRMSE(dat.te[,1], pred.svm)
    predictions.te[,j] <- pred.svm
    j <- j+1
    # plot(dat.te[,1], pred.svm,  main = "RF")
    # abline(0,1)
    
    
    
    #############################################################
    set.seed(1951) 
    ppr.fit <- ppr(dat.tr.full[, -1], dat.tr.full[, 1], nterms = 2, max.terms = 5)
    pred.ppr <- ppr.fit %>% predict(dat.te[, -1])              
    RMSE[i,j] <- calcRMSE(dat.te[,1], pred.ppr)
    predictions.te[,j] <- pred.ppr
    j <- j + 1
    
    #############################################################
    set.seed(1951)    
    brnn.fit <- brnn(as.matrix(dat.tr.full[,-1]), dat.tr.full[, 1])
    pred.brnn <- brnn.fit %>% predict(dat.te[, -1])                        
    RMSE[i,j] <- calcRMSE(dat.te[,1], pred.brnn)
    predictions.te[,j] <- pred.brnn
    j <- j + 1
    # plot(dat.tv[,1], pred.brnn,  main = "BRNN")
    # abline(0,1)
    
    #############################################################
    set.seed(1951)     
    gbm.fit <- gbm(
      ys ~ .,
      data = dat.tr.full,
      distribution = "gaussian",
      cv.folds = 5,
      n.trees = 200)
    
    best.iter <- gbm.perf(gbm.fit, method = "cv")
    
    pred.gbm <- gbm.fit %>% predict(dat.te[, -1],  n.trees = best.iter)                  
    RMSE[i,j] <- calcRMSE(dat.te[,1], pred.gbm)
    predictions.te[,j] <- pred.gbm
    j <- j + 1
    # plot(dat.tv[,1], pred.gbm,  main = "GBM")
    # abline(0,1)
    
    RMSE[i,j] <-  calcRMSE(dat.te[,1],   rowMeans(predictions.te[, 1:15]))
    
    
    colnames(predictions.te) <- c("lasso", "elastic net", "pca", "pls",  "LM+15", "RF", "RF+VS", "Bart", "kPCA","RF+kPCA", "RF+kernel", "svm", "ppr", "brnn", "gbm")
    # pairs(predictions.te)
    
    
    pred.stack.lm <- as.matrix(predictions.te) %*% as.matrix(lm.fit.SG[[i]][[yind]]$coef)
    
    pred.stack.nonneg <- predict(nonneg.fit.SG[[i]][[yind]], predictions.te)
    
    
    pred.stack.rf <- predict(rf.fit.SG[[i]][[yind]], data = cbind(predictions.te, dat.te[, -1]))$prediction
    
    pred.stack.lasso <- predict(lasso.fit.SG[[i]][[yind]], predictions.te)
    # plot(pred.stack.nonneg, pred.stack.lm)
    
    
    RMSE[i,j+1] <- calcRMSE(dat.te[,1],  pred.stack.lm)
    RMSE[i,j+2] <- calcRMSE(dat.te[,1],  pred.stack.nonneg)
    RMSE[i,j+3] <- calcRMSE(dat.te[,1],  pred.stack.rf)
    RMSE[i,j+4] <- calcRMSE(dat.te[,1],  pred.stack.lasso)
  }
  
  colnames(RMSE) <- c("lasso", "elastic net", "pca", "pls", "LM+15", "RF", "RF+VS", "Bart", "kPCA","RF+kPCA", "RF+kernel", "svm", "ppr", "brnn", "gbm", "ensamble MA",  "ensLM",  "ensNonNeg", "ensRF", "enslasso")
  
  RMSEFULL[[yind]] <- RMSE
  
  RMSE1 <- RMSE %>% as.data.frame() %>% mutate(split = as.factor(1:N)) %>% pivot_longer(1:20,"ALGORITHM")
  # RMSE1 %>% ggplot(aes(x = ALGORITHM, y = value)) + geom_boxplot()
  
  p1 <- RMSE1 %>% group_by(ALGORITHM) %>% mutate(meanv = mean(value)) %>% ggplot(aes(x= ALGORITHM, y = value, color = split, group = split)) +
    geom_point()+ geom_line() +geom_line(aes(x = ALGORITHM, y = meanv), col = 1) + theme(legend.position = "none")+
    ylab("RMSE")+ ggtitle(names(ys_tr)[yind])+
    theme(axis.text.x=element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  
  print(p1)
  
  t1[[yind]] <- RMSE1 %>% group_by(ALGORITHM) %>% summarise(meanv = mean(value), sdv=sd(value)) %>% arrange(desc(meanv))
  
}ß

