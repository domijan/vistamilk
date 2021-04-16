set.seed(1979)
library(tidyverse)
# library(rgl)
library(dendextend)
library(glmnet) # for the ridge regression
library(pls)
library(ranger)
library(bartMachine)
library(e1071)

tr <- read.csv("tr.csv") 
te <- read.csv("te.csv") 

# tr <- tr %>% na.omit()
calcRMSE <- function(y, yhat){sqrt(mean((yhat - y)^2))}
# dim(tr)
# dim(te)

y <- tr[,1:3]
tr <- tr[, -c(1:3)]


a <- apply(tr,2,mean)
b <- apply(tr,2,sd)

pairs(log(y)) # all non-normal
# points3d(y[,1], y[,2], y[,3])
matplot(t(tr), type = "l", col = 9, lty = 1)
matplot(t(te), type = "l", col = 9, lty = 1)


# pairs(tr[,1:10])

# 
# plot(cor(tr)[1,], col = 0)
# lines(cor(tr)[1,], col = 9)
# lines(cor(tr)[2,], col = 9)
# matplot(cor(te), type = "l", col = 9, lty = 1)
# 
# 
# library(corrplot)
# corrplot(cor(te), method="color")

# library(corrr)
# 
# tr %>% correlate() %>% network_plot(min_cor=0.6)

tr <- scale(tr,a,b)
te <- scale(te,a,b)
matplot(t(tr), type = "l", col = 9, lty = 1)
matplot(t(te), type = "l", col = 9, lty = 1)



############################################
# Look for clusters in the training data


D <- as.matrix(dist(tr,method="euclidean"))
image(D)
pairs(log(y))

h<- hclust(as.dist(D), method = "ward.D2")
d <- as.dendrogram(h)
plot(d)

heatmap(D,
        Colv=as.dendrogram(h),     
        Rowv=as.dendrogram(h))


d2 <- color_branches(d,k=7, col=c(2,3,5,4, 6,7,8)) 
# auto-coloring 4 clusters of branches.
# plot(d2)

# plot(reorder(d2, D, method = "OLO")) 

labl.Ave <- cutree(d2,7)
matplot(t(tr), type = "l", col = labl.Ave, lty = 1) # see the 7 clusters


#############################################################

library(BKPC)
library(kernlab)
kfunc <-  rbfdot(sigma = 0.0005)
# kfunc <- laplacedot(sigma = 0.0001)
# kfunc <- anovadot(sigma = 1, degree = 1)
Ktrain2 <- kernelMatrix(kfunc, tr)
image(Ktrain2)

Ktest2 <- kernelMatrix(kfunc, te, tr)



kpcKern2 <- kPCA(Ktrain2)


# plot the data projection on the principal components

pairs(kpcKern2$KPCs[ , 1 : 6], col = labl.Ave)
KPCpred3 <- predict(kpcKern2, Ktest2)
pairs(KPCpred3[ , 1 : 6])
#############################################################

# correlation with responses

indexcor <- abs(cor(cbind(y,tr), use = "pairwise.complete.obs"))[1:3, -c(1:3)]
matplot(t(indexcor), type = "l", col = 9, lty = 1)

matplot(t(tr), type = "l", col = labl.Ave, lty = 1)
lines(indexcor[1,]*20,lty = 3)


# clustering the columns (variables) rather than observations
D2 <- as.matrix(1-abs(cor(tr)))

# image(D)


h2<- hclust(as.dist(D2), method = "ward.D2")
dnew <- as.dendrogram(h2)
plot(dnew)

heatmap(D2,
        Colv=as.dendrogram(h2),     
        Rowv=as.dendrogram(h2))


d2 <- color_branches(dnew,k=7, col = 2:8) 
# auto-coloring 4 clusters of branches.
# plot(d2)

plot(reorder(d2, D2, method = "OLO")) 

labl.cols <- cutree(d2,15) # 15 variables
matplot(t(tr), type = "l", col = labl.Ave, lty = 1)
points(1:1060, rep(-5.8, 1060), col = labl.cols)

# Pick 15 top features
nonz1 <- NULL
nonz2 <- NULL
nonz3 <- NULL
for (j in unique(labl.cols)){
  nonz1[j] <- (which(indexcor[1,]==max(indexcor[1,labl.cols==j])))
  nonz2[j] <- (which(indexcor[2,]==max(indexcor[2,labl.cols==j])))
  nonz3[j] <- (which(indexcor[3,]==max(indexcor[3,labl.cols==j])))
  
  
  }

matplot(t(tr), type = "l", col = labl.Ave, lty = 1)
abline(v = nonz1, col = 1:15) # pick values in a sensible region
points(1:1060, rep(-5.8, 1060), col = labl.cols)
# max(indexcor[1,])

pairs(cbind(y[,1], tr[,nonz1]), col = labl.Ave)

#############################################################
#############################################################

# First response:

dat <- cbind((y[,1]), tr)
dat <- as.data.frame(dat)
dat <- dat %>% mutate(class = as.factor(labl.Ave))
dat %>% ggplot(aes(x= V1)) + geom_histogram()
dat <- dat %>% na.omit()
nonz <- nonz1
# dat <- dat %>% filter(V1<15)
dim(dat)

N <- 50
ttt <- replicate(N, sample(nrow(dat), 200))




RMSE <- matrix(0, N, 15)
for (i in 1:N){
  j <- 1
  dat.tr <- dat[ttt[,i],]
  dat.te <- dat[-ttt[,i],]
  
  # pairs(dat[,1:4])
  # dat %>% ggplot(aes(x= (V1)))+geom_histogram()
  # dat %>% ggplot(aes(x= log(V1)))+geom_histogram()
  # qqnorm(y=(dat$V1[dat$V1<15]))
  # qqnorm(y=log(dat$V1))
  
  x <- model.matrix(V1 ~. , data = dat.tr[, -1062])
  
  y1 <- dat.tr[,1]
  
  grid <- 10^seq(-3, 3, length = 100)
  
  
  
  
  lasso.fit <- glmnet(x,y1,alpha=1, lambda = grid) # for lasso
  
  # plot(lasso.fit)
  
  cv.out <- cv.glmnet(x,y1,alpha=1)
  print(cv.out$lambda.min)
  
  lasso.fit <- glmnet(x,y1,alpha=1, 
                      lambda = cv.out$lambda.min)
  # head(coef(lasso.fit))
  nonzero <- which(coef(lasso.fit)[3:1062]!=0) #intercept there twice?!
  
  # matplot(t(tr), type = "l", col = labl.Ave, lty = 1)
  # abline(v = nonzero) # pick values in a sensible region
  
  
  
  # pairs(cbind(y1, dat.tr[,1+nonzero])) # still perfectly correlated predictors
  
  lasso.pred <- predict(lasso.fit, newx = model.matrix(V1 ~. , data = dat.te[, -1062]))
  RMSE[i,j] <- calcRMSE(dat.te[,1], lasso.pred)
  j <- j + 1
  # plot(dat.te[,1], lasso.pred, col = as.numeric(dat.te$class), main = "lasso")
  # abline(a = 0, b = 1)
  # plot(lasso.fit)
  #############################################################
  
  # pc <- prcomp(tr[!is.na(y[,1])&y[,1]<15,], scale = TRUE)
  # pairs(cbind(y1, pc$x[,1:9]), col = as.numeric(dat$class))
  
  # dim(pc$x)
  # summary(lm(dat[, 1]~pc$x[,1:10]))
  # library(car)
  # avPlots(lm(dat[, 1]~pc$x[,1:10]))
  # 
  # 
  # lmfit<- lm(y[!is.na(y[,1])&y[,1]<15,1]~pc$x[,1:397])
  # plot(lmfit,1)
  # pca.pred <- predict(lmfit)
  # calcRMSE(dat[,1], pca.pred)
  # plot(dat[,1], pca.pred, col = as.numeric(dat$class))
  
  
  
  
  pc <- prcomp(dat.tr[, 2:1061], scale = TRUE)
  # pairs(cbind(dat.tr[,1], pc$x[,1:9]), col = as.numeric(dat.tr$class))
  # plot(pc)
  # dim(pc$x)
  # summary(lm(dat.tr[, 1]~pc$x[,1:10]))
  # library(car)
  # avPlots(lm(dat.tr[, 1]~pc$x[,1:10]))
  
  
  lmfit<- lm(dat.tr[, 1]~pc$x[,1:5])
  # plot(lmfit,1)
  # anova(lmfit)
  # summary(lmfit)
  x.te <- predict(pc, dat.te[, 2:1061])
  # pca.pred.poor <- as.matrix(cbind(rep(1, 200), pc$x[,1:5])) %*% as.matrix(lmfit$coef[1:6])
  pca.pred.poor.te <- as.matrix(cbind(rep(1, nrow(x.te)), x.te[,1:5])) %*% as.matrix(lmfit$coef[1:6])
  # plot(dat.tr[,1], pca.pred.poor, col = as.numeric(dat.tr$class))
  # abline(0,1)
  # pca.pred <- predict(lmfit)
  # calcRMSE(dat.tr[,1], pca.pred.poor)
  RMSE[i,j] <- calcRMSE(dat.te[,1], pca.pred.poor.te)
  j <- j + 1
  # plot(dat.te[,1], pca.pred.poor.te, col = as.numeric(dat.te$class), main = "PCA")
  # abline(a = 0, b = 1)
  
  #############################################################
  
  # # library(pls)
  # # pcr.fit <- pcr(y[!is.na(y[,1])&y[,1]<15,1]~tr[!is.na(y[,1])&y[,1]<15,])
  # pcr.fit <- pcr(V1~., data=dat.tr, scale= TRUE, validation = "CV")
  # 
  # summary(pcr.fit)
  # validationplot(pcr.fit, val.type = "MSEP")
  # pcr.pred <- predict(pcr.fit, dat.te[,-1], ncomp = 5)
  # # pcr.pred <- predict(pcr.fit, tr[!is.na(y[,1])&y[,1]<15,], ncomp = 356)
  # # calcRMSE(y[!is.na(y[,1])&y[,1]<15,1], pcr.pred)
  # # plot(y[!is.na(y[,1])&y[,1]<15,1], pcr.pred,  col = as.numeric(dat$class))
  # 
  # calcRMSE(dat.te[,1], pcr.pred)
  # plot(dat.te[,1], pcr.pred,  col = as.numeric(dat.te$class))
  # abline(a = 0, b = 1)
  #############################################################
  
  pls.fit <- plsr(V1 ~., data=dat.tr[, -1062], scale= TRUE, validation = "CV")
  # summary(pls.fit)
  # validationplot(pls.fit, val.type = "MSEP")
  
  pls.pred <- predict(pls.fit, dat.te[,-1], ncomp = 3)
  
  RMSE[i,j] <- calcRMSE(dat.te[,1], pls.pred)
  j <- j+ 1
  
  pls.pred <- predict(pls.fit, dat.te[,-1], ncomp = 4)
  
  RMSE[i,j] <- calcRMSE(dat.te[,1], pls.pred)
  j <- j+ 1
  
  pls.pred <- predict(pls.fit, dat.te[,-1], ncomp = 5)
  
  RMSE[i,j] <- calcRMSE(dat.te[,1], pls.pred)
  j <- j+ 1
  # plot(dat.te[,1],pls.pred,  col = as.numeric(dat.te$class), main = "PLS")
  # abline(a = 0, b = 1)
  
  # summary(pls.fit)
  # validationplot(pls.fit, val.type = "MSEP")
  pls.pred <- predict(pls.fit, dat.te[,-1], ncomp = 6)
  
  RMSE[i,j] <- calcRMSE(dat.te[,1], pls.pred)
  j <- j+ 1
 
  pls.pred <- predict(pls.fit, dat.te[,-1], ncomp = 7)
  
  RMSE[i,j] <- calcRMSE(dat.te[,1], pls.pred)
  j <- j+ 1
  pls.pred <- predict(pls.fit, dat.te[,-1], ncomp = 9)
  
  RMSE[i,j] <- calcRMSE(dat.te[,1], pls.pred)
  j <- j+ 1

  #############################################################
  
  lmfit <- lm(dat.tr[, 1]~ as.matrix(dat.tr[, nonz+1]) )
  # plot(lmfit,1)
  # anova(lmfit)
  # summary(lmfit)
  
  x <- model.matrix(V1 ~. , data = dat.te[, c(1, nonz+1)])
  
  
  
  pca.pred.poor.te <- x %*% as.matrix(lmfit$coef)
  # plot(dat.te[,1], pca.pred.poor.te , col = as.numeric(dat.te$class), main = "Uncorrelated variables")
  # abline(0,1)
  
  
  RMSE[i,j] <- calcRMSE(dat.te[,1], pca.pred.poor.te)
  j <- j+ 1
  # lmfit <- lm(V1~ wave_530+ wave_170+wave_134+wave_141+wave_162+
  #             wave_193+wave_221+wave_719+wave_493+wave_596+
  #             wave_602+wave_656+wave_662+wave_688+wave_713, data = dat.tr)
  # library(condvis2)
  # condvis(dat.tr, lmfit,  response = "V1", sectionvars="wave_530", 
  #         conditionvars=c("wave_170", "wave_134", "wave_141", "wave_162",
  #                         "wave_193", "wave_221", "wave_719", "wave_493", "wave_596",
  #                         "wave_602", "wave_656", "wave_662", "wave_688", "wave_713"), pointColor= "class")
  
  #############################################################
  
  
  
  rf.fit <- ranger(V1 ~ ., data = dat.tr, importance = "impurity")
  # plot(rf.fit$variable.importance)
  pred.rf <- predict(rf.fit , data = dat.te)
  
  RMSE[i,j] <- calcRMSE(dat.te[,1], pred.rf$predictions)
  j <- j+1
  # plot(dat.te[,1], pred.rf$predictions, col = as.numeric(dat.te$class), main = "RF")
  # abline(0,1)
  
  
  rf.fit <- ranger(V1 ~ ., data = dat.tr, importance = "impurity",regularization.factor = 0.2, regularization.usedepth=FALSE)
  # plot(rf.fit$variable.importance)
  pred.rf <- predict(rf.fit , data = dat.te)
  
  RMSE[i,j] <- calcRMSE(dat.te[,1], pred.rf$predictions)
  j <- j+1
  # plot(dat.te[,1], pred.rf$predictions, col = as.numeric(dat.te$class))
  # abline(0,1)
  
  
  
  bart_machine = bartMachine(as.data.frame(dat.tr[,-1]), dat.tr[,1])
  # summary(bart_machine)
  pred.bart <- predict(bart_machine, as.data.frame(dat.te[,-1]))
  RMSE[i,j] <- calcRMSE(dat.te[,1], pred.bart)
  j <-j+1
  # plot(dat.te[,1], pred.bart, col = as.numeric(dat.te$class), main = "BART")
  # abline(0,1)
  
  
  Ktrain2 <- kernelMatrix(kfunc, as.matrix(dat.tr[,2:1061]))
  Ktest2 <- kernelMatrix(kfunc, as.matrix(dat.te[,2:1061]), as.matrix(dat.tr[,2:1061]))
 kpcKern2 <- kPCA(Ktrain2)
  
  
  # plot the data projection on the principal components
  
  # pairs(cbind(dat.tr[,1], kpcKern2$KPCs[ , 1 : 6]), col = dat.tr$class)
  KPCpred3 <- predict(kpcKern2, Ktest2)
  # pairs(KPCpred3[ , 1 : 6], col = dat.te$class)
  lmfit <- lm(dat.tr[, 1]~  kpcKern2$KPCs[ , 1 : 6] )
  # plot(lmfit,1)
  # anova(lmfit)
  # summary(lmfit)
  # predict(lmfit, KPCpred3[ , 1 : 6])
 
  
  
  
  pca.pred.poor.te <- cbind(rep(1,nrow(KPCpred3)), KPCpred3[ , 1 : 6]) %*% as.matrix(lmfit$coef)
  # plot(dat.te[,1], pca.pred.poor.te , col = as.numeric(dat.te$class), main = "KPCs")
  # abline(0,1)
  
  
  RMSE[i,j] <- calcRMSE(dat.te[,1], pca.pred.poor.te)
  j <- j+ 1
  
  
  rf.fit <- ranger(y = dat.tr[, 1], x =  Ktrain2, importance = "impurity")
  # plot(rf.fit$variable.importance)
  pred.rf <- predict(rf.fit , data = Ktest2)
  
  RMSE[i,j] <- calcRMSE(dat.te[,1], pred.rf$predictions)
  j <- j+1
  # plot(dat.te[,1], pred.rf$predictions, col = as.numeric(dat.te$class), main = "RF")
  # abline(0,1)


  modelsvm <- svm(V1~.,dat.tr[, -1062])
  
  #Predict using SVM regression
  predYsvm = predict(modelsvm, dat.te[, -1062])
  RMSE[i,j] <- calcRMSE(dat.te[,1], predYsvm)
  j <- j+1
  # plot(dat.te[,1], predYsvm, col = as.numeric(dat.te$class), main = "RF")
  # abline(0,1)

}
i
j
colnames(RMSE) <- c("lasso", "pca", "pls3", "pls4","pls5", "pls6", "pls7", "pls9", "uncorr", "RF", "varRF", "Bart", "kpca", "kernelRF", "svm")
RMSE1 <- RMSE  %>%as.data.frame()%>% mutate(split = as.factor(1:N))%>% pivot_longer(1:15,"ALGORITHM")
# RMSE1 %>% ggplot(aes(x = ALGORITHM, y = value)) + geom_boxplot()
RMSE1 %>% group_by(ALGORITHM) %>% mutate(meanv = mean(value)) %>% ggplot(aes(x= ALGORITHM, y = value, color = split, group = split))+
  geom_point()+ geom_line() +geom_line(aes(x = ALGORITHM, y = meanv), col = 1)+ theme(legend.position = "none")
