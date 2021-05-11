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
a <- read.csv("recleaneddataset1.csv") 

y_names <- c("kappa_casein",
"alpha_s2_casein",
"alpha_s1_casein",
"beta_casein",
"alpha_lactalbumin",
"beta_lactoglobulin_a",
"beta_lactoglobulin_b",
"Casein_micelle_size",
"Heat_stability",
"Native_pH",
"RCT",
"k20",
"a30",
"a60")

ys <- a %>% select(y_names)
x <- a[,60:590] 
matplot(t(x), type = "l", col = as.numeric(as.factor(a$Breed))+1)
matplot(t(x), type = "l", col = a$Farm+1)
# matplot(t(log(1/tr[,4:1060],10)), type = "l")
calcRMSE <- function(y, yhat){sqrt(mean((yhat - y)^2))}

#names(a)
id_te <- sample(nrow(x), 69)
ys_te <- ys[id_te,]
te <- x[id_te, ]
ys_tr <- ys[-id_te,]
tr <- x[-id_te, ]

me <- apply(tr,2,mean)
sd <- apply(tr,2,sd)

# pairs(log(ys)) # all non-normal
ys %>% pivot_longer(1:14, "key", "value") %>% ggplot(aes(x = log(value+1)))+ 
  geom_histogram()+facet_wrap(vars(key), scales = "free")

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

tr <- scale(tr,me,sd)
te <- scale(te,me,sd)
matplot(t(tr), type = "l", col = a$Farm[-id_te]+1, lty = 1)
matplot(t(te), type = "l", col = as.numeric(as.factor(a$Breed[-id_te]))+1, lty = 1)



############################################
# Look for clusters in the training data


D <- as.matrix(dist(tr,method="euclidean"))
# image(D)
pairs(log(y))

h <- hclust(as.dist(D), method = "ward.D2")
d <- as.dendrogram(h)
plot(d)

heatmap(D,
        Colv=as.dendrogram(h),     
        Rowv=as.dendrogram(h))


d2 <- color_branches(d,k=7, col=c(2,3,5,4, 6,7,8)) 
# auto-coloring 4 clusters of branches.
plot(d2)

# plot(reorder(d2, D, method = "OLO")) 

labl.Ave <- cutree(d2,7)
matplot(t(tr), type = "l", col = labl.Ave, lty = 1, ylab = "training data") # see the 7 clusters

table(a$Farm[-id_te], labl.Ave)
table(a$Breed[-id_te], labl.Ave)

#############################################################

kfunc <-  rbfdot(sigma = 0.0005)
# kfunc <- laplacedot(sigma = 0.0001)
# kfunc <- anovadot(sigma = 1, degree = 1)
Ktrain2 <- kernelMatrix(kfunc, tr)
image(Ktrain2)

Ktest2 <- kernelMatrix(kfunc, te, tr)



kpcKern2 <- kPCA(Ktrain2)


# plot the data projection on the principal components

pairs(kpcKern2$KPCs[ , 1 : 6], col = labl.Ave)

newdat <- cbind(log(ys_tr), kpcKern2$KPCs[ , 1 : 6], labl.Ave)

newdat <- newdat %>% as.data.frame() %>% na.omit() %>% mutate(labl.Ave = as.factor(labl.Ave))
pairs(newdat[,1:20],  col = newdat[,21])
KPCpred3 <- predict(kpcKern2, Ktest2)
pairs(KPCpred3[ , 1 : 6])
#############################################################

# correlation with responses

indexcor <- abs(cor(cbind(ys_tr,tr), use = "pairwise.complete.obs"))[1:14, -c(1:14)]
matplot(t(indexcor), type = "l", col = 9, lty = 1)

# matplot(t(tr), type = "l", col = labl.Ave, lty = 1)
# lines(indexcor[1,]*5,lty = 3)


# clustering the columns (variables) rather than observations
D2 <- as.matrix(1-abs(cor(tr)))

image(D2)


h2<- hclust(as.dist(D2), method = "ward.D2")
dnew <- as.dendrogram(h2)
plot(dnew)
# library(pheatmap)
# pheatmap(tr, cluster_cols = TRUE, clustering_distance_cols = "correlation",  cutree_cols  = 15)
heatmap(D2,
        Colv=as.dendrogram(h2),     
        Rowv=as.dendrogram(h2))


d2 <- color_branches(dnew,k=15, col = 2:16) 
# auto-coloring 4 clusters of branches.
 plot(d2)

plot(reorder(d2, D2, method = "OLO")) 

labl.cols <- cutree(d2,15) # 15 variables
matplot(t(tr), type = "l", col = labl.Ave, lty = 1)
points(1:1060, rep(-5.8, 1060), col = labl.cols)

image(t(indexcor))
# Pick 15 top features
nonz <- matrix(0,14,15)


for (j in unique(labl.cols)){
  for (i in 1:14){nonz[i, j] <- (which(indexcor[i,]==max(indexcor[i,labl.cols==j])))}
  }
  
  

matplot(t(tr), type = "l", col = labl.Ave, lty = 1)
abline(v = nonz1, col = 1:15) # pick values in a sensible region

# max(indexcor[1,])

pairs(cbind(ys_tr[,1], tr[,nonz[1,]]), col = labl.Ave)

#############################################################
#############################################################
# Each response: to pick the model

N <- 50
grid <- 10^seq(-3, 3, length = 100)
RMSEFULL <- vector(mode = "list", length = 14)
t1 <- vector(mode = "list", length = 14)

for(yind in 1:14){
  if (yind ==8) dat <- cbind(log(ys_tr[,yind]), tr)else dat <- cbind(ys_tr[,yind], tr)
  dat <- as.data.frame(dat)
  
  dat %>% ggplot(aes(x= V1)) + geom_histogram()
  dat <- dat %>% na.omit()
  
  nonzy <- nonz[yind,]
  ttt <- replicate(N, sample(nrow(dat), 200))
 
  
  RMSE <- matrix(0, N, 14)
  predictions <- matrix(0,nrow(dat)-200, 13)
  for (i in 1:N){
    j <- 1
    dat.tr <- dat[ttt[,i],]
    dat.te <- dat[-ttt[,i],]

    x <- model.matrix(V1 ~. , data = dat.tr)
    
    y1 <- dat.tr[,1]
    
  
    
    lasso.fit <- glmnet(x,y1,alpha=1, lambda = grid) # for lasso
    
    # plot(lasso.fit)
    
    cv.out <- cv.glmnet(x,y1,alpha=1)
    # print(cv.out$lambda.min)
    
    lasso.fit <- glmnet(x,y1,alpha=1, 
                        lambda = cv.out$lambda.min)
    
    
    lasso.pred <- predict(lasso.fit, newx = model.matrix(V1 ~. , data = dat.te))
    RMSE[i,j] <- calcRMSE(dat.te[,1], lasso.pred)
    predictions[,j] <- lasso.pred 
    j <- j + 1
    
    pc <- prcomp(dat.tr[, 2:ncol(dat.tr)], scale = TRUE)
    
    lmfit<- lm(dat.tr[, 1]~pc$x[,1:5])

    x.te <- predict(pc, dat.te[, 2:ncol(dat.tr)])
    # pca.pred.poor <- as.matrix(cbind(rep(1, 200), pc$x[,1:5])) %*% as.matrix(lmfit$coef[1:6])
    pca.pred.poor.te <- as.matrix(cbind(rep(1, nrow(x.te)), x.te[,1:5])) %*% as.matrix(lmfit$coef[1:6])

    RMSE[i,j] <- calcRMSE(dat.te[,1], pca.pred.poor.te)
    predictions[,j] <- pca.pred.poor.te 
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
    
    pls.fit <- plsr(V1 ~., data=dat.tr, scale= TRUE, validation = "CV")
    # summary(pls.fit)
    # validationplot(pls.fit, val.type = "MSEP")
    
    pls.pred <- predict(pls.fit, dat.te[,-1], ncomp = 3)
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pls.pred)
    predictions[,j] <- pls.pred 
    j <- j+ 1
    
    pls.pred <- predict(pls.fit, dat.te[,-1], ncomp = 4)
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pls.pred)
    predictions[,j] <- pls.pred 
    j <- j+ 1
    
    pls.pred <- predict(pls.fit, dat.te[,-1], ncomp = 5)
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pls.pred)
    predictions[,j] <- pls.pred 
    j <- j+ 1
  
    #############################################################
    
    lmfit <- lm(dat.tr[, 1]~ as.matrix(dat.tr[, nonzy+1]) )
    x <- model.matrix(V1 ~. , data = dat.te[, c(1, nonzy+1)])
    pca.pred.poor.te <- x %*% as.matrix(lmfit$coef)
    
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pca.pred.poor.te)
    predictions[,j] <- pca.pred.poor.te 
    j <- j+ 1
  
    #############################################################
    
    
    
    rf.fit <- ranger(V1 ~ ., data = dat.tr, importance = "impurity")
    # plot(rf.fit$variable.importance)
    pred.rf <- predict(rf.fit , data = dat.te)
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pred.rf$predictions)
    predictions[,j] <- pred.rf$predictions 
    j <- j+1
  
    
    rf.fit <- ranger(V1 ~ ., data = dat.tr, importance = "impurity",regularization.factor = 0.2, regularization.usedepth=FALSE)
    # plot(rf.fit$variable.importance)
    pred.rf <- predict(rf.fit , data = dat.te)
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pred.rf$predictions)
    predictions[,j] <- pred.rf$predictions
    j <- j+1
    # plot(dat.te[,1], pred.rf$predictions, col = as.numeric(dat.te$class))
    # abline(0,1)
    
    
    
    bart_machine = bartMachine(as.data.frame(dat.tr[,-1]), dat.tr[,1], verbose = FALSE)
    # summary(bart_machine)
    pred.bart <- predict(bart_machine, as.data.frame(dat.te[,-1]))
    RMSE[i,j] <- calcRMSE(dat.te[,1], pred.bart)
    predictions[,j] <- pred.bart
    j <-j+1
    # plot(dat.te[,1], pred.bart, col = as.numeric(dat.te$class), main = "BART")
    # abline(0,1)
    
    
    Ktrain2 <- kernelMatrix(kfunc, as.matrix(dat.tr[,2:ncol(dat.tr)]))
    Ktest2 <- kernelMatrix(kfunc, as.matrix(dat.te[,2:ncol(dat.tr)]), as.matrix(dat.tr[,2:ncol(dat.tr)]))
    kpcKern2 <- kPCA(Ktrain2)
    
    
    # plot the data projection on the principal components
    
    # pairs(cbind(dat.tr[,1], kpcKern2$KPCs[ , 1 : 6]), col = dat.tr$class)
    KPCpred3 <- predict(kpcKern2, Ktest2)
    # pairs(KPCpred3[ , 1 : 6], col = dat.te$class)
    lmfit <- lm(dat.tr[, 1]~  kpcKern2$KPCs[ , 1 : 6] )
    
    pca.pred.poor.te <- cbind(rep(1,nrow(KPCpred3)), KPCpred3[ , 1 : 6]) %*% as.matrix(lmfit$coef)
    # plot(dat.te[,1], pca.pred.poor.te , col = as.numeric(dat.te$class), main = "KPCs")
    # abline(0,1)
    
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pca.pred.poor.te)
    predictions[,j] <- pca.pred.poor.te
    j <- j+ 1
    
    
    
    
    rf.fit <- ranger(y = dat.tr[, 1], x =  kpcKern2$KPCs[ , 1 : 6], importance = "impurity")
    # plot(rf.fit$variable.importance)
    pred.rf <- predict(rf.fit , data = KPCpred3[ , 1 : 6])
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pred.rf$predictions)
    predictions[,j] <- pred.rf$predictions
    j <- j+1
    
    
    
    rf.fit <- ranger(y = dat.tr[, 1], x =  Ktrain2, importance = "impurity")
    # plot(rf.fit$variable.importance)
    pred.rf <- predict(rf.fit , data = Ktest2)
    
    RMSE[i,j] <- calcRMSE(dat.te[,1], pred.rf$predictions)
    predictions[,j] <- pred.rf$predictions
    j <- j+1
    # plot(dat.te[,1], pred.rf$predictions, col = as.numeric(dat.te$class), main = "RF")
    # abline(0,1)
    
    
    modelsvm <- svm(V1~.,dat.tr)
    
    #Predict using SVM regression
    predYsvm = predict(modelsvm, dat.te)
    RMSE[i,j] <- calcRMSE(dat.te[,1], predYsvm)
    predictions[,j] <- predYsvm
    j <- j+1
    # plot(dat.te[,1], predYsvm,  main = "RF")
    # abline(0,1)
    
    
    
    RMSE[i,j] <-  calcRMSE(dat.te[,1],   rowMeans(predictions[, 1:13]))
  }
  colnames(RMSE) <- c("lasso", "pca", "pls3", "pls4","pls5",  "uncorr", "RF", "varRF", "Bart", "kpca","kpcaRF", "kernelRF", "svm", "ensamble")
  RMSEFULL[[yind]] <- RMSE
  
  RMSE1 <- RMSE  %>%as.data.frame()%>% mutate(split = as.factor(1:N))%>% pivot_longer(1:14,"ALGORITHM")
  # RMSE1 %>% ggplot(aes(x = ALGORITHM, y = value)) + geom_boxplot()
  p1<- RMSE1 %>% group_by(ALGORITHM) %>% mutate(meanv = mean(value)) %>% ggplot(aes(x= ALGORITHM, y = value, color = split, group = split))+
    geom_point()+ geom_line() +geom_line(aes(x = ALGORITHM, y = meanv), col = 1)+ theme(legend.position = "none")+
    ylab("RMSE")+ ggtitle(names(ys_tr)[yind])
  print(p1)
  
  t1[[yind]] <- RMSE1 %>% group_by(ALGORITHM) %>% summarise(meanv = mean(value), sdv=sd(value)) %>% arrange(desc(meanv))
  
  
  
}



