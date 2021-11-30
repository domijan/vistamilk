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

tr <- read.csv("data/tr.csv")
te <- read.csv("data/te.csv")

calcRMSE <- function(y, yhat) {
  sqrt(mean((yhat - y) ^ 2))
}


y <- tr[, 1:3]
tr <- tr[,-c(1:3)]


a <- apply(tr, 2, mean)
b <- apply(tr, 2, sd)

pairs(log(y)) # all non-normal
matplot(t(tr),
        type = "l",
        col = 9,
        lty = 1)
matplot(t(te),
        type = "l",
        col = 9,
        lty = 1)



tr <- scale(tr, a, b)
te <- scale(te, a, b)
matplot(t(tr),
        type = "l",
        col = 9,
        lty = 1)
matplot(t(te),
        type = "l",
        col = 9,
        lty = 1)



#########################################################
# Look for clusters in the training data


D <- as.matrix(dist(tr, method = "euclidean"))



h <- hclust(as.dist(D), method = "ward.D2")
d <- as.dendrogram(h)
plot(d)

heatmap(D,
        Colv = as.dendrogram(h),
        Rowv = as.dendrogram(h))


d2 <- color_branches(d, k = 7, col = c(2, 3, 5, 4, 6, 7, 8))



labl.Ave <- cutree(d2, 7)
matplot(
  t(tr),
  type = "l",
  col = labl.Ave,
  lty = 1,
  ylab = "training data"
) # see the 7 clusters


#############################################################

kfunc <-  rbfdot(sigma = 0.0005)
# kfunc <- laplacedot(sigma = 0.0001)
# kfunc <- anovadot(sigma = 1, degree = 1)

Ktrain2 <- kernelMatrix(kfunc, tr)
image(Ktrain2)

Ktest2 <- kernelMatrix(kfunc, te, tr)



kpcKern2 <- kPCA(Ktrain2)


# plot the data projection on the principal components

pairs(kpcKern2$KPCs[, 1:6], col = labl.Ave)

newdat <- cbind(log(y), kpcKern2$KPCs[, 1:6], labl.Ave)

newdat <-
  newdat %>% as.data.frame() %>% na.omit() %>% mutate(labl.Ave = as.factor(labl.Ave))
pairs(newdat[, 1:9],  col = newdat[, 10])
KPCpred3 <- predict(kpcKern2, Ktest2)
pairs(KPCpred3[, 1:6])

#############################################################

# correlation with responses

indexcor <-
  abs(cor(cbind(y, tr), use = "pairwise.complete.obs"))[1:3,-c(1:3)]
matplot(t(indexcor),
        type = "l",
        col = 9,
        lty = 1)

matplot(t(tr),
        type = "l",
        col = labl.Ave,
        lty = 1)
lines(indexcor[1, ] * 20, lty = 3)


# clustering the columns (variables) rather than observations
D2 <- as.matrix(1 - abs(cor(tr)))

# image(D)

h2 <- hclust(as.dist(D2), method = "ward.D2")
dnew <- as.dendrogram(h2)
plot(dnew)

# library(pheatmap)
# pheatmap(tr, cluster_cols = TRUE, clustering_distance_cols = "correlation",  cutree_cols  = 15)
heatmap(D2,
        Colv = as.dendrogram(h2),
        Rowv = as.dendrogram(h2))


d2 <- color_branches(dnew, k = 7, col = 2:8)



plot(reorder(d2, D2, method = "OLO"))

labl.cols <- cutree(d2, 15) # 15 variables
matplot(t(tr),
        type = "l",
        col = labl.Ave,
        lty = 1)
points(1:1060, rep(-5.8, 1060), col = labl.cols)

# Pick 15 top features
nonz1 <- NULL
nonz2 <- NULL
nonz3 <- NULL

for (j in unique(labl.cols)) {
  nonz1[j] <-
    (which(indexcor[1,] == max(indexcor[1, labl.cols == j])))
  nonz2[j] <-
    (which(indexcor[2,] == max(indexcor[2, labl.cols == j])))
  nonz3[j] <-
    (which(indexcor[3,] == max(indexcor[3, labl.cols == j])))
}

matplot(t(tr),
        type = "l",
        col = labl.Ave,
        lty = 1)
abline(v = nonz1, col = 1:15) # pick values in a sensible region
points(1:1060, rep(-5.8, 1060), col = labl.cols)
# max(indexcor[1,])

pairs(cbind(y[, 1], tr[, nonz1]), col = labl.Ave)

#############################################################
# Each response: to pick the model

N <- 50
random_split <- replicate(N, sample(nrow(dat), 200))

dat <- cbind((y[, 1]), tr)
dat <- as.data.frame(dat)
dat <- dat %>% mutate(class = as.factor(labl.Ave))
dat %>% ggplot(aes(x = V1)) + geom_histogram()
dat <- dat %>% na.omit()
nonz <- nonz1

RMSE <- matrix(0, N, 15)
predictions <- matrix(0,   nrow(dat) - 200, 14)

for (i in 1:N) {
  j <- 1
  dat.tr <- dat[random_split[, i],]
  dat.te <- dat[-random_split[, i],]
  
  #############################################################
  # lasso
  #############################################################
  
  x <- model.matrix(V1 ~ . , data = dat.tr[,-1062])
  
  y1 <- dat.tr[, 1]
  
  grid <- 10 ^ seq(-3, 3, length = 100)
  
  lasso.fit <- glmnet(x, y1, alpha = 1, lambda = grid) # for lasso
  
  # plot(lasso.fit)
  cv.out <- cv.glmnet(x, y1, alpha = 1)
  # print(cv.out$lambda.min)
  
  lasso.fit <- glmnet(x, y1, alpha = 1,
                      lambda = cv.out$lambda.min)
  # head(coef(lasso.fit))
  # nonzero <- which(coef(lasso.fit)[3:1062]!=0) #intercept there twice?!
  
  # matplot(t(tr), type = "l", col = labl.Ave, lty = 1)
  # abline(v = nonzero) # pick values in a sensible region
  
  # pairs(cbind(y1, dat.tr[,1 + nonzero])) # still perfectly correlated predictors
  
  lasso.pred <-
    predict(lasso.fit, newx = model.matrix(V1 ~ . , data = dat.te[,-1062]))
  
  RMSE[i, j] <- calcRMSE(dat.te[, 1], lasso.pred)
  predictions[, j] <- lasso.pred
  j <- j + 1
  
  # plot(dat.te[,1], lasso.pred, col = as.numeric(dat.te$class), main = "lasso")
  # abline(a = 0, b = 1)
  # plot(lasso.fit)
  #############################################################
  # pca
  #############################################################
  

  pc <- prcomp(dat.tr[, 2:1061], scale = TRUE)
  # pairs(cbind(dat.tr[,1], pc$x[,1:9]), col = as.numeric(dat.tr$class))
  # plot(pc)
  # dim(pc$x)
  # summary(lm(dat.tr[, 1]~pc$x[,1:10]))
  # library(car)
  # avPlots(lm(dat.tr[, 1]~pc$x[,1:10]))
  
  lm.fit <- lm(dat.tr[, 1] ~ pc$x[, 1:5])
  # plot(lm.fit,1)
  # anova(lm.fit)
  # summary(lm.fit)
  x.te <- predict(pc, dat.te[, 2:1061])
  # pca.pred.poor <- as.matrix(cbind(rep(1, 200), pc$x[,1:5])) %*% as.matrix(lm.fit$coef[1:6])
  lm.pred <-
    as.matrix(cbind(rep(1, nrow(x.te)), x.te[, 1:5])) %*% as.matrix(lm.fit$coef[1:6])
  # plot(dat.tr[,1], pca.pred.poor, col = as.numeric(dat.tr$class))
  # abline(0,1)
  # pca.pred <- predict(lm.fit)
  # calcRMSE(dat.tr[,1], pca.pred.poor)
  RMSE[i, j] <- calcRMSE(dat.te[, 1], lm.pred)
  predictions[, j] <- lm.pred
  j <- j + 1
  # plot(dat.te[,1], lm.pred, col = as.numeric(dat.te$class), main = "PCA")
  # abline(a = 0, b = 1)
  
  #############################################################
  # pls
  #############################################################
  

  pls.fit <-
    plsr(V1 ~ .,
         data = dat.tr[,-1062],
         scale = TRUE,
         validation = "CV")
  # summary(pls.fit)
  # validationplot(pls.fit, val.type = "MSEP")
  
  pls.pred <- predict(pls.fit, dat.te[, -1], ncomp = 3)
  
  RMSE[i, j] <- calcRMSE(dat.te[, 1], pls.pred)
  predictions[, j] <- pls.pred
  j <- j + 1
  
  pls.pred <- predict(pls.fit, dat.te[, -1], ncomp = 4)
  
  RMSE[i, j] <- calcRMSE(dat.te[, 1], pls.pred)
  predictions[, j] <- pls.pred
  j <- j + 1
  
  pls.pred <- predict(pls.fit, dat.te[, -1], ncomp = 5)
  
  RMSE[i, j] <- calcRMSE(dat.te[, 1], pls.pred)
  predictions[, j] <- pls.pred
  j <- j + 1
  # plot(dat.te[,1],pls.pred,  col = as.numeric(dat.te$class), main = "PLS")
  # abline(a = 0, b = 1)
  
  # summary(pls.fit)
  # validationplot(pls.fit, val.type = "MSEP")

  #############################################################
  # lm +15
  #############################################################
  
  lm.fit <- lm(dat.tr[, 1] ~ as.matrix(dat.tr[, nonz + 1]))
  # plot(lm.fit,1)
  # anova(lm.fit)
  # summary(lm.fit)
  
  x <- model.matrix(V1 ~ . , data = dat.te[, c(1, nonz + 1)])
  
  lm.pred <- x %*% as.matrix(lm.fit$coef)
  # plot(dat.te[,1], lm.pred , col = as.numeric(dat.te$class), main = "Uncorrelated variables")
  # abline(0,1)
  
  RMSE[i, j] <- calcRMSE(dat.te[, 1], lm.pred)
  predictions[, j] <- lm.pred
  j <- j + 1
  
  #############################################################
  # rf
  #############################################################
  
  
  rf.fit <- ranger(V1 ~ ., data = dat.tr, importance = "impurity")
  # plot(rf.fit$variable.importance)
  rf.pred <- predict(rf.fit , data = dat.te)
  
  RMSE[i, j] <- calcRMSE(dat.te[, 1], rf.pred$predictions)
  predictions[, j] <- rf.pred$predictions
  j <- j + 1
  # plot(dat.te[,1], rf.pred$predictions, col = as.numeric(dat.te$class), main = "RF")
  # abline(0,1)
  
  
  #############################################################
  # rf +vs
  #############################################################
  
  rf.fit <-
    ranger(
      V1 ~ .,
      data = dat.tr,
      importance = "impurity",
      regularization.factor = 0.2,
      regularization.usedepth = FALSE
    )
  # plot(rf.fit$variable.importance)
  rf.pred <- predict(rf.fit , data = dat.te)
  
  RMSE[i, j] <- calcRMSE(dat.te[, 1], rf.pred$predictions)
  predictions[, j] <- rf.pred$predictions
  j <- j + 1
  # plot(dat.te[,1], rf.pred$predictions, col = as.numeric(dat.te$class))
  # abline(0,1)
  
  #############################################################
  # bart
  #############################################################
  
  
  bart.fit <- bartMachine(as.data.frame(dat.tr[, -1]), dat.tr[, 1])
  # summary(bart.fit)
  bart.pred <- predict(bart.fit, as.data.frame(dat.te[, -1]))
  RMSE[i, j] <- calcRMSE(dat.te[, 1], bart.pred)
  predictions[, j] <- bart.pred
  j <- j + 1
  # plot(dat.te[,1], bart.pred, col = as.numeric(dat.te$class), main = "BART")
  # abline(0,1)
  
  
  #############################################################
  # lm + kPCA
  #############################################################
  
  Ktrain2 <- kernelMatrix(kfunc, as.matrix(dat.tr[, 2:1061]))
  Ktest2 <-
    kernelMatrix(kfunc, as.matrix(dat.te[, 2:1061]), as.matrix(dat.tr[, 2:1061]))
  kpcKern2 <- kPCA(Ktrain2)
  
  
  # plot the data projection on the principal components
  
  # pairs(cbind(dat.tr[,1], kpcKern2$KPCs[ , 1 : 6]), col = dat.tr$class)
  KPCpred3 <- predict(kpcKern2, Ktest2)
  # pairs(KPCpred3[ , 1 : 6], col = dat.te$class)
  lm.fit <- lm(dat.tr[, 1] ~  kpcKern2$KPCs[, 1:6])
  # plot(lm.fit,1)
  # anova(lm.fit)
  # summary(lm.fit)
  # predict(lm.fit, KPCpred3[ , 1 : 6])
  
  lm.pred <-
    cbind(rep(1, nrow(KPCpred3)), KPCpred3[, 1:6]) %*% as.matrix(lm.fit$coef)
  # plot(dat.te[,1], lm.pred , col = as.numeric(dat.te$class), main = "KPCs")
  # abline(0,1)
  
  
  RMSE[i, j] <- calcRMSE(dat.te[, 1], lm.pred)
  predictions[, j] <- lm.pred
  j <- j + 1
  
  
  #############################################################
  # rf + kPCA
  #############################################################
  
  
  rf.fit <-
    ranger(y = dat.tr[, 1],
           x =  kpcKern2$KPCs[, 1:6],
           importance = "impurity")
  # plot(rf.fit$variable.importance)
  rf.pred <- predict(rf.fit , data = KPCpred3[, 1:6])
  
  RMSE[i, j] <- calcRMSE(dat.te[, 1], rf.pred$predictions)
  predictions[, j] <- rf.pred$predictions
  j <- j + 1
  
  #############################################################
  # rf + kernel
  #############################################################
  
  
  rf.fit <-
    ranger(y = dat.tr[, 1],
           x =  Ktrain2,
           importance = "impurity")
  # plot(rf.fit$variable.importance)
  rf.pred <- predict(rf.fit , data = Ktest2)
  
  RMSE[i, j] <- calcRMSE(dat.te[, 1], rf.pred$predictions)
  predictions[, j] <- rf.pred$predictions
  j <- j + 1
  # plot(dat.te[,1], rf.pred$predictions, col = as.numeric(dat.te$class), main = "RF")
  # abline(0,1)
  
  #############################################################
  # svm
  #############################################################
  
  svm.fit <- svm(V1 ~ ., dat.tr[,-1062])
  
  #Predict using SVM regression
  svm.pred <-  predict(svm.fit, dat.te[,-1062])
  RMSE[i, j] <- calcRMSE(dat.te[, 1], svm.pred)
  predictions[, j] <- svm.pred
  j <- j + 1
  # plot(dat.te[,1], svm.pred, col = as.numeric(dat.te$class), main = "RF")
  # abline(0,1)
  
  #############################################################
  # rf + kPCA + class
  #############################################################
  
  rf.fit <-
    ranger(
      y = dat.tr[, 1],
      x =  cbind(kpcKern2$KPCs[, 1:6], dat.tr$class),
      importance = "impurity"
    )
  # plot(rf.fit$variable.importance)
  rf.pred <-
    predict(rf.fit , data = cbind(KPCpred3[, 1:6], dat.te$class))
  
  RMSE[i, j] <- calcRMSE(dat.te[, 1], rf.pred$predictions)
  predictions[, j] <- rf.pred$predictions
  j <- j + 1
  
  
  RMSE[i, j] <-
    calcRMSE(dat.te[, 1],   rowMeans(predictions[, 1:13]))
}
#############################################################

colnames(RMSE) <-
  c(
    "lasso",
    "pca",
    "pls3",
    "pls4",
    "pls5",
    "uncorr",
    "RF",
    "varRF",
    "Bart",
    "kpca",
    "kpcaRF",
    "kernelRF",
    "svm",
    "kpcaClassRF",
    "ensemble"
  )

#############################################################


RMSE1 <-
  RMSE  %>% 
  as.data.frame() %>% 
  mutate(split = as.factor(1:N)) %>% 
  pivot_longer(1:15, "ALGORITHM")

# RMSE1 %>% ggplot(aes(x = ALGORITHM, y = value)) + geom_boxplot()

RMSE1 %>% 
  group_by(ALGORITHM) %>% 
  mutate(meanv = mean(value)) %>% 
  ggplot(aes(
  x = ALGORITHM,
  y = value,
  color = split,
  group = split
)) +
  geom_point() + 
  geom_line() + 
  geom_line(aes(x = ALGORITHM, y = meanv), col = 1) + 
  theme(legend.position = "none") +
  ylab("RMSE")

t1 <-
  RMSE1 %>% 
  group_by(ALGORITHM) %>% 
  summarise(meanv = mean(value), sdv = sd(value)) %>% 
  arrange(desc(meanv))


RMSE2 <-
  RMSE  %>% 
  as.data.frame() %>% 
  mutate(split = as.factor(1:N)) %>% 
  pivot_longer(1:15, "ALGORITHM")

# RMSE1 %>% ggplot(aes(x = ALGORITHM, y = value)) + geom_boxplot()
RMSE2 %>% 
  group_by(ALGORITHM) %>% 
  mutate(meanv = mean(value)) %>% 
  ggplot(aes(
  x = ALGORITHM,
  y = value,
  color = split,
  group = split
)) +
  geom_point() + 
  geom_line() + 
  geom_line(aes(x = ALGORITHM, y = meanv), col = 1) + 
  theme(legend.position = "none") +
  ylab("RMSE")

t2 <-
  RMSE2 %>% 
  group_by(ALGORITHM) %>% 
  summarise(meanv = mean(value), sdv = sd(value)) %>% 
  arrange(desc(meanv))


RMSE3 <-
  RMSE  %>% 
  as.data.frame() %>% 
  mutate(split = as.factor(1:N)) %>% 
  pivot_longer(1:15, "ALGORITHM")

# RMSE1 %>% ggplot(aes(x = ALGORITHM, y = value)) + geom_boxplot()

RMSE3 %>% 
  group_by(ALGORITHM) %>% 
  mutate(meanv = mean(value)) %>% 
  ggplot(aes(
  x = ALGORITHM,
  y = value,
  color = split,
  group = split
)) +
  geom_point() + 
  geom_line() + 
  geom_line(aes(x = ALGORITHM, y = meanv), col = 1) + 
  theme(legend.position = "none") +
  ylab("RMSE")

t3 <-
  RMSE3 %>% 
  group_by(ALGORITHM) %>% 
  summarise(meanv = mean(value), sdv = sd(value)) %>% 
  arrange(desc(meanv))


# USE RESULTS FROM ENSEMBLE

#############################################################
#############################################################
# train on full training set
#############################################################
#############################################################



dat <- cbind(y[,3], tr)
dat <- as.data.frame(dat)
dat <- dat %>% mutate(class = as.factor(labl.Ave))
dat %>% ggplot(aes(x= V1)) + geom_histogram()
dat <- dat %>% na.omit()
nonz <- nonz3
# dat <- dat %>% filter(V1<15)
# dim(dat)


predictions <- matrix(0,nrow(te), 14)


#############################################################

j <- 1
dat.tr <- dat
dat.te <- as.data.frame(te)

#############################################################
# lasso
#############################################################

x <- model.matrix(V1 ~ . , data = dat.tr[,-1062])
y1 <- dat.tr[, 1]
grid <- 10 ^ seq(-3, 3, length = 100)
lasso.fit <- glmnet(x, y1, alpha = 1, lambda = grid) # for lasso

cv.out <- cv.glmnet(x, y1, alpha = 1)

lasso.fit <- glmnet(x, y1, alpha = 1,
                    lambda = cv.out$lambda.min)
# head(coef(lasso.fit))

lasso.pred <- predict(lasso.fit, newx = model.matrix( ~ ., dat.te))
predictions[, j] <- lasso.pred
j <- j + 1

#############################################################
# lm +pca
#############################################################

pc <- prcomp(dat.tr[, 2:1061], scale = TRUE)
pairs(cbind(dat.tr[, 1], pc$x[, 1:6]), col = as.numeric(dat.tr[, 1062]))

lm.fit <- lm(dat.tr[, 1] ~ pc$x[, 1:6])

x.te <- predict(pc, dat.te)

lm.pred <-
  as.matrix(cbind(rep(1, nrow(x.te)), x.te[, 1:6])) %*% as.matrix(lm.fit$coef[1:7])
predictions[, j] <- lm.pred
j <- j + 1

#############################################################
# pls
#############################################################

pls.fit <-
  plsr(V1 ~ .,
       data = dat.tr[,-1062],
       scale = TRUE,
       validation = "CV")


pls.pred <- predict(pls.fit, dat.te, ncomp = 3)
predictions[, j] <- pls.pred
j <- j + 1

pls.pred <- predict(pls.fit, dat.te, ncomp = 4)
predictions[, j] <- pls.pred
j <- j + 1

pls.pred <- predict(pls.fit, dat.te, ncomp = 5)

predictions[, j] <- pls.pred
j <- j + 1

#############################################################
# lm +15
#############################################################

lm.fit <- lm(dat.tr[, 1] ~ as.matrix(dat.tr[, nonz + 1]))
# plot(lm.fit,1)
# anova(lm.fit)
# summary(lm.fit)

x <- model.matrix( ~ . , data = dat.te[, nonz])

lm.pred <- x %*% as.matrix(lm.fit$coef)

predictions[, j] <- lm.pred
j <- j + 1

#############################################################
# rf
#############################################################


rf.fit <-
  ranger(V1 ~ ., data = dat.tr[, -1062], importance = "impurity")
# plot(rf.fit$variable.importance)
rf.pred <- predict(rf.fit , data = dat.te)

predictions[, j] <- rf.pred$predictions
j <- j + 1

#############################################################
# rf +vs
#############################################################


rf.fit <-
  ranger(
    V1 ~ .,
    data = dat.tr[, -1062],
    importance = "impurity",
    regularization.factor = 0.2,
    regularization.usedepth = FALSE
  )
# plot(rf.fit$variable.importance)
rf.pred <- predict(rf.fit , data = dat.te)

predictions[, j] <- rf.pred$predictions
j <- j + 1


#############################################################
# bart
#############################################################


bart.fit <-
  bartMachine(as.data.frame(dat.tr[, -c(1, 1062)]), dat.tr[, 1])
# summary(bart.fit)
bart.pred <- predict(bart.fit, as.data.frame(dat.te))
predictions[, j] <- bart.pred
j <- j + 1


#############################################################
# lm + kPCA
#############################################################

Ktrain2 <- kernelMatrix(kfunc, as.matrix(dat.tr[, 2:1061]))
# image(Ktrain2)
Ktest2 <-
  kernelMatrix(kfunc, as.matrix(dat.te), as.matrix(dat.tr[, 2:1061]))
# image(Ktest2)
kpcKern2 <- kPCA(Ktrain2)

# plot the data projection on the principal components

# pairs(cbind(dat.tr[,1], kpcKern2$KPCs[ , 1 : 6]), col = dat.tr$class)
KPCpred3 <- predict(kpcKern2, Ktest2)
# pairs(KPCpred3[ , 1 : 6], col = dat.te$class)

#############################################################

lm.fit <- lm(dat.tr[, 1] ~  kpcKern2$KPCs[, 1:6])
# plot(lm.fit,1)
# anova(lm.fit)
# summary(lm.fit)
# predict(lm.fit, KPCpred3[ , 1 : 6])


lm.pred <-
  cbind(rep(1, nrow(KPCpred3)), KPCpred3[, 1:6]) %*% as.matrix(lm.fit$coef)
predictions[, j] <- lm.pred
j <- j + 1

#############################################################
# rf + kPCA
#############################################################


rf.fit <-
  ranger(y = dat.tr[, 1],
         x =  kpcKern2$KPCs[, 1:6],
         importance = "impurity")
# plot(rf.fit$variable.importance)
rf.pred <- predict(rf.fit , data = KPCpred3[, 1:6])
predictions[, j] <- rf.pred$predictions
j <- j + 1

#############################################################
# rf kernel
#############################################################


rf.fit <-
  ranger(y = dat.tr[, 1], x =  Ktrain2, importance = "impurity")
# plot(rf.fit$variable.importance)
rf.pred <- predict(rf.fit , data = Ktest2)

predictions[, j] <- rf.pred$predictions
j <- j + 1
# plot(dat.te[,1], rf.pred$predictions, col = as.numeric(dat.te$class), main = "RF")
# abline(0,1)

#############################################################
# svm 
#############################################################


svm.fit <- svm(V1 ~ ., dat.tr[,-1062])

#Predict using SVM regression
svm.pred <-  predict(svm.fit, dat.te)
predictions[, j] <- svm.pred
j <- j + 1

#############################################################
# ensemble
#############################################################

predictions[, j] <- rowMeans(predictions[, 1:13])
dim(predictions)
colnames(predictions) <-
  c(
    "lasso",
    "pca",
    "pls3",
    "pls4",
    "pls5",
    "uncorr",
    "RF",
    "varRF",
    "Bart",
    "kpca",
    "kpcaRF",
    "kernelRF",
    "svm",
    "ensemble"
  )


#############################################################

predictions1 <- predictions  %>%
  as.data.frame() %>% 
  mutate(obs = as.factor(1:nrow(te))) %>% 
  pivot_longer(1:14,"ALGORITHM")
# RMSE1 %>% ggplot(aes(x = ALGORITHM, y = value)) + geom_boxplot()
predictions1 %>%  
  group_by(ALGORITHM) %>% 
  mutate(meanv = mean(value)) %>% 
  ggplot(aes(x= ALGORITHM, y = value, color = obs, group = obs)) +
  geom_point() + 
  geom_line() + 
  theme(legend.position = "none") +
  ylab("predictions")

summary(y1)

write.csv(predictions[,14], "kappa_casein.csv", row.names = FALSE, col.names = FALSE)




predictions2 <- predictions  %>%
  as.data.frame() %>% 
  mutate(obs = as.factor(1:nrow(te))) %>% 
  pivot_longer(1:14,"ALGORITHM")

# RMSE1 %>% ggplot(aes(x = ALGORITHM, y = value)) + geom_boxplot()
predictions2 %>% 
  group_by(ALGORITHM) %>% 
  mutate(meanv = mean(value)) %>% 
  ggplot(aes(x= ALGORITHM, y = value, color = obs, group = obs)) +
  geom_point() + 
  geom_line() + 
  theme(legend.position = "none") +
  ylab("predictions")

summary(exp(predictions))
summary(y[,2])

write.csv(exp(predictions[,14]), "Casein_micelle_size.csv", row.names = FALSE, col.names = FALSE)



predictions3 <- predictions %>% 
  as.data.frame() %>% 
  mutate(obs = as.factor(1:nrow(te))) %>% 
  pivot_longer(1:14, "ALGORITHM")
# RMSE1 %>% ggplot(aes(x = ALGORITHM, y = value)) + geom_boxplot()
predictions3 %>% group_by(ALGORITHM) %>% mutate(meanv = mean(value)) %>% ggplot(aes(
  x = ALGORITHM,
  y = value,
  color = obs,
  group = obs
)) +
  geom_point() + 
  geom_line() + 
  theme(legend.position = "none") +
  ylab("predictions")

summary((predictions))
summary(y[, 3])

write.csv(predictions[, 14],
          "Native_pH.csv",
          row.names = FALSE,
          col.names = FALSE)
