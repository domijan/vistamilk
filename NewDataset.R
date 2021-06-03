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
glimpse(x)
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
Ktrain <- kernelMatrix(kfunc, x.sc)
# image(Ktrain)

Ktest <- kernelMatrix(kfunc, te, tr)



kpcTrain <- kPCA(Ktrain)


# plot the data projection on the principal components

# pairs(kpcTrain$KPCs[ , 1 : 6], col = labl.Ave)

# newdat <- cbind(log(ys_tr), kpcTrain$KPCs[ , 1 : 6], labl.Ave)
# 
# newdat <- newdat %>% as.data.frame() %>% na.omit() %>% mutate(labl.Ave = as.factor(labl.Ave))
# pairs(newdat[,1:20],  col = newdat[,21])
kpcTest <- predict(kpcTrain, Ktest)


