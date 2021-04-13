library(tidyverse)
library(rgl)
library(dendextend)
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

pairs(y) # all non-normal
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


D <- as.matrix(dist(tr,method="euclidean"))
# image(D)

h<- hclust(as.dist(D), method = "ward.D2")
d <- as.dendrogram(h)
plot(d)

heatmap(D,
        Colv=as.dendrogram(h),     
        Rowv=as.dendrogram(h))


d2 <- color_branches(d,k=7, col=c(2,3,5,4, 6,7,8)) 
# auto-coloring 4 clusters of branches.
# plot(d2)

plot(reorder(d2, D, method = "OLO")) 

labl.Ave <- cutree(d2,7)
matplot(t(tr), type = "l", col = labl.Ave, lty = 1)


#############################################################

library(BKPC)
?bkpc

#############################################################

indexcor <- abs(cor(cbind(y,tr), use = "pairwise.complete.obs"))[1:3, -c(1:3)]
matplot(t(indexcor), type = "l", col = 9, lty = 1)

matplot(t(tr), type = "l", col = labl.Ave, lty = 1)
lines(indexcor[1,]*20,lty = 3)


dim(D2)
D2 <- as.matrix(1-abs(cor(tr)))

# image(D)
D2[2,1]

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

labl.cols <- cutree(d2,7)
matplot(t(tr), type = "l", col = labl.Ave, lty = 1)
points(1:1060, rep(-5.8, 1060), col = labl.cols)

#############################################################

# First response:

dat <- cbind(y[,1], tr)
dat <- as.data.frame(dat)
dat <- dat %>% na.omit()

dat <- dat %>% filter(V1<15)
dim(dat)
# pairs(dat[,1:4])
# dat %>% ggplot(aes(x= (V1)))+geom_histogram()
# dat %>% ggplot(aes(x= log(V1)))+geom_histogram()
# qqnorm(y=(dat$V1[dat$V1<15]))
# qqnorm(y=log(dat$V1))

x <- model.matrix(V1 ~. , data = dat)

y1 <- dat[,1]

grid <- 10^seq(-3, 3, length = 100)


library(glmnet) # for the ridge regression


lasso.fit <- glmnet(x,y1,alpha=1, lambda = grid) # for lasso

plot(lasso.fit)

cv.out <- cv.glmnet(x,y1,alpha=1)
cv.out$lambda.min

lasso.fit <- glmnet(x,y1,alpha=1, 
                    lambda = cv.out$lambda.min)
nonzero <- which(coef(lasso.fit)[3:1062]!=0) #intercept there twice?!

matplot(t(tr), type = "l", col = labl.Ave, lty = 1)
abline(v = nonzero) # pick values in a sensible region


pairs(cbind(y1, dat[,1+nonzero[1:10]]),) # still perfectly correlated predictors

lasso.pred <- predict(lasso.fit, newx = model.matrix(V1 ~. , data = dat))
calcRMSE(y1, lasso.pred)
plot(y1, lasso.pred)

pc <- prcomp(dat[,-1], scale = TRUE)
pairs(cbind(y[,1], pc$x[,1:9]))

dim(pc$x)
summary(lm(dat[, 1]~pc$x[,1:10]))
library(car)
avPlots(lm(dat[, 1]~pc$x[,1:10]))


lmfit<- lm(dat[, 1]~pc$x[,1:397])
plot(lmfit,1)
pca.pred <- predict(lmfit)
calcRMSE(dat[,1], pca.pred)
plot(dat[,1], pca.pred)



library(pls)
pcr.fit <- pcr(V1~., data=dat, scale= TRUE, validation = "CV")
summary(pcr.fit)
validationplot(pcr.fit, val.type = "MSEP")
pcr.pred <- predict(pcr.fit, dat[,-1], ncomp = 356)
calcRMSE(dat[,1], pcr.pred)
plot(dat[,1], pcr.pred)


pls.fit <- plsr(V1 ~., data=dat, scale= TRUE, validation = "CV")
summary(pls.fit)
validationplot(pls.fit, val.type = "MSEP")
pls.pred <- predict(pls.fit, dat[,-1], ncomp = 250)
calcRMSE(dat[,1], pls.pred)
plot(dat[,1],pls.pred)
