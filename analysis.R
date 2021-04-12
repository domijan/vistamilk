library(tidyverse)
library(rgl)
tr <- read.csv("tr.csv") 
te <- read.csv("te.csv") 
# dim(tr)
# dim(te)

y <- tr[,1:3]
tr <- tr[, -c(1:3)]


pairs(y)
points3d(y[,1], y[,2], y[,3])


matplot(t(tr), type = "l", col = 9, lty = 1)
matplot(t(te), type = "l", col = 9, lty = 1)
