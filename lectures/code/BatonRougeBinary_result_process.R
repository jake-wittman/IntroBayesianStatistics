#####################################################################
## the file reads the posterior samples of w obtained from WinBUGS ##
## and plot the image and contours of their posterior medians in R ##
#####################################################################


library(MBA)
library(fields) ## For using the image.plot function

setwd('D:\\4. Teaching\\UM_SPH8472\\Notes_Sp2017\\5. Univariate_analysis\\lab\\Baton-data')

w.post <- read.table('RatonRougeBinary_w_postsample.txt',header=FALSE)
w.post.mat <- matrix(w.post[,2],5000,50)
w.post.median <- apply(w.post.mat,2,median)


baton.dat <- read.table('BatonRouge.dat',header=T)
coords <- cbind(baton.dat$Longitude,baton.dat$Latitude)


x.res <- 200
y.res <- 200
surf <- mba.surf(cbind(coords, w.post.median), no.X = x.res, no.Y = y.res, h = 10, m = 2, extend = TRUE)$xyz.est
image(surf, xaxs = "r", yaxs = "r", xlab = "Longitude", ylab = "Lantitude")
contour(surf, add=T) ## (Optional) Adds contour lines to the plot


