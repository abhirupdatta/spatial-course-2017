library(ggplot2)
library(RColorBrewer)
library(akima)
library(fields)
library(ggmap)
library(maps)
library(spBayes)
library(MBA)
library(classInt)
library(plotrix)
library(geoR)

### simulated datasets ####
data1=read.csv("../data/dataset1.csv")
data2=read.csv("../data/dataset2.csv")
data3=read.csv("../data/dataset3.csv")

col.br <- colorRampPalette(c("midnightblue", "cyan", "yellow", "red"))
col.pal <- col.br(5)

### function for plotting interpolated surface of a column of a data table
myplot=function(tab,colname){
    
    surf <- mba.surf(tab[,c("sx","sy",colname)], no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
    dev.new()
    image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))
    
}

### data and covariate surface plots
dev.new()
myplot(data1,"y")
myplot(data2,"y")
myplot(data3,"y")
myplot(data3,"x")

### linear regression on x
lmobj1=lm(y~x,data=data1)
lmobj2=lm(y~x,data=data2)
lmobj3=lm(y~x,data=data3)

summary(lmobj1)
summary(lmobj2)
summary(lmobj3)

data1$res=lmobj1$residuals
data2$res=lmobj2$residuals
data3$res=lmobj3$residuals

### residual surface plots
myplot(data1,"res")
myplot(data2,"res")
myplot(data3,"res")

### empirical variograms ###
max.dist <- 0.75*max(iDist(data1[,1:2]))
bins <- 20

n=nrow(data1)
ydiff=diff=rep(0,n^2)
for(i in 2:n) for(j in 1:(i-1)) {
  ydiff[(i-1)*n+j]=(data1$y[i]-data1$y[j])^2
  diff[(i-1)*n+j]=sqrt(sum(data1[i,1:2]-data1[j,1:2])^2)
  }

### variogram cloud 
dev.new()
plot(diff,ydiff,xlab="",ylab="")

### binned variograms
vario1raw <- variog(coords=data1[,1:2], data=data1$y, uvec=(seq(0, max.dist, length=bins)))
plot(vario1raw,pch=16)

vario1 <- variog(coords=data1[,1:2], data=data1$res, uvec=(seq(0, max.dist, length=bins)))
plot(vario1,pch=16)

vario2 <- variog(coords=data1[,1:2], data=data2$res, uvec=(seq(0, max.dist, length=bins)))
plot(vario2,pch=16)

vario3 <- variog(coords=data1[,1:2], data=data3$res, uvec=(seq(0, max.dist, length=bins)))
plot(vario3,pch=16)

#### adding the coordinates into the regression
lmobj2s=lm(y~x+sx+sy,data=data2)
data2$res2=lmobj2s$residuals
myplot(data2,"res2")

vario2s <- variog(coords=data1[,1:2], data=data2$res2, uvec=(seq(0, max.dist, length=bins)))
plot(vario2s,pch=16)

lmobj3s=lm(y~x+sx+sy,data=data3)
data3$res2=lmobj3s$residuals
myplot(data3,"res2")

vario3s <- variog(coords=data1[,1:2], data=data3$res2, uvec=(seq(0, max.dist, length=bins)))
plot(vario3s,pch=16)

### variograms for WEF data ###
data(WEF.dat)

WEF.dat <- WEF.dat[!apply(WEF.dat[,c("East_m","North_m","DBH_cm","Tree_height_m","ELEV_m")], 1, function(x)any(is.na(x))),]

DBH <- WEF.dat$DBH_cm

coords <- as.matrix(WEF.dat[,c("East_m","North_m")])

col.br <- colorRampPalette(c("midnightblue", "cyan", "yellow", "red"))
col.pal <- col.br(5)

surf <- mba.surf(cbind(coords,DBH), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))

max.dist=0.25*max(iDist(coords))
bins=50

vario.DBH <- variog(coords=coords, data=DBH, uvec=(seq(0, max.dist, length=bins)))
plot(vario.DBH)

lm.DBH <- lm(DBH~Species, data=WEF.dat)
summary(lm.DBH)
DBH.resid <- resid(lm.DBH)

surf <- mba.surf(cbind(coords,DBH.resid), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))

vario.DBH.resid <- variog(coords=coords, data=DBH.resid, uvec=(seq(0, max.dist, length=bins)))
plot(vario.DBH.resid)
