library("ggplot2")
library("RColorBrewer")
library("akima")
library("fields")
library("ggmap")
library("maps")
library("spBayes")
library("MBA")
library("classInt")
library("plotrix")
library("geoR")
library("sp")
library("maptools")
library("rgdal")
library("classInt")
library("lattice")
library("raster")
library("sf")

### simulated dataset 3 from lecture 1 ####
data3=read.csv("../data/dataset3.csv")

### function for plotting interpolated surface of a column of a data table
myplot=function(tab,colname){
  
  surf <- mba.surf(tab[,c("sx","sy",colname)], no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
  dev.new()
  image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))
  
}

### data and covariate surface plots
col.br <- colorRampPalette(c("midnightblue", "cyan", "yellow", "red"))
col.pal <- col.br(5)

dev.new()
myplot(data3,"y")
myplot(data3,"x")

### holdout data for model testing ###
set.seed(04152017)
N=nrow(data3)
ind=sample(1:N,N/2,replace=FALSE)
datain=data3[ind,]  ## in-sample data
dataout=data3[-ind,]

### linear regression ###
ols=lm(y~x,data=datain)
datain$res=ols$residuals
beta_ols=ols$coefficients
beta_ols_CI=ols$coefficients+1.96*
  t(t(sqrt(diag(summary(ols)$cov.unscaled))*summary(ols)$sigma))%*%c(-1,1)
myplot(datain,"y")
myplot(datain,"x")
myplot(datain,"res")

### empirical variograms ###
max.dist <- 0.75*max(iDist(datain[,1:2]))
bins <- 20
vario3 <- variog(coords=datain[,1:2], data=datain$res, 
  uvec=(seq(0, max.dist, length=bins)))
plot(vario3,pch=16)

vfit3 <-variofit(vario3, ini.cov.pars=c(0.1,1), ##sigma^2 and 1/phi 
  cov.model="exponential", minimisation.function="optim",
  nugget=0.01, weights="equal")

lines(vfit3,col="red",lwd=2)

### ML estimation ##
mle3 <- likfit(coords=datain[,1:2],data=datain$y,
  trend = trend.spatial(~x, datain), ini.cov.pars=c(0.15,0.25),
  nugget = 0.009,cov.model="exponential",nospatial=TRUE)
lines(mle3,col="blue",lwd=2)

beta_mle=mle3$beta
beta_mle_CI=mle3$beta+1.96*t(t(sqrt(diag(mle3$beta.var))))%*%c(-1,1)

### verifying that the nospatial settings gives the same estimates as ols ###
mle3$nospatial$beta

### AIC and BIC ###
mle3$AIC
mle3$nospatial$AIC.ns

mle3$BIC
mle3$nospatial$BIC.ns

### predictions on the holdout sample ###
Xout=cbind(1,dataout$x)
pred_ols=as.vector(Xout%*%beta_ols)

krig_variofit=krige.conv(coords=datain[,1:2],data=datain$res,
    locations=dataout[,1:2],krige=krige.control(type.krige="SK",obj.model=vfit3,
    beta=0))
pred_variofit=pred_ols+krig_variofit$predict

krig_mlefit=krige.conv(coords=datain[,1:2],data=datain$y,
  locations=dataout[,1:2],krige=krige.control(type.krige="OK",obj.model=mle3,
  trend.d=trend.spatial(~x, datain),trend.l=trend.spatial(~x, dataout)))
pred_mlefit=krig_mlefit$predict
var_mlefit=krig_mlefit$krige.var

### RMSPE ###
sqrt(mean((dataout$y-pred_ols)^2))
sqrt(mean((dataout$y-pred_variofit)^2))
sqrt(mean((dataout$y-pred_mlefit)^2))

### Coverage Prob and CI length ###
CI_ols=pred_ols+1.96*summary(ols)$sigma*cbind(-rep(1,N/2),rep(1,N/2))
mean(CI_ols[,1]<dataout$y & CI_ols[,2]>dataout$y)
mean(CI_ols[,2]-CI_ols[,1])

CI_vario=pred_variofit+1.96*sqrt(krig_variofit$krige.var)%*%t(c(-1,1))
mean(CI_vario[,1]<dataout$y & CI_vario[,2]>dataout$y)
mean(CI_vario[,2]-CI_vario[,1])

CI_mle=pred_mlefit+1.96*sqrt(krig_mlefit$krige.var)%*%t(c(-1,1))
mean(CI_mle[,1]<dataout$y & CI_mle[,2]>dataout$y)
mean(CI_mle[,2]-CI_mle[,1])

### plotting kriged surface ###
xo=yo=seq(0,1,0.02)
s=expand.grid(xo,yo)
surface_X=0.5*sin(10*s[,1]*s[,2])+1*(0.5-s[,1])^2

krig_surface_mlefit=krige.conv(coords=datain[,1:2],data=datain$y,
  locations=s,krige=krige.control(type.krige="OK",obj.model=mle3,
    trend.d=trend.spatial(~x, datain),trend.l=trend.spatial(~surface_X)))

pred_surface_mlefit=krig_surface_mlefit$predict
var_surface_mlefit=krig_surface_mlefit$krige.var

surface_krig_tab=cbind(s,surface_X,pred_surface_mlefit,var_surface_mlefit)
colnames(surface_krig_tab)=c("sx","sy","x","yhat","vyhat")

myplot(surface_krig_tab,"yhat")
myplot(data3,"y")

myplot(surface_krig_tab,"vyhat")
points(datain[,1],datain[,2],pch=16)



######################################################
####  Bartlett Experimental Forestry (BEF) data   ####
######################################################

## Data preliminaries
data(BEF.dat)
BEF.dat <- BEF.dat[BEF.dat$ALLBIO02_KGH>0,]
bio <- BEF.dat$ALLBIO02_KGH*0.001;
log.bio <- log(bio)
## Extract the coordinates
coords <- as.matrix(BEF.dat[,c("XUTM","YUTM")])

## Make a surface plot
x.res <- 100; y.res <- 100

surf <- mba.surf(cbind(coords, log.bio), no.X=x.res, no.Y=y.res, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)")
points(coords)

BEF.small=BEF.dat[,c("XUTM","YUTM","ELEV","SLOPE","SUM_02_TC1","SUM_02_TC2","SUM_02_TC3")]
BEF.small$logbio=log.bio


### linear regression ###
ols=lm(logbio~ELEV+SLOPE+SUM_02_TC1+SUM_02_TC2+SUM_02_TC3,data=BEF.small)
BEF.small$res=ols$residuals
beta_ols=ols$coefficients
beta_ols_CI=ols$coefficients+1.96*
  t(t(sqrt(diag(summary(ols)$cov.unscaled))*summary(ols)$sigma))%*%c(-1,1)

res_OLS=resid(ols)

surf <- mba.surf(cbind(coords, res_OLS), no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
z.lim <- range(surf[[3]], na.rm=TRUE)
dev.new()
image.plot(surf, xaxs = "r", yaxs = "r", zlim=z.lim, main="LM residuals")

### mle ###
max.dist <- 0.25*max(iDist(BEF.small[,1:2]))
bins <- 20
vario <- variog(coords=BEF.small[,1:2], data=BEF.small$res, 
  uvec=(seq(0, max.dist, length=bins)))
plot(vario,pch=16)

mle <- likfit(coords=BEF.small[,1:2],data=BEF.small$logbio,
  trend = trend.spatial(~ELEV+SLOPE+SUM_02_TC1+SUM_02_TC2+SUM_02_TC3,
    BEF.small), ini.cov.pars=c(0.05,100),
  nugget = 0.05,cov.model="exponential",nospatial=TRUE)

mle$AIC
mle$BIC

mle$nospatial$AIC
mle$nospatial$BIC

## Predictions
BEF.shp <- readOGR("BEF-data/BEF_bound.shp")
shp2poly <- BEF.shp@polygons[[1]]@Polygons[[1]]@coords
BEF.poly <- as.matrix(shp2poly)
BEF.grids <- readGDAL("BEF-data/dem_slope_lolosptc_clip_60.img")

## Construct the prediction design matrix for the entire grid extent.
pred.covars <- cbind(BEF.grids[["band1"]], BEF.grids[["band2"]], BEF.grids[["band3"]], BEF.grids[["band4"]], BEF.grids[["band5"]])
colnames(pred.covars)=c("ELEV","SLOPE","SUM_02_TC1","SUM_02_TC2","SUM_02_TC3")
pred.covars=as.data.frame(pred.covars)
#pred.covars <- cbind(rep(1, nrow(pred.covars)), pred.covars)

## Extract the coordinates of the BEF bounding polygon vertices and use the pointsInPoly (spBayes) function to obtain the desired subset of the prediction design matrix and associated prediction coordinates (i.e., pixel centroids).
pred.coords <- SpatialPoints(BEF.grids)@coords
pointsInPolyOut <- pointsInPoly(BEF.poly, pred.coords)
pred.covars <- pred.covars[pointsInPolyOut,]
pred.coords <- pred.coords[pointsInPolyOut,]

krig_mlefit=krige.conv(coords=BEF.small[,1:2],data=BEF.small$logbio,
  locations=pred.coords,krige=krige.control(type.krige="OK",obj.model=mle,
    trend.d=trend.spatial(~ELEV+SLOPE+SUM_02_TC1+SUM_02_TC2+SUM_02_TC3,BEF.small),
    trend.l=trend.spatial(~ELEV+SLOPE+SUM_02_TC1+SUM_02_TC2+SUM_02_TC3, pred.covars)))

pred_mlefit=krig_mlefit$predict
var_mlefit=krig_mlefit$krige.var

## Mapping the predicted values
surf <- mba.surf(cbind(coords, log.bio), no.X=x.res, no.Y=x.res, extend=TRUE, sp=TRUE)$xyz.est
surf <- surf [!is.na((over(surf, BEF.shp)))[,1],]
surf <- as.image.SpatialGridDataFrame(surf)
z.lim <- range(surf[["z"]], na.rm=TRUE)

pred.grid <- as.data.frame(list(pred.coords, pred.mu=pred_mlefit, pred.sd=var_mlefit))
coordinates(pred.grid) = c("x", "y")
gridded(pred.grid) <- TRUE
pred.mu.image <- as.image.SpatialGridDataFrame(pred.grid["pred.mu"])
pred.sd.image <- as.image.SpatialGridDataFrame(pred.grid["pred.sd"])

dev.new()
image.plot(surf, axes=TRUE, zlim=z.lim, col=tim.colors(25), xaxs = "r", yaxs = "r", main="Log metric tons of biomass")
plot(BEF.shp, add=TRUE)
points(coords)
dev.new()
image.plot(pred.mu.image, axes=TRUE, zlim=z.lim, col=tim.colors(25), xaxs = "r", yaxs = "r", main="Log metric tons of biomass")
plot(BEF.shp, add=TRUE)
dev.new()
image.plot(pred.sd.image, col=tim.colors(25), xaxs = "r", yaxs = "r", main="Mean predicted log metric tons of biomass")
plot(BEF.shp, add=TRUE)
points(coords,cex=1)

beta_mle=mle$beta
beta_mle_CI=mle$beta+1.96*t(t(sqrt(diag(mle$beta.var))))%*%c(-1,1)






