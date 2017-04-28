
## R code from BCG Second Edition, pages 140--150

## Load the libraries to be used

library(spBayes)
library(MBA)
library(geoR)
library(fields)
library(sp)
library(maptools)
library(rgdal)
library(classInt)
library(lattice)


## BEF data ###
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

p <- 6 ## This is the number of columns in the design matrix
## Set the prior mean and precision for the regression
beta.prior.mean <- as.matrix(rep(0, times=p))
beta.prior.precision <- matrix(0, nrow=p, ncol=p)

## For use with bayesGeostatExact, do the following
phi <- 0.014 ## Set the spatial range (from the variogram)
alpha <- 0.016/0.08 ## Set the nugget/partial-sill ratio
sigma.sq.prior.shape <- 2.0 ## Set IG shape for sigma.sq (partial sill)
sigma.sq.prior.rate <- 0.08 ## Set IG scale for sigma.sq (partial sill)

## Run bayesGeostatExact to deliver exact posterior samples
sp.exact <- bayesGeostatExact(
log.bio~ELEV+SLOPE+SUM_02_TC1+SUM_02_TC2+SUM_02_TC3,
data=BEF.dat, coords=coords, n.samples=1000,
beta.prior.mean=beta.prior.mean,
beta.prior.precision=beta.prior.precision,
cov.model="exponential",
phi=phi, alpha=alpha,
sigma.sq.prior.shape=sigma.sq.prior.shape,
sigma.sq.prior.rate=sigma.sq.prior.rate,
sp.effects=FALSE)

##Produce the posterior summaries
round(summary(sp.exact$p.samples)$quantiles,3)


## Run spLM to deliver MCMC samples from marginal posterior distributions
n.samples <- 1000
bef.sp <- spLM(log.bio~ELEV+SLOPE+SUM_02_TC1+SUM_02_TC2+SUM_02_TC3,
data=BEF.dat, coords=coords, starting=list("phi"=3/200,"sigma.sq"=0.08,
"tau.sq"=0.02), tuning=list("phi"=0.1, "sigma.sq"=0.05, "tau.sq"=0.05),
               priors=list("phi.Unif"=c(3/1500, 3/50), "sigma.sq.IG"=c(2, 0.08),
                 "tau.sq.IG"=c(2, 0.02)), cov.model="exponential",n.samples=n.samples)

round(summary(mcmc(bef.sp$p.theta.samples))$quantiles,3)

## Recover spatial residuals using spRecover
burn.in <- floor(0.75*n.samples)
bef.sp <- spRecover(bef.sp, start=burn.in, thin=2)

## The posterior samples of the regression coefficients and the spatial effects can then be obtained as
beta.samples = bef.sp$p.beta.recover.samples
w.samples = bef.sp$p.w.recover.samples

## Obtain trace plots for regression coefficients
dev.new()
par(mfrow=c(3,2))
plot(beta.samples, auto.layout=TRUE, density=FALSE)

round(summary(mcmc(bef.sp$p.beta.recover.samples))$quantiles,3)

## Obtain posterior means and sd's of spatial residuals for each location
w.hat.mu <- apply(w.samples,1,mean)
w.hat.sd <- apply(w.samples,1,sd)

## Obtain OLS residuals
lm.bio = lm(log.bio~ELEV+SLOPE+SUM_02_TC1+SUM_02_TC2+SUM_02_TC3, data=BEF.dat)
bio.resid = resid(lm.bio)

## Plot the spatial residual mean surface and a map of sd's
par(mfrow=c(1,2))
surf <- mba.surf(cbind(coords, bio.resid), no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
z.lim <- range(surf[[3]], na.rm=TRUE)
image.plot(surf, xaxs = "r", yaxs = "r", zlim=z.lim, main="LM residuals")
surf <- mba.surf(cbind(coords, w.hat.mu), no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", zlim=z.lim, main="Mean spatial effects (w(s))")

## Predictions
BEF.shp <- readOGR("BEF-data/BEF_bound.shp")
shp2poly <- BEF.shp@polygons[[1]]@Polygons[[1]]@coords
BEF.poly <- as.matrix(shp2poly)
BEF.grids <- readGDAL("BEF-data/dem_slope_lolosptc_clip_60.img")

## Construct the prediction design matrix for the entire grid extent.
pred.covars <- cbind(BEF.grids[["band1"]], BEF.grids[["band2"]], BEF.grids[["band3"]], BEF.grids[["band4"]], BEF.grids[["band5"]])
pred.covars <- cbind(rep(1, nrow(pred.covars)), pred.covars)


## Extract the coordinates of the BEF bounding polygon vertices and use the pointsInPoly (spBayes) function to obtain the desired subset of the prediction design matrix and associated prediction coordinates (i.e., pixel centroids).
pred.coords <- SpatialPoints(BEF.grids)@coords
pointsInPolyOut <- pointsInPoly(BEF.poly, pred.coords)
pred.covars <- pred.covars[pointsInPolyOut,]
pred.coords <- pred.coords[pointsInPolyOut,]

bef.bio.pred <- spPredict(bef.sp, start=burn.in, thin=2, pred.coords=pred.coords, pred.covars=pred.covars)

## Mapping the predicted values
bef.bio.pred.mu = apply(bef.bio.pred$p.y.predictive.samples,1,mean)
bef.bio.pred.sd = apply(bef.bio.pred$p.y.predictive.samples,1,sd)

surf <- mba.surf(cbind(coords, log.bio), no.X=x.res, no.Y=x.res, extend=TRUE, sp=TRUE)$xyz.est
#surf <- surf [!is.na(over(surf, BEF.shp)),]
surf <- surf [!is.na((over(surf, BEF.shp)))[,1],]
surf <- as.image.SpatialGridDataFrame(surf)
z.lim <- range(surf[["z"]], na.rm=TRUE)

pred.grid <- as.data.frame(list(pred.coords, pred.mu=bef.bio.pred.mu, pred.sd=bef.bio.pred.sd))
coordinates(pred.grid) = c("x", "y")
gridded(pred.grid) <- TRUE
pred.mu.image <- as.image.SpatialGridDataFrame(pred.grid["pred.mu"])

par(mfrow=c(1,2))
image.plot(surf, axes=TRUE, zlim=z.lim, col=tim.colors(25), xaxs = "r", yaxs = "r", main="Log metric tons of biomass")
plot(BEF.shp, add=TRUE)
image.plot(pred.mu.image, zlim=z.lim, col=tim.colors(25), xaxs = "r", yaxs = "r", main="Mean predicted log metric tons of biomass")
plot(BEF.shp, add=TRUE)