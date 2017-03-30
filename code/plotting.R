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

### code for visualizing point referenced spatial data ###

### plotting on the xy-plane 
data1=read.csv("../data/dataset1.csv")

##colors <- brewer.pal(5,"Spectral") ## other choices for color palettes 
##colors=heat.colors(11)
#colors=c("#800000","#FFFFFFFF","#000080")
colors=colorRampPalette(c("midnightblue", "cyan","yellow","red"))(5)
quants=classIntervals(data1$y,n=5,style="quantile")
sc=scale_colour_gradientn(colours = colors,values=rescale(quants$brks,c(0,1)),name="y")

#### data locations ####
dev.new()
plot(data1$sx,data1$sy,xlab="",ylab="",main="Data locations")

#### color coded data locations ####
p0 = ggplot(data=data1,aes(x=sx,y=sy))+sc+labs(x = "",y = "")
p1=p0+geom_point(data=data1,aes(x=sx,y=sy,color=y))+guides(fill=guide_colorbar(title="y"))
plot(p1)

### function to create an equispaced sequence based on the range of a vector 
frange=function(x,n) seq(min(x),max(x),length=n)

#### interpolation ####
intobj1=interp(data1$sx,data1$sy,data1$y,
  xo=frange(data1$sx,200),yo=frange(data1$sy,200))
interpdata1=cbind(expand.grid(intobj1$x,intobj1$y),as.vector(intobj1$z))
colnames(interpdata1)=c("sx","sy","y")

#### plotting the interpolated surface using akima package ####
p1interp=p0+geom_raster(data=interpdata1,aes(x=sx,y=sy,fill=y))+
  scale_fill_gradientn(colours=colors,values=rescale(quants$brks,c(0,1)),name="y")
plot(p1interp)

###################################################
### Similar plots for WEF data (spBayes package) without using ggplot2
###################################################
rm(list=ls())

data(WEF.dat)
head(WEF.dat)

WEF.dat <- WEF.dat[!apply(WEF.dat[,c("East_m","North_m","DBH_cm","Tree_height_m","ELEV_m")], 1, function(x)any(is.na(x))),]

### diameter at breast height for the trees
DBH <- WEF.dat$DBH_cm

### size coded locations ###
coords <- as.matrix(WEF.dat[,c("East_m","North_m")])
plot(coords, pch=1, cex=sqrt(DBH)/10, col="darkgreen", xlab="Easting (m)", ylab="Northing (m)")
leg.vals <- round(quantile(DBH),0)
legend("topleft", pch=1, legend=leg.vals, col="darkgreen", pt.cex=sqrt(leg.vals)/10, bty="n", title="DBH (cm)")

### setting the colors ####
col.br <- colorRampPalette(c("midnightblue", "cyan", "yellow", "red"))
col.pal <- col.br(5)

### DBH quantile based color coding of the locations
quant <- classIntervals(DBH, n=5, style="quantile")

## plot showing how the data was divided based on quantiles
plot(quant, pal=col.pal, xlab="DBH", main="Quantile") 

quant.col <- findColours(quant, col.pal)

plot(coords, col=quant.col, pch=19, cex=0.5, main="Quantile", xlab="Easting (m)", ylab="Northing (m)")
legend("topleft", fill=attr(quant.col, "palette"), 
       legend=names(attr(quant.col, "table")), bty="n")

### DBH classes based color coding of the locations
fixed <- classIntervals(DBH, n=4, style="fixed", fixedBreaks=c(0,12.7,30.48,60,max(DBH)+1))

fixed.col <- findColours(fixed, col.pal)

plot(coords, col=fixed.col, pch=19, cex=0.5, main="Forestry tree size classes", xlab="Easting (m)", ylab="Northing (m)")
legend("topleft", fill=attr(fixed.col, "palette"), legend=c("sapling","poletimber","sawtimber","large sawtimber"), bty="n")

###################################################
### interpolation using mBA package
###################################################
x.res <- 100; y.res <- 100

surf <- mba.surf(cbind(coords, DBH), no.X=x.res, no.Y=y.res, h=5, m=2, extend=FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))




