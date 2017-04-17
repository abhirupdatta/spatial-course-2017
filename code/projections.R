### projections ###

library(sp)
library(rgdal)

###### geodesic naiveEuclid ####
lltoxy=function(loc)
{
  i=nchar(loc)-4
  part1=strtoi(substring(loc,1,i))
  part2=strtoi(substring(loc,i+2,i+3))
  part3=substring(loc,i+4,i+4)
  deg=(part1+part2/100)*pi/180
  if(((part3=="W") || (part3=="S")) && (part1<180)) deg=-deg
  deg
}

naiveEuclid=function(lat1,long1,lat2,long2)
{
  lambda1=lltoxy(long1)
  lambda2=lltoxy(long2)
  theta1=lltoxy(lat1)
  theta2=lltoxy(lat2)
  naivedist=6371*sqrt((lambda1-lambda2)^2+(theta1-theta2)^2)
  naivedist
}

naiveEuclid("41.88N","87.63W","44.89N","93.22W")
naiveEuclid("40.78N","73.97W","29.98N","90.25W")


geodesic=function(lat1,long1,lat2,long2)
{
  lambda1=lltoxy(long1)
  lambda2=lltoxy(long2)
  theta1=lltoxy(lat1)
  theta2=lltoxy(lat2)
  geo=6371*acos(round(sin(theta1)*sin(theta2)+cos(theta1)*cos(theta2)*cos(lambda2-lambda1),5))
  geo
}

geodesic("41.88N","87.63W","44.89N","93.22W")
geodesic("40.78N","73.97W","29.98N","90.25W")


### UTM ###
LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}

# Example
y<-c( 41.8,44.89,40.78,29.98)
x<--c( 87.63,93.22,73.97,90.25)
s=LongLatToUTM(x,y,17)
row.names(s)=c("Chicago","Minneapolis","New York","New Orleans")
round(dist(s)/1000,2)
