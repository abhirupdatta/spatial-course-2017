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
library("mvtnorm")
library("MCMCpack")

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

n=nrow(data3)

X=cbind(1,data3$x)
y=data3$y

S=data3[,c("sx","sy")]
dmat=as.matrix(dist(S))


#### Chain length and burn in ####
N=10000
Nb=5001

#### true values and posterior distribution and quantiles ####
truesigs=0.25
truetaus=0.01
truebeta=c(0.2,-0.3)
truephi=2

####### Gibbs sampler (assuming phi is known) ##########

### linear regression and empirical variogram to choose phi  ###
ols=lm(y~x,data=data3)
data3$res=ols$residuals

max.dist <- 0.75*max(iDist(data3[,1:2]))
bins <- 20
vario3 <- variog(coords=data3[,1:2], data=data3$res, 
  uvec=(seq(0, max.dist, length=bins)))
vfit3 <-variofit(vario3, ini.cov.pars=c(0.1,1), ##sigma^2 and 1/phi 
  cov.model="exponential", minimisation.function="optim",
  nugget=0.01, weights="equal")
phi=1/vfit3$cov.pars[2]

R=exp(-phi*dmat)
Rinv=chol2inv(chol(R))
gram=chol2inv(chol(t(X)%*%X))
M=gram%*%t(X)

### creating matrices and vectors to store the samples ###
betamat=matrix(0,N,2)
wmat=matrix(0,N,n)
sigsvec=tausvec=rep(0,N)

### initial values ###
beta=rep(0,2)
sigs=taus=1

### hyper-params ###
asig=1
atau=0.5
bsig=btau=2
res=y-X%*%beta

### Gibbs sampler ###
set.seed(1)
for(i in 1:N){
    Vw=chol2inv(chol(Rinv/sigs+diag(n)/taus))
    w=as.vector(rmvnorm(1,Vw%*%res/taus,Vw))
    beta=as.vector(rmvnorm(1,M%*%(y-w),taus*gram))
    res=y-X%*%beta
    sigs=rinvgamma(1,asig+n/2,bsig+t(w)%*%Rinv%*%w/2)
    taus=rinvgamma(1,atau+n/2,btau+sum((res-y)^2)/2)
    betamat[i,]=beta
    wmat[i,]=w
    sigsvec[i]=sigs
    tausvec[i]=taus
    if(i %% 100 == 0) print(i)
    }

MCMCmat=cbind(betamat,wmat,sigsvec,tausvec)
colnames(MCMCmat)= c(paste0("beta",1:length(beta)),paste0("w",1:length(w)),"sigs","taus")
#write.csv(MCMCmat,"data3MCMCmat.csv",row.names=F,quote=F)
samples=MCMCmat

### loading the pre-saved samples as the MCMC takes a while
### samples=read.csv("data3MCMCmat.csv")

pbsample <- samples[Nb:N,] ## post burn-in samples

qtls <- apply(pbsample,2,quantile,c(0.025,0.5,0.95))
qtls[,c(1:2,n+3,n+4)]

dev.new()
plot(density(pbsample[,1]),lwd=2,main="beta0",xlab="",ylab="")
abline(v=qtls[,1],lwd=2)
abline(v=truebeta[1],col="red",lwd=2)

dev.new()
plot(density(pbsample[,2]),lwd=2,main="beta1",xlab="",ylab="")
abline(v=qtls[,2],lwd=2)
abline(v=truebeta[2],col="red",lwd=2)

dev.new()
plot(density(pbsample[,n+3]),lwd=2,main="sigs",xlab="",ylab="")
abline(v=qtls[,n+3],lwd=2)
abline(v=truesigs,col="red",lwd=2)

dev.new()
plot(density(pbsample[,n+4]),lwd=2,main="taus",xlab="",ylab="")
abline(v=qtls[,n+4],lwd=2)
abline(v=truetaus,col="red",lwd=2)

what=qtls[2,3:(n+2)]
data3$what=what
myplot(data3,"what")

### plotting kriged surface ###
set.seed(1)
subsample=pbsample[sample(1:nrow(pbsample),100),]

xo=yo=seq(0,1,0.02)
so=expand.grid(xo,yo)

Do=rdist(so,S)

### we can calculate this as phi is fixed here 
c=exp(-phi*Do)
weights=c%*%Rinv

### we will use only 100 posterior samples 
subsample=as.matrix(subsample)
wpredmean=subsample[,3:(n+2)]%*%t(weights)

Xo=cbind(1,0.5*sin(10*so[,1]*so[,2])+1*(0.5-so[,1])^2)

### kriging using composition sampling wo|w,params,y then yo|wo,w,params,y
set.seed(1)
predmat=sapply(1:nrow(so),function(i,c,weights,wpredmean,subsample,Xo){
  wpredvar=pbsample[,n+3]*(1-sum(c[i,]*weights[i,]))
  wo=rnorm(100,wpredmean[,i],sqrt(wpredvar))
  yo=as.vector(subsample[,1:2]%*%Xo[i,]) + wo + rnorm(100,rep(0,100),sqrt(subsample[,n+4]))
  yo
},c,weights,wpredmean,subsample,Xo)

surface_krig_tab=cbind(so,Xo[,2],apply(predmat,2,median),apply(predmat,2,var))
colnames(surface_krig_tab)=c("sx","sy","x","yhat","vyhat")

myplot(surface_krig_tab,"yhat")
myplot(surface_krig_tab,"vyhat")

############ MH for the marginalized model using Nimble ##############
library(nimble)

gpcov <- nimbleFunction(
  run = function(dmat=double(2),phi=double(0),n=double(0)){
    returnType(double(2))
    M=matrix(0,n,n)
    for(i in 1:n) 
      for(j in 1:n)
        M[i,j]=exp( -phi * dmat[i,j])
    #M=exp( -ph * dmat)
    return(M)
  })


bayesSpatialLM <- nimbleCode({
  sigs <- 1/invsigmasq
  taus <- 1/invtausq
  G[,] <- gpcov(dmat[,],phi,n)
  V[,] <- sigs*G[,]+ taus * I[,]
  C[,] <- chol(V[,])
  mu[] <- X[,]%*%beta[]
  y[] ~ dmnorm(mu[], cholesky=C[,], prec_param=0)  ## marginalized model
  
  ### parameter priors ###
  beta[] ~ dmnorm(mub[], cholesky=Cb[,], prec_param=0)
  invsigmasq ~ dgamma(2, rate=1)
  invtausq ~ dgamma(0.01, 0.01)
  phi ~  dunif(0,10)
})

constants <- list(n = n, dmat=dmat,
  ones = rep(1,n), I=diag(n), Cb=1000*diag(2), mub=rep(0,2)
)

data <- list(y = y, X=X)

dimensions = list(G = c(n, n),
  V = c(n, n),
  C = c(n, n),
  Cb = c(2, 2),
  mub = c(2),
  mu = c(n),
  ones = c(n),
  I = c(n,n),
  beta = c(2),
  y = c(n),
  X= c(n,2),
  G = c(n,n))

model <- nimbleModel(code=bayesSpatialLM, constants=constants, data=data,
  dimensions = dimensions,check = FALSE)

Cmodel <- compileNimble(model)  ## compiling model 
modelconf <- configureMCMC(Cmodel,print=TRUE) 
modelconf$addMonitors(c('beta','sigs','taus','phi'))

Rmcmc <- buildMCMC(modelconf)
nimbleOptions(showCompilerOutput = FALSE) 
Cmcmc <- compileNimble(Rmcmc, project = model) ## compiling MCMC
set.seed(1)
Cmcmc$run(100) ## running MCMC

samples <- as.matrix(Cmcmc$mvSamples)
#write.csv(samples,"data3nimble.csv",quote=F,row.names=F)

### loading the pre-saved samples
### samples = read.csv("../data/data3nimble.csv")

pbsample <- samples[Nb:N,] ## post burn-in samples

qtls <- apply(pbsample,2,quantile,c(0.025,0.5,0.95))
qtls

#### true values and posterior distribution and quantiles ####
dev.new()
plot(density(pbsample[,"beta.1."]),lwd=2,main="beta0",xlab="",ylab="")
abline(v=qtls[,"beta.1."],lwd=2)
abline(v=truebeta[1],col="red",lwd=2)

dev.new()
plot(density(pbsample[,"beta.2."]),lwd=2,main="beta1",xlab="",ylab="")
abline(v=qtls[,"beta.2."],lwd=2)
abline(v=truebeta[2],col="red",lwd=2)

dev.new()
plot(density(pbsample[,"sigs"]),lwd=2,main="sigs",xlab="",ylab="")
abline(v=qtls[,"sigs"],lwd=2)
abline(v=truesigs,col="red",lwd=2)

dev.new()
plot(density(pbsample[,"taus"]),lwd=2,main="taus",xlab="",ylab="")
abline(v=qtls[,"taus"],lwd=2)
abline(v=truetaus,col="red",lwd=2)

dev.new()
plot(density(pbsample[,"phi"]),lwd=2,main="phi",xlab="",ylab="")
abline(v=qtls[,"phi"],lwd=2)
abline(v=truephi,col="red",lwd=2)

#### Posterior recovery of w ####
set.seed(1)
subsample=pbsample[sample(1:nrow(pbsample),100),]
wsamples=apply(subsample,1,function(params,y,X,dmat){
  beta=params[1:2]
  taus=params[7]
  sigs=params[6]
  phi=params[5]
  Rinv=chol2inv(chol(exp(-phi*dmat)))
  Vw=chol2inv(chol(Rinv/sigs+diag(n)/taus))
  w=as.vector(rmvnorm(1,Vw%*%(y-X%*%beta)/taus,Vw))
  },y,X,dmat  )

what=apply(wsamples,1,median)
data3$what=what
myplot(data3,"what")
