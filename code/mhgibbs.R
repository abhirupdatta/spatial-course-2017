########## Example: Gibbs' sampling ###########
library(mvtnorm)
n=100
set.seed(1)
X=rnorm(n,0,1)
Y=2+5*X+rnorm(n,0,1)
gamma0=100
gamma1=10

N=10000
Nb=5001

beta0_sample_gibbs=beta1_sample_gibbs=rep(0,N)
beta1init=10
for(i in 1:N){
    beta0_sample_gibbs[i]=rnorm(1,sum(Y-X*beta1init)/(n+1/gamma0),sqrt(1/(n+1/gamma0)))
    beta0init=beta0_sample_gibbs[i]
    beta1_sample_gibbs[i]=rnorm(1,sum(X*(Y-beta0init))/(sum(X^2)+1/gamma1),sqrt(1/(sum(X^2)+1/gamma1)))
    beta1init=beta1_sample_gibbs[i]
}

Xmat=cbind(rep(1,n),X)
betavar=solve(t(Xmat)%*%Xmat+diag(1/c(gamma0,gamma1)))
beta_sample=rmvnorm(N-Nb+1, betavar%*%t(Xmat)%*%Y, betavar)

dev.new()
plot(density(beta_sample[,1]),xlab="theta",ylab="p(theta)",main="Posterior density beta_0")
lines(density(beta0_sample_gibbs[Nb:N]), col="red")
legend("topright",legend=c("Direct","Gibbs'"), col=c("black","red"), lwd=rep(1,2))

dev.new()
plot(density(beta_sample[,2]),xlab="theta",ylab="p(theta)",main="Posterior density beta_1")
lines(density(beta1_sample_gibbs[Nb:N]), col="red")
legend("topright",legend=c("Direct","Gibbs'"), col=c("black","red"), lwd=rep(1,2))

############ Example: Metropolis Algorithm #############
n=100
theta=3
sigmasq=1
set.seed(1)
Y=rnorm(n,theta,sqrt(sigmasq))

tausq=5
mu=0
posterior_loglikelihood=function(theta,Y,n,sigmasq,mu,tausq){
	-0.5*n*(mean(Y)-theta)^2/sigmasq-0.5*(theta-mu)^2/tausq
	}
N=10000
Nb=5001
lambda=0.1
theta_sample_metropolis=flag=rep(0,N)
thetainit=0
for(i in 1:N){
	thetastar=rnorm(1,thetainit,sqrt(lambda))
	logr=posterior_loglikelihood(thetastar,Y,n,sigmasq,mu,tausq)-posterior_loglikelihood(thetainit,Y,n,sigmasq,mu,tausq)
	logu=log(runif(1,0,1))
	if(logu <= logr){
		theta_sample_metropolis[i]=thetastar
		flag[i]=1
		}	else	
			{
				theta_sample_metropolis[i]=thetainit
			}
	thetainit=theta_sample_metropolis[i]
	}

acceptance_ratio=sum(flag)/length(flag)
acceptance_ratio

### second chain ###
theta_sample_metropolis2=flag2=rep(0,N)
thetainit2=10
for(i in 1:N){
	thetastar=rnorm(1,thetainit,sqrt(lambda))
	logr=posterior_loglikelihood(thetastar,Y,n,sigmasq,mu,tausq)-posterior_loglikelihood(thetainit2,Y,n,sigmasq,mu,tausq)
	logu=log(runif(1,0,1))
	if(logu <= logr){
		theta_sample_metropolis2[i]=thetastar
		flag2[i]=1
		}	else	
			{
				theta_sample_metropolis2[i]=thetainit2
			}
	thetainit2=theta_sample_metropolis2[i]
	}

dev.new()
plot(theta_sample_metropolis, type="l", col="red", xlab="Iteration i", ylab="theta_i", main="Trace plots")
lines(theta_sample_metropolis2, col="blue")

### direct sample ###
theta_sample=rnorm(N-Nb+1,(n*mean(Y)/sigmasq+mu/tausq)/(n/sigmasq+1/tausq),sqrt(1/(n/sigmasq+1/tausq)))

dev.new()
plot(density(theta_sample),xlab="theta",ylab="p(theta)",main="Posterior density")
lines(density(theta_sample_metropolis[Nb:N]), col="red")
lines(density(theta_sample_metropolis2[Nb:N]), col="blue")
legend("topright",legend=c("Direct","MA 1","MA 2"), col=c("black","red","blue"), lwd=rep(1,3))

####### Tuning ########
lambda=10
theta_sample_metropolis_low=flaglow=rep(0,N)
thetainit=0
for(i in 1:N){
	thetastar=rnorm(1,thetainit,sqrt(lambda))
	logr=posterior_loglikelihood(thetastar,Y,n,sigmasq,mu,tausq)-posterior_loglikelihood(thetainit,Y,n,sigmasq,mu,tausq)
	logu=log(runif(1,0,1))
	if(logu <= logr){
		theta_sample_metropolis_low[i]=thetastar
		flaglow[i]=1
		}	else	
			{
				theta_sample_metropolis_low[i]=thetainit
			}
	thetainit=theta_sample_metropolis_low[i]
	}

acceptance_ratio_low=sum(flaglow)/length(flaglow)
acceptance_ratio_low

plot(theta_sample_metropolis_low, type="l", col="red", xlab="Iteration i", ylab="theta_i", main="lambda=10")

lambda=0.0001
theta_sample_metropolis_high=flaghigh=rep(0,N)
thetainit=0
for(i in 1:N){
	thetastar=rnorm(1,thetainit,sqrt(lambda))
	logr=posterior_loglikelihood(thetastar,Y,n,sigmasq,mu,tausq)-posterior_loglikelihood(thetainit,Y,n,sigmasq,mu,tausq)
	logu=log(runif(1,0,1))
	if(logu <= logr){
		theta_sample_metropolis_high[i]=thetastar
		flaghigh[i]=1
		}	else	
			{
				theta_sample_metropolis_high[i]=thetainit
			}
	thetainit=theta_sample_metropolis_high[i]
	}

acceptance_ratio_high=sum(flaghigh)/length(flaghigh)
acceptance_ratio_high

plot(theta_sample_metropolis_high, type="l", col="red", xlab="Iteration i", ylab="theta_i", main="lambda=0.0001")

######  Example 2: Jacobian ######
library(MCMCpack)
n=100
sigmasq=4
set.seed(1)
Y=rnorm(n,0,sqrt(sigmasq))

alpha=2
beta=1
posterior_loglikelihood=function(sigmasq,Y,n,alpha,beta){
	-(alpha+n/2+1)*log(sigmasq)-(sum(Y^2)/2+beta)/sigmasq
	}
N=10000
Nb=5001
lambda=0.1
sigmasq_sample_metropolis=flag=rep(0,N)
sigmasqinit=1
for(i in 1:N){
	sigmasqstar=exp(rnorm(1,log(sigmasqinit),sqrt(lambda)))
	logr=posterior_loglikelihood(sigmasqstar,Y,n,alpha,beta)+log(sigmasqstar)-posterior_loglikelihood(sigmasqinit,Y,n,alpha,beta)-log(sigmasqinit)
	logu=log(runif(1,0,1))
	if(logu <= logr){
		sigmasq_sample_metropolis[i]=sigmasqstar
		flag[i]=1
		}	else	
			{
				sigmasq_sample_metropolis[i]=sigmasqinit
			}
	sigmasqinit=sigmasq_sample_metropolis[i]
	}

acceptance_ratio=sum(flag)/length(flag)
acceptance_ratio

sigmasq_sample=rinvgamma(N-Nb+1,alpha+n/2,beta+sum(Y^2)/2)

dev.new()
plot(density(sigmasq_sample),xlab="theta",ylab="p(theta)",main="Posterior density")
lines(density(sigmasq_sample_metropolis[Nb:N]), col="red")
legend("topright",legend=c("Direct","MA"), col=c("black","red"), lwd=rep(1,2))
